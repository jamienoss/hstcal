#include <stdio.h>
#include <math.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"
#include "pcte.h"
#include "acs.h"

/*This is the workhorse subroutine; it simulates the readout
  of one column currentColumn[].

  JDIM == RAZ_ROWS
  WDIM == TRAPS  Ws is the input traps number < 999999
  NITs == cte_pars->n_par

  These are already in the parameter structure CTEParams
  int     Ws              the number of traps < 999999
  float     q_w[TRAPS];     the run of charge with level  == qlevq_data
  float   dpde_w[TRAPS];  the run of charge loss with level == dpdew_data
  float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel == rprof fits image
  float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)  == cprof fits image

  W = wcol_data = trap id

  q_w[TRAP] = qlev_q from QPROF  traps as function of packet size = cte->qlevq_data[TRAP]

  currentColumn (read), pixf(cteff) are passed and are 1d arrays which have values for a particular column

  the ttrap reference to the image array has to be -1 for C
  */

int sim_colreadout_l(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows)
{
    extern int status;

    double pixAdd;
    double pixAddMore;//rethink name?
    double pixRemove;
    double pix;
    double releasedFlux;
    double trappedFlux;
    int nTransfersFromTrap;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE*/
    //Look into whether this really has to be computed each iteration?
    //Since this is simulating the readout and thus moving pixels down and out, pmax can only get smaller with
    //each pixel transfer, never greater.
    double pmax = 10;
    for (unsigned i = 0; i < nRows; ++i)
        pmax = (pixelColumn[i] < pmax) ? pmax : pixelColumn[i]; //check assembly (before & after op) to see if actually implemented differently

    //Find highest charge trap to not exceed i.e. map pmax to an index
    unsigned maxChargeTrapIndex = cte->cte_traps-1;
    for (int w = cte->cte_traps-1; w >= 0; --w)//go up or down? (if swap, change below condition)
    {
        if (cte->qlevq_data[w] <= pmax)//is any of this even needed or can we just directly map?
        {
            maxChargeTrapIndex = w;
            break;
        }
    }

    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    for (int w = maxChargeTrapIndex; w >= 0; --w)
    {
        trappedFlux = 0;
        nTransfersFromTrap = cte->cte_len; /*for referencing the image at 0*/
        releasedFlux = 0;

        /*GO UP THE COLUMN PIXEL BY PIXEL*/
        for(unsigned i = 0; i < nRows; ++i)
        {
            pix = pixelColumn[i];
            //What is this doing - for pix_1 >= cte->qlevq_data[w] - 1.???
            if ( (nTransfersFromTrap < cte->cte_len) || ( pix >= cte->qlevq_data[w] - 1. ) )
            {
                if (pixelColumn[i] >= 0 )//seems a shame to need check this every iteration
                {
                    pix = pixelColumn[i] + releasedFlux; /*shuffle charge in*/
                    releasedFlux = pix - floor(pix); /*carry the charge remainder*/
                    pix = floor(pix); /*reset pixel*/
                }

                /*HAPPENS AFTER FIRST PASS*/
                /*SHUFFLE CHARGE IN*/
                //move out of loop to separate instance?
                if (i > 0)
                {
                    if (pixf[i] < pixf[i-1])
                        trappedFlux *= (pixf[i] / pixf[i-1]);
                }

                /*RELEASE THE CHARGE*/
                pixAdd=0.0;
                if (nTransfersFromTrap < cte->cte_len)
                {
                    pixAdd = rprof->data[w*rprof->ny + nTransfersFromTrap] * trappedFlux;
                    ++nTransfersFromTrap;
                }

                pixAddMore = 0;
                pixRemove = 0;
                if ( pix >= cte->qlevq_data[w])
                {
                    pixRemove =  cte->dpdew_data[w] / cte->n_par * pixf[i];  /*dpdew is 1 in file */
                    if (nTransfersFromTrap < cte->cte_len)
                        pixAddMore = cprof->data[w*cprof->ny + nTransfersFromTrap-1] * trappedFlux; //ttrap-1 may not be the same index as ref'd in rprof???
                    nTransfersFromTrap=0;
                    trappedFlux=pixRemove;
                }

                pixelColumn[i] += pixAdd + pixAddMore - pixRemove;
            } //replaces trap continue
        } //end for i
    } //end for w
    return(status);
}

Bool correctCROverSubtraction(double * const pix_ctef, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold, const int doCorrectionFlag)
{
    /*LOOK FOR AND DOWNSCALE THE CTE MODEL IF WE FIND
      THE TELL-TALE SIGN OF READOUT CRS BEING OVERSUBTRACTED;
      IF WE FIND ANY THEN GO BACK UP AND RERUN THIS COLUMN

      THE WFC3 UVIS MODEL SEARCHES FOR OVERSUBTRACTED TRAILS.
      WHICH ARE  DEFINED AS EITHER:

          - A SINGLE PIXEL VALUE BELOW -10E-
          - TWO CONSECUTIVE PIXELS TOTALING -12 E-
          - THREE TOTALLING -15 E-

      WHEN WE DETECT SUCH AN OVER-SUBTRACTED TAIL, WE ITERATIVELY REDUCE
      THE LOCAL CTE SCALING BY 25% UNTIL THE TRAIL IS
      NO LONGER NEGATIVE  THIS DOES NOT IDENTIFY ALL READOUT-CRS, BUT IT DOES
      DEAL WITH MANY OF THEM. FOR IMAGES THAT HAVE BACKGROUND GREATER THAN 10 OR SO,
      THIS WILL STILL END UP OVERSUBTRACTING CRS A BIT, SINCE WE ALLOW
      THEIR TRAILS TO BE SUBTRACTED DOWN TO -10 RATHER THAN 0.
    */

    if (!doCorrectionFlag || !pix_model || !pix_observed || !pix_ctef)
        return False;

    Bool redo = False;
    for (unsigned j = 10; j < nRows-2; ++j)
    {
        if ( (( threshHold > pix_model[j] ) &&
                    ( threshHold > (pix_model[j] - pix_observed[j]))) ||

                (((pix_model[j] + pix_model[j+1]) < -12.) &&
                 (pix_model[j] + pix_model[j+1] - pix_observed[j] - pix_observed[j+1] < -12.)) ||

                (((pix_model[j] + pix_model[j+1] + pix_model[j+2]) < -15.) &&
                 ((pix_model[j] + pix_model[j+1] + pix_model[j+2] -pix_observed[j] -
                   pix_observed[j+1] - pix_observed[j+2]) <-15.))  )
        {
            redo = True;
            unsigned jmax = j;

            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
            for (unsigned jj = j-10; jj <= j; ++jj)
            {
                if ( (pix_model[jj] - pix_observed[jj]) > (pix_model[jmax] - pix_observed[jmax]) )
                    jmax = jj;
                    //should there be a break here?
            }
            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
            for (unsigned jj = jmax; jj <= j; ++jj)
                pix_ctef[jj] *= 0.75;
        }
    } /*end for  j*/
    return redo;
}
