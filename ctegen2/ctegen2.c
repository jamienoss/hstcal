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

int sim_colreadout_l(double *currentColumn, const double *pixf, const CTEParams *cte, const unsigned nRows, const Bool dampen)
{
    extern int status;

    //for (unsigned nthIter = 0; nthIter < nIterations; ++nthIter)//look into moving this further down
    //{
        double padd_3=0.0;
        double prem_3=0.0;
        double padd_2=0.0;
        double fcarry=0.0;
        double pix_1=0.0;
        double ftrap=0.0;
        int ttrap=0;

        /*from the reference table*/
        const FloatHdrData * rprof = cte->rprof;
        const FloatHdrData * cprof = cte->cprof;
        const double rnAmp2 = cte->rn_amp*cte->rn_amp;

        /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
          PMAX SHOULD ALWAYS BE POSITIVE HERE*/
        double pmax = 10.;
        //Look into whether this really has to be computed each iteration?
        for (unsigned i = 0; i < nRows; ++i)
        {
            //check assembly (before & after op) to see if actually implemented differently
            pmax = (currentColumn[i] < pmax) ? pmax : currentColumn[i];
            /*if (pixo[j] > pmax)
                pmax=pixo[j];
                */
        }

        /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
          AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
        for (int w = cte->cte_traps-1; w >= 0; --w){
            if (cte->qlevq_data[w] <= pmax)
            {
                ftrap = 0.0e0;
                ttrap = cte->cte_len; /*for referencing the image at 0*/
                fcarry = 0.0e0;

                /*GO UP THE COLUMN PIXEL BY PIXEL*/
                for(unsigned i = 0; i < nRows; ++i)
                {
                    pix_1 = currentColumn[i];

                    if ( (ttrap < cte->cte_len) || ( pix_1 >= cte->qlevq_data[w] - 1. ) )
                    {
                        if (currentColumn[i] >= 0 )
                        {
                            pix_1 = currentColumn[i] + fcarry; /*shuffle charge in*/
                            fcarry = pix_1 - floor(pix_1); /*carry the charge remainder*/
                            pix_1 = floor(pix_1); /*reset pixel*/
                        }

                        /*HAPPENS AFTER FIRST PASS*/
                        /*SHUFFLE CHARGE IN*/
                        //move out of loop to separate instance
                        if (i > 0)
                        {
                            if (pixf[i] < pixf[i-1])
                                ftrap *= (pixf[i] /  pixf[i-1]);
                        }

                        /*RELEASE THE CHARGE*/
                        padd_2=0.0;
                        if (ttrap <cte->cte_len){
                            ttrap += 1;
                            padd_2 = Pix(rprof->data,w,ttrap-1) *ftrap;
                        }

                        padd_3 = 0.0;
                        prem_3 = 0.0;
                        if ( pix_1 >= cte->qlevq_data[w]){
                            prem_3 =  cte->dpdew_data[w] / cte->n_par * pixf[i];  /*dpdew is 1 in file */
                            if (ttrap < cte->cte_len)
                                padd_3 = Pix(cprof->data,w,ttrap-1)*ftrap;
                            ttrap=0;
                            ftrap=prem_3;
                        }

                        /* Is it more efficient to dampen here rather than external and have to recompute the delta?
                         * the dampening loop is over all rows, unnecessarily, rather than just those computed here.
                         * does it add conditional to this loop now? Don't think compiler will optimize out
                         * conditional but since const, the proc should. (could make inline?)
                         */
                        double correction = padd_2 + padd_3 - prem_3;
                        double correction2 = correction * correction;
                        if (dampen)
                            correction = correction2 / (correction2 + rnAmp2);

                        currentColumn[i] += correction;
                    } /*replaces trap continue*/
                }/*end if i>0*/
            }/* end if qlevq > pmax, replaces continue*/

        }/*end for w*/
    //}
    return(status);
}
