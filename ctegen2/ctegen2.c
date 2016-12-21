#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"
#include "pcte.h"
#include "acs.h"

/*** this routine does the inverse CTE blurring... it takes an observed
  image and generates the image that would be pushed through the readout
  algorithm to generate the observation

  CTE_FF is found using the observation date of the data
  FIX_ROCRs is cte->fix_rocr
  Ws is the number of TRAPS that are < 999999

  this is sub_wfc3uv_raz2rac_par in jays code

  floor rounds to negative infinity
  ceiling rounds to positive infinity
  truncate rounds up or down to zero
  round goes to the nearest integer

  fff is the input cte scaling array calculated over all pixels
  This is a big old time sink function
 ***/

int inverse_cte_blur(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * fff, CTEParams * cte, const int verbose, const double expstart, const unsigned nRows)
{
    extern int status;

    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
      APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
      INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
      CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
      OVER 3 YEARS AFTER INSTALLATION*/

    /*cte scaling based on observation date*/
    const double cte_ff=  (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);
    cte->scale_frac=cte_ff;   /*save to param structure for header update*/
    const double rnAmp2 = cte->rn_amp * cte->rn_amp;

    if(verbose){
        sprintf(MsgText,"CTE_FF (scaling fraction by date) = %g",cte_ff);
        trlmessage(MsgText);
    }

    FloatTwoDArray cteRprof;
    FloatTwoDArray cteCprof;
    initFloatData(&cteRprof);
    initFloatData(&cteCprof);
    allocFloatData(&cteRprof, cte->rprof->data.nx, cte->rprof->data.ny, False);
    allocFloatData(&cteCprof, cte->cprof->data.nx, cte->cprof->data.ny, False);
    copyFloatDataToColumnMajor(&cteRprof, &cte->rprof->data);
    copyFloatDataToColumnMajor(&cteCprof, &cte->cprof->data);

#pragma omp parallel shared(rsc, rsz, cte, cteRprof, cteCprof, fff)
{
    double * observed    = NULL;
    double * model       = NULL;
    double * tempModel   = NULL;
    double * pix_ctef    = NULL;
    double * observedAll = NULL;

    //get rid of asserts
    assert(observedAll = (double *) malloc(sizeof(*observed)*nRows));
    assert(observed = (double *) malloc(sizeof(*observed)*nRows));
    assert(model = (double *) malloc(sizeof(*model)*nRows));
    assert(tempModel = (double *) malloc(sizeof(*tempModel)*nRows));
    assert(pix_ctef = (double *) malloc(sizeof(*pix_ctef)*nRows));

    //Due to rsc->sci.data being used as a mask for subarrays, keep this as 'dynamic'. If we can ensure no empty columns then
    //remove if(!hasFlux) and change schedule to 'static'
#ifdef _OPENMP
    unsigned chunkSize = (nRows / omp_get_num_procs())*0.1; //Have each thread take 10% of its share of the queue at a time
#endif

    //Experiment with chuckSize (ideally static would be best)
    #pragma omp for schedule (dynamic, chunkSize)
    for (unsigned i = 0; i < nRows; ++i)
    {
        Bool hasFlux = False;
        for (unsigned j = 0; j < nRows; ++j)
        {
            if (Pix(rsz->dq.data,i,j))
            {
                hasFlux = True;
                break;
            }
        }

        if (!hasFlux)
        {
            for (unsigned j = 0; j < nRows; ++j){
                Pix(rsc->sci.data, i, j) = 0; //check to see if even needed
            }
            continue;
        }

        //rsz->dq.data is being used as a mask to differentiate the actual pixels in a sub-array, from the entire full frame.
        //Eventually this will not be needed was will only be passed subarray and not the subarray padded to the full image size!
        for (unsigned j = 0; j < nRows; ++j)
        {
            observedAll[j] = Pix(rsz->sci.data,i,j); //Only left in to match master implementation
            observed[j] = Pix(rsz->dq.data,i,j) ? observedAll[j] : 0;
            pix_ctef[j] =  cte_ff * Pix(fff->sci.data, i, j);
        }

        unsigned NREDO = 0;
        Bool REDO;
        do
        {
            REDO = False; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
            memcpy(model, observedAll, sizeof(*observedAll)*nRows);

            /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
              IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
            for (unsigned NITINV = 1; NITINV <= cte->n_forward - 1; ++NITINV)
            {
                memcpy(tempModel, model, sizeof(*model)*nRows);
                simulateColumnReadout(model, pix_ctef, cte, &cteRprof, &cteCprof, nRows, cte->n_par);
                //memcpy(model, tempModel, sizeof(*model)*nRows);

                //Now that the updated readout has been simulated, subtract this from the model
                //to reproduce the actual image, without the CTE trails.
                //Whilst doing so, DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE, THIS IS
                //AN ADDITIONAL AID IN MITIGATING THE IMPACT OF READNOISE
                for (unsigned j = 0; j < nRows; ++j)
                {
                    double delta = model[j] - observed[j];
                    double delta2 = delta * delta;

                    //DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE
                    delta *= delta2 / (delta2 + rnAmp2);

                    //Now subtract the simulated readout
                    model[j] = tempModel[j] - delta;
                }
            }

            //Do the last forward iteration but don't dampen... no idea why???
            memcpy(tempModel, model, sizeof(*model)*nRows);
            simulateColumnReadout(model, pix_ctef, cte, &cteRprof, &cteCprof, nRows, cte->n_par);
            //Now subtract the simulated readout
            for (unsigned j = 0; j < nRows; ++j)
                model[j] = tempModel[j] - (model[j] - observed[j]);

            REDO =  cte->fix_rocr ? correctCROverSubtraction(pix_ctef, model, observed, nRows,
                    cte->thresh) : False;

        } while (REDO && ++NREDO < 5); //If really wanting 5 re-runs then need NREDO++

        //Update source array
        for (unsigned j = 0; j < nRows; ++j)
        {
            Pix(rsc->sci.data, i, j) = Pix(rsz->dq.data,i,j) ? model[j] : 0;
        }
    } //end loop over columns


    if (observedAll)
    {
        free(observedAll);
        observedAll = NULL;
    }
    if (tempModel)
    {
        free(tempModel);
        tempModel = NULL;
    }
    if (observed)
    {
        free(observed);
        observed = NULL;
    }
    if (model)
    {
        free(model);
        model = NULL;
    }
    if (pix_ctef)
    {
        free(pix_ctef);
        pix_ctef = NULL;
    }
}// close scope for #pragma omp parallel
    freeFloatData(&cteRprof);
    freeFloatData(&cteCprof);

    return(status);
}

/*This is the workhorse subroutine; it simulates the readout
  of one column currentColumn[].

  JDIM == nRows
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

int simulatePixelReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows)
{
    extern int status;

    double chargeToAdd;
    double extraChargeToAdd;//rethink name?
    double chargeToRemove;
    double pixel;
    double releasedFlux;
    double trappedFlux;
    int nTransfersFromTrap;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE*/
    //Look into whether this really has to be computed each iteration?
    //Since this is simulating the readout and thus moving pixels down and out, pmax can only get smaller with
    //each pixel transfer, never greater.
    double maxPixel = 10;
    for (unsigned i = 0; i < nRows; ++i)
        maxPixel = pixelColumn[i] > maxPixel ? pixelColumn[i] : maxPixel; //check assembly (before & after op) to see if actually implemented differently

    //Find highest charge trap to not exceed i.e. map pmax to an index
    unsigned maxChargeTrapIndex = cte->cte_traps-1;
    for (int w = maxChargeTrapIndex; w >= 0; --w)//go up or down? (if swap, change below condition)
    {
        if (cte->qlevq_data[w] <= maxPixel)//is any of this even needed or can we just directly map?
        {
            maxChargeTrapIndex = w;
            break;
        }
    }

    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    for (int w = maxChargeTrapIndex; w >= 0; --w)
    {
        nTransfersFromTrap = cte->cte_len; //for referencing the image at 0
        trappedFlux = 0;
        releasedFlux = 0;

        /*GO UP THE COLUMN PIXEL BY PIXEL*/
        for(unsigned i = 0; i < nRows; ++i)
        {
            //if (!pixelColumn[i])
              //  continue;

            pixel = pixelColumn[i];
            Bool isInsideTrailLength = nTransfersFromTrap < cte->cte_len;
            Bool isAboveChargeThreshold = pixel >= cte->qlevq_data[w] - 1.;
            if (isInsideTrailLength || isAboveChargeThreshold)
            {
                if (pixelColumn[i] >= 0 )//seems a shame to need check this every iteration
                {
                    pixel = pixelColumn[i] + releasedFlux; /*shuffle charge in*/
                    releasedFlux = pixel - floor(pixel); /*carry the charge remainder*/
                    pixel = floor(pixel); /*reset pixel*/
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
                chargeToAdd=0;
                if (isInsideTrailLength)
                {
                    ++nTransfersFromTrap;
                    chargeToAdd = rprof->data[w*rprof->ny + nTransfersFromTrap-1] * trappedFlux;
                }

                extraChargeToAdd = 0;
                chargeToRemove = 0;
                if (pixel >= cte->qlevq_data[w])
                {
                    chargeToRemove =  cte->dpdew_data[w] / cte->n_par * pixf[i];  /*dpdew is 1 in file */
                    if (nTransfersFromTrap < cte->cte_len)
                        extraChargeToAdd = cprof->data[w*cprof->ny + nTransfersFromTrap-1] * trappedFlux; //ttrap-1 may not be the same index as ref'd in rprof???
                    nTransfersFromTrap = 0;
                    trappedFlux = chargeToRemove;
                }

                pixelColumn[i] += chargeToAdd + extraChargeToAdd - chargeToRemove;
            } //replaces trap continue
        } //end for i
    } //end for w
    return(status);
}

int simulateColumnReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts)
{
    extern int status;

    //Take each pixel down the detector
    for (unsigned shift = 1; shift <= nPixelShifts; ++shift)
        simulatePixelReadout(pixelColumn, pixf, cte, rprof, cprof, nRows);

    return status;
}

Bool correctCROverSubtraction(double * const pix_ctef, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold)
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

    if (!pix_model || !pix_observed || !pix_ctef)
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
            double deltaFromOffendingCR = pix_model[jmax] - pix_observed[jmax];
            for (unsigned jj = j-10; jj <= j; ++jj)
            {
                if (pix_model[jj] - pix_observed[jj] > deltaFromOffendingCR)
                {
                    jmax = jj;
                    deltaFromOffendingCR = pix_model[jmax] - pix_observed[jmax];
                }
            }
            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
            for (unsigned jj = jmax; jj <= j; ++jj)
                pix_ctef[jj] *= 0.75;
        }
    } /*end for  j*/
    return redo;
}
