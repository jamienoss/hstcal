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



int inverseCTEBlurWithRowMajorInput(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart)
{
    //Convert all arrays to column major for efficiency.
    SingleGroup rszColumnMajor;
    SingleGroup rscColumnMajor;
    SingleGroup trapPixelMapColumnMajor;

    copySingleGroup(&rszColumnMajor, rsz, ROWMAJOR);
    copySingleGroup(&rscColumnMajor, rsc, ROWMAJOR);
    copySingleGroup(&trapPixelMapColumnMajor, trapPixelMap, ROWMAJOR);

    return inverseCTEBlurWithColumnMajorInput(rszColumnMajor, rscColumnMajor, trapPixelMapColumnMajor, cte, verbose, expstart);
}

int inverseCTEBlurWithColumnMajorInput(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart)
{
    return inverseCTEBlur(rsz, rsc, trapPixelMap, cte, verbose, expstart);
}

int inverseCTEBlur(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart)
{
    extern int status;

    const unsigned nRows = rsc->sci.data.ny;
    const unsigned nColumns = rsc->sci.data.nx;

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

#ifdef _OPENMP
    const unsigned nThreads = omp_get_num_procs();
#pragma omp parallel num_threads(nThreads) shared(rsc, rsz, cte, cteRprof, cteCprof, trapPixelMap)
#endif
{
    //get rid of asserts
    double * observed    = malloc(sizeof(*observed)*nRows);
    assert(observed);
    double * model       = malloc(sizeof(*model)*nRows);
    assert(model);
    double * tempModel   = malloc(sizeof(*tempModel)*nRows);
    assert(tempModel);
    double * pix_ctef    = malloc(sizeof(*pix_ctef)*nRows);
    assert(pix_ctef);
    double * observedAll = malloc(sizeof(*observedAll)*nRows);
    assert(observedAll);

    //Due to rsc->sci.data being used as a mask for subarrays, keep this as 'dynamic'. If we can ensure no empty columns then
    //remove if(!hasFlux) and change schedule to 'static'
#ifdef _OPENMP
    unsigned chunkSize = (nColumns / nThreads)*0.1; //Have each thread take 10% of its share of the queue at a time
#endif

    //Experiment with chuckSize (ideally static would be best)
#ifdef _OPENMP
    #pragma omp for schedule (dynamic, chunkSize)
#endif
    for (unsigned j = 0; j < nColumns; ++j)
    {
        Bool hasFlux = False;
        for (unsigned i = 0; i < nRows; ++i)
        {
            if (Pix(rsz->dq.data,j,i))
            {
                hasFlux = True;
                break;
            }
        }

        if (!hasFlux)
        {
            for (unsigned i = 0; i < nRows; ++i){
                Pix(rsc->sci.data, j, i) = 0; //check to see if even needed
            }
            continue;
        }

        //rsz->dq.data is being used as a mask to differentiate the actual pixels in a sub-array, from the entire full frame.
        //Eventually this will not be needed was will only be passed subarray and not the subarray padded to the full image size!
        for (unsigned i = 0; i < nRows; ++i)
        {
            observedAll[i] = Pix(rsz->sci.data,j,i); //Only left in to match master implementation
            observed[i] = Pix(rsz->dq.data,j,i) ? observedAll[i] : 0;
            pix_ctef[i] =  cte_ff * Pix(trapPixelMap->sci.data, j, i);
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
                for (unsigned i = 0; i < nRows; ++i)
                {
                    double delta = model[i] - observed[i];
                    double delta2 = delta * delta;

                    //DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE
                    delta *= delta2 / (delta2 + rnAmp2);

                    //Now subtract the simulated readout
                    model[i] = tempModel[i] - delta;
                }
            }

            //Do the last forward iteration but don't dampen... no idea why???
            memcpy(tempModel, model, sizeof(*model)*nRows);
            simulateColumnReadout(model, pix_ctef, cte, &cteRprof, &cteCprof, nRows, cte->n_par);
            //Now subtract the simulated readout
            for (unsigned i = 0; i < nRows; ++i)
                model[i] = tempModel[i] - (model[i] - observed[i]);

            REDO =  cte->fix_rocr ? correctCROverSubtraction(pix_ctef, model, observed, nRows,
                    cte->thresh) : False;

        } while (REDO && ++NREDO < 5); //If really wanting 5 re-runs then need NREDO++

        //Update source array
        for (unsigned i = 0; i < nRows; ++i)
        {
            Pix(rsc->sci.data, j, i) = Pix(rsz->dq.data, j, i) ? model[i] : 0;
        }
    } //end loop over columns


    delete((void*)&observedAll);
    delete((void*)&tempModel);
    delete((void*)&observed);
    delete((void*)&model);
    delete((void*)&pix_ctef);
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
                    double floored = floor(pixel);
                    releasedFlux = pixel - floored; /*carry the charge remainder*/
                    pixel = floored; /*reset pixel*/
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
    return status;
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
    for (unsigned i = 10; i < nRows-2; ++i)
    {
        if ( (( threshHold > pix_model[i] ) &&
                    ( threshHold > (pix_model[i] - pix_observed[i]))) ||

                (((pix_model[i] + pix_model[i+1]) < -12.) &&
                 (pix_model[i] + pix_model[i+1] - pix_observed[i] - pix_observed[i+1] < -12.)) ||

                (((pix_model[i] + pix_model[i+1] + pix_model[i+2]) < -15.) &&
                 ((pix_model[i] + pix_model[i+1] + pix_model[i+2] -pix_observed[i] -
                   pix_observed[i+1] - pix_observed[i+2]) <-15.))  )
        {
            redo = True;
            unsigned iMax = i;

            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
            double deltaFromOffendingCR = pix_model[iMax] - pix_observed[iMax];
            for (unsigned ii = i-10; ii <= i; ++ii)
            {
                if (pix_model[ii] - pix_observed[ii] > deltaFromOffendingCR)
                {
                    iMax = ii;
                    deltaFromOffendingCR = pix_model[iMax] - pix_observed[iMax];
                }
            }
            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
            for (unsigned ii = iMax; ii <= i; ++ii)
                pix_ctef[ii] *= 0.75;
        }
    } /*end for  j*/
    return redo;
}



int populateTrapPixelMap(SingleGroup * trapPixelMap, CTEParams * cte)
{
    extern int status;

    const unsigned nRows = trapPixelMap->sci.data.ny;
    const unsigned nColumns = trapPixelMap->sci.data.nx;

    int j;
    double cte_i=0.0;
    double cte_j=0.0;
    double ro=0;
    int io=0;
    double ff_by_col[nColumns][4];
    //float hardset=0.0;

    /*These are already in the parameter structure
      int     Ws              the number of traps < 999999, taken from pctetab read
      int     q_w[TRAPS];     the run of charge with level  cte->qlevq_data[]
      float   dpde_w[TRAPS];  the run of charge loss with level cte->dpdew_data[]

      float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel, TRAPS=999
      float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)

      The rprof array gives the fraction of charge that comes out of every parallel serial-shift
      the cummulative distribution in cprof then tells you what's left

*/
    /*SCALE BY 1 UNLESS THE PCTETAB SAYS OTHERWISE, I IS THE PACKET NUM
      THIS IS A SAFETY LOOP INCASE NOT ALL THE COLUMNS ARE POPULATED
      IN THE REFERENCE FILE*/

    for (unsigned i = 0; i < nColumns; ++i){
        ff_by_col[i][0]=1.;
        ff_by_col[i][1]=1.;
        ff_by_col[i][2]=1.;
        ff_by_col[i][3]=1.;
        j= cte->iz_data[i]; /*which column to scale*/
        ff_by_col[j][0]=cte->scale512[i];
        ff_by_col[j][1]=cte->scale1024[i];
        ff_by_col[j][2]=cte->scale1536[i];
        ff_by_col[j][3]=cte->scale2048[i];

        /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
          Index is figured on the final size of the image
          not the current size. Moved above
          */

        for(j=0; j<nRows; j++){//non contig swap with iter over i
            //Pix(trapPixelMap->sci.data,i,j)=hardset; //remove
            ro = j/512.0; /*ro can be zero, it's an index*/
            if (ro <0 ) ro=0.;
            if (ro > 2.999) ro=2.999; /*only 4 quads, 0 to 3*/
            io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
            cte_j= (j+1) / 2048.0;
            cte_i= ff_by_col[i][io] + (ff_by_col[i][io+1] -ff_by_col[i][io]) * (ro-io);
            Pix(trapPixelMap->sci.data,i,j) =  (cte_i*cte_j);
        }
    }

    /*FOR REFERENCE TO JAY ANDERSON'S CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

      int   iz_data[RAZ_ROWS];  column number in raz format
      double scale512[RAZ_ROWS];      scaling appropriate at row 512
      double scale1024[RAZ_ROWS];     scaling appropriate at row 1024
      double scale1536[RAZ_ROWS];     scaling appropriate at row 1536
      double scale2048[RAZ_ROWS];     scaling appropriate at row 2048
      */

    return(status);
}


