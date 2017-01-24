#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"
//#include "pcte.h"
//#include "acs.h"

char MsgText[256];
extern void trlmessage (char *message);


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

//#define PixColumnMajor(a,j,i) ( (a).storageOrder == ROWMAJOR ? (a).data[(j)*(a).tot_nx + (i)] : (a).data[(i)*(a).tot_ny + (j)] )
//#define PixColumnMajor(a,i,j) (a).data[(j)*(a).tot_nx + (i)] // := Pix

int inverseCTEBlurWithRowMajorInput(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart)
{
    clock_t t1 = clock();

    const unsigned nx = rsz->sci.data.nx;
    const unsigned ny = rsz->sci.data.ny;

    //Convert all arrays to column major for efficiency.
    SingleGroup rszColumnMajor;
    SingleGroup rscColumnMajor;
    SingleGroup trapPixelMapColumnMajor;

    //Initialize
    initSingleGroup(&rszColumnMajor);
    initSingleGroup(&rscColumnMajor);
    initSingleGroup(&trapPixelMapColumnMajor);

    //Allocate
    allocSingleGroup(&rszColumnMajor, nx, ny, False);
    allocSingleGroup(&rscColumnMajor, nx, ny, False);
    allocSingleGroup(&trapPixelMapColumnMajor, nx, ny, False);

    //Copy
    assert(!copySingleGroup(&rszColumnMajor, rsz, COLUMNMAJOR));
    assert(!copySingleGroup(&rscColumnMajor, rsc, COLUMNMAJOR));
    assert(!copySingleGroup(&trapPixelMapColumnMajor, trapPixelMap, COLUMNMAJOR));


/*  assert(!copySingleGroup(rsz, &rszColumnMajor, ROWMAJOR));
    assert(!copySingleGroup(rsc, &rscColumnMajor, ROWMAJOR));
    assert(!copySingleGroup(trapPixelMap, &trapPixelMapColumnMajor, ROWMAJOR));
    int ret = inverseCTEBlur(rsz, rsc, trapPixelMap, cte, verbose, expstart);
*/

    int ret = inverseCTEBlur(&rszColumnMajor, &rscColumnMajor, &trapPixelMapColumnMajor, cte, verbose, expstart);

    //copy data back
    copySingleGroup(rsc, &rscColumnMajor, ROWMAJOR);

    freeSingleGroup(&rszColumnMajor);
    freeSingleGroup(&rscColumnMajor);
    freeSingleGroup(&trapPixelMapColumnMajor);

    printf("\nTime taken to swap storage order: %f (secs)\n", ((float)(clock() - t1))/CLOCKS_PER_SEC);
    return ret;
}

int inverseCTEBlur(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart)
{
	//WARNING: This function assumes column major storage for 'rsz', 'rsc', & 'trapPixelMap'
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
    swapFloatStorageOrder(&cteRprof, &cte->rprof->data, COLUMNMAJOR);
    swapFloatStorageOrder(&cteCprof, &cte->cprof->data, COLUMNMAJOR);

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
            if (PixColumnMajor(rsz->dq.data,i,j))
            {
                hasFlux = True;
                break;
            }
        }
        if (!hasFlux)
        	continue;

        //allocated with calloc
       /* if (!hasFlux)
        {
            for (unsigned i = 0; i < nRows; ++i){
                PixColumnMajor(rsc->sci.data, i, j) = 0; //check to see if even needed
            }
            continue;
        }
        */

        //rsz->dq.data is being used as a mask to differentiate the actual pixels in a sub-array, from the entire full frame.
        //Eventually this will not be needed was will only be passed subarray and not the subarray padded to the full image size!
        for (unsigned i = 0; i < nRows; ++i)
        {
            observedAll[i] = PixColumnMajor(rsz->sci.data,i,j); //Only left in to match master implementation
            observed[i] = PixColumnMajor(rsz->dq.data,i,j) ? observedAll[i] : 0;
            pix_ctef[i] =  cte_ff * PixColumnMajor(trapPixelMap->sci.data, i, j);
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
            PixColumnMajor(rsc->sci.data, i, j) = PixColumnMajor(rsz->dq.data, i, j) ? model[i] : 0;
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
        unsigned column = cte->iz_data[i]; /*which column to scale*/
        ff_by_col[column][0]=cte->scale512[i];
        ff_by_col[column][1]=cte->scale1024[i];
        ff_by_col[column][2]=cte->scale1536[i];
        ff_by_col[column][3]=cte->scale2048[i];

        /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
          Index is figured on the final size of the image
          not the current size. Moved above
          */

        for (unsigned j = 0; j < nRows; ++j)
        {
            //Pix(trapPixelMap->sci.data,i,j)=hardset;
            ro = (double)j / 512.0; /*ro can be zero, it's an index*/
            if (ro > 2.999)
            	ro = 2.999; // only 4 quads, 0 to 3
            else if (ro < 0)
            	ro = 0;
            io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
            cte_j = (j+1) / 2048.0;
            cte_i = ff_by_col[i][io] + (ff_by_col[i][io+1] - ff_by_col[i][io]) * (ro-io);
            PixColumnMajor(trapPixelMap->sci.data,j,i) = (cte_i * cte_j);
        }
    }

    /*FOR REFERENCE TO JAY ANDERSON'S CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

      int   iz_data[nRows];  column number in raz format
      double scale512[nRows];      scaling appropriate at row 512
      double scale1024[nRows];     scaling appropriate at row 1024
      double scale1536[nRows];     scaling appropriate at row 1536
      double scale2048[nRows];     scaling appropriate at row 2048
      */

    return(status);
}

int cteSmoothImage(const SingleGroup * input, SingleGroup * output, double readNoiseAmp, unsigned maxThreads, int verbose)
{
    /*
       This routine will read in a RAZ image and will output the smoothest
       image that is consistent with being the observed image plus readnoise. (RSZ image)
       This is necessary because we want the CTE-correction algorithm to produce the smoothest
       possible reconstruction, consistent with the original image and the
       known readnoise.  This algorithm constructs a model that is smooth
       where the pixel-to-pixel variations can be thought of as being related
       to readnoise, but if the variations are too large, then it respects
       the pixel values.  Basically... it uses a 2-sigma threshold.

       This is strategy #1 in a two-pronged strategy to mitigate the readnoise
       amplification.  Strategy #2 will be to not iterate when the deblurring
       is less than the readnoise.

*/

    extern int status;

    const unsigned nRows = input->sci.data.ny;
    const unsigned nColumns = input->sci.data.nx;


    double rms=0;
    double nrms=0;

    clock_t begin = clock();



    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    memcpy(output->sci.data.data, input->sci.data.data, nRows*nColumns*sizeof(*input->sci.data.data));
    memcpy(output->dq.data.data, input->dq.data.data, nRows*nColumns*sizeof(*input->dq.data.data));
    return status;

    /*THE RSZ IMAGE JUST GETS UPDATED AS THE RAZ IMAGE IN THIS CASE*/
    if (readNoiseAmp < 0.1){
        trlmessage("rnsig < 0.1, No read-noise mitigation needed");
        return(status);
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
      */

    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, nColumns, nRows, True);

    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, nColumns, nRows, True);

#ifdef _OPENMP
    const unsigned nThreads = omp_get_num_procs(); //need to change this (and all others) to use maxThreads passed into main()
#pragma omp parallel num_threads(nThreads) shared(input, output, readNoiseAmp, rms, nrms, zadj, rnz)
#endif
    {
    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING nColumns*/
    FloatTwoDArray obs_loc;
    FloatTwoDArray rsz_loc;
    initFloatData(&obs_loc);
    initFloatData(&rsz_loc);
    allocFloatData(&obs_loc, nRows, 3, False);
    allocFloatData(&rsz_loc, nRows, 3, False);

    for(unsigned NIT=1; NIT<=100; ++NIT)
    {
        double rmsLocal = 0;
        double nrmsLocal = 0;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned i = 0; i < nColumns; ++i)
        {
            unsigned imid=i;
            /*RESET TO MIDDLE nColumns AT ENDPOINTS*/
            if (imid < 1)
                imid=1;
            if (imid == nColumns-1)
                imid = nColumns-2;

            //is copy needed?
            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            for(unsigned j = 0; j < nRows; ++j)
            {
                Pix(obs_loc, j, 0) = PixColumnMajor(input->sci.data, j, imid-1);
                Pix(obs_loc, j, 1) = PixColumnMajor(input->sci.data, j, imid);
                Pix(obs_loc, j, 2) = PixColumnMajor(input->sci.data, j, imid+1);

                Pix(rsz_loc, j, 0) = PixColumnMajor(output->sci.data,j,imid-1);
                Pix(rsz_loc, j, 1) = PixColumnMajor(output->sci.data,j,imid);
                Pix(rsz_loc, j, 2) = PixColumnMajor(output->sci.data,j,imid+1);
            }
            for (unsigned j = 0; j < nRows; ++j)
            {
             if(PixColumnMajor(input->dq.data, j, imid))
                PixColumnMajor(zadj.sci.data,j,i) = find_dadj(j, 1+i-imid, &obs_loc, &rsz_loc, nRows, readNoiseAmp);
            }
        } /*end the parallel for*/ //implicit omp barrier



        /*NOW GO OVER ALL THE nColumns AND nRows AGAIN TO SCALE THE PIXELS
        */
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned i = 0; i < nColumns; ++i){
            for(unsigned j = 0; j < nRows; ++j){
                if (PixColumnMajor(input->dq.data,j,i)){
                    PixColumnMajor(output->sci.data,j,i) += (PixColumnMajor(zadj.sci.data,j, i)*0.75);
                    PixColumnMajor(rnz.sci.data,j,i) = (PixColumnMajor(input->sci.data,j,i) - PixColumnMajor(output->sci.data,j,i));
                }
            }
        }//implicit omp barrier

        //we don't care if each thread does this to these shared variables
        rms=0;
        nrms=0;

        /*This is probably a time sink because the arrays are being
          accessed out of storage order, careful of page faults */
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned j = 0; j < nRows; ++j)
        {
            nrmsLocal=0;
            rmsLocal=0;
            for(unsigned i = 0; i < nColumns; ++i)
            {
                if ( (fabs(PixColumnMajor(input->sci.data,j, i)) > 0.1) ||
                        (fabs(PixColumnMajor(output->sci.data,j,i)) > 0.1))
                {
                    unsigned tmp = PixColumnMajor(rnz.sci.data,j, i);
                    rmsLocal  +=  tmp*tmp;
                    nrmsLocal += 1.0;
                }
            }
#ifdef _OPENMP
            #pragma omp critical (accumulate)
#endif
            {
                rms  += rmsLocal;
                nrms += nrmsLocal;
            }
        }//implicit omp barrier

        //we don't care if each thread does this to these shared variables (above crit sec & barrier make this one ok)
        rms = sqrt(rms/nrms);

        // if it is true that one breaks then it is will be true for all
        /*epsilon type comparison*/
        if ( (readNoiseAmp-rms) < 0.00001) break; /*this exits the NIT for loop*/
#ifdef _OPENMP
    #pragma omp barrier
#endif
    } /*end NIT*/
    freeFloatData(&obs_loc);
    freeFloatData(&rsz_loc);
    } // close parallel block
    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);

    if (verbose)
    {
    	double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
    	sprintf(MsgText,"Time taken to smooth image: %.2f(s) with %i procs/threads\n",timeSpent/maxThreads,maxThreads);
    	trlmessage(MsgText);
    }

    return (status);
}

double find_dadj(const unsigned i, const unsigned j, const FloatTwoDArray * obsloc, const FloatTwoDArray * rszloc, const unsigned nRows, const double readNoiseAmp)
{
    /*
       This function determines for a given pixel how it can
       adjust in a way that is not inconsistent with its being
       readnoise.  To do this, it looks at its upper and lower
       neighbors and sees whether it is consistent with either
       (modulo readnoise).  To the extent that it is consistent
       then move it towards them.  But also bear in mind that
       that we don't want it to be more than 2 RN sigmas away
       from its original value.  This is pretty much a tug of
       war... with readnoise considerations pushing pixels to
       be closer to their neighbors, but the original pixel
       values also pull to keep the pixel where it was.  Some
       accommodation is made for both considerations.
       */

    const double mval = Pix(*rszloc, i, j);
    const double dval0  = Pix(*obsloc, i, j) - mval;
    double dval0u = dval0;

    if (dval0u >1.0)
        dval0u =  1.0;
    if (dval0u <-1.0)
        dval0u = -1.0;

    /*COMPARE THE SURROUNDING PIXELS*/
    double dval9 = 0;
    if (j == 1 &&  nRows-1>=i  && i>0 )
    {
        dval9 = Pix(*obsloc, i, j-1)  - Pix(*rszloc, i, j-1) +
                Pix(*obsloc, i, j)    - Pix(*rszloc, i, j)  +
                Pix(*obsloc, i, j+1)  - Pix(*rszloc, i, j+1) +
                Pix(*obsloc, i-1, j-1)- Pix(*rszloc, i-1, j-1) +
                Pix(*obsloc, i-1, j)  - Pix(*rszloc, i-1, j) +
                Pix(*obsloc, i-1, j+1)- Pix(*rszloc, i-1, j+1) +
                Pix(*obsloc, i+1, j-1)- Pix(*rszloc, i+1, j-1) +
                Pix(*obsloc, i+1, j)  - Pix(*rszloc, i+1, j) +
                Pix(*obsloc, i+1, j+1)- Pix(*rszloc, i+1, j+1);
    }

    dval9 =dval9 / 9.;
    double dval9u = dval9;

    if (dval9u > (readNoiseAmp*0.33))
        dval9u =  readNoiseAmp*0.33;
    else if (dval9u <  readNoiseAmp*-0.33)
        dval9u = readNoiseAmp*-0.33;

    double dmod1 = 0;
    if (i>0)
        dmod1 = Pix(*rszloc, i, j-1) - mval;

    double dmod1u = dmod1;
    if (dmod1u > readNoiseAmp*0.33)
        dmod1u =  readNoiseAmp*0.33;
    else if (dmod1u < readNoiseAmp*-0.33)
        dmod1u = readNoiseAmp*-0.33;

    double dmod2 = 0;
    if (i < nRows-1)
        dmod2 = Pix(*rszloc, i, j+1) - mval;

    double dmod2u = dmod2;
    if (dmod2u > readNoiseAmp*0.33)
        dmod2u =  readNoiseAmp*0.33;
    else if (dmod2u < readNoiseAmp*-0.33)
        dmod2u = readNoiseAmp*-0.33;


    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    const double readNoiseAmp2 = readNoiseAmp*readNoiseAmp;
    const double w0 =   (dval0*dval0) / ((dval0*dval0)+ 4.0*(readNoiseAmp2));
    const double w9 =   (dval9*dval9) / ((dval9*dval9)+ 18.0*(readNoiseAmp2));
    const double w1 = (4*readNoiseAmp2) / ((dmod1*dmod1)+4.0*(readNoiseAmp2));
    const double w2 = (4*readNoiseAmp2) / ((dmod2*dmod2)+4.0*(readNoiseAmp2));

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    return  (dval0u * w0 * 0.25f) + /* desire to keep the original pixel value */
            (dval9u * w9 * 0.25f) + /* desire to keep the original sum over 3x3*/
            (dmod1u * w1 * 0.25f) + /* desire to get closer to the pixel below*/
            (dmod2u * w2 * 0.25f) ; /* desire to get closer to the pixel above*/
}

# define RAZ_COLS 8412
# define RAZ_ROWS 2070
double find_dadj_old(int i ,int j, double obsloc[][RAZ_ROWS], double rszloc[][RAZ_ROWS], double rnsig){
    /*
       This function determines for a given pixel how it can
       adjust in a way that is not inconsistent with its being
       readnoise.  To do this, it looks at its upper and lower
       neighbors and sees whether it is consistent with either
       (modulo readnoise).  To the extent that it is consistent
       then move it towards them.  But also bear in mind that
       that we don't want it to be more than 2 RN sigmas away
       from its original value.  This is pretty much a tug of
       war... with readnoise considerations pushing pixels to
       be closer to their neighbors, but the original pixel
       values also pull to keep the pixel where it was.  Some
       accommodation is made for both considerations.
       */

    extern int status;

    double mval=0.0;
    double    dval0, dval0u, w0;
    double    dval9, dval9u, w9;
    double    dmod1, dmod1u, w1;
    double    dmod2, dmod2u, w2;

    dval0=0.;
    dval0u=0.;
    w0=0.;
    dval9=0.;
    dval9u=0.;
    w9=0.;
    dmod1=0.;
    dmod1u=0.;
    w1=0.;
    dmod2=0.;
    dmod2u=0.;
    w2=0.;

    mval = rszloc[i][j];
    dval0  = obsloc[i][j] - mval;
    dval0u = dval0;

    if (dval0u >1.0)
        dval0u =  1.0;
    if (dval0u <-1.0)
        dval0u = -1.0;

    dval9 = 0.;

    /*COMPARE THE SURROUNDING PIXELS*/
    if (i==1 &&  RAZ_ROWS-1>=j  && j>0 ) {

        dval9 = obsloc[i][j-1]  - rszloc[i][j-1] +
            obsloc[i][j]    - rszloc[i][j]  +
            obsloc[i][j+1]  - rszloc[i][j+1] +
            obsloc[i-1][j-1]- rszloc[i-1][j-1] +
            obsloc[i-1][j]  - rszloc[i-1][j] +
            obsloc[i-1][j+1]- rszloc[i-1][j+1] +
            obsloc[i+1][j-1]- rszloc[i+1][j-1] +
            obsloc[i+1][j]  - rszloc[i+1][j] +
            obsloc[i+1][j+1]- rszloc[i+1][j+1];
    }

    dval9 =dval9 / 9.;
    dval9u = dval9;

    if (dval9u > (rnsig*0.33))
        dval9u =  rnsig*0.33;
    if (dval9u <  rnsig*-0.33)
        dval9u = rnsig*-0.33;

    dmod1 = 0.;
    if (j>0)
        dmod1 = rszloc[i][j-1] - mval;

    dmod1u = dmod1;
    if (dmod1u > rnsig*0.33)
        dmod1u =  rnsig*0.33;
    if (dmod1u < rnsig*-0.33)
        dmod1u = rnsig*-0.33;

    dmod2 = 0.;
    if (j < RAZ_ROWS-1)
        dmod2 =  rszloc[i][j+1] - mval;

    dmod2u = dmod2;
    if (dmod2u > rnsig*0.33)
        dmod2u =  rnsig*0.33;
    if (dmod2u < rnsig*-0.33)
        dmod2u = rnsig*-0.33;


    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    w0 =   (dval0*dval0) / ((dval0*dval0)+ 4.0*(rnsig*rnsig));
    w9 =   (dval9*dval9) / ((dval9*dval9)+ 18.0*(rnsig*rnsig));
    w1 = (4*rnsig*rnsig) / ((dmod1*dmod1)+4.0*(rnsig*rnsig));
    w2 = (4*rnsig*rnsig) / ((dmod2*dmod2)+4.0*(rnsig*rnsig));

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    return ((dval0u * w0 * 0.25f) + /* desire to keep the original pixel value */
            (dval9u*w9*0.25f) + /* desire to keep the original sum over 3x3*/
            (dmod1u*w1*0.25f) + /*desire to get closer to the pixel below*/
            (dmod2u*w2*0.25f)) ; /*desire to get closer to the pixel above*/

    //return(status);
}
int cteSmoothImage_old(SingleGroup *raz, SingleGroup *rsz, double rnsig, int max_threads){
    /*
       This routine will read in a RAZ image and will output the smoothest
       image that is consistent with being the observed image plus readnoise. (RSZ image)
       This is necessary because we want the CTE-correction algorithm to produce the smoothest
       possible reconstruction, consistent with the original image and the
       known readnoise.  This algorithm constructs a model that is smooth
       where the pixel-to-pixel variations can be thought of as being related
       to readnoise, but if the variations are too large, then it respects
       the pixel values.  Basically... it uses a 2-sigma threshold.
       This is strategy #1 in a two-pronged strategy to mitigate the readnoise
       amplification.  Strategy #2 will be to not iterate when the deblurring
       is less than the readnoise.
*/

    extern int status;

    int i, j, NIT; /*loop variables*/
    int imid;
    double dptr=0.0;
    double  rms=0.0;
    double  rmsu=0.0;
    double nrms=0.0;
    double nrmsu=0.0;
    float hardset=0.0f;
    double setdbl=0.0;


    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING RAZ_COLS*/
    double obs_loc[3][RAZ_ROWS] ;
    double rsz_loc[3][RAZ_ROWS] ;

    NIT=1;

    /*ALL ELEMENTS TO FLAG*/
    for(i=0;i<3;i++){
        for (j=0; j<RAZ_ROWS; j++){
            obs_loc[i][j]=setdbl;
            rsz_loc[i][j]=setdbl;
        }
    }

    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, RAZ_COLS, RAZ_ROWS, True);


    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    for(i=0;i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS;j++){
            Pix(rsz->sci.data,i,j) = Pix(raz->sci.data,i,j);
            Pix(rsz->dq.data,i,j) = Pix(raz->dq.data,i,j);
            Pix(rnz.sci.data,i,j) = hardset;
            Pix(zadj.sci.data,i,j) = hardset;
        }
    }


    /*THE RSZ IMAGE JUST GETS UPDATED AS THE RAZ IMAGE IN THIS CASE*/
    if (rnsig < 0.1){
        trlmessage("rnsig < 0.1, No read-noise mitigation needed");
        return(status);
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
      */

    rms=setdbl;

    for(NIT=1; NIT<=100; NIT++){
        #pragma omp parallel for schedule(dynamic) \
        private(i,j,imid,obs_loc,rsz_loc,dptr)\
        shared(raz, rsz, rnsig,rms,nrms, zadj)
        for(i=0; i<RAZ_COLS; i++){
            imid=i;
            /*RESET TO MIDDLE RAZ_COLS AT ENDPOINTS*/
            if (imid < 1)
                imid=1;
            if (imid == RAZ_COLS-1)
                imid = RAZ_COLS-2;

            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            for(j=0; j<RAZ_ROWS; j++){
                obs_loc[0][j] = Pix(raz->sci.data,imid-1,j);
                obs_loc[1][j] = Pix(raz->sci.data,imid,j);
                obs_loc[2][j] = Pix(raz->sci.data,imid+1,j);

                rsz_loc[0][j] = Pix(rsz->sci.data,imid-1,j);
                rsz_loc[1][j] = Pix(rsz->sci.data,imid,j);
                rsz_loc[2][j] = Pix(rsz->sci.data,imid+1,j);
            }
            for (j=0; j<RAZ_ROWS; j++){
             if(Pix(raz->dq.data,imid,j)) {
                Pix(zadj.sci.data,i,j) = find_dadj_old(1+i-imid,j, obs_loc, rsz_loc, rnsig);
              }
            }
        } /*end the parallel for*/

        /*NOW GO OVER ALL THE RAZ_COLS AND RAZ_ROWS AGAIN TO SCALE THE PIXELS
        */
        for(i=0; i<RAZ_COLS;i++){
            for(j=0; j<RAZ_ROWS; j++){
                if (Pix(raz->dq.data,i,j)){
                    Pix(rsz->sci.data,i,j) +=  (Pix(zadj.sci.data,i,j)*0.75);
                    Pix(rnz.sci.data,i,j) = (Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j));
                }
            }
        }

        rms=setdbl;
        nrms=setdbl;

        /*This is probably a time sink because the arrays are being
          accessed out of storage order, careful of page faults */
        #pragma omp parallel for schedule(dynamic,1)\
        private(i,j,rmsu,nrmsu) \
        shared(raz,rsz,rms,rnsig,nrms)
        for(j=0; j<RAZ_ROWS; j++){
            nrmsu=setdbl;
            rmsu=setdbl;
            for(i = 0;i<RAZ_COLS; i++){
                if ( (fabs(Pix(raz->sci.data,i,j)) > 0.1) ||
                        (fabs(Pix(rsz->sci.data,i,j)) > 0.1) ){
                    rmsu  +=  ( Pix(rnz.sci.data,i,j) * Pix(rnz.sci.data,i,j) );
                    nrmsu += 1.0;
                }
            }
            #pragma omp critical (rms)
            {   rms  += rmsu;
                nrms += nrmsu;
            }
        }
        rms = sqrt(rms/nrms);

        /*epsilon type comparison*/
        if ( (rnsig-rms) < 0.00001) break; /*this exits the NIT for loop*/
    } /*end NIT*/

    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);

    return (status);
}
