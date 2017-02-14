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

//#define PixColumnMajor(a,j,i) ( (a).storageOrder == ROWMAJOR ? (a).data[(j)*(a).tot_nx + (i)] : (a).data[(i)*(a).tot_ny + (j)] )
//#define PixColumnMajor(a,i,j) (a).data[(j)*(a).tot_nx + (i)] // := Pix

void initPtrRegister(PtrRegister * reg)
{
    reg->cursor = 0; //points to last ptr NOT next slot
    reg->length = PTR_REGISTER_LENGTH;
    reg->ptrs = calloc(reg->length+1, sizeof(*reg->ptrs));
    assert(reg->ptrs);
    reg->freeFunctions = calloc(reg->length+1, sizeof(*reg->freeFunctions));
    if (!reg->freeFunctions)
    {
        free(reg->ptrs);
        assert(0);
    }
    reg->ptrs[0] = reg; //this ptr
    reg->freeFunctions[0] = &free;
}
void addPtr(PtrRegister * reg, void * ptr, void * freeFunc)
{
    if (!reg || !ptr || !freeFunc)
        return;

    if (++reg->cursor >= reg->length)
    {
        reg->length += 10;
        assert(realloc(&reg->ptrs, reg->length*sizeof(*reg->ptrs)));
        assert(realloc(reg->freeFunctions, reg->length*sizeof(*reg->freeFunctions)));
    }
    reg->ptrs[reg->cursor] = ptr;
    reg->freeFunctions[reg->cursor] = freeFunc;
}
void freePtr(PtrRegister * reg, void * ptr)
{
    if (!reg || !ptr)
        return;

    unsigned i;
    for (i = reg->cursor; i > 0 ; --i)
    {
        if (reg->ptrs[i] == ptr)
            break;
    }

    //call function to free ptr
    reg->freeFunctions[i](ptr);

    if (i == reg->cursor)
    {
        reg->ptrs[i] = NULL;
        reg->freeFunctions[i] = NULL;
    }
    else
    {
        //move last one into gap to close - not a stack so who cares
        reg->ptrs[i] = reg->ptrs[reg->cursor];
        reg->ptrs[reg->cursor] = NULL;
        reg->freeFunctions[i] = reg->freeFunctions[reg->cursor];
        reg->freeFunctions[reg->cursor] = NULL;
    }
    --reg->cursor;
}
void freeAll(PtrRegister * reg)
{
    if (!reg || reg->length == 0)
        return;

    for (unsigned i = 1; i < reg->cursor; ++i)
    {
        if (reg->freeFunctions[i] && reg->ptrs[i])
        {
            reg->freeFunctions[i](reg->ptrs[i]);
            reg->ptrs[i] = NULL;
            reg->freeFunctions[i] = NULL;
        }
    }
    reg->cursor = 0;
    reg->length = 0;
    // free 'itself'
    reg->freeFunctions[0](reg->ptrs);
    reg->ptrs[0] = NULL;
    reg->freeFunctions[0](reg->freeFunctions);
    reg->freeFunctions[0] = NULL;
}


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

    int ret = inverseCTEBlur(&rszColumnMajor, &rscColumnMajor, &trapPixelMapColumnMajor, cte);

    //copy data back
    copySingleGroup(rsc, &rscColumnMajor, ROWMAJOR);

    freeSingleGroup(&rszColumnMajor);
    freeSingleGroup(&rscColumnMajor);
    freeSingleGroup(&trapPixelMapColumnMajor);

    printf("\nTime taken to swap storage order: %f (secs)\n", ((float)(clock() - t1))/CLOCKS_PER_SEC);
    return ret;
}

int inverseCTEBlur(const SingleGroup * input, SingleGroup * output, SingleGroup * trapPixelMap, CTEParams * cte)
{
    //WARNING: This function assumes column major storage for 'rsz', 'rsc', & 'trapPixelMap'
    extern int status;

    const unsigned nRows = output->sci.data.ny;
    const unsigned nColumns = output->sci.data.nx;
    const double rnAmp2 = cte->rn_amp * cte->rn_amp;
    const FloatTwoDArray * cteRprof  = &cte->rprof->data;
    const FloatTwoDArray * cteCprof = &cte->cprof->data;

#ifdef _OPENMP
    unsigned nThreads = omp_get_num_procs();
#pragma omp parallel shared(input, output, cte, cteRprof, cteCprof, trapPixelMap)
#endif
{
    //get rid of asserts
    double * observed    = NULL;//malloc(sizeof(*observed)*nRows);
    //assert(observed);
    double * model       = malloc(sizeof(*model)*nRows);
    assert(model);
    double * tempModel   = malloc(sizeof(*tempModel)*nRows);
    assert(tempModel);
    double * observedAll = malloc(sizeof(*observedAll)*nRows);
    assert(observedAll);
    float * traps = NULL;

    //Due to input->sci.data being used as a mask for subarrays, keep this as 'dynamic'. If we can ensure no empty columns then
    //remove if(!hasFlux) and change schedule to 'static'
#ifdef _OPENMP
    unsigned chunkSize = 1;//(nColumns / nThreads)*0.1; //Have each thread take 10% of its share of the queue at a time
#endif

    //Experiment with chuckSize (ideally static would be best)
#ifdef _OPENMP
    #pragma omp for schedule (dynamic, chunkSize)
#endif
    for (unsigned j = 0; j < nColumns; ++j)
    {
        /*Bool hasFlux = False;
        for (unsigned i = 0; i < nRows; ++i)
        {
            if (PixColumnMajor(input->dq.data,i,j))
            {
                hasFlux = True;
                break;
            }
        }
        if (!hasFlux)
            continue;
*/
        //input->dq.data is being used as a mask to differentiate the actual pixels in a sub-array, from the entire full frame.
        //Eventually this will not be needed was will only be passed subarray and not the subarray padded to the full image size!
        //memcpy(observedAll, input->sci.data.data+j*nRows, nRows*sizeof(*input->sci.data.data));
        for (unsigned i = 0; i < nRows; ++i)
        {
            observedAll[i] = PixColumnMajor(input->sci.data,i,j); //Only left in to match master implementation
            //observed[i] = PixColumnMajor(input->dq.data,i,j) ? observedAll[i] : 0;
        }

        traps = &(PixColumnMajor(trapPixelMap->sci.data, 0, j));
        observed = observedAll;
        unsigned NREDO = 0;
        Bool REDO;
        do
        {
            REDO = False; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
            memcpy(model, observedAll, nRows*sizeof(*observedAll));

            /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
              IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
            for (unsigned NITINV = 1; NITINV <= cte->n_forward - 1; ++NITINV)
            {
                memcpy(tempModel, model, nRows*sizeof(*model));
                simulateColumnReadout(model, traps, cte, cteRprof, cteCprof, nRows, cte->n_par);

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
            simulateColumnReadout(model, traps, cte, cteRprof, cteCprof, nRows, cte->n_par);
            //Now subtract the simulated readout
            for (unsigned i = 0; i < nRows; ++i)
                model[i] = tempModel[i] - (model[i] - observed[i]);

            REDO = cte->fix_rocr ? correctCROverSubtraction(traps, model, observed, nRows,
                    cte->thresh) : False;

           // if (REDO)
             //   printf("Readout cosmic ray detected, redoing computation\n");
        } while (REDO && ++NREDO < 5); //If really wanting 5 re-runs then use NREDO++

        //Update source array
        //can't use memcpy as arrays of diff types
        for (unsigned i = 0; i < nRows; ++i)
            PixColumnMajor(output->sci.data, i, j) = model[i];//PixColumnMajor(input->dq.data, i, j) ? model[i] : 0;
    } //end loop over columns


    delete((void*)&observedAll);
    delete((void*)&tempModel);
    //delete((void*)&observed);
    delete((void*)&model);
}// close scope for #pragma omp parallel

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

int simulatePixelReadout(double * const pixelColumn, const float * const traps, const CTEParams * const cte,
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

            if (!isInsideTrailLength && !isAboveChargeThreshold)
                continue;

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
                if ((double)traps[i] < (double)traps[i-1])
                    trappedFlux *= ((double)traps[i] / (double)traps[i-1]);
            }

            /*RELEASE THE CHARGE*/
            chargeToAdd = 0;
            if (isInsideTrailLength)
            {
                ++nTransfersFromTrap;
                chargeToAdd = rprof->data[w*rprof->ny + nTransfersFromTrap-1] * trappedFlux;
            }

            extraChargeToAdd = 0;
            chargeToRemove = 0;
            if (pixel >= cte->qlevq_data[w])
            {
                chargeToRemove =  cte->dpdew_data[w] / cte->n_par * (double)traps[i];  /*dpdew is 1 in file */
                if (nTransfersFromTrap < cte->cte_len)
                    extraChargeToAdd = cprof->data[w*cprof->ny + nTransfersFromTrap-1] * trappedFlux; //ttrap-1 may not be the same index as ref'd in rprof???
                nTransfersFromTrap = 0;
                trappedFlux = chargeToRemove;
            }

            pixelColumn[i] += chargeToAdd + extraChargeToAdd - chargeToRemove;
        } //end for i
    } //end for w
    return status;
}

int simulateColumnReadout(double * const pixelColumn, const float * const traps, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts)
{
    extern int status;

    //Take each pixel down the detector
    for (unsigned shift = 1; shift <= nPixelShifts; ++shift)
        simulatePixelReadout(pixelColumn, traps, cte, rprof, cprof, nRows);

    return status;
}

Bool correctCROverSubtraction(float * const traps, const double * const pix_model, const double * const pix_observed,
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

    if (!pix_model || !pix_observed || !traps)
        return False;

    Bool redo = False;
    for (unsigned i = 10; i < nRows-2; ++i)
    {
        if ( (( threshHold > pix_model[i] ) &&
                    ( threshHold > (pix_model[i] - pix_observed[i]))) ||

                (((pix_model[i] + pix_model[i+1]) < -12.) &&
                 (pix_model[i] + pix_model[i+1] - pix_observed[i] - pix_observed[i+1] < -12.)) ||

                (((pix_model[i] + pix_model[i+1] + pix_model[i+2]) < -15.) &&
                 ((pix_model[i] + pix_model[i+1] + pix_model[i+2] - pix_observed[i] -
                   pix_observed[i+1] - pix_observed[i+2]) < -15.))  )
        {
            redo = True;
            unsigned iMax = i;

            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
            for (unsigned ii = i-10; ii <= i; ++ii)
            {
                if ( (pix_model[ii] - pix_observed[ii]) > (pix_model[iMax] - pix_observed[iMax]) )
                    iMax = ii;
            }
            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
            for (unsigned ii = iMax; ii <= i; ++ii)
                traps[ii] *= 0.75;
        }
    } /*end for  j*/
    return redo;
}

int populateTrapPixelMap(SingleGroup * trapPixelMap, CTEParams * cte, const int verbose, const double expstart)
{
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
    /*FOR REFERENCE TO JAY ANDERSON'S CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

          int   iz_data[nRows];  column number in raz format
          double scale512[nRows];      scaling appropriate at row 512
          double scale1024[nRows];     scaling appropriate at row 1024
          double scale1536[nRows];     scaling appropriate at row 1536
          double scale2048[nRows];     scaling appropriate at row 2048
     */
    extern int status;

    const unsigned nRows = trapPixelMap->sci.data.ny;
    const unsigned nColumns = trapPixelMap->sci.data.nx;

    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
          APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
          INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
          CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
          OVER 3 YEARS AFTER INSTALLATION*/
    // cte scaling based on observation date
    const double cteScale =  (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);
    cte->scale_frac = cteScale; // save to param structure for header update
    if (verbose)
    {
        sprintf(MsgText,"CTE_FF (scaling fraction by date) = %g",cteScale);
        trlmessage(MsgText);
    }
    //double ff_by_col[nColumns][4];

   /* for (unsigned i = 0; i < nColumns; ++i){
        ff_by_col[i][0]=1;
        ff_by_col[i][1]=1;
        ff_by_col[i][2]=1;
        ff_by_col[i][3]=1;
        unsigned column = cte->iz_data[i]; //which column to scale
        ff_by_col[column][0]=cte->scale512[i];
        ff_by_col[column][1]=cte->scale1024[i];
        ff_by_col[column][2]=cte->scale1536[i];
        ff_by_col[column][3]=cte->scale2048[i];
    }*/

#ifdef _OPENMP
    #pragma omp parallel shared(trapPixelMap, cte)
#endif
    {
    double trapColumnScale[4];
    double cte_i;
    double cte_j;
    double ro;
    int io;

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 1)
#endif
    for (unsigned i = 0; i < nColumns; ++i)
    {
        //this is only my guess of what should be happening
        unsigned column = cte->iz_data[i]; /*which column to scale*/
        //should check 'column' within bounds(?)
        trapColumnScale[0]=1;//cte->scale512[column];
        trapColumnScale[1]=1;//cte->scale1024[column];
        trapColumnScale[2]=1;//cte->scale1536[column];
        trapColumnScale[3]=1;//cte->scale2048[column];
        /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
          Index is figured on the final size of the image
          not the current size. Moved above
         */
        for (unsigned j = 0; j < nRows; ++j)
        {
            unsigned jj = cte->subarrayColumnOffset + j;
            ro = jj / 512.0; /*ro can be zero, it's an index*/
            if (ro > 2.999)
            	ro = 2.999; // only 4 quads, 0 to 3
            else if (ro < 0)
            	ro = 0;
            io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
            cte_j = (jj+1) / 2048.0;
            cte_i = trapColumnScale[io] + (trapColumnScale[io+1] - trapColumnScale[io]) * (ro - io);
            //cte_i = ff_by_col[i][io] + (ff_by_col[i][io+1] - ff_by_col[i][io]) * (ro-io);
            PixColumnMajor(trapPixelMap->sci.data,j,i) = cte_i * cte_j * cteScale;
        }
    }
    } // end parallel block
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
    copySingleGroup(output, input, input->sci.data.storageOrder);
    //memcpy(output->sci.data.data, input->sci.data.data, nRows*nColumns*sizeof(*input->sci.data.data));
    //memcpy(output->dq.data.data, input->dq.data.data, nRows*nColumns*sizeof(*input->dq.data.data));

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

    //To get rid of this and adjust in place i.e. using only output:
    //Don't use pointers to output for obs_loc & rsz_loc
    //Copy columns and then just shift these copies (boundary case might be annoying)
    //Use schedule(static) and pre (inner for loop) copy boundary columns to avoid race conditions
    SingleGroup adjustment;
    initSingleGroup(&adjustment);
    allocSingleGroup(&adjustment, nColumns, nRows, False);

    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, nColumns, nRows, False);

#ifdef _OPENMP
    #pragma omp parallel shared(input, output, readNoiseAmp, rms, nrms, rnz)
#endif
    {
    const float * obs_loc[3];
    const float * rsz_loc[3];

    double rmsLocal;
    double nrmsLocal;
    for(unsigned iter = 0; iter < 100; ++iter)
    {
        rmsLocal = 0;
        nrmsLocal = 0;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned i = 0; i < nColumns; ++i)
        {
            unsigned imid = i;
            /*RESET TO MIDDLE nColumns AT ENDPOINTS*/
            // This seems odd, the edge columns get accounted for twice?
            if (i == 0)
                imid = 1;
            else if (i == nColumns-1) // NOTE: use of elseif breaks if nColumns = 1
                imid = nColumns-2;

            /*LOCATE THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            obs_loc[0] = input->sci.data.data + (imid-1)*nRows;
            obs_loc[1] = obs_loc[0] + nRows;
            obs_loc[2] = obs_loc[1] + nRows;

            rsz_loc[0] = output->sci.data.data + (imid-1)*nRows;
            rsz_loc[1] = rsz_loc[0] + nRows;
            rsz_loc[2] = rsz_loc[1] + nRows;

            for (unsigned j = 0; j < nRows; ++j)
            {
                //if(PixColumnMajor(input->dq.data, j, imid))//I think imid should be i here?
                    PixColumnMajor(adjustment.sci.data,j,i) = find_dadj(1+i-imid, j, nRows, obs_loc, rsz_loc, readNoiseAmp);
            }
        } /*end the parallel for*/ //implicit omp barrier

        //move this entire section into the above and remove need for 'adjustment' group!!!
        //NOW GO OVER ALL THE nColumns AND nRows AGAIN TO SCALE THE PIXELS
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned i = 0; i < nColumns; ++i)
        {
            for(unsigned j = 0; j < nRows; ++j)
            {
                //if (PixColumnMajor(input->dq.data,j,i))
                //{
                    PixColumnMajor(output->sci.data,j,i) += (PixColumnMajor(adjustment.sci.data,j, i)*0.75);
                    PixColumnMajor(rnz.sci.data,j,i) = (PixColumnMajor(input->sci.data,j,i) - PixColumnMajor(output->sci.data,j,i));
                //}
            }
        }//implicit omp barrier

#ifdef _OPENMP
        #pragma omp single
        {
            rms=0;
            nrms=0;
        }
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for(unsigned j = 0; j < nColumns; ++j)
        {
            for(unsigned i = 0; i < nRows; ++i)
            {
                if ( /*PixColumnMajor(input->dq.data, i, j) &&*/
                     (fabs(PixColumnMajor(input->sci.data, i, j)) > 0.1 ||
                     fabs(PixColumnMajor(output->sci.data, i, j)) > 0.1))
                {
                    double tmp = PixColumnMajor(rnz.sci.data, i, j);
                    rmsLocal  +=  tmp*tmp;
                    ++nrmsLocal;
                }
            }
        }//implicit omp barrier

#ifdef _OPENMP
        #pragma omp critical (aggregate)
#endif
        {
            rms  += rmsLocal;
            nrms += nrmsLocal;
        }
#ifdef _OPENMP
        #pragma omp barrier
#endif

#ifdef _OPENMP
        #pragma omp single
        {
            rms = sqrt(rms/nrms);
        }
#endif

        // if it is true that one breaks then it is true for all
        /*epsilon type comparison*/
        if ( (readNoiseAmp - rms) < 0.00001)
            break; // this exits loop over iter
#ifdef _OPENMP
        #pragma omp barrier
#endif
    } // end loop over iter
    } // close parallel block
    freeSingleGroup(&adjustment);
    freeSingleGroup(&rnz);

    if (verbose)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        sprintf(MsgText,"Time taken to smooth image: %.2f(s) with %i procs/threads\n",timeSpent/maxThreads,maxThreads);
        trlmessage(MsgText);
    }

    return (status);
}

double find_dadj(const unsigned i, const unsigned j, const unsigned nRows, const float * obsloc[3], const float * rszloc[3], const double readNoiseAmp)
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

    const double mval = (double)*(rszloc[i] + j);
    const double dval0  = (double)*(obsloc[i] + j) - mval;
    double dval0u = dval0;

    if (dval0u > 1)
        dval0u = 1;
    else if (dval0u < -1)
        dval0u = -1;

    /*COMPARE THE SURROUNDING PIXELS*/
    double dval9 = 0.0;
    if (i == 1 &&  j < nRows-1 && j > 0)
    {
        dval9 = (double)*(obsloc[i]   + j-1)   - (double)*(rszloc[i]   + j-1) +
                (double)*(obsloc[i]   + j)     - (double)*(rszloc[i]   + j)   +
                (double)*(obsloc[i]   + j+1)   - (double)*(rszloc[i]   + j+1) +
                (double)*(obsloc[i-1] + j-1)   - (double)*(rszloc[i-1] + j-1) +
                (double)*(obsloc[i-1] + j)     - (double)*(rszloc[i-1] + j)   +
                (double)*(obsloc[i-1] + j+1)   - (double)*(rszloc[i-1] + j+1) +
                (double)*(obsloc[i+1] + j-1)   - (double)*(rszloc[i+1] + j-1) +
                (double)*(obsloc[i+1] + j)     - (double)*(rszloc[i+1] + j)  +
                (double)*(obsloc[i+1] + j+1)   - (double)*(rszloc[i+1] + j+1);

        dval9 = dval9 / 9.0;
    }
    double dval9u = dval9;

    if (dval9u > readNoiseAmp*0.33)
        dval9u = readNoiseAmp*0.33;
    else if (dval9u < readNoiseAmp*-0.33)
        dval9u = readNoiseAmp*-0.33;

    const double dmod1 = j > 0 ? (double)*(rszloc[i] + j-1) - mval : 0;
    const double dmod2 = j < nRows-1 ? (double)*(rszloc[i] + j+1) - mval : 0;

    double dmod1u = dmod1;
    if (dmod1u > readNoiseAmp*0.33)
        dmod1u = readNoiseAmp*0.33;
    else if (dmod1u < readNoiseAmp*-0.33)
        dmod1u = readNoiseAmp*-0.33;

    double dmod2u = dmod2;
    if (dmod2u > readNoiseAmp*0.33)
        dmod2u = readNoiseAmp*0.33;
    else if (dmod2u < readNoiseAmp*-0.33)
        dmod2u = readNoiseAmp*-0.33;

    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    const double readNoiseAmp2 = readNoiseAmp*readNoiseAmp;
    const double w0 =     dval0 * dval0 / (dval0 * dval0 + 4.0 * readNoiseAmp2);
    const double w9 =     dval9 * dval9 / (dval9 * dval9 + 18.0 * readNoiseAmp2);
    const double w1 = 4 * readNoiseAmp2 / (dmod1 * dmod1 + 4.0 * readNoiseAmp2);
    const double w2 = 4 * readNoiseAmp2 / (dmod2 * dmod2 + 4.0 * readNoiseAmp2);

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    return  dval0u * w0 * 0.25f + /* desire to keep the original pixel value */
            dval9u * w9 * 0.25f + /* desire to keep the original sum over 3x3*/
            dmod1u * w1 * 0.25f + /* desire to get closer to the pixel below*/
            dmod2u * w2 * 0.25f; /* desire to get closer to the pixel above*/
}
