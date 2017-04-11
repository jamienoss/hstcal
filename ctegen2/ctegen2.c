#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"

char MsgText[256];
extern void trlmessage (char *message);

//move these elsewhere
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

    //check ptr isn't already registered?

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
void freeReg(PtrRegister * reg)
{
    if (!reg || reg->length == 0)
        return;

    reg->cursor = 0;
    reg->length = 0;
    // free 'itself'
    reg->freeFunctions[0](reg->ptrs);
    reg->ptrs[0] = NULL;
    reg->freeFunctions[0](reg->freeFunctions);
    reg->freeFunctions[0] = NULL;
}


int inverseCTEBlur(const SingleGroup * input, SingleGroup * output, SingleGroup * trapPixelMap, CTEParamsFast * cte)
{
    //WARNING - assumes column major storage order
    assert(trapPixelMap->sci.data.storageOrder == COLUMNMAJOR);
    assert(input->sci.data.storageOrder == COLUMNMAJOR);
    output->sci.data.storageOrder = COLUMNMAJOR;

    extern int status;

    const unsigned nRows = output->sci.data.ny;
    const unsigned nColumns = output->sci.data.nx;
    const double rnAmp2 = cte->rn_amp * cte->rn_amp;
    const FloatTwoDArray * cteRprof  = &cte->rprof->data;
    const FloatTwoDArray * cteCprof = &cte->cprof->data;

#ifdef _OPENMP
#pragma omp parallel shared(input, output, cte, cteRprof, cteCprof, trapPixelMap)
#endif
{
    //get rid of asserts
    double * model       = malloc(sizeof(*model)*nRows);
    assert(model);
    double * tempModel   = malloc(sizeof(*tempModel)*nRows);
    assert(tempModel);
    double * observed = malloc(sizeof(*observed)*nRows);
    assert(observed);
    float * traps = NULL;

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (unsigned j = 0; j < nColumns; ++j)
    {
        // Can't use memcpy as diff types
        for (unsigned i = 0; i < nRows; ++i)
            observed[i] = PixColumnMajor(input->sci.data,i,j);

        traps = &(PixColumnMajor(trapPixelMap->sci.data, 0, j));
        unsigned NREDO = 0;
        Bool REDO;
        do
        {
            REDO = False; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
            memcpy(model, observed, nRows*sizeof(*observed));

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

        // Update source array
        // Can't use memcpy as arrays of diff types
        for (unsigned i = 0; i < nRows; ++i)
            PixColumnMajor(output->sci.data, i, j) = model[i];
    } //end loop over columns


    delete((void*)&tempModel);
    delete((void*)&observed);
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

int simulatePixelReadout(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
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

int simulateColumnReadout(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
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

int populateTrapPixelMap(SingleGroup * trapPixelMap, CTEParamsFast * cte, const int verbose)
{
   /*SCALE BY 1 UNLESS THE PCTETAB SAYS OTHERWISE, I IS THE PACKET NUM
     THIS IS A SAFETY LOOP INCASE NOT ALL THE COLUMNS ARE POPULATED
     IN THE REFERENCE FILE*/
    /*FOR REFERENCE TO JAY ANDERSON'S CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

          int   iz_data[<cte->nScaleTableColumns>];  column number in raz format
          double scale512[<cte->nScaleTableColumns>];      scaling appropriate at row 512
          double scale1024[<cte->nScaleTableColumns>];     scaling appropriate at row 1024
          double scale1536[<cte->nScaleTableColumns>];     scaling appropriate at row 1536
          double scale2048[<cte->nScaleTableColumns>];     scaling appropriate at row 2048
     */

    //WARNING - assumes column major storage order
    trapPixelMap->sci.data.storageOrder = COLUMNMAJOR;

    clock_t begin = clock();

    extern int status;

    const unsigned nRows = trapPixelMap->sci.data.ny;
    const unsigned nColumns = trapPixelMap->sci.data.nx;
    const double cteScale = cte->scale_frac;


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
    #pragma omp for schedule(static)
#endif

    for (unsigned i = 0; i < cte->nScaleTableColumns; ++i)
    {
        unsigned column = cte->iz_data[i] - cte->razColumnOffset; //which column to scale
        if (column < 0 || column >= nColumns)
            continue;
        trapColumnScale[0] = cte->scale512[i];
        trapColumnScale[1] = cte->scale1024[i];
        trapColumnScale[2] = cte->scale1536[i];
        trapColumnScale[3] = cte->scale2048[i];
        //CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
        //  Index is figured on the final size of the image
        //  not the current size.
        for (unsigned j = 0; j < nRows; ++j)
        {
            ro = j / 512.0; //ro can be zero, it's an index
            if (ro > 2.999)
                ro = 2.999; // only 4 quads, 0 to 3
            else if (ro < 0)
                ro = 0;
            io = (int) floor(ro); //force truncation towards 0 for pos numbers
            cte_j = (j+1) / 2048.0;
            cte_i = trapColumnScale[io] + (trapColumnScale[io+1] - trapColumnScale[io]) * (ro - io);
            PixColumnMajor(trapPixelMap->sci.data, j, column) = cte_i * cte_j * cteScale;
        }
    }
    } // end parallel block

    if (verbose)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        sprintf(MsgText,"Time taken to populate pixel trap map image: %.2f(s) with %i threads\n",timeSpent/cte->maxThreads, cte->maxThreads);
        trlmessage(MsgText);
    }

    return(status);
}

int cteSmoothImage(const SingleGroup * input, SingleGroup * output, CTEParamsFast * ctePars, double ampReadNoise, int verbose)
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

    //WARNING - assumes column major storage order
    assert(input->sci.data.storageOrder == COLUMNMAJOR);
    output->sci.data.storageOrder = COLUMNMAJOR;

    extern int status;

    const unsigned nRows = input->sci.data.ny;
    const unsigned nColumns = input->sci.data.nx;

    double rmsGlobal=0;
    double nrmsGlobal=0;

    clock_t begin = clock();

    copySingleGroup(output, input, input->sci.data.storageOrder);

    //Is the readnoise diff per amp? Current method assumes not.
    if (ampReadNoise < 0.1){
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

    SingleGroup readNoise;
    initSingleGroup(&readNoise);
    allocSingleGroup(&readNoise, nColumns, nRows, False);

#ifdef _OPENMP
    #pragma omp parallel shared(input, output, ampReadNoise, rmsGlobal, nrmsGlobal, readNoise)
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
        #pragma omp for schedule(static)
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
                PixColumnMajor(adjustment.sci.data, j, i) = find_dadjFast(1+i-imid, j, nRows, obs_loc, rsz_loc, ampReadNoise);
        } /*end the parallel for*/ //implicit omp barrier

        //NOW GO OVER ALL THE nColumns AND nRows AGAIN TO SCALE THE PIXELS
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for(unsigned i = 0; i < nColumns; ++i)
        {
            for(unsigned j = 0; j < nRows; ++j)
            {
                PixColumnMajor(output->sci.data,j,i) += (PixColumnMajor(adjustment.sci.data,j, i)*0.75);
                PixColumnMajor(readNoise.sci.data,j,i) = (PixColumnMajor(input->sci.data,j,i) - PixColumnMajor(output->sci.data,j,i));
            }
        }//implicit omp barrier

#ifdef _OPENMP
        #pragma omp single
        {
            rmsGlobal=0;
            nrmsGlobal=0;
        }
#endif

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for(unsigned j = 0; j < nColumns; ++j)
        {
            for(unsigned i = 0; i < nRows; ++i)
            {
                if ( (fabs(PixColumnMajor(input->sci.data, i, j)) > 0.1 ||
                     fabs(PixColumnMajor(output->sci.data, i, j)) > 0.1))
                {
                    double tmp = PixColumnMajor(readNoise.sci.data, i, j);
                    rmsLocal  +=  tmp*tmp;
                    ++nrmsLocal;
                }
            }
        }//implicit omp barrier

#ifdef _OPENMP
        #pragma omp critical (aggregate)
#endif
        {
            rmsGlobal  += rmsLocal;
            nrmsGlobal += nrmsLocal;
        }
#ifdef _OPENMP
        #pragma omp barrier
#endif

#ifdef _OPENMP
        #pragma omp single
        {
            rmsGlobal = sqrt(rmsGlobal/nrmsGlobal);
        }
#endif

        // if it is true that one breaks then it is true for all
        /*epsilon type comparison*/
        //WARNING: remove fabs at some point! Rethink what this even means when computed over the entire image
        //instead of per amp-quad
        if ( fabs(ampReadNoise - rmsGlobal) < 0.00001)
            break; // this exits loop over iter

#ifdef _OPENMP
        #pragma omp barrier
#endif
    } // end loop over iter
    } // close parallel block
    freeSingleGroup(&adjustment);
    freeSingleGroup(&readNoise);

    if (verbose)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        sprintf(MsgText,"Time taken to smooth image: %.2f(s) with %i threads\n", timeSpent/ctePars->maxThreads, ctePars->maxThreads);
        trlmessage(MsgText);
    }

    return (status);
}

double find_dadjFast(const unsigned i, const unsigned j, const unsigned nRows, const float * obsloc[3], const float * rszloc[3], const double readNoiseAmp)
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
    const double readNoiseAmpFraction = 0.33;
    double dval9u = dval9;
    if (dval9u > readNoiseAmp*readNoiseAmpFraction)
        dval9u = readNoiseAmp*readNoiseAmpFraction;
    else if (dval9u < readNoiseAmp*-readNoiseAmpFraction)
        dval9u = readNoiseAmp*-readNoiseAmpFraction;

    const double dmod1 = j > 0 ? (double)*(rszloc[i] + j-1) - mval : 0;
    const double dmod2 = j < nRows-1 ? (double)*(rszloc[i] + j+1) - mval : 0;

    double dmod1u = dmod1;
    if (dmod1u > readNoiseAmp*readNoiseAmpFraction)
        dmod1u = readNoiseAmp*readNoiseAmpFraction;
    else if (dmod1u < readNoiseAmp*-readNoiseAmpFraction)
        dmod1u = readNoiseAmp*-readNoiseAmpFraction;

    double dmod2u = dmod2;
    if (dmod2u > readNoiseAmp*readNoiseAmpFraction)
        dmod2u = readNoiseAmp*readNoiseAmpFraction;
    else if (dmod2u < readNoiseAmp*-readNoiseAmpFraction)
        dmod2u = readNoiseAmp*-readNoiseAmpFraction;

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
