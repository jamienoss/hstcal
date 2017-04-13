/* WFC3 -- CTE loss correction for UVIS

   M. sosey  Aug-2014  Adapted for the pipeline from Jay Andersons CTE correction code for wfc3 UVIS
   raw2raz_wfc3uv.F , an edited file was delivered december 2014, and both are different from the
   fortran code currently served on the wfc3 website.

   M. Sosey Aug-2016 Adapted to be used with Subarrays as well as full frame arrays,
   as long as the subarray contains physical overscan pixels, which don't include the science team subarrays
   which can span quads.
*/

# include <time.h>
# include <string.h>
# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <float.h>
#include <assert.h>

# ifdef _OPENMP
#  include <omp.h>
# endif

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"
# include "wf3corr.h"
# include "cte-fast.h"
# include "../../../../ctegen2/ctegen2.h"


int WF3cteFast (char *input, char *output, CCD_Switch *cte_sw,
        RefFileInfo *refnames, int printtime, int verbose, unsigned nThreads) {

    /*
    input: filename
    output: filename
    cte_sw: the calibration flags
    refnames: the names of the calibration reference files
    onecpu: use parallel processing?

    The following are new primary header keywords which will be added to the data
    so that they can be updated by the code. They are also specified in the PCTETAB
    reference file.

    These are taken from the PCTETAB
    CTE_NAME - name of cte algorithm
    CTE_VER - version number of cte algorithm
    CTEDATE0 - date of wfc3/uvis installation in HST, in MJD
    CTEDATE1 - reference date of CTE model pinning, in MJD

    PCTETLEN - max length of CTE trail
    PCTERNOI - readnoise amplitude for clipping
    PCTESMIT - number of iterations used in CTE forward modeling
    PCTESHFT - number of iterations used in the parallel transfer
    PCTENSMD - readnoise mitigation algorithm
    PCTETRSH - over-subtraction threshold
    PCTEFRAC - cte scaling frac calculated from expstart
    PCTERNOI - the readnoise clipping level to use

    #These are taken from getreffiles.c
    DRKCFILE is a new dark reference file used only in the CTE branch *_DRC.fits
    BIACFILE is a new super-bias reference file used only in the CTE branch *_BIC.fits
    PCTETAB is a new reference file FITS table which will contain the software parameter switches for the CTE correction *_CTE.fit

    This is the main workhorse function for removing the CTE from WFC3 UVIS images

    Unfortunately this happens before anything else in wfc3, so there's a lot of reading files
    at the beginning in order to populate needed information. The rest of the pipeline works
    on one chip at a time and the structures are all defined to support that. None of these
    structures are defined until the code enters the single chip loops. This differs from the
    CTE correction in ACS which occurs later in the process after basic structures are defined.
*/

    extern int status;

    WF3Info wf3; /*structure with calibration switches and reference files for passing*/
    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    clock_t begin = (double)clock();

#ifdef _OPENMP
    trlmessage("Using parallel processing provided by OpenMP inside CTE routine");
    omp_set_dynamic(0);
    unsigned ompMaxThreads = omp_get_num_procs();
    if (nThreads > ompMaxThreads)
    {
        sprintf(MsgText, "System env limiting nThreads from %d to %d", nThreads, ompMaxThreads);
        nThreads = ompMaxThreads;
    }
    else
        sprintf(MsgText,"Setting max threads to %d out of %d available", nThreads, ompMaxThreads);

    omp_set_num_threads(nThreads);
    trlmessage(MsgText);
#endif


    /* COPY COMMAND-LINE ARGUMENTS INTO WF3. */
    WF3Init (&wf3); /*sets default information*/
    strcpy (wf3.input, input);
    strcpy (wf3.output, output);

    PrBegin ("WFC3CTE");
    if (wf3.printtime)
        TimeStamp("WFC3CTE Started: ",wf3.rootname);

    /* CHECK WHETHER THE OUTPUT FILE ALREADY EXISTS. */
    if (FileExists (wf3.output)){
        WhichError(status);
        freeOnExit(&ptrReg);
        return (ERROR_RETURN);
    }

    wf3.pctecorr = cte_sw->pctecorr;
    wf3.darkcorr = cte_sw->darkcorr;
    wf3.biascorr = cte_sw->biascorr;
    wf3.blevcorr = cte_sw->blevcorr;
    wf3.printtime = printtime;
    wf3.verbose = verbose;
    wf3.refnames = refnames;

    PrFileName ("input", wf3.input);
    PrFileName ("output", wf3.output);

    if (wf3.biascorr == COMPLETE){
        trlmessage("BIASCORR complete for input image, CTE can't be performed");
        freeOnExit(&ptrReg);
        return(ERROR_RETURN);
    }
    if (wf3.darkcorr == COMPLETE){
        trlmessage("DARKCORR complete for input image, CTE can't be performed");
        freeOnExit(&ptrReg);
        return(ERROR_RETURN);
    }
    if (wf3.blevcorr == COMPLETE){
        trlmessage("BLEVCORR complete for input image, CTE can't be performed");
        freeOnExit(&ptrReg);
        return(ERROR_RETURN);
    }

    /* DETERMINE THE NAMES OF THE TRAILER FILES BASED ON THE INPUT
       AND OUTPUT FILE NAMES, THEN INITIALIZE THE TRAILER FILE BUFFER
       WITH THOSE NAMES.
       */
    if (initCTETrl (input, output))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
    Hdr phdr; /*primary header for input image, all output information saved here*/
    initHdr(&phdr);
    addPtr(&ptrReg, &phdr, &freeHdr);
    if (LoadHdr (wf3.input, &phdr) ){
        WhichError(status);
        freeOnExit(&ptrReg);
        return (ERROR_RETURN);
    }

    /* GET KEYWORD VALUES FROM PRIMARY HEADER. */
    if ((status = GetKeys (&wf3, &phdr))) {
        freeOnExit(&ptrReg);
        return (status);
    }

    if ((status = GetCTEFlags (&wf3, &phdr))) {
        freeOnExit(&ptrReg);
        return (status);
    }

    /*READ IN THE CTE PARAMETER TABLE*/
    CTEParamsFast cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    initCTEParamsFast(&cte_pars, TRAPS, RAZ_ROWS, RAZ_COLS, RAZ_COLS, nThreads);
    addPtr(&ptrReg, &cte_pars, &freeCTEParamsFast);
    if ((status = allocateCTEParamsFast(&cte_pars)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }
    if ((status = GetCTEParsFast(wf3.pctetab.name, &cte_pars)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }
    //Compute scale fraction
    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
          APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
          INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
          CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
          OVER 3 YEARS AFTER INSTALLATION*/
    cte_pars.scale_frac = (wf3.expstart - cte_pars.cte_date0)/(cte_pars.cte_date1 - cte_pars.cte_date0);

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);

        sprintf(MsgText,"CTE_FF (scaling fraction by date) = %g", cte_pars.scale_frac);
        trlmessage(MsgText);
    }

    //Store sizes - these are corrected for subarrays in getSubarray()
    cte_pars.nRowsPerFullFrame = RAZ_ROWS;
    cte_pars.nColumnsPerFullFrame = RAZ_COLS;
    cte_pars.nColumnsPerChip = cte_pars.nColumnsPerFullFrame/2;
    cte_pars.nRowsPerChip = cte_pars.nRowsPerFullFrame;
    cte_pars.nColumnsPerQuad = cte_pars.nColumnsPerFullFrame/4;
    cte_pars.nRowsPerQuad = cte_pars.nRowsPerFullFrame;
    cte_pars.isSubarray = wf3.subarray;
    cte_pars.refAndIamgeBinsIdenticle = True;

    unsigned nChips = wf3.subarray ? 1 : 2;
    PtrRegister chipLoopReg;
    addPtr(&ptrReg, &chipLoopReg, &freeOnExit);
    for(unsigned chip = 1; chip <= nChips; ++chip)
    {
        //This is used for the final output
        SingleGroup raw;
        initSingleGroup(&raw);
        addPtr(&chipLoopReg, &raw, &freeSingleGroup);

        //Load image into 'raw' one chip at a time
        if (wf3.subarray)
        {
            if ((status = getCCDChipId(&wf3.chip, wf3.input, "SCI", 1)))
            {
                freeOnExit(&ptrReg);
                return status;
            }
            if ((status = getSubarray(&raw, &cte_pars, &wf3)))
            {
                freeOnExit(&ptrReg);
                return status;
            }
        }
        else
        {
            cte_pars.columnOffset = 0;
            cte_pars.rowOffset = 0;
            if (chip == 1)
                cte_pars.razColumnOffset = 0;
            else if (chip == 2)
                cte_pars.razColumnOffset = cte_pars.nColumnsPerChip;
            else
                assert(0);//unimplemented

            getSingleGroup(wf3.input, chip, &raw);
            if (hstio_err())
            {
                freeOnExit(&ptrReg);
                return (status = OPEN_FAILED);
            }
        }

        cte_pars.nRows = raw.sci.data.ny;
        cte_pars.nColumns = raw.sci.data.nx;
        findAlignedQuadImageBoundaries(&cte_pars, 25, 30, 19); //25 prescan, 30 postscan, & 19 parallel overscan

        //leave raw as pre-biased image, clone and use copy from here on out
        SingleGroup rowMajorImage;
        initSingleGroup(&rowMajorImage);
        addPtr(&chipLoopReg, &rowMajorImage, &freeSingleGroup);
        if ((status = allocSingleGroupExts(&rowMajorImage, cte_pars.nColumns, cte_pars.nRows, SCIEXT, False)))
        {
            sprintf(MsgText, "Allocation problem with 'rowMajorImage' in 'WF3cteFast()'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return status;
        }
        SingleGroup * image = &rowMajorImage;
        if ((status = copySingleGroup(image, &raw, raw.sci.data.storageOrder)))
        {
            freeOnExit(&ptrReg);
            return status;
        }
        //align raw image for later comparison with aligned corrected image
        alignAmps(&raw, &cte_pars);

        //biac bias subtraction
        if (doCteBias(&wf3, image))
        {
            freeOnExit(&ptrReg);
            return(status);
        }
        alignAmps(image, &cte_pars);

        //CTE correction not sensitive enough to work without amp bias and gain correction
        //which require vertical overscan
        if (cte_pars.isSubarray)
        {
            for (unsigned quad = 0; quad < 2; ++quad)
            {
                if (cte_pars.quadExists[quad] && !cte_pars.hasPrescan[quad])
                {
                    sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",
                            cte_pars.columnOffset, cte_pars.columnOffset + cte_pars.nColumns);
                    trlmessage(MsgText);
                    freeOnExit(&ptrReg);
                    return(ERROR_RETURN);
                }
            }
        }

        //c++ ref would be good here
        const unsigned nRows = cte_pars.nRows;
        const unsigned nColumns = cte_pars.nColumns;

        //copy to column major storage
        SingleGroup columnMajorImage;
        initSingleGroup(&columnMajorImage);
        addPtr(&chipLoopReg, &columnMajorImage, &freeSingleGroup);
        if ((status = allocSingleGroupExts(&columnMajorImage, nColumns, nRows, SCIEXT, False)))
        {
            sprintf(MsgText, "Allocation problem with 'columnMajorImage' in 'WF3cteFast()'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return status;
        }
        if ((status = copySingleGroup(&columnMajorImage, image, COLUMNMAJOR)))
        {
            sprintf(MsgText, "Allocation problem with 'copySingleGroup(columnMajorImage)' in 'WF3cteFast()'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return status;
        }
        image = &columnMajorImage;

        //SUBTRACT AMP BIAS AND CORRECT FOR GAIN
        if (correctAmpBiasAndGain(image, wf3.ccdgain, &cte_pars))
        {
            freeOnExit(&ptrReg);
            return (status);
        }

        /*//Subtract these now - not before alignAmps()
                if (wf3.subarray)
                {
                    cte_pars.nColumnsPerChip -= 2*cte_pars.postscanWidth;
                    cte_pars.nColumnsPerQuad -= cte_pars.postscanWidth;
                }
*/
        /***CALCULATE THE SMOOTH READNOISE IMAGE***/
        trlmessage("CTE: Calculating smooth readnoise image");
        SingleGroup * smoothedImage = &rowMajorImage; //reuse rowMajorImage memory space
        setStorageOrder(smoothedImage, COLUMNMAJOR);




        //WARNING TEST CODE ONLY
       // copySingleGroup(smoothedImage, image, COLUMNMAJOR);



        /***CREATE THE NOISE MITIGATION MODEL ***/
        if (cte_pars.noise_mit == 0)
        {
            //printf("smooth clipping level: %f\n",cte_pars.rn_amp);
            if (cteSmoothImage(image, smoothedImage, &cte_pars, cte_pars.rn_amp, wf3.verbose))
            {
                freeOnExit(&ptrReg);
                return (status);
            }
        }
        else
        {
            trlmessage("Only noise model 0 implemented!");
            freeOnExit(&ptrReg);
            return (status=ERROR_RETURN);
        }

        SingleGroup trapPixelMap;
        initSingleGroup(&trapPixelMap);
        addPtr(&chipLoopReg, &trapPixelMap, &freeSingleGroup);
        if ((status = allocSingleGroupExts(&trapPixelMap, nColumns, nRows, SCIEXT, False)))
        {
           sprintf(MsgText, "Allocation problem with 'trapPixelMap' in 'WF3cteFast()'");
           trlerror(MsgText);
           freeOnExit(&ptrReg);
           return status;
       }
        setStorageOrder(&trapPixelMap, COLUMNMAJOR);
        if (populateTrapPixelMap(&trapPixelMap, &cte_pars, wf3.verbose))
        {
            freeOnExit(&ptrReg);
            return status;
        }

        SingleGroup * cteCorrectedImage = image; // reuse columnMajorImage
        image = NULL;



        //TESTcode only
       // copySingleGroup(cteCorrectedImage, smoothedImage, COLUMNMAJOR);




        // MAIN CORRECTION LOOP IN HERE
        if (inverseCTEBlur(smoothedImage, cteCorrectedImage, &trapPixelMap, &cte_pars))
        {
            freeOnExit(&ptrReg);
            return status;
        }
        freePtr(&chipLoopReg, &trapPixelMap);

        const double scaleFraction = cte_pars.scale_frac;
        //freePtr(&ptrReg, &cte_pars);

        // CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT
        const float ccdgain = wf3.ccdgain;
        float totalCounts = 0;
        float totalRawCounts = 0;
    #ifdef _OPENMP
        #pragma omp parallel shared(cteCorrectedImage, smoothedImage, raw)
    #endif
        {
        float threadCounts = 0;
        float threadRawCounts = 0;
        float delta;
    #ifdef _OPENMP
        #pragma omp for schedule(static)
    #endif
        for (unsigned i = 0; i < nColumns; ++i)
        {
            for(unsigned j = 0; j < nRows; ++j)
            {
                delta = (PixColumnMajor(cteCorrectedImage->sci.data,j,i) - PixColumnMajor(smoothedImage->sci.data,j,i))/ccdgain;
                threadCounts += delta;
                threadRawCounts += Pix(raw.sci.data, i, j);
                Pix(raw.sci.data, i, j) += delta;
            }
        }

    #ifdef _OPENMP
        #pragma omp critical(deltaAggregate)
    #endif
        {
            totalCounts += threadCounts;
            totalRawCounts += threadRawCounts;
        }
        }//end omp threads
        printf("\nTotal count difference (rac-raw) incurred from correction: %f (%f%%)\n\n", totalCounts, totalCounts/totalRawCounts*100);

        cteCorrectedImage = NULL;
        smoothedImage = NULL;
        freePtr(&chipLoopReg, &rowMajorImage);
        freePtr(&chipLoopReg, &columnMajorImage);

     /*   //Put back this alteration before calling alignAmps again (in outputImage())
        if (wf3.subarray)
        {
            cte_pars.nColumnsPerChip += 2*cte_pars.postscanWidth;
            cte_pars.nColumnsPerQuad += cte_pars.postscanWidth;
        }
        */
        if (outputImage(output, &raw, &cte_pars))//needs to pop status
        {
            freeOnExit(&ptrReg);
            return status;
        }

        // SAVE USEFUL HEADER INFORMATION
        if (chip == 2)
        {
            if (cteHistory(&wf3, raw.globalhdr))
            {
                freeOnExit(&ptrReg);
                return status;
            }
            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(raw.globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");
        }

        double time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
        if (verbose){
            sprintf(MsgText,"CTE run time: %.2f(s) with %i procs/threads\n",time_spent/nThreads,nThreads);
            trlmessage(MsgText);
        }
        freeAll(&chipLoopReg);
    }//end of chip for loop

    PrSwitch("pctecorr", COMPLETE);
    if (wf3.printtime)
        TimeStamp("PCTECORR Finished", wf3.rootname);

    freeOnExit(&ptrReg);
    return (status);
}


/********************* SUPPORTING SUBROUTINES *****************************/
int correctAmpBiasAndGain(SingleGroup * image, const float ccdGain, CTEParamsFast * ctePars)
{
    /* Do an additional bias correction using the residual bias level measured for each amplifier from the
     * steadiest pixels in the horizontal overscan and subtracted from the pixels for that amplifier.
     */

    //WARNING - assumes column major storage order
    assert(image->sci.data.storageOrder == COLUMNMAJOR);

    extern int status;

    float biasMean[2] = {0, 0};
    float biasSigma[2]= {0, 0}; // This is not actually used (dummy to pass to findOverScanBias)

    enum OverscanType overscanType = ctePars->isSubarray ? PRESCAN : POSTSCAN;
    unsigned nOverscanColumnsToIgnore = ctePars->isSubarray ? 5 : 3;
    findOverscanBias(image, biasMean, biasSigma, overscanType, nOverscanColumnsToIgnore, ctePars);

    //used to vary for dev purposes
    unsigned rowsStart = 0;//ctePars->imageRowsStart;
    unsigned rowsEnd = image->sci.data.ny;//ctePars->imageRowsEnd;
    unsigned columnsStart[2] = {0, 2103};//{ctePars->imageColumnsStart[0], ctePars->imageColumnsStart[1]};
    unsigned columnsEnd[2] = {2103, 2103*2};//{ctePars->imageColumnsEnd[0], ctePars->imageColumnsEnd[1]};

#ifdef _OPENMP
    #pragma omp parallel shared(image, biasMean, biasSigma)
#endif
    {
    //Image will only ever be at most 2 amps i.e. a single chip
    // SUBTRACT THE EXTRA BIAS CALCULATED, AND MULTIPLY BY THE GAIN
    for (unsigned nthAmp = 0; nthAmp < 2; ++nthAmp)
    {
        //Check if amp even present
        if (!ctePars->quadExists[nthAmp])
            continue;
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        //NOTE: this nolonger overwrites overscan regions!!!
        for (unsigned i = columnsStart[nthAmp]; i < columnsEnd[nthAmp]; ++i)
        {
            for (unsigned j = rowsStart; j < rowsEnd; ++j)
            {
                PixColumnMajor(image->sci.data, j, i) -= biasMean[nthAmp];
                PixColumnMajor(image->sci.data, j, i) *= ccdGain;
            }
        }
    }
    } // close parallel block

    return(status);
}

int findOverscanBias(SingleGroup *image, float *mean, float *sigma, enum OverscanType overscanType, unsigned nOverscanColumnsToIgnore, CTEParamsFast * ctePars)
{
    //WARNING - assumes column major storage order
    assert(image->sci.data.storageOrder == COLUMNMAJOR);

    /*Calculate the post scan and bias after the biac file has been subtracted.
      This calls resistmean, which does a better job clipping outlying pixels
      that just a standard stddev clip single pass.

      The first 25 columns in a quad are serial physical overscan (prescan)
      The last 30 columns in a quad are parallel virtual overscan and are only present in full frame quads (postscan)
      The last 19 rows in a quad are serial
      Actual image size = 2051 (rows) x 2048 (columns)
     */

    extern int status;

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    for (unsigned nthAmp = 0; nthAmp < 2; ++nthAmp)
    {
        if (!ctePars->quadExists[nthAmp])
            continue;

        float * imageOverscanPixels = NULL;
        unsigned nOverscanPixels = 0;
        unsigned overscanWidth = 0;
        unsigned nOverscanRows = ctePars->imageRowsEnd;

        //Find overscan columns
        if (overscanType == PRESCAN)//subarray
        {
            if (!ctePars->hasPrescan[nthAmp])
                continue;
            unsigned overscanStart = ctePars->columnOffset < nOverscanColumnsToIgnore ? nOverscanColumnsToIgnore - ctePars->columnOffset : 0;
            overscanWidth = ctePars->imageColumnsStart[nthAmp] - overscanStart;
            nOverscanPixels = overscanWidth*nOverscanRows;
            imageOverscanPixels = image->sci.data.data + overscanStart*ctePars->nRows;
        }
        else if (overscanType == POSTSCAN)//full frame
        {
            if (!ctePars->hasPostscan[nthAmp])
                continue;
            overscanWidth = ctePars->postscanWidth - nOverscanColumnsToIgnore;
            nOverscanPixels = overscanWidth*nOverscanRows;
            imageOverscanPixels = image->sci.data.data + (ctePars->imageColumnsEnd[nthAmp] + nOverscanColumnsToIgnore-1)*ctePars->nRows;
        }
        else
            assert(0);

        //If we didn't need to skip the 19 rows of parallel virtual overscan, this array would
        //not be needed and a pointer to the data could be passed directly to resistmean instead.
        float * overscanPixels = malloc(nOverscanPixels*sizeof(*overscanPixels));
        addPtr(&ptrReg, overscanPixels, &free);
        if (!overscanPixels)
        {
            sprintf(MsgText, "Out of memory for 'overscanPixels' in 'findOverscanBias'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return (status = OUT_OF_MEMORY);
        }

        assert(overscanPixels);
        for (unsigned column = 0; column < overscanWidth; ++column)
            memcpy(overscanPixels + column*nOverscanRows, imageOverscanPixels + column*ctePars->nRows, nOverscanRows*sizeof(*overscanPixels));

        float rmean = 0;
        float rsigma = 0;
        float min = 0;
        float max = 0;
        resistmean(overscanPixels, nOverscanPixels, 7.5, &rmean, &rsigma, &min, &max);
        mean[nthAmp] = rmean;
        sigma[nthAmp] = rsigma;

        printf("npix=%i\nmean[%i]=%f\nsigma[%i] = %f\n", nOverscanPixels, nthAmp+1, rmean, nthAmp+1, rsigma);

        freePtr(&ptrReg, overscanPixels);
    }

    freeOnExit(&ptrReg);
    return status;
}

//NOTE: wf3 * ctePars should be const however this is too large a refactor for this PR
int getSubarray(SingleGroup * image, CTEParamsFast * ctePars, WF3Info * wf3)
{
    extern int status;

    if (!ctePars->isSubarray)
        return status;

    // get subarray from first extension
    getSingleGroup(wf3->input, 1, image);
    if (hstio_err()){
        freeSingleGroup(image);
        return (status = OPEN_FAILED);
    }
    // getSingleGroup will set image->group_num to 1 here so now correct
    image->group_num = wf3->chip == 2 ? 1 : 2;

    ctePars->nRows = image->sci.data.ny;
    ctePars->nColumns = image->sci.data.nx;

    // These are used to find subarrays location within a full chip (2 quads)
    int sci_bin[2];         // bin size of science image
    int sci_corner[2];      // science image corner location
    int ref_bin[2];
    int ref_corner[2];
    int rsize = 1;          // reference pixel size

    if (GetCorner(&image->sci.hdr, rsize, sci_bin, sci_corner))
    {
        freeSingleGroup(image);
        return (status);
    }
    //Create a dummy header to represent the full chip (both amps) and populate via initChipMetaData
    {
        Hdr fullChipHdr;
        initChipMetaData(wf3, &fullChipHdr, image->group_num);
        if (GetCorner (&fullChipHdr, rsize, ref_bin, ref_corner))
        {
            freeSingleGroup(image);
            return (status);
        }
    }

    ctePars->refAndIamgeBinsIdenticle = ref_bin[0] == sci_bin[0] && ref_bin[1] == sci_bin[1] ? True : False;
    ctePars->columnOffset = sci_corner[0] - ref_corner[0];
    ctePars->rowOffset = sci_corner[1] - ref_corner[1];

    /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
      AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
      CTE CODE VARIABLES.
      */
    if (CompareCTEParamsFast(image, ctePars))
    {
        freeSingleGroup(image);
        return (status);
    }

    return status;
}

int unalignAmps(SingleGroup * image, CTEParamsFast * ctePars)
{
    //The function alignAmps is symmetric, wrapping as aid when reading logic
    return alignAmps(image, ctePars);
}

int alignAmps(SingleGroup * image, CTEParamsFast * ctePars)
{
    //NOTE: There is a similar version of this in acs - code changes should be reflected in both.

    if (!image || !image->sci.data.data)
        return status;

    //WARNING - assumes row major storage
    assert(image->sci.data.storageOrder == ROWMAJOR);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    //Align amps such that they are at the bottom left
    extern int status;

    const unsigned columnOffset = ctePars->columnOffset;
    const unsigned nColumns = ctePars->nColumns;
    const unsigned nRows = ctePars->nRows;

    const Bool isCDAmp = image->group_num == 1 ? True : False;
    if (ctePars->isSubarray)
    {
        if (isCDAmp)
            printf("subarray from CD amp\n");
        else
            printf("subarray from AB amp\n");
    }

    //If subarray only in c quad do nothing as this is already bottom left aligned
    if (isCDAmp && !ctePars->quadExists[1])
    {
        freeOnExit(&ptrReg);
        return status;
    }

    //Find if subarray extends into either b or d quad and flip right to left
    if (ctePars->quadExists[1])//  &columnOffset + nColumns > ctePars->nColumnsPerQuad)
    {
        if (ctePars->isSubarray)
            printf("subarray extends into amps B or D\n");
        //grab a row, flip it, put it back
        unsigned rowLength = columnOffset + nColumns - ctePars->nColumnsPerQuad;
        unsigned quadBoundary = nColumns - rowLength;
#ifdef _OPENMP
        #pragma omp parallel for shared(image, ctePars) schedule(static)
#endif
        for (unsigned i = 0; i < nRows; ++i)
        {
            //find row
            float * row = image->sci.data.data + i*nColumns + quadBoundary;
            //flip right to left
            float tempPixel;
            for (unsigned j = 0; j < rowLength/2; ++j)
            {
                tempPixel = row[j];
                row[j] = row[rowLength-1-j];
                row[rowLength-1-j] = tempPixel;
            }
        }
    }

    //Only thing left is to flip ab chip upside down
    if (!isCDAmp) // isABAmp
    {
        //either physically align all or propagate throughout a mechanism to work on the array upside down (flip b quad though)
        //we'll just flip all for now. See if there's info in the header specifying amp location rel to pix in file,
        //i.e. a way to know whether the chip is 'upside down'. Could then reverse cte corr trail walk direction
        float * tempRow = NULL;
        size_t rowSize = nColumns*sizeof(*tempRow);
        tempRow = malloc(rowSize);
        addPtr(&ptrReg, tempRow, &free);
        if (!tempRow)
        {
            sprintf(MsgText, "Out of memory for 'tempRow' in 'alignAmpData'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return (status = OUT_OF_MEMORY);
        }

#ifdef _OPENMP
        #pragma omp parallel for shared(image, ctePars) schedule(static)
#endif
        for (unsigned i = 0; i < nRows/2; ++i)
        {
            float * topRow = image->sci.data.data + i*nColumns;
            float * bottomRow = image->sci.data.data + (nRows-1-i)*nColumns;
            memcpy(tempRow, topRow, rowSize);
            memcpy(topRow, bottomRow, rowSize);
            memcpy(bottomRow, tempRow, rowSize);
        }
        free(tempRow);
    }

    freeOnExit(&ptrReg);
    return status;
}


int putChip(char * fileName, SingleGroup * image, WF3Info * wf3, double const scaleFraction)
{
    /*** SAVE USEFUL HEADER INFORMATION ***/
    if (cteHistory(wf3, image->globalhdr))
        return status;

    /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
    PutKeyDbl(image->globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
    trlmessage("PCTEFRAC saved to header");

    //Full2Sub(wf3, subChip, fullChip, 0, 1, 1);
    putSingleGroup(fileName, 1, image, 0);
    //freeSingleGroup(subChip);
    return status;
}

int getCCDChipId(int * value, char * fileName, char * ename, int ever)
{
    extern int status;
    // OPEN INPUT IMAGE IN ORDER TO READ ITS SCIENCE HEADER.
    IODescPtr ip = openInputImage (fileName, "SCI", ever);
    if (hstio_err()) {
        sprintf (MsgText, "Image: \"%s\" is not present", fileName);
        trlerror (MsgText);
        return (status = OPEN_FAILED);
    }
    Hdr scihdr;
    initHdr(&scihdr);
    getHeader(ip, &scihdr);
    if (ip)
        closeImage (ip);
    /* Get CCD-specific parameters. */
    if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, ever, value)){
        freeHdr(&scihdr);
        return (status);
    }
    freeHdr(&scihdr);

    return status;
}

int outputImage(char * fileName, SingleGroup * image, CTEParamsFast * ctePars)
{
    extern int status;

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    SingleGroup temp;
    initSingleGroup(&temp);
    if (image->sci.data.storageOrder == COLUMNMAJOR)
    {
        if ((status = allocSingleGroup(&temp, image->sci.data.nx, image->sci.data.ny, False)))
        {
            sprintf(MsgText, "Allocation problem with 'temp' in 'outputImage()'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return status;
        }
        addPtr(&ptrReg, &temp, &freeSingleGroup);
        if ((status = copySingleGroup(&temp, image, ROWMAJOR)))
        {
            sprintf(MsgText, "Allocation problem with 'copySingleGroup(temp)' in 'outputImage()'");
            trlerror(MsgText);
            freeOnExit(&ptrReg);
            return status;
        }
        image = &temp;
    }

    if ((status = unalignAmps(image, ctePars)))
    {
        sprintf(MsgText, "Allocation problem with 'unalignAmps()' in 'outputImage()'");
        trlerror(MsgText);
        freeOnExit(&ptrReg);
        return status;
    }

    if (putSingleGroup(fileName, image->group_num, image, 0))
    {
        sprintf(MsgText, "IO failure in 'putSingleGroup()' in 'outputImage()'");
        trlerror(MsgText);
        freeOnExit(&ptrReg);
        return (status = OPEN_FAILED);
    }
    freeOnExit(&ptrReg);
    return status;
}

void findAlignedQuadImageBoundaries(CTEParamsFast * ctePars, unsigned const prescanWidth, unsigned const postscanWidth, unsigned const parallelOverscanWidth)
{
    //WARNING this assumes quads have already been aligned

    ctePars->prescanWidth = prescanWidth;
    ctePars->postscanWidth = postscanWidth;
    ctePars->parallelOverscanWidth = parallelOverscanWidth;
    ctePars->imageRowsStart = 0;
    //Ignore last 19 rows of parallel virtual overscan, i.e. the last 19 pixels in a column
    ctePars->imageRowsEnd = ctePars->nRows <= ctePars->nRowsPerQuad - parallelOverscanWidth ? ctePars->nRows : ctePars->nRowsPerQuad - parallelOverscanWidth;

    //find image boundaries
    const unsigned nColumns = ctePars->nColumns;
    const unsigned columnOffset = ctePars->columnOffset;
    const unsigned nColumnsPerQuad = ctePars->nColumnsPerQuad;

    if (ctePars->isSubarray)
    {
        ctePars->hasPostscan[0] = False; //always true
        ctePars->hasPostscan[1] = False; //always true
    }
    else
    {
        ctePars->quadExists[0] = True;
        ctePars->quadExists[1] = True;
        ctePars->hasPrescan[0] = True;
        ctePars->hasPrescan[1] = True;
        ctePars->hasPostscan[0] = True;
        ctePars->hasPostscan[1] = True;
        ctePars->imageColumnsStart[0] = prescanWidth;
        ctePars->imageColumnsEnd[0] = nColumnsPerQuad - postscanWidth;
        ctePars->imageColumnsStart[1] = nColumnsPerQuad + prescanWidth; //skip 30 pixels of previous quads virtual overscan
        ctePars->imageColumnsEnd[1] = ctePars->nColumnsPerChip - postscanWidth;
        return;
    }

    //NOTE: For subarrays nColumnsPerChip & nColumnsPerQuad do NOT include the 60 & 30 extra postscan columns
    if (columnOffset < nColumnsPerQuad) //image starts in 1st quad A or C
    {
        ctePars->quadExists[0] = True;
        if (columnOffset < prescanWidth)
        {
            ctePars->hasPrescan[0] = True;
            ctePars->imageColumnsStart[0] = prescanWidth - columnOffset;
        }
        else
        {
            ctePars->hasPrescan[0] = False;
            ctePars->imageColumnsStart[0] = 0;
        }

        if (columnOffset + nColumns > nColumnsPerQuad) //image also extends into 2nd quad B or D
        {
            //NOTE: This code assumes that the extension into the 2nd quad does so beyond the prescan
            //i.e. there is actually an image in the 2nd quad
            assert(columnOffset + nColumns > nColumnsPerQuad + prescanWidth);

            ctePars->quadExists[1] = True;
            ctePars->hasPrescan[1] = True;
            ctePars->imageColumnsEnd[0] = nColumnsPerQuad - columnOffset;
            ctePars->imageColumnsStart[1] = ctePars->imageColumnsEnd[0] + prescanWidth;
            ctePars->imageColumnsEnd[1] = nColumns;
        }
        else
        {
            ctePars->quadExists[1] = False;
            ctePars->hasPrescan[1] = False;
            ctePars->imageColumnsEnd[0] = nColumns;
            ctePars->imageColumnsStart[1] = 0;
            ctePars->imageColumnsEnd[1] = 0;
        }
   }
   else if (columnOffset > nColumnsPerQuad) //image only exits in 2nd quad B or D
   {
       if (columnOffset - nColumnsPerQuad < prescanWidth)
       {
           ctePars->hasPrescan[0] = True;
           ctePars->imageColumnsStart[0] = prescanWidth - (columnOffset - nColumnsPerQuad);
       }
       else
       {
           ctePars->hasPrescan[0] = False;
           ctePars->imageColumnsStart[0] = 0;
       }
   }
   else
       assert(0); //image is subarray but exists in neither quad???
}
