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
# include "wf3err.h"
# include "wf3corr.h"
# include "cte.h"
# include "../../../../ctegen2/ctegen2.h"


int WF3cte (char *input, char *output, CCD_Switch *cte_sw,
        RefFileInfo *refnames, int printtime, int verbose, int onecpu) {

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
    int max_threads=1;
    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    /*check if this is a subarray image.
      This is necessary because the CTE routine will start with the raw images
      from scratch and read them in so that both chips can be used. CTE is
      outside of the normal processing where one chip goes through the pipeline
      at a time, both chips are used at the same time for the correction.

      For the case of subarrays, a fake second chip needs to be created.
      The subarray is also placed inside the confines of a full size image
      and a mask is created to ignore pixels not associated with the original
      data during the cte correction. This is necessary because the pixel location
      itself is used as part of the correction. A secondary option would be to set
      the looping arrays to variable sizes and make sure all array references were
      consistent with the current data being processed. I decided on masking which
      might allow for other considerations in future updates.

      Only subarrays which were taken with physical overscan pixels are currently valid
      This distinction can be made with the CRDS ruleset for PCTECORR but it
      should also be checked here incase users update the header themselves for
      local runs. In order to check for overscan pixels I'm using the array start
      location instead of the APERTURE keyword information (there are known user
      apertures which do not have overscan pixels, but this gets around string
      comparisons and any future name changes or aperture additions in the future)
     */
    clock_t begin = (double)clock();

    /*CONTAIN PARALLEL PROCESSING TO A SINGLE THREAD AS USER OPTION*/
#   ifdef _OPENMP
    trlmessage("Using parallel processing provided by OpenMP inside CTE routine");
    if (onecpu){
        omp_set_dynamic(0);
        max_threads=1;
        sprintf(MsgText,"onecpu == TRUE, Using only %i threads/cpu", max_threads);
    } else {
        omp_set_dynamic(0);
        max_threads = omp_get_num_procs(); /*be nice, use 1 less than avail?*/
        sprintf(MsgText,"Setting max threads to %i of %i cpus",max_threads, omp_get_num_procs());
    }
    omp_set_num_threads(max_threads);
    trlmessage(MsgText);
#   endif


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
        freeAll(&ptrReg);
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
        freeAll(&ptrReg);
        return(ERROR_RETURN);
    }
    if (wf3.darkcorr == COMPLETE){
        trlmessage("DARKCORR complete for input image, CTE can't be performed");
        freeAll(&ptrReg);
        return(ERROR_RETURN);
    }
    if (wf3.blevcorr == COMPLETE){
        trlmessage("BLEVCORR complete for input image, CTE can't be performed");
        freeAll(&ptrReg);
        return(ERROR_RETURN);
    }

    /* DETERMINE THE NAMES OF THE TRAILER FILES BASED ON THE INPUT
       AND OUTPUT FILE NAMES, THEN INITIALIZE THE TRAILER FILE BUFFER
       WITH THOSE NAMES.
       */
    if (initCTETrl (input, output))
    {
        freeAll(&ptrReg);
        return (status);
    }

    /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
    Hdr phdr; /*primary header for input image, all output information saved here*/
    initHdr(&phdr);
    addPtr(&ptrReg, &phdr, &freeHdr);
    if (LoadHdr (wf3.input, &phdr) ){
        WhichError(status);
        freeAll(&ptrReg);
        return (ERROR_RETURN);
    }

    /* GET KEYWORD VALUES FROM PRIMARY HEADER. */
    if (GetKeys (&wf3, &phdr)) {
        freeAll(&ptrReg);
        return (status);
    }

    if (GetCTEFlags (&wf3, &phdr)) {
        freeAll(&ptrReg);
        return (status);
    }

    /*READ IN THE CTE PARAMETER TABLE*/
    CTEParams cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    initCTEParams(&cte_pars, TRAPS, RAZ_ROWS, RAZ_COLS);
    addPtr(&ptrReg, &cte_pars, &freeCTEParams);
    allocateCTEParams(&cte_pars);
    if (GetCTEPars (wf3.pctetab.name, &cte_pars))
    {
        freeAll(&ptrReg);
        return (status);
    }

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);
    }

    //stash sizes
    cte_pars.nRowsPerFullFrame = RAZ_ROWS;
    cte_pars.nColumnsPerFullFrame = RAZ_COLS;
    cte_pars.nColumnsPerChip = cte_pars.nColumnsPerFullFrame/2;
    cte_pars.nRowsPerChip = cte_pars.nRowsPerFullFrame;
    cte_pars.nColumnsPerQuad = cte_pars.nColumnsPerFullFrame/4;
    cte_pars.nRowsPerQuad = cte_pars.nRowsPerFullFrame;

    cte_pars.nRows = cte_pars.nRowsPerFullFrame; // assume full frame for now
    cte_pars.nColumns = cte_pars.nColumnsPerFullFrame; // assume full frame for now

    //This is used for the final output
    SingleGroup raw;
    initSingleGroup(&raw);
    addPtr(&ptrReg, &raw, &freeSingleGroup);

    SingleGroup rowMajorImage;
    initSingleGroup(&rowMajorImage);
    addPtr(&ptrReg, &rowMajorImage, &freeSingleGroup);
    SingleGroup * image = &rowMajorImage;

    if (wf3.subarray)
    {
        if (getCCDChip(&wf3.chip, wf3.input, "SCI", 1) ||
                getSubarray(&raw, &cte_pars, &wf3))
        {
            freeAll(&ptrReg);
            return status;
        }

        alignAmps(&raw, &cte_pars);

        //leave raw as pre-biased image, bifurcate and use copy from here on out
        allocSingleGroupSciOnly(&raw, cte_pars.nColumns, cte_pars.nRows, False);
        copySingleGroup(image, &raw, raw.sci.data.storageOrder);

        //biac bias
       /* if (doCteBias(&wf3, image))
        {
            freeAll(&ptrReg);
            return(status);
        }*/
    }
   /* else
    {
        // Full frame image, just read in the groups
        //   and init the mask to use all pixels

        getSingleGroup(wf3.input, 1, NULL);
        getSingleGroup(wf3.input, 2, NULL);

        // SAVE A COPY OF THE RAW IMAGE FOR LATER
        makeRAZ(&cd,&ab,&raw);

        //SUBTRACT THE CTE BIAS FROM BOTH CHIPS IN PLACE
        if (doCteBias(&wf3,&cd)){
            freeSingleGroup(&cd);
            return(status);
        }
        if (doCteBias(&wf3,&ab)){
            freeSingleGroup(&ab);
            return(status);
        }
        //SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
        //  AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
        //  CTE CODE VARIABLES.

        if (CompareCTEParams(&cd, &cte_pars))
            return (status);
    }
*/
    const unsigned nRows = cte_pars.nRows;
    const unsigned nColumns = cte_pars.nColumns;

    //copy to column major storage
    SingleGroup columnMajorImage;
    initSingleGroup(&columnMajorImage);
    addPtr(&ptrReg, &columnMajorImage, &freeSingleGroup);
    allocSingleGroupSciOnly(&columnMajorImage, nColumns, nRows, False);
    assert(!copySingleGroup(&columnMajorImage, image, COLUMNMAJOR));
    image = &columnMajorImage;

    //SUBTRACT BIAS AND CORRECT FOR GAIN
    if (biasAndGainCorrect(image, wf3.ccdgain, wf3.subarray))
    {
        freeAll(&ptrReg);
        return (status);
    }

    /***CALCULATE THE SMOOTH READNOISE IMAGE***/
    trlmessage("CTE: Calculating smooth readnoise image");

    SingleGroup * smoothedImage = &rowMajorImage; //reuse rowMajorImage memory space
    setStorageOrder(smoothedImage, COLUMNMAJOR);
    /***CREATE THE NOISE MITIGATION MODEL ***/
    if (cte_pars.noise_mit == 0)
    {
        if (cteSmoothImage(image, smoothedImage, cte_pars.rn_amp, max_threads, wf3.verbose))
        {
            freeAll(&ptrReg);
            return (status);
        }
    }
    else
    {
        trlmessage("Only noise model 0 implemented!");
        freeAll(&ptrReg);
        return (status=ERROR_RETURN);
    }

    SingleGroup trapPixelMap;
    initSingleGroup(&trapPixelMap);
    addPtr(&ptrReg, &trapPixelMap, &freeSingleGroup);
    allocSingleGroupSciOnly(&trapPixelMap, nColumns, nRows, False);
    setStorageOrder(&trapPixelMap, COLUMNMAJOR);
    if (populateTrapPixelMap(&trapPixelMap, &cte_pars, wf3.verbose, wf3.expstart))
    {
        freeAll(&ptrReg);
        return status;
    }

    SingleGroup * cteCorrectedImage = image; // reuse columnMajorImage
    image = NULL;
    // MAIN CORRECTION LOOP IN HERE
    if (inverseCTEBlur(smoothedImage, cteCorrectedImage, &trapPixelMap, &cte_pars))
    {
        freeAll(&ptrReg);
        return status;
    }
    freePtr(&ptrReg, &trapPixelMap);

    const double scaleFraction = cte_pars.scale_frac;
    freePtr(&ptrReg, &cte_pars);

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
    freePtr(&ptrReg, &rowMajorImage);
    freePtr(&ptrReg, &columnMajorImage);

    /* COPY BACK THE SCIENCE SUBARRAYS AND
       SAVE THE NEW RAW FILE WITH UPDATED SCIENCE
       ARRAYS AND PRIMARY HEADER TO RAC
       */
    if (outputImage(output, &raw, &cte_pars))
    {
        freeAll(&ptrReg);
        return status;
    }

   /*
        putSingleGroup(output,cd.group_num, &cd,0);
        putSingleGroup(output,ab.group_num, &ab,0);
    }
*/
    /*** SAVE USEFUL HEADER INFORMATION ***/
    if (cteHistory(&wf3, raw.globalhdr))
    {
        freeAll(&ptrReg);
        return status;
    }
    /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
    PutKeyDbl(raw.globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
    trlmessage("PCTEFRAC saved to header");

    double time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
    if (verbose){
        sprintf(MsgText,"CTE run time: %.2f(s) with %i procs/threads\n",time_spent/max_threads,max_threads);
        trlmessage(MsgText);
    }

    PrSwitch("pctecorr", COMPLETE);
    if (wf3.printtime)
        TimeStamp("PCTECORR Finished", wf3.rootname);

    freeAll(&ptrReg);
    return (status);
}


/********************* SUPPORTING SUBROUTINES *****************************/

int biasAndGainCorrect(SingleGroup *image, const float ccdGain, const Bool isSubarray, CTEParams * ctePars)
{
    extern int status;

    const unsigned nColumns = image->sci.data.nx;
    const unsigned nRows = image->sci.data.ny;
    const unsigned nColumnsPerChip = ctePars->nColumnsPerChip; /* for looping over quads  */
    const unsigned columnOffset = ctePars-> subarrayColumnOffset;

    float bias[4];
    float bsig[4];

    /*INIT THE ARRAYS*/
    for (unsigned i = 0; i < 4; ++i)
    {
        bias[i] = 0;
        bsig[i] = 0;
    }

    // Note that for user subarray the image is in only 1 quad, and only
    // has prescan bias pixels so the regions are different for full and subarrays
    //int (*findScanBias)(SingleGroup *, float *, float * ) = isSubarray ? &findPreScanBias : &findPostScanBias;

    // SUBTRACT THE EXTRA BIAS CALCULATED, AND MULTIPLY BY THE GAIN
    //(*findScanBias)(raz, bias, bsig);
    if (isSubarray)
        findOverScanBias(image, bias, bsig, PRESCAN);
    else
        findOverScanBias(image, bias, bsig, POSTSCAN);

#ifdef _OPENMP
    #pragma omp parallel shared(image, bias, bsig)
#endif
    {
    //Image will only ever be at most 2 amps i.e. a single chip
    //quads only split along a row
    unsigned quadStart[2];
    unsigned quadEnd[2];
    quadStart[0] = 0;
    quadEnd[0] = ctePars->nColumnsPerQuad - columnOffset;
    quadStart[1] = quadEnd[0] + 1;//+1???
    quadEnd[1] = nColumns - quadEnd[0];

    Bool isCDAmp = image->group_num == 1 ? True : False;
    unsigned ampIndexPrefix = isCDAmp ? 0 : 2;

    for (unsigned nthAmp = 0; nthAmp < 2; ++nthAmp)
    {
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned i = quadStart[nthAmp]; i < quadEnd[nthAmp]; ++i)
        {
            for (unsigned j = 0; j < nRows; ++j)
            {
                PixColumnMajor(image->sci.data, j, i*nRows) -= bias[ampIndexPrefix + nthAmp];
                PixColumnMajor(image->sci.data, j, i*nRows) *= ccdGain;
            }
        }
    }
    } // close parallel block

    return(status);
}

int findOverScanBias(SingleGroup *raz, float *mean, float *sigma, enum OverScanType overScanType)
{
	/*calculate the post scan and bias after the biac file has been subtracted
	  add some history information to the header

	  Jay gave no explanation why plist is limited to 55377 for full arrays, his
	  subarray limitation was just 1/4 of this value

	  the serial virtual overscan pixels are also called the trailing-edge pixels
	  these only exist in full frame images
	  */

	/*CALCULATE THE PRE SCAN AND BIAS AFTER THE BIAC FILE HAS BEEN SUBTRACTED

	  The serial physical overscan pixels are also known as the serial prescan,
	  they are the only pixels available for subarrays. For full frame arrays
	  the prescan is not used as part of the correction, instead the virtual
	  overscan pixels are used and modeled in findPostScanBias.

	*/
    /** this calls resistmean, which does a better job clipping outlying pixels
      that just a standard stddev clip single pass*/

    extern int status;
    const unsigned arraySize = 55377; // ????????
    const unsigned nRows = raz->sci.data.ny;
    const unsigned nColumns = raz->sci.data.nx;
    const unsigned nColumnsPerChip = nColumns / 4;

    // post-scan settings (full frame)
    unsigned iBegin = nRows + 5;
    unsigned iEnd = nColumnsPerChip - 1;
    // pre-scan settings (subarray)
    if (overScanType == PRESCAN)
    {
        iBegin = 5;
        iEnd = 25;
    }

    unsigned npix; /*track array size for resistant mean*/
    float min=0;
    float max=0;
    const float sigreg = 7.5; /*sigma clip*/
    float rmean;
    float rsigma;
    float * plist = calloc(arraySize, sizeof(*plist));    /*bias pixels to measure*/
    assert(plist);

    for (unsigned nthChip = 0; nthChip < 4; ++nthChip)
    {  /*for each quadrant, CDAB ordered*/
        npix = 0;
        rmean = 0;
        rsigma = 0;
        for (unsigned i = iBegin; i < iEnd; ++i)
        {
            for (unsigned j = 0; j < 2051; ++j) //why 2051 not 2070 (prob where overscan starts)
            { /*all rows*/
                if (npix < arraySize )
                {
                    if (PixColumnMajor(raz->dq.data, j, i+(nthChip*nColumnsPerChip))) // is this needed for full frame?
                    {
                        plist[npix] = PixColumnMajor(raz->sci.data, j, i+nthChip*nColumnsPerChip);
                        npix++;
                    }
                }
            }
         }

        if (npix > 0)
            resistmean(plist, npix, sigreg, &rmean, &rsigma, &min, &max);

        mean[nthChip] = rmean;
        sigma[nthChip] = rsigma;
        if (overScanType == PRESCAN && npix > 0)
            printf("npix=%i\nmean[%i]=%f\nsigma[%i] = %f\n",npix,nthChip+1,rmean,nthChip+1,rsigma);
    }

    if (plist)
        free(plist);

    return status;
}
int initCTETrl (char *input, char *output) {

    extern int status;

    char trl_in[SZ_LINE+1];     /* trailer filename for input */
    char trl_out[SZ_LINE+1];    /* output trailer filename */
    int exist;


    int MkName (char *, char *, char *, char *, char *, int);
    int TrlExists (char *);
    void SetTrlOverwriteMode (int);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;

    /* Input and output suffixes. */
    char *isuffix[] = {"_raw"};
    char *osuffix[] = {"_rac_tmp"};
    char *trlsuffix[] = {""};

    int nsuffix = 1;


    /* Start by stripping off suffix from input/output filenames */
    if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, SZ_LINE)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't determine trailer filename for %s",
                input);
        trlmessage (MsgText);
    }
    if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, SZ_LINE)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't create trailer filename for %s",
                output);
        trlmessage (MsgText);
    }

    /* NOW, CONVERT TRAILER FILENAME EXTENSIONS FROM '.FITS' TO '.TRL' */

    if (MkNewExtn (trl_in, TRL_EXTN) ) {
        sprintf (MsgText, "Error with input trailer filename %s", trl_in);
        trlerror (MsgText);
        WhichError (status);
    }
    if (MkNewExtn (trl_out, TRL_EXTN) ) {
        sprintf (MsgText, "Error with output trailer filename %s", trl_out);
        trlerror (MsgText);
        WhichError (status);
    }

    /* If we are working with a RAW file, then see if a TRL file
       needs to be overwritten after the generic conversion comments.  */
    if (strstr(input, isuffix[0]) != NULL) {
        /* Test whether the output file already exists */
        exist = TrlExists(trl_out);
        if (exist == EXISTS_YES) {
            /* The output file exists, so we want to add to them
             ** the new trailer comments.  */
            SetTrlOverwriteMode (NO);
        }
    }

    /* Sets up temp trailer file for output and copies input
     ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);

    return(status);
}
//NOTE: wf3 * ctePars should be const however this is too large a refactor for this PR
int getSubarray(SingleGroup * image, CTEParams * ctePars, WF3Info * wf3)
{
    extern int status;
    // get subarray from first extension
    getSingleGroup(wf3->input, 1, image);
    // getSingleGroup will set image->group_num to 1 here so now correct
    image->group_num = wf3->chip == 2 ? 1 : 2;
    if (hstio_err()){
        freeSingleGroup(image);
        return (status = OPEN_FAILED);
    }

    ctePars->nRows = image->sci.data.ny;
    ctePars->nColumns = image->sci.data.nx;

  // These are used to find subarrays with physical overscan
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

    ctePars->subarrayColumnOffset = sci_corner[0] - ref_corner[0];
    ctePars->subarrayRowOffset = sci_corner[1] - ref_corner[1];
    const unsigned start = ctePars->subarrayColumnOffset; //column where the subarray starts
    const unsigned finish = start + image->sci.data.nx; //column where the subarray ends

    if ( start >= 25 &&  finish + 60 <= (ctePars->nColumnsPerChip) - 25){
        sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",start,finish);
        trlmessage(MsgText);
        freeSingleGroup(image);
        return(ERROR_RETURN);
    }

    /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
      AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
      CTE CODE VARIABLES.
      */
    if (CompareCTEParams(image, ctePars))
    {
        freeSingleGroup(image);
        return (status);
    }

    /*Put the subarray data into full frame*/
    Bool virtual = True;
    //Sub2Full(wf3, subChip, fullChip, 0, 1, virtual);
    if (virtual && ctePars->subarrayColumnOffset >= 2072)
        ctePars->subarrayColumnOffset += 60; // image starts in B or D regions and we can just shift the starting pixel

    /* SAVE A COPY OF THE RAW IMAGE BEFORE BIAS FOR LATER */
   /* if (fullChip->group_num == 1)
        makeRAZ(fullChip, otherFullChip, backup);
    else
        makeRAZ(otherFullChip, fullChip, backup);
*/

    //save copy before doCteBias

    /* Subtract the BIAC file from the subarray before continuing
       The bias routine will take care of cutting out the correct
       image location for the subarray.*/
    /*if (doCteBias(wf3, subChip)){
        freeSingleGroup(subChip);
        return(status);
    }
    //set the array with bias subtracted data
    //Sub2Full(wf3, subChip, fullChip, 0, 1, 1);
*/
    return status;
}

int unalignAmps(SingleGroup * image, CTEParams * ctePars)
{
    //The function alignAmps is symmetric, wrapping as aid when reading logic
    return alignAmps(image, ctePars);
}

int alignAmps(SingleGroup * image, CTEParams * ctePars)
{
    //WARNING - assumes row major storage

    //Align amps such that they are at the bottom left
    extern int status;

    unsigned columnOffset = ctePars->subarrayColumnOffset;
    unsigned nColumns = ctePars->nColumns;
    unsigned nRows = ctePars->nRows;

    Bool isCDAmp = image->group_num == 1 ? True : False;

    if (isCDAmp)
        printf("subarray from CD amp\n");
    else
        printf("subarray from AB amp\n");


    //If subarray only in c quad do nothing as this is already bottom left aligned
    if (isCDAmp && columnOffset + nColumns < ctePars->nColumnsPerQuad)
        return status;

    //Find how much of subarray extends into either b or d quad and flip right to left
    if (columnOffset + nColumns > ctePars->nColumnsPerQuad)
    {
        printf("subarray extends into amps B or D\n");
        //grab a row, flip it, put it back
        unsigned rowLength = columnOffset + nColumns - ctePars->nColumnsPerQuad;
        unsigned quadBoundary = nColumns - rowLength;
        float * row = NULL;
        for (unsigned i = 0; i < nRows; ++i)
        {
            //find row
            row = image->sci.data.data + i*nColumns + quadBoundary;
            //flip right to left
            float tempPixel;
            for (unsigned j = 0; j < rowLength/2; ++j)
            {
                tempPixel = row[j];
                row[j] = row[rowLength-j];
                row[rowLength-j] = tempPixel;
            }
        }
    }

    //Only thing left is to flip ab chip upside down
    if (!isCDAmp) // isABAmp
    {
        //either physically align all or propagate throughout a mechanism to work on the array upside down (flip b quad though)
        //we'll just flip all for now
        float * tempRow = malloc(nColumns*sizeof(*tempRow));
        float * topRow = NULL;
        float * bottomRow = NULL;
        size_t rowSize = nColumns*sizeof(*tempRow);
        for (unsigned i = 0; i < nRows/2; ++i)
        {
            topRow = image->sci.data.data + i*nColumns;
            bottomRow = image->sci.data.data + (nRows-i-1)*nColumns;
            memcpy(tempRow, topRow, rowSize);
            memcpy(topRow, bottomRow, rowSize);
            memcpy(bottomRow, tempRow, rowSize);
        }
        free(tempRow);
    }
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

int getCCDChip(int * value, char * fileName, char * ename, int ever)
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

int outputImage(char * fileName, SingleGroup * image, CTEParams * ctePars)
{
    //Use for interim dev testing
    SingleGroup temp;
    initSingleGroup(&temp);
    if (image->sci.data.storageOrder == COLUMNMAJOR)
    {
        allocSingleGroup(&temp, image->sci.data.nx, image->sci.data.ny, False);
        copySingleGroup(&temp, image, ROWMAJOR);
        image = &temp;
    }

    unalignAmps(image, ctePars);

    int ret = putSingleGroup(fileName, image->group_num, image, 0);
    freeSingleGroup(&temp);
    return ret;
}
