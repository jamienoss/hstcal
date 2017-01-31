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
    IODescPtr ip = NULL;

    SingleGroup cd; /*SCI 1, chip 2*/
    SingleGroup ab; /*SCI 2, chip 1*/
    SingleGroup subcd; /*subarray chip*/
    SingleGroup subab; /*subarray chip*/
    SingleGroup raw; /* THE RAW IMAGE IN RAZ FORMAT */

    int i,j; /*loop vars*/
    int max_threads=1;
    clock_t begin;
    double  time_spent;
    float hardset=0.0;

    /* These are used to find subarrays with physical overscan */
    int sci_bin[2];			/* bin size of science image */
    int sci_corner[2];		/* science image corner location */
    int ref_bin[2];
    int ref_corner[2];
    int rsize = 1;          /* reference pixel size */
    int start=0;            /*where the subarray starts*/
    int finish=0;           /*where the subarray ends*/

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
    begin = (double)clock();

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
        return(ERROR_RETURN);
    }
    if (wf3.darkcorr == COMPLETE){
        trlmessage("DARKCORR complete for input image, CTE can't be performed");
        return(ERROR_RETURN);
    }
    if (wf3.blevcorr == COMPLETE){
        trlmessage("BLEVCORR complete for input image, CTE can't be performed");
        return(ERROR_RETURN);
    }

    /* DETERMINE THE NAMES OF THE TRAILER FILES BASED ON THE INPUT
       AND OUTPUT FILE NAMES, THEN INITIALIZE THE TRAILER FILE BUFFER
       WITH THOSE NAMES.
       */
    if (initCTETrl (input, output))
        return (status);

    /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
    Hdr phdr; /*primary header for input image, all output information saved here*/
    initHdr(&phdr);
    if (LoadHdr (wf3.input, &phdr) ){
        WhichError(status);
        freeHdr (&phdr);
        return (ERROR_RETURN);
    }

    /* GET KEYWORD VALUES FROM PRIMARY HEADER. */
    if (GetKeys (&wf3, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }

    if (GetCTEFlags (&wf3, &phdr)) {
        freeHdr(&phdr);
        return (status);
    }


    /*SET UP THE ARRAYS WHICH WILL BE PASSED AROUND*/
    initSingleGroup(&raw);
    allocSingleGroup(&raw, RAZ_COLS, RAZ_ROWS, False);

    /*READ IN THE CTE PARAMETER TABLE*/
    CTEParams cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    initCTEParams(&cte_pars, TRAPS, RAZ_ROWS, RAZ_COLS);
    allocateCTEParams(&cte_pars);
    if (GetCTEPars (wf3.pctetab.name, &cte_pars))
        return (status);

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);
    }
    /* Full frame and subarrays always have group 1
       If it's a subarray, the group can be from either chip
       and will still be labled group 1 because it's the FIRST
       and only group, so look at the ccdchip instead.

       amps ab are in chip1, sci,2
       amps cd are in chip2, sci,1

    */
    Hdr scihdr; /*science header in case of subarray image to detect chip*/
    initHdr(&scihdr);
    if (wf3.subarray) {
        /* OPEN INPUT IMAGE IN ORDER TO READ ITS SCIENCE HEADER. */
        ip = openInputImage (wf3.input, "SCI", 1);
        if (hstio_err()) {
            sprintf (MsgText, "Image: \"%s\" is not present", wf3.input);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }
        getHeader (ip, &scihdr);
        if (ip != NULL)
            closeImage (ip);

        /* Get CCD-specific parameters. */
        if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, 1, &wf3.chip)){
            freeHdr(&scihdr);
            return (status);
        }
        freeHdr(&scihdr);

        if (wf3.chip == 2){ /*sci1,cd*/
            start=0;
            finish=0;
            /*get CD subarray from first extension*/
            initSingleGroup (&subcd);
            getSingleGroup (wf3.input, 1, &subcd);
            if (hstio_err()){
                freeSingleGroup(&subcd);
                return (status = OPEN_FAILED);
            }

            /*create an empty full size chip for pasting*/
            initSingleGroup(&cd);
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS, False); // CreateEmptyChip 0 init
            cd.group_num=1;
            CreateEmptyChip(&wf3, &cd);

            if (GetCorner(&subcd.sci.hdr, rsize, sci_bin, sci_corner))
                return (status);
            if (GetCorner(&cd.sci.hdr, rsize, ref_bin, ref_corner))
                return (status);

            start = sci_corner[0] - ref_corner[0];
            finish = start + subcd.sci.data.nx;
            if ( start >= 25 &&  finish + 60 <= (RAZ_COLS/2) - 25){
                sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",start,finish);
                trlmessage(MsgText);
                return(ERROR_RETURN);
            }

            /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
              AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
              CTE CODE VARIABLES.
              */
            if (CompareCTEParams(&subcd, &cte_pars))
                return (status);

            /*Put the subarray data into full frame*/
            Sub2Full(&wf3, &subcd, &cd, 0, 1, 1);

            /* now create an empty chip 1*/
            initSingleGroup(&ab);
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS, False); // CreateEmptyChip 0 init
            ab.group_num=2;
            CreateEmptyChip(&wf3, &ab);

            /* SAVE A COPY OF THE RAW IMAGE BEFORE BIAS FOR LATER */
            makeRAZ(&cd,&ab,&raw);

            /* Subtract the BIAC file from the subarray before continuing
               The bias routine will take care of cutting out the correct
               image location for the subarray.*/

            if (doCteBias(&wf3,&subcd)){
                freeSingleGroup(&subcd);
                return(status);
            }

            /*reset the array after bias subtraction*/
            Sub2Full(&wf3, &subcd, &cd, 0, 1, 1);


        } else { /*chip is 1, ab, sci2*/
            start=0;
            finish=0;
            initSingleGroup(&subab);
            getSingleGroup(wf3.input, 1, &subab);
            if (hstio_err()){
                freeSingleGroup(&subab);
                return (status = OPEN_FAILED);
            }

            /*make an empty fullsize chip for pasting*/ //Redo only allocate enough for subarray and not full frame!
            initSingleGroup(&ab);
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS, False); // Sub2Full will 0 init
            ab.group_num=2;
            CreateEmptyChip(&wf3, &ab);

            if ( GetCorner(&subab.sci.hdr, rsize, sci_bin, sci_corner))
                return (status);

            if ( GetCorner(&ab.sci.hdr, rsize, ref_bin, ref_corner))
                return (status);

            start = sci_corner[0] - ref_corner[0];
            finish = start + subab.sci.data.nx + 2103;
            finish = start + subab.sci.data.nx;
            if ( start >= 25 &&  finish + 60 <= (RAZ_COLS/2) - 25){
                sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",start,finish);
                trlmessage(MsgText);
                return(ERROR_RETURN);
            }
            /*add subarray to full frame image*/
            Sub2Full(&wf3, &subab, &ab, 0, 1, 1);

            /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
              AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
              CTE CODE VARIABLES.
              */
            if (CompareCTEParams(&subab, &cte_pars))
                return (status);

            /* now create an empty chip 2*/
            initSingleGroup(&cd);
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS, False); // CreateEmptyChip 0 init
            cd.group_num=1;
            CreateEmptyChip(&wf3, &cd);

            /* SAVE A COPY OF THE RAW IMAGE FOR LATER */
            makeRAZ(&cd,&ab,&raw);

            /* Subtract the BIAC file from the subarray before continuing*/
            subab.group_num=2;
            if (doCteBias(&wf3,&subab)){
                freeSingleGroup(&subab);
                return(status);
            }

            /*reset the array after bias subtraction*/
            Sub2Full(&wf3, &subab, &ab, 0, 1, 1);
        }

    } else {
        /* Full frame image, just read in the groups
           and init the mask to use all pixels
        */

        initSingleGroup (&cd);
        getSingleGroup (wf3.input, 1, &cd);
        if (hstio_err()){
            return (status = OPEN_FAILED);
        }

        initSingleGroup (&ab);
        getSingleGroup (wf3.input, 2, &ab);
        if (hstio_err()){
            return (status = OPEN_FAILED);
        }

        /*setup the mask*/
        for (unsigned j = 0; j < ab.dq.data.ny; ++j)
        {
            for (unsigned i = 0; i < ab.dq.data.nx; ++i)
            {
                PPix(&ab.dq.data, i, j) = 1;
                PPix(&cd.dq.data, i, j) = 1;
            }
        }

        /* SAVE A COPY OF THE RAW IMAGE FOR LATER */
        makeRAZ(&cd,&ab,&raw);

        /***SUBTRACT THE CTE BIAS FROM BOTH CHIPS IN PLACE***/
        if (doCteBias(&wf3,&cd)){
            freeSingleGroup(&cd);
            return(status);
        }

        if (doCteBias(&wf3,&ab)){
            freeSingleGroup(&ab);
            return(status);
        }
        /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
          AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
          CTE CODE VARIABLES.
          */
        if (CompareCTEParams(&cd, &cte_pars))
            return (status);

    }

    //CONVERT TO RAZ
    SingleGroup raz; /* THE LARGE FORMAT COMBINATION OF CDAB*/
    initSingleGroup(&raz);
    allocSingleGroup(&raz, RAZ_COLS, RAZ_ROWS, False);
    makeRAZ(&cd, &ab, &raz);

    //copy to column major storage
    SingleGroup tempGroup1;
    initSingleGroup(&tempGroup1);
    allocSingleGroup(&tempGroup1, RAZ_COLS, RAZ_ROWS, False);
    SingleGroup * razColumnMajor = &tempGroup1;
    assert(!copySingleGroup(razColumnMajor, &raz, COLUMNMAJOR));

    //SUBTRACT BIAS AND CORRECT FOR GAIN
    if (biasAndGainCorrect(razColumnMajor, wf3.ccdgain, wf3.subarray))
        return (status);

    /***CALCULATE THE SMOOTH READNOISE IMAGE***/
    trlmessage("CTE: Calculating smooth readnoise image");

    /* LARGE FORMAT READNOISE CORRECTED IMAGE */
    SingleGroup tempGroup2;
    initSingleGroup(&tempGroup2);
    allocSingleGroup(&tempGroup2, RAZ_COLS, RAZ_ROWS, False);
    SingleGroup * smoothedImage = &tempGroup2;

    /***CREATE THE NOISE MITIGATION MODEL ***/
    if (cte_pars.noise_mit == 0) {
        if (cteSmoothImage(razColumnMajor, smoothedImage, cte_pars.rn_amp, max_threads, wf3.verbose))
            return (status);
    } else {
        trlmessage("Only noise model 0 implemented!");
        return (status=ERROR_RETURN);
    }

    razColumnMajor = NULL;
    SingleGroup * trapPixelMap = &tempGroup1;
    if (populateTrapPixelMap(trapPixelMap, &cte_pars))
    {
    	//Geez, do we not want to free everything???
        return status;
    }

    //reuse raz
    SingleGroup * cteCorrectedImage = &raz;
    /*THIS IS RAZ2RAC_PAR IN JAYS CODE - MAIN CORRECTION LOOP IN HERE*/
    if (inverse_cte_blur(smoothedImage, cteCorrectedImage, trapPixelMap, &cte_pars, wf3.verbose, wf3.expstart))
        return status;
    trapPixelMap = NULL;
    freeSingleGroup(&tempGroup1);

    const double scaleFraction = cte_pars.scale_frac;
    freeCTEParams(&cte_pars);

    /*** CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT***/
    const float ccdgain = wf3.ccdgain;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static), shared(cteCorrectedImage, smoothedImage, raw)
#endif
    for (unsigned i = 0; i < RAZ_COLS; ++i)
    {
        for(unsigned j = 0; j < RAZ_ROWS; ++j)
            Pix(raw.sci.data,i,j) += (PixColumnMajor(cteCorrectedImage->sci.data,j,i) - PixColumnMajor(smoothedImage->sci.data,j,i))/ccdgain;
            //Pix(raw.sci.data,i,j) += (Pix(cteCorrectedImage->sci.data,i,j) - PixColumnMajor(smoothedImage->sci.data,j,i))/ccdgain;
    }
    cteCorrectedImage = NULL;
    freeSingleGroup(&raz);
    smoothedImage = NULL;
    freeSingleGroup(&tempGroup2);


    /*BACK TO NORMAL FORMATTING*/
    /*Copies rzc data to cd->sci.data and ab->sci.data */
    undoRAZ(&cd, &ab, &raw);
    freeSingleGroup(&raw);

    /* COPY BACK THE SCIENCE SUBARRAYS AND
       SAVE THE NEW RAW FILE WITH UPDATED SCIENCE
       ARRAYS AND PRIMARY HEADER TO RAC
       */
    if (wf3.subarray) {
        if (wf3.chip == 2) {
            /*** SAVE USEFUL HEADER INFORMATION ***/
            if (cteHistory (&wf3, subcd.globalhdr))
                return (status);

            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(subcd.globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");

            Full2Sub(&wf3, &subcd, &cd, 0, 1, 1);
            putSingleGroup(output, 1, &subcd,0);
            freeSingleGroup(&subcd);
        } else {

            /*** SAVE USEFUL HEADER INFORMATION ***/
            if (cteHistory (&wf3, subab.globalhdr))
                return (status);

            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(subab.globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");

            Full2Sub(&wf3, &subab, &ab, 0, 1, 1);
            putSingleGroup(output, 1, &subab,0);
            freeSingleGroup(&subab);
        }

    } else { /*FUll FRAME*/
        /*** SAVE USEFUL HEADER INFORMATION ***/
        if (cteHistory (&wf3, cd.globalhdr))
            return (status);

        /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
        PutKeyDbl(cd.globalhdr, "PCTEFRAC", scaleFraction,"CTE scaling fraction based on expstart");
        trlmessage("PCTEFRAC saved to header");

        putSingleGroup(output,cd.group_num, &cd,0);
        putSingleGroup(output,ab.group_num, &ab,0);
    }

    /** CLEAN UP ON AISLE 3 **/
    freeSingleGroup(&cd);
    freeSingleGroup(&ab);

    time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
    if (verbose){
        sprintf(MsgText,"CTE run time: %.2f(s) with %i procs/threads\n",time_spent/max_threads,max_threads);
        trlmessage(MsgText);
    }

    PrSwitch("pctecorr", COMPLETE);
    if(wf3.printtime)
        TimeStamp("PCTECORR Finished",wf3.rootname);

    return (status);
}


/********************* SUPPORTING SUBROUTINES *****************************/

int biasAndGainCorrect(SingleGroup *raz, const float ccdGain, const Bool isSubarray)
{
    extern int status;

    const unsigned nColumns = raz->sci.data.nx;
    const unsigned nRows = raz->sci.data.ny;
    const unsigned nColumnsPerChip = nColumns / 4; /* for looping over quads  */

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
        findOverScanBias(raz, bias, bsig, PRESCAN);
    else
        findOverScanBias(raz, bias, bsig, POSTSCAN);

#ifdef _OPENMP
    #pragma omp parallel shared(raz, bias, bsig)
#endif
    {
    for (unsigned nthChip = 0; nthChip < 4; ++nthChip)
    {
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned i = 0; i < nColumnsPerChip; ++i)
        {
            for (unsigned j = 0; j < nRows; ++j)
            {
                if (isSubarray)
                {
                    if(PixColumnMajor(raz->dq.data, j, i+nthChip*nColumnsPerChip))
                    //{
                        PixColumnMajor(raz->sci.data, j, i+nthChip*nColumnsPerChip) -= bias[nthChip];
                        PixColumnMajor(raz->sci.data, j, i+nthChip*nColumnsPerChip) *= ccdGain;
                    //}
                }
                else
                {
                    PixColumnMajor(raz->sci.data, j, i+nthChip*nColumnsPerChip) -= bias[nthChip];
                    PixColumnMajor(raz->sci.data, j, i+nthChip*nColumnsPerChip) *= ccdGain;
                }
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
