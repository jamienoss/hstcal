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
    SingleGroup raz; /* THE LARGE FORMAT COMBINATION OF CDAB*/
    SingleGroup rsz; /* LARGE FORMAT READNOISE CORRECTED IMAGE */
    SingleGroup rsc; /* CTE CORRECTED*/
    SingleGroup rzc; /* FINAL CTE CORRECTED IMAGE */
    SingleGroup chg; /* THE CHANGE DUE TO CTE  */
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
    initSingleGroup(&raz);
    allocSingleGroup(&raz, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rsz);
    allocSingleGroup(&rsz, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rsc);
    allocSingleGroup(&rsc, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rzc);
    allocSingleGroup(&rzc, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&raw);
    allocSingleGroup(&raw, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&chg);
    allocSingleGroup(&chg, RAZ_COLS, RAZ_ROWS);

    /*hardset the science arrays*/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0;j<RAZ_ROWS;j++){
            Pix(raw.sci.data,i,j)=hardset;
            Pix(raz.sci.data,i,j)=hardset;
            Pix(rsz.sci.data,i,j)=hardset;
            Pix(rsc.sci.data,i,j)=hardset;
            Pix(rzc.sci.data,i,j)=hardset;
            Pix(chg.sci.data,i,j)=hardset;
        }
    }

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
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS);
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
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS);
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
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS);
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
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS);
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
        for(i=0; i< ab.dq.data.nx; i++){
            for(j=0; j< ab.dq.data.ny; j++){
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


    /*CONVERT TO RAZ, SUBTRACT BIAS AND CORRECT FOR GAIN*/
    if (raw2raz(&wf3, &cd, &ab, &raz))
        return (status);

    /***CALCULATE THE SMOOTH READNOISE IMAGE***/
    trlmessage("CTE: Calculating smooth readnoise image");


    /***CREATE THE NOISE MITIGATION MODEL ***/
    if (cte_pars.noise_mit == 0) {
        if (raz2rsz(&wf3, &raz, &rsz, cte_pars.rn_amp, max_threads))
            return (status);
    } else {
        trlmessage("Only noise model 0 implemented!");
        return (status=ERROR_RETURN);
    }

    //CONVERT THE READNOISE SMOOTHED IMAGE TO RSC IMAGE
    SingleGroup trapPixelMap;
    initSingleGroup(&trapPixelMap);
    allocSingleGroup(&trapPixelMap, RAZ_COLS, RAZ_ROWS);

    if (populateTrapPixelMap(&trapPixelMap, &cte_pars))
        return status;

    /*THIS IS RAZ2RAC_PAR IN JAYS CODE - MAIN CORRECTION LOOP IN HERE*/
    if (inverseCTEBlur(&rsz, &rsc, &trapPixelMap, &cte_pars, wf3.verbose, wf3.expstart))
        return status;

    freeSingleGroup(&trapPixelMap);

    /*** CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT***/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0; j<RAZ_ROWS; j++){
           Pix(chg.sci.data,i,j) = (Pix(rsc.sci.data,i,j) - Pix(rsz.sci.data,i,j))/wf3.ccdgain;
           Pix(rzc.sci.data,i,j) =  Pix(raw.sci.data,i,j) + Pix(chg.sci.data,i,j);
        }
    }

    /*BACK TO NORMAL FORMATTING*/
    /*Copies rzc data to cd->sci.data and ab->sci.data */
    undoRAZ(&cd,&ab,&rzc);

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
            PutKeyDbl(subcd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");

            Full2Sub(&wf3, &subcd, &cd, 0, 1, 1);
            putSingleGroup(output, 1, &subcd,0);
            freeSingleGroup(&subcd);
        } else {

            /*** SAVE USEFUL HEADER INFORMATION ***/
            if (cteHistory (&wf3, subab.globalhdr))
                return (status);

            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(subab.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
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
        PutKeyDbl(cd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
        trlmessage("PCTEFRAC saved to header");

        putSingleGroup(output,cd.group_num, &cd,0);
        putSingleGroup(output,ab.group_num, &ab,0);
    }

    /** CLEAN UP ON AISLE 3 **/
    freeSingleGroup(&rzc);
    freeSingleGroup(&rsc);
    freeSingleGroup(&chg);
    freeSingleGroup(&raz);
    freeSingleGroup(&rsz);
    freeSingleGroup(&raw);
    freeSingleGroup(&cd);
    freeSingleGroup(&ab);

    freeCTEParams(&cte_pars);

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

int raw2raz(WF3Info *wf3, SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){
    /*

       convert a raw file to raz file: CDAB longwise amps, save data array
       for comparison with what jay has during testing

       -->do an additional bias correction using the  residual bias level measured for each amplifier from the
       steadiest pixels in the horizontal overscan and subtracted fom the pixels for that amplifier.

       ---> convert into electrons at the end
       ---> add supplemental bias info to the header

       allocate contiguous 2d array on the heap
       with pointers and return the pointer to the head of the array

       The Following macros are used to represent 2-d indexing.
       Two dimensional arrays are stored in FITS order.

       ny
       ^
       N | a05   a15   a25   a35
       A | a04   a14   a24   a34
       X | a03   a13   a23   a33
       I | a02   a12   a22   a32
       S | a01   a11   a21   a31
       2 | a00   a10   a20   a30
       ---------------------------> nx
       NAXIS1

       NAXIS1 is 4 and NAXIS2 is 6
       PIX(a,1,4) accesses a14

       In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left.
       The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan),
       with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down,
       then serial-shifted to the left.

*/
    extern int status;

    int i,j,k;              /*loop counters*/
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    extern int status;      /* variable for return status */
    float bias_post[4];
    float bsig_post[4];
    float bias_pre[4];
    float bsig_pre[4];
    float gain;

    /*INIT THE ARRAYS*/
    for(i=0;i<4;i++){
        bias_post[i]=0.;
        bsig_post[i]=0.;
        bias_pre[i]=0.;
        bsig_pre[i]=0.;
    }

    gain=wf3->ccdgain;

    /*REFORMAT TO RAZ*/
    makeRAZ(cd,ab,raz);


    /*SUBTRACT THE EXTRA BIAS CALCULATED, AND MULTIPLY BY THE GAIN
      Note that for user subarray the image is in only 1 quad, and only
      has prescan bias pixels so the regions are different for full and subarrays
    */
    if (wf3->subarray){
        findPreScanBias(raz, bias_pre, bsig_pre);
        for (k=0;k<4;k++){
            for (i=0; i<subcol;i++){
                for (j=0;j<RAZ_ROWS; j++){
                    if(Pix(raz->dq.data,i+k*subcol,j))
                        Pix(raz->sci.data,i+k*subcol,j) -= bias_pre[k];
                        Pix(raz->sci.data,i+k*subcol,j) *= gain;
                }
            }
        }
    } else {
        findPostScanBias(raz, bias_post, bsig_post);
        for (k=0;k<4;k++){
            for (i=0; i<subcol;i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->sci.data,i+k*subcol,j) -= bias_post[k];
                    Pix(raz->sci.data,i+k*subcol,j) *= gain;
                }
            }
        }
    }

    return(status);
}

/*calculate the post scan and bias after the biac file has been subtracted
  add some history information to the header

  Jay gave no explanation why plist is limited to 55377 for full arrays, his
  subarray limitation was just 1/4 of this value

  the serial virtual overscan pixels are also called the trailing-edge pixels
  these only exist in full frame images
  */

int findPostScanBias(SingleGroup *raz, float *mean, float *sigma){

    extern int status;
    int arrsize = 55377;
    int i,j,k;              /*Looping variables */
    float plist[arrsize];  /*bias bpixels to measure*/
    float *plistSub;
    float min=0.0;
    float max=0.0;
    float rmean=0.0;
    float rsigma=0.0;
    float sigreg =7.5; /*sigma clip*/


    int subcol = RAZ_COLS/4;
    int npix=0; /*track array size for resistant mean*/

    /*init plist for full size
      We'll allocate heap memory for smaller arrays
      */
    for (i=0;i<arrsize;i++){
        plist[i]=0.;
    }

    for (k=0;k<4;k++){  /*for each quadrant cdab = 0123*/
        npix=0; /*reset for each quad*/
        rmean=0.;
        rsigma=0.;
        for (i=RAZ_ROWS+5;i<= subcol-1; i++){ /*quad area for post scan bias pixels*/
            for (j=0; j<2051; j++){
                if (npix < arrsize){
                    if ( Pix(raz->dq.data,i+k*subcol,j)) {
                        plist[npix] = Pix(raz->sci.data,i+k*subcol,j);
                        npix+=1;
                    }
                }
            }
        }
        if (npix > 0 ){
            plistSub = (float *) calloc(npix, sizeof(float));
            if (plistSub == NULL){
                trlerror("out of memory for resistmean entrance in findPostScanBias.");
                free(plistSub);
                return (ERROR_RETURN);
            }
            for(i=0; i<npix; i++){
                plistSub[i]=plist[i];
            }
            resistmean(plistSub, npix, sigreg, &rmean, &rsigma, &min, &max);
            free(plistSub);
        }
        mean[k]= rmean;
        sigma[k] = rsigma;
    }
    return status;
}

/*CALCULATE THE PRE SCAN AND BIAS AFTER THE BIAC FILE HAS BEEN SUBTRACTED

  The serial physical overscan pixels are also known as the serial prescan,
  they are the only pixels available for subarrays. For full frame arrays
  the prescan is not used as part of the correction, instead the virtual
  overscan pixels are used and modeled in findPostScanBias.

*/

int findPreScanBias(SingleGroup *raz, float *mean, float *sigma){
    /** this calls resistmean, which does a better job clipping outlying pixels
      that just a standard stddev clip single pass*/

    extern int status;
    int arrsize = 55377;
    int i,j,k;              /*Looping variables */
    float plist[arrsize];    /*bias pixels to measure*/
    float *plistSub; /*heap allocation for variable size plist array*/
    float min=0.0;
    float max=0.0;
    float rmean;
    float rsigma;
    float sigreg =7.5; /*sigma clip*/
    int subcol = RAZ_COLS/4;
    int npix=0; /*track array size for resistant mean*/


    /*init plist*/
    for (i=0;i<arrsize;i++){
        plist[i]=0.;
    }

    for (k=0;k<4;k++){  /*for each quadrant, CDAB ordered*/
        npix=0;
        rmean=0.;
        rsigma=0.;
        for (i=5;i<25; i++){
            for (j=0; j<2051; j++){ /*all rows*/
                if (npix < arrsize ){
                    if (Pix(raz->dq.data,i+(k*subcol),j)){
                        plist[npix] = Pix(raz->sci.data,i+k*subcol,j);
                        npix+=1;
                    }
                }
            }
         }

        if (0 < npix ){
            plistSub = (float *) calloc(npix, sizeof(float));
            if (plistSub == NULL){
                trlerror("out of memory for resistmean entrance in findPostScanBias.");
                free(plistSub);
                return (ERROR_RETURN);
            }
            for(i=0; i<npix; i++){
                plistSub[i]=plist[i];
            }
            resistmean(plistSub, npix, sigreg, &rmean, &rsigma, &min, &max);
            free(plistSub);
        }

        mean[k]= rmean;
        sigma[k] = rsigma;
        if(npix>0)
            printf("npix=%i\nmean[%i]=%f\nsigma[%i] = %f\n",npix,k+1,rmean,k+1,rsigma);
    }
    return status;
}


int raz2rsz(WF3Info *wf3, SingleGroup *raz, SingleGroup *rsz, double rnsig, int max_threads){
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

    clock_t begin = clock();

    /*ALL ELEMENTS TO FLAG*/
    for(unsigned i = 0; i < 3; ++i){
        for (unsigned j = 0; j < RAZ_ROWS; ++j){
            obs_loc[i][j]=setdbl;
            rsz_loc[i][j]=setdbl;
        }
    }

    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, RAZ_COLS, RAZ_ROWS);

    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, RAZ_COLS, RAZ_ROWS);


    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    for(unsigned i = 0; i < RAZ_COLS; ++i){
        for (unsigned j = 0; j < RAZ_ROWS; ++j){
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

    for(unsigned NIT=1; NIT<=100; ++NIT){
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) \
        private(imid,obs_loc,rsz_loc,dptr)\
        shared(raz, rsz, rnsig,rms,nrms, zadj)
#endif
        for(unsigned i = 0; i < RAZ_COLS; ++i){
            imid=i;
            /*RESET TO MIDDLE RAZ_COLS AT ENDPOINTS*/
            if (imid < 1)
                imid=1;
            if (imid == RAZ_COLS-1)
                imid = RAZ_COLS-2;

            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            for(unsigned j = 0; j < RAZ_ROWS; ++j){
                obs_loc[0][j] = Pix(raz->sci.data,imid-1,j);
                obs_loc[1][j] = Pix(raz->sci.data,imid,j);
                obs_loc[2][j] = Pix(raz->sci.data,imid+1,j);

                rsz_loc[0][j] = Pix(rsz->sci.data,imid-1,j);
                rsz_loc[1][j] = Pix(rsz->sci.data,imid,j);
                rsz_loc[2][j] = Pix(rsz->sci.data,imid+1,j);
            }
            for (unsigned j = 0; j < RAZ_ROWS; ++j){
             if(Pix(raz->dq.data,imid,j)) {
                find_dadj(1+i-imid,j, obs_loc, rsz_loc, rnsig, &dptr);
                Pix(zadj.sci.data,i,j) = dptr;
              }
            }
        } /*end the parallel for*/

        /*NOW GO OVER ALL THE RAZ_COLS AND RAZ_ROWS AGAIN TO SCALE THE PIXELS
        */
        for(unsigned i = 0; i < RAZ_COLS; ++i){
            for(unsigned j = 0; j < RAZ_ROWS; ++j){
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
        private(rmsu,nrmsu) \
        shared(raz,rsz,rms,rnsig,nrms)
        for(unsigned j = 0; j < RAZ_ROWS; ++j){
            nrmsu=setdbl;
            rmsu=setdbl;
            for(unsigned i = 0; i < RAZ_COLS; ++i){
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

    if (wf3->verbose)
    {
    	double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
    	sprintf(MsgText,"Time taken to smooth image: %.2f(s) with %i procs/threads\n",timeSpent/max_threads,max_threads);
    	trlmessage(MsgText);
    }

    return (status);
}


int find_dadj(int i ,int j, double obsloc[][RAZ_ROWS], double rszloc[][RAZ_ROWS], double rnsig, double *d){
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

    *d = ((dval0u * w0 * 0.25f) + /* desire to keep the original pixel value */
            (dval9u*w9*0.25f) + /* desire to keep the original sum over 3x3*/
            (dmod1u*w1*0.25f) + /*desire to get closer to the pixel below*/
            (dmod2u*w2*0.25f)) ; /*desire to get closer to the pixel above*/

    return(status);
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
