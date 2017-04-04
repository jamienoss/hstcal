#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "acserr.h"

#include "pcte.h"

# ifdef _OPENMP
#  include <omp.h>
# endif
# include "../../../../ctegen2/ctegen2.h"
#include <assert.h>

static int make_amp_array(const ACSInfo *acs, const SingleGroup *input,
                          const int amp,
                          const int xbeg, const int ybeg,
                          SingleGroup * output
                          );

static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *input,
                            const int amp,
                            const int xbeg, const int ybeg,
                            SingleGroup * output
                            );

int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *image,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);


int doPCTEGen2 (ACSInfo *acs, SingleGroup * image)
{

    /* arguments:
       acs     i: calibration switches, etc
       x      io: image to be calibrated; written to in-place
    */

    extern int status;

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    int max_threads=1;
#   ifdef _OPENMP
    trlmessage("Using parallel processing provided by OpenMP inside CTE routine");
    if (acs->onecpu){
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

    char ccdamp[strlen(AMPSTR1)+1]; /* string to hold amps on current chip */
    int numamps;               /* number of amps on chip */
    int amp;                   /* index amp A:0, B:1, etc. */
    char * amploc;             /* pointer to amp character in AMPSORDER */
    int nRows, nColumns;    /* int for C array dimension sizes */
    int amp_xsize, amp_ysize;  /* int for amp array size (x/y in CCD coords) */
    int amp_xbeg, amp_xend;    /* int for beg and end of amp arrays on chip */
    int amp_ybeg, amp_yend;

    /* functions from calacs/lib */
    void parseWFCamps (char *acsamps, int chip, char *ccdamp);
    void TimeStamp (char *message, char *rootname);
    int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);

    /* test whether this we've been given ACS WFC data, since the CTE algorithm
       is currently only valid for WFC data. */
    if (acs->detector != WFC_CCD_DETECTOR) {
        trlerror("(pctecorr) only valid for WFC CCD data, PCTECORR should be OMIT for all others.");
        return (status = ERROR_RETURN);
    }

    if (acs->printtime) {
        TimeStamp("Starting CTE correction...","");
    }

    /**************** read and calculate parameters of CTE model ************/
    /* structure to hold CTE parameters from file */
    ACSCTEParams pars;
    if (PixCteParams(acs->pcte.name, acs->expstart, &pars)) {
        return (status);
    }

    if (CompareCteParams(image, &pars)) {
        return (status);
    }

    sprintf(MsgText, "(pctecorr) Read noise level PCTERNCL: %f", pars.rn_clip);
    trlmessage(MsgText);
    sprintf(MsgText, "(pctecorr) Readout simulation iterations PCTESMIT: %i",
            pars.sim_nit);
    trlmessage(MsgText);
    sprintf(MsgText, "(pctecorr) Number of readout shifts PCTESHFT: %i",
            pars.shft_nit);
    trlmessage(MsgText);
    sprintf(MsgText, "(pctecorr) CTE_FRAC: %f", pars.cte_frac);
    trlmessage(MsgText);

    /* also add cte_frac as header keyword */
    if ((acs->chip == 2) || (acs->subarray == YES)) {
        if (PutKeyDbl(x->globalhdr, "PCTEFRAC", pars.cte_frac,
                      "CTE scaling factor")) {
            trlerror("(pctecorr) Error writing PCTEFRAC to image header");
            return (status = HEADER_PROBLEM);
        }
    }

 /*   if (InterpolatePsi(pars.chg_leak, pars.psi_node, chg_leak, chg_open)) {
        return (status);
    }
    if (InterpolatePhi(pars.dtde_l, pars.q_dtde, pars.shft_nit, dtde_q)) {
        return (status);
    }
    if (FillLevelArrays(chg_leak, chg_open, dtde_q, pars.levels, chg_leak_lt,
                        chg_open_lt, dpde_l)) {
        return (status);
    }
    */
    /********* done reading and calculating CTE model parameters ************/

    /* need to figure out which amps are on this chip */
    ccdamp[0] = '\0'; /* "reset" the string for reuse */
    parseWFCamps(acs->ccdamp, acs->chip, ccdamp);

    //Store sizes - these are corrected for subarrays in getSubarray()

    /* loop over amps on this chip and do CTE correction */
    numamps = strlen(ccdamp);
    for (unsigned i = 0; i < numamps; ++i) {
        sprintf(MsgText, "(pctecorr) Performing CTE correction for amp %c",
                ccdamp[i]);
        trlmessage(MsgText);

        /* get the amp letter and number where A:0, B:1, etc. */
        amploc = strchr(AMPSORDER, ccdamp[i]);
        amp = *amploc - AMPSORDER[0];

        /* get amp array size */
        if (get_amp_array_size_acs_cte(acs, image, amp, amploc, ccdamp,
                               &amp_xsize, &amp_ysize, &amp_xbeg,
                               &amp_xend, &amp_ybeg, &amp_yend)) {
            freeAll(&ptrReg);
            return (status);
        }

        nRows = amp_ysize;
        nColumns = amp_xsize;

        initCTEParamsFast(&pars.baseParams, TRAPS, nRows, nColumns);
        addPtr(&ptrReg, &pars.baseParams, &freeCTEParamsFast);
        allocateCTEParamsFast(&pars.baseParams);
        if (GetCTEParsFast (acs->pcteTabNameFromCmd, &pars.baseParams))
        {
            freeAll(&ptrReg);
            return (status);
        }

        pars.baseParams.nRows = nRows;
        pars.baseParams.nColumns = nColumns;
        pars.baseParams.nRowsPerFullFrame = nRows;
        pars.baseParams.nColumnsPerFullFrame = nColumns;
        pars.baseParams.nColumnsPerChip = pars.baseParams.nColumnsPerFullFrame/2;
        pars.baseParams.nRowsPerChip = pars.baseParams.nRowsPerFullFrame;
        pars.baseParams.nColumnsPerQuad = pars.baseParams.nColumnsPerFullFrame/4;
        pars.baseParams.nRowsPerQuad = pars.baseParams.nRowsPerFullFrame;
        pars.baseParams.isSubarray = acs->subarray;
        pars.baseParams.refAndIamgeBinsIdenticle = True;
        pars.baseParams.columnOffset = 0;
        pars.baseParams.rowOffset = 0;

        //This is used for the final output
        SingleGroup raw;
        initSingleGroup(&raw);
        addPtr(&ptrReg, &raw, &freeSingleGroup);
        allocSingleGroup/*SciOnly*/(&raw, nColumns, nRows, False);

        /* read data from the SingleGroup into an array containing data from
           just one amp */
        if (make_amp_array(acs, x, amp, amp_xbeg, amp_ybeg, &raw))
        {
            freeAll(&ptrReg);
            return (status);
        }

        //copy to column major storage
        SingleGroup columnMajorImage;
        initSingleGroup(&columnMajorImage);
        addPtr(&ptrReg, &columnMajorImage, &freeSingleGroup);
        allocSingleGroup/*SciOnly*/(&columnMajorImage, nColumns, nRows, False);
        assert(!copySingleGroup(&columnMajorImage, &raw, COLUMNMAJOR));

        /***CALCULATE THE SMOOTH READNOISE IMAGE***/
        trlmessage("CTE: Calculating smooth readnoise image");
        SingleGroup smoothedImage;
        initSingleGroup(&smoothedImage);
        addPtr(&ptrReg, &smoothedImage, &freeSingleGroup);
        allocSingleGroup/*SciOnly*/(&smoothedImage, nColumns, nRows, False);
        setStorageOrder(&smoothedImage, COLUMNMAJOR);

        /* do some smoothing on the data so we don't amplify the read noise.
           data should be in electrons. */
        if (cteSmoothImage(&columnMajorImage, &smoothedImage, &pars.baseParams, pars.baseParams.rn_amp /*pars.rn_clip*/, max_threads, acs->verbose))
        {
            freeAll(&ptrReg);
            return (status);
       }

       SingleGroup trapPixelMap;
       initSingleGroup(&trapPixelMap);
       addPtr(&ptrReg, &trapPixelMap, &freeSingleGroup);
       allocSingleGroup/*SciOnly*/(&trapPixelMap, nColumns, nRows, False);
       setStorageOrder(&trapPixelMap, COLUMNMAJOR);
       if (populateTrapPixelMap(&trapPixelMap, &pars.baseParams, acs->verbose, acs->expstart))
       {
           freeAll(&ptrReg);
           return status;
       }

       /* perform CTE correction */
       SingleGroup * cteCorrectedImage = &columnMajorImage;
       if (inverseCTEBlur(&smoothedImage, cteCorrectedImage, &trapPixelMap, &pars.baseParams))
       {
           freeAll(&ptrReg);
           return status;
       }
       freePtr(&ptrReg, &trapPixelMap);

       /* add readout noise back and convert corrected data back to DN.
           add 10% correction to error in quadrature. */
        double temp_err;
        for (unsigned k = 0; k < nRows; ++k)
        {
            for (unsigned m = 0; m < nColumns; ++m)
            {

                float delta = (PixColumnMajor(cteCorrectedImage->sci.data,k,m) - PixColumnMajor(smoothedImage.sci.data,k,m));

                temp_err = 0.1 * fabs(delta - Pix(raw.sci.data, m, k));
                Pix(raw.sci.data,m,k) += delta;

                float err2 = Pix(raw.err.data,m, k);
                err2 *= err2;
                Pix(raw.err.data,m, k) = sqrt(err2 + temp_err*temp_err);
            }
        }

        /* put the CTE corrected data back into the SingleGroup structure */
        if (unmake_amp_array(acs, &raw, amp, amp_xbeg, amp_ybeg, x))
        {
            freeAll(&ptrReg);
            return (status);
        }

        /* free space used by our amp arrays */
        freePtr(&ptrReg, &raw);
        freePtr(&ptrReg, &columnMajorImage);
        freePtr(&ptrReg, &smoothedImage);
        freePtr(&ptrReg, &pars.baseParams);
    }

    if (acs->printtime)
        TimeStamp("CTE corrections complete...","");

    freeAll(&ptrReg);
    return (status);
}


static int make_amp_array(const ACSInfo *acs, const SingleGroup *input,
                          const int amp,
                          const int xbeg, const int ybeg,
                          SingleGroup * output) {

    extern int status;

    const unsigned arr1 = output->sci.data.tot_ny;
    const unsigned arr2 = output->sci.data.tot_nx;

    if (acs->detector == WFC_CCD_DETECTOR) {

    if (amp != AMP_A || amp != AMP_B || amp != AMP_C || amp != AMP_D)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        status = ERROR_RETURN;
        return status;
    }

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
        int r;
        int c;
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned i = 0; i < arr1; ++i) {
            for (unsigned j = 0; j < arr2; ++j) {
                if (amp == AMP_A) {
                    r = ybeg + arr1 - i - 1;
                    c = xbeg + j;
                } else if (amp == AMP_B) {
                    r = ybeg + arr1 - i - 1;
                    c = xbeg + arr2 - j - 1;
                } else if (amp == AMP_C) {
                    r = ybeg + i;
                    c = xbeg + j;
                } else if (amp == AMP_D) {
                    r = ybeg + i;
                    c = xbeg + arr2 - j -1;
                }

                if (input->sci.data.data && output->sci.data.data)
                    Pix(output->sci.data, j, i) = Pix(input->sci.data, c, r);
                if (input->err.data.data && output->err.data.data)
                    Pix(output->err.data, j, i) = Pix(input->err.data, c, r);
            }
        }
        }
    } else {
        sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
        trlerror(MsgText);
        status = ERROR_RETURN;
        return status;
    }

    return status;
}


/* unmake_amp_array does the opposite of make_amp_array, it takes amp array
   views and puts them back into the single group in the right order.
*/
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *input,
                            const int amp,
                            const int xbeg, const int ybeg,
                            SingleGroup * output) {

    return make_amp_array(acs, input, amp, xbeg, ybeg, output);
}


