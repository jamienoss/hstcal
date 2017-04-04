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

static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int xbeg, const int ybeg,
                          SingleGroup * output
                          );

static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                            const int amp,
                            const int xbeg, const int ybeg,
                            SingleGroup * output
                            );

int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *x,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);


int doPCTEGen2 (ACSInfo *acs, SingleGroup *x)
{

    /* arguments:
       acs     i: calibration switches, etc
       x      io: image to be calibrated; written to in-place
    */

    extern int status;

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    /* interpolated cte profile shape parameters */
    double chg_leak[MAX_TAIL_LEN*NUM_LOGQ];
    double chg_open[MAX_TAIL_LEN*NUM_LOGQ];

    /* interpolated profile charge parameters */
    double dtde_q[MAX_PHI];

    /* arrays interpolated at each parameterized charge level */
    double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV];
    double chg_open_lt[MAX_TAIL_LEN*NUM_LEV];
    double dpde_l[NUM_LEV];



    /* temporary variable used during final error calculation */
    double temp_err;

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

    /* iteration variable */
    int i, k, m;

    char ccdamp[strlen(AMPSTR1)+1]; /* string to hold amps on current chip */
    int numamps;               /* number of amps on chip */
    int amp;                   /* index amp A:0, B:1, etc. */
    char * amploc;             /* pointer to amp character in AMPSORDER */
    int nRows, nColumns;    /* int for C array dimension sizes */
    int amp_xsize, amp_ysize;  /* int for amp array size (x/y in CCD coords) */
    int amp_xbeg, amp_xend;    /* int for beg and end of amp arrays on chip */
    int amp_ybeg, amp_yend;

    /* make arrays to hold data amp by amp.
       these can be large arrays so it's best to declare them as pointers and
       get the space for the ararys using malloc */
    double * amp_sci_arr; /* original sci data */
    double * amp_err_arr; /* original err data */
    double * amp_sig_arr; /* decomposed signal */
    double * amp_nse_arr; /* decomposed readout error */
    double * amp_cor_arr; /* cte corrected data */

    /* in this algorithm each pixel has it's own CTE scaling,
       so we need an array for that. */
    double * cte_frac_arr;

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
/*     initCTEParams(&pars.baseParams, TRAPS, nRows, nColumns);
     addPtr(&ptrReg, &pars.baseParams, &freeCTEParams);
     allocateCTEParams(&pars.baseParams);
     if (GetCTEPars (acs->pcte.name, &pars.baseParams))
     {
         freeAll(&ptrReg);
         return (status);
     }
     */
    if (PixCteParams(acs->pcte.name, acs->expstart, &pars)) {
        return (status);
    }

    if (CompareCteParams(x, &pars)) {
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
    for (i = 0; i < numamps; i++) {
        sprintf(MsgText, "(pctecorr) Performing CTE correction for amp %c",
                ccdamp[i]);
        trlmessage(MsgText);

        /* get the amp letter and number where A:0, B:1, etc. */
        amploc = strchr(AMPSORDER, ccdamp[i]);
        amp = *amploc - AMPSORDER[0];

        /* get amp array size */
        if (get_amp_array_size_acs_cte(acs, x, amp, amploc, ccdamp,
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


        /* allocate space to hold this amp's data in its various forms */
        //amp_sci_arr = (double *) malloc(nRows * nColumns * sizeof(double));
        /*amp_err_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
        amp_sig_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
        amp_nse_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
        amp_cor_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
        cte_frac_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
*/

        //This is used for the final output
       SingleGroup raw;
       initSingleGroup(&raw);
       addPtr(&ptrReg, &raw, &freeSingleGroup);
       allocSingleGroup/*SciOnly*/(&raw, nColumns, nRows, False);

        /* read data from the SingleGroup into an array containing data from
           just one amp */
        if (make_amp_array(acs, x, amp, amp_xbeg, amp_ybeg,
                          &raw)) {
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


        /* fill cte_frac_arr

           To match cte_frac_arr: v8.1.3  ->  v8.2
               k*amp_arr2+m  -> (k+20)*(amp_arr2+24) + (m+24)

           offsety (LTV2) always 0 for fullframe and -1 for 2K subarray.
        */
       /* for (k = 0; k < amp_arr1; k++) {
            for (m = 0; m < amp_arr2; m++) {
                cte_frac_arr[k*amp_arr2 + m] = pars.cte_frac *
                    pars.col_scale[m*NAMPS + amp] *
                    ((double) k - acs->offsety) / CTE_REF_ROW;
            }
        }
*/
        /* do some smoothing on the data so we don't amplify the read noise.
           data should be in electrons. */
       if (cteSmoothImage(&columnMajorImage, &smoothedImage, &pars.baseParams, pars.baseParams.rn_amp /*pars.rn_clip*/, max_threads, acs->verbose))
       {
           freeAll(&ptrReg);
           return (status);
       }
     /*   if (DecomposeRN(nRows, nColumns, amp_sci_arr,
                        pars.rn_clip, pars.noise_model,
                        amp_sig_arr, amp_nse_arr)) {
            return (status);
        }
*/

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
     /*   if (FixYCte(nRows, nColumns, amp_sig_arr, amp_cor_arr, pars.sim_nit,
                    pars.shft_nit, pars.sub_thresh, cte_frac_arr, pars.levels,
                    dpde_l, chg_leak_lt, chg_open_lt, acs->onecpu)) {
            return (status);
        }
*/


        /* add readout noise back and convert corrected data back to DN.
           add 10% correction to error in quadrature. */
        for (k = 0; k < nRows; k++) {
            for (m = 0; m < nColumns; m++) {

                float delta = (PixColumnMajor(cteCorrectedImage->sci.data,k,m) - PixColumnMajor(smoothedImage.sci.data,k,m));///ccdgain;
                //amp_cor_arr[k*nColumns + m] = amp_cor_arr[k*nColumns + m] +
                  //  amp_nse_arr[k*nColumns + m];

                temp_err = 0.1 * fabs(delta - Pix(raw.sci.data, m, k));
                Pix(raw.sci.data,m,k) += delta;

                //temp_err = 0.1 * fabs(amp_cor_arr[k*nColumns + m] -
                                   //   amp_sci_arr[k*nColumns + m]);
                float err2 = Pix(raw.err.data,m, k);
                err2 *= err2;
                Pix(raw.err.data,m, k) = sqrt(err2 + temp_err*temp_err);
               // amp_err_arr[k*nColumns + m] = sqrt(
                //    pow(amp_err_arr[k*nColumns + m],2) + pow(temp_err,2));
            }
        }

        /* put the CTE corrected data back into the SingleGroup structure */
        if (unmake_amp_array(acs, &raw, amp, amp_xbeg, amp_ybeg,
                x)) {
            return (status);
        }

        /* free space used by our amp arrays */
       /* free(amp_sci_arr);
        free(amp_err_arr);
        free(amp_sig_arr);
        free(amp_nse_arr);
        free(amp_cor_arr);
        free(cte_frac_arr);
        */
        freePtr(&ptrReg, &raw);
        freePtr(&ptrReg, &columnMajorImage);
        freePtr(&ptrReg, &smoothedImage);

        freePtr(&ptrReg, &pars.baseParams);
        //freeCTEParams(CTEParams * &pars.baseParams);
    }

    if (acs->printtime) {
        TimeStamp("CTE corrections complete...","");
    }

    freeAll(&ptrReg);
    return (status);
}


static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int xbeg, const int ybeg,
                          SingleGroup * output) {

    extern int status;

    /* iteration variables */
    int i, j;

    /* variables for the image row/column we want */
    int r, c;

    const unsigned arr1 = output->sci.data.tot_ny;
    const unsigned arr2 = output->sci.data.tot_nx;

    if (acs->detector == WFC_CCD_DETECTOR) {
        for (i = 0; i < arr1; i++) {
            for (j = 0; j < arr2; j++) {
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
                } else {
                    trlerror("Amp number not recognized, must be 0-3.");
                    status = ERROR_RETURN;
                    return status;
                }

                if (im->sci.data.data && output->sci.data.data)
                    Pix(output->sci.data, j, i) = Pix(im->sci.data, c, r);
                if (im->err.data.data && output->err.data.data)
                    Pix(output->err.data, j, i) = Pix(im->err.data, c, r);
                //amp_sci_array[i*arr2 + j] = Pix(im->sci.data, c, r);
                //amp_err_array[i*arr2 + j] = Pix(im->err.data, c, r);
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
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                            const int amp,
                            const int xbeg, const int ybeg,
                            SingleGroup * output) {

    extern int status;

    /* iteration variables */
    int i, j;

    /* variables for the image row/column we want */
    int r, c;
    const unsigned arr1 = im->sci.data.tot_ny;
    const unsigned arr2 = im->sci.data.tot_nx;

    if (acs->detector == WFC_CCD_DETECTOR) {
        for (i = 0; i < arr1; i++) {
            for (j = 0; j < arr2; j++) {
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
                } else {
                    trlerror("Amp number not recognized, must be 0-3.");
                    status = ERROR_RETURN;
                    return status;
                }

                if (im->sci.data.data && output->sci.data.data)
                    Pix(output->sci.data, c, r) = Pix(im->sci.data, j, i);

                if (im->err.data.data && output->err.data.data)
                    Pix(output->err.data, c, r) = Pix(im->err.data, j, i);
                //Pix(im->sci.data, c, r) = (float) amp_sci_array[i*arr2 + j];
                //Pix(im->err.data, c, r) = (float) amp_err_array[i*arr2 + j];
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


