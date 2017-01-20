#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "acserr.h"

#include "pcte.h"
#include "../../../../ctegen2/ctegen2.h"
//#include "ctegen2.h"

static int get_amp_array_size(const ACSInfo *acs, SingleGroup *x,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);

static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2],
                          const Bool targetColumnMajor);

static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                            const int amp,
                            const int arr1, const int arr2,
                            const int xbeg, const int ybeg,
                            double amp_sci_array[arr1*arr2],
                            double amp_err_array[arr1*arr2],
                            const Bool sourceColumnMajor);

static int arrayToOrFromSingleGroup(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2],
                          const Bool targetColumnMajor,
                          const Bool toSingleGroup);

/* Perform a pixel based CTE correction on the SCI data extension of ACS CCD
   data. Parameters of the CTE characterization are read from the PCTE reference
   file.

   Originally adapted from code written by Jay Anderson and rewritten by
   Pey Lian Lim as a standalone application that operated on FLT files.

   MRD 10 Mar 2011
*/
int doPCTE (ACSInfo *acs, SingleGroup *x) {

    /* arguments:
       acs     i: calibration switches, etc
       x      io: image to be calibrated; written to in-place
    */

    extern int status;

    /* interpolated cte profile shape parameters */
    double chg_leak[MAX_TAIL_LEN*NUM_LOGQ];
    double chg_open[MAX_TAIL_LEN*NUM_LOGQ];

    /* interpolated profile charge parameters */
    double dtde_q[MAX_PHI];

    /* arrays interpolated at each parameterized charge level */
    double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV];
    double chg_open_lt[MAX_TAIL_LEN*NUM_LEV];
    double dpde_l[NUM_LEV];

    /* structure to hold CTE parameters from file */
    ACSCTEParams pars;

    char ccdamp[strlen(AMPSTR1)+1]; /* string to hold amps on current chip */
    int numamps;               /* number of amps on chip */
    int amp;                   /* index amp A:0, B:1, etc. */
    char * amploc = NULL;      /* pointer to amp character in AMPSORDER */
    int amp_arr1, amp_arr2;    /* int for C array dimension sizes */
    int amp_xsize, amp_ysize;  /* int for amp array size (x/y in CCD coords) */
    int amp_xbeg, amp_xend;    /* int for beg and end of amp arrays on chip */
    int amp_ybeg, amp_yend;

    /* make arrays to hold data amp by amp.
       these can be large arrays so it's best to declare them as pointers and
       get the space for the ararys using malloc */
    double * amp_sci_arr = NULL; /* original sci data */
    double * amp_err_arr = NULL; /* original err data */
    double * amp_sig_arr = NULL; /* decomposed signal */
    double * amp_nse_arr = NULL; /* decomposed readout error */
    double * amp_cor_arr = NULL; /* cte corrected data */

    /* in this algorithm each pixel has it's own CTE scaling,
       so we need an array for that. */
    double * cte_frac_arr = NULL;

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

    if (InterpolatePsi(pars.chg_leak, pars.psi_node, chg_leak, chg_open)) {
        return (status);
    }
    if (InterpolatePhi(pars.dtde_l, pars.q_dtde, pars.shft_nit, dtde_q)) {
        return (status);
    }
    if (FillLevelArrays(chg_leak, chg_open, dtde_q, pars.levels, chg_leak_lt,
                        chg_open_lt, dpde_l)) {
        return (status);
    }
    /********* done reading and calculating CTE model parameters ************/

    /* need to figure out which amps are on this chip */
    ccdamp[0] = '\0'; /* "reset" the string for reuse */
    parseWFCamps(acs->ccdamp, acs->chip, ccdamp);

    /* loop over amps on this chip and do CTE correction */
    numamps = strlen(ccdamp);
    for (unsigned ampIterator = 0; ampIterator < numamps; ++ampIterator) {
        sprintf(MsgText, "(pctecorr) Performing CTE correction for amp %c",
                ccdamp[ampIterator]);
        trlmessage(MsgText);

        /* get the amp letter and number where A:0, B:1, etc. */
        amploc = strchr(AMPSORDER, ccdamp[ampIterator]);
        amp = *amploc - AMPSORDER[0];

        /* get amp array size */
        if (get_amp_array_size(acs, x, amp, amploc, ccdamp,
                               &amp_xsize, &amp_ysize, &amp_xbeg,
                               &amp_xend, &amp_ybeg, &amp_yend)) {
            return (status);
        }

        amp_arr1 = amp_ysize;
        amp_arr2 = amp_xsize;

        //NOTE TO SELF - remove asserts here
        /* allocate space to hold this amp's data in its various forms */
        //Check whether amps can have different array sizes, if not move these out of amp iterator
        assert(amp_sci_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double)));
        assert(amp_err_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double)));
        assert(amp_sig_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double)));
        assert(amp_nse_arr = (double *) calloc((size_t)(amp_arr1*amp_arr2), sizeof(double))); //malloc(amp_arr1 * amp_arr2 * sizeof(double)));
        assert(amp_cor_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double)));
        assert(cte_frac_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double)));

        /* read data from the SingleGroup into an array containing data from
           just one amp */
        if (make_amp_array(acs, x, amp, amp_arr1, amp_arr2, amp_xbeg, amp_ybeg,
                           amp_sci_arr, amp_err_arr, True/*store column major*/)) {
            return (status);
        }

        /* fill cte_frac_arr

           To match cte_frac_arr: v8.1.3  ->  v8.2
               k*amp_arr2+m  -> (k+20)*(amp_arr2+24) + (m+24)

           offsety (LTV2) always 0 for fullframe and -1 for 2K subarray.
        */
        //loop over columns
        for (unsigned j = 0; j < amp_arr2; ++j)
        {
            //loop over rows
            for (unsigned i = 0; i < amp_arr1; ++i)
            {
                /*
                cte_frac_arr[k*amp_arr2 + m] = pars.cte_frac *
                    pars.col_scale[m*NAMPS + amp] *
                    ((double) k - acs->offsety) / CTE_REF_ROW;
                 */
                //cte_frac_arr[k*amp_arr2 + m] = pars.cte_frac *

                cte_frac_arr[i + j*amp_arr1] = pars.cte_frac *
                    pars.col_scale[j*NAMPS + amp] *
                    ((double) i - acs->offsety) / CTE_REF_ROW;
            }
        }

        /* do some smoothing on the data so we don't amplify the read noise.
           data should be in electrons. */
        if (DecomposeRN(amp_arr1, amp_arr2, amp_sci_arr,
                        pars.rn_clip, pars.noise_model,
                        amp_sig_arr, amp_nse_arr)) {
            return (status);
        }

        /* perform CTE correction */
        /* OLD ACS CTE correction
        if (FixYCte(amp_arr1, amp_arr2, amp_sig_arr, amp_cor_arr, pars.sim_nit,
                    pars.shft_nit, pars.sub_thresh, cte_frac_arr, pars.levels,
                    dpde_l, chg_leak_lt, chg_open_lt, acs->onecpu))
        */
        //double * pixDontKnowWhatThisIsYet = NULL; //For now, assuming cte_frac_arr???

        FloatTwoDArray cteRprof;
        FloatTwoDArray cteCprof;
        initFloatData(&cteRprof);
        initFloatData(&cteCprof);
        allocFloatData(&cteRprof, pars.baseParams.rprof->data.ny, pars.baseParams.rprof->data.nx, False);
        allocFloatData(&cteCprof, pars.baseParams.cprof->data.ny, pars.baseParams.cprof->data.nx, False);
        //Transpose arrays to column major
        //copyAndTransposeFloatData(&cteRprof, &pars.baseParams.rprof->data);
        //copyAndTransposeFloatData(&cteCprof, &pars.baseParams.cprof->data);

        //correctedColumn isn't needed - can be done in place (when all mem transposed to column major)
        double * correctedColumn = NULL;
        assert(correctedColumn = (double*)malloc(amp_arr1 * sizeof(*amp_sig_arr)));

        //Correct image one column at a time, amp_arr2 = nColumns
        //Before, with FixYCTE, the entire 2D array was passed in a split up internally by the function
        //loop over columns
        for (unsigned j = 0; j < amp_arr2; ++j)
        {

            // HORIZONTAL PRE/POST SCAN POPULATION
            //is this applicable in acs?
            Bool hasFlux = False;
            for (unsigned i = 0; i < amp_arr1; ++i)
            {
                if (correctedColumn[i] > 0)
                {
                    hasFlux = True;
                    break;
                }
            }
            if (!hasFlux)
                continue;

            unsigned NREDO = 0;
            Bool REDO = False; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            do { //while (redo) - post loop eval
                /* perform CTE correction */
                //for (int i = 0; i < pars->n_par; ++i)
                //{
               // if (simulateColumnReadout(correctedColumn, cte_frac_arr, &pars.baseParams, &cteRprof, &cteCprof, amp_arr1, pars.baseParams.n_par))
                 //   return status;
                //add damping

                //Do the last iteration separately so as not to dampen.
                //if (simulateColumnReadout(correctedColumn, cte_frac_arr, &pars.baseParams, &cteRprof, &cteCprof, amp_arr1, pars.baseParams.n_par))
                  //  return status;

               // REDO = pars.baseParams.thresh ? correctCROverSubtraction(cte_frac_arr, correctedColumn, amp_sig_arr, amp_arr1,
                 //       pars.baseParams.fix_rocr) : False;

            } while (REDO && ++NREDO < 5);
            //copy corrected column back into original array
            memcpy(&amp_sig_arr[j*amp_arr1], correctedColumn, sizeof(*amp_sig_arr)*amp_arr1);
        }
        if (correctedColumn)
        {
            free(correctedColumn);
            correctedColumn = NULL;
        }
        if (cte_frac_arr)
        {
            free(cte_frac_arr);
            correctedColumn = NULL;
        }
        freeFloatData(&cteRprof);
        freeFloatData(&cteCprof);

        // add readout noise back
        for (unsigned j = 0; j < amp_arr2; ++j)
        {
            for (unsigned i = 0; i < amp_arr1; ++i)
            {
                unsigned currentPixel = i + j*amp_arr1;
                amp_cor_arr[currentPixel] = amp_sig_arr[currentPixel] + amp_nse_arr[currentPixel];
            }
        }
        free(amp_nse_arr);
        amp_nse_arr = NULL;

        //Compute error and convert corrected data back to DN.
        // add 10% correction to error in quadrature.
        for (unsigned j = 0; j < amp_arr2; ++j)
        {
            for (unsigned i = 0; i < amp_arr1; ++i)
            {
                unsigned currentPixel = i + j*amp_arr1;
                double temp_err = 0.1 * fabs(amp_cor_arr[currentPixel] -
                                              amp_sci_arr[currentPixel]);

                amp_err_arr[currentPixel] = sqrt(
                            pow(amp_err_arr[currentPixel],2) + pow(temp_err,2));
            }
        }
        free(amp_sci_arr);
        amp_sci_arr = NULL;
        free(amp_sig_arr);
        amp_sig_arr = NULL;

        /* put the CTE corrected data back into the SingleGroup structure */
        if (unmake_amp_array(acs, x, amp, amp_arr1, amp_arr2, amp_xbeg, amp_ybeg,
                             amp_cor_arr, amp_err_arr, True /*unmake column major source*/)) {
            return (status);
        }

        /* free space used by our amp arrays */
        free(amp_err_arr);
        amp_err_arr = NULL;
        free(amp_cor_arr);
        amp_cor_arr = NULL;
    }

    if (acs->printtime) {
        TimeStamp("CTE corrections complete...","");
    }

    return (status);
}


/* Returns the x/y dimensions for an array that holds data readout through a
   single amp. currently only works for ACS WFC data.

   the standalone version has the array size hard wired since _flt files will
   always have 2048 x 2048 amp regions starting at pixel 0. Here we want to be
   a bit more careful because the overscan regions are still part of the data.

   the logic for figuring out the amp regions has been copied from doblev.
   - MRD 14 Mar 2011
*/
static int get_amp_array_size(const ACSInfo *acs, SingleGroup *x,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg, int *xend,
                              int *ybeg, int *yend) {
    extern int status;

    int bias_loc;
    int bias_ampx, bias_ampy;
    int bias_orderx[4] = {0,1,0,1};
    int bias_ordery[4] = {0,0,1,1};

    int trimx1, trimx2, trimy1, trimy2;

    if (acs->detector == WFC_CCD_DETECTOR) {
        /* Copy out overscan info for ease of reference in this function*/
        trimx1 = acs->trimx[0];
        trimx2 = acs->trimx[1];
        trimy1 = acs->trimy[0];
        trimy2 = acs->trimy[1];

        bias_loc = *amploc - ccdamp[0];
        bias_ampx = bias_orderx[bias_loc];
        bias_ampy = bias_ordery[bias_loc];

        /* Compute range of pixels affected by each amp */
        *xbeg = (trimx1 + acs->ampx) * bias_ampx;
        *xend = (bias_ampx == 0 && acs->ampx != 0) ? acs->ampx + trimx1 : x->sci.data.nx;
        *ybeg = (trimy1 + acs->ampy) * bias_ampy;
        *yend = (bias_ampy == 0 && acs->ampy != 0) ? acs->ampy + trimy1 : x->sci.data.ny;
        /* Make sure that xend and yend do not extend beyond the bounds of the
           image... WJH 8 Sept 2000
        */
        if (*xend > x->sci.data.nx) *xend = x->sci.data.nx;
        if (*yend > x->sci.data.ny) *yend = x->sci.data.ny;
        *xsize = *xend - *xbeg;
        *ysize = *yend - *ybeg;
    } else {
        sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
        trlerror(MsgText);
        status = ERROR_RETURN;
        return status;
    }

    return status;
}

/* Make_amp_array returns an array view of the data readout through the
   specified amp in which the amp is at the lower left hand corner.
*/
static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2],
                          const Bool targetColumnMajor)
{
    return arrayToOrFromSingleGroup(acs, im, amp, arr1, arr2, xbeg, ybeg,
            amp_sci_array, amp_err_array, targetColumnMajor, False);
}

/* unmake_amp_array does the opposite of make_amp_array, it takes amp array
   views and puts them back into the single group in the right order.
*/
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2],
                          const Bool targetColumnMajor)
{
    return arrayToOrFromSingleGroup(acs, im, amp, arr1, arr2, xbeg, ybeg,
            amp_sci_array, amp_err_array, targetColumnMajor, True);
}

static int arrayToOrFromSingleGroup(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2],
                          const Bool targetColumnMajor,
                          const Bool toSingleGroup) {

    extern int status;

    // variables for the image row/column we want
    unsigned rowOffset;
    unsigned columnOffset;
    unsigned row;
    unsigned column;
    int iSig;
    int jSig;

    if (acs->detector == WFC_CCD_DETECTOR)
    {
        switch(amp)
        {
        case AMP_A :
        {
            rowOffset = ybeg + arr1 - 1;
            columnOffset = xbeg;
            iSig = -1;
            jSig = 1;
            break;
        }
        case AMP_B :
        {
            rowOffset = ybeg + arr1 - 1;
            columnOffset = xbeg + arr2 - 1;
            iSig = -1;
            jSig = -1;
            break;
        }
        case AMP_C :
        {
            rowOffset = ybeg;
            columnOffset = xbeg;
            iSig = 1;
            jSig = 1;
            break;
        }
        case AMP_D :
        {
            rowOffset = xbeg + arr2 - 1;
            columnOffset = xbeg + arr2 - 1;
            iSig = 1;
            jSig = -1;
            break;
        }
        default :
        {
            trlerror("Amp number not recognized, must be 0-3.");
            status = ERROR_RETURN;
            return status;
        }
        }//end switch

        /* These loops could be moved into each of the above case stmnts to remove extra op of *Signatures
         * However, it maybe worth ignoring this for the sake of code neatness and maintainability.
         * NOTE: loop order optimized to column major storage for target.
         */
        if (targetColumnMajor)
        {
            if (toSingleGroup)
            {
                //loop over columns
                for (unsigned j = 0; j < arr2; ++j)
                {
                    column = columnOffset + jSig*j;
                    //loop over rows
                    for (unsigned i = 0; i < arr1; ++i)
                    {
                        row = rowOffset + iSig*i;
                        Pix(im->sci.data, column, row) = (float) amp_sci_array[i + j*arr1];
                        Pix(im->err.data, column, row) = (float) amp_err_array[i + j*arr1];
                    }
                }
            }
            else
            {
                //loop over columns
                for (unsigned j = 0; j < arr2; ++j)
                {
                    column = columnOffset + jSig*j;
                    //loop over rows
                    for (unsigned i = 0; i < arr1; ++i)
                    {
                        row = rowOffset + iSig*i;
                        amp_sci_array[i + j*arr1] = Pix(im->sci.data, column, row);
                        amp_err_array[i + j*arr1] = Pix(im->err.data, column, row);
                    }
                }
            }
        }
        else
        {
            if (toSingleGroup)
            {
                //loop over columns
                for (unsigned j = 0; j < arr2; ++j)
                {
                    column = columnOffset + jSig*j;
                    //loop over rows
                    for (unsigned i = 0; i < arr1; ++i)
                    {
                        row = rowOffset + iSig*i;
                        Pix(im->sci.data, column, row) = (float) amp_sci_array[i*arr2 + j];
                        Pix(im->err.data, column, row) = (float) amp_err_array[i*arr2 + j];
                    }
                }
            }
            else
            {
                //loop over columns
                for (unsigned j = 0; j < arr2; ++j)
                {
                    column = columnOffset + jSig*j;
                    //loop over rows
                    for (unsigned i = 0; i < arr1; ++i)
                    {
                        row = rowOffset + iSig*i;
                        amp_sci_array[i*arr2 + j] = Pix(im->sci.data, column, row);
                        amp_err_array[i*arr2 + j] = Pix(im->err.data, column, row);
                    }
                }
            }
        }
    }
    else
    {
        sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
        trlerror(MsgText);
        status = ERROR_RETURN;
        return status;
    }

    return status;
}
