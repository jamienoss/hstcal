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

int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *amp,
                              const int ampID, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);

static void extractAmp(SingleGroup * amp, const const SingleGroup * image, const unsigned ampID);
static void insertAmp(SingleGroup * amp, const const SingleGroup * image, const unsigned ampID);
static int alignAmpData(FloatTwoDArray * amp, const unsigned ampID);
static int alignAmp(SingleGroup * amp, const unsigned ampID);

int doPCTEGen2 (ACSInfo *acs, SingleGroup * chipImage)
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
    int ampID;                   /* index amp A:0, B:1, etc. */
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

    /* need to figure out which amps are on this chip */
    ccdamp[0] = '\0'; /* "reset" the string for reuse */
    parseWFCamps(acs->ccdamp, acs->chip, ccdamp);

    /* loop over amps on this chip and do CTE correction */
    numamps = strlen(ccdamp);
    Bool headerRead = False;
    Bool headerWritten = False;
    CTEParamsFast pars;
    addPtr(&ptrReg, &pars, &freeCTEParamsFast);
    unsigned nScaleTableColumns = 8412;
    initCTEParamsFast(&pars, TRAPS, 0, 0, nScaleTableColumns);
    allocateCTEParamsFast(&pars);
    if (GetCTEParsFast (acs->pcteTabNameFromCmd, &pars))
    {
        freeAll(&ptrReg);
        return (status);
    }

    for (unsigned i = 0; i < numamps; ++i)
    {
        sprintf(MsgText, "(pctecorr) Performing CTE correction for amp %c", ccdamp[i]);
        trlmessage(MsgText);

        /* get the amp letter and number where A:0, B:1, etc. */
        amploc = strchr(AMPSORDER, ccdamp[i]);
        ampID = *amploc - AMPSORDER[0];

        /* get amp array size */
        if (get_amp_array_size_acs_cte(acs, chipImage, ampID, amploc, ccdamp,
                               &amp_xsize, &amp_ysize, &amp_xbeg,
                               &amp_xend, &amp_ybeg, &amp_yend))
        {
            freeAll(&ptrReg);
            return (status);
        }

        nRows = amp_ysize;
        nColumns = amp_xsize;

        /**************** read and calculate parameters of CTE model ************/
        if (!headerRead)
        {
            headerRead = True;
           // initCTEParamsFast(&pars, TRAPS, nRows, nColumns);
           // allocateCTEParamsFast(&pars);
            if (GetCTEParsFast (acs->pcteTabNameFromCmd, &pars))
            {
                freeAll(&ptrReg);
                return (status);
            }

            pars.nRows = nRows;
            pars.nColumns = nColumns;
            pars.nRowsPerFullFrame = nRows;
            pars.nColumnsPerFullFrame = nColumns;
            pars.nColumnsPerChip = pars.nColumnsPerFullFrame/2;
            pars.nRowsPerChip = pars.nRowsPerFullFrame;
            pars.nColumnsPerQuad = pars.nColumnsPerFullFrame/4;
            pars.nRowsPerQuad = pars.nRowsPerFullFrame;
            pars.isSubarray = acs->subarray;
            pars.refAndIamgeBinsIdenticle = True;
            pars.columnOffset = 0;//amp_xbeg;//acs->offsetx;
            pars.rowOffset = 0;//amp_ybeg;//acs->offsety;
        }

        //warning remove fabs
        pars.scale_frac = fabs((acs->expstart - pars.cte_date0) / (pars.cte_date1 - pars.cte_date0));//0.041907;//(acs->expstart - pars.cte_date0) / (pars.cte_date1 - pars.cte_date0);
        if (!headerWritten)// && (acs->chip == 2 || acs->subarray == YES))
        {
            headerWritten = True;
            if (PutKeyDbl(chipImage->globalhdr, "PCTEFRAC", pars.scale_frac, "CTE scaling factor"))
            {
                trlerror("(pctecorr) Error writing PCTEFRAC to image header");
                freeAll(&ptrReg);
                return (status = HEADER_PROBLEM);
            }

            sprintf(MsgText, "(pctecorr) Read noise level PCTERNCL: %f", pars.rn_amp);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) Readout simulation forward modeling iterations PCTENFOR: %i",
                    pars.n_forward);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) Number of iterations used in the parallel transfer PCTENPAR: %i",
                    pars.n_par);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) CTE_FRAC: %f", pars.scale_frac);
            trlmessage(MsgText);
        }

        //This is used for the final output
        SingleGroup ampImage;
        initSingleGroup(&ampImage);
        addPtr(&ptrReg, &ampImage, &freeSingleGroup);
        allocSingleGroupExts(&ampImage, nColumns, nRows, SCIEXT | ERREXT, False);

       // if (ampID != AMP_C)
         //   continue;

        // read data from the SingleGroup into an array containing data from
        // just one amp
        extractAmp(&ampImage, chipImage, ampID);
        if (alignAmp(&ampImage, ampID))
        {
            freeAll(&ptrReg);
            return (status);
        }

        //copy to column major storage
        SingleGroup columnMajorImage;
        initSingleGroup(&columnMajorImage);
        addPtr(&ptrReg, &columnMajorImage, &freeSingleGroup);
        allocSingleGroupExts(&columnMajorImage, nColumns, nRows, SCIEXT | ERREXT, False);
        assert(!copySingleGroup(&columnMajorImage, &ampImage, COLUMNMAJOR));

        //CALCULATE THE SMOOTH READNOISE IMAGE
        trlmessage("CTE: Calculating smooth readnoise image");
        SingleGroup smoothedImage;
        initSingleGroup(&smoothedImage);
        addPtr(&ptrReg, &smoothedImage, &freeSingleGroup);
        allocSingleGroupExts(&smoothedImage, nColumns, nRows, SCIEXT | ERREXT, False);
        setStorageOrder(&smoothedImage, COLUMNMAJOR);

        // do some smoothing on the data so we don't amplify the read noise.
        if (cteSmoothImage(&columnMajorImage, &smoothedImage, &pars, pars.rn_amp, max_threads, acs->verbose))
        {
            freeAll(&ptrReg);
            return (status);
       }

       SingleGroup trapPixelMap;
       initSingleGroup(&trapPixelMap);
       addPtr(&ptrReg, &trapPixelMap, &freeSingleGroup);
       allocSingleGroupExts(&trapPixelMap, nColumns, nRows, SCIEXT | ERREXT, False);
       setStorageOrder(&trapPixelMap, COLUMNMAJOR);
       if (populateTrapPixelMap(&trapPixelMap, &pars, acs->verbose))
       {
           freeAll(&ptrReg);
           return status;
       }

       //perform CTE correction
       SingleGroup * cteCorrectedImage = &columnMajorImage;
       if (inverseCTEBlur(&smoothedImage, cteCorrectedImage, &trapPixelMap, &pars))
       {
           freeAll(&ptrReg);
           return status;
       }
       freePtr(&ptrReg, &trapPixelMap);

        // add 10% correction to error in quadrature.
        double temp_err;
#ifdef _OPENMP
        #pragma omp parallel for shared(nRows, nColumns, cteCorrectedImage, smoothedImage, ampImage) schedule(static)
#endif
        //This loop order is row major as more row major storage is accessed than column.
        //MORE: look into worth splitting out ops - prob needs a order swap (copy) so perhaps not worth it.
        for (unsigned k = 0; k < nRows; ++k)
        {
            for (unsigned m = 0; m < nColumns; ++m)
            {
                float delta = (PixColumnMajor(cteCorrectedImage->sci.data,k,m) - PixColumnMajor(smoothedImage.sci.data,k,m));

                temp_err = 0.1 * fabs(delta - Pix(ampImage.sci.data, m, k));
                Pix(ampImage.sci.data, m, k) += delta;

                float err2 = Pix(ampImage.err.data, m, k);
                err2 *= err2;
                Pix(ampImage.err.data, m, k) = sqrt(err2 + temp_err*temp_err);
            }
        }

        // put the CTE corrected data back into the SingleGroup structure
        if (alignAmp(&ampImage, ampID))//unmake_amp_array(image, &raw, acs, amp, nRows, nColumns, amp_xbeg, amp_ybeg))
        {
            freeAll(&ptrReg);
            return (status);
        }
        insertAmp(chipImage, &ampImage, ampID);

        /* free space used by our amp arrays */
        freePtr(&ptrReg, &ampImage);
        freePtr(&ptrReg, &columnMajorImage);
        freePtr(&ptrReg, &smoothedImage);
    }

    if (acs->printtime)
        TimeStamp("CTE corrections complete...","");

    freeAll(&ptrReg);
    return (status);
}


static void extractAmp(SingleGroup * amp,  const SingleGroup * image, const unsigned ampID)
{
    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return;

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = 0;
    if (ampID == AMP_B || ampID == AMP_D)
    	offset = nColumns;
    else if (ampID != AMP_A && ampID != AMP_C)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        status = ERROR_RETURN;
        return;
    }

    copyOffsetSingleGroup(amp, image, nRows, nColumns, 0, offset, nColumns, rowSkipLength);

    return;
}

static void insertAmp(SingleGroup * image, const SingleGroup * amp, const unsigned ampID)
{
    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return;

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = 0;
    if (ampID == AMP_B || ampID == AMP_D)
    	offset = nColumns;
    else if (ampID != AMP_A && ampID != AMP_C)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        status = ERROR_RETURN;
        return;
    }

    copyOffsetSingleGroup(image, amp, nRows, nColumns, offset, 0, rowSkipLength, nColumns);

    return;
}

static int alignAmp(SingleGroup * amp, const unsigned ampID)
{
    if (!amp)
        return 0;//be silent

    //sci data
    if (amp->sci.data.data)
    {
        if ((status = alignAmpData(&amp->sci.data, ampID)))
            return status;
    }
    //err data
    if (amp->err.data.data)
    {
        if ((status = alignAmpData(&amp->err.data, ampID)))
            return status;
    }
    //dq data
    if (amp->dq.data.data)
     assert(0);//unimplemented

    return status;
}

static int alignAmpData(FloatTwoDArray * amp, const unsigned ampID)
{
    //Align amps such that they are at the bottom left
    extern int status;

    if (!amp || !amp->data)
        return 0;//be silent

    //WARNING - assumes row major storage
    assert(amp->storageOrder == ROWMAJOR);

    unsigned nColumns = amp->nx;
    unsigned nRows = amp->ny;

    if (ampID == AMP_C)
        return status;

    if (ampID == AMP_B || ampID == AMP_D)
    {
        //grab a row, flip it, put it back
        unsigned rowLength = nColumns;
        float * row = NULL;
        for (unsigned i = 0; i < nRows; ++i)
        {
            //find row
            row = amp->data + i*nColumns;
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
    if (ampID == AMP_A || ampID == AMP_B)
    {
        //either physically align all or propagate throughout a mechanism to work on the array upside down (flip b quad though)
        //we'll just flip all for now. See if there's info in the header specifying amp location rel to pix in file,
        //i.e. a way to know whether the chip is 'upside down'. Could then reverse cte corr trail walk direction
        float * tempRow = NULL;
        size_t rowSize = nColumns*sizeof(*tempRow);
        tempRow = malloc(rowSize);
        float * topRow = NULL;
        float * bottomRow = NULL;
        for (unsigned i = 0; i < nRows/2; ++i)
        {
            topRow = amp->data + i*nColumns;
            bottomRow = amp->data + (nRows-1-i)*nColumns;
            memcpy(tempRow, topRow, rowSize);
            memcpy(topRow, bottomRow, rowSize);
            memcpy(bottomRow, tempRow, rowSize);
        }
        free(tempRow);
    }

    return status;
}
