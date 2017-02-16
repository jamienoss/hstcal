/*
MLS 2015: read in the CTE parameters from the PCTETAB file

    Jan 19, 2016: MLS  Updated to check for existence of PCTENSMD in
                  raw science header

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wf3.h" //need to remove this dependency
#include "ctegen2.h"

/*
# include "hstio.h"
# include "ximio.h"
# include "xtables.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
# include "cte.h"
*/

void initCTEParams(CTEParams * pars, const unsigned _nTraps, const unsigned _nRows, const unsigned _nColumns)
{
    pars->nRows = _nRows;
    pars->nColumns = _nColumns;
    pars->nTraps = _nTraps;

    pars->noise_mit=0;
    pars->thresh=0;
    pars->cte_date0=0;
    pars->cte_date1=0;
    pars->cte_traps=0;
    pars->cte_len=0;
    pars->rn_amp=0;
    pars->n_forward=0;
    pars->n_par=0;
    pars->scale_frac=0; /*will be updated during routine run*/
    pars->fix_rocr = 0;

    pars->nRowsPerFullFrame = 0;
    pars->nColumnsPerFullFrame = 0;
    pars->nColumnsPerChip = 0;
    pars->nRowsPerChip = 0;
    pars->nColumnsPerQuad = 0;
    pars->nRowsPerQuad = 0;
    pars->subarrayRowOffset = 0;
    pars->subarrayColumnOffset = 0;

    pars->iz_data = NULL;
    pars->wcol_data = NULL;
    pars->scale512 = NULL;
    pars->scale1024 = NULL;
    pars->scale1536 = NULL;
    pars->scale2048 = NULL;
    pars->qlevq_data = NULL;
    pars->dpdew_data = NULL;

    pars->rprof = NULL; /*differential trail profile as image*/
    pars->cprof = NULL; /*cummulative trail profile as image*/

    *pars->cte_name='\0';
    *pars->cte_ver='\0';
    *pars->descrip2='\0';
}

void delete(void ** ptr)
{
    if (*ptr)
    {
        free(*ptr);
        *ptr = NULL;
    }
}

void * newAndZero(void ** ptr, size_t count, size_t size)
{
    //if (*ptr)
      //  return *ptr;
    delete(ptr);
    *ptr = calloc(count, size);
    return *ptr;
}

void allocateCTEParams(CTEParams * pars)
{
    newAndZero((void*)&pars->iz_data, pars->nColumns, sizeof(*pars->iz_data));
    newAndZero((void*)&pars->scale512, pars->nColumns, sizeof(*pars->scale512));
    newAndZero((void*)&pars->scale1024, pars->nColumns, sizeof(*pars->scale1024));
    newAndZero((void*)&pars->scale1536, pars->nColumns, sizeof(*pars->scale1536));
    newAndZero((void*)&pars->scale2048, pars->nColumns, sizeof(*pars->scale2048));
    newAndZero((void*)&pars->wcol_data, pars->nTraps, sizeof(*pars->wcol_data));
    newAndZero((void*)&pars->qlevq_data, pars->nTraps, sizeof(*pars->qlevq_data));
    newAndZero((void*)& pars->dpdew_data, pars->nTraps, sizeof(*pars->dpdew_data));

/*
    pars->iz_data = calloc(pars->nRows, sizeof(*pars->iz_data));
    pars->scale512 = calloc(pars->nRows, sizeof(*pars->scale512));
    pars->scale1024 = calloc(pars->nRows, sizeof(*pars->scale1024));
    pars->scale1536 = calloc(pars->nRows, sizeof(*pars->scale1536));
    pars->scale2048 = calloc(pars->nRows, sizeof(*pars->scale2048));
    pars->wcol_data = calloc(pars->nTraps, sizeof(*pars->wcol_data));
    pars->qlevq_data = calloc(pars->nTraps, sizeof(*pars->qlevq_data));
    pars->dpdew_data = calloc(pars->nTraps, sizeof(*pars->dpdew_data));
*/
    assert(pars->iz_data);
    assert(pars->wcol_data);
    assert(pars->scale512);
    assert(pars->scale1024);
    assert(pars->scale1536);
    assert(pars->scale2048);
    assert(pars->qlevq_data);
    assert(pars->dpdew_data);
}

void freeCTEParams(CTEParams * pars)
{
    delete((void*)&pars->iz_data);
    delete((void*)&pars->wcol_data);
    delete((void*)&pars->scale512);
    delete((void*)&pars->scale1024);
    delete((void*)&pars->scale1536);
    delete((void*)&pars->scale2048);
    delete((void*)&pars->qlevq_data);
    delete((void*)&pars->dpdew_data);

    assert(pars->iz_data == NULL);
    assert(pars->wcol_data == NULL);
    assert(pars->scale512 == NULL);
    assert(pars->scale1024 == NULL);
    assert(pars->scale1536 == NULL);
    assert(pars->scale2048 == NULL);
    assert(pars->qlevq_data == NULL);
    assert(pars->dpdew_data == NULL);

    freeFloatHdrData(pars->rprof);
    freeFloatHdrData(pars->cprof);
    delete((void*)&pars->rprof);
    delete((void*)&pars->cprof);
    assert(pars->rprof == NULL);
    assert(pars->cprof == NULL);
}

/************ HELPER SUBROUTINES ****************************/

int GetCTEPars (char *filename, CTEParams *pars) {
    /* Read the cte parameters from the reference table PCTETAB

       These are taken from the PCTETAB global header:
       CTE_NAME - name of cte algorithm
       CTE_VER - version number of cte algorithm
       CTEDATE0 - date of wfc3/uvis installation in HST, in fractional years
       CTEDATE1 - reference date of CTE model pinning, in fractional years

       PCTETLEN - max length of CTE trail
       PCTERNCL - readnoise amplitude and clipping level
       PCTESMIT - number of iterations used in CTE forward modeling
       PCTESHFT - number of iterations used in the parallel transfer
       PCTENSMD - readnoise mitigation algorithm
       PCTETRSH - over-subtraction threshold

       The table has 4 extensions:

Filename: wfc3_cte.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU      21   ()
1    QPROF       BinTableHDU     16   999R x 3C    ['i', 'i', 'i']
2    SCLBYCOL    BinTableHDU     20   8412R x 5C   ['i', 'e', 'e', 'e', 'e']
3    RPROF       ImageHDU        12   (999, 100)   float32
4    CPROF       ImageHDU        12   (999, 100)   float32

     */

    extern int status; /* variable for return status */

    /* VARIABLE FOR FILENAME + EXTENSION NUMBER. */
    char filename_wext[strlen(filename) + 4];


    /* NAMES OF DATA COLUMNS WE WANT FROM THE FILE, DATA WILL BE STORED IN THE PARS STRUCTURE */
    const char wcol[] = "W";
    const char qlevq[] = "QLEV_Q";
    const char dpdew[] = "DPDE_W";
    const char iz[]    = "IZ";
    const char sens512[]= "SENS_0512";
    const char sens1024[] = "SENS_1024";
    const char sens1536[] = "SENS_1536";
    const char sens2048[] = "SENS_2048";

    /* HSTIO VARIABLES */
    Hdr hdr_ptr;
    initHdr(&hdr_ptr);
    /* LOAD PRIMARY HEADER */
    if (LoadHdr(filename, &hdr_ptr)) {
        sprintf(MsgText,"(pctecorr) Error loading header from %s",filename);
        cteerror(MsgText);
        status = OPEN_FAILED;
        return status;
    }

    /* GET CTE_NAME KEYWORD */
    if (GetKeyStr (&hdr_ptr, "CTE_NAME", NO_DEFAULT, "", pars->cte_name, SZ_CBUF)) {
        cteerror("(pctecorr) Error reading CTE_NAME keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }
    sprintf(MsgText,"\nCTE_NAME: %s",pars->cte_name);
    trlmessage(MsgText);

    /* GET VERSION NUMBER  */
    if (GetKeyStr(&hdr_ptr, "CTE_VER", NO_DEFAULT, "", pars->cte_ver, SZ_CBUF)) {
        cteerror("(pctecorr) Error reading CTE_VER keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }
    sprintf(MsgText,"CTE_VER: %s",pars->cte_ver);
    trlmessage(MsgText);

    /* GET DATE OF UVIS INSTALLATION IN HST */
    if (GetKeyDbl(&hdr_ptr, "CTEDATE0", NO_DEFAULT, -999, &pars->cte_date0)) {
        cteerror("(pctecorr) Error reading CTEDATE0 keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"CTEDATE0: %g",pars->cte_date0);
    trlmessage(MsgText);

    /* GET REFRENCE DATE OF CTE MODEL PINNING */
    if (GetKeyDbl(&hdr_ptr, "CTEDATE1", NO_DEFAULT, -999, &pars->cte_date1)) {
        cteerror("(pctecorr) Error reading CTEDATE1 keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"CTEDATE1: %g",pars->cte_date1);
    trlmessage(MsgText);

    /* READ MAX LENGTH OF CTE TRAIL */
    if (GetKeyInt(&hdr_ptr, "PCTETLEN", NO_DEFAULT, -999, &pars->cte_len)) {
        cteerror("(pctecorr) Error reading PCTETLEN keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"PCTETLEN: %d",pars->cte_len);
    trlmessage(MsgText);

    /* GET READ NOISE CLIPPING LEVEL */
    if (GetKeyDbl(&hdr_ptr, "PCTERNOI", NO_DEFAULT, -999, &pars->rn_amp)) {
        cteerror("(pctecorr) Error reading PCTERNOI keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"PCTERNOI: %f",pars->rn_amp);
    trlmessage(MsgText);

    /* GET NUMBER OF ITERATIONS USED IN FORWARD MODEL */
    if (GetKeyInt(&hdr_ptr, "PCTENFOR", NO_DEFAULT, -999, &pars->n_forward)) {
        cteerror("(pctecorr) Error reading PCTENFOR keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }
    sprintf(MsgText,"PCTERNFOR: %d",pars->n_forward);
    trlmessage(MsgText);

    /* GET NUMBER OF ITERATIONS USED IN PARALLEL TRANSFER*/
    if (GetKeyInt(&hdr_ptr, "PCTENPAR", NO_DEFAULT, -999, &pars->n_par)) {
        cteerror("(pctecorr) Error reading PCTENPAR keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"PCTERNPAR: %d",pars->n_par);
    trlmessage(MsgText);

    /* GET READ NOISE MITIGATION ALGORITHM*/
    if (GetKeyInt(&hdr_ptr, "PCTENSMD", NO_DEFAULT, -999, &pars->noise_mit)) {
        cteerror("(pctecorr) Error reading PCTENSMD keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }
    sprintf(MsgText,"PCTENSMD: %d",pars->noise_mit);
    trlmessage(MsgText);

    /* GET OVER SUBTRACTION THRESHOLD */
    if (GetKeyDbl(&hdr_ptr, "PCTETRSH", NO_DEFAULT, -999, &pars->thresh)) {
        cteerror("(pctecorr) Error reading PCTETRSH keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    sprintf(MsgText,"PCTETRSH: %g",pars->thresh);
    trlmessage(MsgText);

    /*FIX THE READOUT CR'S? */
    if (GetKeyInt(&hdr_ptr, "FIXROCR", NO_DEFAULT, -999, &pars->fix_rocr)){
        cteerror("(pctecorr) Error reading FIXROCR keyword from PCTETAB");
        status = KEYWORD_MISSING;
        return status;
    }

    /*
    sprintf(MsgText,"FIXROCR: %d",pars->fix_rocr);
    trlmessage(MsgText);
    */


    /* DONE READING STUFF FROM THE PRIMARY HEADER */
    freeHdr(&hdr_ptr);
    /****************************************************************************/
    /* READ  DATA FROM FIRST TABLE EXTENSIONS */
    sprintf(filename_wext, "%s[%i]", filename, 1);

    /* OPEN  PARAMETERS FILE TO EXTENSION NUMBER 1 */
    IRAFPointer tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0); // xtables table pointer
    if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
        cteerror(MsgText);
        status = OPEN_FAILED;
        c_tbtclo(tbl_ptr);
        return status;
    }

    /* READ DATA FROM TABLE */
    /* get column pointer for w */
    IRAFPointer w_ptr = c_tbcfnd1_retPtr(tbl_ptr, wcol);
    if (c_iraferr() || !w_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",wcol);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* GET COLUMN POINTER FOR QLEVQ */
    IRAFPointer qlevq_ptr = c_tbcfnd1_retPtr(tbl_ptr, qlevq);
    if (c_iraferr() || !qlevq_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",qlevq);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* GET COLUMN POINTER FOR DPDEW */
    IRAFPointer dpdew_ptr = c_tbcfnd1_retPtr(tbl_ptr, dpdew);
    if (c_iraferr() || !dpdew_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",dpdew);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }


    // LOOP OVER TABLE ROWS UP TO SIZE TRAPS
    int ctraps = 0; // actual usable traps, i.e. see if more traps were added to reference file
    for (unsigned j = 0; j < pars->nTraps; ++j) {

        /* GET W FROM THIS ROW */
    	pars->wcol_data[j] = c_tbeGetInt(tbl_ptr, w_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, wcol);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }

        /* GET QLEVQ FROM THIS ROW */
        pars->qlevq_data[j] = c_tbeGetDouble(tbl_ptr, qlevq_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, qlevq);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }

        if (pars->qlevq_data[j] < 999999.)
            ctraps+=1;

        /* GET DPDEW FROM THIS ROW */
        pars->dpdew_data[j] = c_tbeGetDouble(tbl_ptr, dpdew_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, dpdew);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }
        if (ctraps > pars->nTraps){
            sprintf(MsgText,"More TRAPS in reference file than available, update TRAPS: %i -> %i",pars->nTraps,(int)ctraps);
            trlmessage(MsgText);
        }
    }

    /*IF CTRAPS EVER OVERFLOWS INT THIS NEEDS TO BE CHANGED*/
    pars->cte_traps = ctraps;

    /*
    sprintf(MsgText,"(pctecorr) data check for PCTETAB QPROF, row %i, %i\t%g\t%g\ttraps=%i\n",20,
            pars->wcol_data[19],pars->qlevq_data[19], pars->dpdew_data[19], pars->cte_traps);
    trlmessage(MsgText);
    */

    /* CLOSE CTE PARAMETERS FILE FOR EXTENSION 1*/
    c_tbtClose((void*)&tbl_ptr);
    assert(!tbl_ptr);
    /****************************************************************************/
    /****************************************************************************/
    /* READ CTE SCALING DATA FROM SECOND TABLE EXTENSION */
    sprintf(filename_wext, "%s[%i]", filename, 2);

    tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
        cteerror(MsgText);
        status = OPEN_FAILED;
        c_tbtclo(tbl_ptr);
        return status;
    }

    /*get column pointer for iz column*/
    IRAFPointer iz_ptr = c_tbcfnd1_retPtr(tbl_ptr, iz);
    if (c_iraferr() || iz_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",iz);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* get column pointer for sens512 */
    IRAFPointer sens512_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens512);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens512);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* get column pointer for sens1024 */
    IRAFPointer sens1024_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens1024);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens1024);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }
    /* get column pointer for sens1536 */
    IRAFPointer sens1536_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens1536);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens1536);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }
    /* get column pointer for sens2048 */
    IRAFPointer sens2048_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens2048);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens2048);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }


    /* read data from table */
    /* loop over table rows */
    for (unsigned j = 0; j < pars->nColumns; ++j)
    {
        /* get trap from this row */
    	pars->iz_data[j] = c_tbeGetInt(tbl_ptr, iz_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, iz);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale512[j] = c_tbeGetDouble(tbl_ptr, sens512_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens512);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }

        pars->scale1024[j] = c_tbeGetDouble(tbl_ptr, sens1024_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens1024);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale1536[j] = c_tbeGetDouble(tbl_ptr, sens1536_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens1536);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale2048[j] = c_tbeGetDouble(tbl_ptr, sens2048_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens2048);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
    }
    // for testing
    /*{
    unsigned j = pars->nColumns;
    sprintf(MsgText,"(pctecorr) data check for PCTETAB SCLBYCOL row %d, %d %g\t%g\t%g\t%g\ntotal traps = %i",
            j,pars->iz_data[j-1],pars->scale512[j-1],pars->scale1024[j-1],pars->scale1536[j-1],pars->scale2048[j-1],pars->cte_traps);
    trlmessage(MsgText);
    }
    */

    /* close CTE parameters file for extension 2*/
    c_tbtClose((void*)&tbl_ptr);
    assert(!tbl_ptr);

    /****************************************************************************/
    /*  extension 3: differential trail profile as image */
    ctemessage("Reading in image from extension 3");

    /* Get the coefficient images from the PCTETAB */
    pars->rprof = malloc(sizeof(*pars->rprof));
    if (pars->rprof == NULL){
        sprintf (MsgText, "Can't allocate memory for RPROF ref data");
        trlerror (MsgText);
        return (status = 1);
    }
    initFloatHdrData(pars->rprof);
    pars->rprof->data.storageOrder = COLUMNMAJOR;
    if (getFloatHD (filename, "RPROF", 1, pars->rprof)){
        return (status=1);
    }


    /****************************************************************************/
    /* ext number 4 : cummulative trail profile as image */
    ctemessage("Reading in image from extension 4");

    pars->cprof  = malloc(sizeof(*pars->cprof));
    if (pars->cprof == NULL){
        sprintf (MsgText, "Can't allocate memory for CPROF ref data");
        trlerror (MsgText);
        return (status = 1);
    }

    /* Get the coefficient images from the PCTETAB */
    initFloatHdrData (pars->cprof);
    pars->cprof->data.storageOrder = COLUMNMAJOR;
    if (getFloatHD (filename, "CPROF", 1, pars->cprof)){
        return (status=1);
    }

    return(status);
}

/*
  Some CTE parameters  will likely end up in the image headers so users can tune them.
  This will first check whether these parameters are set in the image header. If they
  are then the values found there will be used instead of the values read from
  the PCTETAB. If the values are not set they will be populated with the values
  from the PCTETAB, namely these:

        'CTE_NAME':'pixelCTE 2012', #name of cte algorithm
        'CTE_VER':'1.0' ,  #version number of algorithm
        'CTEDATE0':54962.0, #date of uvis installation in HST in MJD
        'CTEDATE1':56173.0, #reference date of cte model pinning in MJd
        'PCTETLEN':60, #max length of CTE trail
        'PCTERNOI':3.25, #read noise amplitude, clipping limit
        'PCTENFOR':5 ,#number of iterations used in cte forward modeling
        'PCTENPAR':7 ,#number of iterations used in parallel transfer
        'PCTENSMD':0 ,#read noise mitigation algorithm
        'PCTETRSH':-10.0 ,#over subtraction threshold, always use reference value
        'FIXROCR' : 1, #set to 1 for true, fix the readout cr's

 */
int CompareCTEParams(SingleGroup *group, CTEParams *pars) {

    extern int status;

    double rn_amp;
    int cte_len;
    int n_forward;
    int n_par;
    int fix_rocr;
    int noise_mit;

    cte_len=0;
    n_forward=0;
    n_par=0;
    fix_rocr=0;
    noise_mit=0;


    /*always put the cte_name  and cte_ver from the reference file in*/
    if (PutKeyStr(group->globalhdr,"CTE_NAME", pars->cte_name, "CTE algorithm name")){
        trlmessage("(pctecorr) Error updating CTE_NAME keyword in image header");
        return (status=HEADER_PROBLEM);
    }

    if(PutKeyStr(group->globalhdr,"CTE_VER", pars->cte_ver, "CTE algorithm version")){
        trlmessage("(pctecorr) Error updating CTE_VER keyword in image header");
        return (status=HEADER_PROBLEM);
    }

    if (PutKeyDbl(group->globalhdr, "CTEDATE0", pars->cte_date0,"Date of UVIS installation")) {
        trlmessage("(pctecorr) Error putting CTEDATE0 keyword in header");
        return (status=HEADER_PROBLEM);
    }

    if (PutKeyDbl(group->globalhdr, "PCTETRSH", pars->thresh,"cte oversubtraction threshold")) {
        trlmessage("(pctecorr) Error putting PCTETRSH keyword in header");
        return (status=HEADER_PROBLEM);
    }


    if (PutKeyDbl(group->globalhdr, "CTEDATE1", pars->cte_date1, "Date of CTE model pinning")) {
        trlmessage("(pctecorr) Error putting CTEDATE1 keyword in header");
        return (status=HEADER_PROBLEM);
    }

    /*check the PCTENSMD keyword in the header*/
    if (GetKeyInt(group->globalhdr, "PCTENSMD", NO_DEFAULT, -999, &noise_mit)){
        trlmessage("(pctecorr) Error reading PCTENSMD keyword from header");
        return (status=HEADER_PROBLEM);
    }
    if (noise_mit != pars->noise_mit){
        pars->noise_mit = noise_mit;
    }

    /*check the PCTEDIM keyword in header*/
    if (GetKeyInt(group->globalhdr, "PCTETLEN", NO_DEFAULT, -999, &cte_len)) {
        trlmessage("(pctecorr) Error reading PCTETLEN keyword from header");
        return (status=HEADER_PROBLEM);
    }

    if ( (cte_len != pars->cte_len) && (cte_len > 1) ){
        pars->cte_len=cte_len;
    } else {
        if (PutKeyInt(group->globalhdr,"PCTETLEN",pars->cte_len,"max length of CTE trail")){
            trlmessage("(pctecorr) Error updating PCTETLEN in header");
        return (status=HEADER_PROBLEM);
        }
    }

    /*check the PCTERNOI keyword in header*/
    if (GetKeyDbl(group->globalhdr, "PCTERNOI", NO_DEFAULT, -999, &rn_amp)) {
        trlmessage("(pctecorr) Error reading PCTERNOI keyword from header");
        return (status=HEADER_PROBLEM);
    }

    if ( (rn_amp >1.) && (rn_amp != pars->rn_amp)){
        pars->rn_amp=rn_amp;
    } else {
        if(PutKeyDbl(group->globalhdr, "PCTERNOI", pars->rn_amp,"read noise amp clip limit")){
            trlmessage("(pctecorr) Error updating PCTERNOI in header");
        return (status=HEADER_PROBLEM);
        }
    }


    /* get number of iterations used in forward model */
    if (GetKeyInt(group->globalhdr, "PCTENFOR", NO_DEFAULT, -999, &n_forward)) {
        trlmessage("(pctecorr) Error reading PCTENFOR keyword from header");
        return (status=HEADER_PROBLEM);
    }

    if (n_forward > 1 && n_forward != pars->n_forward){
        pars->n_forward = n_forward;
    } else {
        if (PutKeyInt(group->globalhdr, "PCTENFOR",pars->n_forward,"Number of iter in forward model")){
            trlmessage("(pctecorr) Error updating PCTENFOR in header");
            return (status=HEADER_PROBLEM);
        }
    }


    /* get number of iterations used in parallel transfer*/
    if (GetKeyInt(group->globalhdr, "PCTENPAR", NO_DEFAULT, -999, &n_par)) {
        trlmessage("(pctecorr) Error reading PCTENPAR keyword from header");
        return (status=HEADER_PROBLEM);
    }


    if( n_par >1 && n_par != pars->n_par){
        pars->n_par = n_par;
    } else {
        if (PutKeyInt(group->globalhdr, "PCTENPAR",pars->n_par,"Number of iter in parallel transfer")){
            trlmessage("(pctecorr) Error updating PCTENPAR in header");
            return (status=HEADER_PROBLEM);
        }
    }



    /*fix the readout Cr's? */
    if (GetKeyInt(group->globalhdr, "FIXROCR", NO_DEFAULT, -999, &fix_rocr)){
        trlmessage("(pctecorr) Error reading FIXROCR keyword from header");
        return(status = KEYWORD_MISSING);
    }
    if (0> fix_rocr && fix_rocr <=  1 && fix_rocr != pars->fix_rocr){
        pars->fix_rocr = fix_rocr;
    } else {
        if (PutKeyInt(group->globalhdr,"FIXROCR",pars->fix_rocr,"fix readout cosmic rays")){
            trlmessage("(pctecorr) Error updating FIXROCR keyword in header");
            return (status = KEYWORD_MISSING);
        }
    }
    return status;
}

void ctewarn (char *message) {

    char line[SZ_LINE+1];

    /* Use macro for beginning of Warning message */
    sprintf(line,"%s",WARN_PREFIX);
    strcat (line,message);

    ctemessage(line);
}

void cteerror (char *message) {

    char line[SZ_LINE+1];

    /* Use macro for beginning of Warning message */
    sprintf(line,"%s",ERR_PREFIX);
    strcat (line,message);

    ctemessage(line);
}

void ctemessage (char *message) {
    printf ("%s\n", message);
        fflush(stdout);
}