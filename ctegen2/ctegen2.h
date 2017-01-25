#ifndef CTEGEN2_INCL
#define CTEGEN2_INCL

#include "hstio.h"

#define PixColumnMajor(a,i,j) (a).data[(j)*(a).tot_ny + (i)]

typedef struct {
    int noise_mit; /*read noise mitigation algorithm*/
    int cte_len; // max length of cte trail
    int fix_rocr; /*make allowance for readout cosmic rays*/
    unsigned nTraps;
    unsigned nRows;
    unsigned nColumns;
    int n_forward; // number of forward modeling iterations
    int n_par; // number of iterations in parallel transfer
    unsigned cte_traps; // number of valid TRAPS in file for reallocation
    double thresh; /*over subtraction threshold*/
    double rn_amp; // read noise amplitude for clipping
    double cte_date0; /*date of uvis install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    double scale_frac; /*scaling of cte model relative to ctedate1*/


    FloatHdrData * rprof; // differential trail profile as image
    FloatHdrData * cprof; // cumulative trail profile as image

    int * wcol_data; //trap number, insync with number of traps
    int    * iz_data; //column number in raz format
    double * scale512;  //scaling appropriate at row 512
    double * scale1024; //scaling appropriate at row 1024
    double * scale1536; //scaling appropriate at row 1536
    double * scale2048; //scaling appropriate at row 2048
    double * qlevq_data; // charge packet size in electrons
    double * dpdew_data; // trap size in electrons


    char descrip2[256]; /*descrip from table row, not read in for cte purposes*/
    char cte_name[256]; /*name of cte algorithm */
    char cte_ver[256]; /*version of algorithm */
} CTEParams;

int inverseCTEBlurWithRowMajorIput(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart);
int inverseCTEBlur(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParams * cte,
        const int verbose, const double expstart);

int simulatePixelReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows);

int simulateColumnReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts);

Bool correctCROverSubtraction(double * const pix_ctef, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold);

int populateTrapPixelMap(SingleGroup * input, CTEParams * params);
int cteSmoothImage(const SingleGroup * input, SingleGroup * output, double readNoiseAmp, unsigned maxThreads , int verbose);
int cteSmoothImageRowMajor(const SingleGroup * input, SingleGroup * output, double readNoiseAmp, unsigned maxThreads, int verbose);
double find_dadj(const unsigned i ,const unsigned j, const unsigned nRows, const float * obsloc[3], const float * rszloc[3], const double readNoiseAmp);

//helpers
void * newAndZero(void ** ptr, const size_t count, const size_t size);
void delete(void ** ptr);
void initCTEParams(CTEParams * pars, const unsigned _nTraps, const unsigned _nRows, const unsigned _nColumns);
void allocateCTEParams(CTEParams * pars);
void freeCTEParams(CTEParams * pars);

int CompareCTEParams(SingleGroup * input, CTEParams * params);
int GetCTEPars (char *filename, CTEParams * params);
void ctewarn (char *message);
void cteerror (char *message);
void ctemessage (char *message);

//These shouldn't be here, they belong in their own header, in a central hstcal/inlclude with the code
//needing to be in hstcal/lib
int LoadHdr (char *, Hdr *);
int GetKeyStr (Hdr *, char *, int, char *, char *, int);
int GetKeyInt (Hdr *, char *, int , int , int *);
int GetKeyDbl (Hdr *, char *, int , double , double *);
int PutKeyInt (Hdr *, char *, int , char *);
int PutKeyDbl (Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);

#endif
