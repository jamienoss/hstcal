#ifndef CTEGEN2_INCL
#define CTEGEN2_INCL

#include "hstio.h"

//This may need to differ between instruments, but oh well for now.
#define TRAPS 999 /*max number of traps per column = rows in pctetab[1], valid traps are < 999999 in qlev*/

typedef struct {
    unsigned n_forward; // number of forward modeling iterations
    unsigned n_par; // number of iterations in parallel transfer
    unsigned cte_traps; // number of valid TRAPS in file for reallocation
    int cte_len; // max length of cte trail
    int fix_rocr; /*make allowance for readout cosmic rays*/
    double thresh; /*over subtraction threshold*/
    double rn_amp; // read noise amplitude for clipping
    double qlevq_data[TRAPS]; // charge packet size in electrons
    double dpdew_data[TRAPS]; // trap size in electrons
    double cte_date0; /*date of uvis install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    double scale_frac; /*scaling of cte model relative to ctedate1*/
    FloatHdrData *rprof; // differential trail profile as image
    FloatHdrData *cprof; // cumulative trail profile as image
} CTEParams;

int inverseCTEBlur(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * fff, CTEParams * cte,
        const int verbose, const double expstart, const unsigned nRows, const unsigned nColumns);

int simulatePixelReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows);

int simulateColumnReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts);

Bool correctCROverSubtraction(double * const pix_ctef, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold);

#endif
