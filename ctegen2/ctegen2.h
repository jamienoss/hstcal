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
    double thresh; /*over subtraction threshold*/
    int fix_rocr; /*make allowance for readout cosmic rays*/
    double rn_amp; // read noise amplitude for clipping
    double qlevq_data[TRAPS]; // charge packet size in electrons
    double dpdew_data[TRAPS]; // trap size in electrons
    FloatHdrData *rprof; // differential trail profile as image
    FloatHdrData *cprof; // cumulative trail profile as image
} CTEParams;

int simulatePixelReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows);

int simulateColumnReadout(double * const pixelColumn, const double * const pixf, const CTEParams * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts);

Bool correctCROverSubtraction(double * const pix_ctef, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold);

#endif
