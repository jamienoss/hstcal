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
    double rn_amp; // read noise amplitude for clipping
    double qlevq_data[TRAPS]; // charge packet size in electrons
    double dpdew_data[TRAPS]; // trap size in electrons
    FloatHdrData *rprof; // differential trail profile as image
    FloatHdrData *cprof; // cumulative trail profile as image
} CTEParams;

int sim_colreadout_l(double * pixo, const double * pixf, const CTEParams * cte,
        const FloatTwoDArray * rprof, const FloatTwoDArray * cprof, const unsigned nRows, const Bool dampen);

#endif
