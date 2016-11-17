#include <stdio.h>
#include <math.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"
#include "pcte.h"
#include "acs.h"

/*This is the workhorse subroutine; it simulates the readout
  of one column pixi() and outputs this to pixo() using a single
  iteration.  It can be called successively to do the transfer
  in steps.


  JDIM == RAZ_ROWS
  WDIM == TRAPS  Ws is the input traps number < 999999
  NITs == cte_pars->n_par

  These are already in the parameter structure CTEParams
  int     Ws              the number of traps < 999999
  float     q_w[TRAPS];     the run of charge with level  == qlevq_data
  float   dpde_w[TRAPS];  the run of charge loss with level == dpdew_data
  float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel == rprof fits image
  float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)  == cprof fits image


  W = wcol_data = trap id

  q_w[TRAP] = qlev_q from QPROF  traps as function of packet size = cte->qlevq_data[TRAP]

  pixi (curr), pixo (read) , pixf(cteff) are passed and are 1d arrays which have values for a particular column

  the ttrap reference to the image array has to be -1 for C
  */

int sim_colreadout_l(double *pixo, const double *pixf, const CTEParams *cte, const unsigned nRows)
{
    extern int status;

    //for (unsigned nthIter = 0; nthIter < nIterations; ++nthIter)//look into moving this further down
    //{
        double padd_3=0.0;
        double prem_3=0.0;
        double padd_2=0.0;
        double fcarry=0.0;
        double pix_1=0.0;
        double ftrap=0.0;
        int ttrap=0;

        /*from the reference table*/
        const FloatHdrData * rprof = cte->rprof;
        const FloatHdrData * cprof = cte->cprof;

        /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
          PMAX SHOULD ALWAYS BE POSITIVE HERE*/
        double pmax = 10.;
        //Look into whether this really has to be computed each iteration?
        for (unsigned j = 0; j < nRows; ++j)
        {
            pmax = (pixo[j] < pmax) ? pmax : pixo[j];
            /*if (pixo[j] > pmax)
                pmax=pixo[j];
                */
        }

        /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
          AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
        for (int w = cte->cte_traps-1; w >= 0; --w){
            if (cte->qlevq_data[w] <= pmax)
            {
                ftrap = 0.0e0;
                ttrap = cte->cte_len; /*for referencing the image at 0*/
                fcarry = 0.0e0;

                /*GO UP THE COLUMN PIXEL BY PIXEL*/
                for(unsigned j = 0; j < nRows; ++j){
                    pix_1 = pixo[j];

                    if ( (ttrap < cte->cte_len) || ( pix_1 >= cte->qlevq_data[w] - 1. ) ){
                        if (pixo[j] >= 0 ){
                            pix_1 = pixo[j] + fcarry; /*shuffle charge in*/
                            fcarry = pix_1 - floor(pix_1); /*carry the charge remainder*/
                            pix_1 = floor(pix_1); /*reset pixel*/
                        }

                        /*HAPPENS AFTER FIRST PASS*/
                        /*SHUFFLE CHARGE IN*/
                        if ( j> 0  ) {
                            if (pixf[j] < pixf[j-1])
                                ftrap *= (pixf[j] /  pixf[j-1]);
                        }

                        /*RELEASE THE CHARGE*/
                        padd_2=0.0;
                        if (ttrap <cte->cte_len){
                            ttrap += 1;
                            padd_2 = Pix(rprof->data,w,ttrap-1) *ftrap;
                        }

                        padd_3 = 0.0;
                        prem_3 = 0.0;
                        if ( pix_1 >= cte->qlevq_data[w]){
                            prem_3 =  cte->dpdew_data[w] / cte->n_par * pixf[j];  /*dpdew is 1 in file */
                            if (ttrap < cte->cte_len)
                                padd_3 = Pix(cprof->data,w,ttrap-1)*ftrap;
                            ttrap=0;
                            ftrap=prem_3;
                        }

                        pixo[j] += padd_2 + padd_3 - prem_3;
                    } /*replaces trap continue*/
                }/*end if j>0*/
            }/* end if qlevq > pmax, replaces continue*/

        }/*end for w*/
    //}
    return(status);
}
