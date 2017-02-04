#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "ctegen2.h"
//#include "pcte.h"
//#include "acs.h"

char MsgText[256];
extern void trlmessage (char *message);

int inverse_cte_blur(const SingleGroup *_rsz, SingleGroup *_rsc, SingleGroup *_fff, CTEParams *cte, int verbose, double expstart){

    SingleGroup __rsz;
    SingleGroup * rsz = &__rsz;
    initSingleGroup(rsz);
    allocSingleGroup(rsz, RAZ_COLS, RAZ_ROWS, False);

    assert(_rsz->sci.data.storageOrder == COLUMNMAJOR);
    assert(_rsz->dq.data.storageOrder == COLUMNMAJOR);

    assert(!copySingleGroup(rsz, _rsz, ROWMAJOR));

    SingleGroup __fff;
    SingleGroup * fff = &__fff;
    initSingleGroup(fff);
    allocSingleGroup(fff, RAZ_COLS, RAZ_ROWS, False);
    assert(_fff->sci.data.storageOrder == COLUMNMAJOR);
    copySingleGroup(fff, _fff, ROWMAJOR);

    extern int status;

    /*looping vars*/
    int NREDO, REDO;
    int NITINV, NITCTE;
    int i;
    int j,jj;
    double dmod;
    int jmax;
    float hardset=0.0f;
    int totflux=0;

    double cte_ff; /*cte scaling based on observation date*/
    double setdbl=0.0;

    /*DEFINE TO MAKE PRIVATE IN PARALLEL RUN*/
    double *pix_obsd=&setdbl;
    double *pix_modl=&setdbl;
    double *pix_curr=&setdbl;
    double *pix_init=&setdbl;
    double *pix_read=&setdbl;
    double *pix_ctef=&setdbl;

    /*STARTING DEFAULTS*/
    NITINV=1;
    NITCTE=1;
    cte_ff=0.0;
    jmax=0;
    dmod=0.0;

    /*LOCAL IMAGES TO PLAY WITH, THEY WILL REPLACE THE INPUTS*/
    SingleGroup rz; /*pixz_raz*/
    initSingleGroup(&rz);
    allocSingleGroup(&rz, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup rc; /*pixz_rac*/
    initSingleGroup(&rc);
    allocSingleGroup(&rc, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup pixz_fff; /*pixz_fff*/
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS, True);


    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
      APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
      INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
      CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
      OVER 3 YEARS AFTER INSTALLATION*/

    cte_ff=  (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);
    cte->scale_frac=cte_ff;   /*save to param structure for header update*/

    if(verbose){
        sprintf(MsgText,"CTE_FF (scaling fraction by date) = %g",cte_ff);
        trlmessage(MsgText);
    }

    /*SET UP THE SCALING ARRAY WITH INPUT DATA, hardset arrays for safety*/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0;j<RAZ_ROWS;j++){
            Pix(rc.sci.data,i,j)=hardset;
            Pix(rz.sci.data,i,j)=hardset;
            Pix(pixz_fff.sci.data,i,j)=hardset;
            Pix(rz.sci.data,i,j) = Pix(rsz->sci.data,i,j);
            Pix(rz.dq.data,i,j) = Pix(rsz->dq.data,i,j);
            Pix(pixz_fff.sci.data,i,j) =  Pix(fff->sci.data,i,j);//*cte_ff
        }
    }

    #pragma omp parallel for schedule (dynamic,1) \
    private(dmod,i,j,jj,jmax,REDO,NREDO,totflux, \
            pix_obsd,pix_modl,pix_curr,pix_init,\
            pix_read,pix_ctef,NITINV,NITCTE)\
    shared(rc,rz,cte,pixz_fff)

    for (i=0; i< RAZ_COLS; i++){
        pix_obsd = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_modl = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_curr = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_init = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_read = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_ctef = (double *) calloc(RAZ_ROWS, sizeof(double));

        totflux=0;
        /*HORIZONTAL PRE/POST SCAN POPULATION */
        for (j=0; j< RAZ_ROWS; j++){
            if(Pix(rz.dq.data,i,j)){
                pix_obsd[j] = Pix(rz.sci.data,i,j); /*starts as input RAZ*/
                totflux += 1;
            }
        }

        if (totflux >= 1) {/*make sure the column has flux in it*/
            NREDO=0; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            REDO=0; /*FALSE*/
            do { /*replacing goto 9999*/
                /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
                for (j=0; j<RAZ_ROWS; j++){
                    pix_modl[j] =  Pix(rz.sci.data,i,j);
                    pix_ctef[j] =  Pix(pixz_fff.sci.data,i,j);
                }
                /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
                  IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
                for (NITINV=1; NITINV<=cte->n_forward; NITINV++){
                    for (j=0; j<RAZ_ROWS; j++){
                        pix_curr[j]=pix_modl[j];
                        pix_read[j]=pix_modl[j];
                        pix_ctef[j]=Pix(pixz_fff.sci.data,i,j);
                    }

                    /*TAKE EACH PIXEL DOWN THE DETECTOR IN NCTENPAR=7*/
                    for (NITCTE=1; NITCTE<=cte->n_par; NITCTE++){
                        //simulatePixelReadout(pix_read, pix_ctef, cte, &cte->rprof->data, &cte->cprof->data, RAZ_ROWS);
                        sim_colreadout_l(pix_curr, pix_read, pix_ctef, cte);

                        /*COPY THE JUST UPDATED READ OUT IMAGE INTO THE INPUT IMAGE*/
                        for (j=0; j< RAZ_ROWS; j++){
                            pix_curr[j]=pix_read[j];
                        }
                    } /* end NITCTE */

                    /*DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE, THIS IS
                      AN ADDITIONAL AID IN MITIGATING THE IMPACT OF READNOISE*/
                    for (j=0; j< RAZ_ROWS; j++){
                        dmod =  (pix_obsd[j] - pix_read[j]);
                        if (NITINV < cte->n_forward){
                            dmod *= (dmod*dmod) /((dmod*dmod) + (cte->rn_amp * cte->rn_amp));
                        }
                        pix_modl[j] += dmod; /*dampen each pixel as the best is determined*/
                    }
                } /*NITINV end*/

                /*LOOK FOR AND DOWNSCALE THE CTE MODEL IF WE FIND
                  THE TELL-TALE SIGN OF READOUT CRS BEING OVERSUBTRACTED;
                  IF WE FIND ANY THEN GO BACK UP AND RERUN THIS COLUMN


                  THE WFC3 UVIS MODEL SEARCHES FOR OVERSUBTRACTED TRAILS.
                  WHICH ARE  DEFINED AS EITHER:

                  - A SINGLE PIXEL VALUE BELOW -10E-
                  - TWO CONSECUTIVE PIXELS TOTALING -12 E-
                  - THREE TOTALLING -15 E-

                  WHEN WE DETECT SUCH AN OVER-SUBTRACTED TAIL, WE ITERATIVELY REDUCE
                  THE LOCAL CTE SCALING BY 25% UNTIL THE TRAIL IS
                  NO LONGER NEGATIVE  THIS DOES NOT IDENTIFY ALL READOUT-CRS, BUT IT DOES
                  DEAL WITH MANY OF THEM. FOR IMAGES THAT HAVE BACKGROUND GREATER THAN 10 OR SO,
                  THIS WILL STILL END UP OVERSUBTRACTING CRS A BIT, SINCE WE ALLOW
                  THEIR TRAILS TO BE SUBTRACTED DOWN TO -10 RATHER THAN 0.

    */
                if (cte->fix_rocr) {
                    for (j=10; j< RAZ_ROWS-2; j++){
                        if (  (( cte->thresh > pix_modl[j] ) &&
                                    ( cte->thresh > (pix_modl[j] - pix_obsd[j]))) ||

                                (((pix_modl[j] + pix_modl[j+1]) < -12.) &&
                                 (pix_modl[j] + pix_modl[j+1] - pix_obsd[j] - pix_obsd[j+1] < -12.)) ||

                                (((pix_modl[j] + pix_modl[j+1] + pix_modl[j+2]) < -15.) &&
                                 ((pix_modl[j] + pix_modl[j+1] + pix_modl[j+2] -pix_obsd[j] -
                                   pix_obsd[j+1] - pix_obsd[j+2]) <-15.))  ){

                            jmax=j;

                            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
                            for (jj=j-10; jj<=j;jj++){
                                if ( (pix_modl[jj] - pix_obsd[jj]) >
                                        (pix_modl[jmax] - pix_obsd[jmax]) ) {
                                    jmax=jj;
                                }
                            }
                            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
                               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
                            for (jj=jmax; jj<=j;jj++){
                                Pix(pixz_fff.sci.data,i,jj) *= 0.75;
                            }
                            REDO=1; /*TRUE*/
                        } /*end if*/
                    } /*end for  j*/
                }/*end fix cr*/

                if (REDO) NREDO +=1;
                if (NREDO == 5)  REDO=0; /*stop*/
            } while (REDO); /*replacing goto 9999*/
        } /*totflux > 1, catch for subarrays*/

    #pragma omp critical (cte)
        for (j=0; j< RAZ_ROWS; j++){
            if (Pix(rz.dq.data,i,j)){
                Pix(rc.sci.data,i,j)= pix_modl[j];
            }
        }

        free(pix_obsd);
        free(pix_modl);
        free(pix_curr);
        free(pix_init);
        free(pix_read);
        free(pix_ctef);

    } /*end i*/

 /*   for (i=0; i< RAZ_COLS; i++){
        for (j=0; j< RAZ_ROWS; j++){
            if(Pix(rsz->dq.data,i,j)){
                Pix(rsz->sci.data,i,j) = Pix(rz.sci.data,i,j);
                //PixColumnMajor(_rsc->sci.data,j,i) = Pix(rc.sci.data,i,j);
                //Pix(rsc->sci.data,i,j) = Pix(rc.sci.data,i,j);
                Pix(fff->sci.data,i,j) = Pix(pixz_fff.sci.data,i,j);
            }
        }
    }
*/
    copySingleGroup(_rsc, &rc, COLUMNMAJOR);

    freeSingleGroup(&rz);
    freeSingleGroup(&rc);
    freeSingleGroup(&pixz_fff);

    //freeSingleGroup(rsc);
    freeSingleGroup(rsz);
    freeSingleGroup(fff);

    return(status);
}

int sim_colreadout_l(double *pixi, double *pixo, double *pixf, CTEParams *cte){


    extern int status;
    int j;
    int ttrap;

    int w;
    double ftrap;
    double pix_1;
    double padd_2;
    double padd_3;
    double prem_3;
    double pmax;
    double fcarry;

    padd_3=0.0;
    prem_3=0.0;
    padd_2=0.0;
    fcarry=0.0;
    pix_1=0.0;
    w=0;
    j=0;
    ftrap=0.0;
    ttrap=0;

    FloatHdrData *rprof;
    FloatHdrData *cprof;

    /*from the reference table*/
    rprof = cte->rprof;
    cprof = cte->cprof;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE  */
    memcpy(pixo, pixi, RAZ_ROWS*sizeof(*pixi));
    pmax=10.;
    for(j=0; j<RAZ_ROWS; j++){
        pmax = pixo[j] > pmax ? pixo[j] : pmax; //check assembly (before & after op) to see if actually implemented differently

        //pixo[j] = pixi[j];
       // if (pixo[j] > pmax)
         //   pmax=pixo[j];
    }


    unsigned maxChargeTrapIndex = cte->cte_traps-1;
        for (int w = maxChargeTrapIndex; w >= 0; --w)//go up or down? (if swap, change below condition)
        {
            if (cte->qlevq_data[w] <= pmax)//is any of this even needed or can we just directly map?
            {
                maxChargeTrapIndex = w;
                break;
            }
        }
    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    for (w = maxChargeTrapIndex; w>=0; w--){
        //if ( cte->qlevq_data[w] <= pmax )
        {

            ftrap = 0.0e0;
            ttrap = cte->cte_len; /*for referencing the image at 0*/
            fcarry = 0.0e0;

            /*GO UP THE COLUMN PIXEL BY PIXEL*/
            for(j=0; j<RAZ_ROWS;j++){
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
                        padd_2 = PixColumnMajor(rprof->data,ttrap-1,w) *ftrap;
                    }

                    padd_3 = 0.0;
                    prem_3 = 0.0;
                    if ( pix_1 >= cte->qlevq_data[w]){
                        prem_3 =  cte->dpdew_data[w] / cte->n_par * pixf[j];  /*dpdew is 1 in file */
                        if (ttrap < cte->cte_len)
                            padd_3 = PixColumnMajor(cprof->data,ttrap-1,w)*ftrap;
                        ttrap=0;
                        ftrap=prem_3;
                    }

                    pixo[j] += padd_2 + padd_3 - prem_3;
                } /*replaces trap continue*/
            }/*end if j>0*/
        }/* end if qlevq > pmax, replaces continue*/

    }/*end for w*/

    return(status);
}



int rsz2rsc(SingleGroup *pixz_fff, CTEParams *cte) {

extern int status;

int i,j;
double cte_i=0.0;
double cte_j=0.0;
double ro=0;
int io=0;
double ff_by_col[RAZ_COLS][4];
float hardset=0.0;

/*These are already in the parameter structure
  int     Ws              the number of traps < 999999, taken from pctetab read
  int     q_w[TRAPS];     the run of charge with level  cte->qlevq_data[]
  float   dpde_w[TRAPS];  the run of charge loss with level cte->dpdew_data[]

  float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel, TRAPS=999
  float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)

  The rprof array gives the fraction of charge that comes out of every parallel serial-shift
  the cummulative distribution in cprof then tells you what's left

*/

//SingleGroup pixz_fff;
//initSingleGroup(&pixz_fff);
//allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS);

/*SCALE BY 1 UNLESS THE PCTETAB SAYS OTHERWISE, I IS THE PACKET NUM
  THIS IS A SAFETY LOOP INCASE NOT ALL THE COLUMNS ARE POPULATED
  IN THE REFERENCE FILE*/

for(i=0; i<RAZ_COLS;i++){
    ff_by_col[i][0]=1.;
    ff_by_col[i][1]=1.;
    ff_by_col[i][2]=1.;
    ff_by_col[i][3]=1.;
    j= cte->iz_data[i]; /*which column to scale*/
   /* ff_by_col[j][0]=cte->scale512[i];
    ff_by_col[j][1]=cte->scale1024[i];
    ff_by_col[j][2]=cte->scale1536[i];
    ff_by_col[j][3]=cte->scale2048[i];
*/
    /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
      Index is figured on the final size of the image
      not the current size. Moved above
      */

    for(j=0; j<RAZ_ROWS; j++){
        Pix(pixz_fff->sci.data,i,j)=hardset;
        ro = j/512.0; /*ro can be zero, it's an index*/
        if (ro <0 ) ro=0.;
        if (ro > 2.999) ro=2.999; /*only 4 quads, 0 to 3*/
        io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
        cte_j= (j+1) / 2048.0;
        cte_i= ff_by_col[i][io] + (ff_by_col[i][io+1] -ff_by_col[i][io]) * (ro-io);
       Pix(pixz_fff->sci.data,i,j) =  (cte_i*cte_j);
    }
}

/*FOR REFERENCE TO JAYS CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

  int   iz_data[RAZ_ROWS];  column number in raz format
  double scale512[RAZ_ROWS];      scaling appropriate at row 512
  double scale1024[RAZ_ROWS];     scaling appropriate at row 1024
  double scale1536[RAZ_ROWS];     scaling appropriate at row 1536
  double scale2048[RAZ_ROWS];     scaling appropriate at row 2048
  */

/*THIS IS RAZ2RAC_PAR IN JAYS CODE - MAIN CORRECTION LOOP IN HERE*/
//inverse_cte_blur(rsz, rsc, &pixz_fff, cte, wf3->verbose,wf3->expstart);
//freeSingleGroup(&pixz_fff);
return(status);
}




int raz2rsz(const SingleGroup *raz, SingleGroup *rsz, double rnsig, int max_threads){
    /*
       This routine will read in a RAZ image and will output the smoothest
       image that is consistent with being the observed image plus readnoise. (RSZ image)
       This is necessary because we want the CTE-correction algorithm to produce the smoothest
       possible reconstruction, consistent with the original image and the
       known readnoise.  This algorithm constructs a model that is smooth
       where the pixel-to-pixel variations can be thought of as being related
       to readnoise, but if the variations are too large, then it respects
       the pixel values.  Basically... it uses a 2-sigma threshold.

       This is strategy #1 in a two-pronged strategy to mitigate the readnoise
       amplification.  Strategy #2 will be to not iterate when the deblurring
       is less than the readnoise.

*/


    extern int status;

    clock_t begin = clock();

    int i, j, NIT; /*loop variables*/
    int imid;
    double dptr=0.0;
    double  rms=0.0;
    double  rmsu=0.0;
    double nrms=0.0;
    double nrmsu=0.0;
    float hardset=0.0f;
    double setdbl=0.0;


    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING RAZ_COLS*/
    double obs_loc[3][RAZ_ROWS] ;
    double rsz_loc[3][RAZ_ROWS] ;

    NIT=1;

    /*ALL ELEMENTS TO FLAG*/
    for(i=0;i<3;i++){
        for (j=0; j<RAZ_ROWS; j++){
            obs_loc[i][j]=setdbl;
            rsz_loc[i][j]=setdbl;
        }
    }

    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, RAZ_COLS, RAZ_ROWS, True);


    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    for(i=0;i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS;j++){
            Pix(rsz->sci.data,i,j) = PixColumnMajor(raz->sci.data,j,i);
            Pix(rsz->dq.data,i,j) = PixColumnMajor(raz->dq.data,j, i);
            Pix(rnz.sci.data,i,j) = hardset;
            Pix(zadj.sci.data,i,j) = hardset;
        }
    }


    /*THE RSZ IMAGE JUST GETS UPDATED AS THE RAZ IMAGE IN THIS CASE*/
    if (rnsig < 0.1){
        trlmessage("rnsig < 0.1, No read-noise mitigation needed");
        return(status);
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
      */

    rms=setdbl;

    for(NIT=1; NIT<=100; NIT++){
        #pragma omp parallel for schedule(dynamic, 1) \
        private(i,j,imid,obs_loc,rsz_loc,dptr)\
        shared(raz, rsz, rnsig,rms,nrms, zadj)
        for(i=0; i<RAZ_COLS; i++){
            imid=i;
            /*RESET TO MIDDLE RAZ_COLS AT ENDPOINTS*/
            if (imid < 1)
                imid=1;
            if (imid == RAZ_COLS-1)
                imid = RAZ_COLS-2;

            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            for(j=0; j<RAZ_ROWS; j++){
                obs_loc[0][j] = PixColumnMajor(raz->sci.data,j, imid-1);
                obs_loc[1][j] = PixColumnMajor(raz->sci.data,j, imid);
                obs_loc[2][j] = PixColumnMajor(raz->sci.data,j, imid+1);

                rsz_loc[0][j] = Pix(rsz->sci.data,imid-1,j);
                rsz_loc[1][j] = Pix(rsz->sci.data,imid,j);
                rsz_loc[2][j] = Pix(rsz->sci.data,imid+1,j);
            }
            for (j=0; j<RAZ_ROWS; j++){
             if(PixColumnMajor(raz->dq.data,j, imid)) {
                 Pix(zadj.sci.data,i,j) = find_dadj_wf3(1+i-imid,j, obs_loc, rsz_loc, rnsig);
              }
            }
        } /*end the parallel for*/

        /*NOW GO OVER ALL THE RAZ_COLS AND RAZ_ROWS AGAIN TO SCALE THE PIXELS
        */
        for(i=0; i<RAZ_COLS;i++){
            for(j=0; j<RAZ_ROWS; j++){
                if (PixColumnMajor(raz->dq.data,j, i)){
                    Pix(rsz->sci.data,i,j) +=  (Pix(zadj.sci.data,i,j)*0.75);
                    Pix(rnz.sci.data,i,j) = (PixColumnMajor(raz->sci.data,j,i) - Pix(rsz->sci.data,i,j));
                }
            }
        }

        rms=setdbl;
        nrms=setdbl;

        /*This is probably a time sink because the arrays are being
          accessed out of storage order, careful of page faults */
        #pragma omp parallel for schedule(dynamic,1)\
        private(i,j,rmsu,nrmsu) \
        shared(raz,rsz,rms,rnsig,nrms)
        for(j=0; j<RAZ_ROWS; j++){
            nrmsu=setdbl;
            rmsu=setdbl;
            for(i = 0;i<RAZ_COLS; i++){
                if ( (fabs(PixColumnMajor(raz->sci.data,j,i)) > 0.1) ||
                        (fabs(Pix(rsz->sci.data,i,j)) > 0.1) ){
                    rmsu  +=  ( Pix(rnz.sci.data,i,j) * Pix(rnz.sci.data,i,j) );
                    nrmsu += 1.0;
                }
            }
            #pragma omp critical (rms)
            {   rms  += rmsu;
                nrms += nrmsu;
            }
        }
        rms = sqrt(rms/nrms);

        /*epsilon type comparison*/
        if ( (rnsig-rms) < 0.00001){
            printf("NIT - %d\n", NIT);

            break; /*this exits the NIT for loop*/
        }
    } /*end NIT*/

    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);

    if (1)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        unsigned maxThreads = 12;
        sprintf(MsgText,"Time taken to smooth image: %.2f(s) with %i procs/threads\n",timeSpent/maxThreads,maxThreads);
        trlmessage(MsgText);
    }
    return (status);
}

double find_dadj_wf3(int i ,int j, double obsloc[][RAZ_ROWS], double rszloc[][RAZ_ROWS], double rnsig)
{

    /*
       This function determines for a given pixel how it can
       adjust in a way that is not inconsistent with its being
       readnoise.  To do this, it looks at its upper and lower
       neighbors and sees whether it is consistent with either
       (modulo readnoise).  To the extent that it is consistent
       then move it towards them.  But also bear in mind that
       that we don't want it to be more than 2 RN sigmas away
       from its original value.  This is pretty much a tug of
       war... with readnoise considerations pushing pixels to
       be closer to their neighbors, but the original pixel
       values also pull to keep the pixel where it was.  Some
       accommodation is made for both considerations.
       */

    extern int status;

    double mval=0.0;
    double    dval0, dval0u, w0;
    double    dval9u, w9;
    double    dmod1, dmod1u, w1;
    double    dmod2, dmod2u, w2;

    dval0=0.;
    dval0u=0.;
    w0=0.;
    volatile double dval9=0.;
    dval9u=0.;
    w9=0.;
    dmod1=0.;
    dmod1u=0.;
    w1=0.;
    dmod2=0.;
    dmod2u=0.;
    w2=0.;

    mval = rszloc[i][j];
    dval0  = obsloc[i][j] - mval;
    dval0u = dval0;

    if (dval0u >1.0)
        dval0u =  1.0;
    if (dval0u <-1.0)
        dval0u = -1.0;

    dval9 = 0.;

    /*COMPARE THE SURROUNDING PIXELS*/
    if (i==1 &&  RAZ_ROWS-1>j  && j>0 ) {

        dval9 = obsloc[i][j-1]  - rszloc[i][j-1] +
            obsloc[i][j]    - rszloc[i][j]  +
            obsloc[i][j+1]  - rszloc[i][j+1] +
            obsloc[i-1][j-1]- rszloc[i-1][j-1] +
            obsloc[i-1][j]  - rszloc[i-1][j] +
            obsloc[i-1][j+1]- rszloc[i-1][j+1] +
            obsloc[i+1][j-1]- rszloc[i+1][j-1] +
            obsloc[i+1][j]  - rszloc[i+1][j] +
            obsloc[i+1][j+1]- rszloc[i+1][j+1];
        dval9 =dval9 / 9.;

    }
    dval9u = dval9;

    if (dval9u > (rnsig*0.33))
        dval9u =  rnsig*0.33;
    if (dval9u <  rnsig*-0.33)
        dval9u = rnsig*-0.33;

    dmod1 = 0.;
    if (j>0)
        dmod1 = rszloc[i][j-1] - mval;

    dmod1u = dmod1;
    if (dmod1u > rnsig*0.33)
        dmod1u =  rnsig*0.33;
    if (dmod1u < rnsig*-0.33)
        dmod1u = rnsig*-0.33;

    dmod2 = 0.;
    if (j < RAZ_ROWS-1)
        dmod2 =  rszloc[i][j+1] - mval;

    dmod2u = dmod2;
    if (dmod2u > rnsig*0.33)
        dmod2u =  rnsig*0.33;
    if (dmod2u < rnsig*-0.33)
        dmod2u = rnsig*-0.33;


    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    w0 =   (dval0*dval0) / ((dval0*dval0)+ 4.0*(rnsig*rnsig));
    w9 =   (dval9*dval9) / ((dval9*dval9)+ 18.0*(rnsig*rnsig));
    w1 = (4*rnsig*rnsig) / ((dmod1*dmod1)+4.0*(rnsig*rnsig));
    w2 = (4*rnsig*rnsig) / ((dmod2*dmod2)+4.0*(rnsig*rnsig));

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    return ((dval0u * w0 * 0.25f) + /* desire to keep the original pixel value */
            (dval9u*w9*0.25f) + /* desire to keep the original sum over 3x3*/
            (dmod1u*w1*0.25f) + /*desire to get closer to the pixel below*/
            (dmod2u*w2*0.25f)) ; /*desire to get closer to the pixel above*/

    //return(status);
}
