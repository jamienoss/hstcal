# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "err.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

/* These routines facilitate moving between the regular WFC3 image structre
   to the RAZ image structure used in the CTE correction and Sink pixel flagging


   convert a raw file to raz file: CDAB longwise amps, save data array
   for comparison with what jay has during testing

   In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left.
   The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan),
   with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down,
   then serial-shifted to the left.

   The code for the DQ arrays and plain float arrays only converts one chip at a time so that it can be run through the regular
   wf3ccd pipeline which operates on 1 chip at a time.

   Megan Sosey, May 2015

*/



/*convert the sci and dq extensions to the long format*/
int makeRAZ(const SingleGroup *cd, const SingleGroup *ab, SingleGroup *raz)
{
    /*

         convert a raw file to raz file: CDAB longwise amps, save data array
         for comparison with what jay has during testing

         -->do an additional bias correction using the  residual bias level measured for each amplifier from the
         steadiest pixels in the horizontal overscan and subtracted fom the pixels for that amplifier.

         ---> convert into electrons at the end
         ---> add supplemental bias info to the header

         allocate contiguous 2d array on the heap
         with pointers and return the pointer to the head of the array

         The Following macros are used to represent 2-d indexing.
         Two dimensional arrays are stored in FITS order.

         ny
         ^
         N | a05   a15   a25   a35
         A | a04   a14   a24   a34
         X | a03   a13   a23   a33
         I | a02   a12   a22   a32
         S | a01   a11   a21   a31
         2 | a00   a10   a20   a30
         ---------------------------> nx
         NAXIS1

         NAXIS1 is 4 and NAXIS2 is 6
         PIX(a,1,4) accesses a14

         In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left.
         The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan),
         with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down,
         then serial-shifted to the left.

  */
    extern int status;

    const unsigned nColumns = raz->sci.data.nx;
    const unsigned nRows = raz->sci.data.ny;
    const unsigned nColumnsPerChip = nColumns / 4; /* for looping over quads  */

#ifdef _OPENMP
    const unsigned nThreads = omp_get_num_procs();
    #pragma omp parallel for num_threads(nThreads) shared(raz, cd, ab) schedule(static)
#endif
	for (unsigned i = 0; i < nRows; ++i)
	{
		for (unsigned j = 0; j < nColumnsPerChip; ++j)
		{
		    //sci data
			Pix(raz->sci.data, j, i) = Pix(cd->sci.data, j, i);
			Pix(raz->sci.data, j+nColumnsPerChip, i) = Pix(cd->sci.data, nColumnsPerChip*2-j-1, i);
			Pix(raz->sci.data, j+2*nColumnsPerChip, i) = Pix(ab->sci.data, j, nRows-i-1);
			Pix(raz->sci.data, j+3*nColumnsPerChip, i) = Pix(ab->sci.data, nColumnsPerChip*2-j-1, nRows-i-1);

			//dq data
			Pix(raz->dq.data, j, i) = Pix(cd->dq.data, j, i);
			Pix(raz->dq.data, j+nColumnsPerChip, i) = Pix(cd->dq.data, nColumnsPerChip*2-j-1, i);
			Pix(raz->dq.data, j+2*nColumnsPerChip, i) = Pix(ab->dq.data, j, nRows-i-1);
			Pix(raz->dq.data, j+3*nColumnsPerChip, i) = Pix(ab->dq.data, nColumnsPerChip*2-j-1, nRows-i-1);
		}
	}

    return status;
}


/* Transform a RAZ format image back into the separate input arrays calwf3 likes*/
int undoRAZ(SingleGroup *cd, SingleGroup *ab, const SingleGroup *raz)
{
    extern int status;

    const unsigned nColumns = raz->sci.data.nx;
    const unsigned nRows = raz->sci.data.ny;
    const unsigned nColumnsPerChip = nColumns / 4; /* for looping over quads  */

    /*REVERSE THE AMPS TO THE RAW FORMAT*/
   #ifdef _OPENMP
       const unsigned nThreads = omp_get_num_procs();
       #pragma omp parallel for num_threads(nThreads) shared(raz, cd, ab) schedule(static)
   #endif
    for (unsigned j = 0; j < nRows; ++j)
    {
        for (unsigned i = 0; i < nColumnsPerChip; ++i)
        {
             Pix(cd->sci.data, i, j) = Pix(raz->sci.data, i, j);
             Pix(cd->sci.data, nColumnsPerChip*2-i-1, j) = Pix(raz->sci.data, i+nColumnsPerChip, j);
             Pix(ab->sci.data, i, nRows-j-1) = Pix(raz->sci.data, i+2*nColumnsPerChip, j);
             Pix(ab->sci.data, nColumnsPerChip*2-i-1, nRows-j-1) = Pix(raz->sci.data, i+3*nColumnsPerChip, j);

             Pix(cd->dq.data, i, j) = Pix(raz->dq.data, i, j);
             Pix(cd->dq.data, nColumnsPerChip*2-i-1, j) = Pix(raz->dq.data, i+nColumnsPerChip, j);
             Pix(ab->dq.data, i, nRows-j-1) = Pix(raz->dq.data, i+2*nColumnsPerChip, j);
             Pix(ab->dq.data, nColumnsPerChip*2-i-1, nRows-j-1) = Pix(raz->dq.data, i+3*nColumnsPerChip, j);
        }
    }

    return status;
}


/***********THE ROUTINES BELOW OPERATE ON A SINGLE GROUP***********/

/*
The DQ versions here are written to work with SINKDETECT
and work on 1 chip at a time, with Half the Columns of the
full size raz image and all the rows
*/

int makedqRAZ(SingleGroup *x, SingleGroup *raz){
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                Pix(raz->dq.data,i,j) = Pix(x->dq.data,i,j);
                Pix(raz->dq.data,i+subcol,j) = Pix(x->dq.data,subcol*2-i-1,j);
            }
        }

    } else {
        if (x->group_num == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->dq.data,i,j) = Pix(x->dq.data,i,RAZ_ROWS-j-1);
                    Pix(raz->dq.data,i+subcol,j) = Pix(x->dq.data,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }

        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

   return(status);

}


/* Transform dq in a  RAZ format image back into the separate input arrays calwf3 likes*/
int undodqRAZ(SingleGroup *x, SingleGroup *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i< subcol; i++){
            for (j=0; j<RAZ_ROWS; j++){
                 Pix(x->dq.data,i,j) = Pix(raz->dq.data,i,j);
                 Pix(x->dq.data,2*subcol-i-1,j) = Pix(raz->dq.data,i+subcol,j);
            }
        }
    } else {
        if (x->group_num == 2){
            for (i=0; i< subcol; i++){
                for (j=0; j<RAZ_ROWS; j++){
                     Pix(x->dq.data,i,RAZ_ROWS-j-1) = Pix(raz->dq.data,i,j);
                     Pix(x->dq.data,subcol*2-i-1,RAZ_ROWS-j-1) = Pix(raz->dq.data,i+subcol,j);
                }
            }
        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

    return (status);

}

/*convert the science image of a single group to RAZ
 for use in the SINK pixel detection*/

int makeSciSingleRAZ(SingleGroup *x, SingleGroup *raz){
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                Pix(raz->sci.data,i,j) = Pix(x->sci.data,i,j);
                Pix(raz->sci.data,i+subcol,j) = Pix(x->sci.data,subcol*2-i-1,j);
            }
        }

    } else {
        if (x->group_num == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->sci.data,i,j) = Pix(x->sci.data,i,RAZ_ROWS-j-1);
                    Pix(raz->sci.data,i+subcol,j) = Pix(x->sci.data,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }

        } else {
            trlmessage("Invalid group number passed to makeSciSingleRAZ");
            return(status=INVALID_VALUE);
        }
    }

   return(status);

}

/*convert  floating point arrays into raz format*/
int makeFloatRaz(FloatTwoDArray *x, FloatTwoDArray  *raz, int group){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (group == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                PPix(raz,i,j) = PPix(x,i,j);
                PPix(raz,i+subcol,j) = PPix(x,subcol*2-i-1,j);
            }
        }
    } else {
        if (group == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    PPix(raz,i,j) = PPix(x,i,RAZ_ROWS-j-1);
                    PPix(raz,i+subcol,j) =PPix(x,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }
        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

    return(status);
}
