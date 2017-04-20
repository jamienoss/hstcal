# include <math.h>
# include <stdio.h>
# include "hstio.h"
# include "err.h"	/* SIZE_MISMATCH */
# include "wf3.h"
# include "assert.h"

/* Subtract the second SingleGroupLine from the first, leaving the
   result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are combined;
   the data quality arrays are combined.
	
    
   M. Sosey, 2013 Sept 09
   Added new routine to deal with uvis subarrays  for postflash correction
   which knows what to do with different reference file / science sizes and the
   uvis overscan region 
*/

int sub1d (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a		io: input data; output difference
int line		 i: line of input data to subtract 1-d data from
SingleGroupLine *b	 i: second input data
*/

    extern int status;

    // These asserts should really be null checks & returns.
    // However, that requires a more substantial refactoring and investigating
    assert(a);
    assert(b);

    if (a->sci.data.nx != b->sci.tot_nx)
        return (status = SIZE_MISMATCH);

    /* science, error, and DQ data */
    unsigned dimx = a->sci.data.nx;

    /* science array */
    if (a->sci.data.data)
    {
        for (unsigned i = 0;  i < dimx;  ++i)
        {
            Pix(a->sci.data, i, line) =
                    Pix (a->sci.data, i, line) - b->sci.line[i];
        }
    }

    /* error array */
    if (a->err.data.data)
    {
        for (unsigned i = 0;  i < dimx;  ++i)
        {
            float da = Pix (a->err.data, i, line);
            float db = b->err.line[i];
            Pix (a->err.data, i, line) = sqrt (da * da + db * db);
        }
    }

    /* data quality */
    if (a->dq.data.data)
    {
        for (unsigned i = 0;  i < dimx;  ++i)
        {
             short dqa = DQPix (a->dq.data, i, line);
             short dqb = b->dq.line[i];
             short dqab = dqa | dqb;
             DQSetPix (a->dq.data, i, line, dqab);
        }
    }

    return (status);
}

int sub1dreform (SingleGroup *a, int line, int overstart, SingleGroupLine *b) {

/* arguments:
SingleGroup *a		io: input data; output difference
int line		 i: line of input data to subtract 1-d data from
SingleGroupLine *b	 i: second input data
int overstart : this is where the overscan starts in the cut image

This is the same as sub1d, except that it is meant to be used with the postflash
correction for subarrays only. This is because the reference file is always a full frame 4-amp readout.
GO users are restricted to subarrays which are entirely contained in one amp (C amp I think). Group
members can place subarrays anywhere for calibration programs and have cases where the subrarray bisects
the virtual overscan region in the postflash reference file. This means that just that striped section needs to
be removed from the reference image and a new stitched line created for subtraction

So in this case, the second single group line is larger than the first and the virtual overscan
strip is removed from the data so that they are the same size before the subtraction

In this case the size of *a does NOT match the size of *b
*/

    extern int status;

    // These asserts should really be null checks & returns.
    // However, that requires a more substantial refactoring and investigating
    assert(a);
    assert(b);

    unsigned sizea=a->sci.data.nx;

    /* science array */
    if (a->sci.data.data)
    {
        for (unsigned i=0, j=0;  i < sizea;  ++i, ++j)
        {
            if (i == overstart)
                j += 60;
            Pix(a->sci.data, i, line) = Pix (a->sci.data, i, line) - b->sci.line[j];
        }
    }

    /* error array */
    if (a->err.data.data)
    {
        for (unsigned i=0, j=0;  i < sizea;  ++i, ++j)
        {
            if (i == overstart)
                j += 60;
             float da = Pix (a->err.data, i, line);
             float db = b->err.line[j];
             Pix (a->err.data, i, line) = sqrt (da * da + db * db);
        }
    }

    /* data quality */
    if (a->dq.data.data)
    {
        for (unsigned i=0, j=0;  i < sizea;  ++i, ++j)
        {
            if (i == overstart)
                j += 60;
             short dqa = DQPix (a->dq.data, i, line);
             short dqb = b->dq.line[j];
             short dqab = dqa | dqb;
             DQSetPix (a->dq.data, i, line, dqab);
        }
    }

    return (status);
}

