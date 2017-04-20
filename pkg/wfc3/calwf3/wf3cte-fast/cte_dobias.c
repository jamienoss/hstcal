# include <stdio.h>
# include <string.h>
# include <assert.h>
# include <math.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"
# include "cte-fast.h"

int doNewCTEBias( SingleGroup * image, char * filename, CTEParamsFast * ctePars)
{
    //WARNING - assumes row major storage
    assert(image->sci.data.storageOrder == ROWMAJOR);

    extern int status;

    sprintf(MsgText,"CTE: Subtracting BIACFILE: %s for imset %d", filename, image->group_num);
    trlmessage(MsgText);

    if (!ctePars->refAndIamgeBinsIdenticle)
    {
       sprintf (MsgText, "BIAC image and input are not binned to the same pixel size!");
       trlerror (MsgText);
       return (status = SIZE_MISMATCH);
    }

    if (ctePars->verbose)
    {
        sprintf(MsgText,"Image has starting location of %d,%d in the reference image", ctePars->columnOffset, ctePars->rowOffset);
        trlmessage(MsgText);
    }

    SingleGroupLine biacLine;
    initSingleGroupLine(&biacLine);
    openSingleGroupLine (filename, image->group_num, &biacLine);
    if (hstio_err())
        return (status = OPEN_FAILED);

    //used to vary for dev purposes
    unsigned rowsStart = 0;//ctePars->imageRowsStart;
    unsigned rowsEnd = ctePars->nRows;//ctePars->imageRowsEnd;
    unsigned columnsStart[2] = {ctePars->imageColumnsStart[0], ctePars->imageColumnsStart[1]};
    unsigned columnsEnd[2] = {ctePars->imageColumnsEnd[0], ctePars->imageColumnsEnd[1]};

    for (unsigned i = rowsStart; i < rowsEnd; ++i)
    {
        status = getSingleGroupLine(filename, i, &biacLine);
        if (status)
        {
            sprintf(MsgText,"Could not read line %d from bias image.", i+1);
            trlerror(MsgText);
            closeSingleGroupLine(&biacLine);
            freeSingleGroupLine(&biacLine);
            return status;
        }
        for (unsigned amp = 0; amp < 2; ++amp)
        {
            for (unsigned j = columnsStart[amp]; j < columnsEnd[amp]; ++j)
            {
                //image
                Pix(image->sci.data, j, i) -= biacLine.sci.line[ctePars->columnOffset + j];
                //err
                float biacErrLine = biacLine.err.line[ctePars->columnOffset + j];
                float imageErrLine = Pix(image->err.data, j, i);
                Pix(image->err.data, j, i) = sqrt(imageErrLine * imageErrLine + biacErrLine * biacErrLine);
                //dq
                Pix(image->dq.data, j, i) = Pix(image->dq.data, j, i) | biacLine.dq.line[ctePars->columnOffset + j];
            }
        }
    }

    closeSingleGroupLine(&biacLine);
    freeSingleGroupLine(&biacLine);
    return status;
}
