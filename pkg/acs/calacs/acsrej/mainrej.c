/* mainrej -- Main program for controlling ACSREJ */

# include <stdio.h>
# include <stdlib.h>        /* calloc */
# include <string.h>

# include <c_iraf.h>        /* for c_irafinit */

#include "hstcal.h"
# include "acs.h"
# include "hstcalerr.h"
# include "acsrej.h"
# include "hstcalversion.h"
#include "trlbuf.h"

static void FreeNames (char *);

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'
struct TrlBuf trlbuf = { 0 };

int main (int argc, char **argv) {

    int status = HSTCAL_OK;
    char    *input, output[CHAR_LINE_LENGTH];    /* file names */
    clpar   par;                                    /* parameters used */
    int     newpar[MAX_PAR+1];          /* user specifiable parameters */
    char    mtype[SZ_STRKWVAL+1];      /* Role of exposure in association */

    int rej_command (int, char **, char **, char *, clpar *, int []);
    int AcsRej (char *, char *, char *, clpar *, int []);
    void WhichError (int);

    c_irafinit (argc, argv);

    /* Initialize mtype to NULL to signal no change in ASN_MTYP for output*/
    mtype[0] = '\0';

    /* Allocate later */
    input = NULL;

    /* Get input and output file names and switches in the command line. */
    if (rej_command (argc, argv, &input, output, &par, newpar)) {
        FreeNames(input);
        exit (ERROR_RETURN);
    }

    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();
    trlGitInfo();

    /* Reject cosmic rays. */
    if (AcsRej (input, output, mtype, &par, newpar)) {
        sprintf (MsgText,"Error processing %s.", input);
        asnerror (MsgText);
    }

    if (status) {
        WhichError (status);
        FreeNames(input);
        CloseTrlBuf ();
        exit (ERROR_RETURN);
    } else {
        FreeNames(input);
        CloseTrlBuf ();
        exit (0);
    }
}


static void FreeNames (char *input) {
    if (input != NULL)
        free(input);
}
