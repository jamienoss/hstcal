#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ximio.h"
#include "trlbuf.h"
#include "hstcalerr.h"

const unsigned initLength = 2;

/* Internal trailer file buffer routines */
static void ResetTrlBuf (void);
static void AddTrlBuf (char *message);
static void WriteTrlBuf (char *message);
static void CatTrlFile (FILE *ip, FILE *op);
static void CatTrlFile_NoEOF (FILE *ip, FILE *op);
static int AppendTrlFile();

struct TrlBuf trlbuf = { 0 };

int InitTrlFile (char *inlist, char *output)
{
    /*
        This initialization function sets up the trailer file to be used
        as output for all the messages from the task it is called from.
        It will concatenate multiple input trailer files
        or copy intial portion of input trailer file into one output
        file, creating the single trailer file for the messages to be
        written to...
    */

/* Parameters:
char *inlist        i: list of input trailer filenames
char *output        i: full filename of output (final) trailer file
*/
    extern int status;

    IRAFPointer tpin;
    FILE *ip, *tp;
    int n, td;
    char trldata[SZ_LINE+1];
    static char uniq_outtemplate[] = "tmp_calacs_XXXXXX";
    char uniq_outname[SZ_LINE+1];

    void SetTrlOverwriteMode (int);
    int unlink(const char *);

    trldata[0] = '\0';

    /* Copy name of output file to trlbuf */
    strcpy (trlbuf.trlfile, output);

    /*  Open the output file, for concatenating
        exposure trailers or sub-product trailers, or simply for
        appending new comments to end of existing trailer file.
    */
    if ( (trlbuf.fp = fopen(trlbuf.trlfile,"a+")) == NULL) {
        trlopenerr(trlbuf.trlfile);
        return(status=INVALID_TEMP_FILE);
    }


    /* open the input file template */
    tpin = c_imtopen (inlist);

    /* Only if we are building a product/sub-product trailer
        file, do we need to concatenate input trailer files
        or append to an input file.
    */
    if (strcmp(inlist, output) != 0) {
        /* Generate temporary output file name.
           This is to avoid infinite loop when output is
           accidentally the same as one of the inputs.
           Not using tmpfile() because it opens binary stream.
           Not using tmpnam() because compiler complains about danger.
        */
        strcpy(uniq_outname, uniq_outtemplate);
        if ( (td = mkstemp(uniq_outname)) < 0 ) {
            trlerror("Failed to create temporary file name.");
            return(status=INVALID_TEMP_FILE);
        }
        if ( unlink(uniq_outname) < 0) {
            trlerror("Failed to unlink temporary file name.");
            return(status=INVALID_TEMP_FILE);
        }
        if ( (tp = fdopen(td, "w+") ) == NULL ) {
            trlopenerr(uniq_outname);
            return(status=INVALID_TEMP_FILE);
        }

        /* Since we have determined that we are indeed combining
            input trailer files to create a new trailer file, ...
        */
        /* loop all input files */
        for (n = 0; n < c_imtlen(tpin); n++) {

            /* read the next input image name in the template list */
            c_imtgetim (tpin, trldata, SZ_FNAME);

            /* open the file (read-only) and add to temp output file... */
            if ((ip = fopen (trldata, "r")) == NULL) {
                /* Report the error, but the processing can continue anyway */
                trlopenerr(trldata);
            }

            /* Do we have an input file to start with? */
            if (ip != NULL) {
                /* If so, append the input to output */
                CatTrlFile(ip, tp);
                /* Done with this input file... */
                fclose(ip);
            }

        }    /* End loop over all input files */

        fflush(tp);

        /* Copy temporary file content to output trailer file */
        CatTrlFile_NoEOF(tp, trlbuf.fp);
        fclose(tp);  /* Also delete the file because of unlink() */

        /* Reset overwrite now to NO */
            SetTrlOverwriteMode(NO);

        /* End if (strcmp) */
    } else if (trlbuf.overwrite == YES) {

        /* We are revising/creating a single trailer file
           We want to overwrite old CALACS comments with new,
           so set the file pointer to just before the first line
           which starts with TRL_PREFIX...
        */
        if (AppendTrlFile ())
            return(status);

        /* We are finished overwriting old comments, so reset this mode. */
        SetTrlOverwriteMode(NO);

    }

    c_imtclose(tpin);
    trlbuf.init = 1;

    return (status);
}
static void CatTrlFile(FILE *ip, FILE *op)
{
    /* This function appends the ENTIRE contents of the input file (ip)
        to those of the output file (op)...
    */

    char buf[SZ_LINE+1];

    buf[0] = '\0';

    /* We need to insure that we start at the beginning of
        the input file (ip)...
    */
    rewind(ip);

    /* Now copy the file into the output trailer file */
    while ( !feof(ip) ) {
        fgets(buf, SZ_LINE, ip);
        fprintf(op,"%s",buf);
    }

    /* Flush the pipe to make sure that nothing gets left behind
        in whatever cache may be in operation.
    */
    fflush (op);
}
static void CatTrlFile_NoEOF(FILE *ip, FILE *op)
{
    // Like CatTrlFile() but ip does not have EOF

    char buf[SZ_LINE+1];

    buf[0] = '\0';

    /* We need to insure that we start at the beginning of
        the input file (ip)...
    */
    rewind(ip);

    /* Now copy the file into the output trailer file */
    while ( fgets(buf, SZ_LINE+1, ip) != NULL ) {
        fprintf(op,"%s",buf);
    }

    /* Flush the pipe to make sure that nothing gets left behind
        in whatever cache may be in operation.
    */
    fflush (op);
}
static int AppendTrlFile()
{
    /* This function sets the file pointer for the input file (ip)
        up to the line just before encountering TRL_PREFIX ...
        It reads in all the comments prior to TRL_PREFIX into memory,
        then opens the trailer file in overwrite('w+' instead of 'a+') mode
        and prints those comments to it.
        WJH 14 Apr 2000

        18 Jan 2002 (WJH) & 21 June 2002 (HAB): Removed extraneous, bug-ridden
        'tmpptr' buffer used during reallocation of oprefix buffer size.
    */

    extern int status;
    char buf[SZ_LINE+1];

    char *oprefix;


    if ( (oprefix = realloc (NULL, initLength)) == NULL){
        trlerror ("Out of memory for trailer file preface.");
        return (status = OUT_OF_MEMORY);
    }
    buf[0] = '\0';
    oprefix[0] = '\0';

    /* Make sure we start searching from the beginning of the file */
    rewind (trlbuf.fp);

    while ( !feof(trlbuf.fp) ) {
        /* Read in a line */
        fgets(buf, SZ_LINE+1, trlbuf.fp);

        /* If we find the prefix, stop searching */
        if (strstr(buf,TRL_PREFIX) !=NULL) {
            break;
        } else {
            /* Store this line in a buffer to be written out when
                the old file is overwritten...
            */
            oprefix = realloc (oprefix, (strlen(oprefix) + strlen(buf) +2 ));

            if (oprefix == NULL) {
                asnmessage ("Out of memory: Couldn't store trailer file comment.");
                status = OUT_OF_MEMORY;
                /* Clean-up... */
                free (trlbuf.buffer);
                free (trlbuf.preface);
                free (oprefix);
                fclose (trlbuf.fp);
                return (status);
            } else {
                strcat (oprefix, buf);
            }
        }
    }
    /* Now we know what needs to be kept, let's close the file... */
    fclose (trlbuf.fp);
    trlbuf.fp = NULL;

    /* ...reopen it with 'w+' instead of 'a+'... */
    if ( (trlbuf.fp = fopen(trlbuf.trlfile,"w+")) == NULL) {
        trlopenerr(trlbuf.trlfile);
        return(status=INVALID_TEMP_FILE);
    }
    /* ... and write out the preface information we wanted to keep. */
    fprintf (trlbuf.fp, "%s\n", oprefix);
    fflush (trlbuf.fp);

    /* Now, clean up the memory used by the temp buffer. */
    oprefix[0] = '\0';
    free (oprefix);

    return (status);
}
int WriteTrlFile (void)
{
    /*
        This function closes the trailer file.
    */

    extern int status;

    /* Now that we have copied the information to the final
        trailer file, we can close it and the temp file...
    */
    if (trlbuf.fp != NULL)
        fclose (trlbuf.fp);
    /* Reset value of fp to NULL so that later checks will be valid */
    trlbuf.fp = NULL;

    return (status);
}
int InitTrlBuf (void)
{
    /* This initialization routine must be called before any others in this
        file.
    */
    extern int status;

    trlbuf.trlfile[0] = '\0';
    trlbuf.fp = NULL;
    trlbuf.overwrite = NO;  /* Initialize with default of append */
    trlbuf.quiet = NO;      /* Initialize to produce STDOUT messages */
    trlbuf.usepref = YES;      /* Switch to specify whether to output preface */

    if ( (trlbuf.buffer = realloc (NULL, initLength)) == NULL){
        trlerror ("Out of memory for trailer file buffer.");
        return (status = OUT_OF_MEMORY);
    }
    if ( (trlbuf.preface = realloc (NULL, initLength)) == NULL){
        trlerror ("Out of memory for trailer file preface.");
        return (status = OUT_OF_MEMORY);
    }

    trlbuf.buffer[0] = '\0';
    trlbuf.preface[0] = '\0';
    trlbuf.init = 1;
    return(status);
}
void SetTrlPrefaceMode (int use)
{
    /*
        This function sets the switch to specify whether the preface should
        be output to the trailer file or not.  This is used when no CRREJ is
        done on a list of images, and you don't want the preface in the trailer
        files of the sub-products again.
    */

    trlbuf.usepref = use;

}
void SetTrlOverwriteMode (int owrite)
{
    /*
        This function sets the overwrite switch.  Once set, all comments
        after the first CALXXXBEG string will be overwritten.
        This will be used strictly for the RAW file trailer files,
        where the generic conversion comments should be kept, but
        the CALXXX comments should be overwritten with the latest results.
    */

    trlbuf.overwrite = owrite;
}
void SetTrlQuietMode (int quiet)
{
    /* This function records the value of the command-line parameter QUIET
        into the structure...
    */
    trlbuf.quiet = quiet;
}
static void AddTrlBuf (char *message)
{
    /* Add a new message line to the buffer.
        Re-allocate space for the buffer and append line to new buffer.
    */

    /* arguments:
    char *message         i: new trailer file line to add to buffer
    */
    void asnmessage (char *);

    extern int status;

    if ( ! trlbuf.init )
        assert(0); //TRLBUF NOT INIT, YOU MAY HAVE PROBLEMS

    trlbuf.buffer = realloc (trlbuf.buffer,
                 (strlen(trlbuf.buffer) + strlen(message) +2));

    if (trlbuf.preface == NULL) {
        asnmessage ("Out of memory: Couldn't store trailer file comment.");
        status = OUT_OF_MEMORY;
    } else {

        strcat (trlbuf.buffer, message);

        /* Append a newline at the end of every output message */
        strcat (trlbuf.buffer, "\n");
    }
}
void InitTrlPreface (void)
{
    /*
        This function will copy contents of the buffer into the preface
    */

    void asnmessage (char *);
    extern int status;

    trlbuf.preface = realloc (trlbuf.preface, (strlen(trlbuf.buffer) +2));
    if (trlbuf.preface == NULL) {
        asnmessage ("Out of memory: Couldn't store trailer file preface.");
        status = OUT_OF_MEMORY;
    } else {
        strcpy (trlbuf.preface, trlbuf.buffer);
    }
}
void ResetTrlPreface (void)
{
    /* This function will reset the preface to a null string, so nothing
        more gets written out to any other trailer files.  It also
        reallocates the space to preface so that it doesn't take any memory.
    */

    /* Start by freeing the initial space used, and setting pointer
        to NULL
    */
    free (trlbuf.preface);
    trlbuf.preface = NULL;

    /* Now, allocate new pointer as before */
    trlbuf.preface = realloc (NULL, initLength);

    /* ...and initialize it to NULL */
    trlbuf.preface[0] = '\0';
}
static void ResetTrlBuf (void)
{
    /*
        Clear out messages from buffer.  Intended for use after
        messages were already written to file.  Assumes trlbuf.buffer
        already exists!
    */

    extern int status;
    void asnmessage (char *);

    free(trlbuf.buffer);
    trlbuf.buffer = realloc(NULL, initLength);

    if (trlbuf.buffer == NULL) {
        asnmessage ("Out of memory: Couldn't resize buffer for trailer file.");
        status = OUT_OF_MEMORY;
    } else {
        trlbuf.buffer[0] = '\0';
    }
}
static void WriteTrlBuf (char *message)
{
    /*
        Write out message to trailer file, if any exist.
        If no file exists, then append message to buffer.
        Once message has been written out to a file,
        reset the trailer buffer to a single newline.
    */

    /* Check to see if trailer file is open for writing... */
    if (trlbuf.fp != NULL) {

        /* Trailer file is open, so write out buffer...*/
        if (trlbuf.preface[0] != '\0' && trlbuf.usepref == YES) {
            fprintf(trlbuf.fp,"%s\n",trlbuf.preface);
            fflush (trlbuf.fp);
            trlbuf.usepref = NO;
        }
        fprintf(trlbuf.fp,"%s\n",message);
        fflush (trlbuf.fp);

        /* ...then reset buffer, so we don't write old messages out
            next time.
        */
        if (trlbuf.buffer[0] != '\0')
            ResetTrlBuf ();

        /* The user must explicitly use 'ResetTrlPreface' in
            calling function in order to reset the preface.
            This allows the preface to be written out to
            multiple trailer files...
        */
    } else {
        /* If not, then we must append message to end of buffer. */
        AddTrlBuf (message);
    }
}
void CloseTrlBuf (void)
{
    /* Free memory allocated for this list, after writing out any final
        messages to the last used trailer file...
        It will then close any open trailer files as well.
    */

    extern int status;
    FILE *ofp;

    /* Do we have any messages which need to be written out? */
    if (trlbuf.buffer[0] != '\0') {
        /* We do, so open last known trailer file and
            append these messages to that file...
        */

        if ( (ofp = fopen(trlbuf.trlfile,"a+")) == NULL) {
            trlopenerr(trlbuf.trlfile);
            status = INVALID_TEMP_FILE;
            goto cleanup;
        }
        fprintf (ofp,"%s",trlbuf.buffer);

        /* Now that we have copied the information to the final
            trailer file, we can close it and the temp file...
        */
        fclose (ofp);
    }

    cleanup: ;

        free (trlbuf.buffer);
        free (trlbuf.preface);

        if (trlbuf.fp != NULL) {
            fclose (trlbuf.fp);
            trlbuf.fp = NULL;
        }

}
void trlmessage (char *message) {

    void asnmessage (char *);

    /* Send output to STDOUT and explicitly flush STDOUT, if desired */
    if (trlbuf.quiet == NO) {
        asnmessage (message);
    }

    /* Send output to (temp) trailer file */
    WriteTrlBuf (message);

}
void trlwarn (char *message) {

    char line[SZ_LINE+1];

    /* Create full warning message, like that output in ASNWARN */
    /* Use macro to add prefix to beginning of Warning message */
    sprintf(line,"%s",WARN_PREFIX);
    strcat (line,message);

    /* Send output to (temp) trailer file */
    trlmessage (line);

}
void trlerror (char *message) {

    char line[SZ_LINE+1];

    /* Create full warning message, like that output in ASNWARN */
    /* Use macro to add prefix to beginning of ERROR message */
    sprintf(line,"%s",ERR_PREFIX);
    strcat (line,message);

    /* Send output to (temp) trailer file */
    trlmessage (line);

}
void trlopenerr (char *filename) {
    sprintf (MsgText, "Can't open file %s", filename);
    trlerror (MsgText);
}
void trlreaderr (char *filename) {
    sprintf (MsgText, "Can't read file %s", filename);
    trlerror (MsgText);
}
void trlkwerr (char *keyword, char *filename) {
    sprintf (MsgText, "Keyword \"%s\" not found in %s", keyword, filename);
    trlerror (MsgText);
}
void trlfilerr (char *filename) {
    sprintf (MsgText, "while trying to read file %s", filename);
    trlerror (MsgText);
}
