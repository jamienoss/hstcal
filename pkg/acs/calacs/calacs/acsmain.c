# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "acs.h"
# include "hstcalerr.h"
# include "acsversion.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

/* CALACS driver: retrieves input table name and calls pipeline
*/

int	status;		/* value of zero indicates OK */

int main(int argc, char **argv) {

	/* Local variables */
	char input[ACS_FNAME+1];
	int printtime = NO;	/* print time after each step? */
	int save_tmp = DUMMY;	/* save temporary files? */
	int verbose = NO;	/* print info during processing? */
	int debug = NO;		/* print debug statements during processing? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int onecpu = NO;		/* suppress OpenMP usage? */
    int gen1cte = NO; //Use gen1cte algorithm rather than gen2 (default)
    char pcteTabNameFromCmd[255];
    *pcteTabNameFromCmd = '\0';
	int too_many = NO;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */
    unsigned nThreads = 0;

	/* Function definitions */
	void c_irafinit (int, char **);
	int CalAcsRun (char *, int, int, int, int, const unsigned nThreads, const int gen1cte, const char * pcteTabNameFromCmd);
    void WhichError (int);
    
	/* Initialize status to OK and MsgText to null */
	status = ACS_OK;
	MsgText[0] = '\0';
	input[0] = '\0';

	/* Initialize IRAF environment */
	c_irafinit(argc, argv);

	/* Command line arguments: 
	**       0. Check for --version option
  **
  **       1. input file name
	**		   2. print time?
	**		   3. save intermediate files?
	**		   4. verbose?
	*/
	for (i = 1;  i < argc;  i++)
	{
		if (!(strcmp(argv[i],"--version")))
		{
		    printf("%s\n",ACS_CAL_VER);
		    exit(0);
		}
        if (argv[i] && strcmp(argv[i], "--gen1cte") == 0)
        {
            gen1cte = YES;
            continue;
        }
        else if (strncmp(argv[i], "--pctetab", 9) == 0)
        {
            if (i + 1 > argc - 1)
            {
                printf("ERROR - no pctetab specified\n");
                exit(1);
            }
            strcpy(pcteTabNameFromCmd, argv[i+1]);
            printf("WARNING: using pcteTab file '%s'\n", pcteTabNameFromCmd);
            ++i;
            continue;
        }
        else if (strncmp(argv[i], "--nthreads", 10) == 0)
        {
            if (i + 1 > argc - 1)
            {
                printf("ERROR - number of threads not specified\n");
                exit(1);
            }
            ++i;
            nThreads = (unsigned)atoi(argv[i]);
            if (nThreads < 1)
                nThreads = 1;
#ifndef _OPENMP
            printf("WARNING: '--nthreads <N>' used but OPENMP not found!\n");
            nThreads = 1;
#endif
            continue;
        }
        if (argv[i][0] == '-')
        {
            if (argv[i][1] == '-')
            {
                printf ("Unrecognized option %s\n", argv[i]);
                exit (ERROR_RETURN);
            }
            for (j = 1;  argv[i][j] != '\0';  j++)
            {
                if (argv[i][j] == 't') {
                    printtime = YES;
                } else if (argv[i][j] == 's') {
                    save_tmp = YES;
                } else if (argv[i][j] == 'v') {
                    verbose = YES;
                } else if (argv[i][j] == 'd') {
                    debug = YES;
                } else if (argv[i][j] == 'q') {
                    quiet = YES;
                } else if (argv[i][j] == '1') {
                    onecpu = YES;
                } else {
                    printf ("Unrecognized option %s\n", argv[i]);
                    exit (ERROR_RETURN);
                }
            }
        }
        else if (input[0] == '\0')
            strcpy (input, argv[i]);
        else
            too_many = YES;
	}
	
	if (input[0] == '\0' || too_many) {
        printf ("CALACS Version %s\n",ACS_CAL_VER_NUM);
	    printf ("syntax:  calacs.e [-t] [-s] [-v] [-q] [-1|--nthreads <N>] [--gen1cte] [--pctetab <path>] [input \n");
	    exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);
	
    if (gen1cte == YES)
    {
        sprintf(MsgText, "WARNING: using older gen1 CTE algorithm");
        trlwarn(MsgText);
    }

    if (*pcteTabNameFromCmd != '\0')
    {
        sprintf(MsgText, "WARNING: using cmd line specified PCTETAB file: '%s'", pcteTabNameFromCmd);
        trlwarn(MsgText);
    }

#ifdef _OPENMP
    unsigned ompMaxThreads = omp_get_num_procs();
#endif
    if (onecpu)
    {
        if (nThreads)
        {
            sprintf(MsgText, "WARNING: option '-1' takes precedence when used in conjunction with '--nthreads <N>'");
            trlwarn(MsgText);
        }
        nThreads = 1;
    }
    else if (!nThreads)//unset
    {
#ifdef _OPENMP
        nThreads = ompMaxThreads;
#else
        nThreads = 1;
#endif
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
    if (nThreads > ompMaxThreads)
    {
        sprintf(MsgText, "System env limiting nThreads from %d to %d", nThreads, ompMaxThreads);
        nThreads = ompMaxThreads;
    }
    else
        sprintf(MsgText,"Setting max threads to %d out of %d available", nThreads, ompMaxThreads);

    omp_set_num_threads(nThreads);
    trlmessage(MsgText);
#endif

	/* Call the CALACS main program */
	if (CalAcsRun (input, printtime, save_tmp, verbose, debug, nThreads, gen1cte, pcteTabNameFromCmd)) {

        if (status == NOTHING_TO_DO){
            /* If there is just nothing to do, 
                as for ACQ images, just quit...     WJH 27Apr1999 */
            status = 0;
	        sprintf (MsgText, "CALACS did NOT process %s", input);
	        trlmessage (MsgText);
            exit(0); 
        } else {
	        /* Error during processing */
	        sprintf (MsgText, "CALACS processing NOT completed for %s", input);
		    trlerror (MsgText);
            /* Added 19 Mar 1999 - provides interpretation of error for user */
            WhichError (status);
		    CloseTrlBuf ();
	        exit (ERROR_RETURN);
        }
	}

	/* Successful completion */
	sprintf (MsgText, "CALACS completion for %s", input);
	trlmessage (MsgText);

	CloseTrlBuf ();
	
	/* Exit the program */
	exit(0);
}
