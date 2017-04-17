# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "wf3.h"
# include "err.h"
# include "wf3version.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

/* H. Bushouse	07-Sep-2011	Implemented new "--version" command line argument. */
/* M. Sosey     added a -r to also print the version (-v is used, so warren chose r for revision */

/* CALWF3 driver: retrieves input table name and calls pipeline
 */

int	status;		/* value of zero indicates OK */

int main (int argc, char **argv) {

	/* Local variables */
	char input[SZ_FNAME+1];
	int printtime = NO;	/* print time after each step? */
	int save_tmp = DUMMY;	/* save temporary files? */
	int verbose = NO;	/* print info during processing? */
	int debug = NO;		/* print debug statements during processing? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int too_many = NO;	/* too many command-line arguments? */
	int onecpu = NO;  /* suppress openmp usage by using only 1 thread?*/
    Bool fastCTE = False; // Use high performance CTE implementation
    unsigned nThreads = 1;
    int i, j;		/* loop indexes */

	/* Function definitions */
	void c_irafinit (int, char **);
	int  CalWf3Run  (char *, int, int, int, int, unsigned, Bool);
	void WhichError (int);

	/* Initialize status to OK and MsgText to null */
	status     = WF3_OK;
	MsgText[0] = '\0';
	input[0]   = '\0';

	/* Initialize IRAF environment */
	c_irafinit(argc, argv);

	/* Command line arguments: 
	 **	0. Check for --version option
	 **	1. input file name
	 **	2. print time?
	 **	3. save intermediate files?
	 **	4. verbose?
	 */
	for (i = 1;  i < argc;  i++) {
		if (!(strcmp(argv[i],"--version"))) {
			printf("%s\n",WF3_CAL_VER_NUM);
			exit(0);
		}
		if (argv[i][0] == '-') {
		    if (strncmp(argv[i], "--fast-cte", 10) == 0)
		    {
		        fastCTE = YES;
		        continue;
		    }
		    else if (strncmp(argv[i], "--nThreads", 10) == 0)
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
		    else
		    {
                for (j = 1;  argv[i][j] != '\0';  j++) {
                    if (argv[i][j] == 't') {
                        printtime = YES;
                    } else if (argv[i][j] == 's') {
                        save_tmp = YES;
                    } else if (argv[i][j] == 'r'){
                        printf ("Current version: %s\n", WF3_CAL_VER);
                        exit(0);
                    } else if (argv[i][j] == 'v') {
                        verbose = YES;
                    } else if (argv[i][j] == 'd') {
                        debug = NO;
                    } else if (argv[i][j] == 'q') {
                        quiet = YES;
                    } else if (argv[i][j] == '1'){
                        onecpu = YES;
                    } else {
                        printf ("Unrecognized option %s\n", argv[i]);
                        exit (ERROR_RETURN);
                    }
                }
		    }
		} else if (input[0] == '\0') {
			strcpy (input, argv[i]);
		} else {
			too_many = YES;
		}
	}

	if (input[0] == '\0' || too_many) {
		printf ("syntax:  calwf3.e [-t] [-s] [-v] [-q] [-r] [-1] [--fast-cte [--nthreads <N>]] input \n");
		exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode (quiet);



    if (fastCTE)
    {
        sprintf (MsgText, "WARNING: using high performance CTE implementation \n"
                "Best results obtained when built with --O3 configure option");
        trlmessage (MsgText);
    }
    else if (nThreads)
    {
        sprintf(MsgText, "cmd options -nThreads requires --fast-cte");
        trlerror(MsgText);
        exit(1);
    }

    if (onecpu)
    {
        if (nThreads != 1)
        {
            sprintf(MsgText, "WARNING: option '-1' takes precedence when used in conjunction with '--nthreads <N>'");
            trlwarn(MsgText);
        }
        nThreads = 1;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
    unsigned ompMaxThreads = omp_get_num_procs();
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


	/* Call the CALWF3 main program */
	if (CalWf3Run (input, printtime, save_tmp, verbose, debug, nThreads, fastCTE)) {

		if (status == NOTHING_TO_DO) {
			/* If there is just nothing to do, 
			 ** as for ACQ images, just quit. */
			status = 0;
			sprintf (MsgText, "CALWF3 did NOT process %s", input);
			trlmessage (MsgText);
			exit(0); 
		} else {
			/* Error during processing */
			sprintf (MsgText, "CALWF3 processing NOT completed for %s",
					input);
			trlerror (MsgText);
			/* Provide interpretation of error for user */
			WhichError (status);
			CloseTrlBuf ();
			exit (ERROR_RETURN);
		}
	}

	/* Successful completion */
	sprintf (MsgText, "CALWF3 completion for %s", input);
	trlmessage (MsgText);

	CloseTrlBuf ();

	/* Exit the program */
	exit (0);
}

