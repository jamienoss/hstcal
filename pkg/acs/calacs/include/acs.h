#ifndef ACS_INCL
#define ACS_INCL

/* acs.h generic header for calacs */
/*      */
# include <stdio.h>             /* To insure that FILE is defined for TrlPtr */
# include "imphttab.h"
# include "hstcal.h" // NOTE: This shouldn't be here

# define ACS_CBUF           24  /* small buffer for e.g. rootname */
# define ACS_FNAME          255 // Use of ACS_FNAME & ACS_LINE are interchanged throughout and should therefore be identical.
# define ACS_LINE           255 // Use of ACS_FNAME & ACS_LINE are interchanged throughout and should therefore be identical.
# define ACS_FITS_REC       82
# define SZ_STRKWVAL        68

#define ACS_WFC_N_COLUMNS_PER_CHIP_INCL_OVERSCAN 4144
#define ACS_WFC_N_COLUMNS_PER_CHIP_EXCL_OVERSCAN 4096

typedef unsigned char Byte;

#define SIZE_BYTE   8

# define MAX_DQ     65535

# define ATOD_SATURATE 65534

/* Number of lines to extract from binned images for unbinning */
# define SECTLINES  2

/* Three extensions per SingleGroup. */
# define EXT_PER_GROUP 3

void errchk ();                 /* HSTIO error check */

/* Integer codes for string-valued header keywords. */


/* Possible values for detector. */
# define UNKNOWN_DETECTOR   (-1)
# define WFC_CCD_DETECTOR   1
# define HRC_CCD_DETECTOR   2
# define MAMA_DETECTOR      3

/* Definitions for interpreting CCDAMP, and using multiamp struct */
# define        NAMPS           4
# define        AMPSTR1         "CD"
# define        AMPSTR2         "AB"
# define        AMPSORDER       "ABCD"
/* Define array indices for each amp for clarity of code. 22Mar99, WJH */
# define        AMP_A           0
# define        AMP_B           1
# define        AMP_C           2
# define        AMP_D           3

# define       DEFAULT_OFFSET   3

# define        SM4MJD          54967

/* A reference image. */
typedef struct {
    char name[ACS_LINE];            /* name of image */
    char pedigree[ACS_FITS_REC];    /* value of pedigree keyword */
    char descrip[ACS_FITS_REC];     /* value of descrip keyword */
    int exists;                     /* does reference image exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefImage;

/* A reference table. */
typedef struct {
    char name[ACS_LINE];            /* name of table */
    char pedigree[ACS_FITS_REC];    /* value of pedigree (header or row) */
    char descrip[ACS_FITS_REC];     /* value of descrip from header */
    char descrip2[ACS_FITS_REC];    /* value of descrip from row */
    int exists;                     /* does reference table exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefTab;

/* This section is for saving the names of reference files. */
# define ACS_KEYWORD        10
typedef struct ref *RefFileInfoPtr;
typedef struct ref {
	char keyword[ACS_KEYWORD];
	char filename[ACS_FITS_REC];
	RefFileInfoPtr next;
} RefFileInfo;

typedef struct {
    float val[NAMPS];
    int colx;       /* If colx == 0, single amp readout */
    int coly;       /* if coly == 0, NOT 4-amp readout */
    int chip;       /* Chip being processed */
    int detector;   /* Which detector was used */
} multiamp;

/* The following function definitions handle the messages created
	by CALACS during operations.  These will have counterparts which
	send output both the STDOUT and a TRL file.
*/
void asnwarn (char *message);
void asnerror (char *message);
void asnmessage (char *message);

#endif
