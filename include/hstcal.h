#ifndef HSTCAL_INCL
#define HSTCAL_INCL

#define YES         1
#define NO          0

/* Macros for dusing GetKey/PutKey functions.... */
# define USE_DEFAULT    1       /* Use default if keyword is missing */
# define NO_DEFAULT     0       /* missing keyword is fatal error */

#define SZ_FNAME          255 // Use of SZ_FNAME & SZ_LINE are interchanged throughout and should therefore be identical.
#define SZ_LINE           255 // Use of SZ_FNAME & SZ_LINE are interchanged throughout and should therefore be identical.
#define MSG_BUFF_LENGTH   SZ_LINE + 1

/*
// Flag values for calibration switches.
#define DUMMY  (-1)
#define OMIT     0
#define PERFORM  1
#define COMPLETE 2
#define SKIPPED  3
#define IGNORED  4         // e.g. SHADCORR when the exposure time is zero
*/

// Flag values for calibration switches.
enum SwitchVals_{
    BLANK=(-2),
    DUMMY,
    OMIT,
    PERFORM,
    COMPLETE,
    SKIPPED,
    IGNORED,
    SKIP,
    OMITTED,
    PERFORMED
};
typedef enum SwitchVals_ SwitchVals;

/* Codes to specify whether a reference file exists or not. */
#define EXISTS_YES       1
#define EXISTS_NO        0
#define EXISTS_UNKNOWN (-1)

/* For a reference file name, this string means that a name was
   intentionally not given.
*/
#define NOT_APPLICABLE   "n/a" /* Upper case in headers, converted to
                                    lower case by CALACS */

/* nearest integer function */
#define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* Codes for goodPedigree in RefImage and RefTab to specify whether
   pedigree is good, dummy, or not specified.
*/
#define GOOD_PEDIGREE      1
#define DUMMY_PEDIGREE     0
#define PEDIGREE_UNKNOWN (-1)

#define WARN_PREFIX    "Warning    "
#define ERR_PREFIX     "ERROR:    "

# define FITS_EXTN  ".fits"     /* default extension */

/* Standard string buffer for use in messages */
char MsgText[SZ_LINE+1]; // NOTE: This should NOT be here!!!

#endif
