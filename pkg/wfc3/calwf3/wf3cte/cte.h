#ifndef CTE_INCL
#define CTE_INCL

#define NUM_SCALE 4 /*number of scaling points, this is the 4 columns in the second table extension*/
#define CTEFLAG 9999999 /*flag to ignore value in array during cte calculation*/
#define TRAPS 999 /*max number of traps per column = rows in pctetab[1], valid traps are < 999999 in qlev*/

#include "wf3.h"
#include "hstio.h"
#include "wf3info.h"
#include "../../../../ctegen2/ctegen2.h"

/*USEFUL LIB FUNCTIONS*/
int GetGrp (WF3Info *, Hdr *);
void PrBegin (char *);
void PrEnd (char *);
void PrFileName (char *, char *);
void PrHdrInfo (char *, char *, char *);
void PrGrpBegin (char *, int);
void PrGrpEnd (char *, int);
int LoadHdr (char *, Hdr *);
void WF3Init (WF3Info *);
int MkName (char *, char *, char *, char *, char *, int);
int GetKeys (WF3Info *, Hdr *);
int GetSwitch (Hdr *, char *, int *);
int resistmean(float *, int, float, float *,float *,float *,float *);
void TimeStamp (char *, char *);
int  FileExists (char *);
void PrRefInfo (char *, char *,char *, char *, char *);
void PrSwitch (char *, int );
void WhichError (int);
int sub1d (SingleGroup *, int, SingleGroupLine *);
int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
int GetKeyInt (Hdr *, char *, int , int , int *);
int GetKeyDbl (Hdr *, char *, int , double , double *);
int GetKeyFlt (Hdr *, char *, int , float , float *);
int PutKeyInt (Hdr *, char *, int , char *);
int PutKeyFlt (Hdr *, char *, float , char *);
int PutKeyDbl (Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);
int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
int GetKeyStr (Hdr *, char *, int, char *, char *, int);
int streq_ic (char *, char *); /* case insensitive string equal */
int  MkOutName (char *, char **, char **, int, char *, int);
int MkNewExtn (char *, char *);
int Sub2Full(WF3Info *, SingleGroup *, SingleGroup *, int, int, int );
int Full2Sub(WF3Info *, SingleGroup *, SingleGroup *, int, int, int );
int CreateEmptyChip(WF3Info *, SingleGroup *);
int GetCorner (Hdr *, int, int *, int *);
int initChipMetaData(WF3Info *wf3, Hdr * hdr, int groupNumber);


/*FUNCTION SIGNATURES FOR CTE SPECIFIC CODE*/

enum OverScanType {
    PRESCAN,
    POSTSCAN
};

int alignAmps(SingleGroup * image, CTEParams * ctePars);
int unalignAmps(SingleGroup * image, CTEParams * ctePars);
int getSubarray(SingleGroup * image, CTEParams * ctePars, WF3Info * wf3);
int getCCDChip(int * value, char * fileName, char * ename, int ever);
int putChip(char * fileName, SingleGroup * image, WF3Info * wf3, double const scaleFraction);
int doCteBias (WF3Info *, SingleGroup *);
int GetCTEFlags (WF3Info *, Hdr *);
int a2d_raz(WF3Info *);
int biasAndGainCorrect(SingleGroup * raz, const float ccdGain, const Bool isSubarray);
int findOverScanBias(SingleGroup *raz, float *mean, float *sigma, enum OverScanType overScanType);
int cteHistory (WF3Info *, Hdr *);
int free_array(float **ptr, int rows, int columns);
int GetCTESwitch (WF3Info *, Hdr *);
int initCTETrl (char *, char *);
int outputImage(char * fileName, SingleGroup * image, CTEParams * ctePars);

int makeRAZ(const SingleGroup *, const SingleGroup *, SingleGroup *);
int undoRAZ(SingleGroup *, SingleGroup *, const SingleGroup *);

#endif
