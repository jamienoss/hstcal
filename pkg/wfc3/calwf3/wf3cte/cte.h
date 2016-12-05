#ifndef CTE_INCL
#define CTE_INCL

#define NUM_SCALE 4 /*number of scaling points, this is the 4 columns in the second table extension*/
#define CTEFLAG 9999999 /*flag to ignore value in array during cte calculation*/

#include "wf3.h"
#include "hstio.h"
#include "wf3info.h"
#include "../../../../ctegen2/ctegen2.h"

/* structure to hold CTE parameters from the reference files */
typedef struct {
    double scale512[RAZ_COLS]; /*scaling appropriate at row 512 */
    double scale1024[RAZ_COLS];/*scaling appropriate at row 1024 */
    double scale1536[RAZ_COLS];/*scaling appropriate at row 1536 */
    double scale2048[RAZ_COLS];/*scaling appropriate at row 2048 */
    double cte_date0; /*date of uvis install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    double scale_frac; /*scaling of cte model relative to ctedate1*/
    int noise_mit; /*read noise mitigation algorithm*/
    int wcol_data[TRAPS]; /*trap number, insync with number of traps*/
    int iz_data[RAZ_COLS]; /*column number in raz format*/
    int fix_rocr; /*make allowance for readout cosmic rays*/
    char descrip2[SZ_LINE+1]; /*descrip from table row, not read in for cte purposes*/
    char cte_name[SZ_LINE+1]; /*name of cte algorithm */
    char cte_ver[SZ_LINE+1]; /*version of algorithm */
    CTEParams baseParams;
} WF3CTEParams;


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

/*FUNCTION SIGNATURES FOR CTE SPECIFIC CODE*/

int GetCTEPars (char *, WF3CTEParams *);
void initCTEParams(WF3CTEParams *);
int doCteBias (WF3Info *, SingleGroup *);
int GetCTEFlags (WF3Info *, Hdr *);
int a2d_raz(WF3Info *);
int raw2raz(WF3Info *, SingleGroup *, SingleGroup *, SingleGroup *);
int raz2rsz(WF3Info *, SingleGroup *, SingleGroup *, double , int );
int findPostScanBias(SingleGroup *, float *, float * );
int findPreScanBias(SingleGroup *, float *, float *);
int find_dadj(int ,int , double [][RAZ_ROWS], double [][RAZ_ROWS], double , double *);
int rsz2rsc(WF3Info *, SingleGroup *, SingleGroup *, WF3CTEParams * );
int inverse_cte_blur(SingleGroup *, SingleGroup *, SingleGroup *, WF3CTEParams *, int, double);
//int sim_colreadout_l(double *, const double *, const CTEParams *, const unsigned nRows);
int CompareCTEParams(SingleGroup *, WF3CTEParams *);
int cteHistory (WF3Info *, Hdr *);
int free_array(float **ptr, int rows, int columns);
int GetCTESwitch (WF3Info *, Hdr *);
int initCTETrl (char *, char *);

int makeRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
int undoRAZ(SingleGroup *, SingleGroup *, SingleGroup *);

#endif
