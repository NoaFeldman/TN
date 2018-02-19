char USAGE[] = 
/* =============================================================
 * cpmat2mat.cc */

"   Usage: rmDMA(<nrg_tag> [, <varnames ...>])\n\
                                                               \n\
      MEX file to delete operator data of dmNRG from NRGWilson \n\
      data. <nrg_tag> specifies a set NRG .mat files (MatLab   \n\
      binaries).                                               \n\
                                                               \n\
      By default, all variables ending with *KK, *KT, *TK, *TT \n\
      will be removed. If only specific variables should be    \n\
      deleted they can be specified as additional arguments to \n\
      this routine. Regular expressions are accepted.          \n\
                                                               \n\
   Flags                                                       \n\
                                                               \n\
      -v verbose                                               \n\
      -? displays this usage.                                  \n\
                                                               \n\
   AWb © Nov 2006                                              \n";

/* This is a MEX-file for MATLAB.
 * ============================================================= */

#define MAIN

#include "wblib.h"
#include "nrgdata.hh"

#include <regex.h>

#define _TQ int
#define _TD double

unsigned NERR;


int matDeleteVarPattern(
   const char *F, int L,
   const Wb::matFile &f, regex_t &re
){
      
   int i,n,ndel=0,e=0;
   char **vn = matGetDir(f.mfp, &n);

   if (!vn && n<0) wblog(F,L,
   "ERR Error reading variable listing from file\n`%s' (%d)", f.fname.data);
   str[0]=0;

   for (i=0; i<n; i++) {
      if (regexec(&re, vn[i], (size_t)0, NULL, 0)==0) {
         e=matDeleteVariable(f.mfp, vn[i]);
         ndel++;
         if (e) sprintf(str+strlen(str), "%s, ", vn[i]);
      }
   }

   if (str[0]) {
      str[strlen(str)-2]=0;
      wblog(F,L,"WRN Error occured while deleting variable(s)\n%s: %s",
      f.fname.data, str); NERR++;
   }

   return ndel;
}


int main(int nargin, const char *argin[]) {

   unsigned verbose=0;
   const char *p, *f;
   int ndel=0, k=1, N=0;
   Wb::matFile F;

   const char *P=".+(KK|KT|TK|TT|_)$";

   NRGData<_TQ,_TD> NRG("A");
   regex_t re;

   NERR=0;

   while (k<nargin) {
      if (isHelpIndicator(argin[k])) { usage(); return 1; }
      if (!strcmp(argin[k],"-v")) { k++; verbose=1; continue; }

      break;
   }

   if (k>=nargin) { fprintf(stderr,
    "\nERR Invalid number of I/O arguments (%d)\n"
      "Type %s -? for more information.\n\n", k-nargin, argin[0]);
      return 1;
   }

   f=argin[k++];

   if (k<nargin) {
       for (str[0]=0; k<nargin; k++)
       sprintf(str+strlen(str), "\\<%s\\>|", argin[k]);
       strcat(str, "_$");
       p=str;
   }
   else p=P;

   if (regcomp(&re, p, REG_EXTENDED | REG_NOSUB))
   wblog(FL,"ERR Invalid variable pattern\n`%s'", p);

   if (F.open(FL,argin[1],"u")) {
      if (verbose) wblog(FL,
      "<i> Removing variables from file\n`%s'",f);

      ndel+=matDeleteVarPattern(FL,F,re);
   }
   else {
      if (verbose) wblog(FL,"<i> Removing operators from NRG data set");

      NRGData<_TQ,_TD>::SetupIO(FL,f);

      NRG.checkVec("K",3); N=NRGData<_TQ,_TD>::N;
      if (verbose) wblog(FL,"\n`%s' (%d)", f, N);

      for (k=0; k<N; k++) {
         F.Open(FL, NRG.getFileName(k),"u");
         ndel+=matDeleteVarPattern(FL,F,re);
         F.close();
      }
   }

   if (NERR)
        wblog(FL,"WRN Errors occured (%d)", NERR);
   else if (N) {
      if (ndel)
           wblog(FL,"<i> Succesfully removed %d variable(s) from %d files.", ndel, N);
      else wblog(FL,"<i> No NRG variables found / removed.");
   }
   else {
      if (ndel)
           wblog(FL,"<i> Succesfully removed %d variable(s) from file.", ndel);
      else wblog(FL,"<i> No variables found.");
   }

   return 0;
}


