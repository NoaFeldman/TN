char c_USAGE[]=
/* ================================================================== */

"  Usage: fdmNRG_QS <NRGdata> <setup.mat> ['C1'] 'C2'       \n\
                                                            \n\
     This is a wrapper function to the MatLab mex routine   \n\
     fdmNRG_QS(). Output will be appended to setup.mat.     \n\
                                                            \n\
     fdmNRG_QS -H displays usage of fdmNRG_QS.              \n\
     See also fdmNRG_QS.cc MatLab mex file.                 \n\
                                                            \n\
  (C) Wb,Oct16,06 ; Jul28,10                                \n\
";

/* This is a C++ wrapper routine for mex file fdmNRG_QS.cc
 * ================================================================== */

#define MAIN

#include <omp.h>
#include "fdmNRG_QS.cc"

int main (int nargin, char **argin) {

   unsigned i,n=nargin-1;
   const mxArray *aa[n+1];

   if (nargin==2) { 
      if (isHelpIndicator(argin[1])) { printf("\n%s\n",c_USAGE); return 0; }
      if (!strcmp(argin[1],"-H")) {
         const mxArray *a=mxCreateString("-?");
         mexFunction(0,NULL,1,&a);
         return 0;
      }
   }

   if (nargin<4 || nargin>5) {
      fprintf(stderr,
       "\nERR invalid number of I/0 arguments (%d).\n"
         "Type %s -? for more information.\n\n", nargin-1, argin[0]);
      exit(1);
   }

   printf("\n  input arguments to %s ...\n",Wb::filename(argin[0]));
   aa[0]=mxCreateString("--mat");
   for (i=0; i<n; i++) {
      aa[i+1]=mxCreateString(argin[i+1]);
      printf("   %2d: '%s'\n",i+1,argin[i+1]);
   }; printf("\n"); fflush(0);

   mexFunction(0,NULL,n,aa);
   return 0;
}

