char c_USAGE[]=
/* ================================================================== */

"  Usage: NRGWilsonQS <parameter file.mat>                  \n\
                                                            \n\
     This is a wrapper function to the MatLab mex routine   \n\
     NRGWilsonQS.                                           \n\
                                                            \n\
     NRGWilsonQS -H displays usage of NRGWilsonQS.          \n\
     See also NRGWilsonQS MatLab mex file.                  \n\
                                                            \n\
  (C) Wb,Aug21,06 ; Wb,Aug18,10                             \n\
";

/* This is a C++ wrapper routine for mex file NRGWilsonQS.cc
 * ================================================================== */

#define MAIN

#include <omp.h>
#include "NRGWilsonQS.cc"

int main (int nargin, char **argin) {

  if (nargin<2 || nargin>3) {
     fprintf(stderr,
      "\nERR Invalid number of I/0 arguments (%d).\n"
        "Type %s -? for more information.\n\n", nargin-1, argin[0]);
     exit(1);
  }

  if (isHelpIndicator(argin[1])) { printf("\n%s\n",c_USAGE); return 0; }

  if (!strcmp(argin[1],"-H")) {
     const mxArray *a=mxCreateString("-?");
     mexFunction(0,NULL,1,&a);
     return 0;
  }
  else {
     const mxArray *aa[2];

     printf("\n  input arguments to %s ...\n",Wb::filename(argin[0]));
     for (int i=1; i<nargin; i++) {
        aa[i-1]=mxCreateString(argin[i]);
        printf("   %2d: '%s'\n",i,argin[i]);
     }; printf("\n"); fflush(0);

     mexFunction(0,NULL,nargin-1, aa);
     return 0;
  }

  return 0;
}

