char USAGE[] = 
/* =============================================================
 * cpmat2mat.cc */

"   Usage: cpmat2mat([FLAGS,] <from.mat>, <to.mat>, <varnames ...>)\n\
                                                               \n\
      MEX file to append specified variables from one mat      \n\
      (MatLab binary) to another. If no variables are specified,\n\
      all variables in from.mat will be copied.                \n\
      Existing variables in to.mat will be overwritten.        \n\
                                                               \n\
      A variable name can be altered by using a double colon   \n\
      separator as in vold1::vnew1.                            \n\
                                                               \n\
   Flags                                                       \n\
                                                               \n\
      -v verbose                                               \n\
      -f force (overwrites target file)                        \n\
      -? displays this usage.                                  \n\
                                                               \n\
   AWb © Oct 2006                                              \n";

/* This is a MEX-file for MATLAB.
 * ============================================================= */

#define MAIN

#include "wblib.h"

int main(int nargin, const char *argin[]) {

   int e, k=1, force=0, verbose=0;
   mxArray *a, *b;

   wbstring vname;
   char *vin, *vout;

   while (k<nargin) {
      if (isHelpIndicator(argin[k])) { usage(); return 1; }
      if (!strcmp(argin[k],"-f")) { k++; force=1; continue; }
      if (!strcmp(argin[k],"-v")) { k++; verbose=1; continue; }

      break;
   }

   if (nargin<k+2) { fprintf(stderr,
    "\nERR Invalid number of I/O arguments (%d)\n"
      "Type %s -? for more information.\n\n", nargin-k, argin[0]);
      return 1;
   }

   if (!strcmp(argin[k],argin[k+1])) wblog(FL,
   "ERR Input and output file are the same!");

// alreayd checks wheather files exist
   Wb::matFile fin (FL,argin[k++],"r"); // read-only
   Wb::matFile fout(FL,argin[k++],"u"); // update

   if (k<nargin) {
      if (verbose) printf("\nCopy variables from `%s' to `%s':\n\n",
      argin[k-2], argin[k-1]);

      for (; k<nargin; k++) {
         vname=argin[k]; vin=vname.data;

         vout=strstr(vin, "::");
         if (!vout) vout=vin;
         else { vout[0]=0; vout+=2; }

         a=matGetVariable(fin.mfp,vin);
         if (!a) {
            wblog(FL,"WRN Input variable `%s' does not exist.", vin);
            continue;
         }

         b=matGetVariableInfo(fout.mfp,vout);
         if (b) { wblog(FL,
           "WRN Output variable `%s' already exists - %s.",
            vout, force ? "overwrite" : "skip");
            mxDestroyArray(b);

            if (!force) { mxDestroyArray(a); continue; }
         }

         e=matPutVariable(fout.mfp,vout,a);
         if (e) wblog(FL,
         "ERR Failed to write variable `%s' (%lX, %d)", vout, a, e);
         else if (verbose) {
            if (vin==vout) printf("  %s\n", vin);
            else printf("  %s -> %s\n", vin, vout);
         }

         mxDestroyArray(a);
      }
      if (verbose) printf("\n");
   }
   else {
      const char *v;

      if (verbose) printf("\nCopy all variables from `%s' to `%s':\n\n",
      argin[k-2], argin[k-1]);

      while ((a=matGetNextVariable(fin.mfp,&v))) {
         b=matGetVariableInfo(fout.mfp,v);
         if (b) { wblog(FL,
           "WRN Output variable `%s' already exists - %s.",
            v, force ? "overwrite" : "skip");
            mxDestroyArray(b);

            if (!force) { mxDestroyArray(a); continue; }
         }

         e=matPutVariable(fout.mfp,v,a);
         if (e) wblog(FL,
         "ERR Failed to write variable `%s' (%lX, %d)", v, a, e);
         else if (verbose) printf("  %s\n", v);

         mxDestroyArray(a);
      }
      if (verbose) printf("\n\n");
   }

   return 0;
}

