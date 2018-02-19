char USAGE[] =
/* ================================================================== */

"   Usage: [var1, var2, ... ] = dealstruct(S,                       \n\
         'fldname1' [,idx1 ], 'fldname2' [,idx2 ], ...)             \n\
                                                                    \n\
      assigns the data of fields with names fldname1, fldname2, ... \n\
      to the variables var1, var2, ...                              \n\
                                                                    \n\
   S may also be a structure array of arbitrary dimension.          \n\
   Index options are as follows:                                    \n\
                                                                    \n\
    * if field is empty for some element, NaN is returned           \n\
    * idx=0     takes first element if array is present (default)   \n\
    * idx=i     takes element i if array is present                 \n\
    * idx='end' takes last element if array is present              \n\
                                                                    \n\
   AWb © Feb 2006                                                   \n";

/* NB! speed gain of about 40% ====================================== *
 * when compared to usual MatLab a=S.a syntax
 *
 * when including index check, e.g. also check whether empty
 * speed gain about factor 200 (!)  [cf. field2arr script]
 * e.g. 1.416259 seconds vs. 0.006320 seconds  Wb,Aug08,06
 * when run from command line and also from script(!!)
 * i.e. no significant effect by JIT acceleration
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

/* ----------------------------------------------------------------------- *
 * THE GATEWAY ROUTINE
 * ----------------------------------------------------------------------- */

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

/* single input -? or -h just displays usage */
   if (isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin<2) wberror(FL,"Invalid number of I/O arguments.");

   const mxArray *S=argin[0];
   const mwSize *ip=mxGetDimensions(S);
   mxArray *a;

   unsigned i,j,k,l,r,s, m, n=mxGetNumberOfDimensions(S);
   char fname[64];
   int fid;

   wbvector<int> idx;
   wbarray<wbcomplex> D; // allow for complex data
   wbvector<unsigned> SIZE;


   if (!mxIsStruct(S)) wberror(FL,
   "First argument must be structure.");

   SIZE.init(n);
   for (i=0; i<n; i++) SIZE[i]=(unsigned)ip[i];

// get number of fields including index specification
   idx.init(nargin-1); idx.set(-1);

   for (k=0,i=1; i<(unsigned)nargin; i++, k++) {
      if (!mxIsChar(argin[i])) wberror(FL,"ERR Invalid field set.");

   // check field index specifications (next in input line)
      if (i+1<(unsigned)nargin) {
         if (mxIsChar(argin[i+1])) {
            if (!mxGetString(argin[i+1],str))
            if (!strcmp(str,"end")) { idx[k]=INT_MAX; i++; }
         }
         else {
            if (mxGetNumber(argin[i+1],r)) wberror(FL,
            "ERR Invalid field subspecs (%s).", str);
            else {
               if (r==0) wberror(FL,"Index must be 1-based.");
               idx[k]=(r-1); i++;
            }
         }
      }
   }

   idx.Resize(k);

   if (nargout!=(int)k) wblog(FL,
   "ERR Output args do not match number of fields (%d)", k);

   D.init(SIZE); // both MatLab and wbarray are col-major
   s=SIZE.prod();

   for (k=0,l=1; l<(unsigned)nargin; k++, l++) {
      D.set(wbcomplex(NAN,0));

      if (mxGetString(argin[l], fname, 60)) wblog(FL,
      "ERR Invalid field name (%d) ???", l+1);

   // printf("... field %d: %s (%d)\n", k+1, fname, idx[k]);

      fid=mxGetFieldNumber(S,fname);
      if (fid<0) wblog(FL,"ERR Invalid field name (%s) ???",fname);

      for (j=i=0; i<s; i++) {
          a=mxGetFieldByNumber(S,i,fid); // MatLab is col-major
          if (!a) continue; // mayhap (shows up as empty in MatLab)

          m=mxGetNumberOfElements(a);
          if (m==0) continue;

          if (mxMatSafeGuard(FL,a,-1,'C','i')) wberror(FL,str);

          if (idx[k]>=0) { 
             j = (idx[k]==INT_MAX ? m-1 : (unsigned)idx[k]);
             if (j>=m) wblog(FL,"ERR Index out of bounds (%d/d)",j+1,m);
          }

          D[i].r=mxGetPr(a)[j]; if (mxIsComplex(a)) 
          D[i].i=mxGetPi(a)[j];
      }

      argout[k]=D.toMx();

      if (idx[k]>=0) l++;
   }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

