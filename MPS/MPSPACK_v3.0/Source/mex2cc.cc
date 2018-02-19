char c_USAGE[]=
/* ================================================================== */
"  Usage: david2SiteQS <args.mat>                                    \n\
                                                             \n\
     Wrapper program for the MatLab mex function david2SiteQS.mexa64 \n\
     Type david2SiteQS -H to display usage of david2SiteQS.                  \n\
                                                             \n\
  Wb,Mar01,10                                                \n\
";
/* ================================================================== */

#define MAIN

#include "david2SiteQS.cc"  // contains mexFunction()

int main (int nargin, char **argin) {

  if (nargin!=2) {
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
     int i, nargin, nargout; mxArray *a;
     OPTS opts; opts.init(argin[1]);

     opts.getOpt("nargout",a);
     if (a==NULL) wblog(FL,"ERR no 'nargout' found in %s",argin[1]); 
     opts.getOpt("nargout",nargout,'!');

     opts.getOpt("varargin",a);
     if (a==NULL) wblog(FL,"ERR no 'varargin' found in %s",argin[1]); 

     if (!mxIsCell(a)) { wblog(FL,"ERR varargin of type cell required."); }
     nargin=mxGetNumberOfElements(a);
     if (nargin<0) { wblog(FL,"ERR invalid varargin (%d)",nargin); }

     wbvector<mxArray*> A(nargin), R(nargout);
     for (i=0; i<nargin; i++) {
        A[i]=mxGetCell(a,i);
        if (A[i]==NULL)  wblog(FL,
        "ERR invalid varargin{%d/%d} (NULL)",i+1,nargin);
     }

     wblog(FL,"TST nargin=%d, nargout=%d",A.len,R.len); 
     mexFunction(R.len,R.data,A.len,(const mxArray**) A.data);

     opts.checkAnyLeft();

     a=mxCreateCellMatrix(1,R.len);
     for (unsigned i=0; i<R.len; i++) mxSetCell(a,i,R[i]);

     Wb::matFile F;
     F.Open(FL,argin[1],"u");
     if (matPutVariable(F.mfp,"varargout",a)) wblog(FL,
     "ERR failed to write varargout to file %s",argin[1]);
     F.close();

     return 0;
  }

  return 0;
}

