char USAGE[] =
/* =================================================================
 * printfc.cc */

"   Usage: printfc(fmt)                                           \n\
                                                                  \n\
      C++ program that emulates C printf                          \n\
      ie. does not terminated automatically with a newline and    \n\
      allows for escape sequences (in contrast to fprintf).       \n\
                                                                  \n\
   NB! If more than one onput argument is specified, all          \n\
   arguments are forwarded to matlab's sprintf.                   \n\
                                                                  \n\
   AWb (C) Jan 2006                                               \n";

/*    e.g. that also allows to use backspace on command line output
 *    for i=1:999, printfc([zeros(1,10)+8, sprintf('i=%03d',i)]); end
 *
 * This is a MEX-file for MATLAB.
 * =============================================================== */

#include <mex.h>

#include "wblib.h"

/* ------------------------------------------------------------- *
 * THE GATEWAY ROUTINE
 * ------------------------------------------------------------- */

void mexFunction(
    int nargout, mxArray *argin[],
    int nargin, const mxArray *argout[]
){
    if (!nargin) return;
    if (nargout) usage(FLINE, "No output arguments available!");

 /* single input -? or -h just displays usage */
    if (nargin) if (isHelpIndicator(argout[0])) { usage(); return; }

    unsigned i,l,m,n; char *s, *s2;
    mxArray *mstr, *aa[nargin];
    
    memcpy(aa,argout,nargin*sizeof(mxArray*));

 // if (nargin>1) { usage(FLINE,
 // "\n   ERR got more than one input argument"
 // "\n   ==> use matlab sprintf or fprintf, instead."); }

    if (nargin==1)
    if (mxIsChar(argout[0])) {
    // NB! the sprintf() command in MatLab would complain about
    // \e escape sequences

       m=1+mxGetNumberOfElements(argout[0])*sizeof(mxChar);

       s = new char[m];
       mxGetString(argout[0],s,m);

    // replace escape sequences
    // NB! in the string "abc\n" the \n is already replaced by the compiler(!)
       s2 = new char[m]; n=strlen(s);
       strcpy(s2,s);

       for (l=i=0; i<n; i++) {
          if (s[i]!='\\') { s2[l++]=s[i]; continue; }
          switch (s[++i]) {
             case 'n': s2[l++]='\n'; break;
             case 't': s2[l++]='\t'; break;
             case 'r': s2[l++]='\r'; break;
             case 'b': s2[l++]='\b'; break;
             case 'e': s2[l++]='\e'; break;

             default : wblog(FL,
             "unknown escape sequence %c%c<%d>", s[i-1], s[i], s[i]);
          }
       }

       s2[l]=0;
       printf("%s",s2);

       delete [] s; delete [] s2; fflush(0);
       return;
    }

    i=mexCallMATLAB(1, &mstr, nargin, aa, "sprintf");
    if (i!=0) return;

    m=1+mxGetNumberOfElements(mstr)*sizeof(mxChar);

    s = new char[m];

    mxGetString(mstr,s,m); /* mexPrintf("%s", istr); */
    mexPrintf(s);

    mxDestroyArray(mstr);
    delete [] s; fflush(0);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

