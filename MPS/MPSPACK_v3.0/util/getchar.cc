char USAGE[] =
/* =============================================================
 * printfc.cc */

"   Usage: getchar()                                          \n\
                                                              \n\
      C++ program that reads single character from keyboard   \n\
                                                              \n\
   Wb,Jan09,09                                                \n";

/*    e.g. that also allows to use backspace on command line output
 *    for i=1:999, printfc([zeros(1,10)+8, sprintf('i=%03d',i)]); end
 *
 * This is a MEX-file for MATLAB.
 * ============================================================= */

#include <mex.h>

#include "wblib.h"

/* ------------------------------------------------------------- *
 * THE GATEWAY ROUTINE
 * ------------------------------------------------------------- */

void got_alarm(int sig) {
   fprintf(stderr,"Got signal %d\n",sig);
}

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin,  const mxArray *argin[]
){
    wbvector<int> cc(16);
    unsigned i=0;

#if 0
 /* does not have any effect when called from MatLab !!*/
    wbSigHandler SIG(FL);
    alarm(1); // send alarm / stop after 1 sec
    signal(SIGALRM,got_alarm);
#endif

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin || nargout>1)
    usage(FL,"Invalid number of I/O arguments");
    
 // getc returns -1 if Ctrl-C is pressed
    cc[i++]=getc(stdin);

    if (cc[0]==27) { // escape sequence

       cc[i++]=getc(stdin);
       cc[i++]=getc(stdin);

       if (cc[1]=='[') {
          while (1) {
             cc[i++]=getc(stdin);
             if (cc[i-1]<0 || cc[i-1]==126) break;
          }
       }
    }

    cc.len=i;
    argout[0]=cc.toMx(); // MX_USING_ANS

 // SIG.call99();

    return;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

