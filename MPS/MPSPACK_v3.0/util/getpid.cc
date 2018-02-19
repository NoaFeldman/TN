/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]=
"getpid - get process id of current MatLab process  \n\
usage: pid=getpid(); \n\n\
Wb,Jun14,06    \n";

#include "wblib.h"

/* ----------------------------------------------------------------------- *
 * THE GATEWAY ROUTINE
 * ----------------------------------------------------------------------- */
void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    int pid;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin) ExitMsg("No input argument expected.");
    if (nargout>1) ExitMsg("No many output arguments specified.");

    pid=getpid();

    argout[0]=numtoMx(pid); // MX_USING_ANS
};

