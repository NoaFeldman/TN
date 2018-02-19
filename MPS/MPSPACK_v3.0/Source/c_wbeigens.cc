/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]=
"   Usage: $ ./c_wbeigens                                        \n\
                                                                 \n\
       rest rougine for parallelization of                       \n\
       eigenvalue decomposition based on BLAS dsyev()            \n\
       for (genearlized) symmetric rank-2 object.                \n\
                                                                 \n\
   Wb,Apr01,13                                                   \n";

#define MAIN
#include "wblib.h"


int main (int nargin, char **argin) {


    wbarray<double> H,U;
    wbvector<double> E;

    int j, D, nrepeat=16;
    double Dx=4, Dmax=2048;
    char istr[128];

    WbClock myclock("timer_wbeig");

    while (1) { Dx*=sqrt(2); D=int(Dx); if (D>Dmax) break;
       myclock.start();
       for (j=0; j<nrepeat; ++j) {
          H.init(D,D).setRand('s'); H.Symmetrize(FL); U.init(); E.init();
          myclock.resume(); wbEigenS(H,U,E);
          myclock.stop();
       }
       sprintf(istr," %4d ",D);
       myclock.info(istr);
    }
}

