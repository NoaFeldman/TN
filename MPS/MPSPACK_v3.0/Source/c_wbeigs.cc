
#include "wbeigs.cc" // includes USAGE[] and mexFunction()

#include "wblib.h"
#include "wbeigs.hh"


int main (int argc, char **argv) {

   unsigned nargin=2, nargout=2, N=256;
   mxArray *argout[nargout];
   const mxArray *argin[nargin];

   wbMatrix<double> H;
   wbvector<double> psi;

   wblog(FL,"<i> init H and |0> to random (D=%d)%N",N);
   H.init(N,N); H.setRand().Symmetrize();
   psi.init(N); psi.setRand().normalize();

   argin[0]=H.toMx();
   argin[1]=psi.toMx();

   mexFunction(nargout, argout, nargin, argin);
}
   

