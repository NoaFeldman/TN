char USAGE[] =
/* =====================================================================
 * getSmoothSpectral.cc
 *
 *  C++ program that takes spectral raw data and smoothes it using
 *  Spectral::getSmoothSpectral()                                      */

"   Usage: [omega, As [,Info]] = getSmoothSpectral(om, Araw [, opts])  \n\
                                                                       \n\
     om         omega binning for Araw                                 \n\
     Araw       raw spectral data; it may contain several set in       \n\
                separate columns, e.g. like data for different spins   \n\
                                                                       \n\
   Options / Flags                                                     \n\
                                                                       \n\
    'om',..     new discretized set of omega values                    \n\
    'sigma',..  sigma used in log-gauss smoothening (0.6)              \n\
    'eps',...   smoothly replaces log-gauss with regular gaussian (1E-6)\n\
                i.e. in absolute units since this routine does not know\n\
                anything about temperature(!). Typically choose eps=T. \n\
    'sigma2',.. Gaussian width sigma2*eps for w<<eps (0.5)             \n\
                                                                       \n\
    'func'      input data is not discrete but already function of omega\n\
    'nlog',..   number of bins per decade of omega scale               \n\
                for logarithmic discretization   (128 = nlog_fdmNRG/2!)\n\
    'emin',..   minimum energy for omega binning (1E-8)                \n\
    'emax',..   maximum energy for omega binning (10.)                 \n\
                                                                       \n\
    'disp',..   verbose for disp>0 (1)                                 \n\
                                                                       \n\
                Note that the first (last) bin contains all data below \n\
                (above) emin (emax).                                   \n\
   Output                                                              \n\
                                                                       \n\
      om        omega corresponding to the energy binning [-emax,+emax]\n\
      As        smoothened spectral function                           \n\
                                                                       \n\
   see also fdmNRG.                                                    \n\
   AWb (C) May 2006                                                    \n\
";

/* This is a MEX-file for MATLAB.
 * ===================================================================== */


#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getSmoothSpectral"
#endif

#include "wblib.h"
#include "spectral.hh"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,j,n, N, disp=1, nlog=128, isfunc;
    double eps=1E-6, sigma=0.6, sigma2=0.5, emin=1E-8, emax=10.;
    Spectral ASpec;

    wbvector<double> om;
    wbMatrix<double> aa;
    wbarray<double> Om, Ar;
    wbvector<unsigned> D;
    wbMatrix<double> A0, *a0=&A0;
    OPTS opts;

    mxArray *a;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    

    if (nargout<2 || nargin<2)
    usage(FL,"Invalid number of I/O arguments.");

    if (!mxIsDblVector(argin[0]) || !mxIsDblMat(argin[1]))
    wberror(FL,"Invalid input arguments (real data required).");

    om.init (FL,argin[0]); emin=om.aMin();
    aa.init0(argin[1]);


    opts.init(argin+2, nargin-2);

    isfunc=opts.getOpt("func");

    opts.getOpt("eps",   eps  );
    opts.getOpt("sigma", sigma);
    opts.getOpt("sigma2",sigma2);

    opts.getOpt("nlog",  nlog );
    opts.getOpt("emin",  emin );
    opts.getOpt("emax",  emax );

    if (eps<emin) wblog(FL,"WRN eps=%g < min(|om|)=%g", eps, emin);

    opts.getOpt("om", a);
    opts.getOpt("disp",disp);
    if (opts.getOpt("-a0")) a0=NULL;

    opts.checkAnyLeft(FL);


    N=aa.dim2;

    if (om.len!=N) wblog(FL,
       "ERR %s() number of cols of aa must match\n"
       "length of om (%d/%d)",PROG,om.len,N);
    if (disp) wblog(FL,
       "<i> %s() %d data set(s) with %d pts",PROG,aa.dim1,aa.dim2);

    for (n=i=0; i<aa.dim1; i++)
    for (  j=0; j<aa.dim2; j++)
    if (isnan(aa(i,j))) { aa(i,j)=0; n++; }

    if (n) wblog(FL, "    Skipping %d NaN's", n);

    if (isfunc) {
       wbvector<double> dom(N); unsigned m;
       wblog(FL,"<i> Smooth -> discrete (raw) data");

       for (i=1; i<N; i++) if (om[i]<=om[i-1]) wblog(FL,
       "ERR option `func' requires increasing set of omega's");

       dom[0]=om[1]-om[0];
       for (m=N-1, i=1; i<m; i++) dom[i]=0.5*(om[i+1]-om[i-1]);
       dom[i]=om[i]-om[i-1];

       for (i=0; i<aa.dim1; i++)
       for (j=0; j<aa.dim2; j++) aa(i,j)*=dom[j];
    }


    ASpec.init("A", aa.dim1, emin, emax, nlog);

    if (a0) {
       wbindex I; wbMatrix<double> ax;
       wbvector<double> a2; A0.init(aa.dim1,2);
       unsigned gotdelta=0;

       om.maxneg(&i);
       if (int(i)>=0) {
          om.findRange(2*om[i],om[i],I,"[[");
          if (!I.isEmpty()) {
             aa.getCols(I,ax).recSumA(a2); a2.setMin(1E-14);
             for (j=0; j<a2.len; j++) {
                if (fabs(aa(j,i))>a2[j]) {
                    A0(j,0)+=aa(j,i); aa(j,i)=0; ++gotdelta;
                }
             }
          }
       }

       om.minpos(&i);
       if (int(i)>=0) {
          om.findRange(om[i],2*om[i],I,"]]");
          if (!I.isEmpty()) {
             aa.getCols(I,ax).recSumA(a2); a2.setMin(1E-14);
             for (j=0; j<a2.len; j++) {
                if (fabs(aa(j,i))>a2[j]) {
                    A0(j,1)+=aa(j,i); aa(j,i)=0; ++gotdelta;
                }
             }
          }
       }

       if (gotdelta) { wblog(FL,
       " *  got delta(0) [%s]",a0->toStr("%.3g").data);
            A0.save2(ASpec.A0);
       }
       else { A0.init(); a0=NULL; }
    }

    D.init(2); D[0]=1; D[1]=N;
    Om.init2ref(om.data,D);

    for (i=0; i<aa.dim1; i++) {
       Ar.init2ref(aa.rec(i),D);
       ASpec.Add(Om,Ar,i);
    }

    if (a) {
         om.init(FL,a);
         ASpec.getSmoothSpectral(aa,om,eps,sigma,sigma2,disp); }
    else ASpec.getSmoothSpectral(om,aa,eps,sigma,sigma2,disp);

    argout[0] = om.toMx('t');
    argout[1] = aa.toMx('r');

    if (nargout>2) {
       MXPut S(FL); S
        .add(emin,"emin").add(emax,"emax").add(nlog,"nlog")
        .add(sigma,"sigma").add(eps,"eps").add(sigma2,"sigma2").add(A0,"a0");
       S.save2(argout[2]);
    }

    if (disp) printf("\n");
};


