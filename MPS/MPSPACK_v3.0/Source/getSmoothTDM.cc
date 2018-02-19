char USAGE[] =
/* =====================================================================
 * getSmoothTDM.cc
 *
 *  C++ program that takes spectral raw data and smoothes it using
 *  TDM::getSmoothTDM()                                                */

"   Usage: [ At [,Info]] = getSmoothTDM(tt, om, Adisc [, opts])        \n\
          [ At [,Info]] = getSmoothTDM(tt, ARAW, [, opts])             \n\
                                                                       \n\
     tt         discretized set of time dependent data                 \n\
     om         omega binning for Araw                                 \n\
     Adisc      raw TDM data with dimensions [ nom, nop [,nnrg]]       \n\
                with nom=length(om), nop=number of operators and,      \n\
                optionally, nnrg the number of NRG iterations          \n\
     ARAW       alternativly, instead of om and Adisc, the RAW data    \n\
                returned by tdmNRG() in partial mode can be specified  \n\
                as input                                               \n\
                                                                       \n\
   Options / Flags                                                     \n\
                                                                       \n\
    'disc'      discrete data (i.e. does not require d(om) weight)     \n\
    'func'      functional data (i.e. DOES require d(om) weight)       \n\
                default: disc for rank(A)==3 and func for rank(A)==2   \n\
                                                                       \n\
    'alpha',..  exponential decay of energies exp(-alpha dE) [0.]      \n\
                this applies alpha in time domain (F. Anders 2005)     \n\
    'sigma',..  applies broadening in frequency domain (log.gauss)     \n\
                NB! either alpha or sigma can be specified but not both\n\
                default: sigma=0.1 (with alpha ignored)                \n\
    'Lambda',.  Lambda from Wilson disretization (required with alpha) \n\
    'eps',...   smoothly replaces log-gauss with regular gaussian (1E-8)\n\
                                                                       \n\
    'nlog',..   number of bins per decade of omega scale               \n\
                for logarithmic discretization   (100)                 \n\
    'emin',..   minimum energy for omega binning (1E-8)                \n\
    'emax',..   maximum energy for omega binning (10.)                 \n\
                                                                       \n\
                Note that the first (last) bin contains all data below \n\
                (above) emin (emax).                                   \n\
   Output                                                              \n\
                                                                       \n\
      At        input spectral information transformed into time       \n\
                dependent data.                                        \n\
      Info      structure with additional information                  \n\
                                                                       \n\
   see also tdmNRG.                                                    \n\
   AWb © Jan 2007                                                      \n\
";

/* This is a MEX-file for MATLAB.
 * ===================================================================== */



#include "wblib.h"
#include "tdmdata.hh"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,k,K,M,N, disc, func, idisp, r=0, nlog=100; int e;
    double eps=1E-8, sigma=0.1, alpha=0., Lambda=-1, emin=1E-8, emax=10.;
    TDMData tdmData;
    char rawflag=0, istr[128]; istr[0]=0;

    wbvector<double> tt,om;
    wbvector<unsigned> D;
    wbMatrix<wbcomplex> az;
    wbMatrix<double> aa,a0;
    OPTS opts;

    wbarray<double> Om,Ar,AR;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    

    if (nargin) {
       if (!mxIsDblVector(argin[0]))
       wberror(FL,"Invalid first argument (real data required).");
       tt.init(FL,argin[0]);
       if (!tt.len) wblog(FL,"ERR No time data specified.");
    }

    if (mxIsStruct(argin[1])) {
       k=2; if (nargout>2 || nargin<(int)k)
       usage(FL,"Invalid number of I/O arguments.");

       if ((e=tdmData.init(argin[1],emin,emax,nlog)))
       wblog(FL,"ERR Invalid input RAW data (%d).",e);

       tdmData.tt=tt;
    }
    else {
       k=3; if (nargout>2 || nargin<(int)k)
       usage(FL,"Invalid number of I/O arguments.");

       if (!mxIsDblVector(0,L,argin[0]) ||
           !mxIsDblVector(0,L,argin[1]) ||
            mxMatSafeGuard(0,L,argin[2],-1) ||
           (r=mxGetNumberOfDimensions(argin[2]))>3 || r<2) 
       wberror(FL,"invalid input arguments (%s)",str);

       om.init(FL,argin[1]);
       AR.init((mxArray*)argin[2]);

       if (om.len!=AR.SIZE[0]) wblog(FL,
          "ERR omega does not match Araw (%d; %d)", om.len, AR.SIZE[0]);

       rawflag=1;
    }


    opts.init(argin+k, nargin-k);

    if (!rawflag) {
       double dbl=-1; unsigned u=0; i=0;
       i+=opts.getOpt("emin",dbl);
       i+=opts.getOpt("emax",dbl);
       i+=opts.getOpt("nlog",u  );
       if (i) wblog(FL,"WRN emin/emax/nlog ignored on input for alpha>0");
    }
    else {
       opts.getOpt("emin", emin ); eps=emin*10;
       opts.getOpt("emax", emax );
       opts.getOpt("nlog", nlog );
    }

    opts.getOpt("alpha", alpha); if (alpha<=0)
    opts.getOpt("sigma", sigma); else { sigma=0;
    opts.getOpt("Lambda",Lambda); }
    opts.getOpt("eps",   eps  );

    idisp=!opts.getOpt("nolog");

    disc=opts.getOpt("disc");
    func=opts.getOpt("func");

    if (disc && func) wblog(FL,
    "ERR Either 'disc' or 'func' can be specified.");

    opts.checkAnyLeft();


    if (rawflag) { if (AR.SIZE.len<3) AR.addSingletons(3);
       N=AR.dim(1);
       M=AR.dim(2);
       K=AR.dim(3);

       if (AR.SIZE.len==2) wblog(FL,
          "<i> getSmoothTDM: %d data set(s) with %d pts",M,N);
       else if (AR.SIZE.len==3) wblog(FL,
          "<i> getSmoothTDM: %d data set(s) with %d pts (x%d)",M,N,K);
       else wblog(FL,"ERR Araw has size [%s] ???", AR.sizeStr().data);

       replRange(FL,AR.data,AR.numel(),double(NAN),0.);


       tdmData.init(FL,"TDM",tt,K,M, emin,emax,nlog);

       D.init(2); D[0]=1; D[1]=N;
       Om.init2ref(om.data,D);

       for (k=0; k<K; k++)
       for (i=0; i<M; i++) {
          Ar.init2ref(AR.ref(i,k),D);
          tdmData.Add(Om,Ar,k,i);
       }
    }


    tdmData.getSpecData(om,a0);

    if (sigma>0) {
       sprintf(istr,"frequency broadened TDM data (sigma=%g)", sigma);
       wblog(FL,"<i> %s",istr);

       tdmData.getSmoothSpec(sigma);
       tdmData.Fourier(om,aa,tt,az);
    }
    else if (alpha>0) {
       if (Lambda<=1) wblog(FL,"ERR alpha requires Lambda (%g)",Lambda);
       sprintf(istr,"time-domain broadened TDM data (alpha=%g)",alpha);
       wblog(FL,"<i> %s",istr);
       tdmData.Fourier(om,az,alpha,Lambda,idisp);
    }
    else {
       sprintf(istr,"plain Fourier transformed data (no broadening)");
       wblog(FL,"<i> %s",istr);
       Wb::Fourier(om,a0,tt,az, func ? 'f' : 0);
    }

    argout[0] = az.toMx('r');

    if (nargout>1) {
       mxArray *S;
       S=mxCreateStructMatrix(1,1,0,NULL);

       mxAddField2Scalar(FL,S,"istr",   wbstring(istr).toMx());
       mxAddField2Scalar(FL,S,"alpha",  numtoMx(alpha));
       mxAddField2Scalar(FL,S,"sigma",  numtoMx(sigma));
       mxAddField2Scalar(FL,S,"eps",    numtoMx(eps));
       if (alpha>0)
       mxAddField2Scalar(FL,S,"Lambda", numtoMx(Lambda));

       mxAddField2Scalar(FL,S,"om",     om.toMx('t'));
       mxAddField2Scalar(FL,S,"aa",     aa.toMx('r'));
       mxAddField2Scalar(FL,S,"a0",     a0.toMx('r'));

       mxAddField2Scalar(FL,S,"emin",   numtoMx(emin));
       mxAddField2Scalar(FL,S,"emax",   numtoMx(emax));
       mxAddField2Scalar(FL,S,"nlog",   numtoMx(nlog));

       argout[1]=S;
    }

    printf("\n");
}


