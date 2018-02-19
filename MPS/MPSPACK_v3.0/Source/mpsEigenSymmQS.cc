char USAGE[]=
/* ===================================================================== */
"   Usage: [E [,I]]= mpsEigenSymmQS(H [,opts])                        \n\
                                                                      \n\
       Obtains eigenspectrum and eigenbasis in QSpace                 \n\
       represenation split into kept (Nkeep) and discarded            \n\
       state space.                                                   \n\
                                                                      \n\
   Options                                                            \n\
                                                                      \n\
      'Nkeep',..  number of states to keep (default: -1 = keep all)   \n\
      'Etrunc',.. keep states according to energy truncation criteria \n\
                  with E<=E0+Etrunc (default: <=0 = ignore; note that \n\
                  Etrunc is interpreted relative to lowest eigenvalue \n\
                  E0 = min[eig(H)]).                                  \n\
      'Rtrunc',.. keep states according to truncation criteria of     \n\
                  reduced density matrix with eig(R)>Rtrunc with R:=H \n\
                  (this is alternative to Etrunc and thus ignored by  \n\
                  default; note, furthermore, that here Rtrunc is     \n\
                  interpreted absolute, i.e. as is).                  \n\
                                                                      \n\
   AWb (C) Sep 2006 ; Mar 2013                                        \n\
";
// This is a MEX wrapper routine for MATLAB.
/* ===================================================================== */


#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsEigenSymmQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    QSpace<gTQ,gTD> A,Ek,Et,Ak,At;
    wbvector<double> Etot;
    wbMatrix<unsigned> D;
    int Nk,Nkeep=-1; char vflag=0, Rflag=0; unsigned r;
    double Etrunc=0;

    if (nargin && isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin<1) wblog(FL,"ERR Invalid number of I/O arguments");

    A.init(FL,argin[0],'r');

    if ((r=A.rank(FL))%2) wblog(FL,"ERR invalid rank-%d QSpace",r);
    A.skipZeroOffDiag(1E-14);

    if (nargin>1) {
       OPTS opts;
       opts.init(argin+1, nargin-1);

       opts.getOpt("Nkeep",Nkeep);
       if (!opts.getOpt("Etrunc",Etrunc)) {
       if ( opts.getOpt("Rtrunc",Etrunc)) Rflag=1; }
       vflag=opts.getOpt("-v");

       opts.checkAnyLeft();
    }

    A.init_gRG_zdim_bare(FL);

    Nk=Nkeep;
    if (Rflag==0)
         A.EigenSymmetric(Ak,At,Ek,Et,Etot,D,Nk,Etrunc);
    else A.EigenSymmetric(Ak,At,Ek,Et,Etot,D,Nk,Etrunc,
         NULL,wbperm(),0,0,-1, "desc"
    );

    if (Nkeep>=0 && vflag)
    wblog(FL,"<i> kept %d/%d", Nk,Etot.len);

    argout[0]=Etot.toMx();

    if (nargout>1) {
       mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

       mxAddField2Scalar(FL,S,"AK",Ak.toMx());
       mxAddField2Scalar(FL,S,"AT",At.toMx());
       mxAddField2Scalar(FL,S,"EK",Ek.toMx());
       mxAddField2Scalar(FL,S,"ET",Et.toMx());

       mxAddField2Scalar(FL,S,"DB",D.toMx());

       mxAddField2Scalar(FL,S,"NK",numtoMx(Nk));
       if (Rflag==0)
            mxAddField2Scalar(FL,S,"Etrunc",numtoMx(Etrunc));
       else mxAddField2Scalar(FL,S,"Rtrunc",numtoMx(Etrunc));
       argout[1]=S;
    }
}

