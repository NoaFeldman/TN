
char USAGE[] =
/* =============================================================
 * clebsch.cc */

"   Usage:                                                    \n\
                                                              \n\
   [cg,I]=clebsch(j1,j2);                                     \n\
                                                              \n\
      Fast evaluation of clebsch gordan coefficients          \n\
      by direct diagonalization of SSZ = S^2 + Sz             \n\
                                                              \n\
   C Wb,Sep22,09                                              \n";


#define _TQ  int
#define QTYPE _TQ

#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    int i,j, k;
    CGstore<_TQ> S;




    int s1=1, s2=2;

    if (nargin==2 && mxIsNumber(argin[0]) && mxIsNumber(argin[1])) {
        double dbl=0;
        mxGetNumber(argin[0],dbl); s1=int(2*dbl);
        mxGetNumber(argin[1],dbl); s2=int(2*dbl);
    }

    int n1=s1+1, n2=s2+1;
    int z1,z2, z1p,z2p, z,s, S1=abs(s1-s2), S2=s1+s2, D=n1*n2;
    wbarray<double> E(D,D);

    wblog(FL,"TST S1=%g, S2=%g -> S=%g..%g",0.5*s1,0.5*s2,0.5*S1,0.5*S2);
    k=0;

    for (z1 =-s1; z1 <=s1; z1 +=2)
    for (z1p=-s1; z1p<=s1; z1p+=2)
    for (z2 =-s2; z2 <=s2; z2 +=2)
    for (z2p=-s2; z2p<=s2; z2p+=2) {
       double c1, c2;
       qset<_TQ> q1 (s1,z1 ), q2 (s2,z2 ),
                 q1p(s1,z1p), q2p(s2,z2p);

       i=(z1 +s1)/2+n1*((z2 +s2)/2);
       j=(z1p+s1)/2+n1*((z2p+s2)/2);

       if (i>=D || j>=D) wblog(FL,"ERR %d,%d; %d",i+1,j+1,D);

       for (s= S1; s<=S2; s+=2)
       for (z=-s; z<=s; z+=2) { qset<_TQ> q(s,z);

          c1=CG.getCoeff(FL,"SU2",q1, q2, q);
          c2=CG.getCoeff(FL,"SU2",q1p,q2p,q);


          E(i,j)+=c1*c2;
       }
    }

    CG.printStatus(FL);

    if (E.isIdentityMatrix())
       wblog(FL," :) SU(2) coefficients complete");
    else {
       E.put("E");
       wblog(FL,"ERR SU(2) coefficients incomplete (see E)");
    }




    if (nargin>2 && mxIsQSpace(argin[0]) && mxIsQSpace(argin[1])) {
       wbvector< QSpace<_TQ,double> > Fl,Fs; 
       wbvector< blockSpaceQS<_TQ,double> > H4;
       wbvector<double> ff,Etot;
       int nk=1024;
       double E0;

       QSpace<_TQ,double> Ak,At,Hk,Ht,C1,C2,C;
       wbMatrix<_TQ> Q1,Q2;
       wbvector<_TQ> qf,qe(2);
       OpMap<_TQ> M4;
       QMap<_TQ> M;
       QVec qt;

       mxInitQSpaceVec(FL,argin[0],Fl); Fl.put("Fl");
       mxInitQSpaceVec(FL,argin[1],Fs); Fs.put("Fs");
       ff.init(FL,argin[2]);            ff.put("f_");

       if (ff.len!=Fl.len) {
          if (ff.len==1) ff.Resize(Fl.len).set(ff[0]);
          else wblog(FL,"ERR length inconsistency in coupling");
       }

       if (Fl.len==0 || Fl.len!=Fs.len || Fl[0].qtype!=Fs[0].qtype)
       wblog(FL,"ERR invalid operator set (%d/%d, %s %s)",Fl.len,Fs.len,
       Fl[0].qtype.toStr().data, Fs[0].qtype.toStr().data);

       qt=Fl[0].qtype;

       getQall(Fl,Q1);
       getQall(Fs,Q2);

       if (Q1.isEmpty() || Q2.isEmpty()) wblog(FL,
       "ERR failed to determine local state space (empty)");

       CG.getQfinal(FL,qt,Q1,Q2,M);
       M.put(FL,"M"); Q1.put("Q1");

       if (Fl[0].qop.len==0) wblog(FL,
       "ERR info.qop field required for f-operator");
       qf=Fl[0].qop;

       CG.getWeightsH4(FL,M,qe,qe,qt, M4, 'N','N'); M4.put(FL,"M1");
                                          
       CG.getWeightsH4(FL,M,qf,qf,qt, M4, 'C','N');
       CG.getWeightsH4(FL,M,qf,qf,qt, M4, 'N','C');
       CG.getWeightsH4(FL,M,qf,qf,qt, M4, 'N','N');
       CG.getWeightsH4(FL,M,qf,qf,qt, M4, 'C','C');
       CG.getWeightsH4(FL,M,qe,qf,qt, M4, 'N','C'); M4.put(FL,"M2");
       CG.getWeightsH4(FL,M,qe,qf,qt, M4, 'C','N'); M4.put(FL,"M3");

       OpProdCG(FL, Fl[0],Fl[0],C1,'C','N');
       OpProdCG(FL, Fl[0],Fl[0],C2,'N','C'); C=C1; C.Plus(C2);

       C1.put("C1"); C2.put("C2"); C.put("C");

       BuildH4(Fl,Fs,ff.data,M,H4);
       nrgTruncateH4(H4,nk,1E-12,Ak,Hk,At,Ht,Etot,&E0);

       UpdateFOps(Fs,M,Ak,Fl);

       Ak.put("AK"); Hk.put("HK");
       At.put("AT"); Ht.put("HT"); H4.put("H4");
    }

    CG.put(FL,"CG");
}


