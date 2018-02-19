char USAGE[] =
/* ================================================================ */
"\n\
   Usage #1: [U,I]=getSymStates(Sp,Sz);                            \n\
                                                                   \n\
     get symmetry eigenstates for symmetry operations given by     \n\
     rising operators Sp and (diagonal) z-operators Sz.            \n\
     The order of the Sz operators determines the order of the     \n\
     z-labels to be determined.                                    \n\
                                                                   \n\
   Usage #2: [U,I]=getSymStates(SOP)                               \n\
                                                                   \n\
     Moreover, Sp and Sz can be grouped according to each symmetry \n\
     by using a single vector structure Is with the mandatory      \n\
     fields .Sp and .Sz. In this case, an empty Sp is interpreted  \n\
     as Abelian quantum number. Eventually Sp and Sz are combined  \n\
     into a single sequence as in usage #1.                        \n\
                                                                   \n\
     Further optional fields                                       \n\
       .type   symmetry type (in QSpace convention)                \n\
       .qfac   applied to all z- and q-labels (scalar or vector)   \n\
       .jmap   subsequently to qfac, applies linear transform      \n\
               of q-labels (square matrix)                         \n\
                                                                   \n\
   Output: symmetry eigenstates written as full matrix in U,       \n\
   and info structure I with the fields                            \n\
                                                                   \n\
     basically reiterating the input                               \n\
       .R0    input .Sp and .Sz data                               \n\
       .type  symmetry type                                        \n\
       .JMap  J-transformation used (identity if not specified)    \n\
       .dz    number of z-operators for each symmetry              \n\
                                                                   \n\
     actual output data                                            \n\
       .dd    dimension of irreduible multiplets generated         \n\
       .RR    irreducible representations with qfac and jmap applied\n\
       .QZ    interleaved q- and z-labels.                         \n\
       .d2    dimensions of setup in QZ (non-abelian number have   \n\
              the same number appearing twice)                     \n\
                                                                   \n\
   (C) Wb,Jul05,10                                                 \n";

// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

// see also getSymmetryStates.m (MatLab script)

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getSymStates"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"


template <class TD>
int init_ops(const mxArray *a, wbvector<wbsparray<TD> > &X);


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
   int ip=-1,iz=-1,it=-1,ix=-1, ij=-1;
   unsigned i,k,np,nz,D, nm=0, nf=0, Np=0, Nz=0;
   int e=0;

   wbvector< wbsparray<gTD> > UK,X;

   wbvector< genRG_base<gTQ,gTD> > RR, R2;
   genRG_base<gTQ,gTD> R;

   wbvector< wbarray<double> > JMap;
   wbvector< wbvector<double> > qfac;
   wbvector<unsigned> dd,dz;
   wbMatrix<unsigned> iOM;
   wbsparray<gTD> U;

   QVec qq;
   QType q;


   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }

   if (nargout>2 || nargin<1 || nargin>2){
       if (nargin || nargout) wblog(FL,"ERR invalid usage");
       else { usage(); return; }
   }

   if (nargin==2) { RR.init(1); qfac.init(RR.len); JMap.init(RR.len);
      e=init_ops(argin[0],RR[0].Sp);
      if (e) wblog(FL,"ERR invalid usage\n%s",str);

      e=init_ops(argin[1],RR[0].Sz);
      if (e) wblog(FL,"ERR invalid usage\n%s",str);
   }
   else {
      int ip=-1, iz=-1;
      if (mxIsStruct(argin[0])) {
         ip=mxGetFieldNumber(argin[0],"Sp"  );
         iz=mxGetFieldNumber(argin[0],"Sz"  );
         it=mxGetFieldNumber(argin[0],"type");
         ix=mxGetFieldNumber(argin[0],"qfac");
         ij=mxGetFieldNumber(argin[0],"jmap");
      }

      if (ip<0 || iz<0) wblog(FL,"ERR invalid usage "
         "(single structure array with field .Sp and .Sz expected");

      RR.init(mxGetNumberOfElements(argin[0]));
      qfac.init(RR.len); JMap.init(RR.len);

      for (k=0; k<RR.len; ++k) {
         e=init_ops(mxGetFieldByNumber(argin[0],k,ip), RR[k].Sp);
         if (e) wblog(FL,
            "ERR %s() got invalid/empty (%d).Sp field\n%s",FCT,k+1,str);

         e=init_ops(mxGetFieldByNumber(argin[0],k,iz), RR[k].Sz);
         if (e) wblog(FL,
            "ERR %s() got invalid/empty (%d).Sz field\n%s",FCT,k+1,str);

         if (it>=0) {
            QVec Q(FL,mxGetFieldByNumber(argin[0],k,it),1);
            RR[k].q=Q[0];
         }

         nz=RR[k].Sz.len;

         if (ix>=0) {
            wbvector<double> &x=qfac[k];
            x.init(mxGetFieldByNumber(argin[0],k,ix));
            if (x.len) {
               if (x.anyEqual(0)) wblog(FL,"ERR qfac contains 0 value");
               if (x.anyUnequal(1)) ++nf;
               if (x.len!=1 && x.len!=nz) wblog(FL,
                  "ERR qfac: length mismatch (%d/%d)",x.len,nz
               );
            }
         }

         if (ij>=0) {
            wbarray<double> &x=JMap[k];

            x.init(mxGetFieldByNumber(argin[0],k,ij));
            if (x.isScalar()) {
               if (nz>1) x.initIdentity(nz,0,x.data[0]);
            }
            else if (!x.isEmpty() && (!x.isSMatrix() || x.SIZE[0]!=nz)) {
               wblog(FL,"ERR invalid qfac (dimension mismatch: %s; %d)",
               x.sizeStr().data, nz);
            }
         }
         if (!JMap[k].isEmpty()) ++nm;
         else JMap[k].initIdentity(nz);
      }
   }

   dz.init(RR.len); qq.init(RR.len);

   for (Np=Nz=k=0; k<RR.len; ++k) {
      Np+=(np=RR[k].Sp.len);
      Nz+=(nz=RR[k].Sz.len); dz[k]=nz; qq[k]=RR[k].q;

      if (np) {
         genRG_struct<gTQ,gTD> cg; cg.initCommRel(FL,RR[k]);
         if (np!=nz) wblog(FL,
            "ERR %s() non-matching Sp/Sz operators (%d/%d)",FCT,np,nz
         );
      }
      else {
         if (!nz) wblog(FL,"ERR %s() got empty Sp/Sz set !?"); else
         if (nz>1) wblog(FL,
            "ERR %s() got %d z-ops with empty Sp (%d)",FCT,nz,k+1
         );
      }

      for (i=0; i<k; ++i) {
         RR[k].checkOrthoCommRel(FL,RR[i]);
      }
   }

   if (!Nz || Np>Nz) wblog(FL,
      "ERR invalid or emtpy number of Sp/Sz operators (%d/%d)",Np,Nz);

   MXPut IO;
   if (nargout>1) { IO
      .add(RR,"R0").addP(qq.toStr().toMx(),"qtype")
      .add(qfac,"qfac").add(JMap,"JMap");
   }


   if (nf) {
      for (k=0; k<RR.len; ++k) { if (qfac[k].len) {
         wbvector< wbsparray<gTD> > &Sz=RR[k].Sz;
         if (qfac[k].len==1) {
            double q=qfac[k].data[0];
            for (i=0; i<Sz.len; ++i) { Sz[i]*=q; }
         }
         else {
            const double *q=qfac[k].data;
            if (qfac[k].len!=Sz.len) wblog(FL,"ERR %s() "
               "got length mismatch %d/%d",FCT,qfac[k].len,Sz.len);
            for (i=0; i<Sz.len; ++i) { Sz[i]*=q[i]; }
         }
      }}
   }

   R.Sp.init(Np);
   R.Sz.init(Nz);
   for (ip=iz=k=0; k<RR.len; ++k) {
      for (Np=RR[k].Sp.len, i=0; i<Np; ++i) R.Sp[ip++]=RR[k].Sp[i];
      for (Nz=RR[k].Sz.len, i=0; i<Nz; ++i) R.Sz[iz++]=RR[k].Sz[i];
   }
   if (Np) {
      genRG_struct<gTQ,gTD> cg; cg.initCommRel(FL,R);

      if (!R.Sp[0].isSMatrix()) wblog(FL,
         "ERR invalid non-square Sp operator (Sp{1}: %s)",
          R.Sp[0].sizeStr().data);
      for (i=1; i<Np; ++i) if (R.Sp[i].SIZE!=R.Sp[0].SIZE) wblog(FL,
         "ERR severe Sp size inconsistency\nSp[1]: %s, Sp[%d]: %s",
          R.Sp[0].sizeStr().data, i+1, R.Sp[i].sizeStr().data);
      for (i=0; i<Nz; ++i) if (R.Sz[i].SIZE!=R.Sp[0].SIZE) wblog(FL,
         "ERR severe Sz size inconsistency\nSz[1]: %s, Sz[%d]: %s",
          R.Sz[0].sizeStr().data, i+1, R.Sz[i].sizeStr().data
      );
   }
   else {
      for (i=1; i<Nz; ++i) if (R.Sz[i].SIZE!=R.Sz[0].SIZE) wblog(FL,
         "ERR severe Sz size inconsistency\nSz[1]: %s, Sz[%d]: %s",
          R.Sz[0].sizeStr().data, i+1, R.Sz[i].sizeStr().data
      );
   }

   getSymmetryStates(FL,q,R.Sp,R.Sz,UK,dd,R2,&iOM,'v');

   D=R.Sz[0].SIZE[1];
   if (dd.sum()!=D) wblog(FL,
      "ERR severe size inconsistency %d/%d !?",dd.sum(),D);
 
   U.initCAT(FL,UK,2);
   argout[0]=U.toMx();

   if (nargout>1) {
      wbMatrix<double> QZ;
      unsigned l,d,m,n,r,ib,iq;
      wbvector<unsigned> d2(2*dz.len);

      IO.add(dz,"dz").add(dd,"dd").add(iOM,"iOM");

      for (l=n=i=0; i<RR.len; ++i) {
         d2[l++]=dz[i];
         if (RR[i].q.type>QTYPE_ASEP) d2[l++]=dz[i];
      }
      d2.len=l; n=d2.sum();
      QZ.init(D,n);

      if (nm) {
         wbarray<double> JM; JM.BlockDiag(JMap);
         for (i=0; i<R2.len; ++i) {
            R2[i].ApplyQFac(wbvector<double>(),&JM);
         }
      }
      IO.add(R2,"RR");

      for (l=ib=0; ib<R2.len; ++ib, l+=d) { d=dd[ib];
         const wbMatrix<double> &Z=R2[ib].Z;
         const qset<gTQ> &J=R2[ib].J;
         if (d!=Z.dim1) wblog(FL,"ERR %d/%d !??",d,Z.dim1);

         for (k=r=iq=0; iq<RR.len; ++iq, k+=m, r+=m) { m=dz[iq];
            if (RR[iq].q.type>QTYPE_ASEP) {
               for (i=0; i<d; ++i) {
                  cpyRange(J.data+r, QZ.ref(l+i,k), m, '!');
               }; k+=m;
            }

            for (i=0; i<d; ++i) {
            MEM_CPY<double>(QZ.ref(l+i,k),m,Z.ref(i,r)); }
         }
      }
      IO.add(QZ,"QZ").add(d2,"d2").save2(argout[1]);
   }

};


template <class TD>
int init_ops(const mxArray *a, wbvector<wbsparray<TD> > &X
){
   if (!a) { sprintf(str,"%s got empy array",shortFL(FL)); return -1; }

   if (mxIsDblMat(0,0,a)) { X.init(1); X[0].init(FL,a); }
   else {
      if (!mxIsCell(a)) { sprintf(str,
         "%s cell array or matrix expected",shortFL(FL)); return 1; }
      X.init(mxGetNumberOfElements(a));
      for (unsigned i=0; i<X.len; ++i) {
         const mxArray *c=mxGetCell(a,i);
         if (!mxIsDblMat(0,0,c)) { sprintf(str,
            "%s cell array of matrizes or matrix expected",shortFL(FL));
            return 2; }
         X[i].init(FL,c);
      }
   }
   return 0;
};


