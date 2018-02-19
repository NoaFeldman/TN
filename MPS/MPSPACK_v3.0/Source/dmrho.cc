#ifndef __WB_DMRHO_CC__
#define __WB_DMRHO_CC__

// NB! make sure SkipZeroData() skips only numerical noise! using plain
// DEPS=DBL_EPSILON can lead to trace(eig(rho)) inconsistency, otherwise.
// with the main difference coming from considering EIG(rho) instead of rho!
// => keep 1E-6 buffer w.r.t. numerical noise for large data blocks.
// Wb,Mar31,11
   const double DEPS2=1E-6*DBL_EPSILON;

   wbvector<double> gES;

   inline double nrgScale(int iter) {
      if (iter<0 || iter>=(int)gES.len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,iter,gES.len); 
      return gES[iter];
   }


template <class TQ, class TD>
double initRHO(
   wbMatrix<double> &rhoNorm,
   const double T,
   NRGData<TQ,TD> &H,
   wbvector<double> E0,
   const unsigned dloc,
   const unsigned NRho=0,
   char vflag=1
);

template <class TQ, class TD>
void getEData(
   const QSpace<TQ,TD> &HK,
   wbvector<TD> &E, wbvector<unsigned> &D,
   wbvector<unsigned> *dz=NULL
);

template <class TQ, class TD>
TD getErange(const QSpace<TQ,TD> &HK, TD* emin=NULL, TD* emax=NULL);

template <class TQ, class TD>
unsigned getGSDegeneracy(const QSpace<TQ,TD> &H);

template <class TQ, class TD>
bool getBoltzman(
   const QSpace<TQ,TD> &H,
   wbvector<TD> &R, wbvector<unsigned> &D, unsigned iter,
   wbvector<unsigned> *DZ=NULL
);

template <class TD>
bool getBoltzman(
   const wbvector<TD>& E, const wbvector<unsigned> &dz,
   wbvector<TD>& R, unsigned iter, double T=0.
);

template <class TQ, class TD>
double initRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &HK,
   const int iter,
   double efac=1.
);
    
template <class TQ, class TD>
void nrgUpdateRHO(
   const wbMatrix<double> &R,
   NRGData<TQ,TD> &A,
   NRGData<TQ,TD> &H,
   NRGData<TQ,TD> &RHO,
   unsigned NRho,
   unsigned dloc,
   wbSigHandler *sig=NULL
);

template <class TQ, class TD>
void updateRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &AK
);


template <class TQ, class TD>
void checkRho (const char *F, int L,
    const QSpace<TQ,TD> &Rho,
    const double traceRho=1.
){
    unsigned i,n=Rho.QIDX.dim1;
    double rsum=0, rmin=0;

    wbvector< wbvector<gTD> > E(n);
    wbvector<gTD> rr;
    wbarray<gTD> U;

    char cgflag=Rho.gotCGS(FL);

    if (cgflag>0 && Rho.qtype.hasMP()>1) {
       cgflag=2;
    }

    if (!Rho.isConsistentR(2) || !Rho.isBlockDiagMatrix() || !Rho.isHConj())
    wblog(F,L,"ERR %s() inconsistency in Rho QSpace (%d;%d;%d)",
    FCT, Rho.isBlockDiagMatrix(), Rho.isHConj());

    for (i=0; i<n; i++) {
       wbEigenS(*(Rho.DATA[i]),U,E[i]);
       if (cgflag>0) {
          unsigned d=Rho.cgsDimScalar(i);
          if (cgflag<=1) {
             unsigned d0=Rho.qtype.QDim(Rho.QIDX.rec(i)); if (d0!=d)
             wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,d,d0); 
          }
          E[i]*=d;
       }
    }
    rr.Cat(E);

    rmin=rr.min();
    rsum=rr.sum();

    if (rmin<-EPS) { wblog(F,L,
       "ERR eig(rho) is in (%.4g .. %.4g; %.4g)",rmin,rr.max(),rsum);
    }
    else if (fabs(rsum-traceRho)>EPS) {
       char e=(fabs(rsum-traceRho)>1E-8 && rsum>1E-8 && traceRho>1E-8);
       if (e) {
          printf("\n\n");
          MXPut(FL).add(Rho,"Rho").add(rr,"rr").add(traceRho,"t");
       }
       sprintf(str,"%s trace(eig(rho)) inconsistency\n"
         "%.4g / %.4g @ %.3g (%.2g; %.2g; %.2g)", e ? "ERR":"WRN",
          rsum,traceRho,rsum-traceRho,DEPS,DEPS2,EPS); printf("\n");
       wblog(F,L,str);
    }
};



template <class TQ, class TD>
double initRHO(
   wbMatrix<double> &rhoNorm,
   const double T,
   NRGData<TQ,TD> &H,
   wbvector<double> E0,
   const unsigned dloc,
   const unsigned NRho,
   char vflag
){
   unsigned i=0,i1,i2, iter, isneg=0, N=NRG_N; char i0flag=0;
   wbvector<double> R;
   wbvector<unsigned> D;
   double dbl, dfac=1;
   double DFAC=DBL_MAX*1E-2;

   if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)");

   i1 = (dloc && NRho ? NRho-1 : 0);
   i2 = (NRho ? NRho-1 : N-1);

   if (int(i1)>int(i2) || int(i1)<0) wblog(FL,
   "ERR rho_i = [%d..%d] !??",i1,i2);

   rhoNorm.init(N,2);


   dbl=E0.min(i); 
   if (i==0) { i0flag=1;
      E0[0]+=E0[1]; dbl=E0.min(i);
      E0[0]-=E0[1];
   }

   if (i<i2) {
      wblog(FL,"WRN Emin set by site %d/%d !??",i+1,N);
      MXPut(FL).add(E0,"E0");
   }

   getBoltzman(E0,D,E0,-1,T);

   for (iter=i1; iter<=i2; iter++) {
      H.init(FL,"T",iter);
      if (getBoltzman(H.T,R,D,iter)) isneg++;
      rhoNorm(iter,1)=R.sum();

      if (dloc && iter<i2) continue;

      H.init(FL,"K",iter);
      if (getBoltzman(H.K,R,D,iter)) isneg++;
      rhoNorm(iter,0)=R.sum();
   }

   if (isneg) wblog(FL,
      "WRN %s() negative (eigen)energies (%d) !??",FCT,isneg);
   if (i0flag && (rhoNorm(0,0)>1 || rhoNorm(0,1)>1)) wblog(FL,
      "WRN %s() density matrix from first site may be unstable (%.3g,%.3g)",
       FCT, rhoNorm(0,0), rhoNorm(0,1));

   if (!H.K.isEmpty())
        i=getGSDegeneracy(H.K);
   else i=getGSDegeneracy(H.T);

   if (vflag) {
      if (i==1)
           wblog(FL," *  non-degenerate ground state");
      else wblog(FL,"<i> ground state degeneracy k=%d (%d)",iter-1,i);
   }



   if (!dloc) {
      if (!NRho) wblog(FL,"ERR NRho=%d !??", NRho);
      wblog(FL,"<i> generate weights for plain NRG (T=%.3g)", T);

      for (iter=0; iter<NRho; iter++) {
         dbl=rhoNorm.recSum(iter);
         if (dbl<1) wblog(FL,"WRN sum_rhoNorm(%d,:)=%.3g",iter+1,dbl);
         dbl=1/dbl;
         rhoNorm(iter,0)*=dbl;
         rhoNorm(iter,1)*=dbl;
      }
      return 0;
   }


   if (!NRho)
   for (dfac=dloc, iter=i2-1; (int)iter>=0; iter--, dfac*=dloc) {
      rhoNorm(iter,1)*=dfac;

      if (dfac>DFAC) wblog(FL,"WRN dfac reaches DBL_MAX (%.3g)", dfac);
   }

   dbl=rhoNorm.sum();
   if (dbl<=0.) wblog(FL,"ERR %s() normRH0=%g ??? DIV/0",FCT,dbl);
   rhoNorm *= (1/dbl);

   if (vflag) {
      double rmax;
      rhoNorm.recSum(R); rmax=R.max(i);
      wblog(FL," *  maxRho=%.4g at iter=%d/%d (T=%.4g)",rmax,i+1,N,T);
   }

   return (-T*log(dbl));
};


template <class TQ, class TD>
void getEData(
   const QSpace<TQ,TD> &HK,
   wbvector<TD> &E,
   wbvector<unsigned> &D,
   wbvector<unsigned> *DZ=NULL
){
   unsigned i,l,m=0,d, e=0, k,K=HK.QIDX.dim1, *dz=NULL;
   char cgflag=HK.gotCGS(FL);

   if (HK.DATA.len!=K) wblog(FL,
      "ERR %s() severe QSpace inconsistency (%d,%d)",
       FCT,HK.DATA.len,HK.QIDX.dim1);
   if (!HK.isBlockDiagMatrix('d')) wblog(FL,
      "ERR %s() Hamiltonian must be (block)diagonal!",FCT);

   D.init(K);
   for (l=k=0; k<K; k++) D[k]=HK.DATA[k]->SIZE.max();

   E.init(D.sum());
   if (DZ && cgflag>0) { DZ->init(E.len); dz=DZ->data; }

   if (cgflag>0 && HK.qtype.hasMP()>1) {
      cgflag=2;
   }

   for (l=k=0; k<K; k++) {
       const wbarray<TD> &Hk = *HK.DATA[k];
       const wbvector<unsigned> &S = Hk.SIZE;
       const TD* const E0=Hk.data;

       d=D[k]; if (dz) {
          m=HK.cgsDimScalar(k);
          if (cgflag<=1) {
             unsigned m0=HK.qtype.QDim(HK.QIDX.rec(k)); if (m0!=m)
             wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,m,m0); 
          }
       }

       if (S[0]!=S[1]) { if (S[0]!=1 && S[1]!=1) e++; else {
          memcpy(E.data+l, E0, d*sizeof(TD));
          if (dz)
             { for (i=0; i<d; ++i,++l) dz[l]=m; }
          else l+=d;
       }}
       else { if (S.len!=2) e++; else {
          if (d>1 && (fabs(E0[1])>1E-12 || fabs(E0[d])>1E-12)) wblog(FL,
          "WRN HK[%d] not E diagonal (%g,%g)", k+1,E0[1], E0[d]);

          for (i=0; i<d; ++i,++l) {
             E[l]=E0[i*d+i]; if (dz) dz[l]=m;
          }
       }}
       if (e) wblog(FL,
          "ERR HK[%d] must be rank-2 object! (%s)",
           k+1, HK.DATA[k]->sizeStr().data
       );
   }
};


template <class TQ, class TD>
unsigned getGSDegeneracy(const QSpace<TQ,TD> &H) {

   wbvector<TD> E; wbvector<unsigned> D,dz;
   unsigned i=0, n=0;
   TD e0, *ee;

   getEData(H,E,D,&dz); e0=E.min(); ee=E.data;
   for (; i<E.len; i++) if (fabs(ee[i]-e0)<EPS) n++;

   return n;
};

template <class TQ, class TD>
TD getErange(const QSpace<TQ,TD> &H, TD* emin, TD* emax){

   if (!H.isBlockDiagMatrix('d')) {
      QSpace<TQ,TD> AK,HK; wbvector<TD> EE;
      H.EigenSymmetric(AK,HK,EE);
      return getErange(HK,emin,emax);
   }

   wbvector<TD> E; wbvector<unsigned> D;
   TD e1,e2;

   getEData(H,E,D);
      e1=E.min(); if (emin) (*emin)=e1;
      e2=E.max(); if (emax) (*emax)=e2;
   return (e2-e1);
};


template <class TQ, class TD>
bool getBoltzman(
   const QSpace<TQ,TD> &H,
   wbvector<TD> &R, wbvector<unsigned> &D, unsigned iter,
   wbvector<unsigned> *DZ
){
   wbvector<unsigned> dz;
   wbvector<TD> E;
   bool e=0;

   getEData(H,E,D,&dz); e=(E<0);

   if (getBoltzman(E,dz,R,iter)) e=1;
   if (DZ) dz.save2(*DZ);

   return e;
};


template <class TD> 
bool getBoltzman(
   const wbvector<TD>& E,
   const wbvector<unsigned>& dz,
   wbvector<TD>& R,
   unsigned iter,
   double T0
){
   static double T=0.;
   static wbvector<double> E0;

   unsigned i,e=0;

   if (int(iter)<0) { T=T0; E0=R; return 0; }

   if (iter>=E0.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,iter+1,E0.len);
   if (dz.len && dz.len!=E.len) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,dz.len,E.len);

   if (E.isEmpty()) { R.init(); return 0; }


   double beta,r, Escale=nrgScale(iter), Eref=E0[iter]/Escale;


   if (E<-Eref) { e++; wblog(FL, 
      "WRN Eref should be lower bound (%g,%g) !??",
       E.min()*Escale, Eref*Escale);
   }

   if (&R!=&E) R.init(E.len);

   if (T>0.) {
      beta=Escale/T;
      for (i=0; i<E.len; i++) {
         r=exp(-beta*(E[i]+Eref)); if (dz.len) r*=dz[i];
         R[i]=r;
      }
   }
   else {
      for (i=0; i<E.len; i++)
      R[i]=( ABS(E[i]+Eref)<EPS ? (dz.len ? dz[i] : 1) : 0. );

      i=unsigned(R.sum());
      if (i) wblog(FL," *  T=0 ground state degeneracy %d", i);
      else   wblog(FL," *  T=0 got unique ground state");
   }

   return e;
};


template <class TQ, class TD>
double initRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &H,
   const int iter,
   double rhofac
){
   unsigned i,k,l,d,n, K=H.QIDX.dim1, has0=0;
   wbvector<unsigned> D,dz;
   wbvector<TD> R;
   TD Z, *r;

   getBoltzman(H,R,D,iter,&dz);

   Z=R.sum();

   if (Z==0.) wblog(FL,"ERR rhoNorm[%d]=%g ??? DIV/0!", iter+1,Z);
   if (dz.len && R.len!=dz.len) wblog(FL,
   "ERR %s() length mismatch %d/%d",FCT,R.len,dz.len);

   R *= (rhofac/Z); r=R.data;


   if (dz.len)
   for (i=0; i<R.len; i++) {
      if (dz.data[i]>1) r[i]/=dz.data[i];
   }

   Rho=H;

   for (l=k=0; k<K; k++, l+=d) {
       wbarray<TD> &Rk = (*Rho.DATA[k]);

       d=D[k]; Rk.init(d,d);

       for (n=i=0; i<d; i++) {
          Rk.data[i*d+i] = r[l+i];
          if (!n && r[l+i]>DEPS) n++;
       }
       if (!n) has0++;
   }

   if (has0) Rho.SkipZeroData(DEPS2,'b');
   return Z;
};


template <class TQ, class TD>
void updateRho(QSpace<TQ,TD> &Rho, const QSpace<TQ,TD> &AK) {

   QSpace<TQ,TD> Xk;

   if (typeid(TD)!=typeid(double)) wblog(FL,
      "WRN %s() AK not real (%s)",FCT,getName(typeid(TD)).data);
   if (!Rho.isBlockDiagMatrix()) wblog(FL,
      "WRN %s() Rho not blockdiag !??",FCT);


   AK.contract(2,Rho,1,Xk);
   Xk.contract("3,2",AK,"2,3;*",Rho);
};


template <class TQ, class TD>
void nrgUpdateRHO(
   const wbMatrix<double> &R,
   NRGData<TQ,TD> &A,
   NRGData<TQ,TD> &H,
   NRGData<TQ,TD> &RHO,
   unsigned NRho,
   unsigned dloc,
   wbSigHandler *sig
){
   unsigned iter, rzero, N=NRG_N; INDEX_T d; char ie;
   QSpace<gTQ,gTD> Rho;
   wbvector<double> w; double w2, wtot=0;

   const size_t n=16; size_t l=0; char tag[n];

   RHO.A.clearQSpace();

   if (R.dim1!=N || R.dim2!=2)
      wblog(FL, "ERR %s() size mismatch of rhoNorm (%dx%d; %d %d)",
      FCT, R.dim1, R.dim2, N, NRho);

   for (iter=N-1; (int)iter>=0; iter--) {
      if (sig) { sig->call99(); }

      R.getRec(iter,w); rzero=(w<DEPS);
      w2=w.sum();

      l=snprintf(tag,n,
         "%s %2d/%d", !NRho ? "FDM":"DMK", iter,N);
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%s; %d/%d)",FCT,tag,l,n); 

      RHO.updatePara(FL,"A",iter);

      if (rzero && RHO.A.isEmpty()) {
         wblog(FL,"%s skip ...%8s\r\\", tag,"");
         continue;
      }

      A.init(FL,"K",iter);
      A.init(FL,"T",iter);
      H.init(FL,"T",iter);

      if (!A.K.isEmpty()) {
         A.K.getDim(2,&d);
         if (d!=dloc && (iter || !A.T.isEmpty())) {
            sprintf(str,"got inconsistent dloc=%d/%d (k=%d)",d,dloc,iter);
            if (A.T.isEmpty())
                 wblog(FL,"WRN %s() %s",FCT,str);
            else wblog(FL,"ERR %s() %s",FCT,str); 
         }
      }

      ie=RHO.A.isEmpty(); if (!ie)
      updateRho(RHO.A, A.K);

      if (!NRho) {
         if (H.T.isEmpty()) wblog(FL,
            "%s no discarded space (%.3g) %8s\r\\",tag,wtot,"");
         else if (w2!=0) wblog(FL,
            "%s dRho=%.3g (%.3g) %8s\r\\",tag,w2,wtot,"");
         else wblog(FL,
            "%s dRho=%.3g %8s\r\\",tag,wtot,"");
      }
      else wblog(FL,"%-30s\r\\",tag);

      if (!rzero) {
         if (w[1]) {
            initRho(Rho, H.T, iter, w[1]);
            updateRho(Rho, A.T);
            RHO.A+=Rho;
         }
      
         if (w[0]) {
            H.init(FL,"K",iter);
            initRho(Rho, H.K, iter, w[0]);
            updateRho(Rho, A.K);
            RHO.A+=Rho;
         }

         RHO.A.SkipZeroData(DEPS2,'b');
      }

      wtot+=w2;

      if (!RHO.A.isEmpty()) { try {
         checkRho(FL, RHO.A, addRange(R.rec(iter),R.dim2*(N-iter)) );
      }
      catch (...) {
         wblog(FL,"ERR %s() iter=%d/%d\n%s",FCT,iter+1,N,
         iter+1==N?"hint: got kept space at last iteration?":"");
      }}
      else {
         if (ie && rzero) wblog(FL,
            "%s empty RHO (will not contribute) %8s\r\\",tag,"");
         else wblog(FL,"WRN %2d/%d got empty RHO !?? %8s",iter,N,"");
      }
      doflush();
   }
};


#endif

