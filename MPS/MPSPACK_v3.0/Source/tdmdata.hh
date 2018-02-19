#ifndef __WB_TDM_DATA_HCC__
#define __WB_TDM_DATA_HCC__

// switched wbArray (deprecated) => wbarray (col-major)
// Wb,Dec20,10

namespace Wb {
   void Fourier(
      const wbvector<double> &om,
      const wbMatrix<double> &A,
      const wbvector<double> &tt,
      wbMatrix<wbcomplex> &az,
      char domflag=0,
      char lflag=1,
      double alpha=0. // Lorenzian broadening between discrete om values
   );

    void dispSumRule(
       const wbvector<double> &om,
       const wbMatrix<double> &A
    );
}


class TDMData {

 public:

    TDMData() : raw(0), l0(0), fac(1.), check_avg_om(0) {};

    void init (const char *F, int L,
       const char *istr0, const wbvector<double> &t0,
       unsigned Nnrg, unsigned Nop,
       double Emin=1E-8, double Emax=10., unsigned nl=512
    );

    void checkAvgOm() { 
       check_avg_om=1; 
       AW.init(2,AP.dim(2),AP.dim(3));
    };
    wbarray<double> getAvgOm();

    template <class T>
    void Add(
       const wbarray<T> &Om, const wbarray<T> &CC,
       unsigned iter, unsigned s
    );

    void add2RAW(unsigned r, const double*, const double*, unsigned n);
    void addPartialFourier();

    void getOmega(wbvector<double> &om, unsigned &linfit) const;
    void getOmega(wbvector<double> &om) const {
       unsigned linfit=0; getOmega(om,linfit);
    }

    void getSpecData(
       wbvector<double> &om, wbarray<double> &A,
       unsigned linfit=0
    ) const;

    void getSpecData(
       wbvector<double> &om, wbMatrix<double> &A,
       unsigned linfit=0
    );

    void mergeNRGData();

    void getSmoothSpec(
        const double sigma,
        const char bswitch=2
    );

    void Fourier(
       wbvector<double> &om, wbMatrix<wbcomplex> &ST,
       double alpha, double Lambda, char lflag=1
    );

    void Fourier(
       wbvector<double> &om,
       wbMatrix<double> &A,
       const wbvector<double> &tt,
       wbMatrix<wbcomplex> &at,
       double alpha=0.
    );

    void dispSumRule();


    wbarray<double> AP,AN;
    wbMatrix<double> Ap,An;
    wbMatrix<wbcomplex> At;
    wbvector<double> tt;
    wbstring istr;

    mxArray *toMx();
    int init(
       const mxArray *S,
       double &emin, double &emax, unsigned &nlog
    );

    unsigned raw;
    wbvector< wbMatrix<double> > RAW;

 protected:
 private:

    double l0,fac;
    unsigned N;

    char check_avg_om;
    wbarray<double> AW;

    double OM2IDX(const double &om) const {
       return ((log(ABS(om))-l0)*fac);
    }

    unsigned om2idx(const double &om) const {
       if (om==0.) return 0;

       int k = (int)floor(OM2IDX(om));
       if (k<0) k=0; else if (k>=(int)N) k=N-1;

       return (unsigned)k;
    }

    double idx2om(const unsigned &k) const {
       return exp( ((double)k+0.5)/fac + l0 );
    }

    void initLG(
       wbMatrix<double> &wlg,
       unsigned n, double sigma, char bswitch=2
    );

    void foldLogGaussian(
       const wbMatrix<double> &wlg,
       double *aa, double a0, unsigned i0
    );
};


void TDMData::init(const char *F, int L, const char *istr0,
   const wbvector<double> &t0,
   unsigned Nnrg, unsigned Nop,
   double Emin, double Emax, unsigned nl
){
   double l2;

   if (Emin<=0. || Emax<=0. || Emin>=0.1*Emax || nl==0) wblog(F,L,
   "ERR Invalid parameter set [%g %g %d]%s.", Emin, Emax, nl,
   (Emin>0 && Emax>0 && Emin<Emax) ? " (require Emin<Emax/10)":"");

   l0 = log(Emin);
   l2 = log(Emax);

   fac= log10(M_E)*double(nl);

   N = (unsigned)((l2-l0)*fac)+1;

   AP.init(N,Nop,Nnrg);
   AN.init(N,Nop,Nnrg);

   if (check_avg_om) AW.init(2,Nop,Nnrg);

   sprintf(str,"tdmData(%s)", istr0 ? istr0 : "");
   wblog(F,L,"<i> %12s: %.2g .. %.2g (%d/dec; 2x%d)",
   str, Emin, Emax, nl, N);

   istr=istr0; tt=t0;
};


wbarray<double> TDMData::getAvgOm() {

   if (AW.rank()!=3 || AW.dim(1)!=2) wblog(FL,
      "ERR invalid AW data structure [%s] !??",AW.sizeStr().data);

   unsigned i,n=AW.numel();
   wbarray<double> aw(AW); if (aw.isEmpty()) return aw;
   double *a=AW.data, *b=aw.data;

   for (i=0; i<n; i+=2)
   b[i]=((a[i+1]!=0) ? a[i]/a[i+1] : a[i]);

   return aw;
};


template <class T>
void TDMData::Add(
   const wbarray<T> &Om,
   const wbarray<T> &CC,
   unsigned iter,
   unsigned s
){
   static int iter_last=-1;
   unsigned i,k,n=Om.numel();
   T *ap, *an;

   if (!Om.hasSameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR size mismatch between Om and CC!");
   }

   if (s>=AP.dim(2) || iter>=AP.dim(3)) wblog(FL,
      "ERR index out of bounds (op=%d/%d; k=%d/%d)",
       s+1,AP.dim(2),iter+1,AP.dim(3));

   ap=&AP(0,s,iter);
   an=&AN(0,s,iter);

   for (i=0; i<n; i++) {
       k=om2idx(Om.data[i]);
       if (Om.data[i]<0.)
            an[k]+=CC.data[i];
       else ap[k]+=CC.data[i];
   }

   if (raw) {
      if (iter_last!=(int)iter) {
wblog(FL,"%NTST iter=%d",iter);
         if (RAW.len) addPartialFourier();
         RAW.initDef(AP.dim(2));
wblog(FL,"TST done");
      } 

      add2RAW(s, Om.data, CC.data, n);

      iter_last=iter;
   }

   if (check_avg_om) {
      if (AW.rank()!=3 || AW.dim(1)!=2 || 
          AW.dim(2)!=AP.dim(2) || AW.dim(3)!=AP.dim(3)) wblog(FL,
         "ERR size mismatch for AW [%s] ? [%s]",
          AW.sizeStr().data, AP.sizeStr().data
      );

      AW(0,s,iter)+=Om.sumTimesEl(CC);
      AW(1,s,iter)+=CC.sum();
   }
};


inline void TDMData::add2RAW(
   unsigned r, const double* om, const double* cc, unsigned n
){
   unsigned k, s=n*sizeof(double);
   wbMatrix<double> &R=RAW[r];

   k=R.dim2; R.Resize(2,k+n);

   memcpy(R.ref(0,k), om, s);
   memcpy(R.ref(1,k), cc, s);
};


inline void TDMData::addPartialFourier() {

   unsigned s,n,m=RAW.len;
   wbMatrix<wbcomplex> at;
   wbvector<double> om;
   wbMatrix<double> aa;

   if (!RAW.len) return;

   if (m!=AP.dim(2)) wblog(FL,
      "ERR severe size mismatch (%d/%d)",m,AP.dim(2));

   if (At.isEmpty()) At.init(AP.dim(2),tt.len);

   for (s=0; s<m; s++) {
      wbMatrix<double> &R=RAW[s]; n=R.dim2;
      if (R.dim1!=2) {
         wblog(FL,"WRN R.dim1=%d !?? - skip.", R.dim1);
         continue;
      }

      om.init2ref(          n,R.rec(0));
      aa.init2ref(R.dim1-1, n,R.rec(1));

      Wb::Fourier(om,aa,tt,at);
      At.recAddP(s,at.data);
   }

   RAW.init();
}


void TDMData::getOmega(
   wbvector<double> &om, unsigned &linfit
) const {

   unsigned i;
   double dbl;

   if (linfit) linfit=unsigned(fac);
   if (linfit>N) wblog(FL,"ERR linfit>N (%d,%d) ???", N, linfit);

   om.init(2*N);

   dbl = idx2om(linfit) / (linfit+0.5);

   for (i=0; i<linfit; i++) om[N-1-i] = -(om[N+i] = dbl*(i+0.5));
   for (   ; i<N;      i++) om[N-1-i] = -(om[N+i] = idx2om(i)  );
}


void TDMData::getSpecData(
   wbvector<double> &om,
   wbarray<double> &A,
   unsigned linfit
) const {
   unsigned i,r,s, d1=AP.dim(2), d2=AP.dim(3);

   getOmega(om,linfit);
   A.init(2*N,d1,d2);

   if (om.len!=2*N || !AP.hasSameSize(AN)) wblog(FL,
      "ERR severe dimension mismatch (%d,%d).",2*N,om.len);

   for (r=0; r<d1; r++)
   for (s=0; s<d2; s++) {
      for (i=0; i<N; i++) {
         A(N-i-1,r,s) = AN(i,r,s);
         A(N+i  ,r,s) = AP(i,r,s);
      }
   }

   if (linfit) {
      for (r=0; r<d1; r++)
      for (s=0; s<d2; s++) {
          set2avg(A.ref(r,s)+N,        linfit);
          set2avg(A.ref(r,s)+N-linfit, linfit);
      }
   }
};


void TDMData::getSpecData(
   wbvector<double> &om,
   wbMatrix<double> &A,
   unsigned linfit
){

   unsigned i,r;

   getOmega(om,linfit);

   if (Ap.isEmpty()) mergeNRGData();
   A.init(Ap.dim1,2*N);

   if (A.dim2!=om.len || !Ap.hasSameSize(An)) wblog(FL,
      "ERR Severe dimension mismatch (%d,%d).",A.dim2,om.len);

   for (r=0; r<A.dim1; r++) 
   for (i=0; i<N; i++) {
      A(r,N-i-1) = An(r,i);
      A(r,N+i  ) = Ap(r,i);
   }

   if (linfit) {
      for (r=0; r<A.dim1; r++) {
          set2avg(A.rec(r)+N,        linfit);
          set2avg(A.rec(r)+N-linfit, linfit);
      }
   }
}


void TDMData::Fourier(
   wbvector<double> &om,
   wbMatrix<wbcomplex> &ST,
   double alpha, double Lambda,
   char lflag
){
   unsigned k,it,ik,is,io,l=0,no, nt=tt.len, ns=AP.dim(2), nk=AP.dim(3);

   wbarray<wbcomplex> TKS(nt,nk,ns), KST;
   wbMatrix<double>   aTK(nt,nk);
   wbarray<double> AR;
   double ek,*ap,kso,o;
   wbcomplex z,*st;

   if (lflag)
   wblog(FL,"<i> %s (damping alpha=%g, Lambda=%g)", FCT, alpha, Lambda);

   getSpecData(om,AR);
   ap=AR.data; no=om.len;

   if (alpha==0) aTK.set(1); else
   for (ik=0; ik<nk; ik++) {
      ek=pow(Lambda, -double(ik)/2.) * 0.5*(Lambda+1);
      for (it=0; it<nt; it++)
      aTK(it,ik)=exp(-fabs(alpha*ek*tt[it]));
   }

   for (ik=0; ik<nk; ik++)
   for (is=0; is<ns; is++)
   for (io=0; io<no; io++, l++) if (ap[l]!=0) { kso=ap[l]; o=om[io];
       for (it=0; it<nt; it++)
       TKS(it,ik,is)+=kso*expi(o*tt[it]);
   }

   TKS.permute(KST,"2 3 1");
   ST.init(ns,nt);

   for (is=0; is<ns; is++)
   for (it=0; it<nt; it++) {
      st=KST.ref(is,it); z=0;
      if (alpha==0)
           for (k=0; k<nk; k++) z+=st[k];
      else for (ap=aTK.rec(it), k=0; k<nk; k++) z+=ap[k]*st[k];
      ST(is,it)=z;
   }

};


void TDMData::Fourier(
   wbvector<double> &om,
   wbMatrix<double> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &at,
   double alpha
){
   if (Ap.isEmpty()) { mergeNRGData(); dispSumRule(); }
   getSpecData(om,A);

   Wb::Fourier(om,A,tt,at,0,1,alpha);
};


void Wb::Fourier(
   const wbvector<double> &om,
   const wbMatrix<double> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &az,
   char domflag,
   char lflag,
   double alpha
){
   unsigned i,j,l,m,n,b1,b2, nom, nt=tt.len, D=512, stop=0;
   wbMatrix<double> ez1,ez2;
   wbMatrix<double> E1,E2;
   wbcomplex z;
   double f,ot;

   wbvector< wbMatrix<wbcomplex> > AZ;
   wbvector<double> dom(om.len);

   m=A.dim1; nom=om.len;

   if (nom!=A.dim2) wblog(FL,
      "ERR %s - size mismatch (%d/%d)",nom,A.dim2);
   if (!tt.len || !m || !nom) {
      wblog(FL,"WRN Empty data [%d %d %d] - return",tt.len,m,nom);
      return;
   }
   if (nom<2) wblog(FL,"ERR %s - omega too short (%d/2)",nom);

   n=om.len-1;
   for (i=1; i<n; i++) dom[i]=0.5*(om[i+1]-om[i-1]);
   dom[0]=om[1]-om[0]; dom[n]=om[n]-om[n-1];

   for (i=0; i<dom.len; i++) if (dom[i]<=0) { wblog(FL,
   "WRN omega not strictly ascending !??"); break; }

   for (i=0; i<tt.len; i++) if (tt[i]<0) { wblog(FL,
   "WRN negative tt values !??"); break; }

   if (domflag) {
      wbMatrix<double> B(A);

      for (i=0; i<B.dim1; i++)
      for (j=0; j<B.dim2; j++) B(i,j)*=dom[j];

      Wb::Fourier(om,B,tt,az,0,lflag,alpha);
      return;
   }

   if (lflag) wblog(FL,
   "<i> %s of (broadened) spectral data (alpha=%.g)", FCT, alpha);


   AZ.initDef(nt<D ? 1 : nt/D);

   for (l=b1=0; !stop; l++,b1+=D) {
      b2=b1+D-1; if (b2+D>=nt) { b2=nt-1; stop=1; }
      n=b2-b1+1;

      if (n*nom*sizeof(double)/(1<<20)>900) {
         wblog(FL,"WRN trying to allocate %.6g MB !?", 
         n*nom*sizeof(double)/double(1<<20));
      }

      ez1.init(nom,n); ez2.init(nom,n);
      for (i=0; i<nom; i++)
      for (j=0; j<n; j++) { ot=om[i]*tt[b1+j];
         ez1(i,j)=cos(ot);
         ez2(i,j)=sin(ot);

         if (alpha>0) {
            f=exp(-fabs(alpha*dom[i]*tt[b1+j]));
            ez1(i,j)*=f;
            ez2(i,j)*=f;
         }
      }
      wbMatProd(A,ez1,E1);
      wbMatProd(A,ez2,E2); AZ[l].set(E1,E2);
   }

   az.CAT(2,AZ);
}


void TDMData::mergeNRGData() {

   unsigned i,k, d1=AP.dim(1), d2=AP.dim(2), d3=AP.dim(3);
   unsigned s=d1*d2, S=s*d3;
   double *Dp,*Dn,*dp,*dn;

   if (!AP.hasSameSize(AN)) wblog(FL, "ERR severe size mismatch!");

   Ap.init(d2,d1); dp=Ap.data; Dp=AP.data;
   An.init(d2,d1); dn=An.data; Dn=AN.data;

   for (i=0; i<S; i++) { k=i%s;
      dp[k]+=Dp[i];
      dn[k]+=Dn[i];
   }
};


void TDMData::dispSumRule() {

   unsigned i=0; double m=0;
   wbvector<double> sp,sn;
   wbarray<double> x;

   AP.sum("1,3",x); sp.init(x.numel(),x.data);
   AN.sum("1,3",x); sn.init(x.numel(),x.data);

   wblog(FL," *  om>0: %s\n *  om<0: %s",
   sp.toStrf("%8.5f").data, sn.toStrf("%8.5f").data);

   sp+=sn; 
   
   if (sp.len) {
      m=round(sp[0]);
      for (i=0; i<sp.len; i++) if (ABS(sp[i]-m)>1E-2) break;
   }

   if (i==sp.len) { sp-=m; wblog(FL,
   " *  \e[34mtotal %s\e[0m (%g)", sp.toStrf("%8.2g").data, m); }
   else wblog(FL,
   " *  \e[34mtotal %s\e[0m",      sp.toStrf("%8.5f").data);
}

void Wb::dispSumRule(
   const wbvector<double> &om,
   const wbMatrix<double> &A
){
   wbvector<double> s(A.dim1);
   unsigned i; double m=0;

   if (om.len!=A.dim2) wblog(FL,
   "ERR Dimension mismatch (%d,%d)", om.len, A.dim2);

   for (unsigned i=0; i<A.dim1; i++)
   s[i] = wbIntTrapez(om.data, A.rec(i), om.len);

   wblog(FL,"SUM %s", s.toStrf("%8.5f").data);

   if (s.len) {
      for (m=round(s[0]), i=1; i<s.len; i++)
      if (ABS(s[i]-m)>1E-2) return;
   }

   s-=m; wblog(FL,
   "\e[34m ·  total-%g     = %s\e[0m", m, s.toStrf("%8.2g").data);
}


void TDMData::getSmoothSpec(const double sigma, const char bswitch) {

    unsigned i,s;
    wbMatrix<double> Ap0(Ap), An0(An), wlg;

    if (fac<=0. || sigma<=0. || N<4)
    wblog(FL,"ERR fac=%g, sigma=%g, N=%d ???",fac,sigma,N);

    if (fac*sigma<1.) wblog(FL,
    "WRN sigma<discrete log spacing (%.3g<%.3g!)",sigma,1/fac);

    wblog(FL,"<i> %s: %dx%d bins (sigma=%g)", FCT, AP.dim(2), N, sigma);

    if (Ap.isEmpty()) mergeNRGData(); dispSumRule();
    Ap.set(0);
    An.set(0);

    initLG(wlg,N,sigma,bswitch);

    for (s=0; s<Ap.dim1; s++)
    for (i=0; i<Ap.dim2; i++) {
       foldLogGaussian(wlg, Ap.rec(s), Ap0(s,i), i);
       foldLogGaussian(wlg, An.rec(s), An0(s,i), i);
    }

};


void TDMData::initLG(
   wbMatrix<double> &wlg,
   unsigned n, double sigma, char bswitch
){
   unsigned i,j;
   double dbl,efac,s=0;

   if (sigma<=0. || fac<=0.) wblog(FL,
   "ERR sigma=%g, fac=%g ???",sigma,fac);

   efac=1./(fac*sigma); wlg.init(2,n);

   switch (bswitch) {
      case 0:
          strcpy(str,"symmetric in log. space");
          s=0; break;
      case 1:
          strcpy(str,"symmetric log-gauss");
          s=sigma/2.;
          break;
      case 2:
          strcpy(str,"norm preserving log-gauss");
          s=sigma/4.;
          break;
      default:
      wblog(FL,"ERR Invalid bswitch=%d ???",bswitch);
   }
   wblog(FL,"<i> %s - %s", FCT, str);

   for (i=0; i<n; i++) {
       dbl=efac*double(i)+s; wlg(0,i)=exp(-dbl*dbl);
       dbl=efac*double(i)-s; wlg(1,i)=exp(-dbl*dbl);
   }

   wlg*=(1./wlg.sum());

   for (i=0; i<2; i++)
   for (j=0; j<n; j++) if (fabs(wlg(i,j))<1E-20) wlg(i,j)=0;
}


inline void TDMData::foldLogGaussian(
   const wbMatrix<double> &wlg,
   double *aa, double a0, unsigned i0
){
   unsigned i, n=wlg.dim2;
   const double *wp=wlg.rec(0), *wn=wlg.rec(1);
   double wi,nrm=0.;

   for (i=0; i<n; i++)
   nrm += (i0>i ? wp[i0-i] : wn[i-i0]);
   nrm=a0/nrm;

   for (i=0; i<n; i++) {
      wi=(i0>i ? wp[i0-i] : wn[i-i0]);
      if (wi) aa[i]+=nrm*wi;
   }
};


mxArray* TDMData::toMx() {
   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,S,"name",istr.toMx());
   mxAddField2Scalar(FL,S,"AP", AP.toMx());
   mxAddField2Scalar(FL,S,"AN", AN.toMx());
   mxAddField2Scalar(FL,S,"fac",numtoMx(fac));
   mxAddField2Scalar(FL,S,"l0", numtoMx(l0 ));

   if (raw) {
      unsigned i,n=RAW.len;
      mxArray *R=mxCreateCellMatrix(1,n);

      for (i=0; i<n; i++)
      mxSetCell(R,i,RAW[i].toMx('r'));

      mxAddField2Scalar(FL,S,"RAW",R);
   }

   return S;
};

int TDMData::init(const mxArray *S,
   double &emin, double &emax, unsigned &nlog
){
   if (!mxIsStruct(S) || !mxIsScalar(S) ||
        mxGetFieldNumber(S,"name")<0 ||
        mxGetFieldNumber(S,"AP"  )<0 ||
        mxGetFieldNumber(S,"AN"  )<0 ||
        mxGetFieldNumber(S,"fac" )<0 ||
        mxGetFieldNumber(S,"l0"  )<0
   ) return 1;

   istr.init(mxGetField(S,0,"name"  ));
   AP  .init(mxGetField(S,0,"AP"));
   AN  .init(mxGetField(S,0,"AN"));
   mxGetNumber(mxGetField(S,0,"fac"),fac);
   mxGetNumber(mxGetField(S,0,"l0"), l0 );

   if (AP.SIZE!=AN.SIZE || AP.isEmpty() || AN.isEmpty()) return 2;
   if (fac<=0) return 3;

   N=AP.dim(1);

   emin=exp(l0);
   emax=idx2om(N);
   nlog=(unsigned)(fac/log10(M_E) + 0.5);

   wblog(FL,
     "<i> TDMdata `%s' from RAW data: %dx%d bins [x%d]",
      istr.data, AP.dim(1), AP.dim(2), AP.dim(3));

   return 0;
};


#endif

