#ifndef __WB_TDMDATA_HCC__
#define __WB_TDMDATA_HCC__

namespace Wb {
   void Fourier(
      const wbvector<double> &om,
      const wbMatrix<double> &A,
      const wbvector<double> &tt,
      wbMatrix<wbcomplex> &az,
      char domflag=0
   );

    void dispSumRule(
       const wbvector<double> &om,
       const wbMatrix<double> &A
    );
}


class TDMData {

 public:

    TDMData() : l0(0), fac(1.) {};

    void init (const char *F, int L,
       const char *istr0, const wbvector<double> &t0,
       unsigned Nnrg, unsigned Nop,
       double Emin=1E-8, double Emax=10., unsigned nl=512
    );

    template <class T>
    void Add(
       const wbArray<T> &Om, const wbArray<T> &CC,
       unsigned iter, unsigned s
    );

    void getOmega(wbvector<double> &om, unsigned &linfit) const;
    void getOmega(wbvector<double> &om) const {
       unsigned linfit=0; getOmega(om,linfit);
    }

    void getRawData(
       wbvector<double> &om, wbArray<double> &A,
       unsigned linfit=0
    ) const;

    void getRawData(
       wbvector<double> &om, wbMatrix<double> &A,
       unsigned linfit=0
    );

    wbMatrix<wbcomplex>& Fourier(
       wbvector<double> &om, wbMatrix<wbcomplex> &ST,
       double alpha, double Lambda, char lflag=1
    );

    void mergeNRGData();

    void getSmoothTDM(wbMatrix<double> &A, const wbvector<double> &om,
    const double eps, const double sigma=0.6, const char iflag=1);

    void TDMData::getSmoothTDM(
       wbvector<double> &om, wbMatrix<double> &A,
       const double eps, const double sigma=0.6, const char iflag=1
    ){ unsigned linfit=0;
       getOmega(om,linfit);
       getSmoothTDM(A,om,eps,sigma,iflag);
    };

    void dispSumRule();


    wbArray<double> AP,AN;
    wbMatrix<double> Ap,An;
    wbMatrix<wbcomplex> At;
    wbvector<double> tt;
    wbstring istr;

 protected:
 private:

    double l0,fac;
    unsigned N;

    wbMatrix<double> wlgauss;

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

    void initLogGaussian(unsigned n, double sigma);
    double foldLogGaussian(double om, const double *aa);

    double foldGaussian(
       const double om, const double sigma,
       const double *oo, const double *aa, const unsigned n
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

   AP.init(Nnrg,Nop,N);
   AN.init(Nnrg,Nop,N);

   sprintf(str,"tdmData(%s)", istr0 ? istr0 : "");
   wblog(F,L,"<i> %12s: %.2g .. %.2g (%d/dec; 2x%d)",
   str, Emin, Emax, nl, N);

   istr=istr0; tt=t0;
};


template <class T>
void TDMData::Add(
   const wbArray<T> &Om,
   const wbArray<T> &CC,
   unsigned iter,
   unsigned s
){
   unsigned i,k,n=Om.SIZE.prod();
   T *ap, *an;

   if (!Om.hasSameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR Size mismatch between Om and CC!");
   }

   if (iter>=AP.dim(0) || s>=AP.dim(1)) wblog(FL,
   "ERR Index out of bounds (%d/%d; %d/%d)",iter,AP.dim(0),s,AP.dim(1));

   ap=&AP(iter,s,0);
   an=&AN(iter,s,0);

   for (i=0; i<n; i++) {
       k=om2idx(Om.data[i]);
       if (Om.data[i]<0.)
            an[k]+=CC.data[i];
       else ap[k]+=CC.data[i];
   }
};


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


void TDMData::getRawData(
   wbvector<double> &om,
   wbArray<double> &A,
   unsigned linfit
) const {
   unsigned i,r,s, d1=AP.dim(0), d2=AP.dim(1);

   getOmega(om,linfit);

   A.init(d1,d2,2*N);

   if (A.dim(2)!=om.len || !AP.hasSameSize(AN)) wblog(FL,
   "ERR Severe dimension mismatch (%d,%d).",A.dim(2),om.len);

   for (r=0; r<d1; r++)
   for (s=0; s<d2; s++) {
      for (i=0; i<N; i++) {
         A(r,s,N-i-1) = AN(r,s,i);
         A(r,s,N+i  ) = AP(r,s,i);
      }
   }

   if (linfit) {
      for (r=0; r<d1; r++)
      for (s=0; s<d2; s++) {
          set2avg(A.ref(r,s)+N,        linfit);
          set2avg(A.ref(r,s)+N-linfit, linfit);
      }
   }
}


void TDMData::getRawData(
   wbvector<double> &om,
   wbMatrix<double> &A,
   unsigned linfit
){

   unsigned i,r;

   getOmega(om,linfit);

   mergeNRGData();
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


wbMatrix<wbcomplex>& TDMData::Fourier(
   wbvector<double> &om,
   wbMatrix<wbcomplex> &ST,
   double alpha, double Lambda,
   char lflag
){
   unsigned k,it,ik,is,io,l=0,no, nt=tt.len, nk=AP.dim(0), ns=AP.dim(1);

   wbArray<wbcomplex> KST(nk,ns,nt), STK;
   wbMatrix<double>   ETK(nt,nk);
   wbArray<double> AR;
   double ek,*ap,kso,o;
   wbcomplex z,*st;

   if (lflag)
   wblog(FL,"<i> %s (alpha=%g, Lambda=%g)", FCT, alpha, Lambda);

   getRawData(om,AR);
   ap=AR.data; no=om.len;

   if (alpha==0) ETK.set(1); else
   for (ik=0; ik<nk; ik++) {
      ek=pow(Lambda, -double(ik)/2.) * 0.5*(Lambda+1);
      for (it=0; it<nt; it++)
      ETK(it,ik)=exp(-fabs(alpha*ek*tt[it]));
   }

   for (ik=0; ik<nk; ik++)
   for (is=0; is<ns; is++)
   for (io=0; io<no; io++, l++) if (ap[l]!=0) { kso=ap[l]; o=om[io];
       for (it=0; it<nt; it++)
       KST(ik,is,it)+=kso*expi(o*tt[it]);
   }

   KST.permute(STK,"2 3 1");
   ST.init(ns,nt);

   for (is=0; is<ns; is++)
   for (it=0; it<nt; it++) {
      st=STK.ref(is,it); z=0;
      if (alpha==0)
           for (k=0; k<nk; k++) z+=st[k];
      else for (ap=ETK.rec(it), k=0; k<nk; k++) z+=ap[k]*st[k];
      ST(is,it)=z;
   }

wblog(FL,"TST [nk ns nt no]=[%d %d %d %d]", nk, ns, nt, no);
   ST.save2(At);
ST.info("ST"); At.info("At");

   return At;
}


void TDMData::mergeNRGData() {

   unsigned i,k, d1=AP.dim(0), d2=AP.dim(1), d3=AP.dim(2);
   unsigned s=d2*d3, S=d1*d2*d3;
   double *Dp,*Dn,*dp,*dn;

   if (!AP.hasSameSize(AN)) wblog(FL, "ERR Severe size mismatch!");

   Ap.init(d2,d3); dp=Ap.data; Dp=AP.data;
   An.init(d2,d3); dn=An.data; Dn=AN.data;

   for (i=0; i<S; i++) { k=i%s;
      dp[k]+=Dp[i];
      dn[k]+=Dn[i];
   }
}


void TDMData::dispSumRule() {

   unsigned i=0; double m=0;
   wbvector<double> sp,sn;
   wbArray<double> A;

   AP.sum("1 3",A); sp.init(A.SIZE.prod(),A.data);
   AN.sum("1 3",A); sn.init(A.SIZE.prod(),A.data);

   wblog(FL,"SUM %s\nw>0 %s\nw<0 %s",
   istr.data, sp.toStrf("%8.5f").data, sn.toStrf("%8.5f").data);

   sp+=sn; 
   
   if (sp.len) {
      m=round(sp[0]);
      for (i=0; i<sp.len; i++) if (ABS(sp[i]-m)>1E-2) break;
   }

   if (i==sp.len) { sp-=m; wblog(FL,
   "\e[34mtot %s\e[0m (%g)", sp.toStrf("%8.2g").data, m); }
   else wblog(FL,
   "\e[34mtot %s\e[0m",      sp.toStrf("%8.5f").data);
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


void TDMData::getSmoothTDM(
    wbMatrix<double> &A,
    const wbvector<double> &omega,
    const double eps,
    const double sigma,
    const char iflag
){
    unsigned i,s;
    wbvector<double> OM(N),I;
    double x, dbl, w1=0., w2=0., om, aom;

    if (fac<=0. || sigma<=0. || eps<=0. || N<4)
    wblog(FL,"ERR fac=%g, sigma=%g, eps=%g, N=%d ???",fac,sigma,eps,N);

    for (i=0; i<OM.len; i++) OM[i]=idx2om(i);

    if (iflag) {
       wblog(FL,
       "<i> %s: 2x%d bins, %g/dec, sigma=%g",
            FCT, N, fac/log10(M_E), sigma);
       wblog(FL,
       " ·  omega range   = %.3g .. %.3g .. %.3g\n"
       " ·  emin/eps/emax = %.3g .. %g .. %.3g",
            omega.min(), omega.aMin(), omega.max(),
            idx2om(0), eps, idx2om(N));
       dispSumRule();
    }
    else {
       wblog(FL,
         "<i> %s: 2x%d bins\n"
         "    om: %.3g .. %.3g .. %.3g (%d)",
       FCT, N, omega.min(), omega.aMin(), omega.max(), omega.len);
    }

    if (Ap.isEmpty()) mergeNRGData();

    A.init(Ap.dim1, omega.len);

    initLogGaussian(N,sigma);

    for (s=0; s<Ap.dim1; s++)
    for (i=0; i<omega.len; i++) {
        om=omega[i]; aom=ABS(om);

        if (aom>eps) x=1.;
        else {
            if (aom==0.) x=0.;
            else {
               dbl = (log(aom)-log(eps))/sigma;
               x = exp(-dbl*dbl);
               if (x<1E-10) x=0.;
            }
        }

        if (x!=0.)
        w1 = foldLogGaussian(aom, om>0 ? Ap.rec(s):An.rec(s));

        if (x!=1.)
        w2 = foldGaussian( om, 0.8*eps, OM.data, Ap.rec(s), N)
           + foldGaussian(-om, 0.8*eps, OM.data, An.rec(s), N);

        A(s,i) = x*w1 + (1.-x)*w2;
    }

    Wb::dispSumRule(omega,A);
};


void TDMData::initLogGaussian(unsigned n, double sigma){

   unsigned i;
   double dbl,efac,nrm,s;

   if (sigma<=0. || fac<=0.) wblog(FL,
   "ERR sigma=%g, fac=%g ???",sigma,fac);
   wblog(FL,"TST Initializing log gauss");

   efac=1./(fac*sigma); wlgauss.init(2,n);

   if (0) {
      s=sigma/2.;
      nrm = exp(+s*s) / (sigma*sqrt(M_PI));
   }
   else {
      s=sigma/4.;
      nrm = 1./(sigma*sqrt(M_PI));
   }

   for (i=0; i<n; i++) {
       dbl=efac*double(i)-s; wlgauss(0,i)=exp(-dbl*dbl);
       dbl=efac*double(i)+s; wlgauss(1,i)=exp(-dbl*dbl);
   }


   wlgauss *= nrm;

};

inline double TDMData::foldLogGaussian(double om, const double *aa){

   unsigned i,k,n=wlgauss.dim2; double xa=0;
   const double *wp=wlgauss.rec(0), *wn=wlgauss.rec(1);

   for (k=om2idx(om), i=0; i<n; i++)
   xa += aa[i]*(k>i ? wp[k-i] : wn[i-k]);

   return xa/om;
};


double TDMData::foldGaussian(
   const double om, const double sigma,
   const double *oo, const double *aa, const unsigned n
){
   double xa=0., dbl, efac, wtot=0;

   if (sigma<=0.) wblog(FL, "ERR sigma=%g ???", sigma);
   efac=1./sigma;

   for (unsigned i=0; i<n; i++) {
       dbl=efac*(om-oo[i]); dbl=exp(-dbl*dbl); wtot+=dbl;
       xa += aa[i]*dbl;
   }

   xa /= wtot;

   return xa;
};


void Wb::Fourier(
   const wbvector<double> &om,
   const wbMatrix<double> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &az,
   char domflag
){
   unsigned i,j,l,m,n,b1,b2, nom, nt=tt.len, D=128, stop=0;
   wbMatrix<double> ez1,ez2;
   wbMatrix<double> E1,E2;
   wbcomplex z;
   double ot;

   wbvector< wbMatrix<wbcomplex> > AZ;

   m=A.dim1; nom=om.len;

   if (nom!=A.dim2) wblog(FL,
      "ERR %s - size mismatch (%d/%d)",nom,A.dim2);
   if (!tt.len || !m || !nom) {
      wblog(FL,"WRN Empty data [%d %d %d] - return",tt.len,m,nom);
      return;
   }

   if (domflag) {
      wbMatrix<double> B(A);
      wbvector<double> dom(om.len); n=om.len-1;

      for (i=1; i<om.len; i++) dom[i]=0.5*(om[i+1]-om[i-1]);
      dom[0]=dom[1]; dom[n]=dom[n-1];

      for (i=0; i<B.dim1; i++)
      for (j=0; j<B.dim2; j++) B(i,j)*=dom[j];

      Wb::Fourier(om,B,tt,az,0);
      return;
   }


   AZ.initDef(nom/D);

   for (l=b1=0; !stop; l++,b1+=D) {
      b2=b1+D-1; if (b2+D>=nt) { b2=nt-1; stop=1; }
      n=b2-b1+1; ez1.init(nom,n); ez2.init(nom,n);
      for (i=0; i<nom; i++)
      for (j=0; j<n; j++) { ot=om[i]*tt[b1+j];
         ez1(i,j)=cos(ot);
         ez2(i,j)=sin(ot);
      }
      wbMatProd(A,ez1,E1);
      wbMatProd(A,ez2,E2); AZ[l].set(E1,E2);
   }

   az.CAT(2,AZ);
}


#endif

