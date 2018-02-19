#ifndef __WB_SPECTRAL_HCC__
#define __WB_SPECTRAL_HCC__

/* ---------------------------------------------------------------- *
 * class Spectral
 * AWb C May 2006
 * ---------------------------------------------------------------- */

void KKimag(
   const wbMatrix<double> &Rp,
   const wbMatrix<double> &Rn,
   const wbvector<double> &OM,
   wbMatrix<double> &Ip,
   wbMatrix<double> &In
);

class Spectral {
 public:

    Spectral() : l0(0), fac(1.), Delta(0.) {};

    void init () {
        Ap.init(); An.init(); Ap_buf.init(); An_buf.init();
        Ar.init(); A0.init(); mspec.init();
        istr.init();
    };

    void init (const char *i, unsigned Nops,
       double Emin=1E-8, double Emax=10., unsigned nl=512,
       double delta=0., double T_=0, unsigned Niter=0
    );

    void initBuf() {
       if (Ap.dim1!=An.dim1) wblog(FL,
       "ERR Severe dimension mismatch (%d,%d).", Ap.dim1, An.dim1);

       Ap_buf.init(Ap.dim1, Ap.dim2);
       An_buf.init(An.dim1, An.dim2);
    };

    void calcSpectralMoments(const char *F, int L, unsigned n) {
       if (n) {
          if (n>4) wblog(F,L,"WRN %s with n=%d !??",FCT,n);
          else wblog(F,L,"<i> calculate %d spectral moment(s)",n);
          mspec.init(Ap.dim1,n+1);
       }
    };

    template <class T>
    void Add(
       const wbarray<T> &Om, const wbarray<T> &CC, unsigned s,
       char dblock=0,
       char mflag=0
    );

    template <class T>
    void Add2buf(const wbarray<T> &Om, const wbarray<T> &CC, unsigned s);

    void crossAddBuf(const char minor=1, char iflag=1);

    void saveIter(const char *F, int L, unsigned iter);

    void getOmega(wbvector<double> &om, unsigned &linfit) const;
    void getOmega(wbvector<double> &om) const {
       unsigned linfit=0;
       getOmega(om,linfit);
    };

    template <class TQ, class TD>
    void applyIROPfac(
       const wbvector<QSpace<TQ,TD> > &F, wbvector<double> *fac_=NULL
    );

    void getRawData(
       wbvector<double> &om, wbMatrix<double> &A,
       unsigned linfit=0, const char *isbuf=""
    ) const;

    void getPARTIAL(
       wbarray<double> &A1, wbarray<double> &A2) const;

    void pairRawData();

    void detailedBalance(
    const double T, wbvector<double> &om, wbMatrix<double> &aa) const;

    void dispSumRule();
    void dispSumRule(const wbvector<double> &om, const wbMatrix<double> &A);

    void getSmoothSpectral(
       wbvector<double> &om, wbMatrix<double> &A,
       double eps, double sigma=0.6, double sigma2=0.5, char iflag=1
    );

    void getSmoothSpectral(
       wbMatrix<double> &A, const wbvector<double> &om,
       double eps, double sigma=0.6, double sigma2=0.5, char iflag=1
    );

    void KKreal(
       wbMatrix<wbcomplex> &GZ,
       wbvector<double> &OM,
       wbvector<double> &DOM
    ) const;


    wbMatrix<double> Ap, An;
    wbMatrix<double> Ap_buf, An_buf;
    wbvector<double> Ar;
    wbMatrix<double> A0;
    wbMatrix<double> mspec;

    wbarray<double> AP, AN;

    wbstring istr;
    double Temp;

 protected:
 private:

    double l0, fac, Delta;

    double OM2IDX(const double &om) const {
       return ((log(ABS(om))-l0)*fac);
    }

    unsigned om2idx(double om) const {
       if (Delta==0) { if (om==0.) return 0; }
       else {
          if (om> Delta) om-=Delta; else
          if (om<-Delta) om+=Delta; else return 0;
       }

       int k = (int)floor(OM2IDX(om));
       if (k<0) k=0; else if (k>=(int)Ap.dim2) k=Ap.dim2-1;

       return (unsigned)k;
    }

    double idx2om(const unsigned &k) const {
       return exp( ((double)k+0.5)/fac + l0 );
    }

    int crossAddBuf_aux(
       double *a, const double *b, const unsigned M,
       const char minor=1, char iflag=1
    );

    double foldLogGaussian(
       const double om, const double sigma,
       const double *aa, const unsigned n
    );

    double foldGaussian(
       const double om, const double sigma,
       const double *oo, const double *aa, const unsigned n
    );
};


void Spectral::init (const char *istr0,
   unsigned Nops,
   double Emin, double Emax, unsigned nl,
   double d,
   double T,
   unsigned Niter
){
   unsigned Nom;
   double l2;

   if (Emin<=0. || Emax<=0. || Emin>=0.1*Emax || nl==0) wblog(FL,
      "ERR Invalid parameter set [%g %g %d]%s.", Emin, Emax, nl,
      (Emin>0 && Emax>0 && Emin<Emax) ? " (require Emin<Emax/10)":""
   );

   if (d<0) wblog(FL,"ERR Invalid Delta (%g neg!)",d);
   Delta=d;

   l0 = log(Emin);
   l2 = log(Emax);

   fac= log10(M_E)*double(nl);

   Nom = (unsigned)((l2-l0)*fac)+1;

   Ap.init(Nops,Nom);
   An.init(Nops,Nom); Ar.init(Nops);

   if (int(Niter)>0) {
      AP.init(Nom,Nops,Niter);
      AN.init(Nom,Nops,Niter);
   }

   istr=istr0; Temp=T;
};


template <class T>
void Spectral::Add(
   const wbarray<T> &Om,
   const wbarray<T> &CC,
   unsigned s,
   char dblock,
   char mflag
){
   unsigned i,k, n=Om.SIZE.prod();
   double x,q, r=0, emin=0.1*idx2om(0);


   double bfac=(!mflag && Temp>0 ? 10/Temp : 0);

   double zfac=(Temp>0 ? 1/Temp : 0);

   if (!Om.hasSameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR size mismatch between Om and CC!");
   }
   if (dblock) {
      if (Om.SIZE.len!=2 || Om.SIZE[0]!=Om.SIZE[1]) wblog(FL,
         "ERR %s() quadratic block expected (%s)",
         FCT,Om.sizeStr().data);
   }

   if (s>=Ap.dim1 || s>=Ar.len) wblog(FL,
      "ERR index out of bounds (%d/%d;%d)",s, Ap.dim1,Ar.len);

   for (i=0; i<n; ++i) { x=Om.data[i]; k=om2idx(x);
       if (x<0. || (x<=emin && s%2))
            An(s,k)+=CC.data[i];
       else Ap(s,k)+=CC.data[i];

       if (x!=0) {
          q=CC.data[i]/x;
          if (bfac>0) { x=bfac*fabs(x); if (x<10) {
             double w=exp(-x*x);
             q*=(1-w);
          }}
          r+=q;
       }
       else if (mflag && !(s%2)) {"pos." omega
          r+=zfac*CC.data[i];
       }
   }

   Ar[s]+=r;

   if (mspec.dim2) {
      double dbl; unsigned j;

      if (mspec.dim1!=Ap.dim1) {
         wblog(FL,
           "WRN mspec must be initialized AFTER ASpec!\n"
           "resize (%d->%d x %d)", mspec.dim1, Ap.dim1, mspec.dim2
         );
         mspec.Resize(Ap.dim1,mspec.dim2);
      }

      for (i=0; i<n; i++)
      for (dbl=CC.data[i], j=0; j<mspec.dim2; j++) {
         if (j) dbl*=Om.data[i];
         mspec(s,j)+=dbl;
      }
   }
};


template <class T>
void Spectral::Add2buf(
   const wbarray<T> &Om,
   const wbarray<T> &CC,
   unsigned s
){
   unsigned i,k,n=Om.SIZE.prod();

   if (!Om.hasSameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR Size mismatch between Om and CC!");
   }

   if (s>Ap_buf.dim1)
   wblog(FL, "ERR Index out of bounds (%d,%d)",s, Ap_buf.dim1);

   for (i=0; i<n; i++) {
       k=om2idx(Om.data[i]);
       if (Om.data[i]<0.)
            An_buf(s,k)+=CC.data[i];
       else Ap_buf(s,k)+=CC.data[i];
   }
};


void Spectral::crossAddBuf(
   const char minor,
   char iflag
){
   unsigned s, Nom=Ap.dim2;
   int i;

   if (!Ap.hasSameSize(Ap_buf) || !An.hasSameSize(An_buf) ||
       !Ap.hasSameSize(An)
   ){
       Ap.info("Ap"); Ap_buf.info("Ap_buf");
       An.info("An"); An_buf.info("An_buf");
       wblog(FL,"ERR Spectral - severe size mismatch.");
   }




   for (s=0; s<Ap.dim1; s++) { i=0;

      i+=crossAddBuf_aux(Ap.rec(s), Ap_buf.rec(s), Nom,minor,iflag);
      i+=crossAddBuf_aux(An.rec(s), An_buf.rec(s), Nom,minor,iflag);

      if (i>1) {
         wblog(FL,"WRN empty buffer (%d/%d) !??",s,Ap.dim1);
         continue;
      }
   }



};


int Spectral::crossAddBuf_aux(
   double *a,
   const double *b,
   const unsigned Nom,
   const char minor,
   char iflag
){
   unsigned i, i1, i2, amin, amax, bmin, bmax;
   double t, alpha;

   if (minor<1 || minor>2)
   wblog(FL,"ERR Invalid minor flag %d", minor);

   for (i=0;   i<Nom;     ++i) if (a[i]!=0.) break; amin=i;
   for (i=Nom-1; i>=amin; --i) if (a[i]!=0.) break; amax=i;
   for (i=0;   i<Nom;     ++i) if (b[i]!=0.) break; bmin=i;
   for (i=Nom-1; i>=bmin; --i) if (b[i]!=0.) break; bmax=i;


   if (bmin>bmax) {
      return iflag;
   }

   if (amin>amax) {
      for (i=bmin; i<=bmax; i++) a[i]=b[i];
      return 0;
   }

   i1=MAX(amin,bmin); i2=MIN(amax,bmax);

   if (i1==i2) { a[i1]=0.5*(a[i1]+b[i2]); return 0; }


   alpha = 1. / double(i2-i1);

   for (i=0; i<Nom; ++i) {
      if (i<i1 || i>i2) a[i]+=b[i];
      else {
         t = alpha * double(i-i1);
         if (minor==1)
              a[i] = (1.-t)*a[i] + t*b[i];
         else a[i] = t*a[i] + (1.-t)*b[i];
      }
   }

   return 0;
};


void Spectral::saveIter(const char *F, int L, unsigned iter){

   unsigned n=Ap.numel();

   if (AP.SIZE.len!=3 || AP.SIZE!=AN.SIZE || n!=AP.SIZE[0]*AP.SIZE[1])
   wblog(F_L,
     "ERR %s() invalid initialization of AP/AN [%s; %s; %dx%d]",
      FCT,AP.sizeStr().data, AN.sizeStr().data, Ap.dim1, Ap.dim2
   ); 
   if (iter>=AP.SIZE[2]) wblog(F_L,
     "ERR %s() iter out of bounds AP/AN [%s; %d/%d]",
      FCT,AP.sizeStr().data, iter+1, AP.SIZE[2]
   ); 

   memcpy(AP.data+iter*n, Ap.data, n*sizeof(Ap.data[0]));
   memcpy(AN.data+iter*n, An.data, n*sizeof(An.data[0]));
};


void Spectral::getOmega(wbvector<double> &om, unsigned &linfit) const {

   unsigned i, Nom=Ap.dim2;
   double dbl;

   if (linfit) linfit=unsigned(fac);
   if (linfit>Nom) wblog(FL,"ERR Nom<linfit (%d,%d) !??",Nom,linfit);

   om.init(2*Nom);

   dbl = idx2om(linfit) / (linfit+0.5);

   for (i=0; i<linfit; ++i) om[Nom-1-i] = -(om[Nom+i] = dbl*(i+0.5));
   for (   ; i<Nom;    ++i) om[Nom-1-i] = -(om[Nom+i] = idx2om(i)  );
}


template <class TQ, class TD>
void Spectral::applyIROPfac(
   const wbvector<QSpace<TQ,TD> > &F, wbvector<double> *fac_
){
   unsigned i1,i2,s,m,r, n=F.len, gotfac=0; INDEX_T nq; double x;
   wbvector<double> fac(n); fac.set(1);

   if (2*n!=Ap.dim1 || 2*n!=An.dim1 || 2*n!=Ar.len) wblog(FL,
      "ERR %s() size mismatch (%d,%d,%d/2*%d)",
       FCT,Ap.dim1,An.dim1,Ar.len,n); 

   for (s=0; s<n; ++s) { r=F[s].rank(FL); i1=2*s; i2=2*s+1;
      if (r<3) { m=nq=1; }
      else {
         if (r>3) wblog(FL,"ERR %s() invalid rank-%g !?",FCT,r);
         m=F[s].getDim(2,&nq);
      }
      if (nq==1) continue; else
      if (int(nq)<1) wblog(FL,"ERR %s() m=%d, nq=%d !??",FCT,m,nq);

      fac[s]=nq; x=1/fac[s]; ++gotfac;
      wblog(FL," *  applying IROP factor 1/%d to op(%d)",nq,s+1);

      Ap.scaleRow(i1,x); An.scaleRow(i1,x); Ar[i1]*=x;
      Ap.scaleRow(i2,x); An.scaleRow(i2,x); Ar[i2]*=x;
   }

   if (fac_) {
      if (gotfac) (*fac_)=fac;
      else (*fac_).init();
   }
   if (!gotfac) return;

   if (!Ap_buf.isEmpty() || !An_buf.isEmpty()) {
      if (!Ap_buf.hasSameSize(Ap) || !An_buf.hasSameSize(Ap)) wblog(FL,
         "ERR %s() severe size mismatch %dx%d <> %dx%d <> %dx%d !??", FCT,
          Ap_buf.dim1,Ap_buf.dim2,An_buf.dim1,An_buf.dim2,Ap.dim1,Ap.dim2
      );
      for (s=0; s<n; ++s) { if (fac[s]>1) { x=1/fac[s]; i1=2*s; i2=2*s+1;
         Ap_buf.scaleRow(i1,x); An_buf.scaleRow(i1,x);
         Ap_buf.scaleRow(i2,x); An_buf.scaleRow(i2,x);
      }}
   }

   if (!AP.isEmpty() || !AN.isEmpty()) {
      if (AP.SIZE.len!=3 || !AP.hasSameSize(AN) || AP.SIZE[1]!=2*n)
      wblog(FL,
         "ERR %s() size mismatch (%s; %s; %dx%d) !??", FCT,
          AP.sizeStr().data, AN.sizeStr().data, Ap.dim1,Ap.dim2
      );

      double *ap=AP.data, *an=AN.data;
      unsigned n2=2*AP.SIZE[0]; m=n*AP.SIZE[2];

      for (i1=0; i1<m; ++i1, ap+=n2, an+=n2) { s=i1%n;
         if (fac[s]>1) { x=1/fac[s];
            timesRange(ap, x, n2);
            timesRange(an, x, n2);
         }
      }
   }

   if (!A0.isEmpty() || !mspec.isEmpty()) { wblog(FL,
        "WRN IROPfac yet not applied to A0 and mspec (%d: %d,%d)",
         gotfac, A0.isEmpty(), mspec.isEmpty());
      wblog(FL," *  A0: %dx%d",A0.dim1,A0.dim2);
      wblog(FL," *  sp: %dx%d",mspec.dim1,mspec.dim2);
   }
};


void Spectral::getRawData(
   wbvector<double> &om,
   wbMatrix<double> &A,
   unsigned linfit,
   const char *isbuf
) const {
   unsigned i,s,Nom=Ap.dim2;

   getOmega(om, linfit);

   A.init(Ap.dim1, 2*Nom);

   if (A.dim2!=om.len || !Ap.hasSameSize(An)) wblog(FL,
   "ERR Severe dimension mismatch (%d,%d).", A.dim2, om.len);

   if (isbuf && isbuf[0]) {
      if (strcmp(isbuf,"buf"))
      wblog(FL, "ERR Invalid flag `%s'", isbuf);

      if (!Ap.hasSameSize(Ap_buf)) {
         if (Ap_buf.isEmpty())
         wblog(FL, "ERR Buffer not set !?"); else
         wblog(FL, "ERR Severe size mismatch !?");
      }

      for (s=0; s<Ap.dim1; s++) { A(s,Nom)=NAN;
         for (i=0; i<Nom; ++i) {
            A(s,Nom-i-1) = An_buf(s,i);
            A(s,Nom+i  ) = Ap_buf(s,i);
         }
      }
   }
   else {
      for (s=0; s<Ap.dim1; s++) { A(s,Nom)=NAN;
         for (i=0; i<Nom; ++i) {
            A(s,Nom-i-1) = An(s,i);
            A(s,Nom+i  ) = Ap(s,i);
         }
      }
   }

   if (linfit) {
      for (s=0; s<A.dim1; s++) {
          set2avg(A.rec(s)+Nom,        linfit);
          set2avg(A.rec(s)+Nom-linfit, linfit);
      }
   }
};


void Spectral::pairRawData() {
   unsigned i,j,k, m=Ap.dim1, n=Ap.dim2;

   if (!Ap.hasSameSize(An)) wblog(FL, "ERR Severe size mismatch!");
   if (Ap.dim1%2) wblog(FL,
   "ERR Ap has odd number of rows !?? (%d)", Ap.dim1);

   for (i=0; i<m; i+=2)
   for (k=i/2, j=0; j<n; j++) {
      Ap(k,j)=Ap(i,j)+Ap(i+1,j);
      An(k,j)=An(i,j)+An(i+1,j);
   }

   Ap.Resize(m/2,n);
   An.Resize(m/2,n);

   if (mspec.dim1==0) return;


   m=mspec.dim1; n=mspec.dim2;

   if (m%2) wblog(FL,
   "ERR mspec has odd number of rows !?? (%d)", m);

   for (i=0; i<m; i+=2)
   for (k=i/2, j=0; j<n; j++) mspec(k,j)=mspec(i,j)+mspec(i+1,j);

   mspec.Resize(m/2,n);
}


void Spectral::getPARTIAL(
   wbarray<double> &A1, wbarray<double> &A2
 ) const {

   if (AP.isEmpty() && AN.isEmpty()) {
      A1.init(); A2.init(); return;
   }

   if (AP.SIZE.len!=3 || AP.SIZE!=AN.SIZE || AP.SIZE[1]%2) wblog(FL,
     "ERR %s() invalid initialization of AP/AN [%s; %s; %dx%d]",
      FCT,AP.sizeStr().data, AN.sizeStr().data, Ap.dim1, Ap.dim2
   );

   int i,j, m=AP.SIZE[1]/2, n=AP.SIZE[0], K=AP.SIZE[2], n1=n-1, n2=2*n*m;
   const double *ap=AP.data, *an=AN.data;
   double *a1, *a2;

   A1.init(2*n,m,K); a1=A1.data;
   A2.init(2*n,m,K); a2=A2.data;

   for (K*=m, j=0; j<K; ++j) {
      if (j>=m) {
         for (i=0; i<n; ++i) { a1[i]=an[n1-i] - an[n1-i-n2]; }; a1+=n; an+=n;
         for (i=0; i<n; ++i) { a1[i]=ap[   i] - ap[   i-n2]; }; a1+=n; ap+=n;
         for (i=0; i<n; ++i) { a2[i]=an[n1-i] - an[n1-i-n2]; }; a2+=n; an+=n;
         for (i=0; i<n; ++i) { a2[i]=ap[   i] - ap[   i-n2]; }; a2+=n; ap+=n;
      }
      else {
         for (i=0; i<n; ++i) { a1[i]=an[n1-i]; }; a1+=n; an+=n;
         for (i=0; i<n; ++i) { a1[i]=ap[   i]; }; a1+=n; ap+=n;
         for (i=0; i<n; ++i) { a2[i]=an[n1-i]; }; a2+=n; an+=n;
         for (i=0; i<n; ++i) { a2[i]=ap[   i]; }; a2+=n; ap+=n;
      }
   }
};


void Spectral::detailedBalance(
   const double T,
   wbvector<double> &om,
   wbMatrix<double> &aa
) const {

   unsigned i,j,k1,k2,m,n;
   wbvector<double> f1,f2;
   wbMatrix<double> ar;

   if (T<0) wblog(FL,"ERR T=%g !??", T);

   getRawData(om,ar);

   n=Ap.dim2; f1.init(2*n);

   if (T==0) {
      for (i=0; i<n; i++) f1[i]=1;
   }
   else {
      double beta=1./T;
      for (i=0; i<om.len; i++) f1[i]=1./(1.+exp(-beta*om[i]));
   }

   f1.flip(f2);

   m=Ap.dim1; n=ar.dim2;
   aa.init(2*m,n);

   for (i=0; i<m; i++) {
      k1=2*i; k2=k1+1;

      for (j=0; j<n; j++) {
         aa(k1,j)=ar(i,j)*f2[j];
         aa(k2,j)=ar(i,j)*f1[j];
      }
   }
}


void Spectral::dispSumRule() {

   wbvector<double> sp=Ap.recSum(),sn=An.recSum();
   unsigned i=0; double m=0;

   wblog(FL,"\r%60s%N"
       "   sum %sp_raw  : %s%N"
       "    *  %sn_raw  : %s ","",
       istr.data, sp.toStrf("%8.5f").data,
       istr.data, sn.toStrf("%8.5f").data
   );
   sp+=sn; 
   
   if (!A0.isEmpty()) {
      if (A0.dim1!=sp.len) wblog(FL,
      "ERR size mismatch (%d/%dx%d)",sp.len,A0.dim1,A0.dim2);
      sp+=A0.recSum();
   }

   if (sp.len) {
      m=round(sp[0]);
      for (i=0; i<sp.len; i++) if (ABS(sp[i]-m)>1E-2) break;
   }

   if (i==sp.len) { sp-=m; wblog(FL,"\r%60s\r    *  \e[34m"
      "total%-+3g: %s\e[0m ","",-m, sp.toStrf("%8.2g").data); }
   else {
      wblog(FL,"\r%60s\r    *  \e[34mtotal   : \\","");
      for (i=0; i<sp.len; ++i) { m=round(sp[i]);
         if (ABS(sp[i]-m)<1E-5) {
            if (m!=0)
                 printf("%2g%+5.0E ",m,sp[i]-m);
            else printf("%8.1E ",sp[i]-m);
         }
         else { printf("%8.5f ",sp[i]); }
      }
      printf("\e[0m\n");
   }
};

void Spectral::dispSumRule(
   const wbvector<double> &om,
   const wbMatrix<double> &A
){
   wbvector<double> s(A.dim1);
   unsigned i; double m=0;

   if (om.len!=A.dim2) wblog(FL,
   "ERR Dimension mismatch (%d,%d)", om.len, A.dim2);

   for (unsigned i=0; i<A.dim1; i++)
   s[i] = wbIntTrapez(om.data, A.rec(i), om.len);

   wblog(FL,"\r%60s\r   "
   "int AA(om)  : %s","",s.toStrf("%8.5f").data);

   if (!A0.isEmpty()) {
      if (A0.dim1!=s.len) wblog(FL,
      "ERR size mismatch (%d/%dx%d)",s.len,A0.dim1,A0.dim2);
      s+=A0.recSum();
   }


   if (s.len) {
      for (m=round(s[0]), i=1; i<s.len; i++)
      if (ABS(s[i]-m)>1E-2) return;
   }

   s-=m; wblog(FL,"\r%60s\r    *  \e[34m"
   "total%-+3g: %s\e[0m","",-m, s.toStrf("%8.2g").data);
};


void Spectral::getSmoothSpectral(
   wbvector<double> &om,
   wbMatrix<double> &A,
   double eps, double sigma, double sigma2, char iflag
){
   unsigned linfit=0;

   getOmega(om,linfit);
   getSmoothSpectral(A,om,eps,sigma,sigma2,iflag);
};


void Spectral::getSmoothSpectral(
    wbMatrix<double> &A,
    const wbvector<double> &omega,
    double eps,
    double sigma,
    double sigma2,
    char iflag
){
    unsigned i,s,Nom=Ap.dim2;
    wbvector<double> OM(Nom), I;
    double x, dbl, w1=0., w2=0., om, aom;

    if (fac<=0 || sigma<=0 || sigma2<=0 || eps<=0 || Nom<4)
    wblog(FL,"ERR %s()\n"
       "fac=%.3g, sigma=%.3g, sigma2=%.3g, eps=%.3g, Nom=%d !??",
        FCT,fac,sigma,sigma2, eps,Nom
    );

    for (i=0; i<OM.len; i++) OM[i]=idx2om(i);

    if (iflag) {
       double w1=omega.min(), w2=omega.max(), w0=omega.aMin();
       if (w1==-w2)
            sprintf(str,"%.3g .. %.3g",w0,w2);
       else sprintf(str,"%.3g .. %.3g .. %.3g",w1,w0,w2);

       wblog(FL,
       "<i> getSmoothSpec: sigma=%g [eps=%.3g @ %.3g]",
            sigma, eps, sigma2);
       wblog(FL,
       " *  |omega range| = %s (input)\n"
       " *  |emin...emax| = %.3g .. %.3g (2x%d bins, %g/dec)",
            str, idx2om(0), idx2om(Nom-1), Nom, fac/log10(M_E));
       dispSumRule();

    A.init(Ap.dim1, omega.len);

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
        w1 = foldLogGaussian(
        aom, sigma, om>0 ? Ap.rec(s) : An.rec(s), Nom);

        if (x!=1.) {
           double b=MAX(sigma2*eps,om*exp(-sigma*sigma/4));

           w2 = foldGaussian( om, b, OM.data, Ap.rec(s), Nom)
              + foldGaussian(-om, b, OM.data, An.rec(s), Nom);
        }

        A(s,i) = x*w1 + (1.-x)*w2;
    }

    if (iflag) dispSumRule(omega,A);
};


double Spectral::foldLogGaussian(
   const double om, const double sigma,
   const double *aa, const unsigned n
){
   double xa=0., k, dbl, efac;

   if (sigma<=0. || fac<=0.)
   wblog(FL, "ERR sigma=%g, fac=%g ???", sigma, fac);

   k=OM2IDX(om)-0.5;

   efac=1./(fac*sigma);

#if 0

   double s2=sigma/2.;

   for (unsigned i=0; i<n; i++) {
       dbl=efac*(k-double(i))-s2;
       xa += aa[i]*exp(-dbl*dbl);
   }
   xa *= exp(+s2*s2) / (sigma*sqrt(M_PI)*om);

#else

   double s4=sigma/4.;

   for (unsigned i=0; i<n; i++) {
       dbl=efac*(k-double(i))-s4;
       xa += aa[i]*exp(-dbl*dbl);
   }
   xa *= 1 / (sigma*sqrt(M_PI)*om);

#endif

   return xa;
};


double Spectral::foldGaussian(
   const double om, const double sigma,
   const double *oo, const double *aa, const unsigned n
){
   double xa=0., dbl, efac;

   if (sigma<=0.) wblog(FL, "ERR sigma=%g ???", sigma);

   efac = 1./sigma;

   for (unsigned i=0; i<n; i++) {
       dbl=efac*(om-oo[i]);
       xa += aa[i]*exp(-dbl*dbl);
   }

   xa *= 1. / (sigma*sqrt(M_PI));

   return xa;
};


void Spectral::KKreal(
   wbMatrix<wbcomplex> &GZ, wbvector<double> &OM, wbvector<double> &DOM
) const {

   unsigned i,j,s, d=Ap.dim1, Nom;
   wbMatrix<double> II, RR;
   double doi, doj;

   getRawData(OM, II); Nom=OM.len;
   getDiff(OM,DOM);


#if 0

   wblog(FL,"TST Get smooth data ...");

   getSmoothSpectral(OM, II, 1E-6, 0.3, 0); Nom=OM.len;
   getDiff(OM,DOM);

   for (s=0; s<d; s++)
   for (i=0; i<Nom; i++) II(s,i) *= DOM[i];

#endif

   RR.init(d,Nom);


   for (i=0;   i<Nom; ++i)
   for (j=i+1; j<Nom; ++j) {

       doi=DOM[i]/(OM[j]-OM[i]);
       doj=DOM[j]/(OM[i]-OM[j]);

       for (s=0; s<d; s++) {
           RR(s,i) += II(s,j)*doi;
           RR(s,j) += II(s,i)*doj;
       }
   }

#if 0
   wbMatrix<double> Ip2, In2;
   wblog(FL,"TST backtransform real->imag");

   RR*=(1/M_PI); RR.put("RR","base",'r');

   KKimag(Rp,Rn,OM,Ip2,In2);

   Ip2.put("Ip","base",'r');
   In2.put("In","base",'r');

   wblog(FL,"ERR");
#endif


   II *= M_PI;

   GZ.set(RR,II);
}


void KKimag(
   const wbMatrix<double> &RR,
   const wbvector<double> &OM,
   wbMatrix<double> &II
){
   unsigned i,j,s, d=RR.dim1, Nom=RR.dim2;
   double doi, doj;
   wbvector<double> DOM;


   getDiff(OM, DOM); II.init(d,Nom);

   if (OM.len!=Nom) wblog(FL,
      "ERR Severe dimension mismatch (%d,%d) !??", OM.len, Nom);

   for (i=0;   i<Nom; ++i)
   for (j=i+1; j<Nom; ++j) {
       doi=DOM[i]/(OM[i]-OM[j]); doj=DOM[j]/(OM[j]-OM[i]);
       for (s=0; s<d; s++) {
           II(s,i) += RR(s,j)*doi;
           II(s,j) += RR(s,i)*doj;
       }
   }

   II *= (1./M_PI);
}


#endif

