#ifndef __WB_SPECTRAL_HCC__
#define __WB_SPECTRAL_HCC__

/* ---------------------------------------------------------------- *
 * class Spectral
 * AWb © May 2006
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

    Spectral() : l0(0), fac(1.) {};

    void init (unsigned QDIM,
       double Emin=1E-8, double Emax=10., unsigned nl=512
    );

    void initBuf() {
       if (Ap.dim1!=An.dim1) wblog(FL,
       "ERR Severe dimension mismatch (%d,%d).", Ap.dim1, An.dim1);

       Ap_buf.init(Ap.dim1, Ap.dim2);
       An_buf.init(An.dim1, An.dim2);
    };

    template <class T>
    void Add(const wbArray<T> &Om, const wbArray<T> &CC, unsigned s);

    template <class T>
    void Add2buf(const wbArray<T> &Om, const wbArray<T> &CC, unsigned s);

    void crossAddBuf(const char minor=1);

    void getRawData(
    wbvector<double> &om, wbMatrix<double> &A, const char *isbuf="");

    void dispSumRule();
    void dispSumRule(const wbvector<double> &om, const wbMatrix<double> &A);

    void getSmoothSpectral(wbvector<double> &om, wbMatrix<double> &A,
    const double eps, const double sigma=0.6, const char iflag=1);

    void getSmoothSpectral(wbMatrix<double> &A, const wbvector<double> &om,
    const double eps, const double sigma=0.6, const char iflag=1);

    void KKreal(
       wbMatrix<wbcomplex> &GZ,
       wbvector<double> *pO=NULL,
       wbvector<double> *pDO=NULL
    );


    wbMatrix<double> Ap, An;
    wbMatrix<double> Ap_buf, An_buf;

 protected:
 private:

    double l0, fac;

    double OM2IDX(const double &om) const {
       return ((log(ABS(om))-l0)*fac);
    }

    unsigned om2idx(const double &om) const {
       if (om==0.) return 0;

       int k = (int)floor(OM2IDX(om));
       if (k<0) k=0; else if (k>=(int)Ap.dim2) k=Ap.dim2-1;

       return (unsigned)k;
    }

    double idx2om(const unsigned &k) const {
       return exp( ((double)k+0.5)/fac + l0 );
    }

    void crossAddBuf_aux(
    double *a, const double *b, const unsigned N, const char minor=1);

    double foldLogGaussian(
       const double om, const double sigma,
       const double *aa, const unsigned n
    );

    double foldGaussian(
       const double om, const double sigma,
       const double *oo, const double *aa, const unsigned n
    );
};


void Spectral::init (
   unsigned QDIM, double Emin, double Emax, unsigned nl
){
   unsigned N;
   double l2;

   if (Emin<=0. || Emax<=0. || Emin>=0.1*Emax || nl==0) wblog(FL,
   "ERR Invalid parameter set [%g %g %d]%s.", Emin, Emax, nl,
   (Emin>0 && Emax>0 && Emin<Emax) ? " (require Emin<Emax/10)":"");

   l0 = log(Emin);
   l2 = log(Emax);

   fac= log10(M_E)*double(nl);

   N = (unsigned)((l2-l0)*fac)+1;

   Ap.init(QDIM,N);
   An.init(QDIM,N);
};


template <class T>
void Spectral::Add(
   const wbArray<T> &Om,
   const wbArray<T> &CC,
   unsigned s
){
   unsigned i,k,n=Om.SIZE.prod();

   if (!Om.hasSameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR Size mismatch between Om and CC!");
   }

   if (s>Ap.dim1)
   wblog(FL, "ERR Index out of bounds (%d,%d)",s, Ap.dim1);

   for (i=0; i<n; i++) {
       k=om2idx(Om.DATA[i]);
       if (Om.DATA[i]<0.)
            An(s,k)+=CC.DATA[i];
       else Ap(s,k)+=CC.DATA[i];
   }
};


template <class T>
void Spectral::Add2buf(
   const wbArray<T> &Om,
   const wbArray<T> &CC,
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
       k=om2idx(Om.DATA[i]);
       if (Om.DATA[i]<0.)
            An_buf(s,k)+=CC.DATA[i];
       else Ap_buf(s,k)+=CC.DATA[i];
   }
};


void Spectral::crossAddBuf(
   const char minor
){
   unsigned s, N=Ap.dim2;

   if (!Ap.hasSameSize(Ap_buf) || !An.hasSameSize(An_buf) ||
       !Ap.hasSameSize(An)
   ){
       Ap.info("Ap"); Ap_buf.info("Ap_buf");
       An.info("An"); An_buf.info("An_buf");
       wblog(FL,"ERR Spectral - severe size mismatch.");
   }

   for (s=0; s<Ap.dim1; s++) {
      crossAddBuf_aux(Ap.rec(s), Ap_buf.rec(s), N, minor);
      crossAddBuf_aux(An.rec(s), An_buf.rec(s), N, minor);
   }

};


void Spectral::crossAddBuf_aux(
   double *a,
   const double *b,
   const unsigned N,
   const char minor
){
   unsigned i, E1, E2, amin, amax, bmin, bmax;
   double t, alpha;

   if (minor<1 || minor>2)
   wblog(FL,"ERR Invalid minor flag %d", minor);

   for (i=0;   i<N;     i++) if (a[i]!=0.) break; amin=i;
   for (i=N-1; i>=amin; i--) if (a[i]!=0.) break; amax=i;
   for (i=0;   i<N;     i++) if (b[i]!=0.) break; bmin=i;
   for (i=N-1; i>=bmin; i--) if (b[i]!=0.) break; bmax=i;


   if (bmin>bmax) {
      wblog(FL, "ERR buffer contains no data yet (%d,%d; %d) !?",
      bmin, bmax, N);
   }

   if (amin>amax) {
      for (i=bmin; i<=bmax; i++) a[i]=b[i];
      return;
   }

   E1=MAX(amin,bmin); E2=MIN(amax,bmax);
   if (E1==E2) { a[E1]=0.5*(a[E1]+b[E2]); return; }


   alpha = 1. / (exp(double(E2-E1)/fac)-1.);

   for (i=0; i<N; i++) {
      if (i<E1 || i>E2) a[i]+=b[i];
      else {
         t = alpha * (exp(double(i-E1)/fac)-1.);
         if (minor==1)
              a[i] = (1.-t)*a[i] + t*b[i];
         else a[i] = t*a[i] + (1.-t)*b[i];
      }
   }
};


void Spectral::getRawData(
   wbvector<double> &om,
   wbMatrix<double> &A,
   const char *isbuf
){
   unsigned i,s,N=Ap.dim2;
   double omi;

   om.init(2*N+1);
   A.init(Ap.dim1,2*N+1);

   for (i=0; i<N; i++) {
      omi = idx2om(i);
      om[N-1-i]=-omi; om[N+1+i]=omi;
   }

   if (isbuf && isbuf[0]) {
      if (strcmp(isbuf,"buf"))
      wblog(FL, "ERR Invalid flag `%s'", isbuf);

      if (!Ap.hasSameSize(Ap_buf)) {
         if (Ap_buf.isEmpty())
         wblog(FL, "ERR Buffer not set !?"); else
         wblog(FL, "ERR Severe size mismatch !?");
      }

      for (s=0; s<Ap.dim1; s++) {
         A(s,N)=NAN;
         for (i=0; i<N; i++) {
            A(s,N-1-i) = An_buf(s,i);
            A(s,N+1+i) = Ap_buf(s,i);
         }
      }
   }
   else {
      for (s=0; s<Ap.dim1; s++) {
         A(s,N)=NAN;
         for (i=0; i<N; i++) {
            A(s,N-1-i) = An(s,i);
            A(s,N+1+i) = Ap(s,i);
         }
      }
   }
};


void Spectral::dispSumRule(){
   wbvector<double> sp=Ap.recSum(),sn=An.recSum();

   wblog(FL,
       "Sum Ap_raw      = %s\n"
       " ·  An_raw      = %s",
   sp.toStrf("%8.6f").data, sn.toStrf("%8.6f").data);

   sp+=sn; sp-=1;
   wblog(FL,
       " ·  total-1.    = %s", sp.toStrf("%8.2g").data);
}

void Spectral::dispSumRule(
   const wbvector<double> &om,
   const wbMatrix<double> &A
){
   wbvector<double> s(A.dim1);

   if (om.len!=A.dim2) wblog(FL,
   "ERR Dimension mismatch (%d,%d)", om.len, A.dim2);

   for (unsigned i=0; i<A.dim1; i++)
   s[i] = wbIntTrapez(om.data, A.rec(i), om.len);

   wblog(FL,
       "Sum AA_smooth   = %s\n"
       " ·  total-1.    = %s",
   s.toStrf("%8.6f").data, (s-1).toStrf("%8.2g").data);
}


void Spectral::getSmoothSpectral(
   wbvector<double> &om,
   wbMatrix<double> &A,
   const double eps,
   const double sigma,
   const char iflag
){
   unsigned i,N=Ap.dim2;
   double omi;

   om.init(2*N+1);

   for (i=0; i<N; i++) {
      omi = idx2om(i);
      om[N-1-i]=-omi; om[N+1+i]=omi;
   }

   getSmoothSpectral(A,om,eps,sigma,iflag);
};


void Spectral::getSmoothSpectral(
    wbMatrix<double> &A,
    const wbvector<double> &omega,
    const double eps,
    const double sigma,
    const char iflag
){
    unsigned i,s,N=Ap.dim2;
    wbvector<double> OM(N), I;
    double x, dbl, w1=0., w2=0., om, aom;

    if (fac<=0. || sigma<=0. || eps<=0. || N<4)
    wblog(FL,"ERR fac=%g, sigma=%g, eps=%g, N=%d ???",fac,sigma,eps,N);

    for (i=0; i<OM.len; i++) OM[i]=idx2om(i);

    A.init(Ap.dim1, omega.len);

    if (iflag) {
       wblog(FL,
           "<i> getSmoothSpectral ( 2x %d bins, %g/dec)\n"
           " ·  omega range = %.3g .. %.3g .. %.3g (%d)\n"
           " ·  emin/max    = %.3g .. %.3g\n"
           " ·  eps         = %g\n"
           " ·  sigma       = %g",
       N, fac/log10(M_E), omega.min(), omega.aMin(), omega.max(), omega.len,
       idx2om(0), idx2om(N), eps, sigma);

       printf("\n");
       dispSumRule();
       printf("\n");
    }
    else {
       wblog(FL,
         "<i> getSmoothSpectral (2x%d bins)\n"
         "    om: %.3g .. %.3g .. %.3g (%d)",
       N, omega.min(), omega.aMin(), omega.max(), omega.len);
    }

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
        aom, sigma, om>0 ? Ap.rec(s) : An.rec(s), N);

        if (x!=1.)
        w2 = foldGaussian( om, 0.8*eps, OM.data, Ap.rec(s), N)
           + foldGaussian(-om, 0.8*eps, OM.data, An.rec(s), N);

        A(s,i) = x*w1 + (1.-x)*w2;
    }

    dispSumRule(omega,A);
    printf("\n");
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
   wbMatrix<wbcomplex> &GZ, wbvector<double> *pO, wbvector<double> *pDO
){
   unsigned i,j,s, d=Ap.dim1, N=Ap.dim2;
   double doi, doj, oi, oj;
   wbMatrix<double> Rp, Rn;
   wbvector<double> OM, DOM;

   initBuf(); OM.init(N); DOM.init(N);

   for (i=0; i<OM.len; i++) {
       OM[i]=idx2om(i); if (i==0) { DOM[0]=2*OM[0]; continue; }
       DOM[i-1] += (DOM[i] = OM[i]-OM[i-1]);
       DOM[i-1] *= 0.5;
   }


   Rp.init(d,N); Rn.init(d,N);

   const wbMatrix<double> &Ip=Ap, &In=An;

#if 0
   wbMatrix<double> Ip, In;
   if (0) {
      wblog(FL,"TST");
      wbvector<double> om;
      wbMatrix<double> AA;
      getSmoothSpectral(om, AA, 1E-6, 0.6, 1);

      if (AA.dim2!=2*N+1) wblog(FL,"ERR ???");

      AA.getCols(0, N-1, Ip);
      AA.getCols(N+1,2*N,In);

      for (s=0; s<d; s++)
      for (i=0; i<N/2; i++) SWAP(Ip(s,i), Ip(s,N-i-1));

      for (s=0; s<d; s++)
      for (i=0; i<N; i++) { Ip(s,i)*=DOM[i]; In(s,i)*=DOM[i]; }
   }
   else {
      wblog(FL,"TST");
      Ip=Ap; In=An;
      for (i=0; i<N; i++) 
      In(0,i) = Ip(0,i) = ABS(OM[i])<=0.1 ? DOM[i] : 0;
   }
#endif

OM.put("OM"); DOM.put("DOM");
Ip.put("Ip0","base",'r');
In.put("In0","base",'r');

double dbl=0.; i=20;
for (j=0; j<N; j++) {
   dbl+=In(0,j)/(-OM[j]-OM[i]); if (i!=j)
   dbl+=Ip(0,j)/( OM[j]-OM[i]);
}  dbl/=M_PI;
mxPut(dbl,"dbl");

wblog(FL,"TST %g %g; %d: %.8g", In(0,0), Ip(0,0), i+1, dbl);

   for (i=0; i<N; i++) {
       oi=OM[i]; doi=DOM[i]/(oi+oi);
       for (s=0; s<d; s++) {
           Rp(s,i) -= In(s,i)*doi;
           Rn(s,i) += Ip(s,i)*doi;
       }

       for (j=i+1; j<N; j++) {
           oj=OM[j]; doi=DOM[i]/(oi+oj); doj=DOM[j]/(oi+oj);
           for (s=0; s<d; s++) {
               Rp(s,i) -= In(s,j)*doi;
               Rn(s,i) += Ip(s,j)*doi;

               Rp(s,j) -= In(s,i)*doj;
               Rn(s,j) += Ip(s,i)*doj;
           }

           doi=DOM[i]/(oj-oi); doj=DOM[j]/(oi-oj);
           for (s=0; s<d; s++) {
               Rp(s,i) += Ip(s,j)*doi;
               Rn(s,i) -= In(s,j)*doi;

               Rp(s,j) += Ip(s,i)*doj;
               Rn(s,j) -= In(s,i)*doj;
           }
       }
   }


wbMatrix<double> Ip2, In2;
wblog(FL,"TST backtransform real->imag");

Rp*=(1/M_PI); Rp.put("Rp","base",'r');
Rn*=(1/M_PI); Rn.put("Rn","base",'r');

KKimag(Rp,Rn,OM,Ip2,In2);

Ip2.put("Ip","base",'r');
In2.put("In","base",'r');

wblog(FL,"ERR");




   if (pO) {
      (*pO).init(2*N+1);
      for (i=0; i<N; i++) { ( *pO)[N+i+1]=OM[i]; (*pO)[N-i-1]=-OM[i]; }
   }
   if (pDO) {
      (*pDO).init(2*N+1);
      for (i=0; i<N; i++) { (*pDO)[N+i+1] = (*pDO)[N-i-1] = DOM[i]; }
   }

   GZ.init(d,2*N+1);
   for (s=0; s<d; s++) {
      GZ(s,N)=NAN;

      for (i=0; i<N; i++) {
         GZ(s,N-1-i).set(Rn(s,i), M_PI*In(s,i));
         GZ(s,N+1+i).set(Rp(s,i), M_PI*Ip(s,i));
      }
   }
}


void KKimag(
   const wbMatrix<double> &Rp,
   const wbMatrix<double> &Rn,
   const wbvector<double> &OM,
   wbMatrix<double> &Ip,
   wbMatrix<double> &In
){
   unsigned i,j,s, d=Rp.dim1, N=Rp.dim2;
   double doi, doj, oi, oj;
   wbvector<double> DOM;


   Ip.init(d,N); In.init(d,N); DOM.init(N);

   for (i=0; i<OM.len; i++) {
       if (i==0) { DOM[0]=2*OM[0]; continue; }
       DOM[i-1] += (DOM[i] = OM[i]-OM[i-1]);
       DOM[i-1] *= 0.5;
   }

   for (i=0; i<N; i++) {
       oi=OM[i]; doi=DOM[i]/(oi+oi);
       for (s=0; s<d; s++) {
           Ip(s,i) += Rn(s,i)*doi;
           In(s,i) -= Rp(s,i)*doi;
       }

       for (j=i+1; j<N; j++) {
           oj=OM[j]; doi=DOM[i]/(oi+oj); doj=DOM[j]/(oi+oj);
           for (s=0; s<d; s++) {
               Ip(s,i) += Rn(s,j)*doi;
               In(s,i) -= Rp(s,j)*doi;

               Ip(s,j) += Rn(s,i)*doj;
               In(s,j) -= Rp(s,i)*doj;
           }

           doi=DOM[i]/(oj-oi); doj=DOM[j]/(oi-oj);
           for (s=0; s<d; s++) {
               Ip(s,i) -= Rp(s,j)*doi;
               In(s,i) += Rn(s,j)*doi;

               Ip(s,j) -= Rp(s,i)*doj;
               In(s,j) += Rn(s,i)*doj;
           }
       }
   }

   Ip *= (1./M_PI);
   In *= (1./M_PI);
}


#endif

