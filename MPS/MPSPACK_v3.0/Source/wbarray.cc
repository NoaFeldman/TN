#ifndef __WB_ARRAY_COL_MAJOR_CC__
#define __WB_ARRAY_COL_MAJOR_CC__

//-------------------------------------------------------------------//
// wbarray - N-D Array class
//
//    generic N-dimensional class that allows permutations
//    and reshaping; index order is column-major.
//
// Wb,Sep14,09 ;  Wb,Aug11,05
//-------------------------------------------------------------------//

template<class T>
void wbarray<T>::setRand(const char pnflag) {
    static int firstcall=1; size_t i,s=numel(); double fac;

    if (WbUtil<T>::isFloat())
         fac = 1.  /(double)RAND_MAX;
    else fac = 100./(double)RAND_MAX;

    if (firstcall) {
        ::srand((unsigned)time((time_t*)NULL)); firstcall=0;
    }
    for (i=0; i<s; ++i) data[i] = (T)(fac*::rand());

    if (pnflag) {
       fac *= (0.5*(double)RAND_MAX);
       for (i=0; i<s; ++i) data[i] -= fac;
    }
}


template<> inline
double wbarray<wbcomplex>::SkipTiny(const double eps) {
   double a,x2=0;
   for (size_t s=numel(), i=0; i<s; ++i) {
      a=fabs(data[i].r); if (a>0 && a<eps) { x2+=a*a; data[i].r=0; }
      a=fabs(data[i].i); if (a>0 && a<eps) { x2+=a*a; data[i].i=0; }
   }
   return std::sqrt(x2);
}

template<> inline
wbarray<double>& wbarray<double>::SkipTiny_float(double ref) {
   Wb::chopTiny_float(data,numel(),ref);
   return *this;
};

template<> inline
wbarray<wbcomplex>& wbarray<wbcomplex>::SkipTiny_float(wbcomplex ref) {
   Wb::chopTiny_float((double*)data,2*numel(),ref.r);
   return *this;
};


template<class T> inline
void wbarray<T>::operator+=(const wbarray<T> &B) { Plus(B); }

template<class T> inline
void wbarray<T>::operator-=(const wbarray<T> &B) { Plus(B,-1); }

template<class T> inline wbarray<T>&
wbarray<T>::plus(
   const wbarray<T> &B, wbarray<T> &C, T bfac, char iflag) const {
   C=*this; C.Plus(B,bfac,iflag); return C;
}

template<class T> inline wbarray<T>&
wbarray<T>::minus(
   const wbarray<T> &B, wbarray<T> &C, T bfac, char iflag) const {
   C=*this; C.Plus(B,-bfac,iflag); return C;
}

template<class T> inline
wbarray<T>& wbarray<T>::Minus(
    const wbarray<T> &B, T bfac, char iflag) {
    Plus(B,-bfac,iflag); return *this;
}

template<class T>
wbarray<T>& wbarray<T>::Plus(const wbarray<T> &B, T bfac, char iflag) {

    if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
    size_t i, s=numel();

    if (isEmpty() && iflag) { *this=B; (*this)*=bfac; return *this; }

    if (SIZE!=B.SIZE) wblog(FL,
       "ERR %s() severe size mismatch\n--> %s + %s ???",FCT,
        sizeStr().data, B.sizeStr().data);

    if (bfac==T(+1)) for (i=0; i<s; ++i) data[i] += B.data[i]; else
    if (bfac==T(-1)) for (i=0; i<s; ++i) data[i] -= B.data[i];
    else             for (i=0; i<s; ++i) data[i] += bfac*B.data[i];

    return *this;
};


template<class T> inline
char wbarray<T>::sameUptoFac(
  const wbarray<T> &B, T* fac_, double eps
) const {

   size_t i, n=numel(); T a,fac=1;

   if (SIZE!=B.SIZE) return 1;
   if (!n) return 0;

   a=this->aMax(&i);
   if (a<=eps) {
      if (a==0 || B.aMax(&i)>eps) { if (fac_) (*fac_)=0; return 0; }
      else return 1;
   }
   if (ABS(B.data[i])<=eps) return 2;

   fac=data[i]/B.data[i]; if (fac_) (*fac_)=fac;

   for (i=0; i<n; ++i) {
      a=ABS(data[i]-fac*B.data[i]);
      if (a>eps) return 3;
   }
   return 0;
};


template<class T> inline
bool wbarray<T>::hasSameSize (const wbarray<T> &B) const {

   if (SIZE.len==B.SIZE.len) { return (SIZE==B.SIZE); }
   else {
      size_t i,r, r1=SIZE.len, r2=B.SIZE.len, i1=0, i2=0;

      for (; i1<r1; ++i1) if (  SIZE[i1]!=1) break;
      for (; i2<r2; ++i2) if (B.SIZE[i2]!=1) break;

      r=MIN(r1-i1, r2-i2);
      for (i=0; i<r; ++i) { if (SIZE[i1+i]!=B.SIZE[i2+i]) return 0; }

      for (i=i1+r; i<r1; ++i) { if (  SIZE[i]!=1) return 0; }
      for (i=i2+r; i<r2; ++i) { if (B.SIZE[i]!=1) return 0; }

      return 1;
   }
};

template<class T> inline
bool wbarray<T>::hasSameSize(const size_t *s, size_t n) const {

   size_t i=0, m=MIN(n,SIZE.len);
   for (i=0; i<m; ++i) { if (s[i]!=SIZE.data[i]) return 0; }
   for (   ; i<n; ++i) { if (s[i]!=1) return 0; }
   return 1;
};


template<class T> inline
bool wbarray<T>::equal (const wbarray<T> &B, const T& eps) const {
    size_t i,s=numel(); if (SIZE!=B.SIZE) return 0;

    if (eps==0) { if (memcmp(data,B.data,s*sizeof(T))) return 0; }
    else { for (i=0; i<s; i++) if (fabs(data[i]-B.data[i])>eps) return 0; }

    return 1;
}

template<class T> inline
bool wbarray<T>::unequal (const wbarray<T> &B, const T& eps) const {
    return !equal(B,eps);
}

template<class T> inline
bool wbarray<T>::equal (
    const wbarray<T> &B, const double& eps, double &maxdiff) const
{
    if (SIZE!=B.SIZE) return 0;

    size_t i,s=numel(); double d;

    for (maxdiff=0, i=0; i<s; i++) {
        d = ABS( data[i] - B.data[i] );
        if (maxdiff<d) maxdiff=d;
    }

    return (maxdiff<eps);
}

template<class T> inline
bool wbarray<T>::unequal (
    const wbarray<T> &B, const double& eps, double &maxdiff) const {
    return !(*this).equal(B,eps,maxdiff);
}


template<class T>
bool wbarray<T>::isDiag_aux(const T eps, const char* task) const {
 
   size_t i=0;
   wbIndex I(SIZE);

   char isDiag=!strcmp(task,"isDiag");
   char isIdty=!strcmp(task,"isIdty");

   if (!isDiag && !isIdty) wblog(FL,"ERR %s() invalid task >%s<",FCT,task);
   if (SIZE.len%2) wblog(FL,"ERR %s() for even-rank objects only (%s)",
     task,sizeStr().data);

   while (++I) {
       if (!I.isdiag()) {
          if (ABS(data[i])>eps) return 0;
       }
       else {
          if (isIdty && (ABS(data[i]-T(1)))>eps) return 0;
       }
       ++i;
   }

   return 1;
};


template<class T>
bool wbarray<T>::isDiag_aux(double *epsp, const char* task) const {
 
   size_t i,k, s=numel(), r=SIZE.len, l=r-1;
   WBINDEX I(r); bool is=1; double a, eps=0, eref=*epsp;

   char isDiag=!strcmp(task,"isDiag");
   char isIdty=!strcmp(task,"isIdty");
   if (!isDiag && !isIdty) wblog(FL,"ERR Invalid task >%s<", task);

   if (!isRank(2)) { info("this"); wblog(FL,
   "ERR %s only appies to rank-2 objects (%s).",task,sizeStr().data); }

   for (i=0; i<s; i++) {
       if (I[0]!=I[1]) { a=ABS(data[i]);
          if (is && eref<a) is=0;
          if (eps<a) eps=a;
       }
       else if (isIdty) { a=ABS(data[i]-T(1));
          if (is && eref<a) is=0;
          if (eps<a) eps=a;
       }

       k=0; I[0]++;
       while(I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   *epsp=eps; return is;
};


template<class T>
bool wbarray<T>::isProptoId(T &x, const T eps) const {
 
   size_t i=0;
   wbIndex I(SIZE);

   if (SIZE.len%2) wblog(FL,
      "ERR %s() for even-rank objects only (%s)",FCT,sizeStr().data);
   if (x==T(0) && data) x=data[0];

   while (++I) {
      if (!I.isdiag()) 
           { if (ABS(data[i]  )>eps) return 0; }
      else { if (ABS(data[i]-x)>eps) return 0; }; ++i;
   }

   if (i && !data) wblog(FL,
      "ERR %s() invalid empty array (len=%d)",FCT,i);

   return 1;
};


template<class T>
bool wbarray<T>::isOpS(size_t *n) const {

   size_t i,r=SIZE.len, r2=r/2;

   if (r==0) return 1;
   if (r%2) return 0;

   for (i=0; i<r2; i++) if (SIZE[i]!=SIZE[i+r2]) return 0;

   if (n) { (*n)=prodRange(SIZE.data,r2); }

   return 1;
}


template<class T>
wbarray<T>& wbarray<T>::balanceOp(
   const char *F, int L, T &dref, double &dscale
){
   size_t i,j,l,n=0;

   if (!isOpS(&n)) wblog(F_L,
   "ERR %s() requires operator (%s)",FCT,sizeStr().data);

   if (n<=1) { dref=0; dscale=1; return *this; }
   dref=0; dscale=0;

   for (l=i=0; i<n; i++, l+=(n+1)) dref+=data[l];
   dref /= n;

   for (l=i=0; i<n; i++)
   for (j=0; j<n; j++, l++) { if (i==j) data[l]-=dref;
       dscale += ABS2(data[l]);
   }
   dscale = std::sqrt(dscale)/n;
   
   if (dscale!=0) (*this)*=(1./dscale);
   else dscale=1;

   return *this;
}


template<class T>
wbarray<T>& wbarray<T>::Symmetrize(
   const char *F, int L,
   double *delta,
   char cflag,
   char tflag,
   char fflag
){
   size_t i,j,n=0; T d1,d2; double e,E=0,dmax=0;

   if (!isOpS(&n)) wblog(F,L,
      "ERR %s() (generalized) square matrix required (%s)",
       FCT,sizeStr().data);
   if (tflag && !delta) return *this;

   if (cflag && typeid(T)!=typeid(wbcomplex)) cflag=0;

   for (i=0; i<n; i++)
   for (j=i; j<n; j++) { d1=data[i+j*n]; d2=data[j+i*n];
      dmax=MAX(MAX(dmax,ABS(d1)),ABS(d2));
      if (i!=j) {
         e=(cflag ? ABS(d1-CONJ(d2)) : ABS(d1-d2));
         E=MAX(E,e);
      }
   }

   dmax=MAX(1.,dmax); 
   e=E; if (dmax!=0) e/=dmax;

   if (delta) {
      if (dmax!=0) {
         if (e>(*delta)) { tflag=1;
#ifdef MATLAB_MEX_FILE
         this->put(FL,"M_");
#endif
            sprintf(str,"%s() non-symmetric operator (%d)\n"
              "eps=%.3g/%.3g @ dmax=%.3g",FCT,cflag,e,*delta,dmax);

            if (fflag) wblog(F,L,"ERR %s",str);
            else wblog(F,L,"WRN %s",str);
         }
         *delta=(e);
      }
      else if (e==0) *delta=0; else *delta=1;
   }

   if (tflag==0) {
      if (cflag) {
         for (i=0; i<n; i++)
         for (j=i+1; j<n; j++) {
            d1=0.5*(data[i+j*n]+CONJ(data[j+i*n]));
            data[i+j*n]=d1; data[j+i*n] = CONJ(d1);
         }
      }
      else {
         for (i=0; i<n; i++)
         for (j=i+1; j<n; j++) {
            d1=0.5*(data[i+j*n]+data[j+i*n]);
            data[i+j*n] = data[j+i*n] = d1;
         }
      }
   }

   return *this;
};

   



         

template<class T> 
int wbarray<T>::Householder(
   const char *F, int L, size_t k,
   T *u, WBPERM *P,
   T eps
){
   if (rank()!=2) wblog(F_L,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   size_t i=0, j=0, k2=k, dim1=SIZE[0], dim2=SIZE[1];
   if (k>=MIN(dim1,dim2)) wblog(F_L,"ERR %s() "
      "index out of bounds (%d; %s)",FCT,k+1,sizeStr().data);

   const T zero=0;
   T ua, nrm, x2=zero, x0max=zero, *x;

   for (; i<k; ++i) u[i]=zero;

   if (P) {
      if (k==0) P->init(dim2); else
      if (P->len!=dim2) wblog(FL,"ERR %s() "
         "got invalid P for (k=%d: %d/%d) !??",FCT,k+1,P->len,dim2);
      T x2max=zero;

      for (x=data, j=0; j<dim2; ++j, x+=dim1) {
         if (!(P->data[j])) {
            for (x2=zero, i=k; i<dim1; ++i) x2 += (CONJ(x[i])*x[i]);
            if (x2max<x2) { x2max=x2; k2=j; }
         }
         else {
            if (P->data[j]>k) wblog(FL, "ERR %s() "
               "got invalid P (%d: %d/%d) !?",FCT,j+1,P->data[j],k+1);
            if (x0max<x[P->data[j]-1]) x0max=x[P->data[j]-1];
         }
      }; x2=x2max;

      x=data+k2*dim1;

      for (i=k; i<dim1; ++i) { u[i]=x[i]; };
   }
   else { x=data;
      for (i=0; i<k; ++i, x+=dim1) { if (x0max<x[i]) x0max=x[i]; }
      for (; i<dim1; ++i, x+=dim1) { u[i]=x[i]; x2+=(CONJ(x[i])*x[i]); }
   }

   eps*=(x0max<x2 ? x2 : x0max);

   nrm=SQRT(x2); if (nrm<eps) { return -1; }

   u[k] += (x[k]>0 ? +nrm : -nrm);

   nrm=T(1)/SQRT(x2 + CONJ(u[k])*u[k] - CONJ(x[k])*x[k]);
   for (i=k; i<dim1; ++i) { u[i]*=nrm; }


   if (P) {
      for (x=data, j=0; j<dim2; ++j, x+=dim1) {
         if (!(P->data[j])) {
            householder(u+k,x+k,dim1-k,eps);
         }
      }
      P->data[k2]=(k+1);"+1")
   }
   else {
      for (x=data+k*dim1, j=k; j<dim2; ++j, x+=dim1) {
         householder(u+k,x+k,dim1-k,eps);
      }
   }

   return k2;
};


template<class T> inline
void householder(const T *u, T *x, size_t n, const T& eps) {

    size_t i=0; T zero=0, ua=zero;
    for (; i<n; ++i) ua += (CONJ(u[i])*x[i]); 

    if (ua!=zero) { ua*=2;
       for (i=0; i<n; ++i) x[i] -= (u[i]*ua);
       if (eps>zero) { ua=zero;
          for (i=0; i<n; ++i) ua += CONJ(x[i])*x[i];
          if (ua<eps)
          for (i=0; i<n; ++i) x[i]=zero;
       }
    }
};


template<class T>
wbarray<T>& wbarray<T>::swapRows(size_t i1, size_t i2) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (i1>=SIZE[0] || i2>=SIZE[0]) wblog(FL,"ERR %s() "
      "index out of bounds (%d,%d; %s)",FCT,i1+1,i2+1,sizeStr().data);

   if (i1!=i2) { size_t dim1=SIZE[0], dim2=SIZE[1]; if (dim2) {
      T x, *d1=data+i1, *d2=data+i2;
      for (size_t j=0; j<dim2; ++j, d1+=dim1, d2+=dim1) {
         x=(*d1); *d1=(*d2); *d2=x;
      }
   }}
   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::swapCols(size_t j1, size_t j2) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (j1>=SIZE[1] || j2>=SIZE[1]) wblog(FL,"ERR %s() "
      "index out of bounds (%d,%d; %s)",FCT,j1+1,j2+1,sizeStr().data);

   if (j1!=j2) { size_t dim1=SIZE[0]; if (dim1) {
      T x, *d1=data+dim1*j1, *d2=data+dim1*j2;
      for (size_t i=0; i<dim1; ++i) {
         x=d1[i]; d1[i]=d2[i]; d2[i]=x;
      }
   }}
   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::setCol(
   size_t k, const T* d0, const char *F, int L
){
   if (rank()!=2) wblog(F_L,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (k>=SIZE[1]) wblog(F_L,
      "ERR %s() index out of bounds (%d; %s)",FCT,k+1,sizeStr().data);

   size_t dim1=SIZE[0];
   if (dim1) {
      T* d=ref(size_t(0),k);
      if (d!=d0) MEM_CPY<T>(d,dim1,d0);
   }
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setCol(
   size_t k, size_t k0, const char *F, int L
){
   if (rank()!=2) wblog(F_L,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (k>=SIZE[1] || k0>=SIZE[1]) wblog(F_L,
      "ERR %s() index out of bounds (%d,%d; %s)",
       FCT,k+1,k0+1,sizeStr().data);

   size_t dim1=SIZE[0];
   if (k!=k0 && dim1) {
      MEM_CPY<T>(ref(0U,k),dim1,ref(0U,k0));
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setRow(
   size_t k, const T* d0, const char *F, int L
){
   if (rank()!=2) wblog(F_L,
   "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (k>=SIZE[0]) wblog(F_L,
   "ERR %s() index out of bounds (%d; %s)",FCT,k+1,sizeStr().data);

   size_t i=0, dim1=SIZE[0], dim2=SIZE[1];
   T* d=data+k;
   for (; i<dim2; i++, d+=dim1) d[i]=d0[i];

   return *this;
};


template<class T>
wbvector<T>& wbarray<T>::colNorm2(wbvector<T> &a) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);

   size_t i,j, dim1=SIZE[0], dim2=SIZE[1]; 
   const T* d0=data; T x2;
   
   a.init(dim2);
   for (j=0; j<dim2; j++, d0+=dim1) {
      for (x2=T(), i=0; i<dim1; i++) x2+=(CONJ(d0[i])*d0[i]);
      a[j]=x2;
   }
   return a;
};


template<class T>
T wbarray<T>::colNorm2(size_t k) const {

   size_t j=0, dim1=SIZE[0], dim2=SIZE[1]; 
   const T* d=data+k*dim1;
   T a=0;
   
   if (rank()!=2) wblog(FL,
      "ERR %s() for matrizes only (%s)",FCT,sizeStr().data);
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d; %s)",FCT,k,dim2,sizeStr().data);

   for (; j<dim2; ++j) { a+=CONJ(d[j])*d[j]; }
   return a;
};


template<class T>
T wbarray<T>::NormalizeCol(size_t k, char tnorm, char qflag) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2); 

   T *d=data+dim1*k, x2=overlap(d,d,dim1,1,tnorm);
   x2=SQRT(x2);
   
   if (x2!=0) timesRange(d,1/x2,dim1);
   else if (!qflag) wblog(FL,"ERR %s() got vector with norm 0!",FCT); 

   return x2;
};

template<class T>
wbarray<T>& wbarray<T>::ColProject(
   size_t k1, size_t k2, char nflag, char tnorm
){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k1>=dim2 || k2>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d,%d/%d)",FCT,k1,k2,dim2); 
   if (k1==k2) wblog(FL,"WRN %s() got k1=k2=%d",FCT,k1+1);

   gs_project_range(
      data+dim1*k1, data+dim1*k2, dim1, 1,
      nflag, tnorm
   );
   return *this;
};

template<class T>
bool wbarray<T>::isOrthogonalCol(size_t k0, char tnorm) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k0>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k0,dim2); 
   if (!dim1) return 0;

   const T *d=data, *d0=data+dim1*k0;
   T nrm2=overlap(d0,d0,dim1,1,tnorm);

   for (size_t k=0; k<dim2; k++, d+=dim1) {
      if (k!=k0 && ABS(overlap(d,d0,dim1,1,tnorm)/nrm2)>1E-10)
      return 0;
   }

   return 1;
};



template<class T>
unsigned wbarray<T>::QRdecomp(
   const char *F, int L, wbarray<T> &Q, wbarray<T> &R, T eps
){
   if (isEmpty()) { Q.init(); R.init(); return 0; }
   if (rank()!=2 || !numel()) wblog(F_L,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   size_t dim1=SIZE[0], dim2=SIZE[1], n=MIN(dim1,dim2);

   if (dim1<=1 || dim2<=1) {
      Q.init(dim1,n); R.init(n,dim2);

      if (!dim1 || !dim2) return n;
      if (dim1==1) { Q[0]=1; R=*this; return n; }
      if (dim2==1) { Q=*this; 
         if ((R[0]=this->norm())!=0) { Q*=(T(1)/R[0]); }
         else { Q[0]=1; }
         return n;
      }
   }

   if (allEqual(0)) {
      Q.init(dim1,1); R.init(1,dim2); Q[0]=1;
      return 1;
   }


   size_t i,j,k=0;
   WBPERM P;
   wbarray<T> U(dim1,n);

   R=*this;

   for (k=0; k<n; ++k) {
      if (R.Householder(F_L,k,U.col(k),&P,eps)<0) break;
   }

   if (k<n) {
      if (!k) wblog(FL,"ERR %s() resulted in k=%d/%d !??",FCT,k,n);
      R.Resize(k,dim2); n=k;
   }
   else if (n<dim1) {
      R.Resize(n,dim2);
   }

   Q.initIdentityB(dim1,n);


   T zero=0, ua=zero, *x;
   const T *u=U.col(n-1);

   for (k=n; k>0; u-=dim1) { --k; x=Q.data;
      for (j=0; j<n; ++j, x+=dim1) { 
         householder(u+k,x+k,dim1-k,eps);
      }
   }

   
   x=Q.data;
   for (j=0; j<n; ++j, x+=dim1) {
      for (i=0; i<dim1; ++i) { if (ABS(x[i])>eps) {
         if (x[i]<zero) { T *y=R.data+j;
            for (i=0; i<dim1; ++i) { x[i]=-x[i]; }
            for (k=0; k<dim2; ++k, y+=n) { y[0]=-y[0]; }
         }
         break;
      }}
   }

   return n;
};


template<class T>
wbarray<T>& wbarray<T>::ColPermute(const wbperm &P, char iflag){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);
   if (P.len!=SIZE[1]) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P.len,SIZE[1]);

   size_t i=0, dim1=SIZE[0], dim2=SIZE[1];
   wbarray<T> X(*this);

   for (; i<P.len; i++) {
      if (P[i]>=dim2) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,P[i],dim2); 
      if (iflag==0)
           MEM_CPY<T>(data+i*dim1, dim1, X.data+P[i]*dim1);
      else MEM_CPY<T>(data+P[i]*dim1, dim1, X.data+i*dim1);
   }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::SignCol(size_t k, const T *d0, char tnorm) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2); 

   T *d=data+dim1*k, x=overlap(d0,d,dim1,1,tnorm);
   if (x<0) {
      for (size_t i=0; i<dim1; i++) { d[i]=-d[i]; }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::FlipSignCol(size_t k) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2); 

   T *d=data+dim1*k;
   for (size_t i=0; i<dim1; i++) { d[i]=-d[i]; }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::SignConventionCol(
   size_t k0,
   double deps
){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,sizeStr().data);

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   size_t i,k,k1,k2; T *d, eps=
   ((typeid(T)!=typeid(double) && typeid(T)!=typeid(float)) ? 0 : deps);

   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2); 

   if (int(k0)>=0) { k1=k0; k2=k0+1; }
   else { k1=0; k2=dim2; }

   for (k=k1; k<k2; k++) { d=data+dim1*k;
      for (i=0; i<dim1; i++) {
         if (d[i]> eps) break;
         if (d[i]<-eps) {
             for (i=0; i<dim1; i++) { d[i]=-d[i]; }
             break;
         }
      }
   }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::NormalizeCols(
   const char *F, int L, double *amin, double *amax, char tnorm
){
   if (!SIZE.len) return *this;

   size_t i,j,m, n=SIZE[0];
   T *d=data, z=0; double a;
   wbvector<T> nn;

   if (SIZE.len>1) for (m=SIZE[1], i=2; i<SIZE.len; i++) m*=SIZE[i];
   else m=1;

   nn.init(m);

   for (j=0; j<m; j++, d+=n, z=0) {
      z=overlap(d,d,n,1,tnorm); z=SQRT(z);

      a=ABS(z); if (a==0) {
         wblog(F,L,"ERR |sum^2| of column returns %.3g!%s",a,
         tnorm && (typeid(T)==typeid(wbcomplex)) ? " (using t-norm)":"");
      }
      nn[j]=z;
   }

   if (amin) {
      a=nn.aMin(); if ((*amin)>0 && a<(*amin))
      wblog(F,L,"WRN %s() nMin=%.3g (%.3g)",FCT,a,*amin);
      (*amin)=a;
   }
   if (amax) {
      a=nn.aMax(); if ((*amax)>0 && a>(*amax))
      wblog(F,L,"WRN %s() nMax=%.3g (%.3g)",FCT,a,*amax);
      (*amax)=a;
   }

   for (d=data, j=0; j<m; j++, d+=n) 
   for (z=T(1)/nn[j], i=0; i<n; i++) { d[i]*=z; }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::OrthoNormalizeCols(
   const char *F, int L,
   char tnorm,
   char qflag,
   double eps,
   unsigned np
){
   if (!SIZE.len) return *this;

   size_t i,j,k,l,p,m, n=SIZE[0];
   T z, *d=0, *v=data; double a;

   char cflag=(typeid(T)==typeid(wbcomplex));
   char xflag=(qflag=='x' || qflag=='X');

   if (SIZE.len<2) m=1; else
   for (m=SIZE[1], i=2; i<SIZE.len; i++) m*=SIZE[i];

   if (!m || !n) return *this;

   for (l=k=0; k<m; ++k, ++l, v+=n) {
      for (p=0; p<np; ++p) {
         for (d=data, j=0; j<l; ++j, d+=n) {
            z=overlap(d,v,n,1,tnorm); if (z!=T(0)) {
            for (i=0; i<n; ++i) v[i]-=z*d[i]; }
         }
      }

      z=overlap(v,v,n,1,tnorm);
      a=ABS(z);
      
      if (a>eps) { a=SQRT(1/a);
         if (xflag && l<k)
              { for (i=0; i<n; ++i) d[i]=v[i]*a; }
         else { for (i=0; i<n; ++i) v[i]*=a; }
      }
      else {
         if (xflag) --l;
         else if (qflag) {
            for (i=0; i<n; ++i) v[i]=0;
         }
         else {
            wblog(F_L,"ERR |sum^2| of column returns %.3g!%s",
            a, tnorm  && cflag ? " (using t-norm)" : "");
         }
      }
   }

   if (xflag && l<k) {
      if (SIZE.len!=2) wblog(F_L,"ERR %s() qflag='%c' must be "
         "used with rank-2 arrays (%d)",FCT,qflag,qflag,SIZE.len);
      if (!l) wblog(F_L,"ERR %s() got all null-vectors (%g)",FCT,k);
      SIZE[1]=l;
   }

   return *this;
};


template<class T>
bool wbarray<T>::isSym_aux(
  const char *F, int L, const char *fct,
  const wbarray<T> &B, double eps, double* xref,
  const char symflag,
  const char lflag
) const {

   WBINDEX S1,S2;
   size_t i,j,k,dim1,dim2, r=SIZE.len, r2=r/2;
   double x;

   const wbarray<T> &A = (*this);
   bool issame = (&A==&B);
   wbarray<T> A2,B2;

   if (r%2 || r!=B.SIZE.len) {
      if (lflag) sprintf(str,
         "%s %s() only applies to even-rank objects (%ld;%ld)",
          shortFL(F,L),fct, r, B.SIZE.len);
      return 0;
   }
   if (isEmpty() && B.isEmpty()) return 1;

   dim1=dim2=1;
   for (i=0; i<r2; i++) { j=i+r2; dim1*=SIZE[i]; dim2*=SIZE[j];
      if (SIZE[i]!=B.SIZE[j] || SIZE[j]!=B.SIZE[i]) {
         if (lflag) sprintf(str,
           "%s %s() size mismatch [%s; %s].", shortFL(FL), fct,
            sizeStr().data, B.sizeStr().data);
         return 0;
      }
   }

   A.groupIndizes_DREF(A2,r2,r2);
   B.groupIndizes_DREF(B2,r2,r2);

   if (xref) *xref=0;
   if (eps<0) {
      eps=-eps; x=A2.aMax(); if (x>1) eps*=x;
   }

   if (symflag=='s') {
      if (eps==0) {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            if (A2(i,j)!=CONJ(B2(j,i))) {
               if (xref) {
                  x=ABS(A2(i,j)-CONJ(B2(j,i)));
                  *xref = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            x=ABS(A2(i,j)-CONJ(B2(j,i)));
            if (xref) *xref=MAX(*xref,x); else if (x>eps) return 0;
         }}
      }
   }
   else if (symflag=='a') {
      if (eps==0) {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            if (A2(i,j)!=-CONJ(B2(j,i))) {
               if (xref) {
                  x = ABS(A2(i,j) + CONJ(B2(j,i)));
                  *xref = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            x = ABS(A2(i,j) + CONJ(B2(j,i)));
            if (x>eps) { if (xref) *xref=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else wblog(FL,"ERR Invalid symflag=%c<%d>",symflag,symflag);

   if (xref && (*xref)>eps) return 0;

   return 1;
}


template<class T>
char wbarray<T>::hasGroupSize(size_t s1, size_t s2) const {

   size_t i,s=1,s0=1;

   for (i=0; i<SIZE.len; i++) {
      s*=SIZE[i]; if (s0==s1 && (!s || s!=s1)) break;
      s0=s;
   }
   for (s=1; i<SIZE.len; i++) s*=SIZE[i];

   return ((s0!=s1 || s!=s2) ? 0 : 1);
};


template<class T> inline
bool wbarray<T>::isComplex() const { return 0; };

template<> inline
bool wbarray<wbcomplex>::isComplex() const {
   for (size_t s=numel(), i=0; i<s; i++)
   if (data[i].i!=0.) return 1;

   return 0;
};


template<class T> inline
bool wbarray<T>::isZero(double eps, char bflag) const {

   size_t i, s=numel(); if (!s) return 1;

   if (eps==0.) {
      for (i=0; i<s; i++) if (data[i]!=0.) return 0;
   }
   else {
      if (bflag) {
         double e2=eps*eps, n2=0.;
         for (i=0; i<s; i++) {
         n2+=NORM2(data[i]); if (n2>e2) return 0; }
      }
      else {
         for (i=0; i<s; i++) if (ABS(data[i])>eps) return 0;
      }
   }

   return 1;
}


template<class T> inline
size_t wbarray<T>::nnz(
   const T eps, WBINDEX *I
) const {
   size_t i,n, s=numel(); T a;

   for (i=n=0; i<s; i++) {
      a=data[i]; if (a<0) a=-a;
      if (a>eps) n++;
   }

   if (I) {
      I->init(n);
      for (i=n=0; i<s; i++) {
         a=data[i]; if (a<0) a=-a;
         if (a>eps) I->data[n++]=i;
      }
   }

   return n;
}


template<class T> inline 
void wbarray<T>::SetRef (
   const INDEX_T *I, size_t len, const T* d0
){
   size_t i,k,s, S=numel(), r=SIZE.len, l=len-1;
   size_t nz=r-len;

   if (len>r || r==0 || (len && data==NULL)) wblog(FL,
      "ERR %s() index out of bounds ([%s]) for %s ???",FCT,
      (WBINDEX(len,I)+1).toStr().data, sizeStr().data);
   for (i=0; i<len; i++) if (I[i]>=SIZE[i+nz]) wblog(FL,
      "ERR %s() index out of bounds (%s; %s)",FCT,
      (WBINDEX(len,I)+1).toStr().data, sizeStr().data);

   if (!nz || !len) wblog(FL,"ERR invalid reference (len=%d/%d)",len,r);

   for (s=SIZE[0], i=1; i<nz; i++) s*=SIZE[i];
   for (k=I[l], i=l-1; i<l; i--) { k = k*SIZE[i+nz] + I[i]; }
   k*=s;

   if (k+s>S || d0==NULL) wblog(FL,
      "ERR index out of bounds ( %d+%d / %d; 0x%lX)\n"
      "having I=[%s] for %s array",k,s,S,d0,
      (WBINDEX(len,I)+1).toStr().data, sizeStr().data
   );

   MEM_CPY<T>(data+k,s,d0);
};


template<class T> inline
void wbarray<T>::squeeze() {

   size_t i,k=0,r=SIZE.len; if (r==0) return;

   for (i=0; i<r; i++) { if (SIZE[i]==1) continue;
      if (k!=i) SIZE[k]=SIZE[i];
      k++;
   }
   if (k<r) {
      if (k==0) k=1;
      SIZE.len=k;
   }
}

template<class T> inline
void wbarray<T>::skipSingletons() {

   size_t l,r=SIZE.len; if (r<=1) return;

   for (l=r-1; l>0; --l) {
      if (SIZE[l]==1) continue;
      else break;
   }

   if (l<SIZE.len) SIZE.len=l+1;
};

template<class T> inline
void wbarray<T>::addSingletons(
   unsigned r, const char *file, int line
){
   size_t i, r0=SIZE.len;

   if (!SIZE.len) wblog(FL,"ERR %s() called on empty array",FCT);

   if (r>SIZE.len) {
      SIZE.Resize(r);
      for (i=r0; i<r; i++) SIZE[i]=1;

      if (r0==0 && r) {
         wblog(FL,"WRN adding singletons to emtpy space (%d/%d)",r,r0);
         SIZE[0]=0;
      }
   }
   else if (r<SIZE.len) wblog(__FL__,
   "ERR %s() invalid singleton rank adjustment (%d/%d)",FCT,r,SIZE.len);
};

template<class T> inline
void wbarray<T>::prependSingletons(unsigned r) {

   if (!SIZE.len) wblog(FL,"ERR %s() called on empty array",FCT);

   if (r>SIZE.len) {
      size_t i, l=SIZE.len-1, m=r-SIZE.len;
      SIZE.Resize(r);
      for (i=l; i<r; i--) SIZE[i+m]=SIZE[i];
      for (i=0; i<m; i++) SIZE[i]=1;
   }
   else if (r<SIZE.len) wblog(FL,
   "ERR %s() invalid extended rank !?? (%d->%d)",FCT,SIZE.len,r);
}


template<class T> inline
double wbarray<T>::aMax(size_t *k) const {

   size_t i,s=numel(); double d, maxd=0; if (k) (*k)=0;

   if (!s && k) { 
      wblog(FL,"WRN %s() got empty array",FCT);
      (*k)=-1; return maxd;
   }

   for (i=0; i<s; i++) {
      d=(double)ABS(data[i]);
      if (d>maxd) { maxd=d; if (k) (*k)=i; }
   }
   return maxd;
};

template<class T> inline
double wbarray<T>::maxDiff (const wbarray<T> &B) const {
   size_t i,s=numel();
   double dmax=0;

   if (SIZE!=B.SIZE) {
      sprintf(str,"%s:%d data SIZE mismatch ([%s, %s])", FL,
      sizeStr().data, B.sizeStr().data); return NAN;
   }

   for (i=0; i<s; i++)
   dmax=MAX(dmax, ABS(data[i]-B.data[i]));

   return dmax;
}


template<class T> inline
void wbarray<T>::TimesEl(const wbarray<T> &B) {
   size_t i=0, n=numel();

   if (!hasSameSize(B)) {
      wblog(FL,"ERR %s() size mismatch: %s <> %s",FCT,
      sizeStr().data, B.sizeStr().data);
   }
   
   for (i=0; i<n; ++i) data[i] *= B.data[i];
};


template<class T>
wbarray<T>& wbarray<T>::timesEl(
   const wbarray &B, wbarray &C, T bfac, T cfac) const {

   size_t i=0, n=numel();

   if (!cfac || C.isEmpty()) {
      if (!hasSameSize(B)) wblog(FL,"ERR %s() size mismatch\n"
         "%s <> %s",FCT, sizeStr().data, B.sizeStr().data
      );
   
      if (bfac==0) C.init(SIZE);
      else { C=*this;
         if (bfac!=1)
              { for (; i<n; ++i) C.data[i] *= (B.data[i] * bfac); }
         else { for (; i<n; ++i) C.data[i] *= (B.data[i]); }
      }
   }
   else {
      if (bfac==0) return C;

      if (!hasSameSize(B) || !hasSameSize(C)) wblog(FL,
         "ERR %s() size mismatch\n%s <> %s <> %s",FCT,
         sizeStr().data, B.sizeStr().data, C.sizeStr().data
      );

      if (cfac!=1) C*=cfac;
      if (bfac!=1)
           { for (; i<n; ++i) C.data[i] += (data[i]*B.data[i] * bfac); }
      else { for (; i<n; ++i) C.data[i] += (data[i]*B.data[i]); }
   }

   return C;
};


template<class T> inline
T wbarray<T>::sumTimesEl(const wbarray<T> &B) const {
   size_t i,s=numel(); T x=0;

   if (!hasSameSize(B)) {
      wblog(FL, "ERR %s() size mismatch: %s <> %s",FCT,
      sizeStr().data, B.sizeStr().data);
   }
   
   for (i=0; i<s; i++) x += (data[i]*B.data[i]);
   return x;
}


template<class T> inline
T wbarray<T>::weightedAvg(const wbarray<T> &B) const {
   size_t i,s=numel(); T w,x=0,n2=0;

   if (!hasSameSize(B)) {
      wblog(FL, "ERR %s() size mismatch: %s <> %s",FCT,
      sizeStr().data, B.sizeStr().data);
   }
   
   for (i=0; i<s; i++) { w=B.data[i]; x+=(data[i]*w); n2+=w*CONJ(w); }

   if (n2<=0) { wblog(FL,
     "ERR weighted avg with total weight %.4g !??",n2);
      return x;
   }

   return x/std::sqrt(n2);
}


template<class T>
void wbarray<T>::tensorProd(
   const wbarray<T> &B0, wbarray<T> &C,
   const char aflag0, const char bflag0,
   const char kflag
) const {
   
   if (this==&C) {
      wbarray<T> X(*this); X.tensorProd(B0,C,aflag0,bflag0,kflag);
      return;
   }

   wbarray<T> A,B;
   opFlags<T> aflag(aflag0), bflag(bflag0);
   size_t ra=SIZE.len, rb=B0.SIZE.len;

   aflag.applyOrRef(FL,*this,A);
   bflag.applyOrRef(FL, B0,  B);


   if (ra!=rb) {
      if (ra>2 || rb>2) wblog(FL,
         "ERR %s() rank mismatch (%d,%d)",FCT,ra,rb);
      if (!ra || !rb) { C.init(); return; }

      if (ra==1) A.ExpandDiagonal(); else
      if (rb==1) B.ExpandDiagonal();
   }
   else if (ra!=rb) wblog(FL,
   "ERR %s() rank mismatch (%d,%d)",FCT,ra,rb);

   const size_t r=A.SIZE.len, l=r-1, *sa=A.SIZE.data, *sb=B.SIZE.data;
   WBINDEX Ia(r),Ib(r),Ja(r),Jb(r),S;

   size_t i,k,s,ia,ib;

   if (kflag) {
      S.init(r);
      for (i=0; i<r; i++) { S[i]=sa[i]*sb[i]; }
   }
   else {
      S.init(2*r);
      for (k=i=0; i<r; i++) { S[k++]=sa[i]; S[k++]=sb[i]; }
   }

   C.init(S); s=S.prod(); Ja[l]=Jb[l]=0;

   for (k=l,i=0; i<s; i++) {

       for (--k; k<r; --k) {
           Ja[k] = (Ja[k+1]+Ia[k+1])*sa[k];
           Jb[k] = (Jb[k+1]+Ib[k+1])*sb[k];
       }

       ia=Ja[0]+Ia[0];
       ib=Jb[0]+Ib[0];

       C.data[i] = A.data[ia] * B.data[ib];

       k=0; if ((++Ia[k])>=sa[k]) { Ia[k]=0; ++Ib[k]; }
       while (Ib[k]>=sb[k] && k<l) {
          Ib[k]=Ia[k]=0; k++;
          if ((++Ia[k])>=sa[k]) { Ia[k]=0; ++Ib[k]; }
       }
   }
}


template<class T>
wbarray<T>& wbarray<T>::Cat(
   unsigned dim,
   const wbvector< wbarray<T> > &aa
){
   wbvector< const wbarray<T>* > ap(aa.len);
   for (unsigned i=0; i<aa.len; i++) ap[i]=&aa[i];
   Cat(dim,ap); return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Cat(
   unsigned dim,
   const wbvector< wbarray<T>* > &ap
){
   if (ap.len==0) { init(); return *this; }
   if (ap.len==1) { *this=(*ap[0]); return *this; }

   size_t i,j=0,k,n, e=0, N, *s2=0;
   const size_t r=ap[0]->SIZE.len, l=r-1, *sref=ap[0]->SIZE.data;
   const WBINDEX &S0=ap[0]->SIZE;
   WBINDEX S(S0), I(r), I0(ap.len), idx;

   if (dim<1) wblog(FL,"ERR invalid index (not 1-based; %d/%d)",dim,r);
   dim--;

   if (dim<r) N=sref[dim];

   for (i=1; i<ap.len; i++) {
      if (ap[i]==NULL || ap[i]->data==NULL) continue;

      if (r != ap[i]->SIZE.len) e++; else {
         s2=ap[i]->SIZE.data;
         for (j=0; j<r; j++) if (j!=dim && s2[j]!=sref[j]) e++;
      }
      if (e) wblog(FL,"ERR size mismatch %d/%d: %s <> %s (%d)", i+1,ap.len,
      ap[0]->sizeStr().data,ap[i]->sizeStr().data,dim+1);

      I0[i]=N; if (dim<r) N+=s2[dim];
   }

   if (dim>=r) { T *dd;
      S.Resize(dim+1); S[dim]=ap.len; init(S); dd=data; n=ap[0]->numel();
      for (i=0; i<ap.len; i++, dd+=n) {
         MEM_CPY<T>(dd,n,ap[i]->data);
      }
      return *this;
   }

   idx.init(N);
   for (k=i=0; i<ap.len; i++) {
      if (ap[i]==NULL || ap[i]->data==NULL) continue;
      for (n=ap[i]->SIZE.data[dim], j=0; j<n; j++, k++) idx[k]=i;
   }

   S[dim]=N; init(S); n=numel();

   for (i=0; i<n; i++) {
       INDEX_T q=idx[I[dim]];
       const size_t *s=ap[q]->SIZE.data;

       for (j=0,k=l; k<r; k--) {
          if (j) j*=s[k];
          j += (k!=dim ? I[k] : (I[k]-I0[q]));
       }

       data[i]=ap[q]->data[j];

       k=0; I[0]++;
       while(I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::BlockCat(
   const wbarray< const wbarray<T>* > &ap,
   WBINDEX *D1_, WBINDEX *D2_
){
   size_t i,j,m,n,i0,j0,d1,d2,M,N;
   WBINDEX D1,D2;

   if (ap.isEmpty()) {
      init(); if (D1_) D1_->init(); if (D2_) D2_->init();
      return *this;
   }

   if (ap.rank()!=2) wblog(FL,
      "ERR %s() invalid input matrix (%s)",FCT,ap.sizeStr().data);

   m=ap.SIZE[0]; n=ap.SIZE[1]; D1.init(m); D2.init(n);
   for (i=0; i<m; i++)
   for (j=0; j<n; j++) if (ap(i,j) && !ap(i,j)->isEmpty()) {
      const wbarray<T> &a=(*ap(i,j));
      if (a.rank()!=2) wblog(FL,
         "ERR %s() invalid block rank (%d,%d): %s",
          FCT,i+1,j+1,a.sizeStr().data);

      if (D1[i]) {
         if (D1[i]!=a.SIZE[0]) wblog(FL,
         "ERR %s() block-size mismatch (%d,%d: %s <> %dx..",
          FCT, i+1, j+1, a.sizeStr().data, D1[i]);
      }
      else D1[i]=a.SIZE[0];

      if (D2[j]) {
         if (D2[j]!=a.SIZE[0]) wblog(FL,
         "ERR %s() block-size mismatch (%d,%d: %s <> ..x%d",
          FCT, i+1, j+1, a.sizeStr().data, D2[j]);
      }
      else D2[j]=a.SIZE[0];
   }

   if (D1_) D1_->init(); if (D2_) D2_->init();

   M=D1.sum(); N=D2.sum();
   init(M,N); if (isEmpty()) return *this;

   for (j0=j=0; j<n; j++, j0+=d2) { d2=D2[j];
   for (i0=i=0; i<m; i++, i0+=d1) { d1=D1[i];
       ap(i,j)->copyStride(ref(i0,j0),M);
   }}

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::blockTrace(
   size_t D1,
   wbarray<T> &a, size_t D2
) const {

   if (rank()!=2) wblog(FL,
   "ERR %s() only applies to rank-2 objects (%d)",FCT,rank());

   size_t i, dim1=SIZE[0], dim2=SIZE[1], m=(D1 ? dim1/D1 : 0);
   if (int(D2)<0) D2=D1;

   if (!D1 || !D2 || dim1%D1 || dim2%D2 || (dim2/D2)!=m) wblog(FL,
      "ERR %s() invalid block specs (%d/%d; %d/%d; %d)",
       D1,dim1, D2,dim2, m);

   a.init(D1,D2);
   for (i=0; i<m; i++) addBlock(i*D1, i*D2, D1, D2, a);
   
   return a;
};


template<class T>
wbarray<T>& wbarray<T>::BlockDiag(const wbvector< wbarray<T> > &D) {

   size_t i,d1=0,d2=0, i0=0, j0=0, D1=0, D2=0;

   for (i=0; i<D.len; i++) {
      if (!D[i].isMatrix()) wblog(FL,
         "ERR %s() matrizes expected (%s)",FCT,D[i].sizeStr().data);
      D1+=D[i].SIZE[0];
      D2+=D[i].SIZE[1];
   }

   init(D1,D2);
   for (i=0; i<D.len; i++, i0+=d1,j0+=d2) {
      d1=D[i].SIZE[0];
      d2=D[i].SIZE[1];
      Wb::cpyStride(ref(i0,j0), D[i].data,d1,d2,D1);
   }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::addBlock(
   size_t i, size_t j,
   size_t n, size_t m,
   wbarray<T> &a
) const {

   if (rank()!=2) wblog(FL,
   "ERR %s() only applies to rank-2 objects (%d)",FCT,rank());

   size_t dim1=SIZE[0], dim2=SIZE[1];

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d..%d/%d; %d..%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2); 

   if (a.isEmpty()) a.init(n,m);
   Wb::cpyStride(a.data, ref(i,j), n,m, -1,dim1,'+');

   return a;
};


template<class T>
wbarray<T>& wbarray<T>::getBlock(
   size_t i, size_t j,
   size_t n, size_t m,
   wbarray<T> &a
) const {

   if (rank()!=2) wblog(FL,
   "ERR %s() only applies to rank-2 objects (%d)",FCT,rank());

   size_t dim1=SIZE[0], dim2=SIZE[1];

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d..%d/%d; %d..%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2); 

   a.init(n,m); Wb::cpyStride(a.data, ref(i,j),n,m,-1,dim1);
   return a;
};


template<class T>
void wbarray<T>::copyStride(T* dd, size_t stride) const {

   if (rank()!=2) wblog(FL,
   "ERR %s() only applies to rank-2 objects (%d)",FCT,rank());

   size_t i=0, dim1=SIZE[0], dim2=SIZE[1];
   T *d0=data;

   for (; i<dim2; i++, d0+=dim1, dd+=stride) {
      MEM_CPY<T>(dd,dim1,d0);
   }
};


template<class T> inline
wbarray<T>& wbarray<T>::Reshape(const wbvector<size_t> &S) {

    if (numel() != S.prod(0)) {
       wblog(FL,"ERR %s() size not preserved (%s => %s)",
       FCT,sizeStr().data, S.toStrD().data);
       throw (char*)"err::size_mismatch";
    }

    SIZE=S; return *this;
};


template<class T>
void wbarray<T>::GroupIndizes(size_t K) {
    size_t i,n;

    if (!K || SIZE.len%K) wblog(FL,
    "ERR %s() invalid block specs: mod(%d,%d)",FCT,SIZE.len,K);

    n=SIZE.len/K;

    for (i=0; i<n; i++)
    SIZE[i]=prodRange(SIZE.data+i*K, K);

    SIZE.Resize(n);
}


template<class T> inline
void wbarray<T>::groupIndizes_P(
   const WBINDEX &I, int pos, wbperm &P,
   size_t *s1, size_t *s2
) const {

   size_t i,j,j1=0,j2=0,l, m=I.len, n=SIZE.len;
   PERM_T *p; INDEX_T *idx;
   char mark[n];
   
   P.init(n); p=P.data; idx=I.data;

   if (m>n) wblog(FL,
   "ERR %s() group index out of bounds (%d,%d)",FCT,m,n);

   if (pos==1) { j1=0;   j2=m; } else
   if (pos==2) { j1=n-m; j2=0; }
   else wblog(FL,"ERR %s() invalid pos=%d",FCT,pos);

   memset(mark,0,n*sizeof(char));

   for (j=j1, i=0; i<m; i++) { l=p[j++]=idx[i];
      if (l>=n || mark[l]) wblog(FL,
         "ERR %s() index out of bounds or not unique (%d/%d; %d)",
          FCT,l,n,mark[l]);
      mark[l]++;
   }
   for (j=j2, i=0; i<n; i++) { if (!mark[i]) p[j++]=i; }

   if (s1 || s2) {
      for (j1=j2=1, i=0; i<n; i++) {
      if (mark[i]) j1*=SIZE[i]; else j2*=SIZE[i]; }

      if (s1) *s1 = (pos==1 ? j1 : j2);
      if (s2) *s2 = (pos==1 ? j2 : j1);
   }
}


template<class T>
wbarray<T>& wbarray<T>::transpose(const char *F, int L, wbarray<T> &A) const {

    unsigned r=rank(); wbperm P;

    if (r%2) wblog(F_L,"ERR %s() "
       "applies to even-rank arrays only! (%s)",FCT,sizeStr().data);

    P.initTranspose(r);
    permute(A,P);
};


template<class T>
wbarray<T>& wbarray<T>::MatPermute(const wbperm &P, char iflag){

   if (isIdentityPerm(P)) return *this;
   if (!isSMatrix()) wblog(FL,
      "ERR %s() expects square matrix (%s)",FCT,sizeStr().data);

   wbarray<T> X(*this);
   const PERM_T *p=P.data;
   size_t i,j, dim1=SIZE[0], dim2=SIZE[1];

   if (!iflag) {
      for (i=0; i<dim1; i++)
      for (j=0; j<dim2; j++) { data[i+dim1*j]=X.data[p[i]+dim1*p[j]]; }
   }
   else {
      for (i=0; i<dim1; i++)
      for (j=0; j<dim2; j++) { data[p[i]+dim1*p[j]]=X.data[i+dim1*j]; }
   }

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::Permute(const WBPERM &P, char iflag){
    if (isIdentityPerm(P)) return *this;

    wbarray<T> A(*this); A.permute(*this,P,iflag);
    return *this;
}

template<class T>
wbarray<T>& wbarray<T>::Permute(const char* sidx, char iflag) {

    WBPERM P = Str2Idx(sidx,1);
    if (P.isEmpty()) {
        wblog (FL, "ERR Invalid index set %s", sidx);
        return *this;
    }

    return Permute(P,iflag);
}

template<class T>
wbarray<T>& wbarray<T>::permute(
    wbarray<T> &A,
    const char* sidx,
    char iflag
) const {

    WBPERM P = Str2Idx(sidx,1);

    if (P.isEmpty()) wblog(FL,"ERR Invalid index set %s", sidx);

    return permute(A,P,iflag);
}


#ifdef WB_CLOCK
   WbClock wbc_arr_perm("arr::permute");
#endif

template<class T>
wbarray<T>& wbarray<T>::permute(
    wbarray<T> &A,
    const WBPERM &P0,
    char iflag
) const {

    size_t i, j, k, *Sp, *S, r=SIZE.len, l=r-1, s, e=0;
    wbvector<size_t> Sp_;
    WBINDEX I_; INDEX_T *I; 
    WBPERM P; PERM_T *p;

#ifdef WB_CLOCK
   wbc_arr_perm.resume();
#endif

    if (&A==this) { wblog(FL,
       "ERR %s() MUST NOT permute onto itself.",FCT);
        e++;
    }

    if (iflag!=0 && iflag !='I') {
        wblog(FL,"ERR %s() invalid iflag '%c'<%d>",
        FCT,iflag,iflag); e++;
    }

    if (P0.len!=SIZE.len || e) {
        wblog(FL,
         "ERR %s() invalid permutation [%s; %d]",
          FCT,P0.toStr().data, SIZE.len);
    }

    if (!iflag) P=P0; else getIPerm(P0,P);

    if (!requiresDataPerm(P)) {
        A=(*this); A.SIZE.Select(P);
        return A;
    }

    p=P.data;


    SIZE.select(P,Sp_); A.init(Sp_); Sp=Sp_.data; s=numel();
    I_.init(SIZE.len); I=I_.data; S=SIZE.data;

    for (i=0; i<s; i++) {
        for (j=I[p[l]],k=l-1; k<r; k--) j = j*Sp[k] + I[p[k]];

        A.data[j]=data[i];

        k=0; I[0]++;
        while(I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
    }

#ifdef WB_CLOCK
   wbc_arr_perm.stop();
#endif

    return A;
}


template<class T>
wbarray<T>& wbarray<T>::select0(
    const WBPERM &P,
    unsigned dim,
    wbarray<T> &A
) const {

    if (&A==this) {
       wbarray<T> X(A); X.select0(P,dim,A);
       return A;
    }

    const unsigned r=SIZE.len, l=r-1; size_t i,j,k,d;
    WBINDEX S,I;

    if (dim>=SIZE.len) wblog(FL,
       "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

    S=SIZE; d=SIZE[dim]; S[dim]=P.len;

    for (i=0; i<P.len; i++) if (P[i]>=d) wblog(FL,
        "ERR %s() index out of bounds (%d/%d)",FCT,P[i],d);

    A.init(S); I.init(SIZE.len);
    INDEX_T *ip=I.data, *s=S.data; const PERM_T *p=P.data;

    for (d=A.numel(), i=0; i<d; i++) {
        for (j=0, k=l; k<r; k--) {
           if (j) j*=SIZE[k];
           j+=(k!=dim ? ip[k] : p[ip[k]]);
        }

        A.data[i]=data[j];

        k=0; ip[0]++;
        while(ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
    }

    return A;
}


template<class T>
wbarray<T>& wbarray<T>::selectSqueeze(
    wbarray<T> &A, unsigned p,
    unsigned dim
) const {

    if (&A==this) {
       wbarray<T> X(A); X.selectSqueeze(A,p,dim);
       return A;
    }

    const unsigned r=SIZE.len, l=r-1; size_t i,j,k,s;
    WBINDEX S,I;

    if (dim>=SIZE.len) wblog(FL,
    "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

    if (p>=SIZE[dim]) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,p,SIZE[dim]);

    if (SIZE.len==1) {
       A.init(1); A.data[0]=data[p];
       return A;
    }

    S=SIZE; S[dim]=1; A.init(S); s=A.numel(); I.init(SIZE.len);

    for (i=0; i<s; i++) {
        for (j=0, k=l; k<r; k--) {
           if (j) j*=SIZE[k];
           j+=(k!=dim ? I[k] : p);
        }

        A.data[i]=data[j];

        k=0; I[0]++;
        while(I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
    }

    for (i=dim+1; i<S.len; i++) A.SIZE[i-1]=A.SIZE[i];
    A.SIZE.len--; A.SIZE[A.SIZE.len]=0;

    return A;
}


template<class T>
void wbarray<T>::setBlock(
    const wbvector< WBINDEX > &DB,
    const wbindex &IB,
    const wbarray<T> &A 
){
    size_t i,j,ik,k,s, r=SIZE.len, l=r-1;
    wbindex I,I1(r),I2(r);

    if (DB.len!=r || IB.len!=r || A.SIZE.len!=r)
        wblog(FL,"ERR %s() size mismatch (%d %d %d; r)",
        FCT, DB.len, IB.len, A.SIZE.len, r);

    for (k=0; k<r; k++) {
        const WBINDEX &Dk=DB[k]; ik=IB[k];

        if (ik>=Dk.len) wblog(FL,
           "ERR %s() block index out of bounds (%d: %d/%d)",
            FCT,k+1, ik+1, Dk.len);
        if (Dk[ik] != A.SIZE[k]) wblog(FL,
           "ERR %s() block size mismatch (%d: %d/%d)",
            FCT,k+1, Dk[ik], A.SIZE[k]);

        I1[k]=addRange(Dk.data,ik);
        I2[k]=I1[k]+Dk[ik];

        if (I2[k]>SIZE[k]) wblog(FL,
           "ERR %s() index out of bounds (%d: %d/%d)",
            FCT,k+1, I2[k], SIZE[k]);
    }

    I=I1; s=A.numel();

    for (i=0; i<s; i++) {
        for (j=I[l],k=l-1; k<r; k--) j = j*SIZE[k] + I[k];

        data[j]=A.data[i];

        k=0; I[0]++;
        while(I[k]>=I2[k] && k<l) { I[k]=I1[k]; ++I[++k]; }
    }
}


template<class T>
void wbarray<T>::addBlock(
    const WBINDEX &I0,
    const wbarray<T> &A,
    const char dflag
){
    size_t i,j,k, s=A.numel(), r=SIZE.len, l=r-1;
    unsigned ra=(dflag ? r/2 : r), la=ra-1;
    const INDEX_T *S=A.SIZE.data;
    WBINDEX I(ra);

    if (I0.len!=r || A.SIZE.len!=ra || dflag && r%2) wblog(FL,
       "ERR %s() size mismatch (%d/%d, %s; %c<%d>)",
        FCT,I0.len, A.SIZE.len, sizeStr().data, dflag
    );
    for (i=0; i<r; i++) if (I0[i]+S[i>=ra ? i-ra : i]>SIZE[i])
         wblog(FL,"ERR %s() index out of bounds (%d: %s - %s %s)",
         FCT,i+1,I0.toStr().data,A.sizeStr().data,sizeStr().data
    );

    for (i=0; i<s; i++) {
       for (j=0,k=l; k<r; k--) { if (j) j*=SIZE[k];
          j+=(I0[k]+I[k<ra ? k : k-ra]);
       }

       data[j] += A.data[i];

       k=0; I[0]++;
       while(I[k]>=S[k] && k<la) { I[k]=0; ++I[++k]; }
    }
};


template<class T>
T wbarray<T>::min() const {
   size_t i,n=numel(); T x;

   if (n==0) {
       wblog(FLINE,"WRN min() of empty wbarray ???");
       memset(&x,0,sizeof(T)); return x;
   }

   for (x=data[0], i=1; i<n; i++)
   if (x>data[i]) x=data[i];

   return x;
}

template<class T>
T wbarray<T>::max() const {
   size_t i,n=numel(); T x;

   if (n==0) {
       wblog(FLINE,"WRN max() of empty wbarray ???");
       memset(&x,0,sizeof(T)); return x;
   }

   for (x=data[0], i=1; i<n; i++)
   if (x<data[i]) x=data[i];

   return x;
}


template<class T>
T wbarray<T>::aMin(char zflag, size_t *k_) const {

   size_t i,k=-1,n=numel(); T a,x=0;

   if (!n) {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got empty array",FCT); }
      return 0;
   }

   if (zflag==0) {
      x=ABS(data[i++]); k=0; if (x!=0) {
         for (; i<n; ++i) { a=ABS(data[i]); if (x>a) {
            x=a; k=i; if (x==0) break;
         }}
      }
   }
   else {
      for (; i<n; ++i) {
         a=ABS(data[i]); if (a) { x=a; k=i; break; }
      }
      for (; i<n; ++i) { a=ABS(data[i]);
         if (a && x>a) { x=a; k=i; }
      }
   }
   if (k_) (*k_)=k;
   return x;
};


template<class T>
void wbarray<T>::print_ref (const char *mark, const char *newl) const {
   int n=nrefs(); char isref=isRef();
   if (isref || n) {
      if (mark && mark[0]) printf("%s",mark);
      if (isref) printf(" ISREF/%d",isref);
      if (n) printf(" %d XREF%s",n, n!=1 ? "S" : "");
      if (mark && mark[0]) printf("%s",mark);
   }; if (newl && newl[0]) printf("%s",newl);
};


template<class T>
void wbarray<T>::info(const char *F, int L, const char* istr) const {

   wbstring tstr; size_t s=numel();

   if (typeid(T)==typeid(double)) tstr="double";
   else tstr.init(getName(typeid(T)).data);

   printf("  %-8s %-8s %s %s",
      istr, sizeStr().data, tstr.data, s>1 ? " array" : ""
   ); print_ref();
}


template<class T>
void wbarray<T>::info(const char* istr) const {

   unsigned l=0, n=16; char ss[n];
   size_t s=numel(), b=s*sizeof(T);

   wbstring tstr;
   if (typeid(T)==typeid(double)) tstr="double";
   else tstr.init(getName(typeid(T)).data);

   if (b<(1<<10)) l=snprintf(ss,n,"%ld  ",b); else
   if (b<(1<<20)) l=snprintf(ss,n,"%.3g kB",b/double(1<<10)); else
                  l=snprintf(ss,n,"%.3g MB",b/double(1<<20));
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   printf("  %-12s %-10s %12s  @ 0x%7pX  %s%s",
      istr, sizeStr().data, ss, data,
      tstr.data, s>1 ? " array" : ""
   );

   print_ref();
};


template<class T>
void wbarray<T>::print(
   const char *F, int L, const char* istr, const char* fmt0
 ) const {

   if (F) {
      wbstring tstr;
      if (typeid(T)==typeid(double)) tstr="double";
      else tstr.init(((char*)(typeid(T).name())));

      wblog(F,L,"%s = [%s] %s array \\", istr[0] ? istr : "ans",
      sizeStr().data, tstr.data); print_ref();
   }
   else info(istr);

   printdata(istr, fmt0);
};


template<class T>
int wbarray<T>::printdata(
   const char *istr,
   const char *fmt0
) const {

   unsigned i,j,k, r=SIZE.len, l=r-1, s=numel();
   wbstring fmt;

   if (!fmt0[0]) {
       if (typeid(T)==typeid(double) || typeid(T)==typeid(float))
            fmt=" %12.5g";
       else fmt.init2Fmt((T)0, 6);
   } else fmt=fmt0;

   if (s>1) {
      if (SIZE.len<2) {
         printf("\n");
         for (i=0; i<s; i++) printf(fmt.data,data[i]);
         if (s) printf("\n\n");
      }
      else {
         WBINDEX I(r);

         wbarray<T> A; wbperm P(r); P[0]=1; P[1]=0;
         permute(A,P);

         for (k=2, i=0; i<s; i++) {
            if (k) { if (i) printf("\n");
               if (k>1) {
                  if (r>2) {
                     printf("\n  %s(:,:",istr);
                     for (j=2; j<r; j++) printf(",%d",I[j]+1);
                     printf(") = [\n\n");
                  }
                  else printf("%s = [\n",istr);
               }
            }

            printf(fmt.data,A.data[i]);

            k=0; I[0]++;
            while(I[k]>=A.SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
         }
         printf("\n];\n");
      }
   }
   else if (s==1) {
      printf("  %s = ",istr);
      printf(fmt.data,data[0]); printf("\n");
   }
   else printf("\n");

   return (int)s;
};


template<>
int wbarray<wbcomplex>::printdata(
   const char *istr, const char *fmt0
) const {

   unsigned i,j,k, r=SIZE.len, l=r-1, s=numel();
   wbstring fmt;

   fmt=((fmt0 && fmt0[0]) ? fmt0 : " %8.4g%+8.4gi");

   if (s>1) {
      if (SIZE.len<2) {
         printf("\n");
         for (i=0; i<s; i++) printf(fmt.data,data[i].r,data[i].i);
         if (s) printf("\n\n");
      }
      else {
         WBINDEX I(r);

         wbarray<wbcomplex> A; wbperm P(r); P[0]=1; P[1]=0;
         permute(A,P);

         for (k=2, i=0; i<s; i++) {
            if (k) { if (i) printf("\n");
               if (k>1) {
                  if (r>2) {
                     printf("\n  %s(:,:",istr);
                     for (j=2; j<r; j++) printf(",%d",I[j]+1);
                     printf(") = [\n\n");
                  }
                  else printf("%s = [\n",istr);
               }
            }

            printf(fmt.data,A.data[i].r,A.data[i].i);

            k=0; I[0]++;
            while(I[k]>=A.SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
         }
         printf("\n];\n");
      }
   }
   else if (s==1) {
      printf("  %s = ",istr);
      printf(fmt.data,data[0].r,data[0].i); printf("\n");
   }
   else printf("\n");

   return (int)s;
};


template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::comm(
  const wbarray<TB> &B, wbarray<TC> &C, char aflag, char bflag
) const {

   if ((void*)&C==(void*)this || (void*)&C==(void*)&B) {
      wbsparray<TC> X; comm(B,X,aflag,bflag);
      return X.save2(C);
   }

   if (!isMatrix() || !B.isMatrix() || SIZE[0]!=B.SIZE[1] || 
       SIZE[0]!=SIZE[1] || B.SIZE[0]!=B.SIZE[1]) wblog(FL,
      "ERR %s() invalid operators for [A,B] (%s; %s)",
       FCT,sizeStr().data, B.sizeStr().data
   );

   wbMatProd(*this,B,C,aflag,bflag);
   wbMatProd(B,*this,C,bflag,aflag,-1.,1.);

   return C;
};

template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::acomm(
  const wbarray<TB> &B, wbarray<TC> &C, char aflag, char bflag
) const {

   if (!isMatrix() || !B.isMatrix() ||
       SIZE[0]!=SIZE[1] || B.SIZE[0]!=B.SIZE[1]) wblog(FL,
      "ERR %s() invalid operators for [A,B] (%s; %s)",
       FCT,sizeStr().data, B.sizeStr().data
   );

   wbMatProd(*this,B,C,aflag,bflag);
   wbMatProd(B,*this,C,bflag,aflag,+1.,1.);

   return C;
};


template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::contractMat(
   unsigned i1, const wbarray<TB> &B, wbarray<TC> &Q,
   unsigned i2
 ) const {

   unsigned r=SIZE.len, l=r-1;
   wbvector<size_t> S; wbperm P;
   wbarray<T> A;
   wbarray<T> MA;

   if (!B.isMatrix()) wblog (FL,
      "ERR %s() rank 2 object required (%d)",FCT,B.SIZE.len);
   if (i1<1 || i1>SIZE.len || i2<1 || i2>2) wblog(FL,
      "ERR %s() invalid contraction indizes [%d %d, %d %d]", i1,i2,
       FCT,rank(), B.rank());

   if (i1==1) {
       toMatrixRef(MA,i1-1,1,P); SIZE.select(P,S); l=0;
       if (i2==2)
            { wbMatProd(B,MA,Q    ); S[l]=B.SIZE[0]; }
       else { wbMatProd(B,MA,Q,'T'); S[l]=B.SIZE[1]; }
   }
   else {
       toMatrixRef(MA,i1-1,2,P); SIZE.select(P,S);
       if (i2==1)
            { wbMatProd(MA,B,Q        ); S[l]=B.SIZE[1]; }
       else { wbMatProd(MA,B,Q,'N','T'); S[l]=B.SIZE[0]; }
   }

   Q.Reshape(S).Permute(P,'I');

   return Q;
}


template<class T>
template<class TB> inline
wbarray<T>& wbarray<T>::Contract(char *idx1,
   const wbarray<TB> &B, char *idx2,
   const wbperm &pfinal
){
   wbarray<T> A(*this);

   WBINDEX i1, i2;
   i1 = Str2Idx(idx1,1);
   i2 = Str2Idx(idx2,1);

   A.contract(i1, B, i2, *this, pfinal);
   return *this;
}

template<class T>
template<class TB, class TC> inline
wbarray<TC>& wbarray<T>::contract(const char *F, int L,
   const char *idx1, const wbarray<TB> &B,
   const char *idx2, wbarray<TC> &Q,
   const wbperm &pfinal) const
{
   WBINDEX i1, i2;
   i1 = Str2Idx(idx1,1);
   i2 = Str2Idx(idx2,1);

   return contract(F,L, i1,B,i2,Q,pfinal);
}

template<class T>
template<class TB> inline
wbarray<T>& wbarray<T>::Contract(const wbvector<unsigned> &i1,
   const wbarray<TB> &B, const wbvector<unsigned> &i2,
   const wbperm &pfinal)
{
   wbarray<T> A(*this); A.contract(i1, B, i2, *this, pfinal);
   return *this;
}


#ifdef WB_CLOCK
   WbClock wbc_arr_ctr("arr::contract");
#endif

template<class TA>
template<class TB, class TC>
wbarray<TC>& wbarray<TA>::contract(
   const char *F, const int L,
   const wbvector<unsigned> &i1, const wbarray<TB> &B0,
   const wbvector<unsigned> &i2,
   wbarray<TC> &C0, const wbperm &P, const TA& afac
) const {

#ifdef WB_CLOCK
   wbc_arr_ctr.resume();
#endif

   if ((void*)(&C0)==(void*)this || (void*)(&C0)==(void*)(&B0))
   wblog(F,L,"ERR contract must have distinct target space");

   unsigned i,s; size_t k;
   char aflag, bflag;
   wbvector<size_t> S1,S2;
   wbarray<TA> MA;
   wbarray<TB> MB;
   wbarray<TC> Q;
   wbperm p1,p2;

   if (i1.isEmpty() || i2.len!=i1.len) wblog(F,L,"ERR %s() "
      "invalid index set {[%s], [%s]}",FCT,i1.toStr().data,i2.toStr().data);
   for (s=SIZE.len, i=0; i<i1.len; ++i) {
      if (i1[i]>=s) wblog(F,L,"ERR %s() "
      "index out of bounds ([%s]/%d)",FCT,(i1+1).toStr().data,s);
   }
   for (s=B0.SIZE.len, i=0; i<i2.len; ++i) {
      if (i2[i]>=s) wblog(F,L,"ERR %s() "
      "index out of bounds ([%s]/%d)",FCT,(i2+1).toStr().data,s);
   }
   for (i=0; i<i1.len; ++i) {
      if (SIZE[i1[i]]!=B0.SIZE[i2[i]]) wblog(F,L,"ERR wbarray::%s() "
      "incompatible object\n[%s] (%d) <> [%s] (%d)",
      FCT,sizeStr().data,i1[i]+1,B0.sizeStr().data,i2[i]+1);
   }

      toMatrixRef(FL,MA,i1,2,aflag);     SIZE.getI(i1,S1);
   B0.toMatrixRef(FL,MB,i2,1,bflag);  B0.SIZE.getI(i2,S2); k=1;

   if (S1.isEmpty() && S2.len<2) S1.init(1,&k);
   if (S2.isEmpty() && S1.len<2) S2.init(1,&k);

   if (afac==TA(0)) {
      WBINDEX S(S1,S2);
      if (!P.isEmpty()) S.Permute(P);

      if (C0.isEmpty()) C0.init(S);
      else if (C0.SIZE!=S) wblog(F,L,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
          FCT, shortFL(FL), S.toStrD().data, C0.sizeStr().data
      );
#ifdef WB_CLOCK
      wbc_arr_ctr.stop();
#endif
      return C0;
   }


   if (!P.isTranspose(S1.len,S2.len)) {
      wbMatProd(MA,MB,Q, (aflag?'T':'N'), (bflag?'T':'N'));

      Q.Reshape(UVEC(S1,S2));"Q");
      if (!P.isEmpty()) Q.Permute(P);
   }
   else {
      wbMatProd(MB,MA,Q, (bflag?'N':'T'), (aflag?'N':'T'));
   }

   if (afac!=TA(1)) Q*=afac;

   if (C0.isEmpty()) Q.save2(C0);
   else {
      if (!Q.hasSameSize(C0)) wblog(F,L,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
          FCT, shortFL(FL), Q.sizeStr().data, C0.sizeStr().data
      );
      C0+=Q;
   }

#ifdef WB_CLOCK
   wbc_arr_ctr.stop();
#endif

   return C0;
}


template<class T>
wbarray<T>& wbarray<T>::contract(
   unsigned i1, unsigned i2,
   wbarray<T> &C
) const {

   size_t i,j,k,d,q,s,r,l;
   WBINDEX S,I0,I;

   if (&C==this) wblog (FL,
   "ERR %s() output space coincides with *this",FCT);

   if (i1==i2 || !i1 || !i2 || i1>SIZE.len || i2>SIZE.len ||
       SIZE[i1-1]!=SIZE[i2-1]) {
       wblog (FL,"ERR %s() invalid trace indizes (%d,%d; %s)",
       FCT,i1, i2, sizeStr().data); return C;
   }

   i1--; i2--;

   d=SIZE[i1];

   S.init(SIZE.len-2);
   for (k=i=0; i<SIZE.len; ++i) if (i!=i1 && i!=i2) S[k++]=SIZE[i];
   if (S.len==0) { S.init(1); S[0]=1; }

   if (C.isEmpty()) C.init(S); else {
      if (S!=C.SIZE) wblog(FL,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
          FCT,shortFL(FL),C.sizeStr().data, S.toStr().data
      );
   }

   I.init(S.len); s=S.prod(0); I0.init(SIZE.len);
   r=SIZE.len-1; l=S.len-1;

   for (i=0; i<s; ++i) {
      for (q=k=0; k<SIZE.len; ++k) if (k!=i1 && k!=i2) I0[k]=I[q++];
      for (q=0; q<d; ++q) {
          I0[i1]=I0[i2]=q;
          for (j=I0[r],k=r-1; k<r; k--) j = j*SIZE[k] + I0[k];

          C.data[i]+=data[j];
       }

       k=0; I[0]++;
       while(I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return C;
}


template<class T>
wbarray<T>& wbarray<T>::contract(
   const wbMatrix<unsigned> &I12,
   wbarray<T> &C
) const {

   size_t i,j,k,l,lt,l0,q,s,st,e=0;
   WBINDEX S, St, I, I0, It; wbvector<int> mark;

   if (&C==this) wblog (FL,
   "ERR %s() output space coincides with *this",FCT);

   if (I12.dim2!=2 || I12.dim1*I12.dim2>SIZE.len) {
       I12.Print("index set");
       wblog (FL,"ERR %s() invalid trace index set (%dx%d; %d)",
       FCT,I12.dim1, I12.dim2, SIZE.len);
   }

   mark.init(SIZE.len);
   St.init(I12.dim1);

   for (i=0; i<I12.dim1; ++i) {
      for (j=0; j<I12.dim2; ++j) {
         k=I12(i,j); if (k>=SIZE.len || mark[k]++) ++e;
      }
      if (!e) {
         St[i]=SIZE[I12(i,0)];
         if (St[i]!=SIZE[I12(i,1)]) ++e;
      }
      if (e) {
          I12.print("I12");
          wblog (FL,"ERR invalid trace index set (%s)",
          sizeStr().data); return C;
      }
   }

   if (I12.dim1==0) {
      if (C.isEmpty()) C=(*this); else C+=(*this);
      return C;
   }

   S.init(SIZE.len-2*I12.dim1);
   for (k=i=0; i<SIZE.len; ++i) if (!mark[i]) S[k++]=SIZE[i];
   if (S.len==0) { S.init(1); S[0]=1; }

   if (C.isEmpty()) C.init(S); else {
      if (S!=C.SIZE) wblog(FL,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
          FCT,shortFL(FL),C.sizeStr().data, S.toStr().data
      );
   }

   I.init(S.len); s=S.prod(0); I0.init(SIZE.len); st=St.prod(0);
   l=S.len-1; lt=St.len-1; l0=SIZE.len-1;

   for (i=0; i<s; ++i) {
       for (q=k=0; k<SIZE.len; ++k) if (!mark[k]) I0[k]=I[q++];

       It.init(St.len);
       for (q=0; q<st; ++q) {
          for (k=0; k<It.len; ++k) I0[I12(k,0)]=I0[I12(k,1)]=It[k];

          for (j=I0[l0],k=l0-1; k<l0; k--) j = j*SIZE[k] + I0[k];

          C.data[i] += data[j];

          k=0; It[0]++;
          while(It[k]>=St[k] && k<lt) { It[k]=0; ++It[++k]; }
       }

       k=0; I[0]++;
       while(I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return C;
}


template<class T>
wbarray<T> TestContract(
   const wbarray<T> &A, const WBINDEX i1,
   const wbarray<T> &B, const WBINDEX i2,
   const wbperm &pfinal
){
   WBINDEX Si,Sj,Sk, I,J,K, Ia,Ib,Ic, Iq; WBPERM Pa,Pb;
   size_t i,j,k, si,sj,sk, li,lj,lk, q, ra=A.SIZE.len, rb=B.SIZE.len;
   wbarray<T> C;

   T dbl;

   if (i1.isEmpty() || i2.len!=i1.len) wblog (FL,
      "ERR %s() invalid index set {[%s], [%s]}",
       FCT,i1.toStr().data, i2.toStr().data);
   if (i1.anyGE(ra) || i2.anyGE(rb)) wblog (FL,
      "ERR %s() index out of range\n[%s; %d], [%s; %d]",
       FCT,i1.toStr().data,ra,i2.toStr().data,rb);

return A;

   for (k=0; k<i1.len; ++k) {
   if (A.SIZE[i1[k]]!=B.SIZE[i2[k]]) wblog (FL,
      "ERR %s() incompatible objects {[%s], [%s]}\n[%s] <> [%s] at %d",
       FCT,i1.toStr().data, i2.toStr().data,
       A.sizeStr().data, B.sizeStr().data, k);
   }

   Si=A.SIZE.Skip(i1); si=Si.len>0 ? Si.prod(0) : 1; li=Si.len-1;
   Sj=A.SIZE [i1];     sj=Sj.len>0 ? Sj.prod(0) : 1; lj=Sj.len-1;
   Sk=B.SIZE.Skip(i2); sk=Sk.len>0 ? Sk.prod(0) : 1; lk=Sk.len-1;

   Pa = Index(0,ra-1).Move2end  (i1);
   Pb = Index(0,rb-1).Move2front(i2);

   Ia.init(ra); Ib.init(rb); C.init(Si.append(Sk));

   if (C.numel()==0) C.init(1);

   I.init(Si.len);
   for (i=0; i<si; ++i) {

       K.init(Sk.len);
       for (k=0; k<sk; ++k) {

           J.init(Sj.len);
           for (dbl=0, j=0; j<sj; ++j) {
               Iq=I.append(J); for (q=0; q<Iq.len; ++q) Ia[Pa[q]]=Iq[q];
               Iq=J.append(K); for (q=0; q<Iq.len; ++q) Ib[Pb[q]]=Iq[q];

               if (j==sj-1) Ic=I.append(K);

               dbl += A(Ia) * B(Ib);

               if (J.len>0) {
                   q=0; J[0]++;
                   while(J[q]>=Sj[q] && q<lj) { J[q]=0; ++J[++q]; }
               }
           }
           C(Ic) = dbl;

           if (K.len>0) {
               q=0; K[0]++;
               while(K[q]>=Sk[q] && q<lk) { K[q]=0; ++K[++q]; }
           }
       }
       if (I.len>0) {
           q=0; I[0]++;
           while(I[q]>=Si[q] && q<li) { I[q]=0; ++I[++q]; }
       }
   }

   if (pfinal) C.Permute(pfinal);

   return C;
}


bool mxIsWbarray(
   const char *F, int L, const mxArray *a,
   const double **dr, const double **di, const mwSize **sz
){
   if (!a) { return 0; }

   if (!mxIsDouble(a) || mxIsSparse(a)) { if (dr) wblog(F_L,
      "ERR %s() invalid input array (got %s; %d/%d)",FCT,
       mxGetClassName(a), mxIsDouble(a), mxIsSparse(a)); return 0; }

   if (dr) (*dr)=mxGetPr(a);
   if (di) (*di)=mxGetPi(a);
   if (sz) { (*sz)=mxGetDimensions(a); 
      if ((*sz)==NULL) { if (dr) wblog(F_L,
         "ERR mxGetDimensions returned (%lX, dim=%d) ???",
          *sz,mxGetNumberOfDimensions(a)); return 0;
      }
   }

   return 1;
};


template<class T> inline
bool wbarray<T>::requiresDataPerm(const wbperm &P) const {

   const PERM_T *p=P.data; const size_t *s=SIZE.data;
   unsigned i,k,l;

   if (P.len!=SIZE.len) wblog(FL,
      "ERR %s() invalid array permutation (%d,%d)",
       FCT,P.len,SIZE.len);

   if (numel()<=1) return 0;

   for (k=i=0; i<P.len; ++i) {
      l=p[i]; if (s[l]==1) continue;
      if (k>l) return 1;
      k=l;
   }

   return 0;
}


template<class T> inline
void wbarray<T>::toMatrixRef(
   wbarray<T> &A,
   const wbvector<unsigned> &I,
   int pos,
   wbperm &P
) const {

   size_t s1,s2;
   groupIndizes_P(I,pos,P,&s1,&s2);

   if (!requiresDataPerm(P)) { A.init2ref(*this);  }
   else { permute(A,P); }

   A.Reshape(s1,s2); 
}


template<class T> inline
bool wbarray<T>::toMatrixRef(const char *F, int L,
   wbarray<T> &A,
   const wbvector<unsigned> &I,
   int pos,
   char &tflag
) const {

   size_t i,s1,s2; wbperm P;

   int k[2]={pos, pos==1 ? 2 : 1};

   for (i=0; i<2; ++i) {
      groupIndizes_P(I,k[i],P,&s1,&s2);

      if (!requiresDataPerm(P)) { 
         tflag = (pos==k[i] ? 0 : 1);
         A.init2ref(*this); 
         A.Reshape(s1,s2);
         return 0;
      }
   }

   if (&A==this) wblog(F_L,"ERR %s() overlaping output space",FCT);

   toMatrixRef(A,I,pos,P); tflag=0;

   return 1;
}


template<class T>
void wbarray<T>::toMatrixRef(
   wbarray<T> &A,
   unsigned K,
   const wbperm &P0
) const {

   size_t i,s1,s2; wbperm P;

   if (K>SIZE.len) wblog(FL,
      "ERR %s() index out of bounds (%d,%d)",
       FCT, K, SIZE.len);

   if (P0.isEmpty()) P.Index(SIZE.len);
   else {
      if (P0.len!=SIZE.len) wblog(FL,
         "ERR %s() invalid permutation (%d,%d)",FCT,P0.len,SIZE.len);
      P=P0;
   }

   for (s1=1,i=0; i<K; ++i) s1*=SIZE[P[i]];
   for (s2=1; i<P.len; ++i) s2*=SIZE[P[i]];

   if (!requiresDataPerm(P)) { A.init2ref(*this);  }
   else { permute(A,P); }

   A.Reshape(s1,s2); 
}







template<class T> inline
void wbarray<T>::save2(wbMatrix<T> &M) {

   if (M.data) M.init();

   if (nrefs()) wblog(FL,
      "ERR %s() array has %d reference(s) !!",FCT,nrefs(),nrefs());

   if (SIZE.len==2) {
      M.data=data; M.dim1=SIZE[1]; M.dim2=SIZE[0]; M.isref=isRef();
      SIZE.init(); data=NULL;
      if (this->sptr_flags) WB_DELETE(this->sptr_flags);
   }
   else if (SIZE.len) wblog(FL,"ERR invalid rank-%d",SIZE.len);
};

template<class T> inline
wbarray<T>& wbarray<T>::save2(wbarray<T> &A) {

   if (this!=&A) {
      SIZE.save2(A.SIZE);
      A.data=data; data=NULL;
      A.sptr_flags=this->sptr_flags; this->sptr_flags=NULL;
   }
   return A;
};














#ifdef USE_WB_MPFR

template<> inline
mxArray* wbarray<Wb::quad>::toMx() const {

   const char *fields[]={"S","data"};
   mxArray *a=mxCreateStructMatrix(1,1,2,fields);

   wbvector<Wb::quad> X(numel(),data,'r');

   mxSetFieldByNumber(a,0,0, SIZE.toMx());
   mxSetFieldByNumber(a,0,1, X.toMx());

   return a;
};

template<> inline
wbarray<Wb::quad>& wbarray<Wb::quad>::init(
   const char *F, int L, const mxArray *S, char ref, char vec) {

   wbvector<size_t> sz(FL, mxGetFieldByNumber(S,0,0));
   wbvector<Wb::quad> X(mxGetFieldByNumber(S,0,1));

   if (X.len!=sz.prod(0)) wblog(FL,
      "ERR %s() size mismatch %d/%d",FCT,(int)X.len,(int)sz.prod(0));
   if (X.len) init(sz,X.data);

   return *this;
};

#endif


template<class T>
void cell2mat(
   const wbMatrix< wbarrRef<T> > &C,
   wbarray<T> &M,
   WBINDEX &D1,
   WBINDEX &D2
){

   size_t i,j,k,m,r,s,d1=0,d2=0,i0,j0,s1,e=0;
   T *d0, *dd;
   double dc,fac;

   D1.init(C.dim1);
   D2.init(C.dim2);

   for (i=0; i<C.dim1; ++i)
   for (j=0; j<C.dim2; ++j) { if (C(i,j).D==NULL) continue;

       WBINDEX &S = C(i,j).D->SIZE;

       if (S.len==2) { d1=S[0]; d2=S[1]; } else
       if (S.len==1) { d1=d2=S[0]; }
       else wblog(FL,
         "ERR %s() requires 1/2-D arrays (got rank-%d).",FCT,S.len);

       if (C(i,j).tflag) SWAP(d1,d2);

       if (D2[j]==0) D2[j]=d2; else
       if (d2!=D2[j]) wblog(FL,
          "ERR %s() dimension mismatch (%d,%d: %d/%d).",FCT,i,j,d2,D2[j]);

       if (D1[i]==0) D1[i]=d1; else
       if (d1!=D1[i]) wblog(FL,
          "ERR %s() dimension mismatch (%d,%d: %d/%d).",FCT,i,j,d1,D1[i]);
   }

   M.init(D1.sum(),D2.sum());


   for (j0=j=0; j<C.dim2; ++j, j0+=d2) { d2=D2[j];
   for (i0=i=0; i<C.dim1; ++i, i0+=d1) { d1=D1[i]; s1=d1*sizeof(T);

      if (C(i,j).D==NULL) continue;

      d0=(C(i,j).D->data); dd=M.ref(i0,j0); 
      fac=C(i,j).fac;
      dc =C(i,j).dc;

      if (C(i,j).D->SIZE.len==1) {
         if (C(i,j).tflag) wblog(FL,
            "WRN %s() tflag irrelevant for diag-matrix (%d,%d)",
             FCT,i+1,j+1);
         if (d1!=d2) wblog(FL,"ERR d1 != d2 !?? (%,%d)",d1,d2);
      
         for (r=0; r<d1; ++r, dd+=(M.dim1+1))
         dd[r]+=(dc+fac*d0[r]);

         continue;
      }

      if (C(i,j).tflag==0) {

         if (fac==1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]+=     d0[r];
         }
         else if (fac==-1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]-=     d0[r];
         }
         else if (fac!=0.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]+= fac*d0[r];
         }
      }
      else {


         if (fac==1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]+=     d0[s+r*d2];
         }
         else if (fac==-1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]-=     d0[s+r*d2];
         }
         else if (fac!=0.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]+= fac*d0[s+r*d2];
         }
      }
      if (dc!=0.) {
         if (d1!=d2) wblog(FL,
            "WRN %s() constant diagonal shift \n"
            "requires square block (%d,%d: %d,%d)",FCT,i+1,j+1,d1,d2);

         dd=M.ref(i0,j0); m=MIN(d1,d2);

         for (r=0; r<m; ++r, dd+=(M.dim1+1))
         dd[r]+=dc;
      }
   }}
}


template<class T>
void wbarray<T>::Conj() { return; }

template<>
void wbarray<wbcomplex>::getReal(wbarray<double> &R) const {
   R.init(SIZE);
   for (size_t s=numel(), i=0; i<s; ++i)
   R.data[i]=(double)data[i].r;
}

template<class T>
void wbarray<T>::getImag(wbarray<double> &I) const {
   wblog(FL,"ERR wbarray::getImag not defined for type %s.",
   getName(typeid(T)));
}

template<class T>
void wbarray<T>::set(const wbarray<double> &R, const wbarray<double> &I) {
   wblog(FL,"ERR wbarray::set(R,I) not defined for type %s.",
   getName(typeid(T)));
}


template<>
void wbarray<wbcomplex>::getImag(wbarray<double> &I) const {
   I.init(SIZE);
   for (size_t s=numel(), i=0; i<s; ++i)
   I.data[i]=data[i].r;
}

template<>
void wbarray<wbcomplex>::Conj() {
   for (size_t s=numel(), i=0; i<s; ++i)
   data[i].i = -data[i].i;
}

template<>
void wbarray<wbcomplex>::set(
   const wbarray<double> &R,
   const wbarray<double> &I
){
   if (!R.hasSameSize(I) && !I.isEmpty()) {
      R.SIZE.print("R.SIZE"); I.SIZE.print("I.SIZE");
      wblog(FL, "ERR Dimension mismatch!");
   }

   init(R.SIZE);

   if (!I.isEmpty()) {
      for (size_t s=numel(), i=0; i<s; ++i)
      data[i].set(R.data[i], I.data[i]);
   }
   else {
      for (size_t s=numel(), i=0; i<s; ++i)
      data[i].set(R.data[i], 0.);
   }
}


#endif

