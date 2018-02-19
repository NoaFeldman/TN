#ifndef __WB_MATRIX_ROW_MAJOR_CC__
#define __WB_MATRIX_ROW_MAJOR_CC__

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
// quick integer data type // Wb,Dec02,11

template<class T> inline
bool wbMatrix<T>::thisIsInt() const { return 0; };

template <> inline
bool wbMatrix<int>::thisIsInt() const { return 1; };

template<> inline
bool wbMatrix<unsigned>::thisIsInt() const { return 1; };

template<> inline
bool wbMatrix<size_t>::thisIsInt() const { return 1; };

template <> inline
bool wbMatrix<char>::thisIsInt() const { return 1; };

template <> inline
bool wbMatrix<unsigned char>::thisIsInt() const { return 1; };


template <class T>
bool wbMatrix<T>::isNormal() const {
   wblog(FL,"ERR %s() not yet defined for type %s",FCT,
   getName(typeid(T)).data); return 0;
};


template <>
bool wbMatrix<double>::isNormal() const {
   for (size_t n=dim1*dim2, i=0; i<n; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <>
bool wbMatrix<float>::isNormal() const {
   for (size_t n=dim1*dim2, i=0; i<n; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <> bool wbMatrix<unsigned>::isNormal() const { return 1; };
template <> bool wbMatrix<size_t  >::isNormal() const { return 1; };
template <> bool wbMatrix<int     >::isNormal() const { return 1; };
template <> bool wbMatrix<char    >::isNormal() const { return 1; };
template <> bool wbMatrix<long    >::isNormal() const { return 1; };


template <class T>
template <class T0>
wbMatrix<T>& wbMatrix<T>::init_transpose(
   const T0* d0, size_t n, size_t m, char check_type
){

   size_t i=0, j=0, k=0;
   RENEW(m,n);

   if (check_type) { T0 x;
      for (i=0; i<m; ++i) {
      for (j=0; j<n; ++j, ++k) { data[k]=T(x=d0[m*j+i]);
         if (T0(data[k])!=x) wblog(FL,
            "ERR %s() got rounding error !? (%g / %g)",
            FCT, double(data[k]), double(x)
         );
      }}
   }
   else {
      for (i=0; i<m; ++i)
      for (j=0; j<n; ++j, ++k) { data[k]=T(d0[m*j+i]); }
   }

   return *this;
};


template <class T> inline
void wbMatrix<T>::init0(const mxArray *a) {
   size_t i,m,n;
   double *d;

   if (!mxIsDblMat(0,0,a))
   wberror(FL,"Need double array on input.");

   m=mxGetM(a); n=mxGetN(a); d=mxGetPr(a);
   RENEW(n,m); 

   if (typeid(T)==typeid(double))
        memcpy(data, d, m*n*sizeof(T));
   else for (n*=m, i=0; i<n; i++) data[i]=(T)d[i];
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::init(
   const char *F, int L, const mxArray *a, char tflag, char tcheck
){
   if (!a || mxIsEmpty(a)) { init(); return *this; }

   if (!mxIsDblMat(0,0,a,0,'*')) {
      wblog(F,L,"ERR numeric array required on input");
   }

   size_t m=mxGetM(a), n=mxGetN(a);

   if (tflag) {
      RENEW(n,m);
      cpyRangeMx(a, data, m*n, tcheck);
   }
   else {

      mxClassID id=mxGetClassID(a);

      if (id==mxDOUBLE_CLASS) {
         init_transpose((const double*)mxGetPr(a),n,m,tcheck);
      }
      else {
         const void *d=mxGetData(a);
         if (!d && n) wblog(FL,"ERR %s() got NULL pointer (%d)",FCT,n);

         switch (id) {
            case mxCHAR_CLASS  : init_transpose(( mxChar *)d,n,m,tcheck); break;
            case mxINT8_CLASS  : init_transpose((  int8_T*)d,n,m,tcheck); break;
            case mxUINT8_CLASS : init_transpose(( uint8_T*)d,n,m,tcheck); break;
            case mxINT16_CLASS : init_transpose(( int16_T*)d,n,m,tcheck); break;
            case mxUINT16_CLASS: init_transpose((uint16_T*)d,n,m,tcheck); break;
            case mxINT32_CLASS : init_transpose(( int32_T*)d,n,m,tcheck); break;
            case mxUINT32_CLASS: init_transpose((uint32_T*)d,n,m,tcheck); break;
            case mxINT64_CLASS : init_transpose(( int64_T*)d,n,m,tcheck); break;
            case mxUINT64_CLASS: init_transpose((uint64_T*)d,n,m,tcheck); break;
            default: wblog(FL,"ERR %s() invalid type %s",FCT,mxGetClassName(a));
         }
      }
   }

   return *this;
};


template<> inline
wbMatrix<char>& wbMatrix<char>::init(
   const char *F, int L, const mxArray *a, char tflag, char tcheck
){
   if (sizeof(mxChar)!=sizeof(short)) wblog(F_L,
      "ERR severe datatype mismatch (%d/%d)",sizeof(mxChar),sizeof(short));

   if (!a || mxIsEmpty(a)) { init(); return *this; }

   size_t m=mxGetM(a), n=mxGetN(a);
   short *d=(short*)mxGetPr(a);

   if (!mxIsChar(a) || mxGetNumberOfDimensions(a)>2) wblog(F_L,
      "ERR need CHAR array on input (%s; %d)",
      mxGetClassName(a), mxGetNumberOfElements(a)
   );

   if (tflag) {
      RENEW(n,m);
      cpyRange(d, data, m*n, tcheck);
   }
   else {
      init_transpose(d,n,m,tcheck);
   }

   return *this;
};

template<> inline
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::init(
   const char *F, int L, const mxArray *a, char tflag,
   char tcheck __attribute__ ((unused))
){
   if (!a || mxIsEmpty(a)) { init(); return *this; }

   if (mxMatSafeGuard(0,L,a,2,'c')) wblog(F_L,"ERR %s()"
      "double/complex input matrix required\n%s",FCT,str);

   size_t m=mxGetM(a), n=mxGetN(a);

   if (tflag) {
      size_t i=0, N=n*m;
      double *d=mxGetPr(a);

      RENEW(n,m);

      for (i=0; i<N; ++i) data[i].r=d[i];
      d=mxGetPi(a); if (d) {
      for (i=0; i<N; ++i) data[i].i=d[i]; }
   }
   else {
      size_t i,j,k=0;
      double *d=mxGetPr(a);

      RENEW(m,n);

      for (i=0; i<m; ++i) 
      for (j=0; j<n; ++j) data[k++].r=d[m*j+i];

      d=mxGetPi(a); if (d) { k=0;
      for (i=0; i<m; ++i) 
      for (j=0; j<n; ++j) data[k++].i=d[m*j+i]; }
   }

   return *this;
};


template <class T> inline 
void wbMatrix<T>::initI(const mxArray *a) {

   size_t i,j,m,n,k=0;
   double *d;

   if (!mxIsDblMat(0,0,a,'C'))
   wberror(FL,"Need double (complex) array on input.");

   m=mxGetM(a); n=mxGetN(a);
   RENEW(m,n);

   d=mxGetPi(a); if (d==NULL) return; 

   for (i=0; i<m; i++) 
   for (j=0; j<n; j++) data[k++]=(T)d[m*j+i];
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(size_t r, size_t c, T* d) {
   if (data && !isref) init(); isref=1;
   dim1=r; dim2=c; data=d;
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(const wbvector<T> &v, const char tflag) {
   if (data && !isref) init(); isref=1; data=v.data;
   if (tflag)
        { dim1=v.len; dim2=1; }
   else { dim1=1; dim2=v.len; }
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(const wbMatrix<T> &M) {
   if (data && !isref) init(); isref=1;
   dim1=M.dim1; dim2=M.dim2; data=M.data;
   return *this;
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::save2(wbMatrix<T> &A, char ref) {

   if (!ref) ref_check(FLF);

   if (this!=&A) {
      if (A.data) A.init();

      A.data=data; A.dim1=dim1; A.dim2=dim2; A.isref=isref;
      data=NULL; dim1=dim2=0; isref=0;
   }

   return A;
};


template <class T>
template<class TA, class TB>
wbMatrix<T>& wbMatrix<T>::init(
   const char *F, int L,
   const wbMatrix<TA> &M, const wbvector<TB> &V, char dim
){
   size_t i,j,n1=M.dim1,n2=M.dim2;
   if (abs(dim)==2) {
      if (n1!=V.len) wblog(F_L,
         "ERR %s() size mismatch %d/%dx%d",FCT,V.len,n1,n2);
      RENEW(n1,n2+1,0,0);
      if (dim>0) {
         for (i=0; i<n1; ++i) {
            for (j=0; j<n2; ++j) (*this)(i,j)=(T)M(i,j);
            (*this)(i,j)=(T)V[i];
         }
      }
      else {
         for (i=0; i<n1; ++i) {
            (*this)(i,0)=(T)V[i];
            for (j=0; j<n2; ++j) (*this)(i,j+1)=(T)M(i,j);
         }
      }
   }
   else if (abs(dim)==1) {
      if (n2!=V.len) wblog(F_L,
         "ERR %s() size mismatch %dx%d/%d",FCT,n1,n2,V.len);
      RENEW(n1+1,n2,0,0);
      if (dim>0) {
         for (i=0; i<n1; ++i)
         for (j=0; j<n2; ++j) (*this)(i+1,j)=(T)M(i,j);
         for (j=0; j<n2; ++j) (*this)(i,j)=(T)V[j];
      }
      else {
         for (j=0; j<n2; ++j) (*this)(0,j)=(T)V[j];
         for (i=0; i<n1; ++i)
         for (j=0; j<n2; ++j) (*this)(i+1,j)=(T)M(i,j);
      }
   }
   else wblog(F_L,"ERR %s() invalid dim=%d",FCT,dim);

   return *this;
};


template <class T> inline
int matchSortedIdxU(
    const char *F, int L,
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    INDEX_T m=-1, char lex=1
){
    INDEX_T m1=0, m2=0; int i=0;
    if (QA.dim2!=QB.dim2) wblog(FL,"ERR %s() dimension mismatch "
       "(%dx%d <> %d/%d)",FCT,QA.dim1,QA.dim2,QB.dim1,QB.dim2);

    try {
       i=Wb::matchSortedIdx(
          QA.data, QA.dim2, QA.dim1,
          QB.data, QB.dim2, QB.dim1, Ia,Ib, m, lex, &m1, &m2
       );
       if (m1>1 || m2>1) wblog(F_L,
         "ERR %s() got non-unique/sorted input records (%d,%d)",FCT,m1,m2);
    }
    catch (...) {
#ifdef MATLAB_MEX_FILE
       MXPut(FL,"q").add(QA,"A").add(QB,"B").add(int(m),"m").add(lex,"lex"); 
       wblog(FL,"ERR %s()",FCT);
#endif
    }

    return i;
};


template<class T> inline
bool wbMatrix<T>::quickCheckSorted(
   char dir,
   char lex,"row-major" (lex>0), else "col-major"
   const char *F, int L,
   size_t n
) const {

   char c; dir=(dir>0 ? +1 : -1);

   for (size_t i=1; i<dim1; i++) {
      c=recCompare(i,i-1,-1,lex);
      if (c) {
         if (c!=dir) {
            if (F) wblog(F,L,
              "ERR %s() input records not sorted %s (%d)",
               FCT, dir>0 ? "ascendingly":"descendingly", i);
            else return 0;
         }
         if (int(--n)<=0) break;
      }
   }

   return 1;
};


template<>
WBINDEX& wbMatrix<size_t>::toIndex(
   WBINDEX &I,
   const WBINDEX *S
 ) const {

   const size_t *idx=data;
   if (S==NULL) {
      WBINDEX SX(dim2); INDEX_T *s=SX.data; size_t i=0,j;

      for (; i<dim1; ++i, idx+=dim2) {
         for (j=0; j<dim2; ++j) {
            if (s[j]<idx[j]) s[j]=idx[j];
         }
      }
      SX+=1; return toIndex(I,&SX);
   }

   I.init(dim1); if (!dim1) return I;
   if (!dim2) wblog(FL,
      "ERR %s() got empty index matrix (%dx%d)",FCT,dim1,dim2);

   size_t i=0,j, l=dim2-1;
   const INDEX_T *s=S->data;

   for (; i<dim1; ++i, idx+=dim2) { INDEX_T &q=I.data[i];
      for (j=l, q=idx[j--]; j<l; --j) { q = q*size_t(s[j]) + size_t(idx[j]); }
   }


   return I;
};


template<>
WBIDXMAT& wbMatrix<size_t>::toIndex2D(
   const WBINDEX &S, const wbindex &ic, 
   WBIDXMAT &IJ, char pos
 ) const {

   wbvector<char> M(S.len);
   size_t i,j; char *m=M.data;

   if (S.len!=dim2) wblog(FL,
      "ERR %s() size mismatch (%dx%d/%d)",FCT,dim1,dim2,S.len);
   if (!dim2) wblog(FL,
      "ERR %s() got empty index matrix (%dx%d)",FCT,dim1,dim2);

   for (i=0; i<ic.len; ++i) {
      if (ic[i]>=S.len) wblog(FL,"ERR %s() index out of bounds "
         "(%s; %d)",FCT,(ic+1).toStr().data,S.len);
      if ((++m[ic[i]])>1) wblog(FL,"ERR %s() index not unique "
         "(%s; %d)",FCT,(ic+1).toStr().data,S.len
      );
   }

   if (pos==2) {
      for (i=0; i<M.len; ++i) { m[i] = (m[i] ? 0 : 1); }
   }
   else if (pos!=1) wblog(FL,"ERR %s() invalid pos=%d",FCT,pos);

   IJ.init(dim1,2); if (!dim1) return IJ;

   size_t k=0, l=dim2-1;
   const size_t *idx=data;
   size_t s[S.len]; for (i=0; i<S.len; ++i) s[i]=S[i];

   for (i=0; i<dim1; ++i, idx+=dim2, k+=2) {
      INDEX_T &I=IJ.data[k], &J=IJ.data[k+1]; (m[l] ? I : J) = idx[l];
      for (j=l-1; j<l; --j) {
         if (idx[j]>=S.data[j]) wblog(FL,"ERR %s() index out of bounds "
            "(%d: %d/%d)",FCT,j, idx[j],S.data[j]);
         if (m[j])
              { if (I) I*=s[j]; I+=size_t(idx[j]); }
         else { if (J) J*=s[j]; J+=size_t(idx[j]); }
      }
   }

   return IJ;
};


template<>
wbMatrix<double>& wbMatrix<double>::SkipTiny_float(double x) {
   Wb::chopTiny_float(data,dim1*dim2,x);
   return *this;
}

template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::SkipTiny_float(wbcomplex x) {
   Wb::chopTiny_float((double*)data,2*dim1*dim2,x.r);
   return *this;
}


template <class T> inline 
size_t wbMatrix<T>::maxRec(char lex, size_t *m_) const {

   if (!dim1 || !dim2) { wblog(FL,
      "WRN %s() got empty matrix (%d,%d)",FCT,dim1,dim2);
      return 0;
   }

   size_t i=1, k=0, m=1; char c;
   for (; i<dim1; i++) {
      c=recCompare(k,i,-1,lex);
      if (c<0) { k=i; m=1; }
      else if (c==0) { m++; }
   }

   if (m_) (*m_)=m;
   return k;
};


template <> inline 
size_t wbMatrix<double>::maxRec_float(
   char lex, size_t *m_, double xref
 ) const {

   wbMatrix<double> M(*this);
   Wb::chopTiny_float(M.data,dim1*dim2,xref);
   return M.maxRec(lex,m_);
};


template <class T> inline 
T wbMatrix<T>::recMax(size_t r, size_t *k) const {

   if (r>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
   if (!dim2) { wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
      return 0;
   }

   T *d=data+r*dim2, x=d[0]; if (k) (*k)=0;

   for (size_t i=1; i<dim2; ++i) {
      if (x<d[i]) { x=d[i]; if (k) (*k)=i; }
   }

   return x;
};

template <class T> inline 
T wbMatrix<T>::colMax(size_t c, size_t *k) const {

   if (c>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,c,dim2);
   if (!dim1) { if (k) (*k)=0; wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
       return 0;
   }

   T *d=data+c, x=d[0]; d+=dim2; if (k) (*k)=0;

   for (size_t i=1; i<dim1; ++i, d+=dim2) {
      if (x<d[0]) { x=d[0]; if (k) (*k)=i; }
   }

   return x;
};


template <class T> inline 
double wbMatrix<T>::recMaxA(size_t r, size_t *k) const {
   if (k) (*k)=0;

   if (r>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
   if (!dim2) { wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
       return 0;
   }

   T *d=data+r*dim2;
   double a=ABS(d[0]), amax=a;

   for (size_t i=1; i<dim2; i++) {
      a=ABS(d[i]); if (amax<a) { amax=a; if (k) (*k)=i; }
   }

   return amax;
};


template <class T>
double wbMatrix<T>::maxRelDiff(const wbMatrix<T> &M) const {

   size_t i, s=dim1*dim2;
   double maxdiff=0., maxval=0.;

   if (!hasSameSize(M)) wblog(FL,
   "ERR %s() incompatible matrix objects",FCT);

   for (i=0; i<s; i++) {
       maxdiff=MAX(maxdiff, fabs(double(M.data[i]-data[i])));
       maxval=MAX(maxval, 
              MAX( fabs(double(data[i])), fabs(double(M.data[i])) ) );
   }

   return ((maxval==0.) ? 0. : maxdiff/maxval);
};


template <class T>
T wbMatrix<T>::normDiff2(const wbMatrix<T> &M, size_t *k) const {

   if (!hasSameSize(M)) wblog(FL,
      "ERR %s() incompatible matrix objects (%dx%d, %dx%d)",
       FCT,dim1,dim2,M.dim1,M.dim2);
   T fac=1;

   return rangeNormDiff2(M.data, data, dim1*dim2, fac, k);
};

template <class T>
template <class T2>
double wbMatrix<T>::normDiff(const wbMatrix<T2> &M) const {

   if (!hasSameSize(M)) wblog(FL,
      "ERR %s() incompatible matrix objects (%dx%d, %dx%d)",
       FCT,dim1,dim2,M.dim1,M.dim2
   );

   double x, x2=0;
   for (size_t n=numel(), i=0; i<n; ++i) {
      x=(data[i]-M.data[i]); x2+=(x*x);
   }

   return sqrt(x);
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::NormalizeCol(size_t k) {

   if (k>=dim2) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
   if (!dim1) return *this;

   T *d=data+k, x2 = overlap(d,d,dim1,dim2);
   timesRange(d,1/SQRT(x2),dim1,dim2);

   return *this;
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::set2Cols(const size_t i1, const size_t i2){
   const wbMatrix<T> M; save2(M);
   M.getCols(i1,i2,*this); return *this;
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::getCols(
   size_t j1, size_t j2,
   wbMatrix<T> &M
 ) const {

   if (&M==this) {
      wbMatrix<T> X; getCols(j1,j2,X);
      return X.save2(M);
   }

   if (j1==0 && j2+1==dim2) return M.init(*this);

   if (j1>=dim2 || j2>=dim2) {
      if (int(j1-j2)==1 && (j1==0 || j1==dim2)) {
         return M.init(dim1,0);
      }
      else wblog(FL,
     "ERR index out of bounds (%dx%d: %d,%d)",dim1,dim2,j1,j2);
   }
   if (j1>j2) return M.init(dim1,0);
   
   size_t i=0, m=j2-j1+1; const T *p=data+j1;

   M.init(dim1,m);

   for (; i<dim1; ++i)
   MEM_CPY<T>(M.data+i*m, m, p+i*dim2);

   return M;
};


template <class T> inline
wbvector<T>& wbMatrix<T>::getCol(
   const size_t k, wbvector<T> &v
 ) const {

   size_t i=0; const T *p=data+k;

   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%dx%d/%d)",FCT,dim1,dim2,k);
   v.init(dim1);
   for (; i<dim1; ++i, p+=dim2) { v[i]=p[0]; }

   return v;
};


template <class T> 
template <class TI> inline
wbMatrix<T>& wbMatrix<T>::getCols(
   size_t n, const TI *I, wbMatrix<T> &M) const {

   wbMatrix<T> X(dim1,n); T *d0=data, *d=X.data;
   size_t i=0, j;

   for (j=0; j<n; ++j) if (I[j]>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d: %d/%d)",FCT,j,I[j],dim2);
   
   for (; i<dim1; ++i, d+=n, d0+=dim2) {
   for (j=0; j<n; ++j) d[j]=d0[I[j]]; }

   X.save2(M); return M;
};


template <class T> inline
wbvector<T>& wbMatrix<T>::getRec(
   const size_t j0, wbvector<T> &v
 ) const {

   if (j0>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,j0,dim1);
   v.init(dim2,rec(j0));
   return v;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::getRecs(
   const WBINDEX &I, wbMatrix<T> &M
 ) const {

   size_t i=0; M.init(I.len,dim2);
   for (; i<I.len; ++i) {
      if (I[i]>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",I[i],dim1);
      MEM_CPY<T>(M.data+i*dim2, dim2, data+I[i]*dim2);
   }
   return M;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::getRecs(
   size_t i1, size_t i2, wbMatrix<T> &M,
   size_t m
 ) const {

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR index out of bounds (%d,%d/%dx%d)",i1,i2,dim1,dim2);

    if (m==0) {
       if (i1>i2) { M.init(0,dim2); return M; }
       M.init(i2-i1+1,dim2);
       MEM_CPY<T>(M.data, M.dim1*dim2, data+i1*dim2);
    }
    else {
       bool lflag=(int(m)<0); m=abs(int(m));
       if (m>dim2) wblog(FL,
          "ERR %s() m=%d/%d out of bounds",FCT,m,dim2);

       M.init(i1<=i2 ? i2-i1+1 : 0, m);
       if (M.dim1) {
          const T* d0=data+i1*dim2; if (lflag) d0+=(dim2-m);

          for (size_t i=0; i<M.dim1; ++i, d0+=dim2) {
             MEM_CPY<T>(M.data+i*m, m, d0);
          }
       }
    }

    return M;
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlocks(
    size_t j1, size_t j2, const size_t D,
    wbMatrix<T> &M
) const {

    size_t i,isnn, D2=2*D, s=D;
    T *p2,*p0;

    if (D*(j1+1)>dim2 || D*(j2+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
    "ERR invalid block parameters (%d,%d,%d,%d).",j1, j2, D, dim2);

    M.init(dim1,D2);

    isnn=(j2==(j1+1)); 
    if (isnn) s+=s;

    p0=data+D*j1; p2=M.data;
    for (i=0; i<dim1; i++) MEM_CPY<T>(p2+i*D2, s, p0+i*dim2);

    if (!isnn) {
       p0=data+D*j2; p2=M.data+D;
       for (i=0; i<dim1; i++) MEM_CPY<T>(p2+i*D2, s, p0+i*dim2);
    }

    return M;
}


template <class T>
void wbMatrix<T>::ColKron(
   const wbMatrix<T> &B,
   int d_
){
   if (&B==this) {
      wbMatrix<T> X(B); ColKron(X,d_);
      return;
   }

   size_t i,j, da=dim1, d2=dim2, db=B.dim1;
   size_t d=( (d_<0) ? B.dim2 : (size_t)d_);

   if (d>B.dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)", FCT,d,B.dim2);

   if (B.isEmpty() || d==0) return;
   if (isEmpty()) { *this = B; return; }

   Repmat(db,1, 0,B.dim2);

   T *dd, *d0;
   for (i=0; i<db; ++i) {
      d0=B.data+i*B.dim2; dd = data + i*da*dim2 + d2;
      for (j=0; j<da; ++j, dd+=dim2) MEM_CPY<T>(dd,d,d0);
   }
};


template <class T>
wbMatrix<T>& wbMatrix<T>::repmat(
   size_t m, size_t n, wbMatrix<T> &Q,
   size_t pad1, size_t pad2
) const {

   if (this==&Q) {
      wbMatrix<T> X(*this);
      return X.repmat(m,n,Q,pad1,pad2);
   }

   size_t r,i,j;
   size_t d1=dim1+pad1, d2=dim2+pad2, D1=m*d1, D2=n*d2;

   Q.init(D1,D2);

   for (i=0; i<m; i++)
   for (j=0; j<n; j++) {
      T *d = Q.data + (i*n*d1+j)*d2;
      for (r=0; r<dim1; r++, d+=D2) MEM_CPY<T>(d, dim2, data+r*dim2);
   }

   return Q;
}


template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlock(
    size_t j1, const size_t D, wbMatrix<T> &M
) const {
    T *d, *d0=data+D*j1;

    if (D*(j1+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
       "ERR %s() invalid block parameters (%d+%d/%d)",FCT,j1,D,dim2);
    M.init(dim1,D); d=M.data;

    for (size_t i=0; i<dim1; i++, d+=D, d0+=dim2) 
    MEM_CPY<T>(d,D,d0);

    return M;
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlock(
   size_t i, size_t j,
   size_t n, size_t m,
   wbMatrix<T> &M
) const {

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d..%d/%d; %d..%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2
   );

   M.init(n,m); Wb::cpyStride(M.data,ref(i,j),m,n,-1,dim2);
   return M;
};


template <class T>
template <class T2>
wbMatrix<T>& wbMatrix<T>::getBlock(
    size_t k,
    size_t D,
    const wbMatrix<T2>& R, size_t k2, size_t D2,
    wbMatrix<T> &M
) const {

    size_t i,j;
    const T *d0; const T2 *d2; T *d;

    if (D*(k+1)>dim2 || !D || dim2%D
     || D2*(k2+1)>R.dim2 || !D2 || R.dim2%D2) wblog(FL,
       "ERR invalid block setting (%d,%d/%d; %d,%d/%d).",
        k,D,dim2, k2,D2,R.dim2);
    if (dim1!=R.dim1) wblog(FL,
       "ERR %s() dimension mismatch (%d/%d)",FCT,dim1,R.dim1);

    M.init(dim1,D+D2); d=M.data; d0=data+D*k; d2=R.data+D2*k2;

    for (i=0; i<dim1; i++) {
        MEM_CPY<T>(d,D,d0); d0+=dim2; d+=D;
        for (j=0; j<D2; j++) d[j]=T(d2[j]);
        d2+=R.dim2; d+=D2; 
    }

    return M;
};


template <class T> inline
void wbMatrix<T>::setBlock(
   size_t k, const size_t D,
   const T* d0
){ 
   T *d = data+k*D;

   if ((k+1)*D>dim2 || (D ? dim2%D : 0)) wblog(FL,
   "ERR invalid block parameters (%d,%d,%d).",k, D, dim2);

   for (size_t i=0; i<dim1; i++, d+=dim2)
   MEM_CPY<T>(d,D,d0);
};

template <class T> inline
void wbMatrix<T>::SetBlock(
   size_t i0, size_t j0,
   const wbMatrix<T> &M
){
   if (i0+M.dim1>dim1 || j0+M.dim2>dim2) wblog(FL,"ERR %s() index "
      "out of bounds (%d/%d)",FCT,i0+1,M.dim1,dim1,j0+1,M.dim2,dim2);
   const T* d0=M.data; T* d=ref(i0,j0);

   for (size_t i=0; i<M.dim1; i++, d+=dim2, d0+=M.dim2)
   MEM_CPY<T>(d,M.dim2,d0);
};


template <class T>
wbMatrix<T>& wbMatrix<T>::blockSum(size_t D, wbMatrix<T> &Q) const {

   size_t i,r,rk;
   const T *q0; T *q;

   if (&Q==this) wblog(FL,"ERR Output space same as input space!");
   if (!D || dim2%D) wblog(FL,"ERR Invalid block size (%d/%d)",D,dim2);

   rk=dim2/D; if (rk==1) { Q=(*this); return Q; }

   Q.init(dim1,D);

   for (i=0; i<dim1; i++) {
      q0=rec(i); q=Q.rec(i);
      for (r=0; r<rk; r++, q0+=D) addRange(q0,q,D);
   }

   return Q;
}

template <class T> inline 
void wbMatrix<T>::blockSum(
   size_t j1,
   size_t j2,
   const size_t D,
   wbMatrix<double> &Q
) const {

   size_t i,k;
   double *q1, *q2, *q;

   if (D*(j1+1)>dim2 || D*(j2+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
   "ERR Invalid block parameters (%d,%d,%d,%d).",j1, j2, D, dim2);

   Q.init(dim1,D);
   for (k=0; k<dim1; k++) {
       q=rec(k); q1=q+D*j1; q2=q+D*j2;
       q=Q.rec(k);
       for (i=0; i<D; i++) q[i]=q1[i]+q2[i];
   }
}

template <class T> inline 
void wbMatrix<T>::blockSum(
   const wbindex &J,
   const size_t D,
   wbMatrix<double> &Q
) const {

   size_t i,j,r=(D ? dim2/D : 0);
   double *qi, *q;

   if (D ? dim2%D : 0) wblog(FL,"ERR Invalid D=%d (%d)", D, dim2);
   for (i=0; i<J.len; i++) if (J[i]>=r) 
   wblog(FL,"ERR Block index out of range (%d; %d,%d)",J[i],D,dim2);

   Q.init(dim1,D);

   for (i=0; i<dim1; i++) {
      q=Q.rec(i);
      for (j=0; j<J.len; j++) addRange(Q.ref(i,J[j]*D), q, D);
   }
}


template <class T> inline 
void wbMatrix<T>::recSet(size_t k,
    const T* a, const wbindex &ia,
    const T* b, const wbindex &ib
){
    if (ia.len+ib.len!=dim2) wblog(FL,"ERR %s() size mismatch "
       "%dx%d/(%d: %d+%d)",FCT,dim1,dim2,k,ia.len,ib.len); 
    if (k>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%dx%d: %d)",FCT,dim1,dim2,k);

    if (dim2) { size_t i,l=0; T *d=rec(k);
       for (i=0; i<ia.len; ++i, ++l) { d[l]=a[ia[i]]; }
       for (i=0; i<ib.len; ++i, ++l) { d[l]=b[ib[i]]; }
    }
};


template <class T> inline 
void wbMatrix<T>::recSetB(
    size_t r, size_t j, size_t D,
    const T *q0, const T *q2, const T a
){
    T *q = data + r*dim2 + j*D;

    if (r>=dim1 || (j+1)*D>dim2 || (D ? dim2%D : 0)) wblog(FL,
       "ERR setRecB() index out of bounds (%d/%d; %d/%d %d)",
        r,dim1, j,dim2, D);

    if (q2==NULL)
       for (size_t i=0; i<D; i++) q[i]=q0[i];
    else if (a) {
       if (a==+1) for (size_t i=0; i<D; i++) q[i]=q0[i]+q2[i]; else
       if (a==-1) for (size_t i=0; i<D; i++) q[i]=q0[i]-q2[i];
       else       for (size_t i=0; i<D; i++) q[i]=q0[i]+a*q2[i];
    }
}


template <class T> inline 
void wbMatrix<T>::recSetB(
    size_t r, const wbindex &J, size_t D, const T *d0
){
    size_t i,j,k; const INDEX_T *jp=J.data;
    T *d = data + r*dim2;
    
    if (r>=dim1) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
    if (J.len*D!=dim2) wblog(FL,
    "ERR %s() size mismatch ([%s] %d/%d)",FCT,(J+1).toStr().data,D,dim2);

    if (D==1) {
       for (j=0; j<J.len; j++) d[j]=d0[jp[j]];
    }
    else {
       for (j=0; j<J.len; j++, d+=D) { k=jp[j]*D;
          for (i=0; i<D; i++) { d[i]=d0[k+i]; }
       }
    }
}


template <class T> inline
void wbMatrix<T>::recSet (size_t i1, size_t i2) {

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR index out of bounds (%d,%d; %d)", i1, i2, dim1);

    if (i1!=i2)
    MEM_CPY<T>(data+i1*dim2, dim2, data+i2*dim2);
}

template <class T> inline
void wbMatrix<T>::recSet (size_t i, const wbvector<T> &v) {

    if (i>=dim1) wblog(FL,
       "ERR index out of bounds (%d/%d)", i, dim1);
    if (v.len!=dim2) wblog(FL,
       "ERR incompatible objects (%d/%d)", v.len, dim2);

    MEM_CPY<T>(data+i*dim2, dim2, v.data);
}

template <class T> inline
void wbMatrix<T>::recSet (
    size_t i, const wbMatrix<T> &M, size_t j
){
    if (i>=dim1 || j>=M.dim1) wblog(FL,
       "ERR index out of bounds (%d/%d; %d/%d)", i, dim1, j, M.dim1);
    if (M.dim2!=dim2) wblog(FL,
       "ERR incompatible objects (%d/%d)", M.dim2, dim2);

    MEM_CPY<T>(data+i*dim2, dim2, M.data+j*M.dim2);
}

template <class T> inline
void wbMatrix<T>::recSet (
    size_t i, const wbvector<T> &v1, const wbvector<T> &v2
){
    if (i>=dim1) wblog(FL,
       "ERR index out of bounds (%d/%d)", i, dim1);
    if (v1.len+v2.len!=dim2) wblog(FL,
       "ERR incompatible objects (%d+%d; %d)", v1.len, v2.len, dim2);

    MEM_CPY<T>(data+i*dim2, dim2, v1.len, v1.data, v2.data);
}

template <class T> inline 
void wbMatrix<T>::recSetP(size_t i1, const T* v, size_t n) {
    if (i1>=dim1) wblog(FL,
       "ERR record index out of bounds (%d/%d)",i1,dim1);
    if (long(n)<0) n=dim2; else if (n>dim2) wblog(FL,
        "ERR %s() rec-length out of bounds (%d/%d)",FCT,n,dim2);
    T *p=data+i1*dim2; if (p!=v) MEM_CPY<T>(p,n,v);
};

template <class T> inline 
void wbMatrix<T>::recSetP(size_t i1,
    const T *v1, size_t n1,
    const T *v2, size_t n2
){
    if (i1>=dim1) {
       if (n1 || n2) wblog(FL,
          "ERR record index out of bounds (%d/%d)",i1,dim1);
       return;
    }
    if (n1+n2>dim2) wblog(FL,
       "ERR %s() size out of bounds (%d+%d/%d)",FCT,n1+n2,dim2);

    size_t i; T *p=data+i1*dim2;
    for (i=0; i<n1; ++i) { p[i]=v1[i]; }; p+=n1;
    for (i=0; i<n2; ++i) { p[i]=v2[i]; };
};

template <class T> inline 
void wbMatrix<T>::recAddP(size_t i1, const T* v) {
    T *d = data+i1*dim2;
    for (size_t i=0; i<dim2; i++) d[i]+=v[i];
}


template <class T>
template<class T2> inline
void wbMatrix<T>::setCol(size_t k, const wbvector<T2> &v) {
    if (v.len!=dim1) wblog(FL,
       "ERR %s() size mismatch (%d/%d)",FCT,v.len,dim1);
    setCol(k,v.data);
};

template <class T>
template<class T2> inline
void wbMatrix<T>::setCol(size_t k, const T2 *v, size_t stride) {

    if (k>=dim2) {
       if (dim2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
       else wblog(FL,
          "ERR %s() got empty matrix (%dx%d/%d)",FCT,dim1,dim2,k);
    }

    T *d=data+k;
    for (size_t i=0; i<dim1; i++, d+=dim2, v+=stride) (*d)=T(v[0]);
};

template <class T> inline
void wbMatrix<T>::setCol(size_t k, const T &x) {

    if (k>=dim2) {
       if (dim2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
       else wblog(FL,
          "ERR %s() got empty matrix (%dx%d/%d)",FCT,dim1,dim2,k);
    }

    T *d=data+k;
    for (size_t i=0; i<dim1; ++i, d+=dim2) (*d)=x;
};

template <class T>
template <class T2> inline
void wbMatrix<T>::setColT (size_t k, const wbvector<T2> &v) {
    if (v.len!=dim1) wblog(FL,
       "ERR %s() size mismatch (%d/%d)",FCT,v.len,dim1);
    setColT(k,v.data);
};

template <class T>
template <class T2> inline
void wbMatrix<T>::setColT (size_t k, const T2 *v, size_t stride) {

    T *d=data+k;
    if (k>=dim2) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

    for (size_t l=0, i=0; i<dim1; i++, d+=dim2, l+=stride) {
       d[0]=T(v[l]);
    }
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::scaleCol(size_t k, const T& fac) {

    size_t i=0; T *d=data+k; 
    if (k>=dim2) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

    if (fac==-1) { for (; i<dim1; ++i, d+=dim2) d[0]=-d[0]; } else
    if (fac== 0) { for (; i<dim1; ++i, d+=dim2) d[0]=0;     } else
    if (fac!=+1) { for (; i<dim1; ++i, d+=dim2) d[0]*=fac;  }

    return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::scaleRow(size_t k, const T& fac) {

    size_t i=0; T *d=data+k*dim2;
    if (k>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim1);

    if (fac==-1) { for (; i<dim2; ++i) d[i]=-d[i]; } else
    if (fac== 0) { for (; i<dim2; ++i) d[i]=0;     } else
    if (fac!=+1) { for (; i<dim2; ++i) d[i]*=fac;  }

    return *this;
};



template <class T>
size_t wbMatrix<T>::skipNanRecs(wbindex &I) {
   wblog(FL,"WRN %s() irrelevant for type <%s>",FCT,getName(typeid(T)));
   return 0;
}

template <>
size_t wbMatrix<double>::skipNanRecs(wbindex &I) {

   size_t i,j,k=0,skipped=0;
   double *d=data;
   I.init(dim1);

   for (i=0; i<dim1; i++, d+=dim2) {
      for (j=0; j<dim2; j++) if (std::isnan(d[j])) break;
      if (j<dim2) { skipped++; continue; }

      if (k<i) recSet(k,i);
      I[k++]=i;
   }

   if (k==0)
        { Resize(0,dim2); I.init(); }
   else { dim1=k; I.len=k; }

   return skipped;
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::recPermute(const wbperm &P, wbMatrix<T> &M
 ) const {

   if (P.len!=dim1 || !P.isValidPerm()) wblog(FL,
      "ERR %s() invalid permutation [%s; %d]",FCT,P.toStr().data,dim1);

   if (M.dim1!=dim1 || M.dim2!=dim2)
   M.init(dim1,dim2);

   T *d=M.data;

   for (size_t i=0; i<P.len; ++i, d+=dim2)
   MEM_CPY<T>(d, dim2, data+P[i]*dim2);

   return M;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::colPermute(
  const wbperm &P, wbMatrix<T> &M,
  char iflag
) const {

   size_t i,j; 
   const PERM_T *p=P.data;
   const T* d0; T* d2;

   if (!P.isValidPerm(0,0,dim2)) wblog(FL,
      "ERR %s() invalid permutation [%s; %d]", FCT,P.toStr().data, dim2);
   M.init(dim1,dim2); d2=M.data; d0=data;

   if (iflag) {
      for (i=0; i<dim1; i++, d0+=dim2, d2+=dim2)
      for (j=0; j<dim2; j++) d2[p[j]]=d0[j];
   }
   else {
      for (i=0; i<dim1; i++, d0+=dim2, d2+=dim2)
      for (j=0; j<dim2; j++) d2[j]=d0[p[j]];
   }

   return M;
};


template <class T> inline
void wbMatrix<T>::blockPermute(const wbperm &P, wbMatrix<T> &M
) const {

   if (this==&M) {
      wbMatrix<T> X(*this); X.blockPermute(P,M);
      return;
   }

   size_t i,j,D; const T *d0=data; T *d;
    
   if (!P.len || !P.isValidPerm()) wblog(FL,
      "ERR invalid permutation (%d;%d)", P.len, dim2);
   if (dim2%P.len) wblog(FL,
      "ERR block dimension mismatch ( %d = %d*?? )", dim2, P.len);

   M.init(dim1,dim2); d=M.data; D=dim2/P.len;

   for (i=0; i<dim1; i++, d0+=dim2)
   for (j=0; j<P.len; j++, d+=D) MEM_CPY<T>(d, D, d0+D*P[j]);
}


template <class T> inline
wbMatrix<T>& wbMatrix<T>::cols2Front(
  const WBINDEX &I1, wbMatrix<T> &M
) const {

   wbperm P; P.init2Front(I1,dim2);
   return colPermute(P,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cols2End(
  const WBINDEX &I1, wbMatrix<T> &M
) const {

   wbperm P; P.init2End(I1,dim2);
   return colPermute(P,M);
};


template <class T>
WBIDXMAT& wbMatrix<T>::toBlockIndex(
   INDEX_T D, WBIDXMAT &II_,
   wbvector< wbMatrix<T> > *QQ_
) const {

   if (!D || dim2%D)
   wblog(FL,"ERR %s() invalid blocksize (%d/%d)",dim2,D);

   size_t i,m=dim2/D;

   wbvector< wbMatrix<T> > QQ;
   wbvector< wbindex > II;
   WBINDEX dd;
   wbperm P;

   wbvector<  WBINDEX const* > Ip;

   QQ.initDef(m);
   II.initDef(m); Ip.init(m);

   for (i=0; i<m; i++) {
      getBlock(i,D,QQ[i]).groupRecs(P,dd);
      II[i].BlockIndex(dd).Permute(P,'i');
      Ip[i]=(&II[i]);
   }

   II_.CAT(2,Ip);
   if (QQ_) QQ.save2(*QQ_);

   return II_;
};



template <class T>
template <class T2>
void wbMatrix<T>::toBlockIndex(INDEX_T D,
   wbMatrix<T2> &R, INDEX_T DR_,
   WBIDXMAT *IB_,
   WBIDXMAT *I2_,
   WBIDXMAT *SS_,
   wbvector< wbMatrix<T> > *QQ_
) const {

   if (dim1!=R.dim1) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,dim1,R.dim1);
   if (!D || dim2%D) wblog(FL,
      "ERR %s() invalid blocksize (%d/%d)",FCT,dim2,D);

   if (dim1==0 || dim2==0) {
      if (IB_) IB_->init();
      if (I2_) I2_->init();
      if (SS_) SS_->init();
      if (QQ_) QQ_->init(); return;
   }

   size_t i,DR, m=dim2/D;

   wbvector< wbMatrix<T> > QQ;
   wbvector< WBINDEX > IB,I2,SB;
   WBINDEX dd;
   wbMatrix<T2> Ri;
   wbperm pp;

   if (int(DR_)<=0) DR=R.dim2/m; else DR=DR_;
   if (m*DR_!=R.dim2) wblog(FL,
      "ERR block size mismatch (expecting %d*%d = %d / %d)",
       m,DR_,m*DR_,R.dim2
   );

   QQ.init(m);
   IB.init(m); I2.init(m); SB.init(m);

   for (i=0; i<m; i++) {
      getBlock(i,D, R,i,DR, QQ[i])
     .groupRecs(pp,dd, D,
        IB[i], I2[i], SB[i],
        'p'
      );

   }

   if (IB_) IB_->CAT(2,IB);
   if (I2_) I2_->CAT(2,I2);
   if (SS_) SS_->CAT(2,SB); if (QQ_) QQ.save2(*QQ_);
};


template <class T>
void wbMatrix<T>::recPrint(
    size_t k,
    const char *istr0,
    char mflag
) const {

    size_t l=0, n=64; char s[n];
    wbvector<T> d;

    if (k>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim1);

    l=snprintf(s,n,"%.32s.rec(%d)", istr0, k);
    if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

    d.init(dim2, (*this).rec(k));
    d.print(s,mflag);

    return;
};

template <class T>
wbstring wbMatrix<T>::rec2Str(
   size_t k,
   const char *fmt,
   const char *sep,
   size_t stride,
   const char *sep2
 ) const {

   if (k>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim1);
   if (!dim2) return "";

   wbvector<T> d; d.init2ref(dim2, data+k*dim2);
   return d.toStrf(fmt,sep,stride,sep2);
};

template <class T>
wbstring wbMatrix<T>::toStr(
    const char *fmt0, const char *sep, const char *rsep,
    size_t stride, const char *sep2
) const {

    size_t i,j;
    wbstring s(MAX(size_t(128),16+dim1*dim2*16)), fmt;

    if (fmt0 && fmt0[0])
         fmt=fmt0;
    else fmt.init2Fmt((T)0);

    for (i=0; i<dim1; i++) { if (i) s.push(FL,rsep);
    for (j=0; j<dim2; j++) { if (j) s.push(FL,sep);
        s.pushf(FL,fmt.data, data[i*dim2+j]);
        if (stride && ((j+1)%stride)==0 && j+1<dim2) {
           s.push(FL,sep2);
        }
    }}

    return s;
};


template <class T>
void wbMatrix<T>::toMxStruct(mxArray* S, const char *vname) const {

    size_t i,j,k; int fid;
    double *dd;
    mxArray *a;

    fid = mxAddField2Scalar(FL,S,vname);

    a=mxCreateDoubleMatrix(dim1,dim2,mxREAL);
    dd=mxGetPr(a);

    for (k=j=0; j<dim2; j++) 
    for (  i=0; i<dim1; i++) dd[k++]=(double)data[i*dim2+j];

    mxSetFieldByNumber(S,0,fid, a);
}

template <class T>
void wbMatrix<T>::mat2mxs(mxArray* a, char fid) const {

    size_t i,n=dim1*dim2;
    double *d;

    if (!mxIsStruct(a))
    wberror(FL, "Need structure input.");

    if (mxGetNumberOfElements(a)>1)
    wberror(FL, "Need SINGLE structure on input.");

    if (fid>=mxGetNumberOfFields(a))
    wberror(FL, "Index of fields out of bounds.");

    mxSetFieldByNumber(a,0,fid, mxCreateDoubleMatrix(dim2,dim1,mxREAL));
    d=mxGetPr(mxGetFieldByNumber(a,0,fid));

    if (typeid(T)==typeid(double))
         memcpy(d, data, n*sizeof(double));
    else for (i=0; i<n; i++) d[i]=(double)data[i];
}


template <class T> inline
wbMatrix<T>& wbMatrix<T>::sortRecs(
   wbperm &P, char dir, char lex
){
   if (!dim1) { P.init(); return *this; }
   static int use_omp=-1;

#ifdef LOAD_CGC_QSPACE
   if (use_omp<0) { use_omp=0;

      int i=Wb::GetEnv(0,0,"CG_USE_OMP",use_omp);
      if (use_omp<0 || use_omp>1) wblog(FL,
         "ERR %s() invalid env CG_USE_OMP (%d; e=%d)",FCT,use_omp,i);

      if (use_omp) {
         int q=0;
         i=Wb::GetEnv(0,0,"ML_DEBUG",q);
         if (!i && q) use_omp=0;
      }

      if (use_omp && CG_VERBOSE>5)
      wblog(FL," * %s() using OMP+STL",FCT); 
   }

   char isLarge=(dim1>(1<<30));

   if (use_omp)
        sprintf(str,"(parallel mode @ %d)",omp_get_thread_num());
   else strcpy (str,"(serial mode)");

   if (isLarge && CG_VERBOSE>5) wblog(FL,
      "TST %s() got %.3fG entries @ len=%d %s",
      FCT,dim1/double(1<<30), (int)dim2, str
   );
#else
   if (use_omp<0) use_omp=0;
#endif

   if (!use_omp || dim1<128) {
      Wb::hpsort(data,dim2,dim1, P,dir,lex);
      return *this;
   }
   else {

      getSortPerm_OMP(*this,P,dir,lex);
      return recPermute(P);
   }

#ifdef LOAD_CGC_QSPACE
   if (isLarge && CG_VERBOSE>5) wblog(FL,"TST %s() done",FCT);
#endif
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::sortRecs(
   wbMatrix<T> &B, wbperm &P, char dir, char lex
){
   B=*this;
   return B.sortRecs(P,dir,lex);
};


template <class T>
template <class T2>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D, const wbMatrix<T2> &R,
   wbvector<INDEX_T> *Ib,
   wbvector<INDEX_T> *I2,
   wbvector<INDEX_T> *Sb,
   char iflag
){
   if (dim1==0) { P.init(); D.init(); return; }

   sortRecs(P);
   groupSortedRecs(D);

   if (R.dim1!=P.len) wblog(FL,
      "ERR %s() dimension mismatch (%d/%d)",FCT,R.dim1,P.len);
   if (R.dim2==0) { wblog(FL,
      "WRN %s() got empty reference space",FCT); return;
   }

   size_t i,j,d,l=0, m=D.len;
   WBINDEX dd,v, I(D.max());
   wbvector< WBINDEX > SS;
   wbvector< wbindex > II;
   wbMatrix<T2> Q;
   wbperm pp;

   if (I2 || Sb) { II.init(m); if (Sb) SS.init(m); }

   for (i=0; i<m; ++i, l+=d) { d=D[i];
       if (d>1) {
          I.len=d; for (j=0; j<d; j++) I[j]=P[l+j];
          R.getRecs(I,Q).groupRecs(pp,dd);
          for (j=0; j<d; j++) P[l+j]=I[pp[j]];

          if (I2 || Sb) {
             II[i].BlockIndex(dd); if (Sb)
             SS[i].init(pp.len).set(dd.len);
          }
       }
       else {
          if (d!=1) wblog(FL,"ERR d=%g !??",d);
          if (I2 || Sb) {
             II[i].init(1).set(0); if (Sb)
             SS[i].init(1).set(1); 
          }
          continue;
       }
   }

   if (Ib) {
      wbindex J; J.BlockIndex(D); J.save2(*Ib);
      if (Ib->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",Ib->len,P.len);
      if (iflag) Ib->Permute(P,iflag);
   }

   if (I2) {
      wbvector< WBINDEX* > J(II.len);
      for (i=0; i<II.len; i++) J[i]=(WBINDEX*)(&II[i]); I2->Cat(J);
      if (I2->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",I2->len,P.len);
      if (iflag) I2->Permute(P,iflag);
   }

   if (Sb) {
      Sb->Cat(SS);
      if (Sb->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",Sb->len,P.len);
      if (iflag) Sb->Permute(P,iflag);
   }
};


template <class T>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D, size_t nc,
   WBINDEX &Ib,
   WBINDEX &I2,
   WBINDEX &Sb,
   char iflag
){
   if (dim1==0) { P.init(); D.init(); return; }

   sortRecs(P);

   if (nc>dim2) wblog(FL,"ERR dimension out of bounds (%d/%d)",nc,dim2);
   groupSortedRecs(D,nc,Ib,I2,Sb);

   if (iflag) {
      Ib.Permute(P,iflag);
      I2.Permute(P,iflag);
      Sb.Permute(P,iflag);
   }
};


template <class T>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D,
   size_t m,
   char lex,
   wbindex *Ig
){
   if (dim1==0) { P.init(); D.init(); return; }
   if (dim2==0) {
      P.init(dim1); D.init(1); D[0]=dim1;
      if (Ig) { Ig->init(dim1).set(0); }

      if (int(m)<0) dim1=1;
      else if (m) wblog(FL,
         "ERR %s() m=%d out of bounds (%dx%d)",FCT,m,dim1,dim2);
      return;
   }

 #ifdef WB_SPARSE_CLOCK
   Wb::UseClock gr1(&wbMat_group1);
 #endif

   if (!recsSorted(+1,lex))
        sortRecs(P,+1,lex);
   else P.init(dim1);

 #ifdef WB_SPARSE_CLOCK
   gr1.done();
   Wb::UseClock gr2(&wbMat_group2);
 #endif

   if (int(m)<0) { groupSortedRecs(D,0,-1,lex); }
   else {
      if (m==0) {
         D.init(1); D[0]=dim1;
      }
      else groupSortedRecs(D,'k',m,lex);
   }

   if (Ig) {
      size_t i,j,l=0,d=0; INDEX_T *ig; PERM_T *p=P.data;
      Ig->init(P.len); ig=Ig->data;
      for (i=0; i<D.len; ++i, l+=d) { d=D.data[i];
      for (j=0; j<d; ++j) ig[p[l+j]]=i; }
   }
};

template <class T>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D,
   const WBINDEX *I, char lex
){
   if (dim1==0) { P.init(); D.init(); return; }

   sortRecs(P,+1,lex);

   if (!I) { groupSortedRecs(D,0,-1,lex); }
   else {
      wbMatrix<T> X; this->cols2Front(*I,X);
      X.groupSortedRecs(D,'k',I->len,lex);
   }
};


template <class T>
void wbMatrix<T>::groupSortedRecs(
   WBINDEX &d,
   char keepall,
   size_t m,
   char lex
){
   size_t i,ig,n=dim1;
   char c, cref=0;

   if (int(m)<0) { m=dim2; }
   else if (m>dim2) wblog(FL,"ERR number out of bounds (%d/%d)",m,dim2);
   if (!dim2) { 
      if (!dim1) wblog(FL,
         "WRN %s() for %dx%d matrix !??",FCT,dim1,dim2);
      d.init(1); d[0]=dim1; init(); return;
   }
   else if (m==0) { 
      if (m<dim2 && !keepall) wblog(FL,
         "WRN %s() keeping all since (m=%d)<%d",FCT,m,dim2);
      d.init(1); d[0]=dim1; return;
   }

   d.init(n); if (n==0) return;

   ig=0; d[ig]++;
   for (i=1; i<n; i++) {
       c=recCompare(i,i-1,m,lex);
       if (c) {
           if ((++ig)!=i) { if (!keepall)
           MEM_CPY<T>(data+ig*dim2, dim2, data+i*dim2); }

           if (c!=cref) {
               if (cref) { MXPut(FL,"a").add(*this,"M").add(m,"m");
                  wblog(FL,"ERR recs not sorted %d/%d (%d, m=%d/%d)",
                  c,cref,i,m,dim2);
               }
               else cref=c;
           }
       }

       d[ig]++;
   }

   d.Resize(ig+1);
   if (!keepall) Resize(ig+1,dim2);
};


template <class T>
void wbMatrix<T>::groupSortedRecs(
   WBINDEX &d,
   size_t m,
   char lex,
   WBINDEX &Ib,
   WBINDEX &I2,
   WBINDEX &Sb
){
   size_t i,j,i2,i0=0, ig=0;
   char c,c2, cref=0;

   if (int(m)<0) { m=dim2; }
   else if ((!m && dim2) || m>dim2) wblog(FL,
      "ERR number out of bounds (%d/%d)",m,dim2);

   if (!dim1 || !dim2) {
      d.init(); init(); Ib.init(); I2.init(); Sb.init();
      return;
   }

   Ib.init(dim1); d.init(dim1);
   I2.init(dim1);
   Sb.init(dim1);
   
   for (d[ig]++, i=1; i<=dim1; i++, d[ig]++) {
       if (i<dim1)
            c=recCompare(i,i-1,m,lex);
       else c=99;

       if (c) {
           for (i2=0, j=i0+1; j<i; j++) {
              c2=recCompare(j,j-1,-1,lex); if (c2) { i2++;
                 if (c2!=cref) {
                    if (cref) wblog(FL,
                       "ERR input recs not sorted (%d: %d/%d).",i,c2,cref);
                    else cref=c2;
                 }
              }
              I2[j]=i2;
           }
           for (i2++, j=i0; j<i; j++) { Ib[j]=ig; Sb[j]=i2; }
           i0=i; if (i>=dim1) break;

           if (c!=cref) {
               if (cref) wblog(FL,
                  "ERR input recs not sorted (%d: %d/%d).",i,c,cref);
               else cref=c;
           }

           if ((++ig)!=i)
           MEM_CPY<T>(data+ig*dim2, dim2, data+i*dim2);
       }
   }

   d.len=ig+1;
   dim1=ig+1;
};


template <> inline
wbMatrix<double>& wbMatrix<double>::sortRecs_float(wbperm &P, char dir) {

   if (dim1==0) { P.init(); return *this; }

   char lex=1;
   wbMatrix<double> X(*this); X.SkipTiny_float();

   X.sortRecs(P,dir,lex);

   return Set2Recs(P);
};


template <class T>
char wbMatrix<T>::findRecsInSet(
    const wbMatrix<T> &S, wbvector<int> &I) const {

    size_t j,l,k,d, e=0; int ns=(int)S.dim1, ifound;
    WBINDEX dg;
    wbMatrix<T> QQ(*this);
    wbperm is;

    QQ.groupRecs(is,dg); I.init(dim1);

    for (l=k=0; k<QQ.dim1; k++) {
        for (ifound=0; ifound<ns; ifound++) {
            if (QQ.recEqual(k, S.rec(ifound)))
            break;
        }
        if (ifound>=ns) { ifound=-1; e++; }

        for (d=dg[k], j=0; j<d; j++)
        I[is[l++]]=ifound;  
    }

    return e;
}


template <class T>
size_t wbMatrix<T>::findUniqueRecSorted1(size_t n, T eps) const {

   if (long(n)>=0 && (!n || n>dim2)) wblog(FL,
      "ERR %s() invalid n=%d/%d",FCT,n,dim2);
   if (dim1<=1) { return (dim1 ? 0:-1); }

   if (eps==0) {
      char c=0, c0=0;
      for (size_t l=dim1-1, i=0; i<l; ++i) {
         c=recCompare(i,i+1,n);
         if (!c || !c0) { c0=c; continue; }
         else { return i; }
      }
      if (dim1==2 && c) return 0;
   }
   else {
      T c=0, c0=0;
      for (size_t l=dim1-1, i=0; i<l; ++i) {
         c=recDiff2(i,i+1,n);
         if (c<eps || c0<eps) { c0=c; continue; }
         else { return i; }
      }
      if (dim1==2 && (c>eps)) return 0;
   }
   return -1;
};


template <class T>
int wbMatrix<T>::findRecSorted(const T* r, size_t n, char lex) const {

   if (int(n)>=0 && (!n || n>dim2)) wblog(FL,
      "ERR %s() invalid n=%d/%d",FCT,n,dim2);

   if (dim1<3) {
      for (size_t i=0; i<dim1; ++i) {
         if (!recCompareP(i,r,n,lex)) return i;
      }; return -1;
   }

   size_t k1=0, k2=dim1-1, k=(k2-k1)/2;
   char c, cref=recCompare(k1,k2,n,lex);

   if (cref==0) {
      wblog(FL,"WRN all recs the same (%dx%d) ???",dim1,dim2);
      return (recCompareP(0,r,n,lex) ? -1 : 0);
   }
   for (size_t m=MAX(size_t(1),(dim1-2)/4), i=1; i<dim1; i+=m)
   if ((c=recCompare(i-1,i,n,lex))==-cref) wblog(FL,
      "ERR records not sorted (%d/%d) !??",c,cref);

   c=recCompareP(k1,r,n,lex); if (c!= cref) return (c ? -1 : 0);
   c=recCompareP(k2,r,n,lex); if (c!=-cref) return (c ? -1 : dim1-1);
   c=recCompareP(k ,r,n,lex);

   while (1) {
      if (c==cref)
             { k1=k; k+=(k2-k1)/2; if (k==k1) return -1; } else
      if (c) { k2=k; k-=(k2-k1)/2; if (k==k2) return -1; } else
      return k;

      if (recCompare(k1,k2,n,lex)==-cref) wblog(FL,
         "ERR records not sorted !??");
      c=recCompareP(k,r,n,lex);
   }

   return -1;
};


template <class T>
int wbMatrix<T>::findRec(const T* r, size_t n, char lex) const {

   for (size_t k=0; k<dim1; k++)
   if (!recCompareP(k,r,n,lex)) return k;

   return -1;
}


template <class T>
size_t wbMatrix<T>::UnionRecs(const wbMatrix<T> &B0, wbindex &Ia){
    size_t l=0,ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B0.dim1;
    char c;

    wbMatrix<T> A(*this), B(B0);
    wbperm P1,P2; wbindex I;

    if (A.isEmpty() || B.isEmpty()) { init(); return na; }

    if (A.dim2 != B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%d/%d)",
        FCT,A.dim2,B.dim2);

    A.sortRecs(P1); I.init(na);
    B.sortRecs(P2);

    while (ia<na && ib<nb) { c=A.recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++;
       else {
          do { I[l++]=ia++; } while (ia<na && A.recEqual(ia-1,ia));
          do {        ib++; } while (ib<nb && B.recEqual(ib-1,ib));
       }
    }

    if (l) I.len=l; else I.init();

    A.Set2Recs(I); P1.get(I,Ia);

    return (na-dim1);
}


template <class T>
bool wbMatrix<T>::gotRecOverlap(const wbMatrix<T> &B0) const {

    size_t ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B0.dim1;
    char c;

    wbMatrix<T> A(*this), B(B0);
    wbperm P1,P2;

    if (A.isEmpty() || B.isEmpty()) return 0;

    if (A.dim2!=B.dim2) wblog(FL,
    "ERR %s() severe dimension mismatch (%d/%d)",FCT,A.dim2,B.dim2);

    A.sortRecs(P1);
    B.sortRecs(P2);

    while (ia<na && ib<nb) { c=A.recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++; else return 1;
    }

    return 0;
};


template <class T>
bool wbMatrix<T>::gotRecOverlapSA(const wbMatrix<T> &B) const {

    size_t ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B.dim1;
    char c;

    if (isEmpty() || B.isEmpty()) return 0;

    if (dim2!=B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%dx%d/%dx%d)",
        FCT,dim1,dim2,B.dim1,B.dim2);

    while (ia<na && ib<nb) { c=recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++; else return 1;
    }

    return 0;
};


template <class T>
size_t wbMatrix<T>::findValsCol(
    const size_t k, const wbvector<T> &v, wbindex &Ia
) const {
    wbvector<T> vk; getCol(k,vk);
    vk.findValues(v,Ia);
    return Ia.len;
}




template <class T> inline
void matchIndexU(const char *F, int L,
    const wbMatrix<T> &QA,
    const wbMatrix<T> &QB,
    wbindex &IB, const char force
){
    size_t ma,mb; wbindex Ia, Ib;

    matchIndex(QA,QB,Ia,Ib,1,&ma,&mb);
    if (mb) wblog(F,L,"ERR index B not unique (%d)",mb);

    if (force) { wbperm iP;
       if (Ia.len!=QA.dim1) wblog(F,L,
       "ERR failed to find all records (%d/%d)",Ia.len,QA.dim1);

       getIPerm(Ia,iP); Ib.select(iP,IB);
    }
    else {
       IB.init(QA.dim1).set(size_t(-1));
       IB.Set(Ia,Ib);
    }
}


template <class T>
int matchIndex(
    const wbMatrix<T> &QA, const wbMatrix<T> &QB,
    wbindex &Ia, wbindex &Ib,
    char lex,
    INDEX_T *ma, INDEX_T *mb,
    T eps
){
    size_t i;

    wbMatrix<T> Q1(QA), Q2(QB);
    wbperm P1, P2;
    wbindex ix1, ix2;

    if (QA.dim2 != QB.dim2) wblog(FL,
       "ERR %s - severe dimension mismatch (%d/%d)",
        FCT,QA.dim2, QB.dim2);

    if (QA.isEmpty() || QB.isEmpty()) {
       Ia.init(); Ib.init();
       return 0;
    }

    Q1.sortRecs(P1);
    Q2.sortRecs(P2);

    i=matchSortedIdx(Q1,Q2,ix1,ix2,-1,lex,ma,mb,eps);

    P1.get(ix1,Ia);
    P2.get(ix2,Ib);

    return (int)i;
};


template <class T>
int wbMatrix<T>::getDiff(
    const wbMatrix<T> &B, wbindex &Ia, wbindex *Ib
) const {

    wbMatrix<T> Q1(*this), Q2(B);
    wbperm P1,P2;
    wbindex ix1,ix2;

    if (dim2!=B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%d/%d)",
        FCT,dim2, B.dim2
    );

    if (isEmpty() || B.isEmpty()) {
       Ia.init(); if (Ib) Ib->init();
       return 0;
    }

    Q1.sortRecs(P1);
    Q2.sortRecs(P2);

    Q1.getDiffSorted(Q2,ix1, Ib ? &ix2 : NULL);
    P1.get(ix1, Ia); if (Ib) {
    P2.get(ix2,*Ib); }

    return Ia.len;
};


template <class T>
int wbMatrix<T>::getDiffSorted(
   const wbMatrix<T> &B, wbindex &Ia, wbindex *IB
) const {

   size_t ia,ib,la,lb; INDEX_T *Ib=NULL; char c;

   if (dim2!=B.dim2) wblog(FL,
      "ERR %s() dimension mismatch (%d/%d)",FCT,dim2,B.dim2); 

   Ia.init(dim1);
   if (IB) { IB->init(B.dim1); Ib=IB->data; }

   for (la=lb=ia=ib=0; ia<dim1 && ib<B.dim1;) {
       c=recCompareP(ia, B.rec(ib));
       if (c<0) { Ia[la++]=ia++; } else
       if (c>0) { if (Ib) Ib[lb++]=ib++; }
       else {
          while ((++ia)<dim1) {
             c=recCompare(ia-1,ia);
             if (c) {
                if (c>0) wblog(FL,
                   "ERR %s() expecting ascending order (%d)",FCT,ia);
                break;
             }
          }
          while ((++ib)<B.dim1) {
             c=B.recCompare(ib-1,ib);
             if (c) {
                if (c>0) wblog(FL,
                   "ERR %s() expecting ascending order (%d)",FCT,ib);
                break;
             }
          }
       }
   }

   while (ia<dim1) { Ia[la++]=ia++; }
   Ia.Resize(la);

   if (Ib) {
      while (ib<B.dim1) { Ib[lb++]=ib++; }
      IB->Resize(lb);
   }
   
   return Ia.len;
};



#endif

