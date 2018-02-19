#ifndef __WB_MATRIX_ROW_MAJOR_HH__
#define __WB_MATRIX_ROW_MAJOR_HH__

#ifdef WB_SPARSE_CLOCK
WbClock wbMat_group3("gr3::sortRecs");
WbClock wbMat_group2("gr2::groupRecs");
WbClock wbMat_group1("gr1::groupRecs");
#endif

/* ---------------------------------------------------------------- *
 * wbMatrix - row-major matrix class
 * AWb C Jan 2006
 * ---------------------------------------------------------------- */

template <class T>
class wbMatrix {

  public:

    wbMatrix(size_t r=0, size_t c=0)
     : data(NULL), dim1(r), dim2(c), isdiag(0), isref(0) { 
       if (r || c) NEW_DATA();
    };

    wbMatrix(size_t r, size_t c, T* d, const char ref=0)
     : data(NULL), dim1(r), dim2(c), isdiag(0), isref(0) {
       if (ref) { data=d; isref=1; }
       else if (r || c) NEW_DATA(d);
    };

    wbMatrix(const wbMatrix &M)
     : data(NULL), dim1(M.dim1), dim2(M.dim2), isdiag(M.isdiag), isref(0) {
       NEW_DATA(M.data);
    };

    template <class T2>
    wbMatrix(const wbMatrix<T2> &M)
     : data(NULL), dim1(0), dim2(0), isdiag(0), isref(0) {
       if (M.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,M.isdiag);
       initT(M);
    };

    wbMatrix(const wbvector<T> &v, const char ref, const char tflag=0)
     : data(NULL), dim1(1), dim2(v.len), isdiag(0), isref(0) {
       if (tflag) SWAP(dim1,dim2);
       if (ref) { data=v.data; isref=1; }
       else if (v.len) NEW_DATA(v.data);
    };

    wbMatrix(const mxArray *a) 
     : data(NULL), dim1(0), dim2(0), isdiag(0),isref(0) {
       init(FL,a);
    };

    wbMatrix(const char *F, int L, const mxArray *a) 
     : data(NULL), dim1(0), dim2(0), isdiag(0),isref(0) {
       init(F,L,a);
    };

   ~wbMatrix() { if (data && !isref) {
       WB_DELETE(data);
    }};

    void copyStride(T* dd, size_t stride) const {
       T *d0=data;
       for (size_t i=0; i<dim1; ++i, d0+=dim2, dd+=stride)
       MEM_CPY<T>(dd, dim2, d0);
    };

    wbMatrix& operator= (const wbMatrix &B) {
       if (this!=&B) {
          RENEW(B.dim1, B.dim2, B.data);
          isdiag=B.isdiag;
       }
       return *this;
    };

    template<class TB>
    wbMatrix& operator= (const wbMatrix<TB> &B) {
       init(B.dim1, B.dim2);
       for (size_t i=0, n=B.numel(); i<n; ++i) { data[i]=T(B.data[i]); }
       isdiag=B.isdiag;
       return *this;
    };

    bool operator==(const T &x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i)
         { if (data[i]!=x) return 0; }
       return 1;
    };

    bool operator==(const wbMatrix &B) const {
       if (this!=&B) {
          if (dim1!=B.dim1 || dim2!=B.dim2) return 0;
          if (data!=B.data) {
             for (size_t n=dim1*dim2, i=0; i<n; ++i)
             if (data[i]!=B.data[i]) return 0;
          }
       }
       return 1;
    };

    bool operator!=(const wbMatrix &B) const {
       return !((*this)==B);
    };

    bool deepEqualP(const wbMatrix &B) const;

    bool anyEqual(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
          if (data[i]==x) return 1; }
       return 0;
    };

    bool allEqual(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
          if (data[i]!=x) return 0; }
       return 1;
    };

    bool anyUnequal(const T& x) const { return !allEqual(x); };

    bool anyGT(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
           if (data[i]>x) return 1; }
       return 0;
    };

    wbMatrix& operator+= (const wbMatrix &B) {
       size_t i,s=dim1*dim2;

       if (isdiag) if (!B.isdiag) isdiag=0;

       if (dim1!=B.dim1 || dim2!=B.dim2) wblog(FL,
       "ERR Dimension mismatch (%dx%d + %dx%d).",dim1,dim2,B.dim1,B.dim2);

       for (i=0; i<s; ++i) data[i]+=B.data[i];
       return *this;
    };

    wbMatrix operator+ (const wbMatrix &B) const {
       wbMatrix Mout(*this); Mout+=B;
       return Mout;
    };

    wbMatrix operator+ (const T &x) const {
       wbMatrix Mout(*this); Mout+=x;
       return Mout;
    };

    void operator-= (const wbMatrix &B) {
       size_t i,s=dim1*dim2;

       if (isdiag) if (!B.isdiag) isdiag=0;

       if (dim1!=B.dim1 || dim2!=B.dim2) wblog(FL,
       "ERR Dimension mismatch (%dx%d + %dx%d).",dim1,dim2,B.dim1,B.dim2);

       for (i=0; i<s; ++i) data[i]-=B.data[i];
    };

    wbMatrix operator- (const wbMatrix &B) const {
       wbMatrix Mout(*this); Mout-=B;
       return Mout;
    };

    void operator*= (const T c) {
        size_t i, s=dim1*dim2;

        if (c==+1) { return; } else
        if (c== 0) { set((T)0); } else
        if (c==-1) { for (i=0; i<s; ++i) data[i]=-data[i]; }
        else       { for (i=0; i<s; ++i) data[i]*=c; }
    };

    wbMatrix& operator+= (const T c) { if (c) {
       for (size_t s=dim1*dim2, i=0; i<s; ++i) data[i]+=c; }
       return *this;
    };

    void operator-= (const T c) { if (c) {
       for (size_t s=dim1*dim2, i=0; i<s; ++i) data[i]-=c; }
    };

    wbMatrix  operator* (const T x) const {
        wbMatrix X(*this); if (x!=1) X*=x;
        return X;
    };

    wbMatrix& add  (const wbMatrix &B, const T& c=1, const char zflag=0);
    wbMatrix& minus(const wbMatrix &B, const T& c=1, const char zflag=0);
    wbMatrix& times(T fac, wbMatrix &X) const {
        X=*this; X*=fac; return X;
    };

    T prod() const {
       size_t i,n=dim1*dim2; T x=0;
        
       if (n==0) {
          wblog(FL,"WRN %s() of null object !??",FCT);
          return x;
       }
        
       for (x=data[0], i=1; i<n; ++i)
       if (x<data[i]) x*=data[i];

       return x;
    };

    T min() const {
       size_t i,n=dim1*dim2; T x=0;
        
       if (n==0) { wblog(FL,
          "WRN %s() of null object !??",FCT); return x; }
       for (x=data[0], i=1; i<n; ++i) { if (x>data[i]) x=data[i]; }

       return x;
    };

    T max() const {
       size_t i,n=dim1*dim2; T x=0;
        
       if (n==0) { wblog(FL,
          "WRN %s() of null object !??",FCT); return x; }
       for (x=data[0], i=1; i<n; ++i) { if (x<data[i]) x=data[i]; }

       return x;
    };

    T amax() const {
       size_t n=dim1*dim2; T x,m=0;
        
       if (n==0) { wblog(FL,
          "WRN %s() of null object !??",FCT); return m; }
       for (size_t i=0; i<n; ++i) { x=ABS(data[i]); if (m<x) m=x; }

       return m;
    };

    T sum() const {
        size_t i, n=dim1*dim2; T x=0;
        for (i=0; i<n; ++i) x+=data[i];
        return x;
    };

    T norm2() const {
        size_t i, n=dim1*dim2; T x=0;
        for (i=0; i<n; ++i) x+=data[i]*data[i];
        return x;
    };

    size_t maxRec(char lex=1, size_t *m=NULL) const;
    size_t maxRec_float(char lex=1, size_t *m=NULL, T xref=-1) const {
       return maxRec(lex,m);
    };

    T recMax(size_t r, size_t *k=NULL) const;
    T colMax(size_t c, size_t *k=NULL) const;

    T colMax(size_t c, unsigned &k) const { T x;
       size_t K; x=colMax(c,&K); k=K;
       if (size_t(k)!=K) wblog(FL,
          "ERR %s() unsigned out of bounds (%d/%ld)",FCT,k,K);
       return x;
    };

    wbMatrix& NormalizeCol(size_t k);

    double recMaxA(size_t r, size_t *k=NULL) const;
    
    size_t skipNanRecs(wbindex &I);

    wbMatrix& init(size_t r=0, size_t c=0, T* d=NULL) {
       RENEW(r,c,d); return *this;
    };

    wbMatrix& initDef(size_t r=0, size_t c=0, T* d=NULL) {
       RENEW(r,c,NULL,0);
       if (d) for (size_t n=r*c, i=0; i<n; ++i) data[i]=d[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(size_t r, size_t c, T2* d) {
       RENEW(r,c,0,0);
       for (size_t n=r*c, i=0; i<n; ++i) data[i]=(T)d[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(const wbMatrix<T2> &B) {
       if (B.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,B.isdiag);
       RENEW(B.dim1,B.dim2,0,0);
       for (size_t n=dim1*dim2, i=0; i<n; ++i) data[i]=(T)B.data[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(const char *F, int L, const wbMatrix<T2> &B) {
       if (B.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,B.isdiag);
       RENEW(B.dim1,B.dim2,0,0);
       for (size_t n=dim1*dim2, i=0; i<n; ++i) { data[i]=(T)B.data[i];
          if (T2(data[i])!=B.data[i]) wblog(F,L,
             "ERR type conversion changes value (%g,%g)",
              double(B.data[i]), double(T2(data[i]))
          );
       }
       return *this;
    };

    template<class TA, class TB>
    wbMatrix<T>& init(const char *F, int L,
       const wbMatrix<TA> &B, const wbvector<TB> &V, char dim=2
    );

    template<class T2>
    wbMatrix<T>& initDim(const wbMatrix<T2> &B) {
       RENEW(B.dim1, B.dim2);
       return *this;
    };

    wbMatrix& init(const wbMatrix &B) {
         RENEW(B.dim1, B.dim2, B.data); isdiag=B.isdiag;
         return *this;
    };

    template<class T0>
    wbMatrix& init_transpose(
       const T0* d0, size_t n, size_t m, char check_type=1
    );

    void init0(const mxArray *a);
    void initI(const mxArray *a);

    wbMatrix& init(
       const char *F, int L, const mxArray *a,
       char tflag=0, char check_type=1);

    void init(const mxArray *a) { init(FL,a); };

    void initTST() {           
        size_t i,j,k=0;      
        for (i=0; i<dim1; ++i) 
        for (j=0; j<dim2; ++j) data[k++]=10*(i+1)+(j+1);
    };

    wbMatrix& init2ref(size_t r, size_t c, T* d);
    wbMatrix& init2ref(const wbvector<T> &v, const char tflag=0);
    wbMatrix& init2ref(const wbMatrix &M);

    wbMatrix& unRef();
    wbMatrix& save2(wbMatrix &A, char ref=0);

    void ref_check(const char *F, int L, const char *fct) {
       if (isref) wblog(F,L,
      "ERR %s() got reference (isref=%d)",fct,isref);
    };

    wbMatrix& RESIZE_RECS(size_t r);
    wbMatrix& Resize(size_t r, size_t c, const T* dx=NULL);
    wbMatrix& resize(size_t m, size_t n, wbMatrix &M) const;

    wbMatrix& Resize2Mult(const WBINDEX &M);



    wbMatrix& transpose(wbMatrix &M) const;
    wbMatrix& Transpose() {
        wbMatrix X; save2(X); return X.transpose(*this);
    };

    wbMatrix& Reshape(size_t r, size_t c);

    template <class T2>
    wbMatrix<T>& AddCols(size_t n, const T2* d0=NULL) {
       Resize(dim1,dim2+n);
       if (d0) Wb::cpyStride(data+dim2-n, d0, n, dim1, dim2);
       return *this;
    };

    template <class T2>
    wbMatrix<T>& AddRows(size_t n, const T2* d0=NULL) {
       Resize(dim1+n,dim2);
       if (d0) MEM_CPY<T>(data+(dim1-n)*dim2, n*dim2, d0);
       return *this;
    };

    bool isEmpty() const { return data==NULL; };
    bool isNormal() const;

    bool isSquare() const { return (dim1==dim2); };
    bool isSquare(size_t d) const { return (dim1==dim2 && d==dim1); };

    bool isScalar() const { return (dim1==1 && dim2==1); };
    bool isVector() const { return (dim1==1 || dim2==1); };
    bool isDiag() const;
    bool isIdentity(const double &eps, double &maxdiff) const;

    bool isHConj(
      const wbMatrix &B, double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(B, eps, xref, 's'); };
        
    bool isHConj(
      double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(*this,eps,xref,'s'); };

    bool isAHerm(
      const wbMatrix &B, double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(B, eps, xref, 'a'); };
        
    bool isAHerm(
      double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(*this,eps,xref,'a'); };

    bool isComplex() const;

    bool isUnique() const;

    bool isUniqueSorted(char dir=0, char lex=1) const;

    char recsSorted(
       char dir=0,
       char lex=1,
       size_t m=-1
    ) const;


    bool quickCheckSorted(
        char dir,
        char lex=1,
        const char *F=NULL, int L=0,
        size_t n=2
    ) const;

    bool thisIsInt() const;

    int getDiff(
    const wbMatrix &B, wbindex &Ia, wbindex *Ib=NULL) const;

    int getDiffSorted(
    const wbMatrix &B, wbindex &Ia, wbindex *Ib=NULL) const;

    size_t length() const { return (dim1>dim2 ? dim1 : dim2); };
    size_t totSize() const { return (dim1*dim2); };
    size_t numel() const { return dim1*dim2; };

    void flipSign() {
        size_t i,s=dim1*dim2;
        for (i=0; i<s; ++i) data[i]=-data[i];
    };

    wbMatrix& Symmetrize() {
        size_t i,j,r,s;
        if (dim1!=dim2) wblog(FL,
        "ERR Symmetrize() called with %dx%d matrix !??",dim1,dim2);

        for (j=0; j<dim2; ++j)
        for (i=j+1; i<dim1; ++i) {
            r=i*dim2+j; s=j*dim2+i;
            data[r]=data[s]=0.5*(data[r]+data[s]);
        }

        return *this;
    };

    void set(const T &x) {
       size_t i,s=dim1*dim2;
       for (i=0; i<s; ++i) data[i]=x;
       isdiag=0;
    };
    wbMatrix& setRand(double fac=1., double shift=0.);
    void setDiagRand(double fac=1., double shift=0.);

    double maxRelDiff(const wbMatrix &M) const;

    template <class T2>
    double normDiff(const wbMatrix<T2> &M) const;

    T normDiff2(const wbMatrix &M, size_t *k=NULL) const;
    T normDiff (const wbMatrix &M, size_t *k=NULL) const {
       return SQRT(normDiff2(M,k));
    };

    bool hasSameSize(const wbMatrix B) const {
       return (dim1==B.dim1 && dim2==B.dim2);
    };

    WBINDEX& toIndex(WBINDEX &I, const WBINDEX *S=NULL) const {"col-major"
       wblog(FL,"ERR %s() not defined for type '%s'",FCT,
       getName(typeid(T)).data); return I;
    };

    WBIDXMAT& toIndex2D(
      const WBINDEX &S, const wbindex &ic, WBIDXMAT &IJ,
      char pos=1
    ) const {"col-major"
       wblog(FL,"ERR %s() not defined for type '%s'",FCT,
       getName(typeid(T)).data); return IJ;
    };

    wbMatrix& SkipTiny_float(T ref __attribute__ ((unused)) =-1) {
       return *this;
    };

    void reset() {
       if (dim1 && dim2) MEM_SET<T>(data,dim1*dim2);
       isdiag=0;
    }

    size_t nnz() const {
       size_t n=0,i,s=dim1*dim2;
       for (i=0; i<s; ++i) if (data[i]!=T(0)) ++n;
       return n;
    }

    wbMatrix& CAT(C_UINT d12, const wbMatrix** M, C_UINT len);
    wbMatrix& CAT(C_UINT d12, const wbvector< wbMatrix > &M0);
    wbMatrix& CAT(C_UINT d12,       wbvector< wbMatrix const* >  M);

    wbMatrix& CAT(C_UINT d12, const wbvector< wbvector<T> > &M);
    wbMatrix& CAT(C_UINT d12, const wbvector< wbvector<T> const* > &M);

    wbMatrix& cat(C_UINT d12, const wbMatrix&);
    wbMatrix& cat(C_UINT d12, const wbMatrix&, const wbMatrix&);
    wbMatrix& cat(C_UINT d12, const wbMatrix&, const wbMatrix&, const wbMatrix&);

    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&);
    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&, const wbMatrix&);
    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&, const wbMatrix&, const wbMatrix&);

    wbMatrix<T>& CAT(const unsigned dim,
        wbvector< wbMatrix<T> const* > M,
        wbvector< WBINDEX const* > I
    );

    wbMatrix& Cat(const unsigned dim,
       const wbMatrix &M1, const wbindex &I1,
       const wbMatrix &M2, const wbindex &I2
    ){
        wbvector< wbMatrix<T> const* > M(2);
        wbvector< WBINDEX const* > I(2);
        M.data[0]=&M1; I.data[0]=&I1;
        M.data[1]=&M2; I.data[1]=&I2;
        return CAT(dim,M,I);
    };

    void split(C_UINT d12, wbvector< wbMatrix* > &M) const;
    void split(C_UINT d12, wbMatrix** M, C_UINT len) const;
    void split(C_UINT d12, wbMatrix &M1, wbMatrix &M2, wbMatrix &M3) const;

    void swap(wbMatrix &B) {
        if (this!=&B) {
           if (isref || B.isref) wblog(FL,
              "ERR %s() got arrays with isref=(%d,%d)",FCT,isref,B.isref);
           SWAP(data,   B.data  );
           SWAP(dim1,   B.dim1  );
           SWAP(dim2,   B.dim2  );
           SWAP(isdiag, B.isdiag);
        }
    };

    const T& operator() (size_t i, size_t j, const T&x) const {
       if (i<dim1 && j<dim2) return data[i*dim2+j];
       else return x;
    };
    T& operator() (size_t i, size_t j, const T&x) {
       if (i<dim1 && j<dim2) return data[i*dim2+j];
       else return x;
    };

    T& operator() (size_t i, size_t j) { return data[i*dim2+j]; };
    const T& operator() (size_t i, size_t j) const {
        return data[i*dim2+j];
    };

    const T& el(size_t i, size_t j) const {
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };
    T& el(size_t i, size_t j) {
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };

    const T& operator() (size_t i) const { return data+i*dim2; };
          T& operator() (size_t i)       { return data+i*dim2; };

    const T& operator[] (size_t i) const { return data[i]; };
          T& operator[] (size_t i)       { return data[i]; };

    const T* ref(size_t i, size_t j=0) const {
       if (int(j)<0) j=dim2+j;
       if (i>=dim1 || j>dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d, %d/%d)",FCT,i,dim1,j,dim2);
       return data+(i*dim2+j);
    };
    T* ref(size_t i, size_t j=0) {
       if (int(j)<0) j=dim2+j;
       if (i>=dim1 || j>dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d, %d/%d)",FCT,i,dim1,j,dim2);
       return data+(i*dim2+j);
    };

    const T* rec(size_t i) const {
       if (i>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",i,dim1);
       return data+i*dim2;
    };
    T* rec(size_t i) {
       if (i>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",i,dim1);
       return data+i*dim2;
    };

    void getDiag(wbvector<T> &D) const {
       size_t i, n=MIN(dim1,dim2); D.init(n);
       for (i=0; i<n; ++i) D[i]=data[i*dim2+i];
    };
    wbvector<T> getDiag() const {
       wbvector<T> D; getDiag(D); return D;
    };

    void appendRows (size_t n, const T* =NULL);
    void appendRow(const wbvector<T> &v);

    mxArray* toMx (char raw=0) const;
    mxArray* toMxP(char flag=0) const;
    mxArray* toMx_base (char raw=0) const;
    mxArray* toMx_Struct() const;
    mxArray* toMx_StructP() const;
    mxArray* toMx_CellP() const;

    mxArray* toMxT(char raw=0) const { return toMx(); }

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    void toMxStruct(mxArray* S, const char *vname) const;
    void mat2mx (mxArray* &a) const { a=toMx(); };
    void mat2mx0(mxArray* &a) const;
    void mat2mxc(mxArray* C, unsigned idx) const;
    void mat2mxs(mxArray* a, char field_nr) const;

    void info(const char *istr="ans") const;
    void print(const char *istr="", char mflag=0) const;
    void Print(const char *istr="", char mflag=0) const;

    void recPrint(size_t i, const char *istr="", char mflag=0) const;

    void printdata(
    const char *istr, const char *dbl_fmt=" %8.4g", const char *rsep="\n") const;

    wbstring toStr(
       const char *fmt="", const char *sep=" ", const char *rsep="; ",
       size_t stride=0, const char *sep2=""
    ) const;

    wbstring rec2Str(
       size_t k, const char *fmt="", const char *sep=" ",
       size_t stride=0, const char *sep2=""
    ) const;

    wbstring sizeStr() const {
       size_t l=0, n=16; char s[n];
       l=snprintf(s,n,"%ldx%ld",dim1,dim2);
       if (l>=n) wblog(FL,"ERR %s() string out of bounds (%ld/%ld)",FCT,l,n);
       return s;
    };

    void put(const char *F, int L,
      const char *vname="ans", const char *ws="caller", const char raw=0
    ) const {
       mxArray *a=toMx(raw);
       mxPutAndDestroy(F_L,a,vname,ws);
    };

    void put(const char *vname="ans", const char *ws="caller",
    const char raw=0) const { put(0,0,vname,ws,raw); }

    void putg(const char *vname="ans", const char raw=0
    ) const { put(vname, "global",raw); return; };

    void putx(const char *vname="ans", const char raw=0
    ) const { put(vname, "caller", raw); return; };

    void put0(const char *vname="ans", const char *ws="caller") const {
       put(vname,ws,'r');
    }

    void getReal(wbMatrix<double> &R) const;
    void getImag(wbMatrix<double> &I) const;
    void set(const wbMatrix<double> &R, const wbMatrix<double> &I);
    wbMatrix& Conj();

    void setRec (size_t r, const T x) {
        T* d=data+r*dim2;
        if (r>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",r,dim1);
        for (size_t i=0; i<dim2; ++i) d[i]=x;
    };

    void setRec (size_t r, const T x1, const T x2) {
        size_t i=0, m=dim2%2, n=dim2-m; T* d=data+r*dim2;

        if (r>=dim1) wblog(FL,
           "ERR index out of bounds (%d/%d)",r,dim1);

        for (; i<n; i+=2) { d[i]=x1; d[i+1]=x2; }
        if (m) d[i]=x1;
    };

    void recSet (size_t idest, size_t isrc);

    void recSet (size_t, const wbvector<T>&);
    void recSet (size_t, const wbMatrix&, size_t);
    void recSet (size_t i, const wbvector<T>&, const wbvector<T>&);
    void recSetP(size_t, const T*, size_t n=-1);
    void recSetP(size_t, const T*, size_t, const T*, size_t);
    void recSetB(size_t, size_t, size_t D, const T*,
                 const T* =NULL, const T bfac=1);
    void recSetB(size_t, const wbindex&, size_t D, const T*);

    void recSet(size_t k,
       const wbvector<T> &a, const wbindex &ia,
       const wbvector<T> &b, const wbindex &ib
    ){ return recSet(a.data,ia,b.data,ib); };

    void recSet(size_t k,
       const T* a, const wbindex &ia,
       const T* b, const wbindex &ib
    );

    template<class T2>
    void setCol (size_t k, const wbvector<T2> &v);

    template<class T2>
    void setCol (size_t k, const T2 *v, size_t stride=1);

    void setCol (size_t k, const T& x);

    template <class T2>
    void setColT(size_t k, const wbvector<T2> &v);
    template <class T2>
    void setColT(size_t k, const T2* v, size_t stride=1);

    void setLastCol(const T& x) { setCol(dim2-1,x); };

    wbMatrix& scaleCol(size_t k, const T& fac);
    wbMatrix& scaleRow(size_t k, const T& fac);

    void recAddP(size_t, const T*);
    int  recLess (size_t, size_t, char lex=1) const;
    int  recLessE(size_t, size_t, char lex=1) const;
    int  recEqual(size_t, size_t) const;
    int  recEqual(size_t j1, const T *d) const;

    T recDiff2(size_t i1, size_t i2, size_t n=-1) const;
    T recDiff2P(size_t i1, const T* b, size_t n=-1) const;

    bool recAllEqual() {
       for (size_t i=1; i<dim1; ++i) if (recCompare(i,0)) return 0;
       return 1;
    };

    bool recIsZero(size_t r) {
       T *d=data+r*dim2;
       if (r>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",r,dim1);
       for (size_t i=0; i<dim2; ++i) if (d[i]!=0) return 0;
       return 1;
    };

    char recCompare(
       size_t, size_t, size_t m=-1, char lex=1, T eps=0) const;
    char recCompare(
       size_t, const wbvector<T> &r2, char lex=1) const;
    char recCompareP(
       size_t, const T*, size_t m=-1, char lex=1, T eps=0) const;

    void set2Rec(size_t i) { recSet(0,i); Resize(1,dim2); };

    wbMatrix& Set2Recs(const WBINDEX &I) {
       wbMatrix X; save2(X);
       return X.getRecs(I,*this);
    };

    wbvector<T>& getRec (size_t j0, wbvector<T> &v) const;
    wbvector<T> getRec (size_t j0) const {
       wbvector<T> v; getRec(j0,v); return v;
    };

    size_t UnionRecs(const wbMatrix &B, wbindex &Ia);

    bool gotRecOverlap  (const wbMatrix &B) const;
    bool gotRecOverlapSA(const wbMatrix &B) const;

    size_t findValsCol(
       const size_t k, const wbvector<T> &v, wbindex &Ia) const;

    wbMatrix& getRecs(const WBINDEX &, wbMatrix &M) const;
    wbMatrix  getRecs(const WBINDEX &I) const {
       wbMatrix M; getRecs(I,M);
       return M;
    };

    wbMatrix& getRecs(
       size_t j1, size_t j2, wbMatrix &M, size_t m=0) const;

    wbMatrix  getRecs(
       size_t i1, size_t i2, size_t m=0) const {
       wbMatrix M; getRecs(i1,i2,M,m);
       return M;
    };

    wbvector<T>& getCol(const size_t c, wbvector<T> &v) const;
    wbvector<T>  getCol(const size_t c) const {
       wbvector<T> v; getCol(c,v); return v;
    };

    wbMatrix getCols(const size_t j1, const size_t j2) const {
       wbMatrix M; getCols(j1,j2,M);
       return M;
    };

    wbMatrix& getCols(size_t j1, size_t j2, wbMatrix &M) const;

    template <class TI>
    wbMatrix& getCols(size_t n, const TI *I, wbMatrix<T> &M) const;

    wbMatrix& getCols(const wbperm &p, wbMatrix &M) const {
       return getCols(p.len,p.data,M);
    };

    wbMatrix& getCols(
       const wbindex &I, wbMatrix &M, char iflag=0) const {

       if (!iflag) { return getCols(I.len,I.data,M); }
       else {
          wbindex I2; I.invert(dim2,I2);
          return getCols(I2.len,I2.data,M);
       }
    };

    wbMatrix& set2Cols(
       const size_t j1, const size_t j2);

    wbMatrix& set2Cols(const wbindex &J, char iflag=0) {
       wbMatrix X; this->save2(X); X.getCols(J,*this,iflag);
       return *this;
    };

    wbMatrix& recPermute(const wbperm &P, wbMatrix &M) const;
    wbMatrix& recPermute(const wbperm &P) {
       wbMatrix M(*this);
       return M.recPermute(P,*this);
    };

    wbMatrix& colPermute(const wbperm &P, wbMatrix &M, char iflag=0) const;
    wbMatrix& colPermute(const wbperm &P, char iflag=0){
       wbMatrix M; save2(M); return M.colPermute(P,*this,iflag);
    };

    wbMatrix& FlipRecs() { wbperm P(dim1,'r'); return recPermute(P); };
    wbMatrix& FlipCols() { wbperm P(dim2,'r'); return colPermute(P); };

    wbMatrix& cols2Front(
       const WBINDEX &I1, wbMatrix &M) const;
    wbMatrix& cols2End(
       const WBINDEX &I1, wbMatrix &M) const;

    void blockPermute(const wbperm &P, wbMatrix &M) const;
    void BlockPermute(const wbperm &P){
       wbMatrix X; save2(X); X.blockPermute(P,*this);
    };

    void recProd(wbvector<T> &p) const;
    T recProd(size_t) const;

    T recSum (size_t) const;
    wbvector<T>& recSum(wbvector<T> &) const;
    wbvector<T>& recSumA(wbvector<T> &) const;
    wbvector<T> recSum() const { wbvector<T> s; return recSum(s); };
    T colSum (size_t) const;

    wbMatrix& repmat(
       size_t m, size_t n, wbMatrix &Q,
       size_t pad1=0, size_t pad2=0
    ) const;

    wbMatrix& Repmat(
       size_t m, size_t n,
       size_t pad1=0, size_t pad2=0
    ){ wbMatrix X(*this); return X.repmat(m,n,*this,pad1,pad2); };

    void ColKron(
       const wbMatrix &Q,
       int d=-1
    );

    wbMatrix& getBlocks(
       size_t j1, size_t j2,
       const size_t D, wbMatrix &M
    ) const;

    wbMatrix& getBlock(
       size_t j1,
       const size_t D, wbMatrix &M) const;

    wbMatrix getBlock(
       size_t j1,
       const size_t D
    ) const {
       wbMatrix M; getBlock(j1,D,M);
       return M;
    };

    template <class T2>
    wbMatrix<T>& getBlock(
       size_t k, size_t D,
       const wbMatrix<T2>& R, size_t k2, size_t D2,
       wbMatrix &M
    ) const;

    wbMatrix& getBlock(
       size_t i, size_t j,
       size_t n, size_t m,
       wbMatrix &M
    ) const;

    wbMatrix& blockSum(size_t D, wbMatrix&) const;
    void blockSum(size_t, size_t, const size_t, wbMatrix<double>&) const;
    void blockSum(const wbindex &J, const size_t, wbMatrix<double>&) const;

    void setBlock(size_t k, const size_t D, const T*);
    void SetBlock(size_t i0, size_t j0, const wbMatrix &M);

    WBIDXMAT& toBlockIndex(
       INDEX_T D, WBIDXMAT &II,
       wbvector< wbMatrix > *QI=NULL
    ) const;

    template <class T2>
    void toBlockIndex(INDEX_T D,
       wbMatrix<T2> &R, INDEX_T DR,
       WBIDXMAT *II=NULL, WBIDXMAT *IS=NULL,
       WBIDXMAT *SS=NULL, wbvector< wbMatrix > *QI=NULL
    ) const;

    void groupRecs(
       wbperm &P, WBINDEX &d,
       size_t m=-1,
       char lex=1,
       wbindex *Ig=NULL
    );

    void groupRecs(wbindex &Ig, INDEX_T m=-1, char lex=1) {
       wbperm P; wbvector<INDEX_T> d;
       groupRecs(P,d,m,lex,&Ig);
    };

    void groupRecs(
       groupIndex<INDEX_T> &IG, size_t m=-1, char lex=1,
       wbindex *Ig=NULL
    ){
       wbperm P; wbvector<INDEX_T> d;
       groupRecs(P,d,m,lex,Ig);
       IG.initX(P,d);
    };

    void groupRecs(
       groupIndex<INDEX_T> &IG, size_t m, char lex, const PERM_T **p
    ){
       wbperm P; wbvector<INDEX_T> d;
       groupRecs(P,d,m,lex);
       IG.initX(P,d,p);
    };

    void groupRecs(
       wbperm &P, WBINDEX &d, wbMatrix &B,
       size_t m=-1, char lex=1, wbindex *Ig=NULL
    ) const { B=*this; B.groupRecs(P,d,m,lex,Ig); };

    void groupRecs(
       wbperm &P, WBINDEX &d,
       const WBINDEX *I, char lex=1
    );

    template <class T2>
    void groupRecs(wbperm &P, WBINDEX &D, const wbMatrix<T2> &R,
       WBINDEX *Ib=NULL,
       WBINDEX *I2=NULL,
       WBINDEX *Sb=NULL,
       char iflag=0
    );

    void groupRecs(
       wbperm &P, WBINDEX &D, size_t nc,
       WBINDEX &Ib,
       WBINDEX &I2,
       WBINDEX &Sb,
       char iflag=0
    );

    void groupSortedRecs(
       WBINDEX &d, char keepall=0, size_t m=-1, char lex=1);

    void groupSortedRecs(
       WBINDEX &d, const WBINDEX *I);

    void groupSortedRecs(
       WBINDEX &d,
       size_t mc,
       char lex,
       WBINDEX &Ib,
       WBINDEX &I2,
       WBINDEX &Sb
    );

    void makeUnique(){ if (dim1<2) return;
       wbperm P; WBINDEX D;
       groupRecs(P,D);
    };
    void makeUnique(wbperm &P, WBINDEX &D) {
       if (dim1<2) {
          P.init(dim1); D.init(dim1); if (dim1) D[0]=1;
          return;
       }
       groupRecs(P,D);
    };

    void makeUnique1(wbindex &I) {
       if (dim1<2) { I.init(dim1); return; }
       wbperm P; WBINDEX D; groupRecs(P,D); I.init(D.len);
       for (size_t d,l=0,i=0; i<D.len; ++i, l+=d) { d=D[i]; I[i]=P[l]; }
    };

    void findUnique(wbindex &I) const {
       if (dim1<2) { I.init(dim1); return; }

       wbMatrix X(*this);
       wbperm P; WBINDEX D; X.groupRecs(P,D); 
       I.init(D.numelFind(T(1)));
       for (size_t d, l=0, j=0, i=0; i<D.len; ++i, l+=d) {
          if ((d=D[i])==1) I[j++]=P[l];
       }
    };

    size_t findUniqueRecSorted1(size_t n=-1, T eps=1E-12) const;

    void makeUnique_ig(wbindex &Iu){
       wbperm P; wbvector<INDEX_T> D; groupRecs(P,D);
       Iu.initGroup(P,D);
    };
    void makeUnique(groupIndex<INDEX_T> &IG){ groupRecs(IG); };

    wbMatrix& sortRecs(wbperm &P,
       char dir=+1,
       char lex=1
    );

    wbMatrix& sortRecs(
       wbMatrix &B, wbperm &P,
       char dir=+1, char lex=1
    );

    wbMatrix& sortRecs(char dir=+1, char lex=1) {
       wbperm P; return sortRecs(P,dir,lex);
    };

    wbMatrix& sortRecs_float(wbperm &P, char dir=+1) {
       P.init(dim1); return *this;
    };
    wbMatrix& sortRecs_float() { wbperm P; return sortRecs_float(P); };

    char findRecsInSet(const wbMatrix &S, wbvector<int> &I) const;
    int findRecSorted(const T* r, size_t n=-1, char lex=1) const;
    int findRec(const T* r, size_t n=-1, char lex=1) const;

    T *data;

    size_t dim1, dim2;

    char isdiag, isref;

  protected:

    void NEW_DATA(T* d=NULL) { size_t s=dim1*dim2;

        if (data) { WB_DELETE(data); }
        if (s) {
            WB_NEW(data,s);
            MEM_CPY<T>(data,s,d);
        }
        else data=NULL;
    };

    void RENEW(
        const size_t &d1, const size_t &d2, T* dd=NULL, char iflag=1
    ){
        const size_t s=d1*d2; isdiag=0;
        
        if (isref) {
           if (d1==0 || d2==0) {
              if (d1 || d2) wblog(FL,
                 "WRN wbMatrix is declared as reference "
                 "(%dx%d => %dx%d; %d)",dim1,dim2,d1,d2,isref
              );

              dim1=d1; dim2=d2; isref=0;

              data=NULL;

              return;
           }
           else {
              wblog(FL,"ERR wbMatrix is declared as reference (%dx%d; %d)",
              d1,d2, isref);
           }
        }

        if (d1!=dim1 || d2!=dim2) {

            dim1=d1; dim2=d2; if (data) {
                if (dd==data) wblog(FL,
                   "ERR init space equals *this (use Resize instead)");
                WB_DELETE(data);
            }
            if (d1==0 || d2==0) return;
            WB_NEW(data,s);
        }

        if (!iflag && !dd) return;
        if (s) MEM_CPY<T>(data,s,dd);
    };

  private:

    bool isSym_aux(
       const wbMatrix &B, double eps, double *xref,
       const char symflag='s'
    ) const;
};


template <class T>
class wbRecs {

    const T* data;
    size_t n,lda; bool asc, lex;

  public:

    wbRecs(const wbMatrix<T> &A, wbperm &P, char dir_=1, char lex_=1)
     : data(A.data), n(A.dim1), lda(A.dim2),
       asc(dir_> 0 ? 1 : 0), lex(lex_> 0 ? 1 : 0)
    {
       P.init(n);
    };

   ~wbRecs() {};

    bool operator()(const size_t &pa, const size_t &pb) const {

       if (pa>=n || pb>=n) wblog(FL,
          "ERR %s() index out of bounds (%d,%d; %d)",FCT,pa,pb,n);

       const T *a=data+pa*lda, *b=data+pb*lda;
       if (lex) { size_t i=0;
          if (asc) {
             for (; i<lda; ++i) {
                 if (a[i]<b[i]) return 1;
                 if (a[i]>b[i]) return 0;
             }
          }
          else {
             for (; i<lda; ++i) {
                 if (a[i]<b[i]) return 0;
                 if (a[i]>b[i]) return 1;
             }
          }
       }
       else { size_t i=lda-1;
          if (asc) {
             for (; i<lda; --i) {
                 if (a[i]<b[i]) return 1;
                 if (a[i]>b[i]) return 0;
             }
          }
          else {
             for (; i<lda; --i) {
                 if (a[i]<b[i]) return 0;
                 if (a[i]>b[i]) return 1;
             }
          }
       }
       if (asc) return (pa<pb);
       else return (pb<pa);
    };
};


template <class T>
void getSortPerm_OMP(
   const wbMatrix<T> &A, wbperm &P, char dir, char lex
){
   if (!A.dim2) { P.init(); return; }

   wbRecs<T> R(A,P,dir,lex);

 #ifdef WB_SPARSE_CLOCK
   Wb::UseClock gr3(&wbMat_group3);
 #endif

   if (A.dim1<128) {
      sort(P.data, P.data+A.dim1, R);
      return;
   }

   int id=0, i,l, im=0, mp=MIN(
      1 << unsigned(floor(log2(double(A.dim1)))-6),
      omp_get_max_threads()
   );

   double nsub=double(A.dim1)/mp;

   wbindex idx(mp+1); wbperm PX(P.len);
   for (i=0; i<mp; ++i) { idx[i]=size_t(i*nsub+0.5); }
   idx[i]=A.dim1;

#pragma omp parallel for shared (id)
   for (i=0; i<mp; ++i) {
      size_t i1=idx.data[i], i2=idx.data[i+1];
      sort(P.data+i1, P.data+i2, R);
      id=MAX(id,omp_get_thread_num());
   }
#pragma omp barrier

   while (idx.len>2) { int m2=mp/2; P.swap(PX); ++im;

#pragma omp parallel for shared (id)
      for (i=0; i<m2; ++i) {"i+2<=mp"
         size_t k=2*i, i1=idx.data[k], i2=idx.data[k+1], i3=idx.data[k+2];
         merge(
            PX.data+i1, PX.data+i2,
            PX.data+i2, PX.data+i3, P.data+i1, R
         );
         id=MAX(id,omp_get_thread_num());
      }
#pragma omp barrier

      if (mp%2) {
         size_t i1=idx.data[mp-3], i2=idx.data[mp-1], i3=idx.data[mp];
         memcpy(PX.data+i1,P.data+i1,(i2-i1)*sizeof(size_t));
         merge(
            PX.data+i1, PX.data+i2,
            PX.data+i2, PX.data+i3, P.data+i1, R
         );
         idx[mp-1]=idx[mp]; idx.len=(mp--);
      }

      for (l=1, i=2; i<mp; i+=2, ++l) { idx[l]=idx[i]; }
      idx[l]=idx[mp]; mp=l; idx.len=l+1;
   }

};


template <class T>
bool wbMatrix<T>::deepEqualP(const wbMatrix<T> &B) const {
    if (dim1!=B.dim1 || dim2!=B.dim2) return 0;

    if (data!=B.data) {
       size_t i=0, n=numel();
       for (; i<n; ++i) { if (data[i]==B.data[i]) continue;
          if ((!data[i] || !B.data[i])) return 0; else
          if (!((*data[i])==(*B.data[i]))) return 0;
       }
    }

    return 1;
};


template <class T> inline
bool wbMatrix<T>::isDiag() const {
   size_t i,j;

   for (i=0; i<dim1; ++i)
   for (j=0; j<dim2; ++j) if (i!=j) if (data[i*dim2+j]!=0.) return 0;

   return 1;
}

template <class T>
bool wbMatrix<T>::isIdentity(const double &eps, double &maxdiff) const {

   if (dim1!=dim2) return 0;

   size_t i,j,l=0; double d; maxdiff=0;

   for (j=0; j<dim2; ++j)
   for (i=0; i<dim1; ++i, ++l) {
       d = ABS( data[l] - (i!=j ? 0 : 1) );
       if (maxdiff<d) maxdiff=d;
   }

   return (maxdiff<eps);
}


template <class T> inline
bool wbMatrix<T>::isSym_aux(
   const wbMatrix<T> &B, double eps, double *xref,
   const char symflag
) const {

   const wbMatrix<T> &A=(*this);
   char issame = (&A==&B);
   size_t i,j,k;
   double x;

   if (A.dim1!=B.dim2 || A.dim2!=B.dim1) return 0;
   if (xref) *xref=0;


   if (symflag=='s') {
      if (eps==0) {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            if (A(i,j)!=CONJ(B(j,i))) {
               if (xref) {
                  x = ABS(A(i,j) - CONJ(B(j,i)));
                  xref[0] = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            x = ABS(A(i,j) - CONJ(B(j,i)));
            if (x>eps) { if (xref) xref[0]=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else if (symflag=='a') {
      if (eps==0) {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            if (A(i,j)!=-CONJ(B(j,i))) {
               if (xref) {
                  x = ABS(A(i,j) + CONJ(B(j,i)));
                  xref[0] = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            x = ABS(A(i,j) + CONJ(B(j,i)));
            if (x>eps) { if (xref) xref[0]=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else wblog(FL,"ERR invalid symflag=%c<%d>",symflag,symflag);

   if (xref && (*xref)>eps) return 0;

   return 1;
}


template <class T> inline
bool wbMatrix<T>::isComplex() const { return 0; };

template<> inline
bool wbMatrix<wbcomplex>::isComplex() const {
   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   if (data[i].i!=0.) return 1;

   return 0;
}


template <class T> inline
bool wbMatrix<T>::isUnique() const {

   wbMatrix<T> Q(*this);
   wbperm P;

   if (dim2==0) return (dim1<=1);

   Q.sortRecs(P);
   return Q.isUniqueSorted();
};


template <class T> inline
bool wbMatrix<T>::isUniqueSorted(char dir, char lex) const {

   if (dim1<=1 || !dim2) return 1;
   size_t i=1; char c=0;

   if (!dir) { if (dim1<=2) return 1;
      c=recCompare(1,0,-1,lex); if (!c) { return 0; }
      i=2;
   }
   else if (abs(int(dir))>1) {
      if (dir=='A') c=+1; else
      if (dir=='D') c=-1; else
      wblog(FL,"ERR %s() invalid c='%c'<%d> !??",FCT,c,c); 
   }
   else { c=dir; }

   for (; i<dim1; ++i) {
      if (recCompare(i,i-1,-1,lex)!=c) return 0;
   }

   return 1;
};


template <class T> inline
char wbMatrix<T>::recsSorted(char dir, char lex, size_t m) const {

   if (dim1<2 || !dim2) return 1;
   size_t j=1, ndeg=0; char c;

   if (!dir) {
      for (; j<dim1; ++j) {
         c=recCompare(j,j-1,m,lex);
         if (c) { dir=c; ++j; break; } else ++ndeg;
      }
   }
   else if (abs(int(dir))>1) {
      wblog(FL,"ERR %s() invalid direction s=%d",FCT,dir);
   }

   for (; j<dim1; ++j) {
      c=recCompare(j,j-1,m,lex);
      if (c) { if (c!=dir) return 0; } else ++ndeg;
   }

   return (ndeg ? 1 : 2);
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::RESIZE_RECS(size_t nr) {
    
    ref_check(FLF);

    if (nr) {
       if (nr<=dim1 && nr/double(dim1)>0.10) dim1=nr;
       else Resize(nr,dim2);
    }
    else init(0,dim2);

    return *this;
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::resize(
   size_t m, size_t n, wbMatrix<T> &M
 ) const {

   if (this==&M) return M.Resize(m,n);

   M.init(m,n);
   m=MIN(m,dim1); n=MIN(n,dim2); if (m==0 || n==0) return M;

   for (size_t i=0; i<m; ++i)
   MEM_CPY<T>(M.data+i*M.dim2, n, data+i*dim2);

   return M;
};

template <class T> inline 
wbMatrix<T>& wbMatrix<T>::Resize(size_t nr, size_t nc, const T* dx) {
   
   size_t d1,d2,s=nr*nc; T *d0=data;
   ref_check(FLF);

   d1=MIN(dim1,nr); d2=MIN(dim2,nc); 

   if (dim1==nr && dim2==nc) return *this;
   if (d1==0 || d2==0) { RENEW(nr,nc); return *this; }

   WB_NEW(data,s);

   if (nc==dim2) {
      MEM_CPY<T>(data,d1*dim2,d0);
      if (nr>dim1) {
         if (dx)
              MEM_CPY<T>(data+dim1*dim2,(nr-dim1)*dim2,dx);
         else MEM_SET<T>(data+dim1*dim2,(nr-dim1)*dim2);
      }
   }
   else {
      size_t i; T *p=data, *p0=d0;
      for (i=0; i<d1; ++i, p+=nc, p0+=dim2) {
         MEM_CPY<T>(p,d2,p0); if (d2<nc) {
         MEM_SET<T>(p+d2,nc-d2); }
      }
      if (nr>dim1) {
         if (dx) wblog(FL,"ERR %s() "
            "allows no reference data if dim2=%d->%d",FCT,dim2,nc);
         MEM_SET<T>(data+dim1*dim2,(nr-dim1)*dim2);
      }
   }

   dim1=nr; dim2=nc;
   WB_DELETE(d0);

   return *this;
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::Resize2Mult(const WBINDEX &M) {

   size_t m=M.sum();

   if (M.len!=dim1) wblog(FL,
      "ERR %s() size mismatch (%d/%dx%d)",FCT,M.len,dim1,dim2);
   if (!m || !dim2) { return init(m,dim2); } else
   if (int(m)<0) wblog(FL,"ERR %s() got m=%d !??",FCT,m);

   size_t i,j;
   wbMatrix<T> X(m,dim2);
   const T *x0=data; T *x=X.data;

   for (i=0; i<M.len; ++i, x0+=dim2) { m=M.data[i];
   for (j=0; j<m; ++j, x+=dim2) MEM_CPY<T>(x, dim2, x0); }

   return X.save2(*this);
};
   

template <class T> inline
wbMatrix<T>& wbMatrix<T>::unRef() {

    if (!isref) return *this;
    if (!dim1 || !dim2) {
       if (data) wblog(FL,"ERR data=%lX (%dx%d) !??",data,dim1,dim2);
       isref=0; return *this;
    }
    
    size_t s=dim1*dim2;
    T const* const d0=data;
    
    WB_NEW(data,s);

    MEM_CPY<T>(data,s,d0); 
    isref=0;

    return *this;
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::transpose(wbMatrix<T> &M) const {

    if (this==&M) return M.Transpose();

    size_t i,j,k=0; M.init(dim2,dim1);

    if (isdiag) {
       M.isdiag=1; k=MIN(dim1,dim2);
       for (i=0; i<k; ++i) M.data[i*dim1+i]=data[i*dim2+i];
       return M;
    }

    for (j=0; j<dim2; ++j)
    for (i=0; i<dim1; ++i) M.data[k++]=data[i*dim2+j]; 

    return M;
};


template <class T> inline 
wbMatrix<T>& wbMatrix<T>::Reshape(size_t nr, size_t nc) {
    
    if (dim1!=nr || dim2!=nc) {
       if (dim1*dim2!=nr*nc) wblog(FL,
          "ERR %s() size not preserved (%dx%d => %dx%d)",
           FCT,dim1,dim2,nr,nc
       );
       dim1=nr; dim2=nc;
    }
    return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::add(
    const wbMatrix<T> &M,
    const T& c,
    const char zflag
){
    size_t i,s=dim1*dim2;

    if (data==NULL && zflag) {
       RENEW(M.dim1, M.dim2, M.data); (*this)*=c;
       isdiag=M.isdiag;
       return *this;
    }

    if (dim1!=M.dim1 || dim2!=M.dim2)
    wberror(FL,"Dimension mismatch.");

    if (isdiag) if (!M.isdiag) isdiag=0;

    if (c==+1) { for (i=0; i<s; ++i) data[i]+=M.data[i]; } else
    if (c== 0) { return *this; } else
    if (c==-1) { for (i=0; i<s; ++i) data[i]-=M.data[i]; }
    else       { for (i=0; i<s; ++i) data[i]+=(c * M.data[i]); }

    return *this;
}

template <class T> inline
wbMatrix<T>& wbMatrix<T>::minus(
    const wbMatrix<T> &M,
    const T& c,
    const char zflag
){
    size_t i,s=dim1*dim2;

    if (data==NULL && zflag) {
       RENEW(M.dim1, M.dim2, M.data); (*this)*=(-c);
       isdiag=M.isdiag;
       return *this;
    }

    if (dim1!=M.dim1 || dim2!=M.dim2)
    wberror(FL,"Dimension mismatch.");

    if (isdiag) if (!M.isdiag) isdiag=0;

    if (c==+1) { for (i=0; i<s; ++i) data[i]-=M.data[i]; } else
    if (c== 0) { return *this; } else
    if (c==-1) { for (i=0; i<s; ++i) data[i]+=M.data[i]; }
    else       { for (i=0; i<s; ++i) data[i]-=(c * M.data[i]); }

    return *this;
}


template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2
){
    const wbMatrix<T>* M0[] = {this, &M2};
    wbvector< wbMatrix<T> const* > M; M.init(2,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2, const wbMatrix<T> &M3
){
    const wbMatrix<T>* M0[] = {this, &M2, &M3};
    wbvector< wbMatrix<T> const* > M; M.init(3,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2, const wbMatrix<T> &M3, const wbMatrix<T> &M4
){
    const wbMatrix<T>* M0[] = {this, &M2, &M3, &M4};
    wbvector< wbMatrix<T> const* > M; M.init(4,M0);
    return CAT(dim,M);
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2
){
    wbMatrix<T> const* M0[] = {&M1, &M2};
    wbvector< wbMatrix<T> const* > M; M.init(2,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2,
    const wbMatrix<T> &M3
){
    wbMatrix<T> const* M0[] = {&M1, &M2, &M3};
    wbvector< wbMatrix<T> const* > M; M.init(3,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2,
    const wbMatrix<T> &M3, const wbMatrix<T> &M4
){
    wbMatrix<T> const* M0[] = {&M1, &M2, &M3, &M4};
    wbvector< wbMatrix<T> const* > M; M.init(4,M0);
    return CAT(dim,M);
};


template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    wbMatrix<T> const* M[], const unsigned len
){
    wbvector< wbMatrix<T> const* > MM;
    MM.init(len,M);

    return CAT(dim,MM);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbMatrix<T> > &M0
){
    wbvector< wbMatrix<T> const* > M(M0.len);
    for (unsigned i=0; i<M0.len; ++i) M[i]=&M0[i];

    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbvector<T> > &V0
){
    wbvector< wbMatrix<T> > M0(V0.len);
    wbvector< wbMatrix<T> const* > M(V0.len);

    for (size_t i=0; i<V0.len; ++i) {
       if (dim==1)
            M0[i].init2ref(1,V0[i].len,V0[i].data);
       else M0[i].init2ref(V0[i].len,1,V0[i].data);

       M[i]=&M0[i];
    }

    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbvector<T> const* > &V0
){
    wbvector< wbMatrix<T> > M0(V0.len);
    wbvector< wbMatrix<T> const* > M(V0.len);

    for (size_t i=0; i<V0.len; ++i) {
       if (dim==1)
            M0[i].init2ref(1,V0[i]->len,V0[i]->data);
       else M0[i].init2ref(V0[i]->len,1,V0[i]->data);

       M[i]=&M0[i];
    }

    return CAT(dim,M);
};


template <class T>
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    wbvector< wbMatrix<T> const* > M
){
    wbMatrix<T> X;
    size_t i,k;

    if (isdiag) { wblog(FL,"WRN unsetting isdiag"); isdiag=0; }

    if (M.len==0) { this->init(); return *this; }

    for (k=i=0; i<M.len; ++i) {
       if (M[i]==this) {
          if (!(k++)) { this->save2(X); }
          M[i]=&X;
       }
    }

    for (k=i=0; i<M.len; ++i) {
       if (M[i] && (M[i]->dim1 || M[i]->dim2)) {
          if (k<i) { M[k]=M[i]; }; ++k;
          if (M[i]->isdiag) wblog(FL,
             "ERR %s() got diagonal matrix (i=%d)",FCT, i+1
          );
       }
    }
    if (k<M.len) {
       if (k) M.len=k;
       else { this->init(); return *this; }
    }

    if (dim==1) {
       size_t m,d1,D1,d2=0;

       for (D1=k=0; k<M.len; ++k) { D1+=(M[k]->dim1);
          if (k==0) d2=(M[k]->dim2); else
          if (d2!=(M[k]->dim2)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d2,M[k]->dim2
          );
       }

       this->init(D1,d2);

       for (m=k=0; k<M.len; ++k) { d1=(M[k]->dim1);
          MEM_CPY<T>(data+d2*m, d1*d2, M[k]->data);
          m+=d1;
       }
    }
    else if (dim==2) {
       size_t i,m,d1=0,d2,D2;
       T *p, *p2;

       for (D2=k=0; k<M.len; ++k) { D2+=(M[k]->dim2);
          if (k==0) d1=(M[k]->dim1); else
          if (d1!=(M[k]->dim1)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d1,M[k]->dim1
          );
       }

       this->init(d1,D2);

       for (m=k=0; k<M.len; ++k) { d2=(M[k]->dim2);
          p=data+m; p2=(M[k]->data);

          for (i=0; i<d1; ++i)
          MEM_CPY<T>(p+D2*i, d2, p2+d2*i);

          m+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat dim=%d",dim);

    return *this;
};


template <class T>
wbMatrix<T>& wbMatrix<T>::CAT(
    const unsigned dim,
    wbvector< wbMatrix<T> const* > M,
    wbvector< WBINDEX const* > I
){
    wbMatrix<T> X;
    size_t i,k;

    if (M.len!=I.len) wblog(FL,
       "ERR %s() size mismatch %d/%d",FCT,M.len,I.len);
    if (!M.len) { this->init(); return *this; }

    for (k=i=0; i<M.len; ++i) {
       if (M[i]==this) {
          if (!(k++)) { this->save2(X); }
          M[i]=&X;
       }
    }

    for (k=i=0; i<M.len; ++i) {
       if (M[i] && (M[i]->dim1 || M[i]->dim2 || I[i]->len)) {
          if (k<i) { M[k]=M[i]; I[k]=I[i]; }; ++k;
          if (M[i]->isdiag) wblog(FL,
             "ERR %s() got diagonal matrix (i=%d)",FCT, i+1
          );
       }
    }
    if (k<M.len) {
       if (k) { M.len=k; I.len=k; }
       else { this->init(); return *this; }
    }

    if (isdiag) isdiag=0;

    if (dim==1) {
       size_t d1,d2=0, D1=0;

       for (k=0; k<M.len; ++k) { D1+=(I[k]->len);
          if (k==0) d2=(M[k]->dim2); else
          if (d2!=(M[k]->dim2)) wblog(FL,
             "ERR dimension mismatch (vcat: dim2=%d/%d)",d2,M[k]->dim2
          );
       }

       this->init(D1,d2);
       T *x=data;

       for (k=0; k<M.len; ++k) {
          const T* x0=M[k]->data;
          const INDEX_T *Ik=I[k]->data; D1=M[k]->dim1; d1=I[k]->len;

          for (i=0; i<d1; ++i, x+=d2) {
              if (Ik[i]>=D1) wblog(FL,
                 "ERR %s() index out of bounds (%d/%d)",FCT,Ik[i]+1,D1);
              MEM_CPY<T>(x, d2, x0+Ik[i]*d2);
          }
       }
    }
    else if (dim==2) {
       size_t i,d1=0,d2,D1,D2=0;

       for (k=0; k<M.len; ++k) { D2+=(d2=M[k]->dim2);
          if (k==0) d1=(I[k]->len); else
          if (d1!=(I[k]->len)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d1,I[k]->len
          );
       }

       this->init(d1,D2);
       T *x=data;

       for (k=0; k<M.len; ++k) {
          const T* x0=M[k]->data;
          const INDEX_T *Ik=I[k]->data; D1=M[k]->dim1; d2=(M[k]->dim2);

          for (i=0; i<d1; ++i) {
             if (Ik[i]>=D1) wblog(FL,
                "ERR %s() index out of bounds (%d/%d)",FCT,Ik[i]+1,D1);
             MEM_CPY<T>(x+i*D2, d2, x0+Ik[i]*d2);
          }; x+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat dim=%d",dim);

    return *this;
};


template <class T> inline
void wbMatrix<T>::split(
    const unsigned dim,
    wbvector< wbMatrix<T>* > &M) const {
    
    split(dim,M.data,M.len);
    return;
}


template <class T> inline
void wbMatrix<T>::split(
    const unsigned dim,
    wbMatrix<T> &M1,
    wbMatrix<T> &M2,
    wbMatrix<T> &M3
) const {

    wbMatrix<T>* M[] = {&M1,&M2,&M3};
    split(dim,M,3);
    return;
}

template <class T>
void wbMatrix<T>::split(
    const unsigned dim, wbMatrix<T>** M, const unsigned len
)const {

    unsigned k;

    for (k=0; k<len; ++k) {
       if (M[k]==this)
       wberror(FL, "This matrix itself is included in split list!");
    }

    if (dim==1) {
       size_t m,d1,D1;

       for (D1=k=0; k<len; ++k) D1+=(*(M[k])).dim1;

       if (D1!=dim1)
       wberror(FL, "Dimension mismatch.");

       for (k=0; k<len; ++k) {
          if (M[k]->dim2!=dim2) M[k]->init(M[k]->dim1, dim2);
       }

       for (m=k=0; k<len; ++k) { d1=M[k]->dim1;
          MEM_CPY<T>(M[k]->data, d1*dim2, data+dim2*m);
          m+=d1;
       }
    }
    else if (dim==2) {
       size_t i,m,d2,D2;
       T *p, *p2;

       for (D2=k=0; k<len; ++k) D2+=M[k]->dim2;

       if (D2!=dim2)
       wblog(FL,"ERR dimension mismatch (%d/%d)",D2,dim2); 

       for (k=0; k<len; ++k) {
          if (M[k]->dim1!=dim1) M[k]->init(dim1,M[k]->dim2);
       }

       for (m=k=0; k<len; ++k) {
          d2=M[k]->dim2; p=data+m; p2=M[k]->data;

          for (i=0; i<dim1; ++i)
          MEM_CPY<T>(p2+d2*i, d2, p+D2*i);

          m+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat index (%d).", dim);
};


template <class T> inline
int wbMatrix<T>::recLess(size_t i1, size_t i2, char lex) const {

    size_t i;
    T *r1=data+i1*dim2, *r2=data+i2*dim2;

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);

    if (lex>0) {
       for (i=0; i<dim2; ++i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }
    else {
       for (i=dim2-1; i<dim2; --i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }

    return 0;
}


template <class T> inline
int wbMatrix<T>::recLessE(size_t i1, size_t i2, char lex) const {

    size_t i;
    T *r1=data+i1*dim2, *r2=data+i2*dim2;

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);

    if (lex) {
       for (i=0; i<dim2; ++i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }
    else {
       for (i=dim2-1; i<dim2; --i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }

    return 1;
}


template <class T> inline
int wbMatrix<T>::recEqual(size_t i1, const T *d) const {
    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,i1,dim1);
    if (data==d) return 1;
    else return (!memcmp(data+i1*dim2, d, dim2*sizeof(T)));
};

template <class T> inline
int wbMatrix<T>::recEqual(size_t i1, size_t i2) const {
    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);
    if (i1==i2) return 1;
    else return (!memcmp(data+i1*dim2, data+i2*dim2, dim2*sizeof(T)));
};

template <class T> inline
T wbMatrix<T>::recDiff2(size_t i1, size_t i2, size_t n) const {
    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",i1,i2,dim1);
    if (i1==i2) return 0;
    return rangeNormDiff2(
       data+i1*dim2, data+i2*dim2, int(n)<0 ? dim2 : n
    );
};

template <class T> inline
T wbMatrix<T>::recDiff2P(size_t i1, const T* b, size_t n) const {
    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",i1,dim1);
    return rangeNormDiff2(
       data+i1*dim2, b, int(n)<0 ? dim2 : n
    );
};


template <class T> inline
char wbMatrix<T>::recCompare(
   size_t i1, size_t i2,
   size_t m,
   char lex,"row-major" (lex>0), else "col-major"
   T eps
 ) const {

   if (i1>=dim1 || i2>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);
   if (!dim2) wblog(FL,
      "ERR %s() got nothing to compare (dim2=%d)",FCT,dim2);

   if (int(m)<0) m=dim2; else
   if (!m || m>dim2) wblog(FL,"ERR size out of bounds (%d/%d)",m,dim2); 

   i1*=dim2; i2*=dim2;
   if (lex<=0) {
      i1+=(dim2-m);
      i2+=(dim2-m);
   }

   if (eps<=0) 
        return Wb::recCompare(data+i1, data+i2, m, lex);
   else return Wb::recCompare(data+i1, data+i2, m, lex, eps);
};


template <class T> inline
char wbMatrix<T>::recCompare(
   size_t i1, const wbvector<T> &v, char lex
 ) const {

   if (v.len!=dim2) wblog(FL,
      "ERR %s() size mismatch (%d,%d)",FCT,v.len,dim2);
   return Wb::recCompare(data+i1*dim2, v.data, dim2, lex);
};


template <class T> inline
char wbMatrix<T>::recCompareP(
  size_t i1, const T* v,
  size_t m,
  char lex,"row-major" (lex>0), else "col-major"
  T eps
) const {

    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,i1,dim1);
   if (!dim2) wblog(FL,
      "ERR %s() got nothing to compare (dim2=%d)",FCT,dim2);

    if (int(m)<0) m=dim2; else
    if (!m || m>dim2) wblog(FL,"ERR size out of bounds (%d/%d)",m,dim2);

    i1*=dim2;
    if (lex<=0) i1+=(dim2-m);

    return Wb::recCompare(data+i1, v, m, lex, eps);
};


template <class T>
inline void wbMatrix<T>::recProd(wbvector<T> &pp) const {
   size_t i,j; T *d=data, p;

   pp.init(dim1);

   if (!dim2) wblog(FL,"ERR %s() records of length 0",FCT);

   for (i=0; i<dim1; ++i, d+=dim2) {
      for (p=d[0], j=1; j<dim2; ++j) p*=d[j];
      pp[i]=p;
   }
};

template <class T>
inline T wbMatrix<T>::recProd(size_t i) const {

   if (i>=dim1 || dim2==0) { wblog(FL,
      "ERR %s() got empty record length (%dx%d; i)",FCT,dim1,dim2,i);
      return 0;
   }

   const T *d=data+i*dim2; T p=d[0];
   for (size_t j=1; j<dim2; ++j) { p*=d[j]; }
   return p;
};

template <class T>
inline T wbMatrix<T>::colSum(size_t c) const {
   T *d=data+c, x=0;

   if (c>=dim2) wblog(FL,"ERR index out of bounds (%d/%d)", c, dim2);
   if (dim1==0) wblog(FL,
   "WRN Sum over column of length zero (%d/%d).", c, dim2);

   for (size_t i=0; i<dim1; ++i, d+=dim2) x+=(*d);
   return x;
}

template <class T> inline
T wbMatrix<T>::recSum(size_t i) const {

   if (i>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d,%dx%d)",FCT,i+1,dim1,dim2);
   if (!dim2) wblog(FL,
      "ERR %s() got empty object (%d,%dx%d)",i,dim1,dim2);

   return addRange(rec(i),dim2);
};

template <class T> inline
wbvector<T>& wbMatrix<T>::recSum(wbvector<T> &s) const {

   s.init(dim1);
   if (!dim2) { if (dim1) wblog(FL,
      "WRN %s() got empty object (%dx%d)",FCT,dim1,dim2);
      return s;
   }

   size_t i=0; const T* d=data;
   for (; i<dim1; ++i, d+=dim2) { s[i]=addRange(d,dim2); }

   return s;
};

template <class T> inline
wbvector<T>& wbMatrix<T>::recSumA(wbvector<T> &s) const {

   s.init(dim1);
   if (!dim2) { if (dim1) wblog(FL,
      "WRN %s() got empty object (%dx%d)",FCT,dim1,dim2);
      return s;
   }

   size_t i=0; const T* d=data; double q;
   for (; i<dim1; ++i, d+=dim2) {
      q=addRange2(d,dim2);
      s[i]=std::sqrt(q);
   }

   return s;
};


template <class T>
wbMatrix<T>& wbMatrix<T>::setRand(double fac, double shift) {
   size_t i,s=dim1*dim2;
   static char first_call=1;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX; isdiag=0;

   if (shift==0.)
        for (i=0; i<s; ++i) data[i]=(T)(fac*rand());
   else for (i=0; i<s; ++i) data[i]=(T)(fac*rand()+shift);

   return *this;
}


template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::setRand(double fac, double shift) {
   size_t i,s=dim1*dim2;
   static char first_call=1;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX; isdiag=0;

   if (shift==0.) {
      for (i=0; i<s; ++i)
      data[i].set(fac*rand(), fac*rand());
   }
   else {
      for (i=0; i<s; ++i)
      data[i].set(fac*rand()+shift, fac*rand()+shift);
   }

   return *this;
}

template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::Symmetrize() {
    size_t i,j,r,s;
    if (dim1!=dim2) wblog(FL,
    "ERR Symmetrize() called with %dx%d matrix !??",dim1,dim2);

    for (j=0; j<dim2; ++j)
    for (i=j+1; i<dim1; ++i) { r=i*dim2+j; s=j*dim2+i;
        data[r]=data[s]=wbcomplex(
            0.5*(data[r].r+data[s].r),
            0.5*(data[r].i-data[s].i)
        );
        data[s].Conj();
    }

    return *this;
};



template <class T>
inline void wbMatrix<T>::setDiagRand(double fac, double shift) {
   size_t i, s=MIN(dim1,dim2);
   static char first_call=1;

   if (data==NULL) return;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX;

   MEM_SET<T>(data,dim1*dim2);

   if (shift==0.)
        for (i=0; i<s; ++i) data[i*dim2+i]=(T)(fac*rand());
   else for (i=0; i<s; ++i) data[i*dim2+i]=(T)(fac*rand()+shift);

   isdiag=1;
}


template <class T>
class mxMatIO {

  public:

    mxMatIO() : ax(0), ar(0), ai(0) {};

    mxMatIO(size_t d1, size_t d2, char cflag=0)
     : ax(0), ar(0), ai(0) { if (d1 && d2) init(d1,d2,cflag); }

    mxMatIO& init(size_t d1, size_t d2, char cflag=0) {
       if (ax || ar || ai) wblog(FL,
          "ERR %s() already got initialzed mxArray !?",FCT);
       if (!d1 || !d2) return *this;

       mwSize dims[]={d1,d2}; Create(2,dims,cflag);
       if (!ax) wblog(FL,
          "ERR %s() failed to allocate mxArray (%dx%d)",FCT,d1,d2);
       ar=(T*)mxGetPr(ax);
       if (!ar) wblog(FL,"ERR %s() got null pointer (%dx%d)",FCT,d1,d2);

       if (cflag) {
          ai=(T*)mxGetPi(ax); if (!ai) wblog(FL,
         "ERR %s() got null pointer (%dx%d; imag)",FCT,d1,d2);
       } else ai=0;

       return *this;
    };

    void Create(size_t r, mwSize *dims, char cflag=0);

    template<class T2>
    mxArray* cpyRange(
       size_t d1, size_t d2, const T2* b,
       size_t stride=1, char raw=0, char iflag=0);

    mxArray *ax;
    T *ar, *ai;

  protected:
  private:
};


template<class T> inline
void mxMatIO<T>::Create(size_t r, mwSize *dims, char cflag){
   ax=mxCreateNumericArray(
      r, dims, mxDOUBLE_CLASS, (cflag ? mxCOMPLEX : mxREAL)
   );
};


template<> inline
void mxMatIO<mxChar>::Create(size_t r, mwSize *dims, char cflag){
   if (cflag) wblog(FL,"ERR %s() cflag=%d will be ignored",FCT,cflag);
   ax=mxCreateCharArray(r,dims);
};

template<> inline
void mxMatIO<int32_T>::Create(size_t r, mwSize *dims, char cflag){
   if (cflag) wblog(FL,"ERR %s() cflag=%d will be ignored",FCT,cflag);
   ax=mxCreateNumericArray(r,dims,mxINT32_CLASS,mxREAL);
};

template<> inline
void mxMatIO<uint32_T>::Create(size_t r, mwSize *dims, char cflag){
   if (cflag) wblog(FL,"ERR %s() cflag=%d will be ignored",FCT,cflag);
   ax=mxCreateNumericArray(r,dims,mxUINT32_CLASS,mxREAL);
};

template<> inline
void mxMatIO<int64_T>::Create(size_t r, mwSize *dims, char cflag){
   if (cflag) wblog(FL,"ERR %s() cflag=%d will be ignored",FCT,cflag);
   ax=mxCreateNumericArray(r,dims,mxINT64_CLASS,mxREAL);
};

template<> inline
void mxMatIO<uint64_T>::Create(size_t r, mwSize *dims, char cflag){
   if (cflag) wblog(FL,"ERR %s() cflag=%d will be ignored",FCT,cflag);
   ax=mxCreateNumericArray(r,dims,mxUINT64_CLASS,mxREAL);
};


template<class T>
template<class T2> inline
mxArray* mxMatIO<T>::cpyRange(
    size_t d1, size_t d2, const T2* b, size_t stride,
    char raw, char iflag
){
    if (!ar && !ai && !ax) {
       if (!raw)
            init(d1,d2);
       else init(d2,d1);
    }

    if (d1 && d2) {
       T *a=ar;
       if (iflag) { a=ai; }
       if (!a) wblog(FL,"ERR %s() uninitialized pointer",FCT); 

       if (!raw) { size_t i,j,k=0;
          for (i=0; i<d1; ++i)
          for (j=0; j<d2; ++j, k+=stride) a[i+j*d1]=(T)b[k];
       }
       else {
          size_t i,k=0, n=d1*d2;
          for (i=0; i<n; ++i, k+=stride) a[i]=(T)b[k];
       }
    }

    return ax;
};


template <class T>
mxArray* wbMatrix<T>::toMx_base(char raw) const {
    return mxMatIO<double>().cpyRange(dim1,dim2,data,1,raw);
};


template <>
mxArray* wbMatrix<char>::toMx_base(char raw) const {
    return mxMatIO<mxChar>().cpyRange(dim1,dim2,data,1,raw);
};

template<>
mxArray* wbMatrix<wbcomplex>::toMx_base(char raw) const {

    char isC = isComplex();
    if (!raw) {
       mxMatIO<double> M(dim1,dim2,isC);

       M.cpyRange(dim1,dim2,&(data[0].r),2,raw); if (isC) {
       M.cpyRange(dim1,dim2,&(data[0].i),2,raw,'i'); }

       return M.ax;
    }
    else {
       mxMatIO<double> M(dim2,dim1,isC);

       M.cpyRange(dim1,dim2,&(data[0].r),2,raw); if (isC) {
       M.cpyRange(dim1,dim2,&(data[0].i),2,raw,'i'); }

       return M.ax;
    }
};


template<class T>
mxArray* wbMatrix<T>::toMx_Struct() const {

   mxArray *S=data->mxCreateStruct(dim1,dim2);

   size_t i,j;
   for (i=0; i<dim1; ++i)
   for (j=0; j<dim2; ++j) data[j+i*dim2].add2MxStruct(S,i+j*dim1);

   return S;
};

template<class T>
mxArray* wbMatrix<T>::toMx_StructP() const {

   if (data==NULL) {
      return mxCreateStructMatrix(dim1,dim2,0,NULL);
   }
   else {
      mxArray *S=(*data)->mxCreateStruct(dim1,dim2);

      size_t i,j;
      for (i=0; i<dim1; ++i)
      for (j=0; j<dim2; ++j) data[j+i*dim2]->add2MxStruct(S,i+j*dim1);

      return S;
   }
};

template<class T>
mxArray* wbMatrix<T>::toMx_CellP() const {

   if (data==NULL) {
      return mxCreateCellMatrix(dim1,dim2);
   }
   else {
      mxArray *S=(*data)->mxCreateCell(dim1,dim2);

      size_t i,j;
      for (i=0; i<dim1; ++i)
      for (j=0; j<dim2; ++j) data[j+i*dim2]->add2MxCell(S,i+j*dim1);

      return S;
   }
};


template<class T>
mxArray* wbMatrix<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbMatrix<T>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      size_t s=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s) wblog(FL,
      "ERR %s() must follow mxCreateCell()\n%lx, %d/%d",
       FCT,S,i+1,s);
   }

   mxSetCell(S,i,toMx());
}



template <class T> inline 
mxArray* wbMatrix<T>::toMx(char raw) const {
   if (raw) wblog(FL,
      "WRN raw flag will be ignored for typeid=%s",getName(typeid(T)).data);

   return toMx_Struct();
};

template <class T> inline 
mxArray* wbMatrix<T>::toMxP(char flag) const {

   if (!flag || flag=='s' || flag=='S') return toMx_StructP();
   if (         flag=='c' || flag=='C') return toMx_CellP();

   wblog(FL,"ERR %s() invalid flag (%s)",FCT,char2Str(flag).data);
   return 0;
};


template <> inline 
mxArray* wbMatrix<unsigned>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<unsigned>::toMxT(char raw) const {
   return mxMatIO<uint32_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<int>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<int>::toMxT(char raw) const {
   return mxMatIO<int32_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<unsigned long>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<unsigned long>::toMxT(char raw) const {
   return mxMatIO<uint64_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<long>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<long>::toMxT(char raw) const {
   return mxMatIO<int64_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<unsigned char>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<unsigned char>::toMxT(char raw) const {
   return mxMatIO<uint8_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<char>::toMx(char raw) const {
   return toMx_base(raw);
};
template <> inline 
mxArray* wbMatrix<char>::toMxT(char raw) const {
   return mxMatIO<int8_T>().cpyRange(dim1,dim2,data,1,raw);
};

template <> inline 
mxArray* wbMatrix<double>::toMx(char raw) const {
   return toMx_base(raw);
};

template <> inline 
mxArray* wbMatrix<wbcomplex>::toMx(char raw) const {
   return toMx_base(raw);
};


template <class T>
void wbMatrix<T>::mat2mx0(mxArray* &a) const {

    double *d=mxGetPr(a);
    a=mxCreateDoubleMatrix(dim2,dim1,mxREAL);

    if (typeid(T)==typeid(double))
       memcpy(d, data, dim1*dim2*sizeof(double));
    else {
        size_t i,n=dim1*dim2;
        for (i=0; i<n; ++i) d[i]=(double)data[i];
    }
}


template <class T>
void wbMatrix<T>::mat2mxc(mxArray* C, unsigned idx) const {

    size_t i,j,k;
    double *d;
    mxArray *a;

    if (!mxIsCell(C)) wblog(FL,"ERR need cell structure on input.");

    if (idx>=mxGetNumberOfElements(C)) wblog(FL,
    "ERR index out of bounds (%d/%d)",idx,mxGetNumberOfElements(C));

    a=mxCreateDoubleMatrix(dim1,dim2,mxREAL);
    mxSetCell(C,idx,a);

    d=mxGetPr(a); k=0;

    for (j=0; j<dim2; ++j) 
    for (i=0; i<dim1; ++i) d[k++]=(double)data[i*dim2+j];
}


template <class T>
void wbMatrix<T>::info(const char *istr) const {

    size_t l=0, n=16; char s[n];
    snprintf(s,n,"%dx%d",dim1,dim2);
    if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
    mexPrintf("  %-12s %-10s @ 0x%lX  double array\n", istr, s, data);
};

template <class T>
void wbMatrix<T>::print(const char *istr, char mflag) const {

    size_t i,n=dim1*dim2;
    mxArray *a;
    double *dd;

    wbMatrix<T> M; (*this).transpose(M); 

    a=mxCreateDoubleMatrix(dim1, dim2,mxREAL);
    dd=mxGetPr(a);

    for (i=0; i<n; ++i)
    dd[i]=(double)M.data[i];

    if (!mflag) {
        if (istr[0])
        mexPrintf("\n%s = [%dx%d double]\n\n", istr, dim1, dim2);
        else mexPrintf("\n");
    }
    else
    mexPrintf("\n%s = [\n", istr[0] ? istr : "ans");

    mexCallMATLAB(0,NULL,1,&a, "disp");

    if (mflag) mexPrintf("];\n");

    mxDestroyArray(a);
}

template <class T>
void wbMatrix<T>::Print(const char *istr, char mflag) const {

    size_t i,j;

    if (!mflag) {
        if (istr[0])
        printf("\n%s = [%dx%d double]\n\n", istr, dim1, dim2);
        else printf("\n");
    }
    else
    printf("\n%s = [\n", istr[0] ? istr : "ans");

    for (i=0; i<dim1; ++i) {
        for (j=0; j<dim2; ++j) printf(" %8.3g", (double)data[i*dim2+j]);
        printf("\n");
    }

    if (mflag) printf("];\n");
}

template <class T>
void wbMatrix<T>::printdata(
    const char *istr, const char *dbl_fmt, const char *rsep
) const {

    size_t i=0,j;
    if (istr && istr[0]) printf("%s = [\n ",istr); else printf("["); 

    for (; i<dim1; ++i) { if (i) printf("%s ",rsep);
    for (j=0; j<dim2; ++j) { printf(dbl_fmt, double(data[i*dim2+j])); }}

    printf(" ]\n");
}


template <class T>
void wbMatrix<T>::appendRow(const wbvector<T> &v) {

    if (dim1) {
       if (dim2!=v.len) wblog(FL,
       "ERR Dimension mismatch (%dx%d : %d)", dim1,dim2, v.len);

       appendRows(1, v.data);
    }
    else RENEW(1,v.len,v.data);
};


template <class T>
void wbMatrix<T>::appendRows(size_t n, const T* v) {

    size_t s=dim1*dim2, ds=n*dim2;
    T *d0=data;

    dim1+=n; if (!n || !dim2) return;

    WB_NEW(data,s+ds);
    MEM_CPY<T>(data,s+ds,s,d0,v);
    WB_DELETE(d0);
};


template<class T>
void wbMatrix<T>::getReal(wbMatrix<double> &R) const {
   wblog(FL,"ERR wbMatrix::getReal not defined for type %s.",
   getName(typeid(T)));
}

template<class T>
void wbMatrix<T>::getImag(wbMatrix<double> &I) const {
   wblog(FL,"ERR wbMatrix::getImag not defined for type %s.",
   getName(typeid(T)));
}

template<class T>
void wbMatrix<T>::set(const wbMatrix<double> &R, const wbMatrix<double> &I) {
   wblog(FL,"ERR wbMatrix::set(R,I) not defined for type %s.",
   getName(typeid(T)));
}

template<class T>
wbMatrix<T>& wbMatrix<T>::Conj() { return *this; }


template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::Conj() {
   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   data[i].i = -data[i].i;
   return *this;
}

template<>
void wbMatrix<wbcomplex>::getReal(wbMatrix<double> &R) const {
   R.init(dim1, dim2);
   
   if (isdiag) {
      size_t i,k,s=MIN(dim1,dim2);
      for (i=0; i<s; ++i) { k=i*dim2+i; R.data[k]=data[k].r; }
      R.isdiag=isdiag;
   }
   else {
      for (size_t s=dim1*dim2, i=0; i<s; ++i)
      R.data[i]=data[i].r;
   }
}

template<>
void wbMatrix<wbcomplex>::getImag(wbMatrix<double> &I) const {
   I.init(dim1, dim2);
   
   if (isdiag) {
      size_t i,k,s=MIN(dim1,dim2);
      for (i=0; i<s; ++i) { k=i*dim2+i; I.data[k]=data[k].i; }
      I.isdiag=isdiag;
   }
   else {
      for (size_t s=dim1*dim2, i=0; i<s; ++i)
      I.data[i]=data[i].i;
   }
}

template<>
void wbMatrix<wbcomplex>::set(
   const wbMatrix<double> &R,
   const wbMatrix<double> &I
){
   if (!R.hasSameSize(I)) wblog(FL,
   "ERR Dimension mismatch (%dx%d; %dx%d)!",R.dim1,R.dim2,I.dim1,I.dim2);

   init(R.dim1,R.dim2);

   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   data[i].set(R.data[i], I.data[i]);
}


template <class T> inline
void matchIndexU(const char *F, int L,
    const wbMatrix<T> &QA, const wbMatrix<T> &QB,
    wbindex &Ib, const char force=1
);

template<class T>
int matchIndex(C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    char lex=1,"row-major" (lex>0), else "col-major"
    INDEX_T *ma=NULL, INDEX_T *mb=NULL,
    T eps=0
);

template <class T> inline
int matchIndex(
    const wbvector<T> &qA, const wbvector<T> &qB,
    wbindex &Ia, wbindex &Ib,
    char lex=1,"row-major" (lex>0), else "col-major"
    INDEX_T *ma=NULL, INDEX_T *mb=NULL,
    T eps=0
){
    wbMatrix<T> QA, QB;  int r;
    QA.init2ref(qA,'t');
    QB.init2ref(qB,'t'); if (QA.dim2!=1 || QB.dim2!=1) wblog(FL,"ERR");
    
    r=matchIndex(QA,QB,Ia,Ib,lex,ma,mb,eps);
    return r;
};


template <class T> inline
int matchSortedIdx(
    const T *da, size_t lda, size_t na,
    const T *db, size_t ldb, size_t nb, wbindex &Ia, wbindex &Ib,
    INDEX_T m=-1,
    char lex=1
){
    INDEX_T m1=0,m2=0;
    return Wb::matchSortedIdx(da,lda,na, db,ldb,nb, Ia,Ib,m,lex,&m1,&m2);
};

template <class T> inline
int matchSortedIdx(C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    INDEX_T m,
    char lex,"row-major" (lex>0), else "col-major"
    INDEX_T *ma=NULL, INDEX_T *mb=NULL,
    T eps=0
){
    return Wb::matchSortedIdx(
       QA.data, QA.dim2, QA.dim1,
       QB.data, QB.dim2, QB.dim1, Ia,Ib, m, lex, ma, mb, eps
    );
};

template <class T> inline
int matchSortedIdx(
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    INDEX_T m=-1,
    char lex=1
){
    INDEX_T m1=0,m2=0;
    return Wb::matchSortedIdx(
       QA.data, QA.dim2, QA.dim1,
       QB.data, QB.dim2, QB.dim1, Ia,Ib, m,lex,&m1,&m2
    );
};

template <class T> inline
int matchSortedIdxU(
    const char *F, int L,
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    INDEX_T m=-1, char lex=1
);

template <class T> inline
int matchSortedIdxU(
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib, INDEX_T m=-1
){  return matchSortedIdxU(FL,QA,QB,Ia,Ib,m); };


#endif

