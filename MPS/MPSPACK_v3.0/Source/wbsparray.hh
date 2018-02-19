#ifndef __WB_SPARRAY_HH__
#define __WB_SPARRAY_HH__


//===================================================================//
// wbsparray - N-D sparse array class                                //
//                                                                   //
//    generic N-dimensional class that allows permutations           //
//    and reshaping. Matrix elements are stored as follows:          //
//                                                                   //
//    [ full index element 1, value1                                 //
//      full index element 2, value2                                 //
//      ...                                                          //
//    ];                                                             //
//                                                                   //
// Wb,Jan05,11                                                       //
//===================================================================//

bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k=0);

bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k,
   const mxArray **as, const mxArray **ad, const mxArray **ai
);

bool mxIsWbsparray(const mxArray *a) { return mxIsWbsparray(0,0,a); }

#ifdef WB_SPARSE_CLOCK
WbClock wbc_sparse_cont("sp::contract");
WbClock wbc_sparse_cntf("spF::contract");
WbClock wbc_sparse_cnt2("sp2::contract");
WbClock wbc_sparse_cnt3("sp3::contract");
WbClock wbc_sparse_cnt4("sp4::contract");
WbClock wbc_sparse_cmpr("sp::compress");
#endif



#define wbSparray    wbsparray<double>
#define wbSparrayTD  wbsparray<TD>

template<class TD>
class wbsparray {

  public:

    wbsparray() : isref(0) {};

    wbsparray(SPIDX_T d1) : isref(0) { if (d1) init(d1); };
    wbsparray(SPIDX_T d1, SPIDX_T d2) : isref(0) { init(d1,d2); };
    wbsparray(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3) : isref(0) { init(d1,d2,d3); };
    wbsparray(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T d4)
     : isref(0) { init(d1,d2,d3,d4); };

    wbsparray(const wbvector<SPIDX_T> &d, SPIDX_T l=0)
     : isref(0) { init(d,l); };

    wbsparray(char *sidx, SPIDX_T l=0)
     : isref(0) { init(Str2Idx(sidx),l); };

    wbsparray(const wbsparray &a)
     : isref(0) { SIZE=a.SIZE; IDX=a.IDX; D=a.D; };

    wbsparray(const char *F, int L, const mxArray *a)
     : isref(0) { init(F_L,a); };

    wbsparray(const mxArray *a) : isref(0) { init(FL,a); };


    wbsparray& init(
       const char *F, int L, const wbarray<TD> &A, TD eps=0
    );

    wbarray<TD>& toFull(wbarray<TD> &A) const;

    template<class T2>
    wbvector<T2>& toFull(wbvector<T2> &A) const;

    template<class DB>
    wbsparray& init(const wbsparray<DB>& a) {
       if ((void*)this!=(void*)&a) {
          SIZE=a.SIZE; IDX=a.IDX; D=a.D;
          isref=0;
       }
       return *this;
    };

    SPIDX_T catRecs(SPIDX_T i0, const wbsparray &a, SPIDX_T id);

    wbsparray& init(const wbvector<SPIDX_T> &S0, SPIDX_T l=0) {
       SIZE=S0; D.init(l);
       if (S0.len) IDX.init(l,S0.len);
       else { 
          if (l>1) wblog(FL,
             "WRN %s() use initDiag() with S=[] and nz=%d",FCT,l);
          IDX.init();
       }
       return *this;
    };

    template<class I2>
    wbsparray& init(const wbvector<I2> &S_, SPIDX_T l=0) {
       wbvector<SPIDX_T> S0(S_);
       return init(S0,l);
    };

    wbsparray& init() {
       SIZE.init(); IDX.init(); D.init(); isref=0;
       return *this;
    };

    wbsparray& init(SPIDX_T d1) {
       if (d1) SIZE.init(1,&d1); else SIZE.init();
       IDX.init(); D.init(); isref=0;
       return *this;
    };

    wbsparray& init(SPIDX_T d1, SPIDX_T d2) {
       SPIDX_T s[]= {d1,d2}; SIZE.init(2,s); IDX.init(); D.init();
       return *this;
    };

    wbsparray& init(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3) {
       SPIDX_T s[]= {d1,d2,d3}; SIZE.init(3,s); IDX.init(); D.init();
       return *this;
    };

    wbsparray& init(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T d4) {
       SPIDX_T s[]= {d1,d2,d3,d4}; SIZE.init(4,s); IDX.init(); D.init();
       return *this;
    };

    wbsparray& initz() {
       IDX.init(0,SIZE.len); D.init();
       return *this;
    };

    wbsparray& initz(SPIDX_T d1, SPIDX_T n) {
       SIZE.init(1,&d1); IDX.init(n,SIZE.len); D.init(n);
       return *this;
    };

    wbsparray& initz(SPIDX_T d1, SPIDX_T d2, SPIDX_T n) {
       SPIDX_T s[]= {d1,d2}; SIZE.init(2,s); IDX.init(n,SIZE.len); D.init(n);
       return *this;
    };

    wbsparray& initz(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T n) {
       SPIDX_T s[]= {d1,d2,d3}; SIZE.init(3,s); IDX.init(n,SIZE.len); D.init(n);
       return *this;
    };

    wbsparray& init_rzd(SPIDX_T r, SPIDX_T n) {
       SIZE.init(r); IDX.init(n,r); D.init(n); isref=0;
       return *this;
    };

    wbsparray& init_nnz(SPIDX_T n=0, char kflag=0) {
       if (!SIZE.len) wblog(FL,
          "ERR sparse::%s() got empty size SIZE!",FCT);
       if (isref) {
          if (SIZE.isref) wblog(FL,"ERR %s() got isref "
             "(%d,%d,%d,%d)",FCT,isref,SIZE.isref,IDX.isref,D.isref);
          isref=0;
       }

       if (kflag) {
          if (SIZE.len%2) wblog(FL,
             "ERR %s() even rank expected (%d)",FCT,SIZE.len); 
          IDX.init(n,SIZE.len/2);
       }
       else { IDX.init(n,SIZE.len); }

       D.init(n); return *this;
    };

    wbsparray& setRand(double x=1, char pnflag=0);
    wbsparray& setRandDist();

    template<class DB>
    wbsparray& init(const char *F, int L,
       SPIDX_T d1, SPIDX_T d2, const wbsparray<DB> &B);

    template<class DB>
    wbsparray& init(const char *F, int L,
       SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, const wbsparray<DB> &B);

    template<class DB>
    wbsparray& init(const char *F, int L,
       SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T d4, const wbsparray<DB> &B);

    wbsparray& initEmpty(const wbsparray &a) {
       if (this!=&a) SIZE=a.SIZE; 
       IDX.init(0,SIZE.len); D.init();
       return *this;
    };

    wbsparray& init(
       const char *F, int L, const mxArray *a, unsigned k=0);

    wbsparray& init2ref(const wbsparray &a) {
       SIZE=a.SIZE; IDX.init2ref(a.IDX); D.init2ref(a.D);
       isref=1; return *this;
    };

    wbsparray& unRef() {
       if (!isref) return *this;
       IDX.unRef(); D.unRef(); isref=0;
       return *this;
    };

    wbsparray& initCAT(const char *F, int L,
       wbvector< wbsparray* > &X,
       wbvector<SPIDX_T> *I0=NULL, char del=0
    );

    wbsparray& initCAT(const char *F, int L,
       wbvector< wbsparray > &X, const wbvector<SPIDX_T> *I0=NULL
    ){
       wbvector< wbsparray* > xp(X.len);
       for (unsigned i=0; i<X.len; ++i) xp[i]=(&X[i]);
       return initCAT(F_L,xp,I0);
    };


    wbsparray& initCAT(
       const char *F, int L, const wbvector< const wbsparray* > &X,
       unsigned dim
    );

    wbsparray& initCAT(
       const char *F, int L, const wbvector< wbsparray > &X,
       unsigned dim
    ){
       wbvector< const wbsparray* > xp(X.len);
       for (unsigned i=0; i<X.len; ++i) xp[i]=(&X[i]);
       return initCAT(F,L,xp,dim);
    };

    wbsparray& Cat(unsigned dim,
       const wbsparray &B, const char *F=NULL, int L=0
    );

    wbsparray& init_kron(
       const wbvector<SPIDX_T> &sa, const wbvector<SPIDX_T> &sb,
       SPIDX_T nab, char kflag
    ){
       setSIZE_kron(FL,sa,sb);
       return init_nnz(nab,kflag);
    };

    void setSIZE_kron(const char *F, int L,
       const wbvector<SPIDX_T> &sa, const wbvector<SPIDX_T> &sb,
       char kflag=0);
    void setSIZE_kron(const char *F, int L, char kflag=0);

    wbsparray& setRec(SPIDX_T i, SPIDX_T j) {
       if (i>=D.len || j>=D.len || IDX.dim1!=D.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d/%d,%d",FCT,i,j,IDX.dim1,D.len);
       if (i!=j) { D[i]=D[j]; IDX.recSet(i,j); }
       return *this;
    };

    template<class T2>
    wbsparray& setRec(SPIDX_T i, SPIDX_T i1, SPIDX_T i2, const T2& x) {
       if (IDX.dim2!=2) wblog(FL,"ERR %s() "
          "rank-2 array required (%d: %s)",FCT,SIZE.len,sizeStr().data);
       if (i>=D.len || IDX.dim1!=D.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d,%d)",FCT,i,IDX.dim1,D.len);
       SPIDX_T *I=IDX.rec(i);
       I[0]=i1; I[1]=i2; D[i]=x;
       return *this;
    };

    wbsparray& setRecP(SPIDX_T i, const SPIDX_T *I, const TD &x) {
       if (i>=D.len && IDX.dim1!=D.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d /%d,%d",FCT,i,IDX.dim1,D.len);
       IDX.recSetP(i,I); D[i]=x;
       return *this;
    };

    wbsparray& setRec_kron(
       SPIDX_T l, TD x, unsigned r, SPIDX_T *i1, SPIDX_T *i2,
       unsigned kflag=0
    );

    wbsparray& Resize(SPIDX_T l) {
       IDX.Resize(l,IDX.dim2); D.Resize(l);
       return *this;
    };

    wbsparray& operator= (const wbsparray &a) {
        if ((void*)this!=(void*)&a) init(a);
        return *this;
    };
    template<class DB>
    wbsparray& operator= (const wbsparray<DB> &a) {
        if ((void*)this!=(void*)&a) init(a);
        return *this;
    };

    wbsparray& save2(wbsparray &a, char ref=0) {
        if (this!=&a) {
           SIZE.save2(a.SIZE); IDX.save2(a.IDX,ref); D.save2(a.D);
           a.isref=isref; isref=0;
        }
        return a;
    };

    SPIDX_T sub2ind(SPIDX_T k) const {
       return Wb::sub2ind(SIZE.data,IDX.rec(k),SIZE.len);
    };

    void ind2sub(SPIDX_T k, SPIDX_T *I) const {
       Wb::ind2sub(k, SIZE.data, I, SIZE.len);
    };

    void ind2sub(SPIDX_T k, wbindex &I) const {
       if (I.len!=SIZE.len) I.init(SIZE.len);
       Wb::ind2sub(k, SIZE.data, I.data, SIZE.len);
    };

    wbsparray<TD>& make2D(const wbindex &ic, wbsparray<TD> &B,
      char pos=1
    ) const {
      if (B.D.len) B.init(); B.D=D;
      IDX.toIndex2D(SIZE,ic,B.IDX,pos);"col-major"

      SPIDX_T s2, s1=SIZE.prod(ic.data,ic.len,s2);

      if (pos==1)
           { SPIDX_T s[2]={s1,s2}; B.SIZE.init(2,s); }
      else { SPIDX_T s[2]={s2,s1}; B.SIZE.init(2,s); }

      return B.Compress();
    };

    
    wbsparray& initDiag(SPIDX_T d) {
       D.init(d); SIZE.init(); IDX.init();
       return *this;
    };

    bool isDiagMatrix(TD eps=1E-14) const;

    char isDiag(const char *F=0, int L=0) const {
       if (!SIZE.len) {
          if (IDX.dim1==0 && IDX.dim2==0) return (D.len>=1);
       }
       checkSize(F_LF);
       return 0;
    };

    wbvector<TD>& getDiag(const char *F, int L, wbvector<TD> &dd);
    wbsparray& getDiag(const char *F, int L, wbsparray &A);

    wbsparray& diag2reg(const char *F=0, int L=0, unsigned r=-1);
    wbsparray& diag2reg(
       const char *F, int L, wbsparray &A, unsigned r=-1) const;

    wbsparray& toDiag(const char *F=0, int L=0) {
       wbsparray X; getDiag(F_L,X);
       return X.save2(*this);
    };

    int toDiagSafe(const char *F, int L, TD eps=1E-14) {
       if (!isDiagMatrix(eps)) { if (F) {
          wblog(F_L, "ERR %s() got non-diagonal matrix",FCT); }
          return  0;
       }
       wbsparray X; getDiag(F_L,X);
       X.save2(*this); return 1;
    };

    wbsparray& initIdentity(SPIDX_T d, TD dval=1){
       initDiag(d); D.set(dval);
       return *this;
    };

    wbsparray& initIdentity(
        const wbvector<SPIDX_T> &S0, int m=-1, sSPIDX_T k=0, char fflag='f');
    wbsparray& initIdentityB3(
        SPIDX_T d1, SPIDX_T d2, SPIDX_T D,
        sSPIDX_T i0=0, char fflag='f'
    );

    unsigned rank(const char *F=NULL, int L=0) const {
       checkSize(F_LF);
       if (!SIZE.len) {
          return (isDiag() ? 2 : 0);
       }
       return SIZE.len;
    };

    SPIDX_T dim() const {
       if (isDiag()) return D.len;
       if (SIZE.len!=2 || SIZE[0]!=SIZE[1]) wblog(FL,
          "ERR %s() got %s array",FCT,sizeStr().data); 
       return SIZE[0];
    };

    SPIDX_T dim(unsigned k) const {
       if ((--k)>9999U) wblog(FL,
          "ERR %s() input dim is 1-based (%d)",FCT,k+1);
       return dim0(k);
    };

    SPIDX_T dim0(unsigned k) const {
       if (isDiag() && k<2) return D.len;
       if (k<SIZE.len) return SIZE[k];

       wblog(FL,"ERR %s() "
         "dim=%d out of bounds (%s)",FCT,k+1,sizeStr().data);
       return 0;
    };

    SPIDX_T matDim(unsigned k=1) {
       if (!k || k>2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k);
       if (isDiag(FL)) { return D.len; } else
       if (k==1) { return prodRange(SIZE.data,SIZE.len/2); }
       else {
          unsigned l=SIZE.len/2;
          return prodRange(SIZE.data+l,SIZE.len-l);
       }
    };

    SPIDX_T maxDim(unsigned r=-1) const {
       if (int(r)<0) {
          if (!SIZE.len)
               { return (isDiag(FL) ? D.len : 0); }
          else { return SIZE.max(); }
       }
       else {
          if (SIZE.len) {
             if (!r || r>SIZE.len) wblog(FL,
                "ERR %s() invalid dim range (%d/%d)",FCT,r,SIZE.len);
             unsigned i=0, d=SIZE[0];
             for (; i<r; ++i) { if (d<SIZE[i]) d=SIZE[i]; }
             return d;
          }
          else {
             if (isDiag(FL)) {
                if (!r || r>2) wblog(FL,
                   "ERR %s() invalid dim range (%d/2)",FCT,r);
                return D.len;
             }
             else {
                if (r) wblog(FL,
                   "ERR %s() invalid dim range (%d/0)",FCT,r);
                return 0;
             }
          }
       }
    };

    TD norm2(char tnorm=0) const {
       return overlap(D.data,D.data,D.len,1,tnorm);
    };

    wbstring sizeStr() const {
       if (SIZE.len) return SIZE.toStrD();" (sparse)");
       else if (D.len) {
          wbstring s(16); s.pushf(FL,"[%dx%d]",D.len,D.len);
          return s;
       }
       else return wbstring("[]");
    };

    double SkipTiny(TD eps=1E-14) {
       SPIDX_T l=0, i=0; TD x2=0;
       Wb::scale_eps(eps,D.data,D.len);

       for (i=0; i<D.len; ++i) {
          if (ABS(D[i])>eps) { if (l<i) { setRec(l,i); }; ++l; }
          else x2+=(D[i]*D[i]);
       }
       if (l<i) {
          IDX.RESIZE_RECS(l);
          D.RESIZE(l);
       }
       return x2;
    };

    bool isVector(const char *F=0, int L=0) const {
       unsigned r=rank(F_L); if (r==1) return 1;
       return (SIZE.count(1)>=(r-1));
    };

    bool isMatrix(const char *F=0, int L=0) const { return rank(F_L)==2; };
    bool isRank(unsigned r) const { return (rank(FL)==r); };

    bool isSMatrix(const char *F=0, int L=0, SPIDX_T *d=NULL) const {
       unsigned r=rank(F_L);
       if (r!=2 || (SIZE.len && SIZE[0]!=SIZE[1])) return 0;
       if (d) {
          SPIDX_T d0=(SIZE.len ? SIZE[0] : D.len);
          if (sSPIDX_T(*d)>=0) { if((*d)!=d0) return 0; }
          else (*d)=d0;
       }
       return 1;
    };
    bool isSMatrix(const char *F, int L, const SPIDX_T &d) const {
       unsigned r=rank(F_L); SPIDX_T d0=(SIZE.len ? SIZE[0] : D.len);
       return (
         r!=2 || (SIZE.len && SIZE[0]!=SIZE[1]) || d!=d0
       ) ? 0 : 1;
    };

    bool isEmpty() const {
       if (SIZE.len || IDX.dim1 || D.len) return 0;
       return 1;
    };

    bool isZero(TD eps=0, char flag=0) const;

    SPIDX_T nnz() const {
       if (isScalar(FL))
            return (D.len==1 && D.data[0]!=0 ? 1 : 0);
       else return D.len;
    };

    SPIDX_T nnz(const TD eps1, const TD eps2=1E9) const {
       SPIDX_T i=0, n=0; TD a,e=0;
       for (i=0; i<D.len; ++i) {
          a=ABS(D[i]); if (a>eps1) ++n;
          else if (a>eps2 && e<a) e=a;
       }
       if (e) wblog(FL,
          "WRN %s() got small value(s) %.3g (%g;%g)",FCT,
           double(e), double(eps1), double(eps2));
       return n;
    };

    SPIDX_T nnz_diag() const {
       if (isDiag(FL)) return D.nnz();
       if (IDX.dim2!=2) wblog(FL,
          "ERR %s() requires rank-2 object (S=%s)",FCT,sizeStr().data);
       SPIDX_T i=0, n=0; const SPIDX_T *I=IDX.data;
       for (; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1] && D[i]!=0) ++n; }
       return n;
    };

    SPIDX_T numel() const {
       return (isDiag(FL) ? D.len*D.len : SIZE.prod(0));
    };

    TD norm(char tnorm=0) const {
       return SQRT(overlap(D.data,D.data,D.len,1,tnorm));
    };

    template<class T2>
    wbvector<T2>& norm2vec(unsigned k,
       wbvector<T2> &A, char tnorm=0
     ) const;

    void swap(wbsparray &B) {
        if (this!=&B) {
           if (isref || B.isref) wblog(FL,
              "ERR %s() got arrays with isref=(%d,%d)",FCT,isref,B.isref);
           SIZE.swap(B.SIZE);
           IDX.swap(B.IDX);
           D.swap(B.D);
        }
    };

    wbsparray& TimesEl(const wbsparray &B, char tnorm=0);

    TD dotProd(
      const char *F, int L, const wbsparray &B, char tnorm=0) const;

    TD froNorm2(const wbsparray &B) const { return dotProd(0,0,B); };

    double normDiff2(const wbsparray &B) const {
       wbsparray X(*this); X-=B; return X.norm2();
    };
    double normDiff (const wbsparray &B) const {
       wbsparray X(*this); X-=B; return X.norm();
    };

    TD Normalize(char tnorm=0, char qflag=0);
    TD NormalizeCol(SPIDX_T k, char tnorm=0, char qflag=0);

    wbsparray& OrthoNormalizeCols(
       const char *F=NULL, int L=0, char tnorm=0,
       char qxflag=0, TD eps=1E-14, unsigned np=1
    );

    wbsparray& setCol(
       const char *F, int L, SPIDX_T k, const wbsparray &v);

    wbsparray& setCol(SPIDX_T k, const wbsparray &v) {
       return setCol(0,0,k,v);
    };

    double sparsity() const {
       if (isScalar(FL))
            return (D.len==1 && D.data[0]!=0 ? 1 : 0);
       else return D.len/double(SIZE.len ? SIZE.prod(0) : 1);
    };

    int checkTrailingSingletons(const char *F=0, int L=0, unsigned *r=NULL) const;

    int skipTrailingSingletons(const char *F=0, int L=0, unsigned r=-1);
    int skipTrailingSingletons(
       const char *F, int L, wbsparray &c, unsigned r=-1) const;

    bool hasSingletons(const wbindex &I) const;

    wbsparray& skipSingletons(
      const char *F, int L, const wbindex &I, wbsparray &A
    ) const;

    wbsparray& addTrailingSingletons(
       const char *F, int L, unsigned r);

    wbsparray& initScalar(const TD x=0, unsigned r=-1) {
       if (int(r)<=0) {
          if (D.len!=1) D.init(1); D[0]=x;
          if (SIZE.len) SIZE.init();
          IDX.init();
       }
       else {
          SIZE.init(r).set(1); 
          if (x) {
             if (D.len!=1) D.init(1); D[0]=x;
             IDX.init(1,r);
          }
          else { IDX.init(0,r); D.init(); }
       }
       isref=0; return *this;
    };

    char isScalar(const char *F=0, int L=0, char zflag=0) const {
       if (!SIZE.len) {
          if (IDX.dim1 || IDX.dim2) sperror_this(F_L);
          return (zflag ? D.len<=1 : D.len==1);
       }
       return SIZE.allEqual(1);
    };

    wbsparray& toScalar(const char *F=0, int L=0, unsigned r=-1);

    TD getScalar(const char *F=0, int L=0) const;

    char checkSize(
      const char *F=0, int L=0, const char *fct=0, const char *istr=0
    ) const;

    char checkSize2(const char *F, int L, const wbsparray &B) const;

    void check_IDX_range(const char *F=0, int L=0) const;

    bool isProptoId(TD &q,TD eps=1E-14) const;
    bool isIdentityMatrix(TD eps=1E-14) const {
       TD q=1; return isProptoId(q,eps);
    };

    char isIdentity(TD eps=1E-14) const;
    char isIdentity(const wbperm &P, TD *dval=NULL, TD eps=1E-14) const;

    bool hasSize(SPIDX_T d1) const;
    bool hasSize(SPIDX_T d1, SPIDX_T d2) const;
    bool hasSize(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3) const;

    template<class T>
    wbvector<T>& getSize(wbvector<T> &S) const {
       if (SIZE.len) { return S.initT(SIZE); }
       if (D.len && isDiag()) { return S.init(2).set(D.len); }
       return S.init();
    };

    char sameSize(const wbsparray &, char strict=0) const;
    char sameSizeUp2Singletons(
       const wbvector<SPIDX_T> &S, char strict=0) const;

    char matchWithSingletons(const char *F, int L,
      const wbvector<SPIDX_T> &S, wbvector<int> &M
    ) const;

    char sameUptoFac(
      const wbsparray &B, TD *fac=NULL, char lex=-1,
      TD eps=1E-12) const;

    bool operator== (const wbsparray &B) const {
       if (this==&B) return 1;
       if (SIZE!=B.SIZE || D.len!=B.D.len) return 0;
       if (IDX==B.IDX && D==B.D) return 1;
       return (normDiff2(B)<1E-12);
    };

    bool operator!= (const wbsparray &B) const {
       return !(*this==B);
    };

    wbsparray& Plus(const char *F, int L,
       const wbsparray &B,
       TD bfac=1, TD afac=1, char strict=0
    );

    TD sum() const { return D.sum(); };
    TD aMin() const { return D.aMin(1); };
    TD aMax() const { return D.aMax(); };

    SPIDX_T findRecSortedP(SPIDX_T *idx, SPIDX_T n,
       char lex=-1
    );

    wbsparray& Eldiv(const char *F, int L, const wbvector<TD> &x_) {
       SPIDX_T i=0, j;
       const SPIDX_T *idx=IDX.data;
       const TD *x=x_.data;

       if (!isVector()) wblog(F_L,
          "ERR %s() invalid vector (%s)",FCT,sizeStr().data);
       for (; i<SIZE.len; ++i) { if (SIZE[i]>1) { idx+=i; break; }}

       for (i=0; i<IDX.dim1; ++i, idx+=IDX.dim2) { j=(*idx);
          if (j>=x_.len) wblog(F_L,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,x_.len);
          if (x[j]==0) wblog(F_L,
             "ERR %s() div/0 !??  (%g)",FCT,double(x[j]));
          D.data[i]/=x[j];
       }
       return *this;
    };

    wbsparray& operator+=(const wbsparray &B) { return Plus(FL,B,+1); }
    wbsparray& operator-=(const wbsparray &B) { return Plus(FL,B,-1); }

    wbsparray& Conj() { D.Conj(); return *this; }

    wbsparray& operator*=(TD x) {
       if (isref) {
          D*=x;
       }
       else {
          if (x!=0) { if (x!=1) D*=x; }
          else init_nnz();
       }
       return *this;
    };

    wbsparray& operator/=(TD x) {
       if (x==0) wblog(FL,"ERR %s() div/0",FCT); 
       if (x!=1) D/=x;
       return *this;
    };

    const TD& operator[](SPIDX_T i) const {
       if (i>=D.len) {
          wblog(FL,"ERR %s() index out of bounds (%d/%d; %s)",
          FCT,i,D.len,sizeStr().data);
       }
       return D.data[i];
    };

    TD& operator[](SPIDX_T i) {
       if (i>=D.len) {
          if (i==0 && IDX.dim1==0) {
             if (SIZE.len && SIZE.allEqual(1)) {
                IDX.init(1,SIZE.len);
                D.init(1); return D.data[i];
             }
             else if (!SIZE.len && IDX.dim2==0) {
                IDX.init();
                D.init(1); return D.data[i];
             }
          }
          wblog(FL,"ERR %s() index out of bounds (%d/%d; %s)",
          FCT,i,D.len,sizeStr().data);
       }
       return D.data[i];
    };

    TD value(const SPIDX_T *idx, SPIDX_T n, char lex=-1
     ) const {
       SPIDX_T i=IDX.findRecSorted(idx,n,lex);
       return (sSPIDX_T(i)<0 ? TD(0.) : D.data[i]);
    };

    wbstring info2Str(const char *istr=0) const {
       wbstring s(32);
       if (istr) {
          s.pushf(FL,"%s(%s; %dx%d; %d)",istr,
          sizeStr().data, IDX.dim1,IDX.dim2,D.len);
       }
       else {
          s.pushf(FL,"%s; %dx%d; %d",
          sizeStr().data, IDX.dim1,IDX.dim2,D.len);
       }
       return s;
    };

    void info(const char* istr="") const { return info(0,0,istr); };
    void info(const char* F, int L, const char* ="") const;
    void print(const char* istr="", const char *fmt="") const;

    void print(
       const char *F, int L,
       const char *istr="", const char *fmt=""
    ) const;

    void printdata(const char *istr="ans", const char *fmt="") const;

    void sperror_this(
       const char *F, int L, const char *fct=0, const char *istr=0
     ) const {
       wblog(F_L,"ERR %s() %s S=%s (%dx%d; %d) %s",
          fct ? fct : "sparse", istr ? istr : "got",
          SIZE.len ? SIZE.toStrD().data : "[]",
          IDX.dim1,IDX.dim2,D.len, isref ? " **REF** ":""
       );
    };

    TD trace() const;

    template<class T2>
    wbvector<T2>& trace(unsigned k, wbvector<T2> &A) const;

    bool isSym_aux(
      const char *file, int line, const char *fct,
      const wbsparray &B, TD eps, TD* xref,
      const char symflag='s',
      const char lflag=1
    ) const;

    bool isHConj(
      const wbsparray &B, TD eps=1E-12, TD* xref=NULL
    ) const { return isSym_aux(FLF, B, eps, xref, 's'); };

    bool isHConj(
      TD eps=1E-12, TD* xref=NULL
    ) const { return isSym_aux(FLF,*this,eps,xref,'s'); };

    bool isAHerm(
      const wbsparray &B, TD eps=1E-12, TD* xref=NULL
    ) const { return isSym_aux(FLF, B, eps, xref, 'a'); };

    bool isAHerm(
      TD eps=1E-12, TD* xref=NULL
    ) const { return isSym_aux(FLF,*this,eps,xref,'a'); };

    wbsparray& Sort(char lex=-1) {
       if (IDX.data) {
          if (IDX.dim1!=D.len) wblog(FL,
             "ERR %s() %d/%d !??",FCT,IDX.dim1,D.len);
          wbperm P; IDX.sortRecs(P,+1,lex); D.Select(P);
       }
       return *this;
    };

    wbsparray& Compress(
       const char *F=0, int L=0, TD eps=0,
       char lex=-1,
       const wbperm *P=NULL
    );

    bool isCompressed(
      const char *F, int L, char lex=-1, char dir=0
    ) const { return IDX.isUniqueSorted(dir,lex); };

    int splitSparseCM(
       const char *F, int L, wbvector< wbsparray > &X,
       wbvector<SPIDX_T> *I0=NULL);

    wbsparray& reshape(
       const wbvector<SPIDX_T> &S, wbsparray &A) const;

    wbsparray& Reshape(const wbvector<SPIDX_T> &S){
       wbsparray X; reshape(S,X); return X.save2(*this);
    };

    wbsparray& Reshape(SPIDX_T s1, SPIDX_T s2){
       wbsparray X; wbvector<SPIDX_T> S(2); S[0]=s1; S[1]=s2;
       reshape(S,X); return X.save2(*this);
    };

    wbsparray& Reshape(SPIDX_T s1, SPIDX_T s2, SPIDX_T s3){
       wbsparray X;  wbvector<SPIDX_T> S(3); S[0]=s1; S[1]=s2; S[2]=s3;
       reshape(S,X); return X.save2(*this);
    };

    wbsparray& select0(
       const WBINDEX &, unsigned dim, wbsparray&) const;
    wbsparray& Select0(const WBINDEX &I, unsigned dim) {
       wbsparray X; save2(X); X.select0(I,dim,*this);
       return *this;
    };

    wbsparray& Permute(const char* s, char iflag=0);
    wbsparray& Permute(const wbperm&, char iflag=0);
    wbsparray& permute(const char* s, wbsparray&, char iflag=0) const;
    wbsparray& permute(const wbperm&, wbsparray&, char iflag=0) const;

    wbsparray& transpose(const char *F, int L, wbsparray &A) const;
    wbsparray& Transpose(const char *F=0, int L=0) {
        wbsparray X;
        transpose(F_L,X); X.save2(*this);
        return *this;
    };

    wbsparray& ColPermute(const wbperm &P, char iflag=0);
    wbsparray& MatPermute(const wbperm &P, char iflag=0);

    wbsparray& mtimes(
        const char *F, int L, const wbsparray&B, wbsparray&C,
        char sflag=0, const TD eps=0
     ){
        const wbsparray &A=*this;
        if (A.rank()!=2 || B.rank()!=2) wblog(FL,
           "ERR %s() requires rank-2 wbsparrays (%s * %s)",
           FCT,A.sizeStr().data,B.sizeStr().data);
        wbindex ia(1), ib(1); ia[0]=1; ib[0]=0;

        return contract(F,L,ia,B,ib,C,wbperm(),1,0,sflag,eps);
    };

    char contract_scalar(
       const char* F, int L,  const wbindex &ica,
       const wbsparray &B, const wbindex &icb, wbsparray &C,
       const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
    ) const;

    wbsparray& contract_diag_diag(
       const char* F, int L,  const wbindex &ica,
       const wbsparray &B, const wbindex &icb, wbsparray &C,
       const wbperm &pfinal, const TD& afac, const TD& cfac
    ) const;

    char contract_check_2full(
       const char* F, int L,  const wbindex &ica,
       const wbsparray &B, const wbindex &icb, wbsparray &C,
       const wbperm &pfinal, const TD& afac, const TD& cfac,
       char flag, double eps, SPIDX_T *nnz_C=NULL
    ) const;

    wbsparray& contract(const char *F, int L,  const wbindex &i1,
        const wbsparray &B, const wbindex &i2, wbsparray &C,
        const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.,
        char sflag=0, const double eps=0
    ) const;

    wbsparray& contract(const char *F, int L,  const char *s1,
        const wbsparray &B, const char *s2, wbsparray &C,
        const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.,
        char sflag=0, const double eps=0
    ) const {
        return contract(F,L,
        wbindex(s1,1),B,wbindex(s2,1),C,pfinal,afac,cfac,sflag,eps);
    };

    wbarray<TD>& contract(const char* F, int L, const wbindex &ica,
       const wbsparray &B, const wbindex &icb, wbarray<TD> &C,
       const wbperm &pfinal, const TD& afac, const TD& cfac
    ) const;

    wbsparray& contract(const char *F, int L,
       unsigned i1, const wbvector<TD> &B, wbsparray &C
    ) const;

    wbsparray& comm(const char *F, int L,
      const wbsparray &B, wbsparray &C,
      char aflag='N', char bflag='N'
    ) const;

    wbsparray& tensorProd(
      const char *F, int L, const wbsparray& B, wbsparray& C,
      char aflag='N', char bflag='N', char kflag=0
    ) const;

    wbsparray& Kron(const wbsparray& B,
      char aflag='N', char bflag='N', char kflag='k'
    ){
       wbsparray X; tensorProd(FL,B,X,aflag,bflag,kflag);
       return X.save2(*this);
    };

    wbsparray& kron(
      const char *F, int L, const wbsparray& B, wbsparray& C,
      char aflag='N', char bflag='N', char kflag='k'
    ) const {
      return tensorProd(F,L,B,C,aflag,bflag,kflag);
    };


    mxArray* toMxSp() const;

    mxArray* toMx(int m=-99) const {
       unsigned r=rank(FL);

       if (m==-99 && r<=2 && isBaseType(typeid(TD))) {
          return toMxSp();
       }

       mxArray *a=mxCreateStruct(1,1);
       add2MxStruct(a,0,m);
       return a;
    };

    mxArray* IDtoMx() const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *a, unsigned i, int ma=-99) const;

    mxArray* mxCreateCell(unsigned m, unsigned n) const;
    void add2MxCell(mxArray *a, unsigned i) const;


    mxArray* toMx0(int m=-99) const {
       mxArray *a=mxCreateStructX(1,1,m);
       add2MxStructX(a,0,m);
       return a;
    };
    mxArray* mxCreateStructX(unsigned m, unsigned n, int ma=-99) const;
    void add2MxStructX(mxArray *a, unsigned i, int ma=-99) const;

    mxArray* toMX() const {
       mxArray *a=mxCreateSTRUCT(1,1);
       return add2MxSTRUCT(a,0);
    };
    mxArray* mxCreateSTRUCT(unsigned m, unsigned n) const;
    mxArray* add2MxSTRUCT(mxArray *a, unsigned i) const;

    void put(const char *vname, unsigned ma=-99, const char *ws="caller"
    ) const { put(0,0,vname,ma,ws); };

    void put0(const char *F, int L,
      const char *vname, const char *ws="caller"
    ) const { put(F,L,vname,-99,ws,'X'); };

    void put(const char *F, int L, const char *vname,
       unsigned ma=-99, const char *ws="caller", char flag=0
     ) const {

       mxArray *a=0; {
          if (!flag) a=toMx(ma); else
          if (flag=='X') a=toMX(); else
          if (flag=='0') a=toMx0(); else
          wblog(FL,"ERR %s() invalid flag=%s",FCT,char2Str(flag).data);
       }
       int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F,L,
          "I/O putting variable '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };


    wbvector<SPIDX_T> SIZE;
    wbMatrix<SPIDX_T> IDX;

    wbvector<TD> D;

    char isref;

  protected:

    unsigned get2DSize(
        const char *F, int L,
        SPIDX_T &d1, SPIDX_T &d2, int m=-1, char fflag=0) const;

    unsigned get2DIndex(
        wbMatrix<SPIDX_T> &IJ, SPIDX_T &d1, SPIDX_T &d2,
        int m=-1) const;

  private:

   TD dotProd_full_diag(const wbsparray &B, char cflag) const;

   wbsparray& timesEl_full_diag(
     const wbsparray &B, wbsparray &X, char cflag
   ) const;

};


template<class TD> inline 
wbsparray<TD> operator+(const wbsparray<TD> &A) { return A; };

template<class TD> inline 
wbsparray<TD> operator-(const wbsparray<TD> &A) {

   wbsparray<TD> X(A); TD *x=X.D.data;

   for (SPIDX_T n=X.D.len, i=0; i<n; ++i) { x[i]=-x[i]; }

   return X;
};


template<class TD>
class iterSparseCMref {

  public:

     iterSparseCMref(const char *F, int L, wbsparray<TD> &A);

     SPIDX_T next(wbsparray<TD> &a);


     SPIDX_T numIter() const {
         if (!IX.len) wblog(FL,
            "ERR %s() got non-initialized/empty object",FCT);
         return IX.len-1;
     };

     bool gotLast() const {
         if (sSPIDX_T(it)<0 || it+1>=IX.len) wblog(FL,
            "ERR %s() not within loop (%d/%d)",FCT,it,IX.len-1);
         return (it+2==IX.len);
     };

     SPIDX_T index() const {
         if (sSPIDX_T(it)<0 || it+1>=IX.len) wblog(FL,
            "ERR %s() not within loop (%d/%d)",FCT,it,IX.len-1);
         return IX[it];
     };

     mxArray* toMx() const;

     SPIDX_T it;

     wbvector<SPIDX_T> I0;

     wbvector<SPIDX_T> IX;

     wbMatrix<SPIDX_T> IDX;

  protected:
  private:

     unsigned r;
     const SPIDX_T *S;
     const TD* data;
};


template <class TD>
mxArray* iterSparseCMref<TD>::toMx() const {

   const char *fields[]={"SIZE","r","it","I0","IX","IDX","data"};
   mxArray *a=mxCreateStructMatrix(1,1,7,fields);

   mxSetFieldByNumber(a,0, 0, wbvector<SPIDX_T>(r,S,'r').toMx());
   mxSetFieldByNumber(a,0, 1, numtoMx(r));
   mxSetFieldByNumber(a,0, 2, numtoMx(it));
   mxSetFieldByNumber(a,0, 3, I0.toMx());
   mxSetFieldByNumber(a,0, 4, IX.toMx());
   mxSetFieldByNumber(a,0, 5, IDX.toMx());
   mxSetFieldByNumber(a,0, 6, wbvector<TD>(IDX.dim1,data,'r').toMx());

   return a;
};


template <class TD>
iterSparseCMref<TD>::iterSparseCMref(
   const char *F, int L, wbsparray<TD> &A) : r(A.SIZE.len) {

   unsigned l=r-1; SPIDX_T D;
   if (r<2 || !(D=A.SIZE.end())) wblog(F_L,
      "ERR %s() got %s (%d)",FCT,A.sizeStr().data,D);

   A.Compress(F_L);
   S=A.SIZE.data; data=A.D.data; it=-1;

   if (A.IDX.dim1!=A.D.len || A.SIZE.len!=A.IDX.dim2)
       A.sperror_this(F_LF);

   if (A.IDX.dim1==0) {
      I0.init(1);
      IX.init(1); IX[0]=D;
      IDX.init(0,r-1); return;
   }
   if (A.IDX.dim1==1) {
      I0.init(2); I0[1]=1;
      IX.init(2); IX[0]=A.IDX(0,l); IX[1]=D;
      IDX.init2ref(1,l,A.IDX.data);
      return;
   }

   SPIDX_T i0=0, i=1, k=1, *I=A.IDX.ref(1,l);
   const int d2=A.IDX.dim2;
   wbvector<SPIDX_T> Dg;

   for (k=1; i<A.IDX.dim1; ++i, I+=d2) { if (I[0]!=I[-d2]) ++k; }
   Dg.init(k); IX.init(k+1);
   
   k=0; I=A.IDX.ref(0,l);
   IX[0]=I[0]; I+=d2;

   for (i=1; i<A.IDX.dim1; ++i, I+=d2) {
      if (I[0]!=I[-d2]) { Dg.data[k++]=i-i0; i0=i; IX[k]=I[0]; }
   }
   Dg.data[k++]=i-i0; IX[k]=D;

   Dg.cumsum0(I0,'N');
   A.IDX.getCols(0,A.IDX.dim2-2,IDX);
};


template <class TD>
SPIDX_T iterSparseCMref<TD>::next(wbsparray<TD> &a) {

   if ((++it)+1>=I0.len) { it=-1; return (it+1); }

   SPIDX_T i0=I0[it], n=I0[it+1]-i0;
   if (sSPIDX_T(n)<=0) wblog(FL,"ERR %s() got I0(%d:%d) !??",FCT,i0,i0+n);

   a.SIZE.init2ref(r-1,S);
   a.D.init2ref(n, data+i0);
   a.IDX.init2ref(n,IDX.dim2,IDX.rec(i0));
   a.isref=1;

   return (it+1);
};


template <class T>
inline int checkContract(const char* F, int L,
   const wbvector<T> &SA, const wbindex &ica,
   const wbvector<T> &SB, const wbindex &icb,
   const wbperm *pfinal=NULL, wbvector<T> *SC=NULL
);

template <class TD>
wbarray<TD>& contract(const char *F, int L,
   const wbsparray<TD>&, const wbindex&,
   const wbarray<TD>&,   const wbindex&, wbarray<TD>&,
   const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
);

template <class TD>
wbarray<TD>& contract(const char *F, int L,
   const wbarray<TD>&,   const wbindex&,
   const wbsparray<TD>&, const wbindex&, wbarray<TD>&,
   const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
);

template <class TD>
TD wbVMatVprod(
   const wbsparray<TD> &v1, const wbsparray<TD> &M, const wbsparray<TD> &v2
);

template <class TD>
wbsparray<TD>& wbMatProd(
  const wbsparray<TD> &A, const wbsparray<TD> &B, wbsparray<TD> &C,
  char aflag='N', char bflag='N', TD afac=TD(1.), TD cfac=TD(0.)
);


template <class TD> inline
wbarray<TD>& contract(const char *F, int L,
   const wbsparray<TD> &A, const char* s1,
   const wbsparray<TD> &B, const char* s2, wbsparray<TD> &C,
   const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
){
   return
   A.contract(F_L, wbindex(s1,1), B,wbindex(s2,1), C, pfinal,afac,cfac);
};

template <class TD> inline
wbarray<TD>& contract(const char *F, int L,
   const wbsparray<TD> &A, const char* s1,
   const wbarray<TD> &B,   const char* s2, wbarray<TD> &C,
   const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
){
   return
   contract(F,L, A,wbindex(s1,1), B,wbindex(s2,1), C, pfinal,afac,cfac);
};

template <class TD> inline
wbarray<TD>& contract(const char *F, int L,
   const wbarray<TD> &A,   const char* s1,
   const wbsparray<TD> &B, const char* s2, wbarray<TD> &C,
   const wbperm &pfinal=wbperm(), const TD& afac=1., const TD& cfac=0.
){
   return
   contract(F,L, A,wbindex(s1,1), B,wbindex(s2,1), C, pfinal,afac,cfac);
};


template<class TD> inline
void wbsparray<TD>::info(const char* F, int L, const char* istr) const {

   SPIDX_T nz=D.len; checkSize(F_L);

   if (!D.len) {
      if (SIZE.len) {
         wblog(F_L,"sparse %3s  %12s  (all zero)",
         istr && istr[0] ? istr:"", SIZE.toStrD().data);
      }
      else {
         wblog(F_L,"sparse %3s  (empty)",
         istr && istr[0] ? istr:"", SIZE.toStrD().data);
      }
      return;
   }

   if (isDiag(F_L)) {
      if (nz>4) {
         wblog(F_L,"sparse %3s  %12s  (diag)",
         istr && istr[0] ? istr:"", SIZE.toStrD().data);
      }
      else if (nz>1) {
         wblog(F_L,"sparse %3s  %8sdiag( [%s] )",
         istr && istr[0] ? istr:"", SIZE.toStrD().data, D.toStr().data);
      }
      else {
         wblog(F_L,"sparse %3s  %5s(scalar)  %-8.4g",
         istr && istr[0] ? istr:"", SIZE.toStrD().data, double(D[0]));
      }
      return;
   }

   if (nz>4) {
      wblog(F_L,"sparse %3s  %12s  (nz=%d @ %-6.3g)",
      istr && istr[0] ? istr:"", SIZE.toStrD().data, nnz(), sparsity());
   }
   else {
      wblog(F_L,"sparse %3s  %12s  {%s}",
      istr && istr[0] ? istr:"", SIZE.toStrD().data, D.toStrf("%.4g",", ").data);
   }
};


template<class TD> inline
void wbsparray<TD>::printdata(const char *istr, const char *fmt0) const {

   SPIDX_T i; const char *fmt; wbstring FMT;

   if (fmt0 && fmt0[0]) fmt=fmt0;
   else { FMT.init2Fmt((TD)0); fmt=FMT.data; }

   checkSize(FLF);

   if (isScalar()) {
      printf("%s = %.8g;\n",istr,D.len ? D[0] : 0);
   }
   else if (isDiag()) {
      printf("%s = [",istr);
      for (i=0; i<D.len; ++i) printf(" %.8g",D[i]);

      printf(" ]'; %s = spdiags(%s,0,%d,%d);\n",
      istr,istr,istr,istr,D.len,D.len);
   }
   else if (SIZE.len==1) {
      printf("%s = [",istr);
      for (i=0; i<IDX.dim1; ++i) {
         printf(" 1 %d %.8g;",IDX(i,1)+1,D[i]);
      }
      printf(" ]; %s = sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d);\n",
      istr,istr,istr,istr,1,SIZE[0]);
   }
   else if (SIZE.len==2) {
      printf("%s = [",istr);
      for (i=0; i<IDX.dim1; ++i) {
         printf(" %d %d %.8g;",IDX(i,0)+1,IDX(i,1)+1,D[i]);
      }
      printf(" ]; %s = sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d);\n",
      istr,istr,istr,istr,SIZE[0],SIZE[1]);
   }
   else {
      wblog(FL,"ERR %s() only applies to rank-2 (%s)",
      FCT,sizeStr().data);
   }

};


template<class TD> inline
unsigned wbsparray<TD>::get2DSize(
   const char *F, int L,
   SPIDX_T &d1, SPIDX_T &d2, int m_,
   char fflag
 ) const {

   const SPIDX_T* const s=SIZE.data;
   unsigned m=0, i=0, r=rank(F_L);

   if (isDiag()) wblog(FL,"ERR %s() got DIAG_REP",FCT);

   if (!r) { d1=d2=0; return 0; }

   if (m_<0) {
      if (r%2 && fflag) wblog(F,L,
         "ERR sparse::%s() even rank or 2nd argument required (%s)",
         FCT,sizeStr().data);
      else {

         SPIDX_T dm=SIZE[0], d=std::sqrt(double(numel()));
         for (m=1; m<SIZE.len; ++m) {
            if (dm>=d) break; else dm*=SIZE[m];
         }
      }
   }
   else { m=m_;
      if (m>r) wblog(F,L,
     "ERR sparse::%s() invalid m=%d (%s)",FCT,m,sizeStr().data);
   }

   if (m>0) { for (d1=s[i++]; i<m; ++i) { d1*=s[i]; }} else d1=1;
   if (m<r) { for (d2=s[i++]; i<r; ++i) { d2*=s[i]; }} else d2=1;

   return m;
};


template<class TD> inline
unsigned wbsparray<TD>::get2DIndex(
   wbMatrix<SPIDX_T> &IJ, SPIDX_T &d1, SPIDX_T &d2, int m_
) const {

   if (isDiag(FL) || !D.len) {
      d1=d2=D.len; IJ.init(d1,2);
      for (SPIDX_T *I=IJ.data, i=0; i<d1; ++i, I+=2) { I[0]=I[1]=i; }
      return 1;
   }
   if (!IDX.dim2) sperror_this(FLF);

   SPIDX_T i,j,k, n=IDX.dim1;
   unsigned r=rank(FL), l=r-1, m=get2DSize(FL,d1,d2,m_);

   IJ.init(IDX.dim1,2);

   if (!d1 || !d2 || !r) {
      if (IJ.data || D.data) wblog(FL,
         "ERR %s() got d1=%d, d2=%d, yet S=[%s] (%d) !??",
         FCT,d1,d2,SIZE.toStr().data,D.len);
      return m;
   }

   SPIDX_T *ij=IJ.data; char e=0;
   const SPIDX_T *idx=IDX.data, *s=SIZE.data;

   if (!m || m==r) {
      for (i=0; i<n; ++i, ij+=2, idx+=r) {
         for (j=l, k=idx[j--]; j<r; --j) { k=k*s[j]+idx[j]; };
         if (m==r)
              { ij[0]=k; if (k>=d1) { e=1; break; }}
         else { ij[1]=k; if (k>=d2) { e=2; break; }}
      }
   }
   else {
      for (i=0; i<n; ++i, ij+=2, idx+=r) {
         j=l;

         for (k=idx[j--]; j>=m;--j) { k=k*s[j]+idx[j]; };
         ij[1]=k;

         for (k=idx[j--]; j<r; --j) { k=k*s[j]+idx[j]; };
         ij[0]=k;

         if (ij[0]>=d1 || ij[1]>=d2) { e=3; break; }
      }
   }

   if (e) wblog(FL,
      "ERR %s() index out of bounds (%d: (%d,%d) %dx%d; %s: %d; %d)",
       FCT, i+1, ij[0]+1, ij[1]+1, d1,d2,sizeStr().data,m, e
   );

   return m;
};


template<class TD> inline
wbsparray<TD>& wbsparray<TD>::initIdentity(
   const wbvector<SPIDX_T> &S0,
   int m,
   sSPIDX_T id,
   char fflag
){
   SPIDX_T d1=1,d2=1;
   wbindex I(S0.len); SIZE=S0; init(S0);

   get2DSize(FL,d1,d2,m);

   if (d1 && d2) { SPIDX_T i=0, n=0;
      if (id==0) {
         n=MIN(d1,d2); init_nnz(n); D.set(1);
         for (; i<n; ++i) ind2sub(i+d1*i,IDX.rec(i));
      }
      else if (id>0) { SPIDX_T k=id; if (k<d2) {
         n=MIN(d1,d2-k); init_nnz(n); D.set(1);
         for (; i<n; ++i) ind2sub(i+d1*(i+k),IDX.rec(i));
      }}
      else if (id<0) { SPIDX_T k=-id; if (k<d1) {
         n=MIN(d1-k,d2); init_nnz(n); D.set(1);
         for (; i<n; ++i) ind2sub((i+k)+d1*i,IDX.rec(i));
      }}
   }

   if (!D.len && fflag)
      wblog(FL,"ERR sparse::%s() to empty (k=%d given %s)",
      FCT, id, sizeStr().data
   );

   return *this;
};


template<class TD> inline
wbsparray<TD>& wbsparray<TD>::initIdentityB3(
   SPIDX_T d1, SPIDX_T d2, SPIDX_T D,
   sSPIDX_T i0,
   char fflag
){
   SPIDX_T s[3]={d1,d2,D}; SIZE.init(3,s);
   return initIdentity(SIZE,2,i0,fflag);
};


template<class TD>
class sparseCollector2D {

 public:

   sparseCollector2D(SPIDX_T m=0, SPIDX_T n=0)
    : M(m), N(n), nnz(0), NIJ(NULL), data(NULL) { SPIDX_T MN=M*N;
      if (MN) {
         WB_NEW(data,MN); MEM_SET<TD>(data,MN);
         WB_NEW(NIJ,MN); memset(NIJ,0,MN*sizeof(SPIDX_T));
      }
   };

  ~sparseCollector2D() { clear(); }

   void clear() {
      if ((data==NULL) ^ (NIJ==NULL)) wblog(FL,
         "ERR %s() got partial initialization (%lX, %lX) !??",FCT,data,NIJ);
      if (data) { WB_DELETE(data); }
      if (NIJ) { WB_DELETE(NIJ); }
      M=N=nnz=0;
   };

   void init(SPIDX_T m=0, SPIDX_T n=0) { SPIDX_T MN=M*N;
      clear(); M=m; N=n; nnz=0;
      if (MN) {
         WB_NEW(data,MN); MEM_SET<TD>(data,MN);
         WB_NEW(NIJ,MN); memset(NIJ,0,MN*sizeof(SPIDX_T));
      }
   };

   sparseCollector2D& add(
      const char *F, int L, SPIDX_T i, SPIDX_T j, const TD& a
   ){
      if (i>=M || j>=N || !data) wblog(F_L,
         "ERR %s() index out of bounds (%d,%d)/(%dx%d) [%lX]",
         FCT,i,j,M,N,data);
      SPIDX_T k=i*N+j;
      data[k]+=a; if ((++NIJ[k])==1) ++nnz;
      return *this;
   };

   void getData(wbMatrix<SPIDX_T> &IJ, wbvector<TD> &X) {
      SPIDX_T i,j, k=0, l=0;
      IJ.init(nnz,2); X.init(nnz);

      for (i=0; i<M; ++i)
      for (j=0; j<N; ++j, ++l) if (data[l]!=0) {
         if (k>=nnz || !NIJ[l]) wblog(FL,"ERR %s() "
            "possibly incorrect nnz (%g/%g; %d) !??",FCT,k,nnz,NIJ[l]);
         X.data[k]=data[l];
         IJ(k,0)=i; IJ(k,1)=j; ++k;
      }
      if (nnz && k<nnz) {
         if (k) { IJ.dim1=k; X.len=k; }
         else {
            wblog(FL,"WRN %s() got all-zero data",FCT);
            IJ.init(0,2); X.init();
         }
      }
   };

   mxArray* toMx() const {
      const char *fields[]={"M","N","nnz","data","NIJ"};
      mxArray *S=mxCreateStructMatrix(1,1,5,fields);

      mxSetFieldByNumber(S,0,0, numtoMx(M));
      mxSetFieldByNumber(S,0,1, numtoMx(N));
      mxSetFieldByNumber(S,0,2, numtoMx(nnz));
      mxSetFieldByNumber(S,0,3, wbMatrix<TD>(M,N,data,'r').toMx());
      mxSetFieldByNumber(S,0,4, wbMatrix<SPIDX_T>(M,N,NIJ ,'r').toMx());

      return S;
   };

 protected:
 private:

   SPIDX_T M,N,nnz,*NIJ;
   TD* data;
};



#endif

