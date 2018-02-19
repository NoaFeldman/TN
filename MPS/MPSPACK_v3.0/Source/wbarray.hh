#ifndef __WB_ARRAY_COL_MAJOR_HH__
#define __WB_ARRAY_COL_MAJOR_HH__


//===================================================================//
// wbarray - N-D Array class                                         //
//                                                                   //
//    generic N-dimensional class that allows permutations           //
//    and reshaping; index order is column-major.                    //
//                                                                   //
// NOTES                                                             //
//                                                                   //
// Wb,Sep14,09 ;  Wb,Aug11,05                                        //
//===================================================================//

// Warning: A = Str2Idx("3,5"); works!!!
//    It is interpreted using the constructor
//    wbarray(const wbvector<size_t> &d);
// The result is an initialized 2D array of dimensions 3x5 !!!




template <class T>
class wbarray : public wb_sptr<T> {

  public:

    wbarray(size_t d1=0) { init(d1); };
    wbarray(size_t d1, size_t d2) { init(d1,d2); };
    wbarray(size_t d1, size_t d2, size_t d3) { init(d1,d2,d3); };
    wbarray(size_t d1, size_t d2, size_t d3, size_t d4) {
       init(d1,d2,d3,d4); };

    wbarray(const wbarray &a, const char ref=0) { init(a,ref); };

    wbarray(char *sidx   ) { init(Str2Idx(sidx)); };
    wbarray(const UVEC &d) { init(d); };

    template<class IT>
    wbarray(const wbvector<IT> &d) { init(d); };

    wbarray(const mxArray* a, char ref=1, char vec=0) {
       init(a,ref,vec); };

    wbarray(const char *F, int L,
       const mxArray* a, char ref=1, char vec=0) { init(F,L,a,ref,vec); };

    wbarray(const wbMatrix<T> &a) { init(a); }

   ~wbarray() {
       DELETE_DATA();
    };

    wbarray& init(size_t d1=0) {
       if (d1) NEW(wbvector<size_t>(1,&d1));
       else if (data) DELETE_DATA();
       return *this;
    };

    wbarray& init(size_t d1, size_t d2) {
       size_t s[2]= {d1,d2}; NEW(wbvector<size_t>(2,s));
       return *this;
    };

    wbarray& init(size_t d1, size_t d2, size_t d3) {
       size_t s[3]= {d1,d2,d3}; NEW(wbvector<size_t>(3,s));
       return *this;
    };

    wbarray& init(size_t d1, size_t d2, size_t d3, size_t d4) {
       size_t s[4]= {d1,d2,d3,d4}; NEW(wbvector<size_t>(4,s));
       return *this;
    };

    wbarray& initp(size_t d1, const T *d) {
       NEW(wbvector<size_t>(1,&d1),d);
       return *this;
    };

    wbarray& initp(size_t d1, size_t d2, const T *d) {
       size_t s[]= {d1,d2};
       NEW(wbvector<size_t>(2,s),d);
       return *this;
    };

    wbarray& initp(size_t d1, size_t d2, size_t d3, const T *d) {
       size_t s[]= {d1,d2,d3};
       NEW(wbvector<size_t>(3,s),d);
       return *this;
    };

    wbarray& init2ref(const wbarray &a) { return init(a,'r'); }

    wbarray& init(const wbarray &a, char ref=0) {
       if (!ref) {
          NEW(a.SIZE, a.data);
       }
       else {
          SIZE=a.SIZE; 
          if (a.sptr_flags)
               { a.add_dref(*this); }
          else { wb_sptr<T>::init2ref(a.data); }
       }
       return *this;
    };


    wbarray& initDef(size_t d1) {
        size_t s[]={d1};
        NEW_DATA_BARE(d1);
        SIZE.init(1,s); return *this;
    };
    wbarray& initDef(size_t d1, size_t d2) {
        size_t s[]={d1,d2};
        NEW_DATA_BARE(d1*d2);
        SIZE.init(2,s); return *this;
    };
    wbarray& initDef(size_t d1, size_t d2, size_t d3) {
        size_t s[]={d1,d2,d3};
        NEW_DATA_BARE(d1*d2*d3);
        SIZE.init(3,s); return *this;
    };

    wbarray& init(const wbMatrix<T> &a, char ref=0) {
       if (ref) {
          size_t s[]={a.dim2,a.dim1}; SIZE.init(2,s);

          if (data!=a.data) {
             if (a.isEmpty()) init();
             else wb_sptr<T>::init2ref(a.data);
          }
       }
       else { 
          wbarray X; save2(X);
          X.permute(*this,"2,1");
       }

       return *this;
    };

    wbindex& ind2sub(size_t k, wbindex &I) const {
       const size_t *s=SIZE.data, n=SIZE.len; 
       size_t i,r, e=0, *idx;

       I.init(n); idx=I.data;
       for (i=0; i<n; ++i) {
           if (!s[i]) wblog(FL,
              "ERR wbarray::%s() got size %s",FCT,sizeStr().data);

           r=k/s[i]; idx[i]=k-r*s[i]; k=r;
           if (idx[i]>=s[i]) e++;
       }

       if (e || k) wblog(FL,
          "ERR wbarray::%s() out of bounds (%d => [%s; %s])",
           FCT, (I+1).toStr().data, sizeStr().data
       ); 

       return I;
    };

    void adjustMMat(
       const char* F, int L,
       char flag, const T fac, wbarray &A
    ) const {
       
       if (!isMatrix()) wblog(F,L,
       "ERR %d - rank-2 object required (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
       "ERR %d - invalid flag %c<%d>",FCT,flag,flag);

       size_t s=numel();

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';

       if (fac!=T(1) || flag!='N') { size_t i;
          if (flag!='N') permute(A,"2 1"); else A=(*this);
          if (flag=='C') for (i=0; i<s; i++) A[i]=CONJ(A[i]);
          if (fac!=T(1)) for (i=0; i<s; i++) A[i]*=fac;
       }
       else { A.init2ref(*this); }
    }

    void adjustDMat(
       const char* F, int L,
       char flag, const T fac, wbarray &A
    ) const {
       
       if (SIZE.len!=1) wblog(F,L,
       "ERR %d - 1D representation required (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
       "ERR %d - invalid flag %c<%d>",FCT,flag,flag);

       size_t n=SIZE[0];

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';
       
       if (fac!=T(1) || flag=='C') { size_t i; A=(*this);
          if (flag=='C') for (i=0; i<n; i++) A[i]=CONJ(A[i]);
          if (fac!=T(1))    for (i=0; i<n; i++) A[i]*=fac;
       }
       else { A.init2ref(*this); }
    }

    void adjustVec(
       const char* F, int L, char &flag, wbarray &A, char pos
     ) const { A=(*this); A.adjustVec(F,L,flag,pos); };
    
    void adjustVec(const char* F, int L, char &flag, char pos) {
       if (SIZE.len!=1) wblog(F,L,
          "ERR %s() vector expected (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
          "ERR %s() invalid flag %c<%d>",FCT,flag,flag);

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';
       if (flag=='C') {
          for (size_t s=numel(), i=0; i<s; ++i)
          data[i]=CONJ(data[i]);
       }
       flag='N';

       SIZE.Resize(2);
       if (pos==1) { SIZE[1]=SIZE[0]; SIZE[0]=1; } else
       if (pos==2) { SIZE[1]=1; } else
       wblog(F,L,"ERR %s() invalid pos=%d",FCT,pos);
    };

    void adjustCMat(const char* F, int L,
       T cfac, size_t s1, size_t s2, char cforce
    ){
       if (cfac!=T(0)) {
           if (data==NULL) {
              if (cforce) wblog(F,L,
             "WRN C = A*B + c*[] with c=%s !??", toStr(cfac).data);
           }
           else if (!isMatrix() || SIZE[0]!=s1 || SIZE[1]!=s2) {
              wblog(F,L,"ERR %s() dimension mismatch: C=(%s) =? (%d,%d).",
              FCT, sizeStr().data, s1, s2); return;
           }
           else {
               if (cfac!=T(1)) {
                  size_t i, s=SIZE.prod(0);
                  for (i=0; i<s; i++) data[i]*=cfac;
               }
               return;
           }
       }
       init();
    };

    void set (const wbarray &a, const T fac=1) {
       if (fac==T(0)) NEW(a.SIZE);
       else {
          NEW(a.SIZE, a.data);
          if (fac!=T(1)) {
             if (fac==T(-1)) {
                for (size_t n=numel(), i=0; i<n; i++)
                data[i]=-data[i];
             }
             else {
                for (size_t n=numel(), i=0; i<n; i++)
                data[i]*=fac;
             }
          }
       }
    };

    void setBlock(
        const wbvector< WBINDEX > &D,
        const wbindex &IB,
        const wbarray &A 
    );

    void addBlock(
        const WBINDEX &I,
        const wbarray &A, const char dflag=0
    );

    wbarray& addBlock(
       size_t i, size_t j,
       size_t n, size_t m,
       wbarray &a
    ) const;

    wbarray& getBlock(
       size_t i, size_t j,
       size_t n, size_t m,
       wbarray &a
    ) const;

    wbarray& BlockDiag(const wbvector< wbarray > &D);

    wbarray& blockTrace(
       size_t D,
       wbarray &a, size_t D2=-1
    ) const;

    void copyStride(T* dd, size_t stride) const;

    template<class T0>
    void initT(const wbarray<T0> &a) { NEW(a.SIZE);
       for (size_t n=SIZE.prod(), i=0; i<n; i++)
       data[i]=T(a.data[i]);
    };

    template <class T0>
    void initT(const char *F, int L, const wbarray<T0> &a) {
       NEW(a.SIZE);
       for (size_t n=SIZE.prod(), i=0; i<n; i++) {
          data[i]=T(a.data[i]);
          if (T0(data[i])!=a.data[i]) wblog(F,L,
             "ERR type conversion changes value (%g,%g)",
              double(a.data[i]), double(T0(data[i]))
          );
       }
    };

    wbarray& init(char *s, T* d=NULL) {
       init(Str2Idx(s), d); return *this;
    };

    wbarray& init(const WBINDEX &S, T* d=NULL) {
       NEW(S,d); return *this;
    };

    template <class IT>
    wbarray& init(const wbvector<IT> &S_, T* d=NULL) {
       WBINDEX S(S_.len);
       for (size_t i=0; i<S.len; ++i) { S.data[i]=size_t(S_.data[i]); }
       NEW(S,d); return *this;
    };


    wbarray& init(const char *F, int L,
    const mxArray *a, char ref=1, char vec=0);

    wbarray& init(const char *F, int L,
    const mxArray *a, wbperm &P0);

    wbarray& init(const mxArray *a, char ref=0, char vec=0) {
    return init(FL,a,ref,vec); }
    wbarray& init(const mxArray *a, wbperm &P0) {
    return init(FL,a,P0); }




    wbarray& init2ref(const T *x, const WBINDEX &S) {
       SIZE=S;
       wb_sptr<T>::init2ref(x);
       return *this;
    };

    wbarray& init2ref(size_t d1, const T *x) {
       SIZE.init(1,&d1);
       wb_sptr<T>::init2ref(x);
       return *this;
    };

    wbarray& init2ref(size_t d1, size_t d2, const T *x) {
       size_t s[]= {d1,d2}; SIZE.init(2,s);
       wb_sptr<T>::init2ref(x);
       return *this;
    };

    wbarray& init2ref(size_t d1, size_t d2, size_t d3, const T *x) {
       size_t s[]= {d1,d2,d3}; SIZE.init(3,s); 
       wb_sptr<T>::init2ref(x);
       return *this;
    };

    void init2ref(const T *d0) {
       if (d0) {
          if (d0==data) return;

          wb_sptr<T>::init2ref(d0);
       }
       else {
          if (SIZE.prod(0)) wblog(FL,
          "ERR data/size inconsistency (%d)",SIZE.prod(0));
       }
    };

    wbarray& Instantiate(const char *F=NULL, int L=0) {
       wb_sptr<T>::Instantiate(F,L,SIZE.prod(0));
       return *this;
    };

    wbarray& unRef() { return Instantiate(); };

    void init2Vec(const T *d0, size_t len, char ref=0) {
       SIZE.init(1,&len);
       if (ref)
            wb_sptr<T>::init2ref(d0);
       else wb_sptr<T>::new_dref(len,d0);
    };

    void init2ref(const wbarray &A, size_t s1, size_t s2) {
        size_t i,s=1; int k=-1;

           for (i=0; i<SIZE.len; i++) { if (s==s1) { k=i; break; }
              s*=SIZE[i];
           }
           for (s=1; i<SIZE.len; i++) s*=SIZE[i];
           if (k<0 || s!=s2) wblog(FL,
             "ERR cannot reshape %s array into %dx%d",
              sizeStr().data, s1, s2
           );

        SIZE.init(2); SIZE[0]=s1; SIZE[1]=s2;
        set_dref(A);
    };

    void initTst();

    void init2Vec(const wbvector<T> &a) {
       WBINDEX S(1); S[0]=a.len;
       NEW(S, a.data);
    };

    wbarray& initIdentity(size_t d, char dflag=0, T dval=1);
    wbarray& initIdentity(const WBINDEX &S, char dflag=0);
    wbarray& initIdentityB(
        size_t d1, size_t d2, size_t k=0);
    wbarray& initIdentityB3(
        size_t d1, size_t d2, size_t D, size_t i0, T one=1);

    void Expand2Projector(const T eps=0);
    void ExpandDiagonal(unsigned i1, unsigned i2);
    void ExpandDiagonal();

    void Diag2Vec();

    wbarray& initDiag(unsigned n, const T *d){
       init(n,n);
          for (size_t i=0; i<n; ++i) data[i+i*n]=d[i];
       return *this;
    };

    wbarray& initDiag(const wbvector<T> &d){
       return initDiag(d.len, d.data); };

    void initVec(size_t d, const T val=0) {
       init(d); if (val) set(val);
    };


    void resize(const wbvector<size_t> &S, wbarray &B) const;

    void Resize(C_UVEC &S, T* i=NULL) { RESIZE(S,i); };
    void Resize(size_t s1, size_t s2) {
       size_t s[2]={s1,s2};
       if (SIZE.len!=2) wblog(FL,"ERR rank-2 required (%d)",SIZE.len);

       wbarray X; this->save2(X);
       X.resize(wbvector<size_t>(2,s), *this);
    };

    wbarray& Enlarge(
       size_t dim,
       size_t D,
       size_t istart=0
    );

    wbarray& save2(wbarray&);
    void save2(wbMatrix<T>&);

    wbarray& Append2(
       const char *F, int L, wbarray &B, unsigned dim);
    wbarray& Append2(wbarray &B, unsigned dim){
       return Append2(0,0,B,dim);
    }

    void Append2( 
       const char *F, int L,
       wbarray &B, unsigned dim,
       wbarray<double> &W, double eps
    );
    void Append2(wbarray &B, unsigned dim, wbarray &W, double eps){
       return Append2(0,0,B,dim,W,eps);
    };

    void rand(const UVEC &d) { init(d); setRand(); };
    void setRand(const char pnflag=0);

    wbarray& set(const T& c) {
        size_t s=SIZE.prod();
        if (c!=T(0)) { for (size_t i=0; i<s; i++) data[i]=c; }
        else { MEM_SET<T>(data,s); }
        return *this;
    };

    wbarray& operator= (const wbarray &a) {
        if (this==&a) return *this;

        NEW(a.SIZE, a.data);
        return *this;
    };

    template <class T0>
    wbarray<T>& operator= (const wbarray<T0> &a) {
        if ((void*)this==(void*)&a) return *this;

        NEW(a.SIZE);
        for (size_t n=a.numel(), i=0; i<n; i++) data[i]=T(a.data[i]);

        return *this;
    };

    void SetRef(const INDEX_T *I, size_t len, const T* d0);
    void SetRef(size_t k, const T* d0) { SetRef(&k,1,d0); };

    const T* ref(const INDEX_T *I, size_t n) const {
       if (!n) return NULL;
       return data+serial_index(I,n);
    };
    T* ref(const INDEX_T *I, size_t n) {
       if (!n) return NULL;
       return data+serial_index(I,n);
    };

    T* col(size_t k) {
       if (SIZE.len!=2 || k>=SIZE[1]) wblog(FL,
          "ERR %s() index out of bounds (%s; %d)",FCT,sizeStr().data,k); 
       return (data+k*SIZE[0]);
    };
    const T* col(size_t k) const {
       if (SIZE.len!=2 || k>=SIZE[1]) wblog(FL,
          "ERR %s() index out of bounds (%s; %d)",FCT,sizeStr().data,k); 
       return (data+k*SIZE[0]);
    };

    wbvector<T>& getCol(size_t k, wbvector<T> &v) const {
       return v.init(SIZE[0],col(k));
    };

    T* ref(const WBINDEX &I) { return ref(I.data,I.len); };
    const T* ref(const WBINDEX &I) const {
        return ref(I.data,I.len);
    };

    const T* ref(INDEX_T i) const { return ref(&i,1); }
    T* ref(INDEX_T i) { return ref(&i,1); }

    const T* ref(INDEX_T i1, INDEX_T i2) const {;
       INDEX_T I[]={i1,i2}; return ref(I,2);
    };
    T* ref(INDEX_T i1, INDEX_T i2) {;
       INDEX_T I[]={i1,i2}; return ref(I,2);
    };

    const T* ref(INDEX_T i1, INDEX_T i2, INDEX_T i3) const {
       INDEX_T I[]={i1,i2,i3}; return ref(I,3);
    };
    T* ref(INDEX_T i1, INDEX_T i2, INDEX_T i3) {
       INDEX_T I[]={i1,i2,i3}; return ref(I,3);
    };

    const T& operator[] (size_t i) const { return data[i]; }
    T& operator[] (size_t i) { return data[i]; }

    const T& operator() (size_t i) const {
        if (SIZE.len!=1) wblog(FL,
           "ERR %s requires vector type (%s)",FCT,sizeStr().data);
        return data[i];
    };
    T& operator() (size_t i) {
        if (SIZE.len!=1) wblog(FL,
           "ERR %s requires vector type (%s)",FCT,sizeStr().data);
        return data[i];
    };

    const T& operator() (size_t i, size_t j) const {
        if (SIZE.len!=2) wblog(FL,
           "ERR %s(i,j) requires matrix type (%s)",FCT,sizeStr().data);
        return data[ i + j*SIZE[0] ];
    };
    T& operator() (size_t i, size_t j) {
        if (SIZE.len!=2) wblog(FL,
           "ERR %s(i,j) requires matrix type (%s)",FCT,sizeStr().data);
        return data[ i + j*SIZE[0] ];
    };

    const T& operator() (size_t i, size_t j, size_t k) const {
        if (SIZE.len!=3) wblog(FL,
           "ERR %s(i,j,k) requires 3D-object (%s)",FCT,sizeStr().data);
        return data[ i + SIZE[0] * (j + SIZE[1]*k) ];
    };
    T& operator() (size_t i, size_t j, size_t k) {
        if (SIZE.len!=3) wblog(FL,
           "ERR %s(i,j,k) requires 3D-object (%s)",FCT,sizeStr().data);
        return data[ i + SIZE[0] * (j + SIZE[1]*k) ];
    };


    const T& operator() (const WBINDEX &I) const {
       if (I.len!=SIZE.len) wblog(FL,
          "ERR %s() invalid index [%s]",FCT,I.toStr().data);
       return data[serial_index(I.data,I.len)];
    };
    T& operator() (const WBINDEX &I) {
       if (I.len!=SIZE.len) wblog(FL,
       "ERR invalid index [%s] having %s",I.toStr().data,sizeStr().data);
       return data[serial_index(I.data,I.len)];
    };

    const T& element(const WBINDEX &I) const {
       if (I.len!=SIZE.len) wblog(FL,"ERR %s() invalid index [%s] "
          "having %s",FCT,I.toStr().data,sizeStr().data);
       for (size_t i=0; i<I.len; i++) if (I.data[i]>=SIZE.data[i])
           wblog(FL,"ERR %s() index out of bounds (%s; %s)",
           FCT,I.toStr().data,sizeStr().data); 
       return data[serial_index(I.data,I.len)];
    };
    T& element(const WBINDEX &I) {
       if (I.len!=SIZE.len) wblog(FL,"ERR %s() invalid index [%s] "
          "having %s",FCT,I.toStr().data,sizeStr().data);
       for (size_t i=0; i<I.len; i++) if (I.data[i]>=SIZE.data[i])
           wblog(FL,"ERR %s() index out of bounds (%s; %s)",
           FCT,I.toStr().data,sizeStr().data); 
       return data[serial_index(I.data,I.len)];
    };

    size_t numel() const { 
       if (SIZE.len) {
          const size_t *s=SIZE.data;
          size_t n=s[0], r=SIZE.len, i=1; for (; i<r; i++) n*=s[i];
          return n;
       }
       return 0;
    };

    wbstring sizeStr() const {
       if (SIZE.len) return SIZE.toStrf("","x");
       else return wbstring("[]");
    };

    void squeeze();
    void skipSingletons();
    void addSingletons(unsigned r, const char *F=NULL, int L=0);
    void prependSingletons(unsigned r);

    size_t size(unsigned i) const { return SIZE[i]; };
    size_t dim (unsigned i) const {
       if (i==0 || i>SIZE.len) wblog(FL,
          "ERR index out of bounds (%d/%d)", i-1, SIZE.len);
       return SIZE[i-1];
    };

    size_t dim0(unsigned i) const {
       if (i>=SIZE.len) wblog(FL,
          "ERR index out of bounds (%d/%d)", i, SIZE.len);
       return SIZE[i];
    };

    void getMatSize(
       const char *F, int L, size_t &dim1, size_t &dim2
    ) const {

       size_t i, r=SIZE.len, r2=r/2;
       const size_t *const &s=SIZE.data;

       if (r%2) wblog(F_L,
          "ERR %s() even rank object required (%s)",FCT,sizeStr().data);
       if (r) {
          dim1=s[0]; dim2=s[r2];
          for (i=1; i<r2; i++) { dim1*=s[i]; dim2*=s[i+r2]; }
       } else { dim1=dim2=0; }
    };

    unsigned rank() const { return SIZE.len; };

    bool isRank(unsigned r) const {
       if (SIZE.len!=r) {
          if (SIZE.len<r) return 0;
          for (unsigned i=r; i<SIZE.len; i++) if (SIZE[i]!=1) return 0;
       }
       return 1;
    };

    bool isOpS(size_t *n=NULL) const;
    bool isOpS(unsigned &n) const { bool q;
       size_t N; q=isOpS(&N); n=N;
       if (size_t(n)!=N) wblog(FL,
          "ERR %s() unsigned out of bounds (%d/%ld)",FCT,n,N);
       return q;
    };


    wbarray& SkipTiny_float(const T eps=1E-14) { return 0; };

    double SkipTiny(const double eps_=1E-14) {
       T zero=T(0), x2=zero, eps=T(eps_);
       if (eps!=zero) { size_t i=0, n=numel(); T a;
          for (; i<n; ++i) { if (data[i]!=zero) {
             a=ABS(data[i]); if (a<eps) { x2+=NORM2(a); data[i]=zero; }
          }}
       }
       return std::sqrt(double(x2));
    };

    bool isMatrix() const { return SIZE.len==2; }
    bool isSMatrix(unsigned d=-1) const {
       if (SIZE.len!=2 || SIZE[0]!=SIZE[1]) return 0;
       if (int(d)>=0 && d!=SIZE[0]) return 0;
       return 1;
    };

    char isEmpty() const { return (data==NULL); };
    bool isComplex() const;
    bool isZero(double eps=0., char flag=0) const;

    wbarray& Symmetrize(
       const char *F, int L,
       double *delta=NULL, char cflag=1, char tflag=0, char fflag=0
    );

    wbarray& swapRows(size_t i1, size_t i2);
    wbarray& swapCols(size_t j1, size_t j2);

    wbarray& setCol(
       size_t k, const T* d, const char *F=NULL, int L=0);
    wbarray& setCol(
       size_t k, size_t k0, const char *F=NULL, int L=0);
    wbarray& setRow(
       size_t k, const T* d, const char *F=NULL, int L=0);

    wbvector<T>& colNorm2(wbvector<T> &a) const;

    T colNorm2(size_t k) const;

    unsigned QRdecomp(
       const char *F, int L, wbarray<T> &Q, wbarray<T> &R,
       T eps=1E-15);

    int Householder(
       const char *F, int L,
       size_t k, T *u,
       WBPERM *P=NULL,
       T eps=1E-15);

    wbarray& Hessenberg(size_t k, wbvector<T> &u, wbarray<T> *U=NULL);

    wbarray& ColProject(
       size_t k1, size_t k2, char nflag=0, char tnorm=0);

    wbarray& ColPermute(const wbperm &P, char iflag=0);

    bool isOrthogonalCol(size_t k, char tnorm=0) const;

    wbarray& SignCol(size_t k, const T *d0, char tnorm=0);
    wbarray& FlipSignCol(size_t k);
    wbarray& SignConventionCol(
       size_t k=-1,
       double eps=1E-12
    );

    T NormalizeCol(size_t k, char tnorm=0, char qflag=0);
    wbarray& NormalizeCols(
       const char *F, int L,
       double *amin=NULL, double *amax=NULL, char tnorm=0
    );

    wbarray& OrthoNormalizeCols(
       const char *F=NULL, int L=0, char tnorm=0,
       char qxflag=0, double eps=1E-14, unsigned np=1
    );

    wbarray& balanceOp(
       const char *F, int L, T &xref, double &xscale
    );

    bool isProptoId(T &x, T eps=1E-14) const;

    bool isDiagMatrix(double eps=1E-14) const {
    return isDiag_aux(eps,"isDiag"); }

    bool isIdentityMatrix(double eps=1E-14) const { T one=1;
       return isProptoId(one,eps);
    }

    bool isDiagMatrix(double *eps) const {
    return isDiag_aux(eps,"isDiag"); }

    bool isIdentityMatrix(double *eps) const {
    return isDiag_aux(eps,"isIdty"); }

    size_t nnz(const T eps=0, WBINDEX *I=NULL) const;

    bool isHConj(
      const wbarray &B, double eps=1E-12, double* xref=NULL
    ) const { return isSym_aux(FLF, B, eps, xref, 's'); };

    bool isHConj(
      double eps=1E-12, double* xref=NULL
    ) const { return isSym_aux(FLF,*this,eps,xref,'s'); };

    bool isAHerm(
      const wbarray &B, double eps=1E-12, double* xref=NULL
    ) const { return isSym_aux(FLF, B, eps, xref, 'a'); };

    bool isAHerm(
      double eps=1E-12, double* xref=NULL
    ) const { return isSym_aux(FLF,*this,eps,xref,'a'); };

    char hasGroupSize(size_t s1, size_t s2) const;

    double aMax(size_t *i=NULL) const;
    double maxDiff (const wbarray &B) const;

    T froNorm2(const wbarray &B) const;
    T norm2() const;
    T norm() const { return SQRT(norm2()); };

    T normCol2(size_t k) const;

    T scalarProd(const wbarray &B) const;
    T sum() const;
    T trace() const;

    wbarray& sum(const WBINDEX &I, wbarray&) const;
    wbarray& sum(const char *sidx, wbarray &A) const {
       sum(Str2Idx(sidx,1), A);
       return A;
    };

    bool hasSameSize(const size_t *s, size_t n) const;

    bool hasSameSize(const wbarray &) const;
    char sameUptoFac(
      const wbarray &B, T *fac=NULL, double eps=1E-12
    ) const;

    bool operator== (const wbarray &B) const {
        if (this==&B) return 1;
        if (SIZE!=B.SIZE) return 0;
        if (data==B.data) return 1;
        return (!memcmp(data, B.data, numel()*sizeof(T)));
    };

    bool operator!= (const wbarray &B) const {
       return !(*this==B);
    };

    double normDiff2(const wbarray &B) const;
    double normDiff2(const wbvector<T> &B) const;

    double normDiff(const wbarray &B) const {
       return std::sqrt(normDiff2(B)); };

    double normDiff(const wbvector<T> &B) const {
       return std::sqrt(normDiff2(B)); };

    void TimesEl(const wbarray &);
    wbarray& timesEl(const wbarray &B, wbarray &C, T afac=1, T cfac=0) const;

    T sumTimesEl(const wbarray &) const;
    T weightedAvg(const wbarray &) const;

    void TensorProd(const wbarray &B,
       const char aflag='N', const char bflag='N', const char kflag=0
    ){
       wbarray X(*this);
       X.tensorProd(B,*this,aflag,bflag,kflag);
    };

    void tensorProd(
       const wbarray &, wbarray &,
       const char aflag='N', const char bflag='N', const char kflag=0
    ) const;

    void kron(
       const wbarray &B, wbarray &X, char aflag='N', char bflag='N'
    ) const { tensorProd(B,X,aflag,bflag,'k'); };

    void Kron(const wbarray &B, char aflag='N', char bflag='N') {
       wbarray X(*this);
       X.tensorProd(B,*this,aflag,bflag,'k');
    };

    wbarray& Cat(unsigned dim, const wbvector< wbarray* > &ap);
    wbarray& Cat(unsigned dim, const wbvector< wbarray  > &ap);

    wbarray& BlockCat(
       const wbarray< wbarray > &aa,
       WBINDEX *D1=NULL, WBINDEX *D2=NULL
    ){
       wbarray< const wbarray* > ap(aa.SIZE);
       for (size_t n=aa.numel(), i=0; i<n; i++) ap[i]=&aa[i];
       return BlockCat(ap,D1,D2);
    };

    wbarray& BlockCat(
       const wbarray< const wbarray* > &ap,
       WBINDEX *D1=NULL, WBINDEX *D2=NULL
    );

    void operator+= (T x) {
       if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
       for (size_t s=numel(), i=0; i<s; i++) data[i]+=x;
    };

    void operator-= (T x) {
       if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
       for (size_t s=numel(), i=0; i<s; i++) data[i]-=x;
    };

    void operator*=(const T a) {
       if (a==T(+1)) return;
       size_t i, s=numel();
       if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
       if (a==T(-1)) for (i=0; i<s; i++) data[i] = -data[i]; else
       if (a==T( 0)) MEM_SET<T>(data,s);
       else     { for (i=0; i<s; i++) data[i] *= a; }
    };

    void operator+= (const wbarray &);
    void operator-= (const wbarray &);

    wbarray& plus (
    const wbarray &B, wbarray &C, T bfac=1, char iflag=0) const;

    wbarray& minus(
    const wbarray &B, wbarray &C, T bfac=1, char iflag=0) const;

    wbarray& Plus (const wbarray &B, T bfac=1, char iflag=0);
    wbarray& Minus(const wbarray &B, T bfac=1, char iflag=0);

    bool equal  (const wbarray&, const T& eps) const;
    bool unequal(const wbarray&, const T& eps) const;
    bool equal  (const wbarray&, const double& eps, double& maxdiff) const;
    bool unequal(const wbarray&, const double& eps, double& maxdiff) const;

    bool allEqual(const T &x) {
       for (size_t n=numel(), i=0; i<n; ++i) { if (data[i]!=x) return 0; }
       return 1;
    };
    bool allLE(const T &x) {
       for (size_t n=numel(), i=0; i<n; ++i) { if (data[i]>x) return 0; }
       return 1;
    };

    bool isVector() const {
        if (SIZE.len==1 || (SIZE.len==2 && (SIZE[0]==1 || SIZE[1]==1)))
        return 1; else return 0;
    };
    bool isScalar() const { return SIZE.allEqual(1); };

    bool requiresDataPerm(const wbperm &P) const;

    void markUnequalZero(const T m=1, const T eps=0) {
        size_t s=numel(), i=0;
        if (eps==0) {
           for (; i<s; i++) data[i] = ((data[i]!=0) ? m : 0);
        }
        else {
           for (; i<s; i++) data[i] = ((ABS(data[i])>eps) ? m : 0);
        }
    };

    void toMatrixRef(
        wbarray &A,
        const UVEC &I,
        int pos,
        wbperm &P
    ) const;

    bool toMatrixRef(const char *F, int L,
        wbarray &A,
        const UVEC &I,
        int pos,
        char &tflag
    ) const;

    void toMatrixRef(
      wbarray &A, unsigned K, const wbperm &P0=wbperm()
    ) const;

    void toMatrixRef(
      wbarray &A, unsigned i, int pos, wbperm &P
    ) const { return toMatrixRef( A, UVEC(1,&i), pos, P); };

    void toMatrixRef(
      wbarray &A, const char* sidx, int pos, wbperm &P
    ) const { return toMatrixRef( A, Str2Idx(sidx,1), pos, P); };

    wbarray& Reshape(const wbvector<size_t> &);

    wbarray& Reshape(size_t s1, size_t s2) {
       WBINDEX S(2); S.data[0]=s1; S.data[1]=s2;
       return Reshape(S);
    };

    wbarray& Reshape(size_t s1, size_t s2, size_t s3) {
       WBINDEX S(3); S.data[0]=s1; S.data[1]=s2;
       S.data[2]=s3; return Reshape(S);
    };

    void GroupIndizes(size_t K);

    wbarray& groupIndizes_DREF(
       wbarray &A, unsigned k1, unsigned k2
    ) const {
       size_t i, s1=1, s2=1;

       if (k1+k2!=SIZE.len) wblog(FL,
       "ERR invalid index group %d+%d=%d ???",k1,k2,SIZE.len);

       k2=SIZE.len;

       for (i=0; i<k1; ++i) { s1*=SIZE[i]; }
       for (   ; i<k2; ++i) { s2*=SIZE[i]; }

       A.SIZE.init(2); A.SIZE[0]=s1; A.SIZE[1]=s2;
       A.wb_sptr<T>::init2ref(data);

       return A;
    };

    void GroupIndizes(unsigned k1, unsigned k2) {
       GroupIndizes_DREF(*this,k1,k2);
    };

    void groupIndizes_P(
       const WBINDEX&, int, wbperm&,
       size_t* =NULL, size_t* =NULL
    ) const;

    wbarray& transpose(const char *F, int L, wbarray &A) const;
    wbarray& Transpose(const char *F=0, int L=0) {
        wbarray A0(*this);
        A0.transpose(F_L,*this);
        return *this;
    };

    wbarray& MatPermute(const wbperm &P, char iflag=0);

    wbarray& Permute(const char* s, char iflag=0);
    wbarray& Permute(const WBPERM&, char iflag=0);
    wbarray& permute(wbarray&, const char* s, char iflag=0) const;
    wbarray& permute(wbarray&, const WBPERM&, char iflag=0) const;


    wbarray& RPermute(C_UVEC &P, char iflag=0);

    wbarray& select0(
       const WBPERM &, unsigned dim, wbarray&) const;
    wbarray& Select0(const UVEC &I, unsigned dim) {
       wbarray A(*this); A.select0(I,dim,*this);
       return *this;
    };


    wbarray& select0(
      size_t i1, size_t i2, unsigned dim, wbarray &Q
    ) const {
       if (dim>=SIZE.len) wblog(FL,
          "ERR dim out of bounds (%d/%d)",dim+1,SIZE.len);
       if (i1>=SIZE[dim] || i2>=SIZE[dim]) wblog(FL,
          "ERR index out of bounds (%d,%d/%d)",i1,i2,SIZE[dim]);

       if (i2<i1) { 
          Q.init(); Q.SIZE=SIZE; Q.SIZE[dim]=0;
          return Q;
       }

       wbindex I; I.Index(i1,i2);
       return select0(I,dim,Q);
    };

    wbarray& selectSqueeze(
       wbarray &A, unsigned p,
       unsigned dim) const;

    T min() const;
    T max() const;

    T aMin(char zflag=0, size_t *k=NULL) const;

    void info(const char* ="") const;
    void info(const char *F, int L, const char* ="") const;

    void print(
       const char *F, int L,
       const char *istr="", const char *fmt=""
     ) const;

    void print(const char *istr="", const char *fmt="") const {
       print(0,0,istr,fmt);
    };

    int  printdata(const char *istr="ans", const char *fmt="") const;
    void print_ref(const char *mark=" ***", const char *newl="\n") const;

    wbstring toStr() const {
       unsigned l=0, n=10*(numel()+SIZE.len); char s[n];
       l=snprintf(s,n,"[%s](%s)",
          wbvector<T>(numel(),data,'r').toStr().data,
          sizeStr().data
       );
       if (l>=n) wblog(FL,
          "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
       return s;
    };

    void contractDiag(
       unsigned ic, const wbarray &B,
       T a, T b, wbarray &C, char Iflag=0, char bflag='N'
    ) const;

    template<class TB, class TC>
    wbarray<TC>& contractMat(
    unsigned, const wbarray<TB>&, wbarray<TC>&, unsigned=1) const;

    wbarray& ContractMat(unsigned i1, const wbarray &B, unsigned i2=1){
        wbarray A; this->save2(A);
        return A.contractMat(i1,B,*this,i2);
    };

    template<class TB, class TC>
    wbarray<TC>& comm(
      const wbarray<TB>&, wbarray<TC>&, char aflag='N', char bflag='N'
    ) const;

    template<class TB, class TC>
    wbarray<TC>& acomm(
      const wbarray<TB>&, wbarray<TC>&, char aflag='N', char bflag='N'
    ) const;

    template<class TB, class TC>
    wbarray<TC>& contract(const char *F, int L,
       const wbvector<unsigned>&, const wbarray<TB>&,
       const wbvector<unsigned>&, wbarray<TC>&,
       const wbperm &pfinal=wbperm(), const T& afac=1.
    ) const;

    template<class TB>
    wbarray<T>& Contract(
       char*, const wbarray<TB>&, char*, 
       const wbperm &pfinal=wbperm()
    );

    template<class TB>
    wbarray<T>& Contract(
       const wbvector<unsigned>&, const wbarray<TB>&,
       const wbvector<unsigned>&,
       const wbperm &pfinal=wbperm()
    );

    template<class TB, class TC>
    wbarray<TC>& contract(const char *F, const int L,
       const char*, const wbarray<TB>&, const char*, wbarray<TC>&,
       const wbperm &pfinal=wbperm()
    ) const;

    template<class TB, class TC>
    wbarray<TC>& contract(
       char *i1, const wbarray<TB> &B, char *i2, wbarray<TC> &Q,
       const wbperm &pfinal=wbperm()) const
    {
       return contract(FL,i1,B,i2,Q,pfinal);
    };

    template<class TB, class TC>
    wbarray<TC>& contract(
       const wbvector<unsigned> &i1, const wbarray<TB> &B0,
       const wbvector<unsigned> &i2,
       wbarray<TC> &Q0, const wbperm &P) const
    {
       return contract(FL,i1,B0,i2,Q0,P);
    };

    template<class TB, class TC>
    wbarray<TC>& contract(
       unsigned i1, const wbarray<TB>&B, unsigned i2, wbarray<TC>&C,
       const wbperm &pfinal=wbperm()) const
    {
       wbvector<unsigned> I1, I2;
       if (!i1 || !i2) wblog(FL,"ERR Index into contract() is 1-based!");
       i1--; i2--; I1.init(1,&i1); I2.init(1,&i2);

       return contract(I1,B,I2,C, pfinal);
    };

    wbarray operator+(const wbarray &B) const {
       wbarray C;
       return plus(B,C);
    };

    wbarray operator-(const wbarray &B) const {
       wbarray C;
       return plus(B,C,-1);
    };

    wbarray operator*(const wbarray &B) const {
       wbarray C;
       wbvector<unsigned> ica(1), icb(1); ica[0]=1; icb[0]=0;
       return contract(FL,ica,B,icb,C);
    };

    wbarray& contract(unsigned i1, unsigned i2, wbarray &C) const;
    wbarray& contract(const wbMatrix<unsigned> &I12, wbarray &C) const;

    mxArray* toMx() const;
    mxArray* toMx_base() const { return cpyRange2Mx(data,SIZE); }
    mxArray* toMx_Struct() const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tflag=0) const;

    void put(const char *vname="ans", const char* ws="base") const {
       mxArray *a=toMx();
       mxPutAndDestroy(FL,a,vname,ws);
    };

    void put(const char *F, int L,
       const char *vname="ans", const char* ws="base"
    ) const {
       mxArray *a=toMx();
       wblog(F,L,"I/O putting array '%s' to %s.",vname,ws);
       mxPutAndDestroy(F,L,a,vname,ws);
    };

    void getReal(wbarray<double> &R) const;
    void getImag(wbarray<double> &I) const;
    void Conj();

    void set(const wbarray<double> &R, const wbarray<double> &I);

    const wbarray& opA(char tflag, wbarray& aux) {
       if (!isMatrix()) wblog(FL,
       "ERR Cannot use opA() with rank-%d object.", SIZE.len);

       if (tflag=='N') { return *this; } else
       if (tflag=='T' || tflag=='C') {
          permute(aux, wbperm("2 1")); if (tflag=='C') aux.Conj();
          return aux;
       }
       wblog(FL, "ERR Invalid flag %c<%d>", tflag, tflag); 
       return *this;
    };

    bool isSym_aux(
      const char *file, int line, const char *fct,
      const wbarray &B, double eps, double* xref,
      const char symflag='s',
      const char lflag=1
    ) const;



    using wb_sptr<T>::data;
    using wb_sptr<T>::isRef;
    using wb_sptr<T>::nrefs;

    wbvector<size_t> SIZE;


  protected:

    bool isDiag_aux(const T eps, const char* task) const;
    bool isDiag_aux(
       double *eps,
       const char* task) const;

    size_t serial_index(const INDEX_T *I, size_t n) const;



  private:

    void print_rec (UVEC&, size_t, const char*, const char*) const;

    void make2D() {
        if (SIZE.len<2) {
            SIZE.Resize(2); if (SIZE[0]) { SIZE[1]=1; }
        }
        else if (SIZE.len>2) wblog(FL,
        "ERR make2D() - invalid usage (got rank-%d)",SIZE.len);
    };

    void MemErrMsg(const char* file, int line, size_t s) {
        wblog(file,line, "ERR Out of memory? (%dx%d)", s, sizeof(T));
    };

    void DELETE_DATA() { wb_sptr<T>::delete_dref(); };

    void NEW_DATA_BARE(size_t s=0) { wb_sptr<T>::new_dref_bare(s); };
    void NEW_DATA(size_t s, const T* d0=NULL) {
       wb_sptr<T>::new_dref(s,d0);
    };

    void NEW(const wbvector<size_t> &S, size_t n, const T* d0=NULL) {
       SIZE.init(n,S);
       wb_sptr<T>::new_dref(SIZE.prod(0), d0);
    };

    void NEW (const wbvector<size_t> &S, const T* d0=NULL) {
       SIZE=S;
       wb_sptr<T>::new_dref(SIZE.prod(0), d0);
    };

    void RESIZE(const wbvector<size_t> &S, const T* d0=NULL) {
       if (d0!=NULL) { NEW(S,d0); }
       else if (S!=SIZE) {
          wbarray<T> X; this->save2(X);
          X.resize(S,*this);
       };
    };

    T& dummy() const {
       static T d; memset(&d,0,sizeof(T));
       return d;
    };

};


bool mxIsWbarray(
   const char *F, int L, const mxArray *a,
   const double **dr=NULL, const double **di=NULL, const mwSize **sz=NULL
);

bool mxIsWbarray(const mxArray *a) {
   return mxIsWbarray(0,0,a);
};

template<class T> inline
void householder(const T *u, T *x, size_t n, const T& eps);


template <class T>
class wbarrRef {
  public:
     
    wbarrRef() : dc(0.), fac(1.), tflag(0), D(NULL) {};

    wbarrRef(const wbarrRef<T> &r) {
       memcpy(this, &r, sizeof(r));
    };
    
    wbarrRef<T>& operator=(const wbarrRef<T> &r) {
       memcpy(this, &r, sizeof(r));
       return *this;
    };

    void print(const char *vname="ans") {
       printf("\nwbarrRef %s\n", vname);
       printf("... wbarray : 0x%lX (%s%s)\n", D,
          D && D->SIZE.len==1 ? "vec of length " : "",
          D ? D->sizeStr().data : "");
       printf("... dc      : %g\n", dc);
       printf("... fac     : %g\n", fac);
       printf("... tflag   : %d\n\n", tflag);
    };

    double dc;
    double fac;
    bool tflag;

    wbarray<T> *D;

  protected:
  private:
};


template<class T>
wbarray<T> TestContract(
    const wbarray<T> &A, C_UVEC i1,
    const wbarray<T> &B, C_UVEC i2, const wbperm &pfinal=wbperm());

template<class T>
void cell2mat(
   const wbMatrix< wbarrRef<T> > &C,
   wbarray<T> &M, WBINDEX &D1, WBINDEX &D2
);


template<class T> inline
size_t wbarray<T>::serial_index(const INDEX_T *I, size_t len) const {

   size_t i,k, r=SIZE.len, l=r-1, nz=r-len;

   if (len==0) return -1;
   if (!r || len>r || (len && data==NULL)) wblog(FL,
      "ERR %s[%s] having %s array",FCT,
      (WBINDEX(len,I)+1).toStr().data, sizeStr().data);

#ifdef CHECK_ELEMENT_RANGE
   for (i=0; i<len; i++) if (I[i]>=SIZE[i+nz]) {
       wblog(FL,"ERR %s[%s] index out of bounds (%s)",FCT,
      (WBINDEX(len,I)+1).toStr().data, sizeStr().data);
   }
#endif

   for (k=I[l-nz], i=l-1; i<r; i--) {
      k = k*SIZE[i] + ((i>=nz) ? I[i-nz] : 0);
   }

#ifdef CHECK_ELEMENT_RANGE
   if (k>=numel()) wblog(FL,
   "ERR %s() index out of bounds (%d; %s)",FCT,k,sizeStr().data);
#endif

   return k;
};


template <class T> inline
mxArray* wbarray<T>::toMx() const { return toMx_Struct(); }


template <> inline
mxArray* wbarray<unsigned> ::toMx() const { return toMx_base(); }
template <> inline
mxArray* wbarray<int>      ::toMx() const { return toMx_base(); }
template <> inline
mxArray* wbarray<char>     ::toMx() const { return toMx_base(); }
template <> inline
mxArray* wbarray<double>   ::toMx() const { return toMx_base(); }
template <> inline 
mxArray* wbarray<wbcomplex>::toMx() const { return toMx_base(); }



template<class T>
mxArray* wbarray<T>::toMx_Struct() const {

   size_t i,dim1,dim2,r=SIZE.len, s=numel();
   WBINDEX I(r);
   T A;

   if (r==2) { dim1=SIZE[0]; dim2=SIZE[1]; } else
   if (r==1) { dim2=SIZE[0]; dim1=(dim2 ? 1 : 0); } else
   if (r==0) { dim1=dim2=0; } else wblog(FL,
   "ERR %s() not implemented yet for rank>2 (%d)",SIZE.len);

   mxArray *S=A.mxCreateStruct(dim1,dim2);

   for (i=0; i<s; i++) data[i].add2MxStruct(S,i);

   return S;
};


template<class T>
mxArray* wbarray<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbarray<T>::add2MxStruct(mxArray *S, unsigned i, char tflag) const {

   if (tflag) { size_t s=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d",FCT,S,i+1,s);
   }

   mxSetCell(S,i,toMx());
}


template<class T>
void wbarray<T>::resize(const wbvector<size_t> &S, wbarray<T> &B) const {

   if (SIZE.len!=S.len) wblog(FL,
      "ERR length mismatch (%d/%d)",SIZE.len,S.len);

   if (S==SIZE) { B=(*this); return; }
   else { B.NEW(S); }

   size_t i,j,k, r=S.len, l=r-1, *s1=SIZE.data, *s2=B.SIZE.data;
   WBINDEX I(r); T *b=B.data;

   size_t s[r]; 
   for (i=0; i<r; i++) s[i]=MIN(S[i],SIZE[i]);

   while (I[l]<s[l]) {
      for (i=j=I[l], k=l-1; k<r; k--) {
          i=i*s1[k]+I[k];
          j=j*s2[k]+I[k];
      }

      b[j]=data[i];

      k=0; I[0]++;
      while(I[k]>=s[k] && k<l) { I[k]=0; ++I[++k]; }
   }
};


template<class T>
wbarray<T>& wbarray<T>::Enlarge(
   size_t dim,
   size_t D,
   size_t istart
){
   WBINDEX S(SIZE);
   wbarray<T> X;
   wbIndex I(S);

   if (!dim || dim>S.len) wblog(FL,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim,S.len); 
   if (D<S[--dim]) wblog(FL,
      "ERR %s() envalid new size (%d: %d/%d)",FCT,dim+1,S[dim],D); 
   if (istart+S[dim]>=D) wblog(FL,
      "ERR %s() index out of bounds (%d+%d/%d)",FCT,istart+1,S[dim],D); 

   S[dim]=D; X.init(S);
   if (istart) {
      for (size_t i=0; ++I; ++i) {
          I.data[dim]+=istart; X(I)=data[i];
          I.data[dim]-=istart;
      }
   } else { for (size_t i=0; ++I; ++i) X(I)=data[i]; }

   X.save2(*this);
   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::Append2(
   const char *F, int L,
   wbarray<T> &B,
   unsigned dim
){
   size_t i=0, D=0, ib=0;

   if (B.isEmpty()) { save2(B); return B; }

   if (SIZE.len!=B.SIZE.len) wblog(F_L,
      "ERR %s() incompatible objects (%s <> %s)",
       FCT,sizeStr().data,B.sizeStr().data);
   if (dim>=SIZE.len) wblog(F_L,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len); 

   for (i=0; i<SIZE.len; ++i) {
      if (i!=dim) { if (SIZE[i]!=B.SIZE[i]) wblog(F_L,
         "ERR %s() size mismatch (%s <> %s @ %d/%d; %d)",
          FCT,sizeStr().data,B.sizeStr().data,i+1,SIZE.len,dim+1);
      }
      else { D=SIZE[i]+B.SIZE[i]; ib=B.SIZE[i]; }
   }

   WBINDEX S(SIZE); S[dim]=D;
   wbIndex Ia(SIZE), Ib(B.SIZE);
   size_t &ia=Ia.data[dim];
   wbarray<T> X(S);

   for (i=0; ++Ib; ++i) X(Ib)=B.data[i];
   for (i=0; ++Ia; ++i) {
       ia+=ib; X(Ia)=data[i];
       ia-=ib;
   }

   X.save2(B); init();
   return B;
};

template<class T>
void wbarray<T>::Append2(
   const char *F, int L,
   wbarray<T> &B, unsigned dim,
   wbarray<double> &W, double eps
){
   size_t i,n, D=0, ib=0;
   const double *w=W.data;

   if (this->isEmpty()) return;
   if (dim>=SIZE.len) wblog(F_L,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len); 

   if (B.isEmpty()) {
      WBINDEX S(SIZE); S[dim]=0;
      B.init(S);
   }

   if (SIZE.len!=B.SIZE.len) wblog(F_L,
      "ERR %s() incompatible objects (%s <> %s)",
       FCT,sizeStr().data,B.sizeStr().data);
   for (i=0; i<SIZE.len; ++i) { if (i!=dim) {
       if (SIZE[i]!=B.SIZE[i]) wblog(F_L,
      "ERR %s() size mismatch (%s <> %s; %d; %d/d)",
       FCT,sizeStr().data,B.sizeStr().data,dim+1,i+1,SIZE.len);
   }}

   wbindex J(SIZE[dim]);
   size_t m=0, *j=J.data; n=W.numel();

   if (!W.isVector() || n!=SIZE[dim]) wblog(F_L,
      "ERR %s() invalid weigths W (%s <> %s @ %d)",
       FCT, W.sizeStr().data, sizeStr().data, dim+1);
   for (i=0; i<n; ++i) { if (w[i]<eps) { j[m++]=i; }}

   if (!m) return;
   if (m==n) { Append2(F_L,B,dim); return; }

   J.len=m;
   D=m+B.SIZE[dim]; ib=B.SIZE[dim];

   WBINDEX S(SIZE); S[dim]=m;
   wbIndex Ix(S), Ib(B.SIZE);  S[dim]=D;
   wbindex Ia(S.len),J2;
   wbarray<T> X(S);

   size_t &ia=Ia.data[dim], &ix=Ix.data[dim];

   for (i=0; ++Ib; ++i) { X(Ib)=B.data[i]; }
   for (i=0; ++Ix; ++i) { Ia.setp(Ix.data); ia=j[ia];
       ix+=ib; X(Ix)=(*this)(Ia);
       ix-=ib;
   }


   X.save2(B); 
   J.invert(SIZE[dim],J2); Select0(J2,dim);
};


template<class T>
wbarray<T>& wbarray<T>::init(
   const char *F, int L,
   const mxArray *a, char ref, char vec
){
   size_t n,r=mxGetNumberOfDimensions(a);
   const mwSize *sz=0;
   const double *dr=0, *di=0;

   if (!mxIsWbarray(F,L,a,&dr,&di,&sz)) { init(); return *this; }

   bool isdbl=(typeid(T)==typeid(double));
   bool isz=(typeid(T)==typeid(wbcomplex));

   if (ref && !isdbl && !isz) wblog(F_L,
      "ERR %s() cannot reference mxArray using wbarray<%s>",
       FCT,getName(typeid(T)).data);
   if (ref && di && isz) { static unsigned char nr_call=0;
      if (++nr_call==1) wblog(F_L,
      "WRN %s() cannot reference complex data (%s; %d,%d)",
       FCT,getName(typeid(T)).data,ref,isz); }
   if (di && !isz) wblog(F_L,
      "ERR %s() got complex data (%s; %d)",
       FCT,getName(typeid(T)).data,isz);

   wbvector<size_t> S;
   S.initT(r,sz); n=S.prod(0); init();

   if (vec && r==2) {
      if (S[1]==1) S.len=1; else
      if (S[0]==1) { S[0]=S[1]; S.len=1; }
   }

   if (isdbl && ref) {
      SIZE=S;
      wb_sptr<T>::init2ref((T*)dr);
   }
   else { init(S); 
      if (isdbl) memcpy(data,dr,n*sizeof(double)); else
      if (isz) cpyZRange(dr,di,data,n);
      else cpyRange(dr,data,n);
   }

   return *this;
};



template<>
wbarray<char>& wbarray<char>::init(
   const char *F, int L,
   const mxArray *a, char ref, char vec
){
   if (!a || mxIsEmpty(a)) { init(); return *this; }

   size_t n, r=mxGetNumberOfDimensions(a);
   const mwSize *sz=mxGetDimensions(a); 

   if (mxIsComplex(a))
      wblog(F_L,"ERR severe datatype mismatch (%s; %d/%d)",
      mxGetClassName(a), sizeof(mxChar), sizeof(short));
   if (ref) wblog(F_L,
      "WRN %s() ref=%d will be ignored (having %s)",FCT,ref,
      mxGetClassName(a)
   );

   WBINDEX S; S.initT(r,sz); n=S.prod(0);

   if (vec && r==2) {
      if (S[1]==1) S.len=1; else
      if (S[0]==1) { S[0]=S[1]; S.len=1; }
   }

   init(S);

   long q=cpyRangeMx(a,data,n);


   if (q<0) wblog(FL,"ERR %s() invalid numeric type "
      "(%s; %d)",FCT,a ? mxGetClassName(a) : "null",q);

   return *this;
};


template<class T>
wbarray<T>& wbarray<T>::init(
   const char *F, int L, const mxArray *a, wbperm &P
){

   size_t i,j,k,s,l,r=mxGetNumberOfDimensions(a);
   const mwSize *ip=mxGetDimensions(a);
   WBINDEX Sp,I,S0;
   double *dr, *di;

   bool isz = (typeid(T)==typeid(wbcomplex));

   if (ip==NULL) wblog(F,L,
      "ERR mxGetDimensions returned (%lX, dim=%d) ???",ip,r);
   if (!mxIsDouble(a) || mxIsSparse(a)) wblog(F,L,
      "ERR %s() invalid input array (%d/%d)",FCT,
       mxIsDouble(a), mxIsSparse(a));
   dr=mxGetPr(a); di=mxGetPi(a);
   if (di && !isz) wblog(F,L,
      "ERR %s() got complex input data (%s; %d,%d)",
       FCT,getName(typeid(T)).data,ref,isz);

   if (P.isEmpty()) P.Index(r);

   S0.init(r); for (i=0; i<r; i++) S0[i]=(size_t)ip[i];
   S0.permute(Sp,P); s=S0.prod(0);

   init(Sp); if (r==0) return;
   I.init(r); l=r-1;

   for (i=0; i<s; i++) {
       for (j=I[P[l]],k=l-1; k<l; k--) j = j*Sp[k] + I[P[k]];

       DSET(data[j], dr, di, i);

       k=0; I[0]++;
       while(I[k]>=S0[k] && k<l) { I[k]=0; ++I[++k]; }
   }
}


template<class T> inline
void wbarray<T>::initTst() {
 
   size_t i,j,k, r=SIZE.len, l=r-1, s=numel(), x;
   WBINDEX I(r);

   for (i=0; i<s; i++) {
      for (j=I[l], k=l-1; k<r; k--) j = j*SIZE[k] +I[k];
      for (x=I[0]+1, k=1; k<r; k++) x = x*10 + (I[k]+1);

      data[j]=(T)x;

      k=0; I[0]++;
      while(I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }
}


template<class T> inline 
wbarray<T>& wbarray<T>::initIdentity(size_t d, char dflag, T dval) {

   if (dflag) { init(d); set(T(1)); }
   else {
      init(d,d);
      for (size_t i=0; i<d; ++i) data[i+d*i]=dval;
   }

   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::initIdentity(
   const WBINDEX &S, char dflag
){
   size_t i,d,m=S.len/2;
   size_t const* const s=S.data;

   if (S.len%2) wblog(FL,
   "ERR %s() requires even rank object (%s)",FCT,sizeStr().data);
   for (i=0; i<m; i++) if (s[i]!=s[i+m]) wblog(FL,
   "ERR %s() requires symmetric object (%s)",FCT,sizeStr().data);
   if (m>1 && dflag) wblog(FL,
   "ERR %s() rank-2 object required with dflag (%s)",FCT,sizeStr().data);

   if (!m) { init(); return *this; }
   if (dflag) { init(s[0]).set(T(1)); return *this; }
   if (m==1 && s[0]==1) { init(S).set(T(1)); return *this; }
  
   init(S); for (d=s[0], i=1; i<m; i++) d*=s[i];
   for (i=0; i<d; i++) data[i+d*i]=1;

   return *this;
};


template<class T> inline
wbarray<T>& wbarray<T>::initIdentityB(
   size_t d1, size_t d2, size_t k
){
   init(d1,d2);

   if (int(k)>=0) {
      if (k<d2) {
         size_t j=k, dx=d1+1; T *x=data+j*d1;
         for (; j<d2 && (j-k)<d1; ++j, x+=dx) { x[0]=1; }
      }
      else wblog(FL,"WRN %s() "
     "got all-zero block (%dx%d w/ k=%d)",FCT,d1,d2,k);
   }
   else { k=-k;
      if (k<d1) {
         size_t i=k, dx=d1+1; T *x=data+i;
         for (; i<d1 && (i-k)<d2; ++i, x+=dx) { x[0]=1; }
      }
      else wblog(FL,"WRN %s() "
     "got all-zero block (%dx%d w/ k=%d)",FCT,d1,d2,-k);
   }

   return *this;
};


template<class T> inline
wbarray<T>& wbarray<T>::initIdentityB3(
   size_t d1, size_t d2, size_t D,
   size_t i0,
   T one
){
   size_t d12=d1*d2, i=0;

   if (i0+d12>D) wblog(FL,
      "ERR %s() index out of bounds\n(%d*%d) x %d starting at %d",
       FCT,d1,d2,D,i0+1);

   init(d1,d2,D); i0*=d12;
   for (; i<d12; i++) data[i0+i+d12*i]=one;

   return *this;
};


template<class T> inline 
void wbarray<T>::Diag2Vec() {
   size_t i,m;
   WBINDEX S(1);

   if (!isSMatrix()) wblog(FL,
      "ERR Calling %s for rank-%d object (%s).",
       FCT, SIZE.len, sizeStr().data);

   m=S[0]=SIZE[0];
   for (i=1; i<m; i++) data[i]=data[i+i*m];

   Resize(S);
}


template<class T>
void wbarray<T>::contractDiag(
   unsigned ic, const wbarray<T> &B,
   T a, T b, wbarray<T> &C, char Iflag, char bflag
) const {

   if (B.SIZE.len!=2) wblog(FL,"ERR invalid rank-%d",SIZE.len);
   if (ic<1 || ic>SIZE.len) wblog(FL,
      "ERR invalid contraction index (ic=%d)",ic);
   if (!strchr("NTCc",bflag)) wblog(FL,
   "ERR invalid flag %c<%d>",bflag,bflag);

   size_t i,k, s=numel(), r=SIZE.len, l=r-1;
   size_t n=B.SIZE.min(), m=B.SIZE[0];
   WBINDEX I_(r);
   wbvector<T> h_(n);
   wbarray<T> X;

   size_t *I=I_.data, *S=SIZE.data;
   T *h=h_.data;

   ic--;
   if (n!=SIZE[ic]) wblog(FL,"ERR dimension mismatch (%d/%d)",n,SIZE[ic]);
   if (n==0) return;

   for (i=0; i<n; i++) h[i]=B.data[i+i*m];
   if (bflag=='C' || bflag=='c') h_.Conj();

   X=*this; if (a!=1) X*=a;

   for (i=0; i<s; i++) {
       if (Iflag)
            X.data[i]/=(b-h[I[ic]]+1E-33);
       else X.data[i]*=(b-h[I[ic]]);

       k=0; I[0]++;
       while(I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   if (C.isEmpty()) X.save2(C); else {
      if (C.SIZE!=SIZE) wblog(FL,"ERR size mismatch in C: %s <> %s",
          sizeStr().data, C.sizeStr().data);
      C+=X;
   }
};


template<class T> 
void wbarray<T>::ExpandDiagonal() {

   if (SIZE.len!=1) { sprintf(str,
      "%s() called for %s object",FCT,sizeStr().data);
       wblog(FL,SIZE.len==2 ? "WRN %s" : "ERR %s",str); return;
   }

   size_t i,n=numel();
   wbarray<T> X; save2(X); init(n,n);

   for (i=0; i<n; i++) data[i+n*i]=X.data[i];
};


template<class T>
void wbarray<T>::ExpandDiagonal(unsigned i1, unsigned i2) {

   if (i1>=SIZE.len || i2>=SIZE.len || i1==i2) wblog(FL,
      "ERR %s() invalid indices (%d %d; %d)",FCT,i1+1,i2+1,SIZE.len);
   if (SIZE[i1]!=1 && SIZE[i2]!=1) wblog(FL,
      "ERR %s() array should be diagonal in (%s; %d %d)",
       FCT, sizeStr().data, i1+1, i2+1);
   if (SIZE[i1]==1 && SIZE[i2]==1) return;

   if (SIZE[i1]!=1) SWAP(i1,i2);

   size_t i,j,k,*s0,*s, r=SIZE.len, l=r-1;
   WBINDEX S=SIZE;
   wbarray<T> X;
   wbindex I(r);
 
   save2(X); s0=X.SIZE.data;
   S[i1]=S[i2]; init(S); s=SIZE.data;

   for (i=0; I[l]<s0[l]; i++) {
       for (j=0,k=l; k<r; k--) {
          if (j) { j*=s[k]; }
          j += (k!=i1 ? I[k] : I[i2]);
       }
       data[j] = X.data[i];

       k=0; I[0]++;
       while(I[k]>=s0[k] && k<l) { I[k]=0; ++I[++k]; }
   }
}


template<class T>
void wbarray<T>::Expand2Projector(T eps) {

   size_t i,j,m,n=numel();
   WBINDEX S(2);
   wbvector<char> mark(n);

   if (n==0) { if (SIZE.len>1) SIZE[0]=0; return; }
   if (!isVector()) wblog(FL,
      "ERR %s() array must be 1-dim vector (got %s)",
       FCT, sizeStr().data
   );

   for (m=i=0; i<n; i++) if (ABS(data[i])>eps) { mark[i]=1; m++; }

   S[0]=n; S[1]=m; init(S);

   for (i=j=0; i<n; i++)
   if (mark[i]) { data[i+(j++)*n]=1; }
}


template<class T> inline
T wbarray<T>::froNorm2(const wbarray<T> &B) const {
   T x=0;

   if (!hasSameSize(B)) wblog(FL,"ERR %s() size mismatch (%s/%s)",
   FCT, sizeStr().data, B.sizeStr().data);

   for (size_t n=numel(), i=0; i<n; i++) x+=CONJ(data[i])*B.data[i];
   return x;
};


template<class T> inline
T wbarray<T>::norm2() const {
   T d=0;
   for (size_t n=numel(), i=0; i<n; i++) d+=data[i]*data[i];
   return d;
};


template<> inline
wbcomplex wbarray<wbcomplex>::norm2() const {
   double d=0;

   for (size_t n=numel(), i=0; i<n; i++)
   d+=data[i].abs2();

   return d;
}


template<class T> inline
T wbarray<T>::normCol2(size_t k) const {

   if (!isMatrix()) wblog(FL,
      "ERR %s() matrix expected (%s)",FCT,sizeStr().data);
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);

   size_t i=0, dim1=SIZE[0];
   const T *d=data+k*dim1; T d2=0;

   for (i=0; i<dim1; i++) d2+=d[i]*d[i];
   return d2;
};


template<> inline
wbcomplex wbarray<wbcomplex>::normCol2(size_t k) const {

   if (!isMatrix()) wblog(FL,
      "ERR %s() matrix expected (%s)",FCT,sizeStr().data);
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);

   size_t i=0, dim1=SIZE[0];
   const wbcomplex *d=data+k*dim1; double d2=0;

   for (i=0; i<dim1; i++) d2+=d[i].abs2();
   return d2;
};


template<class T> inline
T wbarray<T>::scalarProd(const wbarray<T> &B) const {
   T x=0; 

   if (!hasSameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,
       sizeStr().data, B.sizeStr().data);

   for (size_t n=numel(), i=0; i<n; i++)
   x+=data[i]*B.data[i];

   return x;
}


template<> inline 
wbcomplex wbarray<wbcomplex>::scalarProd(
  const wbarray<wbcomplex> &B
) const {

   wbcomplex x=0;

   if (!hasSameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",
       FCT, sizeStr().data, B.sizeStr().data);

   for (size_t n=numel(), i=0; i<n; i++)
   x+=data[i]*B.data[i].conj();

   return x;
}


template<class T> inline 
T wbarray<T>::sum() const {
   T x=0; 

   for (size_t n=numel(), i=0; i<n; i++)
   x+=data[i];

   return x;
}


template<class T> inline 
wbarray<T>& wbarray<T>::sum(
    const WBINDEX &K,
    wbarray<T>&A
) const {

    size_t i,j,k,q, r=SIZE.len, l=r-1, s=numel();
    WBINDEX S=SIZE;
    wbvector<char> mark(r);
    char *m=mark.data;
    wbindex I;

    for (k=i=0; i<K.len; i++) {
       if (K[i]<SIZE.len) { mark[K[i]]++; k++;
          if (mark[K[i]]>1) wblog(FL,
           "ERR %s() index not unique (%s; %d)",
            FCT,(K+1).toStr().data,SIZE.len);
       }
       else wblog(FL,
        "ERR %s() index out of bounds (%s; %d/%d)",
         FCT, (K+1).toStr().data, K[i],SIZE.len);
    }

    for (k=i=0; i<mark.len; i++) if (!mark[i]) S[k++]=SIZE[i];
    if (k) S.len=k;
    else { A.init(1,1); A[0]=addRange(data,s); return A; }

    A.init(S); I.init(SIZE.len);
    for (q=l; q<r; q--) if (!m[q]) break;

    for (i=0; i<s; i++) {
        for (j=I[q], k=q-1; k<r; k--) if (!m[k]) j=j*SIZE[k]+I[k];

        A.data[j]+=data[i];

        k=0; I[0]++;
        while(I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
    }

    return A;
}


template<class T> inline
T wbarray<T>::trace() const {
   T x=0; 
   size_t n;

   if (!isMatrix()) wblog(FL,
   "ERR cannot calculate trace for rank-%d object", SIZE.len);
   if (SIZE[0]!=SIZE[1]) wblog(FL,
   "WRN Calculating trace of non-square array (%s)",sizeStr().data);

   n=SIZE[0];

   for (size_t m=SIZE.min(), i=0; i<m; i++)
   x+=data[i+i*n];

   return x;
}


template<class T> inline 
double wbarray<T>::normDiff2(const wbarray<T> &B) const {
   if (SIZE!=B.SIZE) wblog(FL,
      "ERR size mismatch (%s; %s)",sizeStr().data,B.sizeStr().data);
   return rangeNormDiff2(data, B.data, numel());
};

template<class T> inline 
double wbarray<T>::normDiff2(const wbvector<T> &b) const {
   if (!isVector() || numel()!=b.len) wblog(FL,
      "ERR size mismatch (%s; %d)",sizeStr().data,b.len);
   return rangeNormDiff2(data, b.data, numel());
};


#endif

