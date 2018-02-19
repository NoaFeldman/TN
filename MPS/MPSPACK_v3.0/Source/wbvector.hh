#ifndef __WB_SIMPLE_VECTOR_HCC__
#define __WB_SIMPLE_VECTOR_HCC__

template <class T> class wbvector;

/* ---------------------------------------------------------------- *
 * class wbvector
 * AWb (C) Jan 2006
 * ---------------------------------------------------------------- */

WBINDEX Index(INDEX_T i1, INDEX_T i2);
char isUniqueIdxSet(const WBINDEX &I, INDEX_T imax);
char validPerm(const WBPERM &P, PERM_T d=-1);
void getIPerm (const WBPERM &P, WBPERM &iP);

   inline bool mxIsWbvector(
      const char *F, int L, const mxArray *a,
      size_t *n_=NULL, const char *istr=NULL, char dflag='d');

   inline bool mxIsWbvector(const mxArray *a) {
      return mxIsWbvector(0,1,a);
   }

template <class T>
class wbvector {

  public:

    wbvector() : data(NULL), len(0), isref(0) {};

    wbvector(size_t l)
     : data(NULL), len(0), isref(0) { NEW(l); };

    wbvector(size_t l, const T* d0, char ref_flag=0)
     : data(NULL), len(0), isref(0) {
       if (ref_flag) init2ref(l,d0);
       else NEW(l,d0);
    };

    wbvector(const wbvector &v)
     : data(NULL), len(0), isref(0) { NEW(v.len, v.data); };

    template <class T2>
    wbvector(const wbvector<T2> &v)
     : data(NULL), len(0), isref(0) { initT(v.len,v.data); }

    wbvector(const wbvector &v1, const wbvector &v2)
     : data(NULL), len(0), isref(0) {
       size_t i,k=0;
       NEW(v1.len+v2.len);

       for (i=0; i<v1.len; ++i) data[k++]=v1[i];
       for (i=0; i<v2.len; ++i) data[k++]=v2[i];
    };

    wbvector(const char *F, int L, const mxArray *a, char check_type=1)
     : data(NULL), len(0), isref(0) {
       init(F,L,a,NULL,check_type);
    };

    wbvector(const mxArray *a, char check_type=1)
     : data(NULL), len(0), isref(0) {
       init(FL,a,NULL,check_type);
    };

   ~wbvector() {
       if (!isref && data) { WB_DELETE(data); }
    };

    wbvector& init(size_t l=0) { return RENEW(l); };
    wbvector& init(size_t l, const T* d0) { return RENEW(l,d0); };
       
    wbvector& init(size_t l, const T* d0, char ref) { 
       if (ref)
            return init2ref(l,d0);
       else return RENEW(l,d0); 
    };

    wbvector& init2val(size_t l, const T& x) { return RENEW_VAL(l,x); };

    wbvector& init_stride(
       const T* d0, size_t n, size_t m, size_t D0
    ){
       RENEW(n*m);
       Wb::cpyStride(data, d0, n, m, -1, D0);
       return *this;
    };


    wbvector& init(const wbvector &v) {
       RENEW(v.len);
       for (size_t i=0; i<len; i++) data[i]=v.data[i];
       return *this;
    };

    wbvector& init(size_t n, const wbvector &v) {
       RENEW(n); if (n>v.len) n=v.len;
       for (size_t i=0; i<n; ++i) data[i]=v.data[i];
       return *this;
    };

    template <class T2>
    wbvector<T>& initT(const wbsparray<T2> &S);

    template <class T2>
    wbvector<T>& initT(size_t l, const T2* v) {
       RENEW(l);
       for (size_t i=0; i<len; ++i) data[i]=T(v[i]);
       return *this;
    };

    template <class T2>
    wbvector<T>& initT(const wbvector<T2> &v) {
       RENEW(v.len);
       for (size_t i=0; i<len; ++i) data[i]=T(v.data[i]);
       return *this;
    };

    template <class T2>
    wbvector<T>& initT(const char *F, int L, const wbvector<T2> &v);


    wbvector& init2ref(const wbvector &v) {
       return init2ref(v.len, v.data);
    };

    wbvector& init2ref(size_t l, const T* d0) {
       if (data && !isref) RENEW(0); isref=1;
       len=l; data=(T*)d0;

       if (len && !data) wblog(FL,
          "ERR init2ref with data=0 (%d) !??", len);
       return *this;
    };

    wbvector& unRef();

    wbvector& init(size_t l1, const T* d1, size_t l2, const T* d2) {
       RENEW(l1+l2); MEM_CPY<T>(data,len,l1,d1,d2);
       return *this;
    };

    int init(
       const char*, int, const mxArray *a,
       const char *istr=NULL, char check_type=1
    );

    int init(
       const mxArray *a, const char *istr=NULL, char check_type=1
    ){ return init(0,0,a,istr,check_type); };

    int initMPFR(
       const char *F, int L, const mxArray *a, char base=-1
    );

    void initDef(size_t n) {
        if (data) { WB_DELETE(data); }
        len=n; if (n) { WB_NEW(data,n); }
    };

    void initDef(const wbvector &v) {
       if (len!=v.len) initDef(v.len);
       for (size_t i=0; i<v.len; i++) data[i]=v[i];
    };

    wbvector& operator= (const wbvector &v) {
       if (this!=&v) RENEW(v.len, v.data);
       return *this;
    };

    template<class TB>
    wbvector<T>& operator= (const wbvector<TB> &v) {
       init(v.len);
       for (size_t i=0; i<v.len; ++i) { data[i]=T(v.data[i]); }
       return *this;
    };

    wbvector& set(const T &x) {
       for (size_t i=0; i<len; i++) data[i]=x;
       return *this;
    };

    wbvector& setMin(const T &x) {
       for (size_t i=0; i<len; ++i) { if (data[i]<x) data[i]=x; }
       return *this;
    };
    wbvector& setMax(const T &x) {
       for (size_t i=0; i<len; ++i) { if (data[i]>x) data[i]=x; }
       return *this;
    };

    wbvector& set(const wbvector<size_t> &I, const T &x);
    wbvector& add(const wbvector<size_t> &I, const T &x);

    wbvector& Set(const wbvector<size_t> &I, const wbvector &v);

    void Sort(char dir=+1) { wbperm p; Wb::hpsort(*this,p,dir); };
    void Sort(wbperm &p, char dir=+1) {  Wb::hpsort(*this,p,dir); };
    void sort(wbvector &v2) const { wbperm p;  v2=(*this); Wb::hpsort(v2,p); };
    void sort(wbvector &v2, wbperm &p) const { v2=(*this); Wb::hpsort(v2,p); };

    bool isUnique() const {
       wbvector v(*this); wbperm p; Wb::hpsort(v,p);
       for (size_t i=1; i<v.len; i++) if (v[i]==v[i-1]) return 0;
       return 1;
    };

    int contains(const T& x) const {
       for (size_t i=0; i<len; ++i) { if (data[i]==x) return 1; }
       return 0;
    };

    bool isNormal() const;

    size_t find1(const T& x) const {
       for (size_t i=0; i<len; ++i) { if (data[i]==x) return i; }
       return -1;
    };
    size_t find1() const {
       for (size_t i=0; i<len; ++i) { if (data[i]) return i; }
       return -1;
    };

    size_t count(const T& x) const {
       size_t i=0, n=0;
       for (; i<len; ++i) if (data[i]==x) ++n;
       return n;
    };

    size_t icount(const T& x) const {"inverse" count
       size_t i=0,n=0;
       for (; i<len; ++i) if (data[i]!=x) ++n;
       return n;
    };

    size_t nnz() const { return icount(0); };

    size_t numZeros(T eps=1E-14) const {
       size_t i=0,n=0;
       if (eps>0) 
            { for (; i<len; ++i) if (ABS(data[i])<eps) ++n; }
       else { for (; i<len; ++i) if (data[i]) ++n; }
       return n;
    };

    wbindex& find(const T &x, wbindex &I, char iflag=0) const;
    wbindex& find(const T &x, wbindex &I, wbindex &Ix) const;
    wbindex& findGT(const T &x, wbindex &I, char iflag=0) const;
    wbindex& findGT(const T &x, wbindex &I, wbindex &Ix) const;

    wbindex& findzeros(wbindex &I) const { return find(0,I); }
    wbindex& find(wbindex &I) const { return find(0,I,'i'); }

    wbindex& find(wbindex &I, wbindex &Ix) const {
       find(0,Ix,I); return I; }

    wbvector<size_t> find() const {
       wbvector<size_t> I;
       return find(I);
    };

    size_t numelFind(const T&) const;
    size_t numelFindGT(const T&) const;

    wbindex& findRange(
       const T &x1, const T &x2, wbindex &I, const char *w=NULL
    ) const;

    int findClosestSorted(T r) const;
    int findClosest(T r) const;

    size_t findValues(const wbvector &v, wbindex &I) const;

    wbvector& setp(const T *x) {
       for (size_t i=0; i<len; i++) data[i]=x[i];
       return *this;
    };

    void set(size_t offset, const T *x, size_t n) {
       if (offset+n>=len) wblog(FL,
          "ERR index out of bounds (%d+%d; %d)",offset,n,len);
       MEM_CPY<T>(data+offset,n,x);
    };

    void reset() { if (len) MEM_SET<T>(data,len); };

    void cpy(const T* v) {
       if (len) MEM_CPY<T>(data,len,v);
    };

    void cpy(const T *x, size_t m) {
       if (m>len) { wblog(FL,"ERR index out of bounds (%d/%d)",m,len); m=len; }
       if (m) MEM_CPY<T>(data,m,x);
    };

    wbvector& save2(wbvector &v);

    void swap (wbvector &v) {
       if (this!=&v) {
          SWAP(len,   v.len  );
          SWAP(data,  v.data );
          SWAP(isref, v.isref);
       }
    };

    wbvector& swap(size_t i, size_t j) {
       if (i<len && j<len)
            { T x=data[i]; data[i]=data[j]; data[j]=x; }
       else { wblog(FL,
         "ERR %s() index out of bounds (%d,%d/%d)",FCT,i+1,j+1,len); }
       return *this;
    };

    wbvector& Collect(T& A1, T& A2) {
       init(2); A1.save2(data[0]); A2.save2(data[1]);
       return *this;
    };

    wbvector& Collect(T& A1, T& A2, T& A3) {
       init(3); A1.save2(data[0]); A2.save2(data[1]); A3.save2(data[2]);
       return *this;
    };

    wbvector& Skip(size_t);
    wbvector& Skip(const WBINDEX &, char rflag=0);

    wbvector& Move2_FE  (const WBINDEX &I, char FE_flag);
    wbvector& Move2front(const WBINDEX &I) { return Move2_FE(I,'F'); }
    wbvector& Move2end  (const WBINDEX &I) { return Move2_FE(I,'E'); }

    wbvector& Move2front (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[k];
       for (size_t i=k; i>0; i--) data[i]=data[i-1];
       data[0]=xk; return *this;
    }

    wbvector& Move2end (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[k];
       for (size_t i=k+1; i<len; i--) data[i-1]=data[i];
       data[len-1]=xk; return *this;
    }

    wbvector& mvLastTo (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[len-1];
       for (size_t i=len-1; i>k; i--) data[i]=data[i-1];
       data[k]=xk; return *this;
    }

    wbvector& mvFirstTo (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[0];
       for (size_t i=0; i<k; i++) data[i]=data[i+1];
       data[k]=xk; return *this;
    }


    bool checkSameLength(
       const wbvector &a,
       const char *F=NULL, int L=0, const char *istr=""
    )  const;

    bool checkSameLength(
       const wbvector &a,
       const wbvector &b,
       const char *F=NULL, int L=0, const char *istr=""
    )  const;

    void operator+= (const T c) {
       if (c==T(0)) return;
       for (size_t i=0; i<len; ++i) data[i]+=c;
    };

    void operator-= (const T c) {
       if (c==T(0)) return;
       for (size_t i=0; i<len; ++i) data[i]-=c;
    };

    void operator*= (const T c) const {
       if (c==T(1)) return;
       if (c==T(0)) MEM_SET<T>(data,len); else
       if (c==T(-1))
            for (size_t i=0; i<len; ++i) data[i]=-data[i];
       else for (size_t i=0; i<len; ++i) data[i]*=c;
    };

    void operator/= (const T c) {
       if (c==T(1)) return;
       if (c==T(0)) wblog(FL,"ERR c=%g DIV/0", double(c));
       if (c==T(-1))
            for (size_t i=0; i<len; ++i) data[i]=-data[i];
       else for (size_t i=0; i<len; ++i) data[i]/=c;
    };

    wbvector operator* (const T c) const {
       wbvector u(*this); if (c!=1) u*=c;
       return u;
    };

    wbvector operator+ (const T c) const {
       wbvector u(*this);
       if (c!=0) for (size_t i=0; i<u.len; i++) u.data[i]+=c;
       return u;
    };

    wbvector operator- (const T c) const {
       wbvector u(*this);
       if (c!=0) for (size_t i=0; i<u.len; i++) u.data[i]-=c;
       return u;
    };

    wbvector operator+ (const wbvector &b) const {
       wbvector u(*this); checkSameLength(b,FLF);
       for (size_t i=0; i<u.len; i++) u.data[i]+=b.data[i];
       return u;
    };

    wbvector operator- (const wbvector &b) const {
       wbvector u(*this); checkSameLength(b,FLF);
       for (size_t i=0; i<u.len; i++) u.data[i]-=b.data[i];
       return u;
    };

    void operator+= (const wbvector &b) {
       checkSameLength(b,FLF);
       for (size_t i=0; i<len; i++) data[i]+=b.data[i];
    };

    void operator-= (const wbvector &b) {
       checkSameLength(b,FLF);
       for (size_t i=0; i<len; i++) data[i]-=b.data[i];
    };

    void Plus (const wbvector &b, T bfac=1) const {
       size_t i=0; checkSameLength(b,FLF);

       if (bfac==+1) for (; i<len; i++) data[i]+=b.data[i]; else
       if (bfac==-1) for (; i<len; i++) data[i]-=b.data[i]; else
       if (bfac!= 0) for (; i<len; i++) data[i]+=(bfac*b.data[i]);
    };

    void plus (
       const wbvector &b, wbvector &c, T bfac=1, T cfac=0
    ) const {
       if (cfac) {
              checkSameLength(b,c,FLF); c*=cfac; c.Plus(*this); }
       else { checkSameLength(b,  FLF); c=*this; }
       
       c.Plus(b,bfac);
    };

    wbvector& times(T fac, wbvector &x) const {
        x=*this; x*=fac; return x;
    };

    wbvector& RevertSigns();


    int isum() const {
       int s=0; if (len==0) return s;

       s=data[0]; for (size_t i=1; i<len; i++) s+=data[i];
       return s;
    };

    T sum() const {
       T s=T(); if (len==0) return s;

       s=data[0]; for (size_t i=1; i<len; i++) s+=data[i];
       return s;
    };

    T sum(size_t i1, const char *estr) const {
       if (!estr || strcmp(estr,"end")) wblog(FL,
       "ERR Invalid end string `%s'.", estr);

       return sum(i1,len-1);
    };

    T sum(size_t i1, size_t i2) const {
       T s=T();

       if (len==0) if (i1==0 && i2==0) return s;
       if (i1>=len || i2>=len) wblog(FL,
          "ERR index out of bounds (%d,%d; %d)",i1,i2,len);
       if (i1>i2) return s;

       s=data[i1];
       for (size_t i=i1+1; i<=i2; i++) s+=data[i];

       return s;
    };

    T cumsum0(wbvector &cs, char xflag=0) const {
       size_t i=1, n=(xflag ? len+1 : len);
       cs.init(n); if (!len) { if (xflag) cs[0]=0; return 0; }
       for (cs[0]=0; i<n; ++i) {
           cs[i]=cs[i-1]+data[i-1];
       }
       return (xflag ? cs[len] : cs[len-1]+data[len-1]);
    };

    T Cumsum(char zflag=0) const {
       if (!len) return 0;
       if (zflag) {
          T q,x=data[0]; data[0]=0;
          for (size_t i=1; i<len; ++i) {
             q=data[i]; data[i]=data[i-1]+x; x=q;
          }
          return data[len-1]+x;
       }
       else {
          for (size_t i=1; i<len; ++i) { data[i]+=data[i-1]; }
          return data[len-1];
       }
    };
    
    wbvector& Eldiv( const wbvector &v) const {
       if (len!=v.len) wblog(FL,
          "ERR %s() size mismatch (%d/%d)",FCT,len,v.len);
       for (size_t i=0; i<len; ++i) {
          data[i] /= v.data[i];
          if (v.data[i]==0) wblog(FL,
             "ERR %s() div/0 !??  (%g)",FCT,double(data[i])
          );
       }
       return *this;
    };

    wbvector& tensorProd("col-major")
       const wbvector &v2,
       wbvector &vv
    ) const {

       vv.init(len*v2.len); if (!vv.isEmpty()) {
          size_t i,j,l=0;
          const T *d2=v2.data; T *d=vv.data;

          for (j=0; j<v2.len; ++j)
          for (i=0; i<len; ++i, ++l) d[l]=data[i]*d2[j];
       }

       return vv;
    };

    T scalarProd(const wbvector &v) const {
        if (len!=v.len) wblog(FL,
           "ERR scalarProd - length mismatch (%d/%d)",len,v.len);

        if (len) return scalarProd(v.data);
        else return 0;
    };

    T scalarProd(const T *v) const {
        T s=T(); for (size_t i=0; i<len; i++) s+=(data[i]*CONJ(v[i]));
        return s;
    };

    T norm2() const {
       T s=0; if (len==0) return s;
       s=data[0]*data[0];

       for (size_t i=1; i<len; i++) s+=NORM2(data[i]);
       return s;
    };
    T norm() const { return SQRT(norm2()); }

    double normDiff(const wbvector &v, T fac=1) const {
       return std::sqrt(double(normDiff2(v,fac)));
    };

    T normDiff2(const wbvector &v, T fac=1) const {
       if (len!=v.len) wblog(FL,
          "ERR %s() vector length mismatch (%d/%d)",FCT,len,v.len);
       return rangeNormDiff2(data,v.data,len,fac);
    };

    T GSProject(const T *v, char isnorm=0, char tflag=0) const {
       return gs_project_range(data,v,len,1,isnorm,tflag);
    };

    wbvector& Normalize(const char *F=NULL, int L=0) {
       double n=std::sqrt(norm2());
       if (n==0) {
          if (F) wblog(F,L,"WRN normalize zero vector!? - randomize");
          setRand(1,-0.5); return Normalize(F,L);
       }
       (*this)*=(1./n); return *this;
    };

    wbvector& add(const T* d) {
       for (size_t i=0; i<len; ++i) data[i]+=d[i];
       return *this;
    };

    wbvector& subtract(const T* d) {
       for (size_t i=0; i<len; ++i) data[i]-=d[i];
       return *this;
    };

    wbvector& Conj() { return *this; };
    wbvector& applyFlag(char &flag, const T &fac, wbvector &a);

    bool allEqual() const {
       for (size_t i=1; i<len; ++i) if (data[i]!=data[0]) return 0;
       return 1;
    }

    bool allZero() const { return allEqual(T(0)); };
    bool allOnes() const { return allEqual(T(1)); };

    bool allEqual(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]!=a) return 0;
       return 1;
    }

    bool allEqual2(const T &a, T eps=1E-15) const {
       if (eps==0) return allEqual(a);
       else {
          T x, x2=0; if (a!=0)  eps*=a;
          if (!len) wblog(FL,"WRN %s() got len=%d",FCT,len);
          for (size_t i=0; i<len; ++i) { x=data[i]-a; x2+=(x*x); }
          return (ABS(x2)<eps);
       }
    };

    bool allUnequal(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]==a) return 0;
       return 1;
    }

    bool any() const {
       for (size_t i=0; i<len; ++i) if (data[i]!=0) return 1;
       return 0;
    }

    bool anyUnequal(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]!=a) return 1;
       return 0;
    }

    bool anyEqual(const T &a, size_t *I=NULL) const {
       for (size_t i=0; i<len; ++i) {
           if (data[i]==a) { if (I) (*I)=i; return 1; }
       }
       return 0;
    }

    bool anyGT(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]>a) return 1;
       return 0; }
    bool anyGE(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]>=a) return 1;
       return 0;
    };

    bool anyST(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]<a) return 1;
       return 0; };
    bool anySE(const T &a) const {
       for (size_t i=0; i<len; ++i) if (data[i]<=a) return 1;
       return 0;
    };

    bool allIn(const T &a1, const T &a2) const {
       for (size_t i=0; i<len; ++i) {
          if (data[i]<a1 || data[i]>a2) return 0; }
       return 1;
    };

    bool anySmallVals(const T eps=1E-15) const {
       size_t i=0; if (eps>0)
            { for (; i<len; ++i) if (ABS(data[i])<eps) return 1; }
       else { for (; i<len; ++i) if (data[i]) return 1; }
       return 0;
    };

    bool isEqual(const T* v) const {
       return (!memcmp(data, v, len*sizeof(T)));
    }

    bool anyEqual(const wbvector<size_t> &I, const T &x,
       size_t *k=NULL) const;
    bool anyUnequal(const wbvector<size_t> &I, const T &x,
       size_t *k=NULL) const;

    bool isEmpty() const {
       if (data==NULL) {
          if (len) wblog(FL,"ERR severe data inconsistency (%d,0x%lX)",len,data);
          return 1;
       }
       else if (len==0)
       wblog(FL,"WRN data[] not initialized !?? (%d, 0X%lX)",len,data);

       return 0;
    };

    bool isEqual (const wbvector &v) const {
       if (len!=v.len) return 0;
       for (size_t i=0; i<len; i++) if (data[i]!=v.data[i]) return 0;
       return 1;
    };

    bool operator==(const wbvector &v) const {
       if (this!=&v) {
          if (len!=v.len) return 0;
          if (data!=v.data) {
             for (size_t i=0; i<len; ++i)
             if (data[i]!=v.data[i]) return 0;
          }
       }
       return 1;
    };

    bool operator^(const wbvector &v) const {
       if (this==&v) { return (len ? 0 : 1); }
       if (len!=v.len) return 0;

       for (size_t i=0; i<len; ++i) {
       if (!((data[i]==0) ^ (v.data[i]==0))) return 0; }
       
       return 1;
    };

    bool operator!=(const wbvector &v) const { return !((*this)==v); };

    bool deepEqualP(const wbvector &B) const;

    bool operator< (const wbvector &v) const;
    bool operator<=(const wbvector &v) const;
    bool operator> (const wbvector &v) const { return !((*this)<=v); };
    bool operator>=(const wbvector &v) const { return !((*this)< v); };

    int operator< (const T &r) const {
       for (size_t i=0; i<len; i++) if (data[i]>=r) return 0;
       return (len ? 1 : 0);
    };

    T& el(size_t i) const {
       if (i>=len) wblog(FL,
          "ERR element index out of bounds (%d/%d)",i+1,len);
       return data[i];
    };

    T elx(size_t i, const T &x) const {
       return (i<len ? data[i] : x);
    };

    const T& operator[] (size_t i) const { return data[i]; };
          T& operator[] (size_t i)       { return data[i]; };

    T& last() const { return end(); }
    T& end() const {
       if (len==0) wblog(FL,
          "ERR %s() index out of bounds (len=%d)",FCT,len);
       return data[len-1];
    };

    T& p2last() const {
       if (len==0) wblog(FL,"ERR No last element in [].");
       return data+len-1;
    };

    T* ref(size_t i) const { return (data+i); };

    wbvector& setRand(double fac=1., double shift=0.);
    void setIdx (double shift=0., double fac=1.);

    int set2Group(const wbperm &P,
        const WBINDEX &D, const T* S0=NULL, char iflag=0
    );

    wbvector& RESIZE(size_t n) {
       if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

       if (n) {
          if (n<=len && n/double(len)>0.10) len=n;
          else Resize(n);
       }
       else init();

       return *this;
    };

    wbvector& Resize(size_t n, const T* d=NULL) {
       if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

       if (n==len) return *this;
       if (!n) { len=n;
          if (data) { WB_DELETE(data); }
          return *this;
       }

       T *d0=data; WB_NEW(data,n);
       MEM_CPY<T>(data, n, (n<len ? n : len), d0, d);

       len=n;
       
       WB_DELETE(d0); return *this;
    };

    wbvector& resize(size_t n, wbvector& X) const {
       if (!n) { X.init(); return X; }
       if (n==len) { X=(*this) ; return X; }

       WB_NEW(X.data,n);
       if (len) MEM_CPY<T>(X.data, MIN(n,len), data);
       if (n>len) MEM_SET<T>(X.data+len,n-len);

       return X;
    };

    wbvector& cat(const wbvector &v2) {
       Resize(len+v2.len, v2.data);
       return *this;
    };
    wbvector& Cat(const wbvector &v1, const wbvector &v2) {
       RENEW(v1.len+v2.len);
       MEM_CPY<T>(data,len,v1.len,v1.data,v2.data);
       return *this;
    };
    wbvector& Cat(const wbvector<wbvector  > &vv);
    wbvector& Cat(const wbvector<wbvector* > &vv);

    template <class T1, class T2, class T3>
    wbvector& Cat(
       const wbvector<T1> *v1,const wbvector<T2> *v2,
       const wbvector<T3> *v3=NULL
    ){
       if ((void*)v1==(void*)this || (void*)v2==(void*)this ||
           (void*)v3==(void*)this ){
          wbvector X; X.Cat(v1,v2,v3);
          X.save2(*this); return *this;
       }

       size_t i=0, j=0,
          l1=(v1? v1->len:0), l2=(v2? v2->len:0), l3=(v3? v3->len:0);
       RENEW(l1+l2+l3);

       for (; j<l1; ++j, ++i) { data[i]=T(v1->data[j]); }; j=0;
       for (; j<l2; ++j, ++i) { data[i]=T(v2->data[j]); }; j=0;
       for (; j<l3; ++j, ++i) { data[i]=T(v3->data[j]); }

       return *this;
    };

    wbvector& Append(const wbvector &v2) {
       Resize(len+v2.len, v2.data);
       return *this;
    };
    wbvector& Append(const T &x) {
        Resize(len+1, &x);
        return *this;
    };
    wbvector& Append(size_t d2, const T *x) {
        Resize(len+d2, x);
        return *this;
    };

    wbvector append(const wbvector &v) const {
        wbvector v2; append(v.data,v.len,v2);
        return v2;
    };

    wbvector& append(const wbvector &v2, wbvector &vout) const {
        append(v2.data, v2.len, vout);
        return vout;
    };
    void append(const T *x, size_t d2, wbvector &vout) const {
        vout.RENEW(len+d2);
        MEM_CPY<T>(vout.data,vout.len,len,data,x);
    };

    int recLess (const T* v) const;
    int recLessE(const T* v) const;
    int recCompare(const T* v) const;

    T max(size_t *k=NULL) const;
    T max(unsigned &k) const { T x;
      size_t k_; x=max(&k_); safeConvert(FL,k_,k);
      return x;
    };

    T min(size_t *k=NULL) const;
    T min(unsigned &k) const { T x;
      size_t k_; x=min(&k_); safeConvert(FL,k_,k);
      return x;
    };

    T aMax(size_t *k=NULL) const;

    T aMax(unsigned &k) const { T x;
       size_t k_; x=aMax(&k_); safeConvert(FL,k_,k);
       return x;
    };

    T aMin(char zflag=0, size_t *k=NULL) const;

    T maxneg(size_t *k=NULL) const;
    T minpos(size_t *k=NULL) const;

    T avg() const;
    T std() const;

    void range(T &xmin, T&xmax) const;
    T range() const;

    wbvector& Permute(const wbperm &P, char iflag=0);
    wbvector& permute(wbvector &, const wbperm &P, char iflag=0) const;

    wbvector permute(const wbperm &P, char iflag=0) const {
       wbvector<T> X; permute(X,P,iflag);
       return X;
    };

    wbvector& blockPermute(
       const wbperm &P, wbvector &, char iflag=0) const;

    wbvector& BlockPermute(const wbperm &P, char iflag=0) {
       wbvector<T> x(*this);
       return x.blockPermute(P,*this,iflag);
    };

    int blockCompare(
       const wbperm &pa, const wbvector &b, const wbperm &pb) const;

    void get(const WBINDEX &I0, wbvector &v) const {
       const INDEX_T *I=I0.data;
       size_t i=0, n=I0.len; v.init(n);

       for (; i<n; i++) {
          if (I[i]>=len) wblog(FL,
          "ERR index out of bounds (%d/%d) ???",I[i],len);
          v.data[i]=data[I[i]];
       }
    };

    void get(const wbvector<size_t> &I, T* v) const {
        for (size_t i=0; i<I.len; i++) {
           if (I[i]>=len)
           wblog(FL,"ERR index out of bounds (%d/%d) ???", I[i],len);

           v[i]=data[I[i]];
        }
        return;
    };

    void get(size_t i1, size_t i2, wbvector &v) const {

        if (i1>=len || i2>=len) wberror(FL,"index out of range.");
        if (i2<i1)
             v.init();
        else v.init(i2-i1, data+i1);
    };

    void selectSU(const WBINDEX &I);
    void select(const WBINDEX &I, T* r) const;
    wbvector& select(const WBINDEX &I, wbvector &a) const;
    wbvector& Select(const WBINDEX &I);
    wbvector& BlockSelect(const WBINDEX &I, const WBINDEX &D);

    wbvector select(const WBINDEX &I) const {
       wbvector dI; select(I,dI);
       return dI;
    };

    wbvector& getI(size_t k, wbvector &v) const;
    wbvector& getI(const wbindex &k, wbvector &v, char uflag=1) const;

    wbvector& getI(const WBINDEX &k, wbvector &v, char uflag=1) const {
       return getI((const wbindex&)k,v,uflag);
    };

    wbvector& initI(const wbvector &a, const wbindex &ia){
       return a.getI(ia,*this); };
    wbvector& initI(
       const wbvector &a, const wbindex &ia,
       const wbvector &b, const wbindex &ib,
       char uflag=1
    );

    wbvector& Flip() {
        size_t i, m2=len/2, m=len-1;
        for (i=0; i<m2; i++) SWAP(data[i],data[m-i]);
        return *this;
    };

    void flip(wbvector &a) const {
        if (&a==this) { a.Flip(); return; }
        a.init(len);
        for (size_t m=len-1, i=0; i<len; i++) a.data[i]=data[m-i];
    };

    wbvector flip() const {
       wbvector a; flip(a);
       return a;
    };


    wbvector operator[] (const wbvector<size_t>& I) const {  
        wbvector v(I.len);
        size_t i=0, *const k=I.data;

        for (; i<I.len; ++i) {
           if (k[i]>=len) wblog(FL,
              "ERR index out of bounds (%d: %d/%d)",i,I.max(),len);
           v[i]=data[k[i]];
        }

        return v;
    };

    T prod() const {
        T p=0;

        if (len>0) {
            p=data[0];
            for (size_t i=1; i<len; i++) p*=data[i];
        }
        else wblog(FL,"WRN %s() empty data set",FCT);

        return p;
    }

    T prod(const T& p0) const {
        T p; if (len==0) return p0;

        p=data[0]; for (size_t i=1; i<len; i++) p*=data[i];
        return p;
    }

    T prod(size_t i1, size_t i2) const {
        size_t i; T p=0;

        if (i1<len && i2<len && i1<=i2) {
            for (p=data[i1], i=i1+1; i<len && i<=i2; i++)
            p*=data[i];
        }
        else wblog(FL, "ERR empty index set (%d,%d/%d)",i1,i2,len);
        return p;
    };

    T prod(const INDEX_T *I, size_t n) const {
        T p=0;
        if (n) { p=1;
           for (size_t k,i=0; i<n; ++i) { k=I[i];
              if (k>=len) wblog(FL,
                 "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
              p*=data[k];
           }
        }
        return p;
    };

    template <class T1, class T2>
    T2 prod(const T1 *I, size_t n, T2 &s2) const {
       size_t i,k;
       char m[len]; for (i=0; i<len; ++i) m[i]=0;
       for (i=0; i<n; ++i) { k=I[i];
          if (k>=len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
          if ((++m[k])>1) wblog(FL,"ERR %s() index not unqiue",FCT);
       }
       if (!len) { s2=0; return 0; }
       else {
          T2 s1=1; s2=1;
          for (i=0; i<len; ++i) { (m[i] ? s1 : s2) *= data[i]; }
          return s1;
       }
    };

    void prod2(T &p1, T &p2, const T &p0=0) const {
        size_t l=len/2; if (len%2) wblog(FL,
        "ERR cannot split vector into two equal parts (%d)",len); 
       
        if (len==0) { p1=p0; p2=p0; return; }
        p1=data[0]; p2=data[l];
        
        for (size_t i=1; i<l; i++) { p1*=data[i]; p2*=data[i+l]; }
    }

    void info(const char *istr) const;
    void print(const char *istr, char mflag=0) const;
    void printdata(const char *istr=0) const {
       if (istr && istr[0]) printf("%s = ",istr);
       printf("[ %s ]\n",toStr().data);
    };

    void put(
      const char *F, int L,
      const char *vname="ans", char tflag=0, const char *ws="base"
    ) const;

    void put(const char *vname="ans", char tflag=0, const char *ws="base"
    ) const { put(0,0,vname,tflag,ws); }

    void dispWithMsg(const char* file, int line,
    const char *msg, const char* vname="") const {
        print(vname);
        wblog(file,line,msg);
    }

    void mat2mx(mxArray* &a, char dim=2) const;
    void mat2mxs(mxArray* a, char field_nr, char dim=2) const;

    mxArray* toMx_Struct (const char tflag=0) const;
    mxArray* toMx_StructP(const char tflag=0) const;

    mxArray* toMx (const char tflag=0) const;
    mxArray* toMxP(const char tflag=0) const;
    mxArray* toMx_base  (const char tflag=0) const;

    mxArray* toMx_offset(const T& offset, const char tflag=0) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tflag=0) const;
    void add2MxStruct(mxArray *S, const char *vname, char tflag=0) const;

    mxArray* mxCreateCell(unsigned m, unsigned n) const;
    void add2MxCell(mxArray *S, unsigned i, char tflag=0) const;

    void putg(const char *vname="ans", char tflag=0) const {
        put(vname,tflag,"global");
        return;
    };

    void putx(const char *vname="ans", char tflag=0) const {
        put(vname,tflag,"base");
        return;
    };

    wbstring toStr(int n=-1, const char *sep=" ") const;
    wbstring toStrD() const;
    wbstring toStrf(
       const char *fmt0="",
       const char *sep=" ",
       const char stride=0,
       const char *sep2=""
    ) const;

    T *data;

    size_t len;

    char isref;

  protected:

    void NEW(size_t n, const T* d=NULL) {
        len=n; data=NULL; isref=0;
        if (len) { WB_NEW(data,len); MEM_CPY<T>(data,len,d); }
    }

    wbvector& RENEW(size_t n, const T* d=NULL) {
        if (isref) { data=NULL; isref=0; } else
        if (n!=len && data) { WB_DELETE(data); }

        len=n; if (len) {
           if (!data) { WB_NEW(data,len); }
           MEM_CPY<T>(data,len,d);
        }
        return *this;
    };

    wbvector& RENEW_VAL(size_t n, const T& x) {
        if (isref) { data=NULL; isref=0; } else
        if (n!=len && data) { WB_DELETE(data); }

        len=n;
        if (len) {
           if (!data) { WB_NEW(data,len); }
           for (size_t i=0; i<len; ++i) data[i]=x;
        }
        return *this;
    };

  private:

};


template <class T>
wbvector<T>& wbvector<T>::unRef() {

    if (!isref) return *this;
    if (!len) {
       if (data) wblog(FL,"ERR data=%lX (%d)",data,len);
       isref=0; return *this;
    }

    T const* const d0=data;

    WB_NEW(data,len); MEM_CPY<T>(data,len,d0);
    isref=0;

    return *this;
};


template <class T>
bool wbvector<T>::deepEqualP(const wbvector<T> &B) const {
    if (len!=B.len) return 0;

    if (data!=B.data) {
       for (size_t i=0; i<len; i++) { if (data[i]==B.data[i]) continue;
          if ((!data[i] || !B.data[i])) return 0; else
          if (!((*data[i])==(*B.data[i]))) return 0;
       }
    }

    return 1;
};


template <class T> inline 
bool wbvector<T>::operator<(const wbvector<T> &v) const {

   const T *b=v.data;
   for (size_t n=MIN(len,v.len), i=0; i<n; ++i) {
      if (data[i]<b[i]) return 1;
      if (data[i]>b[i]) return 0;
   }
   if (len<v.len) return 1; else return 0;
};

template <class T> inline 
bool wbvector<T>::operator<=(const wbvector<T> &v) const {

   const T *b=v.data;
   for (size_t n=MIN(len,v.len), i=0; i<n; i++) {
      if (data[i]<b[i]) return 1;
      if (data[i]>b[i]) return 0;
   }
   if (len<=v.len) return 1; else return 0;
};

template <class T> inline 
bool wbvector<T>::checkSameLength(
   const wbvector<T> &a,
   const char *F, int L, const char *istr
) const {

   bool i=(len==a.len);
   if (!i && F) {
      if (istr && istr[0]) sprintf(str,"%s() - ",istr); else str[0]=0;
      wblog(F,L, "ERR %svector length mismatch (%d/%d)",str,len,a.len);
   }
   return i;
};

template <class T> inline 
bool wbvector<T>::checkSameLength(
   const wbvector<T> &a,
   const wbvector<T> &b,
   const char *F, int L, const char *istr
) const {

   bool i=(len==a.len && len==b.len);
   if (!i && F) {
      if (istr && istr[0]) sprintf(str,"%s() - ",istr); else str[0]=0;
      wblog(F,L,"ERR %svector length mismatch (%d,%d/%d)",str,a.len,b.len,len);
   }
   return i;
};


template <> inline 
wbvector<wbcomplex>& wbvector<wbcomplex>::Conj() {
   for (size_t i=0; i<len; i++) data[i].Conj();
   return *this;
};

template <class T> inline 
wbvector<T>& wbvector<T>::applyFlag(
   char &flag, const T &fac, wbvector<T> &a
){
   if (!strchr("NCTc",flag))
   wblog(FL,"ERR %s() invalid flag %c<%d>",__FUNCTION__,flag,flag);

   if (flag=='c') flag='C';
   if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='N';
   
   if (fac!=1 || flag=='C') { a=(*this);
      if (flag=='C') a.Conj();
      if (fac!=1) a*=fac;
   }
   else a<<*this;
};


template <class T> inline 
wbvector<T>& wbvector<T>::save2(wbvector<T> &v) {


   if (this!=&v) {
      if (v.data) v.init();
      v.data=data; v.len=len; v.isref=isref;
      data=NULL; len=0; isref=0;
   }

   return v;
};


template <class T> inline
T wbvector<T>::min(size_t *k) const {
   T x=0; size_t i=1; if (k) (*k)=0;

   if (!len) wblog(FL,"WRN %s() got empty vector",FCT);
   else {
      for (x=data[0]; i<len; ++i) { if (x>data[i]) {
           x=data[i]; if (k) (*k)=i;
      }}
   }
   return x;
};


template <class T> inline
T wbvector<T>::max(size_t *k) const {
   T x=0; size_t i=1; if (k) (*k)=0;

   if (!len) wblog(FL,"WRN %s() got empty vector",FCT);
   else {
      for (x=data[0]; i<len; ++i) { if (x<data[i]) {
           x=data[i]; if (k) (*k)=i; 
      }}
   }
   return x;
};


template <class T> inline
T wbvector<T>::minpos(size_t *k) const {
   T x=0; bool found=0; if (k) (*k)=-1;

   if (len==0) wblog(FL,"WRN %s() of null vector",FCT);
   else {
      for (size_t i=0; i<len; i++)
      if (data[i]>0) {
         if (found) { if (x>data[i]) { x=data[i]; if (k) (*k)=i; }}
         else { x=data[i]; if (k) (*k)=i; found=1; }
      }
   }
   return x;
};


template <class T> inline
T wbvector<T>::maxneg(size_t *k) const {
   T x=0; bool found=0; if (k) (*k)=-1;

   if (len==0) wblog(FL,"WRN %s() of null vector",FCT);
   else {
      for (size_t i=0; i<len; i++)
      if (data[i]<0) {
         if (found) { if (x<data[i]) { x=data[i]; if (k) (*k)=i; }}
         else { x=data[i]; if (k) (*k)=i; found=1; }
      }
   }

   return x;
};


template <class T> inline 
T wbvector<T>::aMax(size_t *k_) const {

   if (len) {
      size_t k=-1, i=0; T a,x=0;
      for (; i<len; ++i) {
         a=ABS(data[i]); if (x<a) { x=a; k=i; }
      }
      if (k_) (*k_)=k;
      return x;
   }
   else {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got null vector",FCT); }
      return 0;
   }
};


template <class T> inline 
T wbvector<T>::aMin(char zflag, size_t *k_) const {

   if (!len) {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got null vector",FCT); }
      return 0;
   }

   size_t i=0, k=-1; double a,x=0;

   if (zflag==0) {
      x=ABS(data[i++]); k=0; if (x!=0) {
         for (; i<len; ++i) { a=ABS(data[i]); if (x>a) {
            x=a; k=i; if (x==0) break;
         }}
      }
   }
   else {
      for (; i<len; ++i) {
         a=ABS(data[i]); if (a) { x=a; k=i; break; }
      }
      for (; i<len; ++i) { a=ABS(data[i]);
         if (a && x>a) { x=a; k=i; }
      }
   }
   if (k_) (*k_)=k;
   return x;
};


template <class T> inline 
void wbvector<T>::range(T &xmin, T&xmax) const {
   if (len==0) {
       xmin=xmax=0;
       wblog(FL,"WRN %s() of null vector",FCT);
   }
   else {
      xmin=xmax=data[0];
      for (size_t i=1; i<len; i++) {
         if (xmin>data[i]) xmin=data[i];
         if (xmax<data[i]) xmax=data[i];
      }
   }
};


template <class T> inline 
T wbvector<T>::range() const {
   T xmin, xmax; range(xmin,xmax);
   return (xmax-xmin);
};


template <class T> inline 
T wbvector<T>::avg() const {

   if (!len) {
      wblog(FL,"WRN %s() of null vector !?",FCT);
      return 0;
   }
   else {
      T x=data[0];
      for (size_t i=1; i<len; ++i) { x+=data[i]; }
      return (x/T(len));
   }
};

template <class T> inline 
T wbvector<T>::std() const {

   if (len<2) {
      wblog(FL,"WRN %s() of vector of length %d !?",FCT,len);
      return 0;
   }
   else {
      T x=data[0], x2=data[0]*data[0];
      for (size_t i=1; i<len; ++i) { x+=data[i]; x2+=data[i]*data[i]; }
      return SQRT( ( x2 - (x*x)/T(len) ) / T(len-1) );
   }
};


template <class T> inline 
int wbvector<T>::findClosest(T r) const {
   if (len==0) return -1;

   T dbl, eps=ABS(data[0]-r); size_t i,k=0;

   for (i=1; i<len; i++) {
      dbl = ABS(data[i]-r);
      if (dbl<eps) { k=i; eps=dbl; if (eps==0) break; }
   }

   return k;
}

template <class T>
int wbvector<T>::findClosestSorted(T r) const {
   if (len==0) return -1;

   size_t i, k, n1=0, n2=len, e=0;
   char c, cref;

   if (len==0) return -1;
   if (len==1) return 0;

   cref=NUMCMP(data[0], data[len-1]);

   if (cref) {
      if (len>6) {
        for (i=1; i<4; i++) if (NUMCMP(data[i-1],data[i])==-cref) e++;
        for (i=len-3; i<len; i++)
        if (NUMCMP(data[i-1],data[i])==-cref) e++;
      }
      else
      for (i=1; i<len; i++) if (NUMCMP(data[i-1],data[i])==-cref) e++;
   }
   else for (i=1; i<len; i++) if (NUMCMP(data[i-1],data[i])) e++;

   if (e) wblog(FL,"ERR Data not sorted! (%d)", e);
   if (cref==0) { wblog(FL,
      "WRN find closest: all data constant!??\n(%g; %g %dx)",
      (double)r, (double)data[0], len); return 0;
   }

   c=NUMCMP(data[len-1],r); if (c!=-cref) return len-1;
   c=NUMCMP(data[0],r);     if (c!= cref) return 0;

   for (k=0;;) {
      if (c==cref) {
          n1=n2/2; n2-=n1; if (n1) k+=n1; else {
             if (k+1<len && (ABS(r-data[k]) > ABS(r-data[k+1]))) k++;
             break;
          }
      }
      else if (c) {
          n2=n1/2; n1-=n2; if (n2) k-=n2; else {
             if (k) if (ABS(r-data[k]) > ABS(r-data[k-1])) k--;
             break;
          }
      }
      else break;

      c=NUMCMP(data[k],r);
   }

   return k;
}


template <class T>
wbvector<T>& wbvector<T>::setRand(double fac, double shift) {
   size_t i;
   static char first_call=1;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX;

   if (shift==0.)
        for (i=0; i<len; ++i) data[i]=(T)(fac*rand());
   else for (i=0; i<len; ++i) data[i]=(T)(fac*rand()+shift);

   return *this;
}


template<>
wbvector<wbcomplex>& wbvector<wbcomplex>::setRand(double fac, double shift) {
   size_t i;
   static char first_call=1;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX;

   if (shift==0.) {
      for (i=0; i<len; ++i)
      data[i].set(fac*rand(), fac*rand());
   }
   else {
      for (i=0; i<len; ++i)
      data[i].set(fac*rand()+shift, fac*rand()+shift);
   }

   return *this;
}


template <class T> inline 
void wbvector<T>::setIdx(double shift, double fac) {
   size_t i;

   if (shift==0.) {
        if (fac==1.)
             for (i=0; i<len; i++) data[i]=(T)(i);
        else for (i=0; i<len; i++) data[i]=(T)(fac*i);
   }
   else {
        if (fac==1.)
             for (i=0; i<len; i++) data[i]=(T)(i+shift);
        else for (i=0; i<len; i++) data[i]=(T)(fac*i+shift);
   }
}


template <class T>
void wbvector<T>::mat2mx(mxArray* &a, char dim) const {

    double *d;

    if (dim==1) a=mxCreateDoubleMatrix(len,1,mxREAL); else
    if (dim==2) a=mxCreateDoubleMatrix(1,len,mxREAL); else
    wberror(FLINE, "Invalid dim parameter.");

    d=mxGetPr(a);
    for (size_t i=0; i<len; i++) d[i]=(double)data[i];
}


template <class T>
void wbvector<T>::mat2mxs(mxArray* a, char fid, char dim) const {

    size_t i, m=1, n=len;
    double *d;

    if (!mxIsStruct(a))
    wberror(FLINE, "Need structure input.");

    if (mxGetNumberOfElements(a)>1)
    wberror(FLINE, "Need SINGLE structure on input.");

    if (fid>=mxGetNumberOfFields(a))
    wberror(FLINE, "index of fields out of bounds.");

    if (dim==1) { m=len; n=1; } else
    if (dim==2) { m=1; n=len; } else
    wberror(FLINE, "Invalid dim parameter.");

    mxSetFieldByNumber(a,0,fid, mxCreateDoubleMatrix(m,n,mxREAL));
    d=mxGetPr(mxGetFieldByNumber(a,0,fid));

    for (i=0; i<len; i++) d[i]=(double)data[i];
}


template <class T> inline 
int wbvector<T>::recLess(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return 1; else
    if (data[i]>v[i]) return 0;

    return 0;
}

template <class T> inline 
int wbvector<T>::recLessE(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return 1; else
    if (data[i]>v[i]) return 0;

    return 1;
}

template <class T> inline 
int wbvector<T>::recCompare(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return -1; else
    if (data[i]>v[i]) return +1;

    return 0;
}


template<class T>
mxArray* wbvector<T>::toMx_base(const char tflag) const {

    size_t i; double *dd;

    mxArray *a;

    if (isEmpty()) { return mxCreateDoubleMatrix(0, len, mxREAL); }

    if (tflag)
         a=mxCreateDoubleMatrix(len, 1, mxREAL);
    else a=mxCreateDoubleMatrix(1, len, mxREAL);

    dd=mxGetPr(a);
    for (i=0; i<len; i++) dd[i]=double(data[i]);

    return a;
}


template<>
mxArray* wbvector<wbcomplex>::toMx_base(const char tflag) const {

    mxArray *a;
    size_t i=0, m=1, n=len, isreal=1;

    if (tflag) { m=n; n=1; }
    for (; i<len; ++i) { if (data[i].i!=0) { isreal=0; break; }}

    a=mxCreateDoubleMatrix(m, n, isreal ? mxREAL : mxCOMPLEX);
    if (!a) wblog(FL,"ERR got a==null !?");

    if (len) {
       double *dd=mxGetPr(a); if (!dd) wblog(FL,"ERR got pr==null !?");
       for (i=0; i<len; ++i) dd[i]=double(data[i].r);

       if (!isreal) {
          dd=mxGetPi(a); if (!dd) wblog(FL,"ERR got pi==null !?");
          for (i=0; i<len; ++i) dd[i]=double(data[i].i);
       }
    }

    return a;
};


template<>
mxArray* wbvector<long>::toMx_base(const char tflag) const {

    size_t i=0, m=1, n=len; if (tflag) { m=n; n=1; }

    mxArray *a=mxCreateNumericMatrix(m,n, mxINT64_CLASS, mxREAL);
    if (!a) wblog(FL,"ERR got a==null !?");

    if (len) {
       int64_T *dd=(int64_T*)mxGetData(a);
          if (!dd) wblog(FL, "ERR got pr==null !?");
       for (; i<len; ++i) dd[i]=data[i];
    }

    return a;
};

template<>
mxArray* wbvector<unsigned long>::toMx_base(const char tflag) const {

    size_t i=0, m=1, n=len; if (tflag) { m=n; n=1; }

    mxArray *a=mxCreateNumericMatrix(m,n,mxUINT64_CLASS, mxREAL);
    if (!a) wblog(FL,"ERR %s() got a=null !?",FCT); 

    if (len) {
       uint64_T *dd=(uint64_T*)mxGetData(a); 
          if (!dd) wblog(FL,"ERR got pr==null !?");
       for (; i<len; ++i) dd[i]=data[i];
    }

    return a;
};


template<class T> inline 
mxArray* wbvector<T>::toMx(const char tflag) const {
   return toMx_Struct(tflag);
};

template<class T> inline 
mxArray* wbvector<T>::toMxP(const char tflag) const {
   return toMx_StructP(tflag);
};


template<> inline 
mxArray* wbvector<unsigned>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<int>::toMx(const char tflag) const
{ return toMx_base(tflag); }

template<> inline 
mxArray* wbvector<short>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<unsigned short>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<long>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<unsigned long>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<long long>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<unsigned long long>::toMx(const char tflag) const
{ return toMx_base(tflag); };

template<> inline 
mxArray* wbvector<double>::toMx(const char tflag) const
{ return toMx_base(tflag); }

template<> inline 
mxArray* wbvector<char>::toMx(const char tflag) const
{ return toMx_base(tflag); }

template<> inline 
mxArray* wbvector<unsigned char>::toMx(const char tflag) const
{ return toMx_base(tflag); }

template<> inline 
mxArray* wbvector<wbcomplex>::toMx(const char tflag) const
{ return toMx_base(tflag); }

template<> inline 
mxArray* wbvector<long double>::toMx(const char tflag) const
{ return toMx_base(tflag); }

#ifdef __WB_MPFR_HH__




#endif


template<class T>
mxArray* wbvector<T>::toMx_Struct(const char tflag) const {

    mxArray *S;

    if (tflag)
         S=data->mxCreateStruct(len,1);
    else S=data->mxCreateStruct(1,len);

    for (size_t i=0; i<len; ++i) data[i].add2MxStruct(S,i);
    return S;
}


template<class T>
mxArray* wbvector<T>::toMx_StructP(const char tflag) const {

    mxArray *S;

    if (tflag)
         S=(*data)->mxCreateStruct(len,1);
    else S=(*data)->mxCreateStruct(1,len);

    for (size_t i=0; i<len; i++) data[i]->add2MxStruct(S,i);
    return S;
};


template<class T>
mxArray* wbvector<T>::toMx_offset(const T& offset, const char tflag)
const {
    if (!len) { return mxCreateDoubleMatrix(0,0,mxREAL); }

    size_t i; double *dd;

    mxArray *a=(tflag ?
       mxCreateDoubleMatrix(len,1,mxREAL) :
       mxCreateDoubleMatrix(1,len,mxREAL)
    );

    dd=mxGetPr(a);

    for (i=0; i<len; i++)
    dd[i]=double(data[i]+offset);

    return a;
}

template <class T>
void wbvector<T>::put(
    const char *file, int line, const char *vname, char tflag,
    const char* ws
) const {

    size_t i;
    mxArray *a=toMx(tflag);

    if (file) wblog(file,line,
    "I/O putting vector '%s' to %s%s",vname,ws,tflag ? " (t)" : ".");

    i=mexPutVariable(ws, vname, a);

    if (i==1) {
       info("this"); wblog(__FL__,
      "ERR failed to put variable %s into workspace %s.",vname,ws);
    }

    mxDestroyArray(a);
    return;
}


inline WBINDEX Index(size_t i1, size_t i2) {
    WBINDEX ii;

    if (i1>i2) return ii;

    ii.init(i2-i1+1);
    for (size_t i=0; i<ii.len; i++) ii[i]=i+i1;

    return ii;
};


inline char isUniqueIdxSet(const WBINDEX &I0, INDEX_T imax) {

    size_t i, n=I0.len, e=0;
    const INDEX_T *I=I0.data;
    wbvector<char> mark(imax); char *m=mark.data;

    for (i=0; i<n; i++) if (I[i]<imax) m[I[i]]++; else e++;
    for (i=0; i<imax;  i++) if (m[i]>1) e++;

    return (e==0);
}


inline char validPerm(const WBPERM &P, PERM_T d) {

    size_t i,n=P.len;
    wbvector<char> f(n);

    if (sPERM_T(d)>=0) if (P.len!=d) return 0;

    for (i=0; i<n; i++) {
        if (P[i]>=n) return 0;
        f[P[i]]++;
    }

    for (i=0; i<n; i++) if (f[i]!=1) return 0;

    return 1;
}


inline void getIPerm(const WBPERM &P, WBPERM &iP) {

    size_t i,n=P.len;

  #ifdef __WBDEBUG__

    if (!validPerm(P))
    wberror(FLINE,"Invalid permutation.");

  #endif

    iP.init(n);
    for (i=0; i<n; i++) iP[P[i]]=i;

    return;
}


char isIdentityPerm(const WBPERM &P) {

    for (size_t i=0; i<P.len; i++) if (P[i]!=i) return 0;
    return 1;
}


#endif

