#ifndef __WB_INDEX_HH__
#define __WB_INDEX_HH__

// ----------------------------------------------------------------- //
// index class
// Wb,Apr19,06
// ----------------------------------------------------------------- //

// inherit class wbvector<unsigned>

class wbindex : public wbvector<INDEX_T> {

  public:

    wbindex(INDEX_T n=0, char iflag=0) {
       if (iflag) Index(n); else init(n);
    };

    wbindex(const ctrIdx& ic) { init(ic); };

    explicit wbindex(INDEX_T n, const INDEX_T* d)
    : wbvector<INDEX_T>(n,d) {};

    explicit wbindex(const char* s, unsigned offset=0) {
       int i=Str2Idx(s,*this,offset);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
    };

    wbindex& init(INDEX_T n=0, const INDEX_T *d=NULL) {
       wbvector<INDEX_T>::init(n,d);
       return *this;
    };

    wbindex& init(const ctrIdx &ic);

    wbindex& init(
       const char *F, int L, const mxArray *a, INDEX_T offset=0);

    INDEX_T init(
       const char *F, int L, const mxArray *a,
       const char *istr, char check_type=1
    ){ return wbvector<INDEX_T>::init(F,L,a,istr,check_type); };

    wbindex& init(const mxArray *a, INDEX_T offset=0) {
       return init(0,0,a,offset);
    };

    wbindex& initGroup(
       const wbvector<INDEX_T> &P, wbvector<INDEX_T> &D
    ){
       INDEX_T i,d, k=0, ig=0, m=D.sum();
       if (m!=P.len) wblog(FL,
          "ERR %s() invalid group index (%d/%d)",FCT,m,P.len);
       init(P.len+D.len);
    
       for (; ig<D.len; ++ig) {
          data[k+ig]=d=D[ig]; m=ig+1;
          for (i=0; i<d; ++i, ++k) data[k+m]=P[k];
       }

       return *this;
    };

    wbindex& initGroup0(
       const wbvector<INDEX_T> &P, wbvector<INDEX_T> &D
    ){
       INDEX_T l=0, i=0; init(D.len);
       for (; i<D.len; ++i) {
          if (l>=P.len) wblog(FL,"ERR %s() index out of bounds",FCT,l,P.len);
          data[i]=P[l]; l+=D[i];
       }
       return *this;
    };


    void initStr(const char* s, INDEX_T offset=0) {
       int i=Str2Idx(s,*this,offset);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
    };

    bool isIndex(INDEX_T n=0) {
       if (n) {
          if (sINDEX_T(n)>0) {
             for (INDEX_T i=0; i<len; ++i) {
             if (data[i]!=(n+i)) return 0; }
          }
          else {
             n=(1-n-len); if (sINDEX_T(n)<0) return 0;
             for (INDEX_T i=0; i<len; ++i) {
             if (data[i]!=(n+i)) return 0; }
          }
       }
       else {
          for (INDEX_T i=0; i<len; ++i) { if (data[i]!=i) return 0; }
       }
       return 1;
    };

    wbindex& Index(INDEX_T n) {
       init(n);
       if (n) { for (INDEX_T i=0; i<n; ++i) data[i]=i; }
       return *this;
    };

    wbindex& Index_ex(INDEX_T n, INDEX_T ix) {
       INDEX_T i,k;

       if (ix>=n) wblog(FL,
       "WRN index to be excluded of no relevance (%d/%d)",ix,n);

       init(ix<n ? n-1 : n); if (n==0) return *this;

       for (k=i=0; i<n; i++) if (i!=ix) data[k++]=i;

       return *this;
    };

    wbindex& Index(INDEX_T i1, INDEX_T i2) {
       if (i1>i2) { init(); return *this; }
       init(i2-i1+1);
       for (INDEX_T i=0; i<len; i++) data[i]=i+i1;

       return *this;
    };

    wbindex& BlockIndex(const wbvector<INDEX_T> &D, char sub=0) {
       INDEX_T i,j,d,l; init(D.sum());

       if (sub==0) {
          for (l=i=0; i<D.len; ++i)
          for (d=D[i], j=0; j<d; ++j,++l) data[l]=i;
       }
       else {
          for (l=i=0; i<D.len; ++i)
          for (d=D[i], j=0; j<d; ++j,++l) data[l]=j;
       }

       return *this;
    };

    wbindex& BlockIndex(INDEX_T n, INDEX_T D) {
       INDEX_T i,j,l; init(n*D);

       for (l=i=0; i<n; ++i)
       for (  j=0; j<D; ++j,++l) data[l]=i;

       return *this;
    };

    wbindex& BlockIndexI(INDEX_T n, INDEX_T D) {
       INDEX_T i,j,l; init(n*D);

       for (l=j=0; j<D; ++j)
       for (  i=0; i<n; ++i,++l) data[l]=i;

       return *this;
    };

    wbindex& IndexTranspose(INDEX_T m, INDEX_T n) {
       INDEX_T i,j,l;
       init(n*m); if (!n || !m) return *this;

       for (l=j=0; j<m; ++j) {
          data[l++]=j;
          for (i=1; i<n; ++i,++l) data[l]=data[l-1]+m;
       }

       return *this;
    };

    wbindex& operator= (const char* s) {
       int i=Str2Idx(s,*this,0);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
       return *this;
    };

    int extend2Perm(INDEX_T N, wbperm &P, const char *pstr=NULL) const;

    void toPerm(wbperm &P) const {
        P.initT(len, data);

        if (!P.isValidPerm()) { print("this");
        wblog(FL, "WRN No valid perm."); }
    };

    wbindex& invert(INDEX_T N, wbindex &I2, char uflag=1) const;

    wbindex& Invert(INDEX_T N, char uflag=1) {
       wbindex X; save2(X);
       return X.invert(N,*this,uflag);
    };

    wbindex flipIdx(INDEX_T ndim) const;
    char isUnique(INDEX_T imax, INDEX_T offset=0) const;

    mxArray* toMx(const char tflag=0) const {
       return toMx_offset(1,tflag);
    };


  protected:
  private:

};


wbperm::wbperm (const wbindex &P) {
   NEW(P.len,P.data);
   if (!isValidPerm()) dispInvalidPerm(FL);
};


class wbIndex : public wbvector<INDEX_T> {

  public:

    wbIndex() : wbvector<INDEX_T>() {};
    wbIndex(const wbvector<INDEX_T> &S) { init(S); };
    wbIndex(INDEX_T n, const INDEX_T *S) { init(n,S); };

    wbIndex(const wbIndex &I) : wbvector<INDEX_T>() {
       if (!I.isEmpty()) {
          wbvector<INDEX_T>::init(I.len); SIZE=I.SIZE;
          if (len) --data[0];
       }
    };


    wbIndex& init() {
       wbvector<INDEX_T>::init(); SIZE.init();
       return *this;
    };

    wbIndex& init(const wbvector<INDEX_T> &S) {
       const INDEX_T* const &s=S.data;
       for (INDEX_T k=0; k<S.len; k++) if (!s[k]) {
          wbvector<INDEX_T>::init(); SIZE.init();
          return *this;
       }

       wbvector<INDEX_T>::init(S.len); SIZE=S;
       if (len) --data[0];
       return *this;
    };

    wbIndex& reset() {
       if (len!=SIZE.len) wblog(FL,
          "ERR %s() severe size mismatch (%d/%d)",FCT,len,SIZE.len);
       if (len) {
          data[0]=-1;
          for (INDEX_T i=1; i<len; ++i) data[i]=0;
       }
       return *this;
    };

    template<class TI>
    wbIndex& init(unsigned n, const TI *S) {
       for (unsigned k=0; k<n; ++k) if (!S[k]) {
          wbvector<INDEX_T>::init(); SIZE.init();
          return *this;
       }

       wbvector<INDEX_T>::init(n); SIZE.initT(n,S);
       if (len) --data[0];
       return *this;
    };

    wbIndex& set(wbvector<INDEX_T> &I) {
       INDEX_T k=0, *const &i=I.data, *const &s=SIZE.data;
       if (I.len!=len) wblog(FL,
          "ERR wbIndex::%s() size mismatch (%d/%d)",FCT,I.len,len);
       for (; k<len; k++) if (i[k]>=s[k]) {
           wblog(FL,"ERR wbIndex::%s() index out of bounds (%d: %d/%d)",
              FCT,k+1,i[k],s[k]);
           data[k]=i[k];
       }
       return *this;
    };

    bool isValid() const {
       if (len!=SIZE.len) return 0;
       for (INDEX_T i=0; i<len; ++i) if (data[i]>=SIZE.data[i]) return 0;
       return 1;
    };

    void checkValid(const char *F, int L) const {
       if (!isValid()) wblog(F_L,
          "ERR %s() invalid hyperindex\n%s",FCT,toStr().data);
    };

    bool operator++() {
       INDEX_T k=0, l=len-1, *const &s=SIZE.data;

       if (len!=SIZE.len) wblog(FL,
          "ERR wbIndex::%s() got corrupted index (%d/%d)",
           FCT, len, SIZE.len
       );
       if (!len) return 0;

       k=0; data[0]++;
       while (data[k]==s[k] && k<l) { data[k]=0; ++data[++k]; }

       if (data[k]>=s[k]) {
          if (k!=l) wblog(FL,
             "WRN wbIndex::%s() premature stop (%d/%d: %s)",
              FCT, k+1,len, toStr().data);
          return 0;
       }
       return 1;
    };

    bool operator--() {
       INDEX_T k=0, l=len-1, *const &s=SIZE.data;

       if (len!=SIZE.len) wblog(FL,
          "ERR wbIndex::%s() got corrupted index (%d/%d)",
           FCT, len, SIZE.len
       );
       if (!len) return 0;

       k=0; data[0]--;
       while (sINDEX_T(data[k])==-1 && k<l) { data[k]=s[k]-1; --data[++k]; }

       if (data[k]>=s[k]) {
          if (k!=l) wblog(FL,
             "WRN wbIndex::%s() premature stop (%d/%d: %s)",
              FCT, k+1,len, toStr().data);
          return 0;
       }
       return 1;
    };

    bool isEmpty() const {
       if (len!=SIZE.len) wblog(FL,"ERR %d/%d",len,SIZE.len);
       if (!len) return 1;
       else {
          const INDEX_T *s=SIZE.data;
          for (INDEX_T i=0; i<len; i++) if (!s[i]) return 1;
       }
       return 0;
    };

    bool gotmore() const {
       if (len!=SIZE.len) wblog(FL,"ERR %d/%d",len,SIZE.len);
       for (INDEX_T i=0; i<len; ++i) {
          if (data[i]+1<SIZE.data[i]) return 1; else
          if (data[i]>=SIZE.data[i]) wblog(FL,
             "WRN %s() index out of bounds (%s)",FCT,toStr().data
          );
       }
       return 0;
    };

    bool isdiag() const {
       if (len%2) wblog(FL,
          "ERR %s() applies to even-rank objects only (%s)",
           FCT,toStrf("","x").data
       );
       for (INDEX_T r=len/2, i=0; i<r; i++) {
          if (data[i]!=data[i+r]) return 0;
       }
       return 1;
    };

    INDEX_T numel() const { return  SIZE.prod(0); };

    INDEX_T serial(
       const wbvector<INDEX_T> *S2=NULL
    ){
       INDEX_T is, k=len-1, *const &s=(S2 ? S2->data : SIZE.data);
       if (!len) wblog(FL,"ERR wbIndex::%s() got empty index",FCT);
       for (is=data[k--]; k<len; --k) { is = is*s[k]+data[k]; }
       return is;
    };

    wbstring sizeStr() const { return SIZE.toStrf("","x"); };
    wbstring toStr() const {
       char s[128]; snprintf(s,128,"[%s] %s",
          wbvector<INDEX_T>::toStr().data,
          SIZE.toStrf("","x").data);
       return s;
    };

    mxArray* toMx() const {
       const char *fields[]={"S","idx"};
       mxArray *S=mxCreateStructMatrix(1,1,2,fields);

       mxSetFieldByNumber(S,0,0, SIZE.toMx());
       mxSetFieldByNumber(S,0,1,((*this)+1).toMx());
       return S;
    };

    wbvector<INDEX_T> SIZE;

  protected:
  private:

};



template <class T=INDEX_T>
class groupIndex {

  public:

    groupIndex() {};

    groupIndex(
       const wbvector<T> &P_, const wbvector<T> &D_,
       const PERM_T **p=NULL
    ){
       P=P_; D=D_;
       setup(); if (p) (*p)=P.data;
    };

    groupIndex(const wbvector<T> &DI, INDEX_T ng=-1) { init(DI,ng); }

    void init() { P.init(); D.init(); I0.init(); };

    groupIndex& init(
       const wbvector<T> &P_, const wbvector<T> &D_
    ){ P=P_; D=D_; setup(); return *this; };

    groupIndex& initX(
       wbvector<PERM_T> &P_, wbvector<T> &D_,
       const PERM_T **p=NULL
    ){
       P_.save2(P); D_.save2(D);
       setup(); if (p) (*p)=P.data;
       return *this;
    };

    groupIndex& init(const wbvector<T> &DI, INDEX_T ng=-1) {
       if (DI.isEmpty()) { init(); return *this; }

       if (sINDEX_T(ng)<0) {
          INDEX_T n=0,l=0;
          for (; l<DI.len; ++n, l+=(1+DI[l]));
          if (l!=DI.len) wblog(FL,
             "ERR %s() invalid group index (%d/%d)",FCT,l,DI.len);
          ng=n;
       }
       P.init(DI.len-ng); D.init(ng);

       INDEX_T i,d,k=0,l=0; T *dp=P.data;
       while (l<DI.len) { d=DI.data[l++];
          for (i=0; i<d; ++i, ++k, ++l) { dp[k]=DI.data[l]; }
       }
       if (k!=P.len || l!=DI.len) wblog(FL,
          "ERR %s() %d/%d %d/%d",FCT,k,P.len,l,DI.len);
       setup(); return *this;
    };

    groupIndex& init(INDEX_T ng) {
       D.init(ng); I0.init(ng+1); P.init();
       return *this;
    };

    groupIndex& init_g(INDEX_T ig, INDEX_T nz) {
       if (ig>=D.len || ig+1>=I0.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d,%d)",FCT,ig,D.len,I0.len-1);
       if (ig==0) I0[ig]=0;
       I0[ig+1]=nz; D[ig]=nz-I0[ig];
       return *this;
    };

    groupIndex& setup() {
       INDEX_T n=D.cumsum0(I0,'N');
       if (n!=P.len) wblog(FL,
          "ERR %s() severe size mismatch (%d/%d)",FCT,n,P.len);
       return *this;
    };

    INDEX_T dim(INDEX_T ig) {
       if (ig>=D.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,ig,D.len);
       return D.data[ig];
    };

    INDEX_T numGroups() const { return D.len; };

    void getGroup(INDEX_T ig, T* &p0, T &d) const {
       if (ig>=D.len || I0.len!=D.len+1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       p0=P.data+I0.data[ig];
       d=D.data[ig];
    };
    INDEX_T getGroup(INDEX_T ig, T &j, T &j2) const {
       if (sINDEX_T(ig)==-1) return 0;
       if (ig>=D.len || I0.len!=D.len+1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       j =I0.data[ig  ];
       j2=I0.data[ig+1]; return ig+1;
    };

    template <class T2>
    void getRecs0(const wbMatrix<T2> &, wbMatrix<T2> &, INDEX_T m) const;

    INDEX_T groupIter(INDEX_T &ig, T &j, T &j2) const {
       if (sINDEX_T(ig)>=0) { ++ig; } else { ig=0;
          if (I0.len!=D.len+1) wblog(FL,"ERR %s() "
         "index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       }
       if (ig<D.len) {
          j =I0.data[ig  ];
          j2=I0.data[ig+1];
       }
       else {
          ig=-1;
       }

       return ig+1;
    };

    mxArray* toMx() const {
       const char *fields[]={"P","D","I0","pp"};
       mxArray *S=mxCreateStructMatrix(1,1,4,fields);

       mxSetFieldByNumber(S,0,0,(P+1).toMx());
       mxSetFieldByNumber(S,0,1, D.toMx());
       mxSetFieldByNumber(S,0,2,(I0+1).toMx());

       wbvector<INDEX_T> p;
       mxArray *c=mxCreateCellMatrix(1,D.len);

       for (INDEX_T i=0; i<D.len; ++i) {
          p.init(D[i],P.data+I0[i]);
          mxSetCell(c,i,p.toMx());
       }
       mxSetFieldByNumber(S,0,3,c);

       return S;
    };

  protected:
  private:

     wbvector<T> D;

     wbvector<T> I0;

     wbvector<T> P;
};


#define CTR_IS_NUM(c)      ((c)>='1' && (c)<='9')

#define CTR_IS_SEP_ANY(c)  ((c)==',' || (c)==' ' || (c)==';')
#define CTR_IS_CONJ(c)     ((c)=='*')
#define CTR_IS_OTHER(c)    (CTR_IS_SEP_ANY(c) || CTR_IS_CONJ(c))

#define CTR_IS_VALID(c)    (CTR_IS_NUM(c) || CTR_IS_OTHER(c))

class ctrIdx: public wbvector<unsigned> {

 public:

    ctrIdx(unsigned l=0, unsigned *d=NULL, char c=0)
     : wbvector<unsigned>(l,d), conj(c) {};

    ctrIdx(const char *F, int L, const char *s, char c=0)
     : conj(c) { init(F,L,s); };

    ctrIdx(const char *F, int L, const mxArray *a)
     : conj(0) { init(F,L,a); };

    ctrIdx(const ctrIdx &I)
     : wbvector<unsigned>(I), conj(I.conj) {};

    ctrIdx(const wbvector<unsigned> &I, char c=0)
     : wbvector<unsigned>(I), conj(c) {};

    ctrIdx& init(unsigned l=0, unsigned *d=NULL, char c=0) {
       wbvector<unsigned>::init(l,d); conj=c;
       return *this;
    };

    ctrIdx& init(const wbvector<unsigned> &I, char c=0) {
       wbvector<unsigned>::init(I); conj=c;
       return *this;
    };

    int init(const char *F, int L, const char *s);

    ctrIdx& init(const char *F, int L, const mxArray *a);

    ctrIdx& Index(unsigned l, char c=0) {
       wbvector<unsigned>::init(l); conj=c;
       for (unsigned i=0; i<l; ++i) data[i]=i;
       return *this;
    };

    ctrIdx& Conj() {
       conj=(conj ? 0 : 1);
       return *this;
    };

    mxArray* toMx() const { return toStr().toMx(); };

    wbstring toStr() const;


    char conj;

 protected:
 private:

};


bool isCtrIdx(const mxArray *a) {

   unsigned i=0, rmax=30;

   if (mxIsDblVector(0,0,a)) { wbvector<double> x(FL,a);
      if (mxGetM(a)!=1 || x.len>rmax) return 0;

      for (; i<x.len; ++i) {
         if (x[i]<0 || double(int(x[i]))!=x[i] || x[i]>rmax) return 0;
      }
      if (!x.wbvector<double>::isUnique()) return 0;
      return 1;
   }
   else if (mxIsChar(a)) { ctrIdx idx; wbstring s(a);
      if (idx.init(0,0,s.data)<=0) return 0;
      for (; i<idx.len; ++i) { if (idx[i]>rmax) return 0; }
      if (!idx.wbvector<unsigned>::isUnique()) return 0;
      return 1;
   }
   return 0;
};


#define IDT unsigned long
#define ITAG_LEN sizeof(IDT)

#define CC_ITAG '*'

#define IS_ITAGS_SEP(c) ((c)==',' || (c)==';' || (c)=='|')
#define IS_CHAR_ITAGS(c) ((c)>32 && (c)<127)
#define IS_CHAR_ITAG(c) (((c)>32 && (c)<127) && !((c)==CC_ITAG || IS_ITAGS_SEP(c)))

#define IT2STR(a) (a).itags.toStr().data
#define IT2STR__      itags.toStr().data


class iTag {

 public:

   iTag() : t(0) {};
   iTag(const char *F, int L, const char *s) { init(F,L,s); };
   iTag(const char *s) { init(FL,s); };
   iTag(const iTag &x) : t(x.t) {};

   iTag& init(const char* F, int L, const char *s, unsigned n=-99);

   iTag& init(const char *s=NULL) { return init(FL,s); };

   iTag& init_qdir(char c) { t=0;
      if (c=='-' || CC_ITAG) Conj();
      else if (c!='+') wblog(FL,"ERR %s() invalid c=%c<%d>",FCT,c,c); 
      return *this;
   };

   template<class TQ, class TD>
   iTag& init(const char *F, int L,
      const wbvector< QSpace<TQ,TD> > &A, wbindex ia);

#ifndef NOMEX
   iTag& init(const char* F, int L, const mxArray *a);
#endif

   iTag& operator=(unsigned n) { t=n; return *this; }
   iTag& operator=(const char *s) { return init(FL,s); }

   bool operator==(const iTag &x) const { return t==x.t; };
   bool operator!=(const iTag &x) const { return t!=x.t; };
   bool operator! () const { return !t; };

   bool sameAs(const iTag &x, char *zflag=NULL) const;

   bool isConj(const iTag &x) const {
      IDT c128=128; return (t^c128)==x.t;
   };

   int isValid() const;

   iTag& Conj() {
      ((char&)t) ^= char(128);
      return *this;
   };

   bool isConj() const { 
      return (((char&)t)<0);
   };

   iTag& deConj() {
      char &c=((char&)t); if (c<0) c+=char(128);
      return *this;
   };

   iTag conj()   { iTag x(*this); return x.Conj();   }
   iTag deconj() { iTag x(*this); return x.deConj(); }

   iTag& SetFlag(unsigned k) {
      if (k>=ITAG_LEN) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,ITAG_LEN);
      ((char*)(&t))[k] ^= char(128);
      return *this;
   };

   iTag& UnsetFlag(unsigned k) {
      if (k>=ITAG_LEN) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,ITAG_LEN);
      ((char*)(&t))[k] &= ~char(128);
      return *this;
   };

   iTag& UnsetFlags() {
      unsigned k=1, n=ITAG_LEN;
      char *q = (char*)(&t), c128=char(128);
      for (; k<n; ++k) { if (q[k]<0) {
         if (k>1) wblog(FL,
            "WRN %s() got iTag flag at %d/%d",FCT,k,n);
         q[k]+=c128;
      }}
      return *this;
   };

   unsigned to_str(char *s, unsigned len) const;

   wbstring toStr() const {
      wbstring q(ITAG_LEN+4); to_str(q.data,q.len);
      return q;
   };

   mxArray* toMx() const { return toStr().toMx(); };

   IDT t;

 protected:
 private:

};


class IndexLabels : public wbvector<iTag> {

 public:

   IndexLabels(unsigned n=0) : wbvector<iTag>(n) {};
   IndexLabels(const char *F, int L, const char *s) { init(F,L,s); };
   IndexLabels(const char *s) { init(FL,s); };

   unsigned init(const char* F, int L, const char *s);
   unsigned init(const char *s=NULL) { return init(FL,s); };

#ifdef LOAD_CGC_QSPACE
   unsigned init_qdir(const char *s);
#endif

#ifndef NOMEX
   unsigned init(const char* F, int L, const mxArray *a);
#endif

   IndexLabels& init(unsigned n) {
      wbvector<iTag>::init(n);
      return *this;
   };

   IndexLabels& init3() {
      wbvector<iTag>::init(3); data[2].Conj();
      return *this;
   };

   unsigned Init(const char* F, int L, const char *s, unsigned);

   void SetConj(const char* F, int L, unsigned k) {
      if (!len) return;
      if (k>=len) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,len); 
      data[k].Conj();
   };

   void SetFlag(const char* F, int L, unsigned k, unsigned l) {
      if (!len) return;
      if (k>=len) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,len); 
      data[k].SetFlag(l);
   };

   void UnsetFlag(const char* F, int L, unsigned k, unsigned l) {
      if (!len) return;
      if (k>=len) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,len); 
      data[k].UnsetFlag(l);
   };

   void UnsetFlags() { 
      for (unsigned k=0; k<len; ++k) data[k].UnsetFlags();
   };

   unsigned match(
      const IndexLabels &b, ctrIdx &ia, ctrIdx &ib,
      const char tflag=0) const;

   bool sameAs(unsigned i, const iTag &q, char *zflag=NULL) const {
      return data[i].sameAs(q,zflag);
   };

   bool sameConj(const IndexLabels &b) {
      if (len!=b.len) return 0;
      for (unsigned i=0; i<len; ++i) {
         if (data[i].isConj() ^ b.data[i].isConj()) return 0;
      }
      return 1;
   };

   bool isOp(unsigned r=-1) {
       if (int(r)>=0) {
          if (r<2 || r>3) wblog(FL,"ERR %s() invalid r=%d !?",FCT,r);
          if (len!=r) return 0;
       }
       else if (len<2 || len>3) return 0;
       if (!data[0].isConj(data[1]) || !data[1].isConj()) return 0;
       return 1;
   };

   wbstring toStr(char vflag=0) const;

   mxArray* toMx() const {
      if (len) {
         mxArray* a=mxCreateCellMatrix(1,len);
         for (unsigned i=0; i<len; ++i) {
            mxSetCell(a,i,data[i].toMx());
         }
         return a;
      }
      return mxCreateCellMatrix(0,0);
   };

   int got_op_labels(const char *F=NULL, int L=0) const {
      int q=0, x=0; unsigned i=0, l=len/2;
      for (; i<l; ++i) { 
         if ((x=data[i].isValid())) {
            if (q>0) { if (q!=x) q=-1; } else
            if (!q) { q=x; }
         }
         else wblog(F_L,
            "ERR %s() got invalid itags (%s)",FCT,toStr().data);
         if (data[i].isConj() || data[i].conj()!=data[i+l]) return 0;
      }
      return q;
   };


 protected:
 private:

   unsigned got_text_labels() const {
      if (len) {
         for (unsigned k=0; k<len; ++k) {
            if (!data[k].isValid()) return 0;
         }
      }; return len;
   };

   unsigned checkUnique(const char *F=NULL, int L=0) const {
      if (len) {
         unsigned i,j;
         for (i=0; i<len; ++i)
         for (j=i+1; j<len; ++j) if (data[i]==data[j]) {
            if (F) wblog(F,L,
               "WRN got non-unique itags! (%s)",toStr().data);
            return 0;
         }
      }
      return len;
   };

};


#endif

