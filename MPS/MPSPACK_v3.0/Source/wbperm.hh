#ifndef __WB_PERM_HH__
#define __WB_PERM_HH__

// ----------------------------------------------------------------- //
// permutation class
// Wb,Apr19,06
// ----------------------------------------------------------------- //

class wbperm : public wbvector<PERM_T> {

  public:

    wbperm(PERM_T n=0, char reverse=0) { Index(n,reverse); };

    wbperm(const wbperm &P, char iflag=0)
     : WBPERM() { init(P,iflag); };

    wbperm(const wbperm &p1, const wbperm &p2) 
     : WBPERM() { init(p1,p2); };

    wbperm(const wbperm &P, size_t r) { init_pad(P,r); };

    template<class T>
    wbperm(const wbvector<T> &b, PERM_T offset=0) {
       init(FL,b,offset);
    };

    wbperm(const wbindex &P);

    wbperm(const char *s, PERM_T offset=1) {
       initStr(s,offset);
    }


    PERM_T el1(PERM_T i) const {
       if (len) {
          if (i>=len) wblog(FL,
             "ERR element index out of bounds (%d/%d)",i+1/len);
          return data[i];
       }
       else return i;
    };

    wbperm& operator= (const char *s) {
       return initStr(s,0);
    };

    wbperm& initStr(const char *s, PERM_T offset=1) {
       Str2Idx(s,*this,offset); isValidPerm(FL);
       return *this;
    };

    wbperm& init(const wbperm &P, char iflag=0) {
       if (&P!=this) {
          if (iflag) { P.getIPerm(*this); }
          else NEW(P.len,P.data);
       }
       else if (iflag) Invert();
       return *this;
    };

    wbperm& init(const wbperm &p1_, const wbperm &p2_);

    wbperm& init_pad(const wbperm &P, PERM_T r) {
       if (r<P.len) wblog(FL,
          "ERR invalid %s contruction (%d/%d)",FCT,r,P.len);
       NEW(r); memcpy(data,P.data,P.len*sizeof(P.data[0]));
       for (PERM_T i=P.len; i<r; i++) data[i]=i;
       if (!isValidPerm()) dispInvalidPerm(FL);
       return *this;
    };


    wbperm& init(PERM_T n=0, char reverse=0) {
       return Index(n,reverse);
    };

    wbperm& init(PERM_T n, const PERM_T* d) {
       RENEW(n,d); isValidPerm(FL);
       return *this;
    };

    wbperm& init(const char *F, int L,
       const mxArray* a, PERM_T offset=0
    ){
       WBPERM::init(F,L,a);
       if (offset) operator-=(offset); isValidPerm(F_L);
       return *this;
    }

    template<class T>
    wbperm& init(
       const char *F, int L, const wbvector<T> &b, PERM_T offset=0
    ){
       initT(b.len,b.data);
       if (offset) { (*this)-=offset; }; isValidPerm(F_L);
       return *this;
    };

    wbperm& initTranspose(PERM_T r) {
       WBPERM::init(r); PERM_T r2=r/2;
       if (r%2) wblog(FL,"ERR %s() even input rank required (%d)",FCT,r);
       for (PERM_T i=0; i<r2; ++i) { data[i]=r2+i; data[r2+i]=i; }
       return *this;
    };

    wbperm& initPerm2(PERM_T r) {
       WBPERM::init(r); PERM_T r2=r/2;
       if (r%2) wblog(FL,"ERR %s() even input rank required (%d)",FCT,r);
       for (PERM_T l=0,i=0; i<r2; ++i) {
          data[l++]=i;
          data[l++]=r2+i;
       }
       return *this;
    };

    wbperm& initFirstTo (PERM_T k, PERM_T N) {
       Index(N); mvFirstTo(k);
       return *this;
    };

    wbperm& initLastTo (PERM_T k, PERM_T N) {
       Index(N); mvLastTo(k);
       return *this;
    };

    wbperm& init2Front(const WBPERM &I, PERM_T N);
    wbperm& init2End  (const WBPERM &I, PERM_T N);

    wbperm& init2FrontB(PERM_T m, PERM_T N);
    wbperm& init2EndB  (PERM_T m, PERM_T N);

    wbperm& init2Front(PERM_T k, PERM_T N);
    wbperm& init2End  (PERM_T k, PERM_T N);

    wbperm& Cycle(PERM_T k1, PERM_T k2);
    wbperm& initCycle(PERM_T n, PERM_T k1, PERM_T k2);

    wbperm& Index(PERM_T n, char reverse=0) {
       RENEW(n); if (n==0) return *this;
       if (reverse)
            for (PERM_T l=n-1, i=0; i<n; ++i) data[i]=l-i;
       else for (PERM_T i=0; i<n; ++i) data[i]=i;
       return *this;
    };

    wbperm& Index(PERM_T i1, PERM_T i2) {
       if (i1>i2) { RENEW(0); return *this; }

       RENEW(i2-i1+1);
       for (PERM_T i=0; i<len; i++) data[i]=i+i1;

       return *this;
    };

    char isValidPerm(const char *F, int L, PERM_T r=-1) const;
    char isValidPerm(PERM_T r=-1) const { return isValidPerm(0,0,r); };

    wbperm& Invert() {
       wbperm P0(*this);
       return P0.getIPerm(*this);
    };

    wbperm& getIPerm(wbperm &iP) const {
       if (&iP==this) {
          wbperm X; getIPerm(X); X.save2(iP);
          return iP;
       }

     #ifdef __WBDEBUG__
       isValidPerm(FL);
     #endif

       iP.RENEW(len);

       for (PERM_T *p=iP.data, i=0; i<len; ++i) {
          if (data[i]>=len) wblog(FL,"ERR %s() "
             "permutation out of bounds (%d/%d)",FCT,data[i],len);
          p[data[i]]=i;
       }
       return iP;
    };

    wbperm& Rotate(sPERM_T l);
    wbperm& initRotate(PERM_T n, sPERM_T k);

    wbperm& Permute(const wbperm &P, char iflag, PERM_T r);
    wbperm& Permute(const wbperm &P, char iflag=0){
       return Permute(P, iflag, len>P.len ? len : P.len);
    };

    bool sameAs(const wbperm &b) const {
       if (!len) {
          if (b.len) return b.isIdentityPerm();
          else return 1;
       }
       if (!b.len) return isIdentityPerm();
       else return (*this)==b;
    };

    bool isIdentityPerm() const {
        bool q=1; PERM_T i=0;
        for (; i<len; ++i) if (data[i]!=i) { q=0; break; }
        for (; i<len; ++i) if (data[i]>=len) wblog(FL,"ERR %s() "
            "permutation out of bounds (%d: %d/%d)",FCT,i,data[i],len);
        return q;
    };

    bool isIdentityPerm(const char *F, int L, PERM_T r) const {
        if (len!=r) { if (F) wblog(F,L,
           "ERR %s() invalid permutation (len=%d/%d)",FCT,len,r);
           return 0;
        }
        bool q=1; char mark[len]; memset(mark,0,len);
        PERM_T i=0;

        for (; i<len; ++i) {
           if (data[i]==i) { ++mark[i]; } else { q=0; break; }
        }
        for (; i<len; ++i) {
           if (data[i]>=r) wblog(F_L, "ERR %s() invalid permutation\n"
              "(index out of range %d/%d)",FCT,data[i],r); 
           if ((++mark[data[i]])>1) wblog(F_L,"ERR %s() "
              "invalid permutation (non-unique index)",FCT); 
        }
        return q;
    };

    bool isReversePerm() const {
        if (len==0) return 0;
        for (PERM_T m=len-1, i=0; i<len; i++) if (data[i]!=m-i) return 0;
        return 1;
    };

    bool isTranspose(PERM_T m, PERM_T n=-1) const {

        if (sPERM_T(n)<0) { n = (len>m ? (len-m) : 0); }
        if (!len) {
           if (m==0 && n==0) return 1;
           else return 0;
        }
        if (len!=m+n) wblog(FL,
           "ERR %s() size mismatch (%d == %d + %d)",FCT,len,m,n);

        PERM_T i=0;
        for (; i<m;   i++) if (data[i]!=i+n) return 0;
        for (; i<len; i++) if (data[i]!=i-m) return 0;
        return 1;
    };

    wbperm& times(const wbperm &p2, wbperm &Pout) const {
       if (!isValidPerm()) wblog(FL,
          "ERR %s() invalid permutation (%s)",__FUNCTION__,toStr().data);
       if (!p2.isValidPerm()) wblog(FL,
          "ERR %s() invalid permutation (%s)",__FUNCTION__,p2.toStr().data);
       if (len!=p2.len) wblog(FL,
          "ERR %s() incompatible permutations (%d/%d)",
         __FUNCTION__,len,p2.len);

       wbperm P(len);
       for (PERM_T i=0; i<len; i++) P.data[i]=data[p2.data[i]];

       P.save2(Pout); return Pout;
    };

    wbperm& flipIdx(wbperm &P) const;
    wbperm& FlipIdx();
    wbperm  flipIdx() const {
       wbperm P; flipIdx(P);
       return P;
    };

    void ToRaw() { ToRaw(len); };
    void ToRaw(PERM_T n) {
        if (n<len) wblog(FL,"WRN ToRaw() with n=%d/%d !??",n,len);
        if (n==0) return; else n--;

        Flip();

        for (PERM_T i=0; i<len; i++) {
            if (i>n) wblog(FL,"ERR wbperm - invalid n=%d (%d)",n,i);
            data[i]=n-data[i];
        }
    }

    mxArray* toMx(const char tflag=0) const {
       return toMx_offset(1,tflag);
    };


  protected:
  private:

    void dispInvalidPerm(const char* F, int L) const {
    wblog(F,L,"WRN invalid permutation [%s]",toStr().data); };

};


wbperm& wbperm::flipIdx(wbperm &P) const {

    if (&P==this) wblog(FL,
       "ERR flipIdx(P) calls itself\nused flipIdx() instead!");

    P.RENEW(len); if (!len) return P;

    for (PERM_T m=len-1, i=0; i<len; i++)
    P.data[i]=m-data[i];

    return P;
};

inline wbperm& wbperm::FlipIdx() {
    if (!len) return *this;

    for (PERM_T m=len-1, i=0; i<len; i++) {
       if (data[i]>m) wblog(FL,"ERR index out of bounds (%d/%d)",data[i],m);
       data[i]=m-data[i];
    }

    return *this;
};


#endif

