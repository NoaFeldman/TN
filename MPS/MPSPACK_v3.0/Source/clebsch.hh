#ifndef __WB_CLEBSCH_GORDAN_HH__
#define __WB_CLEBSCH_GORDAN_HH__

/*====================================================================//
// NOTES FOR LATER REGARDING 6J ETC.

## CGR_TRACE_FOR_ID
   given the normalization |CGC|^2=1 up to OM
   isIdentityCGC() for rank-4,6,... tensors requires extra information
   other than QSet and size info => store trace(CGC) as info.cgr.cgt

   Comments
   The following contraction results in projectors which, by construction,
   are orthogonal (due to orthogonality of input standard CGCs)!
         contract( (1,1|2), 3, (1,1|2), 3 )
     vs. contract( (1,1|0), 3, (1,1|0), 3 )
   and of different (Schmidt) rank 1 and 3, respectively!

   applying 1j symbol to invert out-to-in say, changes the
   interpretation to combining 3 input spins (S=1) into S=1 out.
   Again 1+1=0+3 for the intermediate state space
   which subsequently can be combined with another S=1
   into S=1

    => it appears, this operation results in the *same* rank-4 CGC's
       (having OM=2)
    => applying a 1j symbol (\propto unitary), this does not change
       the underlying eigenvalues / Schmidt spectrum, hence the
       eigenvalues of the original projector *remain intact*!
    => trace(CGC) is *preserved* when applying a 1j symbol
       (\propto unitary) to revert arrows! (in this sense,
       trace(CGC) may be useful also in other contexts later!


*/







#ifndef RTD
  #ifdef __WB_MPFR_HH__


    #define MTI unsigned long
    #define RTD Wb::quad

  #else

    #define MTI unsigned long
    #define RTD double

  #endif
#endif


#ifndef gTQ


   #define gTQ int
#endif

#ifndef gTD
#define gTD double
#endif

#define CDATA_TQ    CData<TQ,RTD>
#define cdata__     cdata<RTD>
#define cgsparray   wbsparray<RTD>

   template <class TQ, class TD> class CData;
   template <class TQ> class CRef;

   class QVec;

   int CG_VERBOSE=0;



   #define CG_SKIP_REPS 1E-17



#ifdef __WB_MPFR_HH__
   #define CG_SKIP_EPS1 1E-14
   #define CG_SKIP_EPS2 1E-20


   #define CG_EPS1 1E-10
   #define CG_EPS2 1E-12

#else

   #define CG_SKIP_EPS1 1E-14
   #define CG_SKIP_EPS2 1E-16



   #define CG_EPS1 1E-10
   #define CG_EPS2 1E-12

#endif


enum QTYPE_SET {
   QTYPE_UNKNOWN,
      QTYPE_P,
      QTYPE_ZN,
      QTYPE_U1,
   QTYPE_ASEP,
      QTYPE_SUN,
      QTYPE_SpN,
   QTYPE_NUM_TYPES
};



const char* QTYPE_STR[QTYPE_NUM_TYPES]=  {
   "???",
   "P",
   "ZN",
   "A",
   "!?!",
   "SUN",
   "SpN"
};

#define CGC_ALL_ABELIAN "A*"

enum CGD_TYPE { CGD_DEFAULT,
   CGD_ABELIAN,
   CGD_IDENTITY,
   CGD_REF_INIT,
   CGD_BSZ_INIT,
   CGD_FROM_CGR,
   CGD_FROM_CTR,
   CGD_FROM_DEC,
   CGD_NEW3,
   CGD_J1SY_CGC,
CGD_NUM_TYPES };

const char* CGD_TYPE_STR[CGD_NUM_TYPES]=  { "",
  "d:iA",
  "d:iE",
  "d:iR"," (bare info only)
  "d:iS",
  "d:cR",
  "d:cC",
  "d:cD",
  "d:c3",
  "d:J1r"
};


enum CGR_TYPE { CGR_DEFAULT,
   CGR_ABELIAN,
   CGR_CTR_SCALAR,
   CGR_CTR_ZERO,
CGR_NUM_TYPES };


const char* CGR_TYPE_STR[CGR_NUM_TYPES]=  { "",
  "r:iA",
  "r:cS",
  "r:c0",
};

#ifdef WB_SPARSE_CLOCK
WbClock wbc_sparse_gss("cgs::getSymStates");
#endif


template <class T> class qset;
template <class T> class CStore;

template <class T> class QMap;
template <class TQ, class TD> class blockSpaceQS;

#include <map>

namespace CG {

   template <class TD> inline
   int signFirstVal(const char *F, int L,
      const TD *d, SPIDX_T n,
      double eps1=CG_EPS1, double eps2=CG_EPS2
   );

   template <class TD> inline
   int rangeSignConvention(const char *F, int L,
      TD *d, SPIDX_T n,
      double eps1=CG_EPS1, double eps2=CG_EPS2
   );

   template<class T>
   double FixRational(const char *F, int L,
      T* d, SPIDX_T n, unsigned niter=0,
      T eps1=-1, T eps2=-1
   );
};


class QType {

  public:

    QType(QTYPE_SET t=QTYPE_UNKNOWN, unsigned s=0)
     : type(t), sub(s) { if (t || s) validType(FL); };

    QType(const char *s ) : type(QTYPE_UNKNOWN), sub(0) { init(0,0,s); };
    QType(const QType &q) : type(q.type), sub(q.sub) { };

    QType(const char *F, int L, const mxArray *a)
      : type(QTYPE_UNKNOWN), sub(0) { init(F,L,a); };

   ~QType() { init(); };

    QType& init() { type=QTYPE_UNKNOWN; sub=0; return *this; };

    QType& initx(QTYPE_SET t, unsigned m=0) {
       if (!t && !m) { init(); return *this; }
       type=t; sub=m; validType(FL);
       return *this;
    };

    QType& init(QTYPE_SET t, unsigned m=0) {

       if (!t && !m) { init(); return *this; }

       type=t;
       if (type==QTYPE_SpN) { sub=m/2;
          if (m%2) wblog(FL,"ERR %s() invalid sub=%d with SpN",FCT,m);
       }
       else { sub=m; }

       validType(FL); return *this;
    };

    int init_s(const char *s);

    QType& init(const char *s) { return init(FL,s); }

    QType& init(const char *F, int L, const char *s) {
       if (!s || !s[0]) { type=QTYPE_UNKNOWN; sub=0; }
       else {
          if (init_s(s)) wblog(FL,"ERR invalid QType '%s'",s);
          else validType(F_L);
       }
       return *this;
    };

    QType& init(const char *F, int L, const mxArray *a) {
       wbstring s(F_L,a);
       return init(F_L,s.data);
    };

    QType& init(const QType &q) {
       type=q.type; sub=q.sub; return *this;
    };

    bool isUnknown() const {
       return (type==QTYPE_UNKNOWN && sub==0);
    };

    bool isKnown() const {
       validType();
       return (type!=QTYPE_UNKNOWN);
    };

    bool validType(const char *F=0, int L=0) const;

    bool isAbelian(const char *F=NULL, int L=0) const {
       if (type<=QTYPE_UNKNOWN || type>=QTYPE_NUM_TYPES) wblog(F_L,
          "ERR invalid type %d, `%s'",type,toStr().data); 
       return ((type && type<QTYPE_ASEP)? 1 : 0);
    };

    bool isNonAbelian(const char *F=NULL, int L=0, char lflag=0) const {
       if (!lflag && !isUnknown()) {
       if (type<=QTYPE_UNKNOWN || type>=QTYPE_NUM_TYPES) wblog(F_L,
          "ERR invalid type %d, `%s'",type,toStr().data); }
       return ((type && type>QTYPE_ASEP)? 1 : 0);
    };

    bool isSU2() const { return (type==QTYPE_SUN && sub==1); };
    bool isSUN() const { return (type==QTYPE_SUN); };
    bool isSU(unsigned N) const { return (type==QTYPE_SUN && sub+1==N); };

    bool isU1()  const { return (type==QTYPE_U1  && sub==0); };

    int hasMP() const {
       return (isAbelian() ? 0 : (sub<=1 ? 1 : 2));
    };

    int hasMP(unsigned r) const {
       return (isAbelian() || r<3 || (sub<2 && r<4) ? 0 : 1);
    };

    unsigned qlen() const {
       if (type && type<QTYPE_ASEP) { return 1; }
       if (type>=QTYPE_NUM_TYPES) wblog(FL,
          "ERR invalid type `%s' (%d)",toStr().data,type);
       return sub;
    };

    unsigned qrank() const {
       if (type && type<QTYPE_ASEP) { return 0; }
       if (type>=QTYPE_NUM_TYPES) wblog(FL,
          "ERR invalid type `%s' (%d)",toStr().data,type);
       return sub;
    };

    bool isLargeD(unsigned d) const {
       if (type && type<QTYPE_ASEP) { return (d>1); }
       switch (type) {
          case QTYPE_SUN: return (d>(sub!=1 ? exp10(qlen()) : 20));
          case QTYPE_SpN: return (d>exp10(qlen())); 

          case QTYPE_UNKNOWN : return (d>100);
          default: wblog(FL,
         "ERR invalid type `%s' (%d)",toStr().data,type); 
       }
       return 1;
    };

    bool isLargeD(const wbvector<unsigned> &dd) const {
       for (unsigned i=0; i<dd.len; ++i) {
          if (isLargeD(dd[i])) return 1;
       }
       return 0;
    };

    template <class T>
    unsigned qdim(const T* q) const;

    template <class T>
    wbvector<unsigned>& QDim(
       const T* q, unsigned r, unsigned stride,
       wbvector<unsigned> &S) const;

    template<class TQ>
    unsigned getDual(const TQ* q0, TQ* q) const {
       switch (type) {
          case QTYPE_P  :
          case QTYPE_U1 : { q[0]=-q0[0]; return 1; }
          case QTYPE_ZN : { q[0]=(sub-q0[0])%sub; return 1; }

          case QTYPE_SUN: {
             unsigned i=0, n=qlen(), l=n-1;
             for (; i<n; ++i) { q[l-i]=q0[i]; }; return n;
          }
          case QTYPE_SpN: {
             unsigned i=0, n=qlen();
             for (; i<n; ++i) { q[i]=q0[i]; }; return n;
          }
          default:
          wblog(FL,"ERR dual not yet defined for '%s'",toStr().data); 
       }

       return 0;
    };

    wbstring toStr(const char tflag=0) const { wbstring s;
       char e=0;

       if (type==QTYPE_U1  ) { s=(tflag ? QTYPE_STR[type] : "U(1)"  ); } else
       if (type==QTYPE_P   ) { s=(tflag ? QTYPE_STR[type] : "Parity"); } else
       if (!type && !sub) { s=(tflag ? "???" : "unknown");}

       if (s.data) {
          if (sub) wblog(FL,
             "ERR invalid QType %s with sub=%d",s.data,sub);
          return s;
       }

       if (sub<1 || sub>99) { e=1; }
       else if (type==QTYPE_ZN) { if (sub<2) e=2;
          if (tflag) 
               { s.init(3); sprintf(s.data,"Z%d",  sub); }
          else { s.init(5); sprintf(s.data,"Z(%d)",sub); }
       }
       else if (type==QTYPE_SUN) {
          if (tflag)
               { s.init(5); sprintf(s.data,"SU%d",  sub+1); }
          else { s.init(7); sprintf(s.data,"SU(%d)",sub+1); }
       }
       else if (type==QTYPE_SpN) {
          if (tflag) 
               { s.init(5); sprintf(s.data,"Sp%d",  2*sub); }
          else { s.init(7); sprintf(s.data,"Sp(%d)",2*sub); }
       }
       else {
          if (type) wblog(FL,"ERR %s() invalid type=%d",FCT,type);
          s=(tflag ? "???" : "unknown");
       }

       if (e) wblog(FL,
          "ERR invalid QType %s with sub=%d", type < QTYPE_NUM_TYPES ?
           QTYPE_STR[type] : "(type out of bounds)", sub
       );

       return s;
    };

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    mxArray* toMx(const char tflag=0) const {
       const size_t n=32; wbstring s(n);
       size_t l=snprintf(s.data,n,"%s",toStr(tflag).data);
       if (l>=n) wblog(FL,
          "WRN %s() QType too long (%s; %d) !??",FCT,s.data,l);
       return s.toMx();
    };

    QType& operator=(const char *s) { return init(FL,s); };
    QType& operator=(const QType &q) { return init(q); };

    QType& operator=(QTYPE_SET t) {
       type=t; sub=0; if (t) validType(FL);
       return *this;
    };

    bool operator<(const QType &q) const {
       if (type!=q.type) return (type<q.type);
       else return (sub<q.sub);
    };

    bool operator> (const QType &q) const { return !((*this)<=q); };
    bool operator<=(const QType &q) const {
       if (type!=q.type) return (type<q.type);
       else return (sub<=q.sub);
    };

    bool operator== (const QType &q) const {
        return (type==q.type && sub==q.sub);
    };
    bool operator!= (const QType &q) const {
        return (type!=q.type || sub!=q.sub);
    };

    QTYPE_SET type;
    unsigned sub;


  private:

    int atoi(const char *s, int &k, unsigned n) {
       if (s[0]!='(') { return Wb::atoi(s,k,n); }
       else {
          unsigned l=strlen(++s);
          if (l && s[l-1]==')') {
             char x[l]; memcpy(x,s,l); x[l-1]=0;
             return Wb::atoi(x,k,n);
          }
       }
       return -6;
    }
};

  QType SU2("SU2");

bool operator!(const QType &q) {
   if (q.type == QTYPE_UNKNOWN) {
      if (q.sub) wblog(FL,
         "ERR %s() invalid QType %s",FCT,q.toStr().data);
      return 1;
   }
   return 0;
};


   map<QType, time_t> load_cstore;

QType load_store_qtype(const char *F, int L, const Wb::matFile &M);
void load_RCStore(const char *F, int L, const QVec &qvec);


class QVec : public wbvector<QType> {
  public:

    QVec() : wbvector<QType>() {};

    QVec(unsigned l, const QType* d0, char ref=0)
     : wbvector<QType>(l,d0,ref) {};

#ifndef NOMEX
    QVec(const char *F, int L, const mxArray *a, unsigned d=0)
     : wbvector<QType>() { init(F,L,a,d); };

     QVec& init(const char *F, int L, const mxArray *a, unsigned d=0){
        if (!a) { init(); return *this; }
        if (mxIsChar(a)) { wbstring s(a); return init(F,L,s.data,d); }
        else wblog(FL,"ERR %s() invalid usage",FCT);
        return *this;
     };
#endif


     QVec& init(unsigned l=0)  {
        wbvector<QType>::init(l); return *this;
     };

     QVec& init(const char *F, int L, const char *s, unsigned d=0);

     unsigned Qlen() const;

     unsigned Qrank() const;

     unsigned Qlen(wbvector<unsigned> &dd) const;
     unsigned Qlen(
        wbvector<unsigned> &dd,
        wbvector<unsigned> &dz
     ) const;

     unsigned Qpos(wbvector<unsigned> &dc) const;

     unsigned Qlenz() const {
        unsigned i=0, n=0;
        for (; i<len; ++i) { n+=(data[i].qlen() + data[i].qrank()); }
        return n;
     }

     unsigned largestRank() const {
        unsigned i=0, q=0, r=0;
        for (; i<len; ++i) { q=data[i].qrank(); if (r<q) { r=q; }}
        return r;
     };

     template <class TQ>
     unsigned QDim(const TQ *q) const;

     template <class TQ>
     unsigned QDim(const TQ *qq, unsigned k) const;

     template <class TQ>
     wbvector<unsigned>& QDim(
        const TQ *qq, wbvector<unsigned> &S
     ) const;

     template <class TQ>
     wbvector<unsigned>& QDim(
        const TQ *qq, unsigned k, unsigned r, wbvector<unsigned> &S
     ) const;

     template <class TQ>
     wbMatrix<TQ>& getQsub(
        const wbMatrix<TQ>& Q, const wbindex &I,
        wbMatrix<TQ>& QI,
        wbMatrix<TQ> *Qx=NULL
     ) const;

     bool isAbelian() const { return allAbelian(); }

     char allAbelian() const { 
        for (unsigned i=0; i<len; ++i) {
           if (!data[i].isAbelian()) return 0;
        }
        return 1;
     };

     bool anyAbelian() const {
        if (!len) return 1;
        for (unsigned i=0; i<len; ++i) {
           if (data[i].isAbelian()) return 1; }
        return 0;
     };

     bool noneAbelian() const { return !anyAbelian(); };

     int hasMP() const {
        int q=0, qi; unsigned i=0;
        if (len)
           for (; i<len; ++i) { if ((qi=data[i].hasMP())>q) q=qi; }
        else wblog(FL,"WRN %s() got empty QVec",FCT);

        return q;
     };

     int hasMP(unsigned r) const {
        if (!len) wblog(FL,"WRN %s() got empty QVec",FCT);
        for (unsigned i=0; i<len; ++i) { if (data[i].hasMP(r)) return 1; }
        return 0;
     };


     bool sameType (const QVec &b) const {
        if (len!=b.len) return 0;
        for (unsigned i=0; i<len; i++)
        if (data[i].type!=b.data[i].type) return 0;
        return 1;
     };

     wbstring toStr(const char vflag=0) const;
};


class QDir: public wbvector<char> {

  public:

    QDir() : wbvector<char>() {};

    QDir(unsigned l, const char* d=NULL, char ref=0)
     : wbvector<char>(l,d,ref) {};

    QDir(const IndexLabels& it) { init(it); };
    QDir(const char *F, int L,const char *s) { init(F,L,s); };
    QDir(const char *s) { init(FL,s); };

#ifndef NOMEX
    QDir(const char *F, int L, const mxArray *a) { init(F,L,a); };
#endif


    QDir& init(unsigned l=0, const char* d=NULL, char ref=0)  {
       wbvector<char>::init(l,d,ref);
       return *this;
    };

    QDir& operator=(const char *s) { return init(FL,s); };

    QDir& init(const char *F, int L,const char *s);
    QDir& init(const char *s) { return init(FL,s); };
    QDir& init(const IndexLabels& it);

    QDir& init(const char *F, int L,const mxArray *a) {
       wbstring s(F,L,a); init(F,L,s.data);
       return *this;
    };

    QDir& init_iout(const char *F, int L,
       unsigned r, unsigned iout
    );

    unsigned isSorted() const {
       unsigned i=0, j=0;

     #ifndef NSAFEGUARDS
       for (; i<len; ++i) if (!data[i]) {
          wblog(FL,"ERR %s() invalid qdir '%s'",FCT,toStr().data);
       }; i=0;
     #endif

       for (; i<len; ++i) { if (data[i]>0) break; };
       if (i==len) return 1;
       if (i) return 0;

       for (; i<len; ++i) { if (data[i]<0) break; }
       if (i==len) return (len+1);

       for (j=i+1; i<len; ++i) { if (data[i]>0) break; }
       if (i==len) return j;

       return 0;
    };

    bool operator!=(const char *s) const { return !((*this)==s); };
    bool operator==(const char *s) const;

    bool operator!=(const QDir &q) const { return !((*this)==q); };
    bool operator==(const QDir &q) const {
       return this->wbvector<char>::operator==(q);
    };

    QDir& Conj() {
       for (unsigned i=0; i<len; ++i) { data[i]=-data[i]; }
       return *this;
    };

    char& Conj(const unsigned &i) {
       if (i>=len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,len);
       return (data[i]=-data[i]);
    };

    unsigned nconj() {
       unsigned i=0, n=0;
       for (; i<len; ++i) { if (data[i]<0) ++n; }
       return n;
    };

    unsigned nconj(const ctrIdx &Ix) {
       unsigned n=0, i=0;
       wbvector<char> q(*this); const unsigned* ix=Ix.data;

       for (; i<Ix.len; ++i) { 
          if (ix[i]>=len) wblog(FL,"ERR %s() index out of bounds"
             "(%d/%d: %d/%d)",FCT,i,Ix.len, ix[i],len);
          q[ix[i]]=0;
       }
       if (Ix.conj)
            { for (i=0; i<len; ++i) { if (q[i]>0) ++n; }}
       else { for (i=0; i<len; ++i) { if (q[i]<0) ++n; }}

       return n;
    };

    wbstring toStr() const;
    wbstring toTag() const;

    mxArray* toMx() const { return toStr().toMx(); };

  protected:
  private:

};


class cgdStatus {

  public:

    cgdStatus(CGD_TYPE t_=CGD_DEFAULT)
     : t(t_), ctime(0), mtime(0), ID(0) { init_time(); };

    cgdStatus(double t_) : t(CGD_TYPE(t_)), ctime(0), mtime(0), ID(0) {
       if (double(t)!=t_ || t>=CGD_NUM_TYPES) wblog(FL,
          "ERR %s() invalid CData type %g",FCT,t_
       ); 
       init_time();
    };

    cgdStatus(const char *F, int L, mxArray *a)
     : t(CGD_DEFAULT), ctime(0), mtime(0), ID(0) { init(F,L,a); };

    template <class TQ>
    cgdStatus(const CDATA_TQ& C) { init(C.cstat); }

    cgdStatus& init(const CGD_TYPE &t_) {
       t=t_; init_time(); return *this;
    };

    cgdStatus& update_m() { init_time('m'); return *this; };
    cgdStatus& update_m(const CGD_TYPE &t_) {
       t=t_; init_time('m'); return *this;
    };

    CGR_TYPE init(const char *F, int L, const mxArray *a);

    cgdStatus& operator=(const cgdStatus &b) {
       t=b.t; ctime=b.ctime; mtime=b.mtime; ID=b.ID;
       return *this;
    };

    int cmp(const char *F, int L, const cgdStatus &b) const {
       if (ID!=b.ID || ctime!=b.ctime) wblog(F_L,
          "ERR %s() got CData status mismatch!%N   %s%N<> %s",
          FCT, toStr('V').data, b.toStr('V').data);
       if (mtime<ctime || b.mtime<b.ctime) wblog(F_L,
          "ERR %s() invalid mtime",FCT);
       return (mtime<b.mtime ? -1 : (mtime==b.mtime ? 0 : +1));
    };

    bool olderThan(const cgdStatus &b) const { return (cmp(FL,b)<0); };

    bool hasID(char lflag=0) const {
       if (!lflag)
            return (ctime && mtime && ID);
       else return (ctime || mtime || ID);
    };

    bool operator==(const CGD_TYPE t_) const { return (t==t_); };
    bool operator!=(const CGD_TYPE t_) const { return (t!=t_); };

    bool operator==(const cgdStatus &b) const {
        return (ctime==b.ctime && mtime==b.mtime && ID==b.ID);
    };

    bool operator!=(const cgdStatus &b) const { return !((*this)==b); };

    bool isEmpty() const {
       return (!ctime && !mtime && !ID && t==CGD_DEFAULT);
    };

    bool sameAs(const cgdStatus &b, char lflag=0) const {
       if (ctime!=b.ctime || ID!=b.ID) return 0;

       if (lflag>2) {
          if (lflag=='l') lflag=1; else
          if (lflag=='L') lflag=2; else
          wblog(FL,"ERR %s() invalid lflag=%d",FCT,lflag);
       }
       if (lflag<1 && t!=b.t) return 0;
       if (lflag<2 && mtime!=b.mtime) return 0;
       return 1;
    };

    const char* tstr() const {
       if (t>=CGD_NUM_TYPES) wblog(FL,"ERR %s() "
          "cgd type out of bounds (%d/%d)",FCT,t,CGD_NUM_TYPES);
       return CGD_TYPE_STR[t];
    };

    wbstring toStr(char vflag=0) const;

    mxArray* toMx() const;
    mxArray* toMx(CGR_TYPE rt) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    mxArray* add2MxStruct(mxArray *S, unsigned k) const;

    CGD_TYPE t;

    double ctime;
    double mtime;
    unsigned ID;

  protected:
  private:

    cgdStatus& init_time(char flag=0);
    wbstring to_vstr(double t, char vflag=0) const;

    void setID() {
       ID=unsigned(rand())>>2;
    };
};

bool operator!(const cgdStatus &q) {
   return (q.t==CGD_DEFAULT && !q.ctime && !q.mtime && !q.ID);
};


class cgrType {

  public:

    cgrType(CGR_TYPE t_=CGR_DEFAULT) : t(t_) {};

    cgrType(double t_) { init(FL,t_); };

    template<class T>
    cgrType& init(const char *F, int L, T t_) {
        t=CGR_TYPE(t_);
        if (T(t)!=t_ || t>=CGR_NUM_TYPES) wblog(F_L,
           "ERR %s() invalid CRef type %g",FCT,double(t_)
        ); 
        return *this;
    };

    cgrType& operator=(const cgrType  &q) { t=q.t; return *this; };
    cgrType& operator=(const CGR_TYPE &t_) { t=t_; return *this; };

    bool operator< (const CGR_TYPE t_) const { return (t< t_); };
    bool operator<=(const CGR_TYPE t_) const { return (t<=t_); };
    bool operator> (const CGR_TYPE t_) const { return (t> t_); };
    bool operator>=(const CGR_TYPE t_) const { return (t>=t_); };
    bool operator==(const CGR_TYPE t_) const { return (t==t_); };
    bool operator!=(const CGR_TYPE t_) const { return (t!=t_); };

    bool operator< (const cgrType &q) const { return (q.t< t); };
    bool operator<=(const cgrType &q) const { return (q.t<=t); };
    bool operator> (const cgrType &q) const { return (q.t> t); };
    bool operator>=(const cgrType &q) const { return (q.t>=t); };
    bool operator==(const cgrType &q) const { return (q.t==t); };
    bool operator!=(const cgrType &q) const { return (q.t!=t); };

    const char* tostr() const {
       if (t>=CGR_NUM_TYPES) wblog(FL,"ERR %s() "
          "cgrType out of bounds (%d/%d)",FCT,t,CGR_NUM_TYPES);
       return CGR_TYPE_STR[t];
    };

    wbstring toStr() const { return tostr(); }

    CGR_TYPE t;

  protected:
  private:
};


template <class TQ>
class qset : public wbvector<TQ> {


  public:
    qset() : wbvector<TQ>() {};

    qset(unsigned l, const TQ* d0=NULL, char ref=0)
     : wbvector<TQ>(l,d0,ref) {};

    qset(const qset<TQ> &q)
     : wbvector<TQ>(q.len,q.data) {};

    qset(const QType &Q, const TQ *qq)
     : wbvector<TQ>() { init(Q,qq); };

    qset(const QVec &Q, const TQ *qq)
     : wbvector<TQ>() { init(Q,qq); };

    qset(const char *F, int L, const mxArray *a,
       const QVec *Q=NULL, char check_type=1) : wbvector<TQ>() {
       init(F,L,a,Q,check_type);
    };

    qset(const qset<TQ> &q1, const qset<TQ> &q2)
     : wbvector<TQ>() { init(q1,q2); };

    qset(const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q3)
     : wbvector<TQ>() { init(q1,q2,q3); };

    qset(const qset<TQ> &q1, char op, const qset<TQ> &q2);

    qset& init(const qset<TQ> &q1) {
       wbvector<TQ>::init(q1);
       return *this;
    };

    qset& init(
       const char *F, int L, const mxArray *a,
       const QVec *Q=NULL, char check_type=1
    ){
       wbvector<TQ>::init(F,L,a,NULL,check_type);
       if (Q) {
          unsigned n=Q->Qlen();
          if ((this->len)%n) wblog(FL,
             "ERR %s() got invalid qset.len=%d @ %d",FCT,this->len,n
          );
       }
       return *this;
    };

    qset& init(const QType &Q, const TQ *qq) {
       wbvector<TQ>::init(Q.qlen(),qq);
       return *this;
    };

    qset& init(const QVec &Q, const TQ *qq) {
       init(Q.Qlen(),qq);
       return *this;
    };

    qset& init(const qset<TQ> &q1, const qset<TQ> &q2) {
       init (q1.len, q1.data, q2.len, q2.data);
       return *this;
    };

    qset& init(const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q3) {
       init(q1.len, q1.data, q2.len, q2.data, q3.len, q3.data);
       return *this;
    };

    qset& init() { wbvector<TQ>::init(); return *this; };

    qset& init(unsigned l, const TQ *d0=NULL) {
       wbvector<TQ>::RENEW(l,d0);
       return *this;
    };

    qset& init(
       unsigned l1, const TQ *d1,
       unsigned l2, const TQ *d2
    );

    qset& init(
       unsigned l1, const TQ* d1,
       unsigned l2, const TQ *d2,
       unsigned l3, const TQ *d3
    );

    qset<TQ>& init(const TQ q1) { return init(1U,&q1); }

    qset<TQ>& init(const TQ q1, const TQ q2) {
       return init(1,&q1,1,&q2);
    };
    qset<TQ>& init(const TQ q1, const TQ q2, const TQ q3) {
       return init(1,&q1,1,&q2,1,&q3);
    };

    qset<TQ>& initStride(
       const TQ* d_, unsigned n, unsigned r, unsigned QDIM
    ){
       init(n*r); if (d_ && this->len) {
          Wb::cpyStride(this->data,d_,n,NULL,r,-1,QDIM);
       }
       return *this;
    };

    bool operator<(const qset<TQ> &b) const;

    int checkNQs(const QVec &qq, const char *F=NULL, int L=0);

    bool isConsistent(
    const QVec &qq, int rank=-1, const char *F=NULL, int L=0) const;

    unsigned QDim(const QVec &q, int rank=-1) const;

    qset<TQ>& times(TQ fac, qset<TQ> &x) const {
        x=*this; x*=fac; return x;
    };

};


template <class TQ>
bool qset<TQ>::operator<(const qset<TQ> &b) const {

   if (!this->len || !b.len) wblog(FL,
      "ERR %s() got empty object (%d/%d)",FCT,this->len,b.len);


   return wbvector<TQ>::operator<(b);
};


template <>
bool qset<double>::operator<(const qset<double> &b) const {

   if (!len || !b.len) wblog(FL,
      "ERR %s() got empty object (%d/%d)",FCT,len,b.len);
   if (len!=b.len) wblog(FL,
      "ERR %s() got length mismatch (%d/%d)",FCT,len,b.len);

   double *d=b.data;
   for (unsigned i=0; i<len; i++) {
      if (float(data[i])<float(d[i])) return 1;
      if (float(data[i])>float(d[i])) return 0;
   }
   return 0;
};

template <>
bool qset<float>::operator<(const qset<float> &b) const {
   wblog(FL,"ERR %s() double precision prefered over '%s' (0x%lX)",
   getName(typeid(float)).data,&b);
   return 0;
};


template <class TQ>
class QSet {

  public:

    QSet() : t(QTYPE_UNKNOWN) {};

    QSet(const QSet &B) { t=B.t; qs=B.qs; qdir=B.qdir; };

    QSet(const QType &t_, const qset<gTQ> &qs_,
       unsigned iout=0,
       char isref=0
    ){
       init(t,qs,iout,isref);
    };

    QSet(const QType &t_, const IndexLabels &it, unsigned r){
       t=t_; qdir.init(it); qs.init(r*t.qlen()); 
    };

    QSet(const CRef<TQ> &A) { init(A); };

    QSet(const char *F, int L,
       const CRef<TQ> &A, const ctrIdx &ica, 
       const CRef<TQ> &B, const ctrIdx &icb, wbperm *P=NULL
    ){ init(F,L,A,ica,B,icb,P); };

    QSet& init(const CRef<TQ> &A) {
       if (!A.cgr) init();
       else { init(*A.cgr);
          if (A.gotPerm()) Permute(A.cgp);
          if (A.gotConj()) Conj();
       }
       return *this;
    };

    QSet& init(const QType t_ = QTYPE_UNKNOWN) { t=QTYPE_UNKNOWN; 
       t=t_; if (qs.len  ) qs.init();
             if (qdir.len) qdir.init();
       return *this;
    };

    QSet& init(const QSet &B) {
       t=B.t; qs=B.qs; qdir=B.qdir; return *this; };

    QSet& init(
       const QType &t_, const qset<gTQ>&qs_,
       unsigned iout=0,
       char isref=0
    );

    QSet& init(
       const QType t_, const TQ* qs_, unsigned r, unsigned N,
       unsigned iout=0
    );

    QSet& init(
       const QType t_, const TQ* qs_, unsigned r, unsigned N,
       const IndexLabels &it = IndexLabels()
    );

    QSet& initZ(const QType &t_, const TQ *q, char iflag=0) {
       unsigned m=t_.qlen();
       t=t_; qdir.init2val(2, iflag ? -1 : +1);

       qs.init(2*m); memcpy(qs.data, q, m*sizeof(TQ));
       t.getDual(q, qs.data+m);

       return *this;
    };

    QSet& init2(const QType &t_, const TQ *q=NULL) {
       unsigned m=t_.qlen();
       t=t_; qdir.init_iout(FL,2,2); qs.init(2*m);
       if (q) { unsigned l=m*sizeof(TQ);
          memcpy(qs.data,   q, l);
          memcpy(qs.data+m, q, l);
       }
       return *this;
    };

    QSet& init2(const QType &t_,
       const qset<TQ> &q1, const qset<TQ> &q2, const char *qds="++") {

       unsigned m=t_.qlen(), l=m*sizeof(TQ);
       if (q1.len<m || q2.len<m) wblog(FL,
          "ERR %s() invalid qset length (%s: %d,%d/%d)",
          FCT,t_.toStr().data,q1.len,q2.len,m
       );
       t=t_; qdir=qds; qs.init(2*m);
       memcpy(qs.data,   q1.data, l);
       memcpy(qs.data+m, q2.data, l);
       return *this;
    };

    QSet& init3(const QType &t_) {
       t=t_; qs.init(3*t_.qlen());
       qdir.init_iout(FL,3,3); return *this;
    };

    QSet& init3(const QType &t_,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J
    ){
       unsigned n=t_.qlen();
       if (J1.len!=n || J1.len!=J2.len || J1.len!=J.len) {
          if (J1.len!=n || J2.len!=n || J.len!=n) wblog(FL,
             "ERR %s() severe length inconsistency [%d %d %d; %d]", FCT,
             J1.len,J2.len,J.len,n
          ); 
       }

       t=t_; qs.init(J1,J2,J);
       qdir.init_iout(FL,3,3);
    };

    QSet& init(const char *F, int L,
       const CRef<TQ> &A, const ctrIdx &ica, 
       const CRef<TQ> &B, const ctrIdx &icb, wbperm *Pcgd=NULL
    );

    QSet& init_str(const char *F, int L, const char *s);

    QSet& Conj() { qdir.Conj(); return *this; };

    unsigned rank(const char *F=NULL, int L=0) const {
       unsigned n=t.qlen();
       if (!n || qs.len%n) wblog(F_L,
          "ERR %s() invalid qs (%d@%d)",FCT,qs.len,n);
       return (qs.len/n);
    };

    QSet& permute(QSet &B, const wbperm &P, char iflag=0) const;

    QSet& Permute(const wbperm &P, char iflag=0) {
       QSet<TQ> X; this->save2(X);
       return X.permute(*this,P,iflag);
    };

    QSet& save2(QSet &X) {
       X.t=t; qdir.save2(X.qdir); qs.save2(X.qs);
       return X;
    };


    QSet& operator=(const QSet &B) { return init(B); };

    bool operator==(const QSet &B) const {
       return (t==B.t && qs==B.qs && qdir==B.qdir);
    };

    bool operator!=(const QSet &B) const { 
       return (t!=B.t || qs!=B.qs || qdir!=B.qdir);
    };

    bool isEmpty(char lflag=0) const {
       if (qs.len || qdir.len) {
          if ((t==QTYPE_UNKNOWN && !lflag) || qs.len!=qdir.len*t.qlen())
             wblog(FL,"ERR %s() got invalid QSet: %s",FCT,toStr().data);
          return 0;
       }
       return (t==QTYPE_UNKNOWN || lflag ? 1 : 0);
    };

    bool isStd3() const {
       if (qdir.len!=3) return 0;
       if (qdir[0]<=0 || qdir[1]<=0 || qdir[2]>=0) return 0;
       if (qs.len%3) wblog(FL,"ERR %s() %s !?",FCT,toStr().data);
       return 1;
    };

    bool isScalar() const;
    bool isJ1Symbol() const;

    bool isAbelian() const { return (this->t.isAbelian()); };

    int checkQ_abelian(const char *F, int L) const {
       if (t.type==QTYPE_U1) return checkQ_U1(F,L);
       if (t.type<QTYPE_ASEP) wblog(F,L,
          "WRN %s() not yet implemented",FCT);
       return -1;
    };

    QSet& Sort(wbperm *cgp=NULL, char *conj=NULL, char iflag=0);
    bool isSorted() const;

    int checkQ_U1(const char *F, int L) const;
    int checkQ_SU2(const char *F, int L) const;

    wbstring QStrS(const wbperm *cgp=NULL, char sep='|') const;
    wbstring QStr() const;

    wbstring toStr(const char *istr=NULL) const;
    wbstring toTag() const;

    mxArray* toMx() const;

    QType t;
    qset<TQ> qs;

    QDir qdir;

 protected:
 private:

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;
};


template <class TQ>
class QHash {

  public: 

    size_t operator()(const QSet<TQ> &S) const {

       unsigned i, j=1, r=S.qdir.len, nq=S.qs.len,
          l=1+(( nq + r + (nq || r ? -1 : 0) )/sizeof(long));
       unsigned long  h=5381U, x[l];
       char *s=(char*)x;

       const TQ *qs=S.qs.data;

       if ((!r && nq) || (r && nq%r)) wblog(FL,
          "ERR %s() got QSet inconsistency (%d/%d)",FCT,nq,r);
       x[l-1]=0;
       
       s[0]=S.t.type;
       s[1]=S.t.sub;

       for (i=0; i<nq; ++i) {
           s[++j]=char(qs[i]);
           if (TQ(s[j])!=qs[i]) wblog(FL,
              "ERR %s() down conversion alters value (%g/%g) !??",
              FCT,double(TQ(s[j])),double(qs[i])
           );
       }

       memcpy(s,S.qdir.data,r);

       for (i=0; i<l; ++i) {
          h ^= ((h<<6) + (h>>2)) + x[i];
       }

       return h;
    };
};







template <class TD>
class cdata : public wbsparray<TD> {


 public:

    cdata& init(SPIDX_T l=0) { 
       wbsparray<TD>::init(l);
       return *this;
    };

    cdata& init(
       SPIDX_T D1, SPIDX_T D2, SPIDX_T D,
       unsigned M=1
    ){ 
       if (M==1)
            wbsparray<TD>::init(D1,D2,D);
       else wbsparray<TD>::init(D1,D2,D,M);
       return *this;
    };

    cdata& init(const wbvector<SPIDX_T> &S) {
       wbsparray<TD>::init(S); return *this;
    };

    cdata& init(const cdata &C) { 
       wbsparray<TD>::init(C);
       return *this;
    };

    template <class TQ>
    cdata& init(const CRef<TQ> &R, char full=1);

    template <class T2>
    cdata& init(wbsparray<T2> &C) { 
       wbsparray<TD>::init(C);
       return *this;
    };

    cdata& init(
       const char *F, int L, const mxArray *a, unsigned k,
       const QType &q, unsigned r=-1
    );

    cdata& operator=(const cdata &C) { 
       wbsparray<TD>::init(C);
       return *this;
    };

    cdata& initScalar(double x=1) { 
       wbsparray<TD>::initScalar(x);
       return *this;
    };

    cdata& initIdentity(const QType &q, SPIDX_T d, TD dval=1);
    cdata& initIdentity(const wbvector<SPIDX_T> &S);

    bool isScalar() const {
       return wbsparray<TD>::isScalar(FL,'z');
    };

    unsigned len() const { return this->SIZE.len; };

    unsigned getOM(
       const char *F=NULL, int L=0, const char *istr=NULL) const;

    cdata& Conj() { return *this; };

    cdata& permute(cdata &B, const wbperm &P0, char iflag=0) const;
    cdata& Permute(const wbperm &P0, char iflag=0) {
       cdata X; this->save2(X);
       return X.permute(*this,P0,iflag);
    };

    cdata& Kron(const cdata &B){
        wbsparray<TD>::Kron(B,'N','N','k');
        return *this;
    };
    
    TD contract(
       const char *F, int L, const wbindex &ica,
       const cdata &B, const wbindex &icb, cdata &C,
       cdata *Cx=NULL, wbperm *P=NULL
    ) const;

    double SkipTiny(const char *F=NULL, int L=0);

    TD NormSignC(
       const char *F=NULL, int L=0, unsigned m=0,
       TD eps=CG_EPS1,
       TD eps2=CG_EPS2
    );

    void info() const {
       printf("\n data: %s\n",this->sizeStr().data);
    };


 protected:
 private:
};



template <class TQ, class TD>
class CData : public QSet<TQ> {

 public:


    CData() : cstat(CGD_DEFAULT) { init(); };

    CData(const CData &C)
     : QSet<TQ>((QSet<TQ>&)C), cgd(C.cgd), cstat(C.cstat) {};

    CData(const QSet<TQ> &Q) : QSet<TQ>(Q) {};

    CData(const CRef<TQ> &R) { init(R); }


    void init() {
       this->t.init(); this->qs.init(); this->qdir.init();
       cgd.init(); cstat=CGD_DEFAULT;
    };

    CData& init(const CRef<TQ> &S);

    CData& initAbelian(const QSet<TQ> &Q) {
       this->t=Q.t; this->qs=Q.qs; this->qdir=Q.qdir;
       cstat.init(CGD_ABELIAN); cgd.init();
       return *this;
    };

    CData& init3(const QType &t_,
       SPIDX_T D1, SPIDX_T D2, SPIDX_T D,
       unsigned M=1
    ){
       QSet<TQ>::init3(t_);
       cgd.init(D1,D2,D, M); cstat=CGD_DEFAULT;
       return *this;
    };

    CData& init3(const QType &t_,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J,
       unsigned M=1
    ){
       SPIDX_T D1,D2,D;
       QSet<TQ>::init3(t_,J1,J2,J);
       cstat.init(CGD_DEFAULT);

       D1=t_.QDim(J1.data);
       D2=t_.QDim(J2.data);
       D =t_.QDim(J .data); cgd.init(D1,D2,D,M);

       return *this;
    };

    template <class T2>
    CData& init(const CData<TQ,T2> &B) {
       this->t=B.t; this->qs.init(B.qs); this->qdir=B.qdir;
       cgd.wbsparray<TD>::init(B.cgd); cstat=B.cstat;
       return *this;
    };

    CData& RefInit(
       const char *F, int L, const CData<TQ,TD> &B, char bare=0);

    CData& RefInit_auxtr(
       const char *F, int L,
       const wbvector<SPIDX_T> &S, const wbvector<RTD> &cgt
    );

    CData& Reduce2Ref(const char *F=0, int L=0);

    CData& LoadRef(const char *F=0, int L=0);

    int init(
       const char *F, int L, const mxArray *S, unsigned k,
       char bflag=0, char check=1
    );

    TQ* qptr(unsigned k) const {
       unsigned n=this->t.qlen();
       if (this->qs.len%n || k>=this->qs.len/n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d <> %d)",
          FCT,this->qs.len,n,k);
       return (this->qs.data+k*n);
    };

    CData& operator=(const CData &B) {
       this->t=B.t; this->qs=B.qs; this->qdir=B.qdir;
       cgd=B.cgd; cstat=B.cstat;
       return *this;
    };

    CData& operator=(const QSet<TQ> &B) {
       this->t=B.t; this->qs=B.qs; this->qdir=B.qdir;
       return *this;
    };

    CData& save2(CData &B) {
       B.t=this->t; this->t.init(); B.cstat=cstat; cstat=CGD_DEFAULT;
       this->qdir.save2(B.qdir);
       this->qs.save2(B.qs); cgd.save2(B.cgd);
       return B;
    };

    bool olderThan(const CData &B) const {
       if ((QSet<TQ>&)*this!=(QSet<TQ>&)B) wblog(FL,
          "ERR %s() got incompatible CData\n   %s\n<> %s",
          FCT, toStr().data, B.toStr().data
       );
       return cstat.olderThan(B.cstat);
    };

    bool operator!=(const CData &B) const { return !((*this)==B); };
    bool operator==(const CData &B) const {
       return (
          this->t==B.t && this->qs==B.qs && this->qdir==B.qdir &&
          cgd==B.cgd && cstat.sameAs(B.cstat,0)
       );
    };

    bool sameAs(const CData &B, char lflag=0) const {
       if (!lflag) { return (*this)==B; }
       if (this->t!=B.t || this->qs!=B.qs || this->qdir!=B.qdir) return 0;

       if (lflag=='l') lflag=1; else
       if (lflag=='L') lflag=2;

       if (!(cstat==CGD_REF_INIT ^ B.cstat.t==CGD_REF_INIT)) {
          if (lflag<2) { if (cgd!=B.cgd) return 0; }
          if (cstat.t!=B.cstat.t) return 0;
       }
       return 1;
    };

    bool operator!=(const QSet<TQ> &B) const {
       return (this->QSet<TQ>::operator!=(B));
    };
    bool operator==(const QSet<TQ> &B) const {
       return (this->QSet<TQ>::operator==(B));
    };

    bool operator!=(const cgdStatus &b) const { return (cstat!=b); };
    bool operator==(const cgdStatus &b) const { return (cstat==b); };

    bool isEmpty() const {
       return (QSet<TQ>::isEmpty() && cgd.isEmpty());
    };

    bool isAbelian(const char *F=NULL, int L=0) const {
       if (this->t.isAbelian()) {
          if (cstat.t==CGD_ABELIAN) {
              if (cgd.D.len || cgd.SIZE.len) wblog(F_L,
                 "ERR %s() invalid abelian CData\n%s",FCT,toStr().data);
              return 1;
          }
          if (!isScalar() || (cgd.D.len && (cgd.D.len>1 || cgd.D[0]!=1)))
          wblog(F_L,"ERR %s() invalid abelian CData\n%s",toStr().data);
          return 1;
       }
       return 0;
    };

    bool isScalar() const {
       if (!cgd.SIZE.len) {
          if (cgd.D.len==1 || (!cgd.D.len && CGD_ABELIAN)) return 1;
       }
       else if (cgd.isScalar()) return 1;
       return 0;
    };

    TD getScalar(const char *F=NULL, int L=0) const {
       if (!isScalar()) wblog(F_L,
          "ERR %s() got non-scalar CData\n%s",FCT,toStr().data);
       return (cgd.D.len ? cgd.D.data[0] : TD(1));
    };

    unsigned rank(const char *F=NULL, int L=0) const;

    unsigned rankS() const { return cgd.rank(); };

    bool isrank(unsigned r, unsigned *r_=NULL) const {
       unsigned l=rank(FL); if (r_) { (*r_)=l; }
       return (l!=r);
    };

    template <class T>
    wbvector<T>& getSize(
       wbvector<T> &S, char bare=0, const wbperm *cgp=NULL) const;

    bool isRefInit() const {
       if (cstat!=CGD_REF_INIT) return 0;
       const unsigned r=cgd.SIZE.len, m=cgd.D.len;
       const SPIDX_T *s=cgd.SIZE.data;
       if (!r || (m>1 && m!=s[r-1])) wblog(FL,
          "ERR %s() got invalid CData size ref (%s @ D.len=%d)",
          FCT, cgd.sizeStr().data, cgd.D.len);
       return 1;
    };

    bool checkValidStat(const char *F=0, int L=0) const {
       if (this->t.isAbelian()) {
          if (cstat.hasID('l') || cstat.t!=CGD_ABELIAN) {
             if (F) wblog(F,L,"ERR %s()\n%s\n%s",
                FCT,cstat.toStr('V').data,toStr().data);
             return 0;
          }
       }
       else if (!cstat.hasID() || cstat.t<=CGD_ABELIAN) {
          if (F) wblog(F,L,"ERR %s()\n%s\n%s",
             FCT,cstat.toStr('V').data,toStr().data);
          return 0;
       }
       return 1;
    };

    bool sameType(const CData &S, const char *F=NULL, int L=0) const;

    template <class DB>
    bool sameSizeR(
       const CData<TQ,DB> &B,
       unsigned *r=NULL, const char *F=NULL, int L=0) const;

    template <class DB>
    bool sameSizeR(
       const cdata<DB> &B,
       unsigned *r=NULL, const char *F=NULL, int L=0) const;

    template <class T2>
    int sameSizeR(const char *F, int L, const wbvector<T2> &S) const;

    bool isSymmmetric() const;

    unsigned checkOM(const char *F=NULL, int L=0) const {
       if (this->qdir.len) {
          unsigned r=this->rank(F_L);
          if (this->t.hasMP(r)) {
             if (this->cgd.SIZE.len>r) {
                SPIDX_T i=this->cgd.SIZE[r]; if (!i)
                   wblog(FL,"ERR %s() got OM=%d !?",FCT,i);
                return i;
             }
             return 1;
          }
          else if (this->cgd.SIZE.len>r) wblog(FL,
          "ERR %s() got CData with invalid OM\n%s",FCT,toStr().data);
       }
       return 0;
    };

    unsigned gotOM(const char *F=NULL, int L=0) const {
       if (this->qdir.len) {
          unsigned l=cgd.SIZE.len, r=this->rank(F_L);
          if (l==r+1) { 
             SPIDX_T i=cgd.SIZE[r];
             if (!i || (i>1 && !this->t.hasMP(r))) wblog(FL,
                "ERR %s() got invalid OM setting (i=%d,r=%d) !?\n%s",
                FCT,i,r,toStr().data);
             return i;
          }
          else if (l && l!=r) wblog(FL,
             "ERR %s() invalid OM setting\n%s",FCT,toStr().data
          );
       }
       return 0;
    };

    unsigned getOM(const char *F=NULL, int L=0) const {

       if (cgd.SIZE.len) {
          unsigned m=gotOM(F,L); return (m ? m : 1);
       }

       if (this->qdir.len==2 && cgd.isDiag()) { return 1; }
       if (cgd.isEmpty() && cstat==CGD_ABELIAN) { return 1; }


       if (!isEmpty()) {
          this->rank(F_L);
          cgd.wbsparray<TD>::checkSize(F_LF);
       }
       return 0;
    };

    CData& initOM(const char *F, int L);

    double normDiff(
       const char *F, int L, const CData &S,
       double eps=1E14
    ) const;

    double norm2(unsigned k=-1) const;

    wbvector<RTD>& trace(
       const char *F, int L, wbvector<RTD> &cgt) const;

    wbvector<RTD> trace(const char *F=NULL, int L=0) const {
       wbvector<RTD> cgt; return trace(F,L,cgt);
    };

    CData& AddMultiplicity(const char *F, int L, const cdata<TD> &c);

    CData& AddMultiplicity(const char *F, int L, const CData &C) {
       if (this->t!=C.t || this->qs!=C.qs) wblog(F_L,
          "ERR %s() incompatible CData\n[%s] <> [%s] (%s, %s)",
          FCT, this->toStr().data, C.toStr().data,
          sizeStr().data, C.sizeStr().data
       );
       return AddMultiplicity(F_L,C.cgd);
    };

    bool isBasicCG(const char *F=NULL, int L=0) const;

    SPIDX_T getBasicCGSize(unsigned *m=NULL) const;

    wbsparray<TD>&
    getBasicCG(unsigned k, wbsparray<TD> &a) const;

    int getCG_set(
       const char *F, int L, wbIndex &Idx, wbsparray<TD> &a) const;

    CData& getJ1Symbol(const char *F, int L, CData<TQ,TD> &C) const;

    int checkConsistency(const char *F=NULL, int L=0) const;
    int checkNormSign(const char *F=NULL, int L=0, char xflag=0) const;

    int checkQ(const char *F, int L,
       const QSet<TQ> &Q, const char* istr=NULL) const;

    bool checkAdditivityZ(const char *F, int L, 
       const wbMatrix<double> &Z1,
       const wbMatrix<double> &Z2, const wbMatrix<double> &Z,
       TD eps=CG_EPS2
    ) const;

    CData& Sort(const char *F=NULL, int L=0);

    double SkipTiny(const char *F=NULL, int L=0);

    CData& Permute(const wbperm &P, char iflag=0) {
       CData<TQ,TD> X; this->save2(X);
       return X.permute(*this,P,iflag);
    };

    CData& permute(CData &B, const wbperm &P, char iflag=0) const {
       QSet<TQ>::permute((QSet<TQ>&)B,P,iflag);

       cgd.permute(B.cgd,P,iflag);
       B.cstat=cstat;

       return B;
    };


    void info(const char *istr=NULL, const char *F=NULL, int L=0) const;

    wbstring dStr() const { return this->qdir.toStr(); };
    wbstring qStr() const { return this->t.toStr(); };


    wbstring toStr(char vflag=0) const;
    wbstring sizeStr() const;

    mxArray* toMx() const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    mxArray* mxCreateCell(unsigned m, unsigned n) const {
       wblog(FL,"ERR %s(%d,%d)",FCT,m,n); return 0; };
    void add2MxCell(mxArray *S, unsigned i, char tst=0) const {
       wblog(FL,"ERR %s(%lX,%d,%d)",FCT,S,i,tst); };

    void put(const char *vname, const char *ws="caller"
    ) const { put(0,0,vname,ws); };

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx();
       int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F_L,
          "I/O putting variable '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };


    cdata<TD> cgd;

    cgdStatus cstat;

 private:
};


template <class TQ>
class CRef {

 public:

    CRef() : cgr(NULL), conj(0) {};

    CRef(const CRef &S) { *this=S; };


    CRef& operator=(const CRef &S) {
       cgr=S.cgr; rtype=S.rtype;
       cgw=S.cgw; cgp=S.cgp; conj=S.conj;
       return *this;
    };

    CRef& init(CGR_TYPE rt_=CGR_DEFAULT) {
       cgr=NULL; cgw.init(); cgp.init(); conj=0; rtype=rt_;
       return *this;
    };

    CRef& init(const CRef &S) { *this=S; return *this; };

    CRef& init(
       const char *F, int L, const mxArray *a, unsigned k,
       char refC=0, const QSet<TQ> *Q=NULL
    );

    CRef<TQ>& initDecompose(
       const char *F, int L, cdata__ X,
       CDATA_TQ *cgr_=NULL
    );

    CRef<TQ>& initIdentityR(const char *F, int L,
       const QType &q, const TQ *qs, unsigned dim=-1, char xflag=0);

    CRef<TQ>& initIdentityZ(const char *F, int L,
       const QType &q, const TQ *qs, unsigned dim=-1, char xflag=0);

    CRef<TQ>& Reduce2Identity(char xflag=0);

    wbvector<RTD>& trace(
       const char *F, int L, wbvector<RTD> &cgt) const;

    wbvector<RTD> trace(const char *F, int L) const {
       wbvector<RTD> cgt; return trace(F,L,cgt);
    };

    RTD trace() const;

    CRef& initAbelian(){
       cgr=NULL; cgw.init(); rtype=CGR_ABELIAN;
       cgp.init(); conj=0;
       return *this;
    };

    CRef& initCtrScalar(RTD x=1){
       cgr=NULL; cgw.init(1,&x); rtype=CGR_CTR_SCALAR;
       cgp.init(); conj=0;
       return *this;
    };

    CRef& initBase(const CDATA_TQ* cgr_, unsigned i, unsigned m);

    const CRef& LoadRef(const char *F=0, int L=0) const {
       if (!isRefInit(F,L)) return *this;
       if (!cgr) wblog(FL,
          "ERR %s() got cgr=null !?\n%s",FCT,toStr().data);
       ((CDATA_TQ*)cgr)->LoadRef(F,L);
       return *this;
    };

    bool isEmpty() const {
        if (!cgr || cgr->isEmpty()) {
           if (cgw.len>1 || cgp.len) wblog(FL,
              "ERR %s() invalid abelian cgref (%d)",FCT,toStr().data);
           return (rtype==CGR_DEFAULT && !cgw.len);
        }
        return 0;
    };

    bool isRefInit(const char *F=0, int L=0) const {
       if (cgr && cgr->isRefInit()) {
          if (rtype!=CGR_DEFAULT) wblog(F_L,"ERR %s() "
             "got REF_INIT mismatch\n%s",FCT,toStr().data);
          return 1;
       }
       return 0;
    };

    bool isAbelian(const char *F=NULL, int L=0) const;

    bool isScalar() const;

    bool isScalarCtr() const {
       if (cgr || rtype!=CGR_CTR_SCALAR || cgw.len!=1) return 0;
       return 1;
    };

    unsigned numel() const;

    unsigned Size(unsigned k) const;
    unsigned Size(unsigned k, unsigned r) const;

    unsigned qdim(unsigned k, const QType *qt=NULL) const;

    int checkAbelian(const char *F=NULL, int L=0) const;

    bool isDiagCGC(RTD eps=1E-14) const;
    bool isIdentityCGC(RTD *nrm=NULL, RTD eps=1E-14) const;

    char sameSize(const CRef &B, char strict=0);
    char sameUptoFac(const CRef &B, RTD *fac=NULL, RTD eps=1E-12) const;
    char sameAs(const CRef &B, RTD eps=1E-12) const;

    bool sameQSet(const QSet<TQ> &Q) const;
    bool sameQSet(const CRef &B) const;
    bool sameQDir(const char *F, int L, const QDir &qd) const;

    bool sameQDir(const char *F, int L, const CRef& B) const {
       QDir qdb; B.get_qdir(qdb);
       return sameQDir(F,L,qdb);
    };

    int gotSameCData(
      const char *F, int L, const CRef& B, char c2flag=0
    ) const;

    RTD NormStd(const char *F, int L, unsigned r=-1);

    RTD normExt(const char *F=NULL, int L=0) const;

    int SortDegQ(const char *F=NULL, int L=0, QSet<TQ> *QS=NULL);

    bool isSortedDQ(QSet<TQ> *QS, wbperm &pxt, char &cxt) const;
    bool isSortedDQ(QSet<TQ> *QS=NULL) const {
       wbperm pxt; char cxt=0;
       return isSortedDQ(QS,pxt,cxt);
    };

    char got3(const char *F, int L,
       const QType &q, const qset<TQ> &J12, const qset<TQ> &J3
     ) const;

    bool operator==(const CRef &B) { return (sameAs(B)==0 ? 1 : 0); };
    bool operator!=(const CRef &B) { return (sameAs(B)==0 ? 0 : 1); };

    RTD scalarProd(const char *F, int L, const CRef &B) const;

    RTD normDiff2(const char *F, int L, const CRef &B) const;
    RTD normDiff (const char *F, int L, const CRef &B) const {
       return SQRT(F,L,normDiff2(B)); };

    RTD sumTimesEl(const char *F, int L, const CRef &B) const;

    int check(const char *F, int L) const;
    int check(const char *F, int L,
      const QType &q, const QDir &qdir) const;

    RTD norm2(char xflag=0) const;

    RTD NormSignR(
       const char *F=NULL, int L=0,
       RTD eps=CG_EPS1,
       RTD eps2=CG_EPS2
    );

    unsigned getOM(const char *F=NULL, int L=0) const;
    unsigned checkOM(const char *F=NULL, int L=0) const {
       return (cgr ? cgr->checkOM(F,L) : 0);
    };

    unsigned rankS(char lflag=0) const;
    unsigned rank(const char *F=NULL, int L=0, char lflag=0) const;

    template <class T>
    wbvector<T>& getSize(wbvector<T> &S, char bare=0) const;

    CRef& Permute(const wbperm &P, char iflag=0, char isnew=0);

    CRef& permute(
       const wbperm &P, CRef &B, char iflag=0, char isnew=0
     ) const {
       B=(*this); return B.Permute(P,iflag,isnew);
    };

    CRef& Conj() {
       if (cgr) { conj=(conj ? 0 : 1); }
       return *this;
    };

    int HConjOpScalar();

    bool gotConj() const { return (conj!=0); };

    bool gotPerm(const char *F=NULL, int L=0) const {
       if (cgp.len) {
          if (!cgr) wblog(F_L,"ERR %s() "
             "got invalid perm (len=%d; cgr=NULL)",FCT, cgp.len);
          if (cgp.len!=cgr->qdir.len) wblog(F_L,"ERR %s() "
             "got invalid perm (len=%d/%d)",FCT, cgp.len, cgr->qdir.len);
          return !cgp.isIdentityPerm(F_L,cgr->qdir.len);
       }
       return 0;
    };

    bool anyTrafo(char pflag=1) const {
       if (conj || cgp.len) return 1;
       if (pflag && !cgp.isIdentityPerm()) return 1;
       return 0;
    };

    bool anyTrafo(char &ip, char &ic, char pflag=1) const {
       if (pflag)
            { ip=(cgp.len && !cgp.isIdentityPerm()); }
       else { ip=(cgp.len ? 1 : 0); }
       ic=(conj!=0);
       return (ip || ic);
    };

    bool gotnoTrafos() const {
       return !anyTrafo(); };

    unsigned getP(unsigned k) const {
       if (cgp.len) {
          if (k>=cgp.len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,k,cgp.len);
          return cgp.data[k];
       }
       return k;
    };

    wbperm& getP(wbperm &p, char iflag=0) const {
       if (!cgp.len)
            return p.init(rank(FL));
       else return p.init(cgp,iflag);
    };

    char stat() const { char i=0;
       if (rtype!=CGR_DEFAULT) i ^= 1;
       if (conj) i ^= 2;
       if (cgp.len && !cgp.isIdentityPerm()) i ^= 4;
       return i;
    };

    bool affectsQDir() const;

    qset<TQ>& adapt(qset<TQ> &qs, char iflag=0) const {
       if (cgp.len) qs.BlockPermute(cgp,iflag ? 0 : 'i');
       return qs;
    };

    QSet<TQ>& adapt(QSet<TQ> &Q, char iflag=0) const {
       if (cgp.len) Q.Permute(cgp,iflag ? 0 : 'i');
       if (conj) Q.qdir.Conj();
       return Q;
    };

    ctrIdx& adapt(ctrIdx &I, char iflag=0) const {
       if (cgp.len && !cgp.isIdentityPerm()) {
          unsigned i=0; wbvector<unsigned> x(I.len,I.data);
          if (I.len>cgp.len) wblog(FL,"ERR %s() "
             "unexpected index length (%d/%d)",FCT,I.len,cgp.len);
          if (iflag==0) {
             for (; i<I.len; ++i) { I[i]= cgp.el(x[i]); }}
          else {
             wbperm cpi(cgp,'i');
             for (; i<I.len; ++i) { I[i]= cpi.el(x[i]); }
          }
       }
       if (conj) I.Conj();
       return I;
    };

    QDir& get_qdir(QDir &qd) const {
       if (!cgr) { qd.init(); return qd; }

       if (cgp.len)
            cgr->qdir.permute(qd,cgp);
       else qd=cgr->qdir;

       if (conj) qd.Conj();

       return qd;
    };

    QDir get_qdir() const { QDir qd;
       return get_qdir(qd);
    };

    RTD safeCpy(const char *F, int L, const CRef &B, CRef &C) const;

    CRef& save2(CRef &S) {
       if (this!=&S) {
          S.cgr=cgr;
          cgw.save2(S.cgw); S.conj=conj; conj=0;
          cgp.save2(S.cgp); S.rtype=rtype; rtype=0;
       }
       return S;
    };

    wbstring qdir2Str(char vflag=0) const {

       if (!cgr || !cgr->qdir.len) return "";

       if (vflag<=0) return get_qdir().toStr();
       if (vflag==1 || vflag=='v') {
          unsigned l, n=2*(cgr->qdir.len)+5; char s[n];
          l=snprintf(s,n,"%s => %s",
             cgr->qdir.toStr().data, get_qdir().toStr().data);
          if (l>=n) wblog(FL,
             "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
          return s;
       }
       else {
          unsigned l, n=2*(cgr->qdir.len)+5; char s[n];
          l=snprintf(s,n,"%s[%s]%s",
             cgr->qdir.toStr().data, (cgp+1).toStr(1,"").data,
             gotConj()?"*":"");
          if (l>=n) wblog(FL,
             "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
          return s;
       }
    };

    wbstring toStr() const;
    wbstring sizeStr() const;

    wbstring statStr(char vflag=0) const;

    wbstring qStr() const { return cgr ? cgr->qStr() : ""; };
    wbstring QStr() const { 
       return cgr ? cgr->QStrS(&cgp) : ""; 
    };


    mxArray* toMx(char flag=0) const;

    mxArray* toMX() const { return toMx('f'); };

    mxArray* mxCreateStruct(unsigned m, unsigned n, char flag=0) const;
    void add2MxStruct(mxArray *S, unsigned i, char flag=0) const;


    const CDATA_TQ *cgr;


    wbvector<RTD> cgw;

    wbperm cgp;

    char conj;

    cgrType rtype;
};



#define MAP32  std::map < qset<TQ>, wbvector< CRef<TQ> > >
#define MAP31  std::map < qset<TQ>, MAP32 >
#define MAP30  std::map < QType, MAP31 >

#define MAPA1 std::map < qset<TQ>, CDATA_TQ >
#define MAPA0 std::map < QType, MAPA1 >

template <class TQ>
class CStore {

  public:

    CStore() {
       BUF.max_load_factor(0.2);
    };


    void clear() {
       map3.clear();
       BUF.clear();
    };

    typedef typename MAP30::iterator iMAP30;
    typedef typename MAP31::iterator iMAP31;
    typedef typename MAP32::iterator iMAP32;

    typedef typename MAP30::const_iterator ciMAP30;
    typedef typename MAP31::const_iterator ciMAP31;
    typedef typename MAP32::const_iterator ciMAP32;

    typedef typename MAPA0::const_iterator ciMAPA0;
    typedef typename MAPA1::const_iterator ciMAPA1;

    void init() {};

    void getCData(const char *F, int L, const QType &type,
       const qset<TQ> &J1, const qset<TQ> &J2,
       wbvector < CDATA_TQ* > &S
    );
    
    const CDATA_TQ& getCData(const char *F, int L,
       const QType &q, const qset<TQ> &J1, const qset<TQ> &J2,
       const qset<TQ> &J, char force=0
    );

    const CDATA_TQ& getCData(const char *F, int L,
       const QVec &q, unsigned k,
       const qset<TQ> &J1, const qset<TQ> &J2,
       const qset<TQ> &J, char force=0
    );

    const CDATA_TQ& getIdentityC(const char *F, int L,
       const QType &t, const TQ *qs, unsigned dim=-1,
       char loadRC='l'
    );
    const CDATA_TQ& getIdentityZ(const char *F, int L,
       const QType &t, const TQ *qs, unsigned dim=-1,
       char loadRC='l'
    );

    void add_CGData_U1(
       const QType &type, const qset<TQ> &J1, const qset<TQ> &J2);
    void add_CGData_ZN(
       const QType &type, const qset<TQ> &J1, const qset<TQ> &J2);
    void add_CGData_P(
       const QType &type, const qset<TQ> &J1, const qset<TQ> &J2);

    inline void add_CGData(
       const QType &q, const qset<TQ> &J1,const qset<TQ> &J2
    );

    void add2BUF(const char *F, int L, const CDATA_TQ &S, char nflag=0);

    void testCG(const char *F=NULL, int L=0);

    unsigned add3(const QType &q, const mxArray* S, unsigned k=-1);
    unsigned add2BUF(const QType &q, const mxArray* S, unsigned k=-1);

    CDATA_TQ& getBUF(
       const char *F, int L, const QSet<TQ> &Q, char loadRC=0
    );

    void getQfinal(const char *F, int L, const QType &type,
       const qset<TQ> &J1, const qset<TQ> &J2, wbMatrix<TQ> &J,
       wbMatrix< const wbvector<CRef<TQ> >* > *ss=NULL
    );

    void getQfinal(const char *F, int L, const QVec &type,
       const qset<TQ> &q1, const qset<TQ> &q2, wbMatrix<TQ> &QQ,
       wbMatrix< const wbvector< CRef<TQ> >* > *SS=NULL
    );

    int getQfinal_zdim(
       const char *F, int L, const QVec &qvec,
       const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q,
       INDEX_T &s1, INDEX_T &s2, INDEX_T &s, unsigned *m=NULL
    );

    int getQfinal_zdim(
       const char *F, int L, const QVec &qvec,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q, 
       wbvector<INDEX_T> &s1, wbvector<INDEX_T> &s2, wbvector<INDEX_T> &s,
       wbvector<unsigned> *m=NULL,
       wbMatrix< const wbvector< CRef<TQ> >* > *S3=NULL
    );

    void getQfinal(const char *F, int L, const QVec &type,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, QMap<TQ> &M,
       const char cgflag=0
    );

    void Info(const char *F, int L);
    void printStatus(const char *F, int L);

    unsigned LoadStore(const char *F, int L, const char *file);

    mxArray* toMx(const QType &q=QTYPE_UNKNOWN) const;

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' to %s (%d)",vname,i,ws);
       if (F) wblog(F_L,
          "I/O putting variable '%s' to %s.",vname,ws);
       mxDestroyArray(a);
    };

    void put(const char *vname, const char *ws="caller") const {
       put(0,0,vname,ws);
    };

    CDATA_TQ* find_cdata(const QSet<TQ> &Q) {
       auto it = BUF.find(Q);
       return (it!=BUF.end() ? &(it->second) : NULL);
    };

    int find_map3(const QType &q, const qset<TQ> &J1, const qset<TQ> &J2
     ) const {
       qset<TQ> J12(J1,J2);
       if (J1.len!=J2.len || J1.len!=q.qlen()) wblog(FL,
          "ERR %s() qset length inconsistency (%d/%d/%d)",
          FCT,J1.len,J2.len,q.qlen());
       return find_map3(q,J12);
    };

    int find_map3(const QType &q, const qset<TQ> &J12) const {
       if (J12.len!=2*q.qlen()) {
          if (!q.isKnown()) wblog(FL,
             "ERR %s() got undefined symmetry",FCT); else
          wblog(FL,"ERR %s() "
             "qset length inconsistency (%d/2*%d)",FCT,J12.len,q.qlen()
          );
       }
       auto it1 = map3.find(q);
       if (it1!=map3.end()) {
          auto it2 = it1->second.find(J12);
          if (it2!=it1->second.end()) { return it2->second.size(); }
       }
       return -1;
    };


    map <QType,
      map <qset<TQ>,
        map <qset<TQ>,
          wbvector< CRef<TQ> >
        >
      >
    > map3;


    unordered_map <
       QSet<TQ>,
       CDATA_TQ,
       QHash<TQ>
    > BUF;




  protected:
  private:


    iMAP31 get_mp3_data(const char *F, int L, 
    const QType &q, const qset<TQ> &J1, const qset<TQ> &J2, char force=0){

       if (J1.len!=J2.len) wblog(FL,
          "ERR qset length inconsistency (%d/%d)",J1.len,J2.len);

       iMAP31 iJ12;
       qset<TQ> J12(J1,J2);

       iMAP30 itype = map3.find(q);
       if (itype!=map3.end()) {
          iJ12 = itype->second.find(J12);
          if (iJ12!=itype->second.end()) {
             if (!iJ12->second.size()) wblog(FL,
                "ERR %s() got empty map3 data !?",FCT);
             return iJ12;
          }
       }

       if (!force) {
          add_CGData(q,J1,J2);
          return get_mp3_data(F,L,q,J1,J2,'f');
       }

       wblog(F,L,
         "ERR failed to find given target space\n[%s] x [%s] => [%s]",
          J1.toStr().data, J2.toStr().data, J12.toStr().data);
       return iJ12;
    };

    void print_Coeff(const char *F, int L, const QType &type,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J,
       double c
    ){
       wblog(F,L," *  CG coeff. %-8s [%s, %s; %s] = %.6g",
       type.toStr(), J1.toStr().data, J2.toStr().data,
       J.toStr().data,c);
    }
};

   CStore<gTQ> gCG;

   CData<gTQ,double> EMPTY_CGSTORE;

   mxArray* map3_CreateStructMatrix(unsigned m, unsigned n);

   template <class TQ>
   void map3_add2MxStruct(mxArray *S, unsigned k,
      const QType &q, const qset<TQ> &J12);


template <class TM>
class cgc_contract_id : public wbvector<TM> {

  public:

    cgc_contract_id() { wbvector<TM>::init(); };


    template <class TQ, class TD>
    cgc_contract_id(
       const CData<TQ,TD> &A, const ctrIdx &ica,
       const CData<TQ,TD> &B, const ctrIdx &icb
    );

    template <class TQ, class TD>
    void extract(
       CData<TQ,TD> &A, ctrIdx &ica,
       CData<TQ,TD> &B, ctrIdx &icb) const;

  protected:
  private:

     template <class T1, class T2>
     unsigned set_val(const char *F, int L, T1* &s, const T2 &x) {
         s[0]=x;
         if (T2(s[0])!=x) wblog(F_L,"ERR got non-%s data entry (%g)",
            getName(typeid(T1)).data, double(x)
         );
         s+=1; return 1;
     };

     template <class T1, class T2>
     unsigned set_vec(const char *F, int L, T1* &s, const wbvector<T2> &x) {
         s[0]=x.len;
         if (unsigned(s[0])!=x.len) wblog(F_L,"ERR unexpected size "
            "(out of bounds %s @ %d)",getName(typeid(T1)).data,x.len
         );
         set_range(F_L,++s,x.len,x.data);
         return (1+x.len);
     };

     template <class T1, class T2>
     unsigned set_range(const char *F, int L, T1* &s, unsigned n, const T2 *x) {
         for (unsigned i=0; i<n; ++i) {
            s[i]=x[i];
            if (T2(s[i])!=x[i]) wblog(F_L,"ERR got non-%s data entry (%g)",
               getName(typeid(T1)).data, double(x[i])
            );
         }
         s+=n; return n;
     };

     template <class TQ, class TD>
     unsigned set_CData(const char *F, int L,
        char* &s, unsigned n, const CData<TQ,TD> &A, char tflag=1
     ){
        char *s0=s;
        unsigned r=A.rank(), l=4+A.qs.len+r+(4+r)*sizeof(unsigned);


        wbvector<size_t> S;

#ifndef NSAFEGUARDS
        { const cdata<TD> &a=A.cgd;
          if (a.IDX.dim1 && a.IDX.dim2!=a.SIZE.len) wblog(F_L,
             "ERR %s() cdata size mismatch !?? (a: %dx%d; %d; %s)",
             FCT, a.IDX.dim1, a.IDX.dim2, a.D.len, a.sizeStr().data);
          if (!a.D.len && A.cstat.t!=CGD_REF_INIT) wblog(FL,
             "ERR %s() got empty CGC data !?? (A: %dx%d; %d)",
             FCT, a.IDX.dim1, a.IDX.dim2, a.D.len
          );
        }
#endif

        if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

        if (tflag) {
           set_val(FL,s, A.t.type);
           set_val(FL,s, A.t.sub );
        }

        set_vec(FL,s, A.qs);
        set_vec(FL,s, A.qdir);

        unsigned *u=(unsigned*)s;

        set_vec(FL,u, A.getSize(S,'b'));

        s+=(u-(unsigned*)s)*sizeof(unsigned); l=s-s0;
        if (l>=n) wblog(FL,
           "ERR %s() string out of bounds (%d/%d)",FCT,l,n
        );

        return l;
     };

     template <class T1, class T2>
     void get_range(
        const char *F, int L, const T1* &s, unsigned n, T2 *x) const {

        for (unsigned i=0; i<n; ++i) {
           x[i]=s[i];
           if (T1(x[i])!=s[i]) wblog(F_L,"ERR got non-%s data entry (%g)",
              getName(typeid(T1)).data, double(x[i])
           );
        }
        s+=n;
     };

     template <class T1, class T2>
     void get_vec(
        const char *F, int L, const T1* &s, wbvector<T2> &x) const {
        if (!s || int(s[0])<0) wblog(FL,
           "ERR %s() got invalid vec.len=%d",FCT,s?s[0]:-1);
        x.init(s[0]); ++s;
        get_range(F,L,s,x.len,x.data);
     };

     template <class T1, class T2>
     void get_val(
        const char *F, int L, const T1* &s, T2 &x) const {
        get_range(F,L,s,1,&x);
     };

     template <class TQ, class TD>
     unsigned get_CData(const char *F, int L,
        const char* &s, CData<TQ,TD> &A, char tflag=1) const {

        A.init();
        if (tflag) {
           A.t.initx(QTYPE_SET(s[0]),unsigned(s[1]));
           s+=2;
        }

        get_vec(F_L,s, A.qs);
        get_vec(F_L,s, A.qdir);

        const unsigned *u = (const unsigned*)s;
        get_vec(F_L,u, A.cgd.SIZE);
        A.cstat.init(CGD_BSZ_INIT);



        s+=(u-(unsigned*)s)*sizeof(unsigned);

        return 0;
     };

};


template <class TM>
class MHash {

  public: 

    size_t operator()(const cgc_contract_id<TM> &S) const {

       MTI h=5381U, *const x=S.data;
       for (unsigned i=0; i<S.len; ++i) {
          h ^= ((h<<6) + (h>>2)) + x[i];
          h ^= ((h<<6) + (h>>2)) + x[i];
       }
       return h;
    };
};


template <class TQ, class TD>
class x3map {

  public:

    x3map() : cgr(NULL), conj(0), rtype(CGR_DEFAULT) {};

    x3map(const char *F, int L,
       const CRef<TQ> &A_, const ctrIdx &ica,
       const CRef<TQ> &B_, const ctrIdx &icb, wbperm &cgp, char xflag=0
     ) : cgr(NULL), conj(0), rtype(CGR_DEFAULT)
     { initCtr(F_L,A_,ica,B_,icb, cgp, xflag); };


    x3map<TQ,TD>& init() {
       a.init(); b.init(); c.init(); pab.init(); x3.init();
       cgr=NULL; conj=0; rtype=CGR_DEFAULT;
       return *this;
    };

    x3map<TQ,TD>& initCtr(const char *F, int L,
       const CRef<TQ> &A_, const ctrIdx &ica_,
       const CRef<TQ> &B_, const ctrIdx &icb_, wbperm &cgp,
       char xflag=1
    );

    bool isEmpty() const {
       return (
          a.isEmpty() && b.isEmpty() && c.isEmpty() &&
          !pab.len && x3.isEmpty()
       );
    };

    wbstring toStr(char vflag=0) const;
    wbstring x3Str(char vflag=0) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;

    void add2MxStruct(
       mxArray *S, unsigned l,
       const cgc_contract_id<MTI> *idc=NULL) const;

    mxArray* toMx(const cgc_contract_id<MTI> *idc=NULL) const {
       mxArray *S=mxCreateStruct(1,1);
       add2MxStruct(S,0,idc); return S;
    };


    CData<TQ,RTD> a, b, c;

    const CDATA_TQ* cgr;

    wbperm pab;


    char conj;

    cgrType rtype;

    wbarray<TD> x3;

  protected:
  private:
};


template <class TQ>
class X3Map {

  public:

    X3Map() {
       XBUF.max_load_factor(0.2);
    };

    void clear() {
       XBUF.clear();
    };
   
    CRef<TQ>& contract(const char *F, int L,
       const CRef<TQ> &A, const ctrIdx &ica,
       const CRef<TQ> &B, const ctrIdx &icb, CRef<TQ> &C
    );

    CRef<TQ>& contractDQ(const char *F, int L,
       const CRef<TQ> &A, CRef<TQ> &B
    );

    mxArray* toMx() const;

    unsigned add(
       const char *F, int L, const cgc_contract_id<MTI> &idc,
       const mxArray* S, unsigned k
    );

   
   
 

    unordered_map <
       cgc_contract_id<MTI>,
       x3map<TQ,RTD>,
       MHash<MTI>
    > XBUF;

  protected:
  private:
};

   X3Map<gTQ> gCX;


template <class TQ>
class QMap {

  public:

    void init() {
       Q1.init(); Q2.init(); Q.init();
       I1.init(); I2.init(); II.init(); D.init(); qvec.init();
    };

    QMap<TQ>& init(const char *F, int L,
       const QVec &qv, const wbMatrix<TQ> &q1, const wbMatrix<TQ> &q2, 
       const wbMatrix< wbMatrix<TQ> > &QQ, 
       const wbMatrix< wbMatrix< const wbvector< CRef<TQ> >* > > *SS=NULL 
    );

    QMap<TQ>& init(const char *F, int L, const mxArray *a);

    template<class TD>
    int getCGZlist(
       const char *F, int L,
       wbMatrix<TQ> &qc,
       wbvector<TD> &dc,
       wbvector<INDEX_T> *IC=NULL,
       const wbvector<INDEX_T>* Ix=NULL,
       const qset<TQ> *q3=NULL,
       wbIndex *m3=NULL,
       const wbvector<double>* cc=NULL,
       double eps1=1E-10,
       double eps2=1E-14
    ) const;

    void Skip(const wbindex &I);

    template<class TD>
    void getIdentityQ(const char *F, int L,
       const wbvector<INDEX_T> &Sa, const wbvector<INDEX_T> &Sb,
       QSpace<TQ,TD> &A, char vflag=1
     ) const;

    bool isConsistent(
       const char *F=NULL, int L=0, const QVec *qv=NULL
     ) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;
    
    mxArray* toMx() const {
       mxArray *S=mxCreateStruct(1,1);
       add2MxStruct(S,0); return S;
    };

    void put(const char *vname, const char *ws="caller") const
    {  put(0,0,vname,ws); }

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);
       
       if (i) wblog(F_L,
          "ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F_L,
          "I/O putting variable '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };


    wbMatrix<TQ> Q1,Q2,Q;
    wbvector<INDEX_T> I1,I2,D;
    wbMatrix<INDEX_T> II;

    wbMatrix< const wbvector< CRef<TQ> >* > cg3;

    QVec qvec;

};

template <class TQ>
bool QMap<TQ>::isConsistent(const char *F, int L, const QVec *qv) const {

   unsigned d=qvec.Qlen();

   if (Q1.dim2!=d || Q2.dim2!=d) {
      if (!F) return 0; else wblog(F,L,
      "ERR QMap: invalid qset record length (%d,%d/%d)",
       Q1.dim2, Q1.dim2, d);
   }
   if (qv && (!qvec.sameType(*qv) || d!=(qv->Qlen()))) {
      if (!F) return 0; else wblog(F,L,
      "ERR QMap: Q type mismatch: %s; %s (%d/%d)",
       qv->toStr().data, qvec.toStr().data, qv->Qlen(), d);
   }
   return 1;
};


template <class T>
class pairPatternCG {
  public:

    pairPatternCG& init(unsigned i0=0, unsigned j0=0, unsigned n=0) {
       i=i0; j=j0; k.init(n); fac.init(n);
       return *this;
    };

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    unsigned i,j;
    wbvector<unsigned> k;
    wbvector<T> fac;

  protected:
  private:
};


template <class TQ, class TD>
class genRG_base{

  public:

     bool isEmpty() const {
        return (!q && !Sp.len && !Sz.len && !J.len && Z.isEmpty());
     };

     double normDiff(const char *F, int L,
        const genRG_base& B,
        TD eps=1E12
      ) const;

     void checkOrthoCommRel(
        const char *F, int L, const genRG_base &B) const;

     genRG_base& Sort();
     double SkipTiny(const char *F=0, int L=0);

     genRG_base& ApplyQFac(
        const wbvector<double> &qfac, const wbarray<double> *JM=NULL);

     genRG_base& operator=(const genRG_base& S) {
        q=S.q; J=S.J; Z=S.Z; Sp=S.Sp; Sz=S.Sz;
        return *this;
     };

     genRG_base& init() {
        q=0; J.init(); Z.init(); Sp.init(); Sz.init();
        return *this;
     };

     void compareStdSU2(const char *F=NULL, int L=0) const;

     genRG_base& save2(genRG_base& S) {
        S.q=q; J.save2(S.J); Z.save2(S.Z); q=0;
        Sp.save2(S.Sp); Sz.save2(S.Sz);
#ifdef CG_CHECK_MW_PERM
        P0.save2(S.P0);
#endif
        return S;
     };

     wbstring info() const {
        unsigned n=128; char s[n]; char sep[2]=" ";
        if (!q.isAbelian() && J.wbvector<TQ>::allIn(0,9)) sep[0]=0;
        unsigned l=snprintf(s,n,"%s [%s]",
            q.toStr().data, J.toStrf("",sep,q.qlen(),";").data);
        if (l>=n) wblog(FL,
           "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
        return s;
     };

     mxArray* toMx() const {
        mxArray *a=mxCreateStruct(1,1);
        add2MxStruct(a,0);
        return a;
     };

     mxArray* mxCreateStruct(unsigned m, unsigned n) const;
     void add2MxStruct(mxArray *S, unsigned i) const;

     QType q;

     wbvector< wbSparrayTD > Sp;
     wbvector< wbSparrayTD > Sz;

     qset<TQ> J;

     wbMatrix<double> Z;

#ifdef CG_CHECK_MW_PERM
     wbperm P0;
#endif

  protected:
  private:
};




template <class TQ, class TD>
class genRG_struct {
  public:

    genRG_struct() {};

    genRG_struct& checkInit(const char *F, int L, const QType &q0) {
       unsigned n=RSet.size();
       if (n>0) {
          if (!q0.type) wblog(FL,"ERR %s() got %d elements "
             "in buf for '%s'",FCT,n,q.toStr().data);
          return *this;
       }
       q=q0; return initBase(F,L);
    };

    genRG_struct& initBase(const char *F, int L, qset<TQ> *qs_=NULL) {
        q.validType(F_L);
        if (qs.isEmpty()) {
           switch (q.type) {
              case QTYPE_SUN: initBase_SUN(F_L,&qs); break;
              case QTYPE_SpN: initBase_SpN(F_L,&qs); break;
              default: wblog(FL,"ERR %s() "
             "type '%s' not implemented yet",FCT,q.toStr().data);
           }
        }
        if (qs_) (*qs_)=qs;
        return *this;
    };

    genRG_struct& initBase_SUN(const char *F, int L, qset<TQ> *qs=NULL);
    genRG_struct& initBase_SpN(const char *F, int L, qset<TQ> *qs=NULL);

    void initSU2(const qset<TQ> &J);

    void initCommRel(const char *F, int L,
       const genRG_base<TQ,TD> &R, char vflag=0);

    double checkCommRel(
      const char *F, int L, const genRG_base<TQ,TD> &R
    ) const;

    void printStatus(const char *F, int L, const char *istr=NULL);

    void getTensorProdReps(const qset<TQ> &q1, const qset<TQ> &q2){
       q.validType(FL);
       int n=getTensorProdReps_gen(q1,q2,'t');
       if (n>=0) return;

       wblog(FL,"ERR %s() "
         "failed to access RSets of input spaces (%d) !?",FCT,n);

    };

    int getTensorProdReps_gen(
       const qset<TQ> &q1, const qset<TQ> &q2, char tflag=0);

    unsigned getTensorProdReps(const qset<TQ> &q);
    unsigned getTensorProdReps(unsigned dmax=-1, char vflag=0);

    void put(const char *vname, const char *ws="caller") const
    {  put(0,0,vname,ws); }

    mxArray* toMx() const;
    void put(const char *F, int L,
      const char *vname, const char *ws="caller"
    ) const;

    QType q;

    wbvector< pairPatternCG<double> > CR;
    wbMatrix<double> DZ;

    qset<TQ> qs;

    typedef typename
       std::map < qset<TQ>, genRG_base<TQ,TD> >::const_iterator
       crMAP;

    map <qset<TQ>,
         genRG_base<TQ,TD>
    > RSet;

  protected:
  private:
};


template <class TQ, class TD>
class RStore {

 public:

   RStore(const char *s=NULL) {
      if (s) {
         if (!strcmp(s,"info")) {
            wbstring s=getName(typeid(TD));
            #ifdef __WB_MPFR_HH__
               wblog(FL,"<i> %s() uses type %s (=mpfr@%d) ",
               FCT,s.data, TD().prec()); 
            #else
               wblog(FL,"<i> %s() uses type %s",FCT,s.data); 
            #endif
         }
         else wblog(FL,"ERR %s() invalid usage (%s)",FCT,s);
      }
   };

   genRG_base<TQ,TD>* find_RSet(const QType &q, const TQ *qptr){
      genRG_struct<TQ,TD> &B=buf[q];
      if (B.RSet.size()) {
         qset<TQ> qs(q.qlen(),qptr,'r');
         auto it = B.RSet.find(qs);
         return (it!=B.RSet.end() ? &(it->second) : NULL);
      }
      else if (!B.q.isKnown()) { B.checkInit(FL,q); }
      return NULL;
   };

   genRG_struct<TQ,TD>& getRSet(
      const char *F, int L, const QType &q,
      qset<TQ> *qvec=NULL
   ){
      genRG_struct<TQ,TD> &B=buf[q];
         if (!B.RSet.size()) { B.q=q; B.initBase(F,L,qvec); } else
         if (qvec) { (*qvec)=B.qs; }
      return B;
   };

   void Info(const char *F=0, int L=0) const;

   genRG_base<TQ,TD>& getR(
      const char *F, int L, const QType &qtype, const TQ *q) const;

   genRG_base<TQ,TD>& getR(
      const char *F, int L, const CDATA_TQ &S,
      unsigned k
   ) const { return getR(F_L, S.t, S.qptr(k)); };

   wbMatrix<TQ> getZ(
      const char *F, int L, const QType &qtype, const TQ *q
   ) const {
      const genRG_base<TQ,TD> &R = getR(F_L,qtype,q);
      return R.Z;
   };

   unsigned qdim(
      const char *F, int L, const QType &qtype, const TQ *q
   ) const {

      if (qtype.isAbelian()) return 1;
      const genRG_base<TQ,TD> &R = getR(F_L,qtype,q);

      if (qtype.isSU2() && TQ(R.Z.dim1)!=q[0]+1) wblog(FL,"ERR %s() "
         "got multiplet dim mismatch (%d/%d)",FCT,R.Z.dim1,int(q[0]+1));

      return R.Z.dim1;
   };

   mxArray* toMx(const QType &q=QTYPE_UNKNOWN) const;

   unsigned LoadStore(const char *F, int L, const char *file);
   unsigned add(const QType &q, const mxArray* S, unsigned k);

   typedef typename 
      std::map < QType, genRG_struct<TQ,TD> >::const_iterator
      ciMAP;

   typedef typename
      std::map < qset<TQ>, genRG_base<TQ,TD> >::const_iterator
      crMAP;

   map <QType,
        genRG_struct<TQ,TD>
   > buf;

   protected:
   private:
};

   RStore<gTQ,RTD> gRG;"info");


class RCStore {

 public:

   RCStore() {
      const char *stag="RC_STORE", *sval=getenv(stag);
      if (!sval || !sval[0]) {
         root="(env RC_STORE not yet defined)";
         return;
      }

      root=sval;
      if (!Wb::fexist(root.data,'d')) wblog(FL,
         "ERR directory in %s does not yet exist\n%s",stag,root.data);

      if (!strncmp(root.data,"/home/",7))
      wblog(FL,"WRN %s() CG_STORE points to HOME",FCT);
   };


   template <class TQ, class TD>
   int save_CData(const char *F, int L, const CData<TQ,TD> &A);

   template <class TQ, class TD>
   int load_CData(
      const char *F, int L, const QSet<TQ> &Q, CData<TQ,TD> &A,
      char bflag=0);

   template <class TQ>
   int save_Std3(
      const char *F, int L, const QType &t, const qset<TQ> &J12) const;

   template <class TQ>
   int load_Std3(const char *F, int L,
       const QType &t, const qset<TQ> &J1, const qset<TQ> &J2) const;

   template <class TQ, class TD>
   int save_RSet(
      const char *F, int L, const genRG_base<TQ,TD> &R, char xflag=0);

   template <class TQ>
   int load_RSet(
      const char *F, int L, const QType &t, const qset<TQ> &J);

   int save_XMap(
      const char *F, int L, const cgc_contract_id<MTI> &idc) const;

   int load_XMap(
      const char *F, int L, const cgc_contract_id<MTI> &idc) const;

   template <class TQ>
   int get_directory(
      const char *F, int L, wbstring &file, const QSet<TQ> &Q,
      const char *ext=NULL, char mflag=0) const;

   template <class TQ, class TD>
   int get_directory(
      const char *F, int L, wbstring &file,
      const genRG_base<TQ,TD> &R, const char *ext=NULL, char mflag=0
    ) const;

   int get_directory(
      const char *F, int L, wbstring &file,
      const cgc_contract_id<MTI> &idc, const char *ext=NULL, char mflag=0
    ) const;

   wbstring root;

 protected:
 private:

   int get_directory(
      const char *F, int L, char *dir, unsigned n,
      const QType &t, const char *tag, const char *sub=NULL,
      char mflag=0) const;
};

   RCStore gStore;


template <class TQ, class TD>
double getSymmetryStates(const char *F, int L, const QType &q,
   const wbvector< wbSparrayTD > &Sp,
   const wbvector< wbSparrayTD > &Sz,
   wbvector< wbSparrayTD > &U,
   wbvector<unsigned> &dd,
   wbvector<genRG_base<TQ,TD> > &RR,
   wbMatrix<unsigned> *iOM=NULL,
   char vflag=1
);

template <class TQ>
SPIDX_T findMaxWeight(
   const QType &q, const wbMatrix<double> &Z, qset<TQ> *J=NULL, wbperm *P_=NULL
);


template <class TQ, class TD>
void get_SU2mat(TQ s2,
   wbvector<double> &sz,
   wbsparray<TD> &Sp, wbsparray<TD> &Sz,
   wbsparray<TD> &S2, wbsparray<TD> &E
);


#endif

