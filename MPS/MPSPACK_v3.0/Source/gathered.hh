#ifndef __WB_GATHERED_HH__
#define __WB_GATHERED_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

namespace Wb {
   inline const char* basename(const char *data, char c='/');

   template<class T>
   int GetEnv (const char *F, int L, const char *name, T& val);

   wbstring repHome(const char *file);

   inline int INT(const double a) {
      int i=::round(a);
      if (std::fabs(i-a)>1E-10) wblog(FL,"ERR %s(%g) !??",FCT,a);
      return i;
   }

   const char* strstri (const char *s1, const char *s2);
   const char* strstrw (const char *s1, const char *s2);
   const char* strstrwi(const char *s1, const char *s2);
   const char* strchri(const char *s, char c);

   bool isword(const char &c) { return (
        (c>='0' && c<='9') ||
        (c>='a' && c<='z') ||
        (c>='A' && c<='Z') || c=='_'
   ); };

   unsigned strnlen(const char *s, unsigned n) {
      for (unsigned i=0; i<n; ++i) if (!s[i]) return i;
      return n;
   };

   template <class T>
   bool isLower(const T* a, const T* b, const size_t n);

   template <class T>
   bool isEqual(const T* a, const T* b, const size_t n);

   template <class T>
   bool allEqual(const T* a, const size_t n, const T &x);

   template <class T>
   bool anyEqual(const T* a, const size_t n, const T &x);

   template <class T>
   bool anyUnqual(const T* a, const size_t n, const T &x);


   template <class T>
   T getdscale(const T* d, size_t n);

   template <class T> inline
   void scale_eps(T &eps, const T* d, size_t n);

}


template <class T, class T2> inline
char GOT_EPS(const T &a, const T2 &eps) {
   if (a==0 || eps!=0 && ABS(a)<eps) return 1;
   return 0;
}


void initSigHandler(const char i='i');
void WBSigHandler(int i);

char *gsh_F; int gsh_L;


class wbSigHandler {

  public:
    wbSigHandler(const char *F, int L) {
       gsh_F = new char[strlen(F)+1]; strcpy(gsh_F, F); gsh_L=L;
#ifndef DBSTOP
       initSigHandler('i');
#endif
    };

   ~wbSigHandler() {
#ifndef DBSTOP
       initSigHandler('r');
#endif
       delete [] gsh_F; gsh_F=NULL; gsh_L=0;
    };

    void call99() { doflush(); 
#ifndef DBSTOP
       WBSigHandler(99);
#endif
    };

  protected:
  private:

};


template <class TX>
class wbREF {

  public:

    wbREF (const TX &A, char flag) : orig(&A), ref(0) {
       if (flag) ref=new TX(A); else ref=&A;
    };

    TX& Ref() { return *ref; }

   ~wbREF() { if (ref!=orig) delete ref; }
   
  protected:
  private:
    const TX *orig, *ref;
};


template <class T>
class num2Fmt {

  public:
   num2Fmt(const char *F, int L,
      char *s, int m0=-9, int p0=-9, char t0=-1) : m(m0), p(p0) {
      t[0]=t0; t[1]=t[2]=0;
      get_fmt(F_L,s);
   };

   int check_paras();
   void get_fmt(const char *F, int L, char *s);

  private:

   char m,p;
   char t[3];
};

template <class T>
char* defaultFmt(char *fmt, const T &x, int n=-1, int p=-1, char t=0);


   template<class T> inline
   wbstring num2Str(const T &x, const char *fmt="");

   template<class T> inline
   int num2int(const char *F, int L, const T &x);

   template<class T> inline
   unsigned checkInt(
      const char *F, int L, const T* x, size_t n,
      double eps=1E-14
   );

   char* strpad(char *s, char p, unsigned n, unsigned w=1) {
      unsigned i,l=strlen(s);
      for (i=0; i<w; i++) s[l++]=' ';
      for (i=l; i<n; i++) s[l++]=p;
      return s;
   };

   template <class T>
   int charGetNumber(const char *F, int L, const char *s, T &x);

   wbstring getName(const std::type_info &type_id);
   bool isBaseType(const std::type_info &type_id);
   bool isComplexType(const std::type_info &type_id);

   wbstring wbTimeStamp();

   void inl(unsigned n) { for (unsigned i=0; i<n; i++) printf("\n"); }

   void quietErrLog(const char *F, int L, const char *msg, char qflag) {
       if (qflag) snprintf(str,256,"%s [%s:%d]",msg,basename(F),L);
       else wblog(F,L,"ERR %s",msg);
   }

   template <class T1, class T2> inline
   void safeConvert(const char *F, int L, const T1 &x1, T2 &x2);


   template <class T1, class T2> inline
   int checkUpdate(T1& a, const T2 &r) {
       if (a!=T1(r)) { a=T1(r); return 1; } else return 0;
   };


   template <class T>
   bool uniformRange(const T* d, const size_t n);

   template <class T>
   T maxRange(const T* d, const size_t n);

   template <class T>
   T minRange(const T* d, const size_t n);

   template <class T>
   T addRange2(const T* d, const size_t n);

   template <class T>
   T addRange(const T* d, const size_t n);

   template <class T>
   void addRange(
      const T* a, const T* b, T* c, const size_t n,
      size_t stride=1
   );

   template <class T>
   void addRange(const T* a, T* c, size_t n);

   template <class T>
   void diffRange(
      const T* a, const T* b, T* c, size_t n,
      size_t stride=1
   );

   template <class T>
   void setRange2avg(T* a, size_t m);

   template <class T>
   size_t nnzRange(const T* a, size_t m);

namespace Wb {

   template <class T, class T2> inline
   void cpyStride(
      T2* d,
      const T* d0,
      unsigned n,
      const unsigned *idx, unsigned m,
      unsigned D2=-1,
      unsigned D0=-1,
      char add_flag=0
   );

   template <class T, class T2> inline
   void cpyStride(
      T2* d,
      const T* d0,
      unsigned n,
      unsigned m,
      unsigned D2=-1,
      unsigned D0=-1,
      char add_flag=0
   ){
      return cpyStride(d,d0,n,(unsigned*)NULL,m,D2,D0,add_flag);
   };
};

   template <class T, class T2> inline
   void cpyRange(const T* a, T2* c, const size_t n, char check_type=1);

   template <class T> inline
   void cpyRangeR(const wbcomplex* a, T* c, const size_t n);
   template <class T> inline
   void cpyRangeI(const wbcomplex* a, T* c, const size_t n);
   template <class T> inline
   void cpyRangeA(const wbcomplex* a, T* c, const size_t n);

   template <class T>
   T prodRange(const T* d, const size_t n);
   template <class T>
   T prodRange(const T* d, const size_t n, T x0);

   template <class T>
   T prodRange(const T* d1, const T* d2, size_t n);
   template <class T>
   T prodRange(const T* d1, const T* d2, size_t n, T x);

   template <class T>
   void timesRange(
      T* d, T x, const size_t n, size_t stride=1
   );

   template <class T>
   T rangeNorm2(const T* d, size_t n, size_t *k=NULL);

   template <class T>
   T rangeNormDiff2(const T* d1, const T* d2, size_t n,
      const T &fac, size_t *k=NULL
   );

   template <class T>
   T rangeNormDiff2(const T* d1, const T* d2, size_t n);

   template <class T>
   T rangeMaxDiff(const T* d1, const T* d2, size_t n, size_t *k=NULL);

   template<class T> inline
   T overlap(
      const T *v1, const T *v2, size_t n,
      size_t stride=1, char tnorm=0
   );

   template<class T> 
   T gs_project_range(
      T *v1, const T *v2,
      size_t n, size_t stride=1, char isnorm=0, char tnorm=0
   );

namespace Wb {

   template<class T>
   void chopTiny_float(T *d, size_t n, T dref=-1);

   template <class T> inline
   int recCompare(
      const T* a, const T* b, const size_t n, char lex=1
   ){
      if (lex>0) {
         for (size_t i=0; i<n; ++i) {
            if (a[i]<b[i]) return -1; else
            if (a[i]>b[i]) return +1;
         }
      }
      else {
         for (size_t i=n-1; i<n; --i) {
            if (a[i]<b[i]) return -1; else
            if (a[i]>b[i]) return +1;
         }
      }
      return 0;
   };

   template <class T> inline
   int recCompare(
      const T* a, const T* b, const size_t n, char lex,
      T eps, int *nwrn_=NULL
   ){
      if (eps<=0) return recCompare(a,b,n,lex);

      int nwrn0=0, *nwrn=(nwrn_ ? nwrn_ : &nwrn0);
      T x=0.2*eps;

      if (lex>0) {
         for (size_t i=0; i<n; ++i) {
            if (a[i]<b[i]-eps) return -1; else
            if (a[i]>b[i]+eps) return +1; else
            if (!nwrn[0] && std::fabs(a[i]-b[i])>x) {
               wblog(FL,"WRN %s() got |a-b|/eps=%.3g (%d/%d)",
               FCT,(a[i]-b[i])/eps,i+1,n); ++nwrn[0];
            }
         }
      }
      else {
         for (size_t i=n-1; i<n; --i) {
            if (a[i]<b[i]-eps) return -1; else
            if (a[i]>b[i]+eps) return +1; else
            if (!nwrn[0] && std::fabs(a[i]-b[i])>x) {
               wblog(FL,"WRN %s() got |a-b|/eps=%.3g (%d/%d)",
               FCT,(a[i]-b[i])/eps,i+1,n); ++nwrn[0];
            }
         }
      }
      return 0;
   };

   template <class T>
   int cmpRange(
      const T *a, const T *b, const size_t n, char lex=+1
   ){
      if (lex>0) {
         for (size_t i=0; i<n; ++i) { if (a[i]!=b[i]) {
            return (a[i]<b[i] ? -1 : +1);
         }}
         return 0;
      }
      else {
         if (!lex) wblog(FL,"WRN %s() got lex=%d",FCT,lex);
         for (size_t i=n-1; i<n; --i) { if (a[i]!=b[i]) {
            return (a[i]<b[i] ? -1 : +1);
         }}
         return 0;
      }
   };

}

   template <class T> void cpyZRange(
   const double *R, const double *I, T *Z, const size_t n);

   template <class T> void splitZRange(
   const T *Z, double *R, double *I, const size_t n);


   void invertIndex(
      const wbvector<unsigned> &i1, unsigned N,
      wbvector<unsigned> &i2, char lflag=1
   );

   void setRand(wbvector<wbcomplex> &zz, double fac=1., double shift=0.);

   template<class T>
   void set2avg(T *dd, const size_t n);

   template<class T>
   void getDiff(const wbvector<T> &xx, wbvector<T> &dx);

   template<class T>
   T wbIntTrapez(const T *xx, const T *yy, const size_t n);

   inline bool contains(const char* s, const char &x, const unsigned &m=16) {
      if (s)
      for (unsigned i=0; s[i] && i<m; i++) { if (s[i]==x) return 1; };

      return 0;
   };


   template<class T>
   void markSet(
      wbvector< wbvector<T>* > &E0,
      wbvector< wbvector<char> > &mark,
      int &Nkeep,
      double Etrunc,
      wbvector<T> &E,
      double eps=0.,
      double b=0.,
      int dmax=-1,
      const char *dir=0"asc" (keep lowest)
   );

   template<class T>
   void markSet(const char *F, int L,
      const wbvector<T> &E,
      wbvector<char> &mark,
      int &Nkeep,
      double Etrunc=0,
      double b=0.,
      int dmax=-1,
      const char *dir=0"asc" (lowest Nkeep), "desc" (largest Nkeep)
   );


#endif

