#ifndef __WB_UTIL_HH__
#define __WB_UTIL_HH__

/* -------------------------------------------------------------------- */
namespace Wb {
/* -------------------------------------------------------------------- */

   int Rational(
      double &x,          // rational number (returns approximation of any)
      long &P, long &Q,   // rational representation P/Q (output)
      double *rerr=NULL,  // relative error
      unsigned *nmax=NULL,// max number of iterations actually used
      wbvector<double>* aa=NULL, // continued fraction sequence
      unsigned niter=6,   // (max) number of iterations
      long pqmax=-1,      // max enumerator/denominator (<0: no limit)
      double eps=1E-6,    // threshold for acceptance within continued fraction
      double eps2=1E-12,  // overall threshold for acceptance
      char vflag=0
   );

   template<class T>
   double FixRational(
      const char *F, int L, T *d, unsigned n, char rflag=0,
      unsigned niter=6, 
      long pqmax=-1,
      double eps=1E-6, double eps2=1E-12,
      double *rz=NULL, double *ra=NULL, char vflag=0
   );

   template<class T>
   double SkipZeros(T *d, unsigned n, double eps=1E-14){ return 0; };

   void CleanUp(const char *F, int L, const char *istr=NULL);

   unsigned long getSize(const char *f, const char *tag);

   wbstring size2Str(double n);
   const char* filename(const char *s);

   size_t getFilesize(const char* f);

   bool fexist(const char *f, char type='*');

   template <class TI>
   TI sub2ind(const TI *s, const TI *I, unsigned n);

   template <class T1, class T2> inline 
   void ind2sub(T1 k, const T1 *s, T2 *I, unsigned n);

   int atoi(const char *s, int &iout, unsigned n, char white=0);

   template <class T>
   size_t findfirst_sorted(const char *F, int L,
      const T* d0, const T* dd, size_t m, size_t N, size_t M,
      char lex=1"row-major")
   );

   template <class T>
   size_t findlast_sorted(const char *F, int L,
      const T* d0, const T* dd, size_t m, size_t N, size_t M,
      char lex=1"row-major")
   );


   template <class T0, class T> inline
   size_t findfirst_sorted(const char *F, int L,
      const T0* d0, const T* dd, size_t m, size_t N, size_t M, char lex=1
   ){ 
      T d2[m];
      for (unsigned i=0; i<m; ++i) d2[i]=T(d0[i]);
      return findfirst_sorted(F,L,d2,dd,m,N,M,lex);
   };

   template <class T0, class T> inline
   size_t findlast_sorted(const char *F, int L,
      const T0* d0, const T* dd, size_t m, size_t N, size_t M, char lex=1
   ){
      T d2[m];
      for (unsigned i=0; i<m; ++i) d2[i]=T(d0[i]);
      return findlast_sorted(F,L,d2,dd,m,N,M,lex);
   };

};


template <class T>
class WbUtil {
  public:


     static void adjust_tnorm(char &tnorm);

     static bool hasConj();
     static bool hasNoConj();

     static bool isPOD();
     static bool isFloat();
     static bool isInt();

     static bool isComplex();
     static bool isReal();

  protected:
  private:
};


template <class T> inline
bool WbUtil<T>::isPOD(){ return 0; };


template <> inline bool WbUtil<char>::isPOD(){ return 1; };
template <> inline bool WbUtil<unsigned char>::isPOD(){ return 1; };

template <> inline bool WbUtil<int>::isPOD(){ return 1; };
template <> inline bool WbUtil<unsigned>::isPOD(){ return 1; };

template <> inline bool WbUtil<short>::isPOD(){ return 1; };
template <> inline bool WbUtil<unsigned short>::isPOD(){ return 1; };

template <> inline bool WbUtil<long>::isPOD(){ return 1; };
template <> inline bool WbUtil<unsigned long>::isPOD(){ return 1; };

template <> inline bool WbUtil<float>::isPOD(){ return 1; };
template <> inline bool WbUtil<double>::isPOD(){ return 1; }; 
template <> inline bool WbUtil<long double>::isPOD(){ return 1; }; 
template <> inline bool WbUtil<wbcomplex>::isPOD(){ return 1; };


template <class T> inline
void WbUtil<T>::adjust_tnorm(char &tnorm){ tnorm=0; };
template <> inline
void WbUtil<wbcomplex>::adjust_tnorm(char &tnorm __attribute__ ((unused))){};

template <class T> inline bool WbUtil<T>::hasConj(){ return 0; };
template <> inline bool WbUtil<wbcomplex>::hasConj(){ return 1; };

template <class T> inline bool WbUtil<T>::hasNoConj(){ return 1; };
template <> inline bool WbUtil<wbcomplex>::hasNoConj(){ return 0; };


template <class T> inline
bool WbUtil<T>::isComplex(){ return 0; };


template <> inline bool WbUtil<wbcomplex>::isComplex(){ return 1; };


template <class T> inline
bool WbUtil<T>::isReal(){ return 1; };


template <> inline bool WbUtil<wbcomplex>::isReal(){ return 0; };


template <class T> inline
bool WbUtil<T>::isFloat(){ return 0; };


template <> inline bool WbUtil<double>::isFloat(){ return 1; }; 
template <> inline bool WbUtil<long double>::isFloat(){ return 1; }; 
template <> inline bool WbUtil<float>::isFloat(){ return 1; };
template <> inline bool WbUtil<wbcomplex>::isFloat(){ return 1; };


template <class T> inline
bool WbUtil<T>::isInt(){ return 0;}


template <> inline bool WbUtil<int>::isInt(){ return 1; }
template <> inline bool WbUtil<char>::isInt(){ return 1; }
template <> inline bool WbUtil<long>::isInt(){ return 1; }

template <> inline bool WbUtil<unsigned>::isInt(){ return 1; }
template <> inline bool WbUtil<unsigned char>::isInt(){ return 1; }
template <> inline bool WbUtil<unsigned long>::isInt(){ return 1; }


#include <sys/time.h>
#include <sys/resource.h>

class cpu_time {

  public:

    cpu_time() { tcpu=time('c'); tsys=time('s'); };
    void set() { tcpu=time('c'); tsys=time('s'); };

    cpu_time since() {
       cpu_time dt; dt.tcpu-=tcpu; dt.tsys-=tsys;
       return dt;
    };

    static double time(char flag) {
       struct rusage result;
       getrusage(RUSAGE_SELF,&result);
       if (flag=='c') {
          return (
             double(result.ru_utime.tv_sec)
           + 1e-6 * result.ru_utime.tv_usec
          );
       } else {
          return (
             double(result.ru_stime.tv_sec)
           + 1e-6 * result.ru_stime.tv_usec
          );
       }
    };

    wbstring toStr(char flag='c');
    double tcpu, tsys;

  protected:
  private:
};

  cpu_time gCPUTime;


class wbsys {

  public:


    wbstring MemTot2Str() const;
    unsigned long getMemTot() const {
       return Wb::getSize("/proc/meminfo","MemTotal"); };

    wbstring MemFree2Str() const;
    unsigned long getMemFree() const {
       return Wb::getSize("/proc/meminfo","MemFree"); };

    wbstring SwapTot2Str() const;
    unsigned long getSwapTot() const {
       return Wb::getSize("/proc/meminfo","SwapTotal"); };

    wbstring SwapFree2Str() const;
    unsigned long getSwapFree() const {
       return Wb::getSize("/proc/meminfo","SwapFree"); };

    char checkSwapSpace(const char *F=NULL, int L=0) const;

  protected:
  private:

};


class wbtop {

  public:

    wbtop() { pid=getpid(); ppid=getppid(); };

    wbstring VmSize2Str() const;
    unsigned long getVmSize() const { return getSize_proc("VmSize"); };

    wbstring VmPeak2Str() const;
    unsigned long getVmPeak() const { return getSize_proc("VmPeak"); };

    char runningLarge(const char *F=NULL, int L=0,
       double th0=0.80,
       double fac=1.2
    );

    pid_t pid, ppid;

  protected:
  private:

    unsigned long getSize_proc(const char *tag) const {
       const size_t n=32; char f[n];
       int i=snprintf(f,n,"/proc/%d/status",pid);
       if (i<0 || unsigned(i)>=n) {
          wblog(FL,"ERR %s() string too short (%s)",FCT,f);
          doflush();
       }
       return Wb::getSize(f,tag);
    };
};



#endif

