#ifndef __WB_WB_LIB_MEX_HCC__
#define __WB_WB_LIB_MEX_HCC__

/* ---------------------------------------------------------------- *
 * library for mex programs
 * Wb (C) Jan 2006
 * ---------------------------------------------------------------- */

#include <time.h>
#include <signal.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>

#include <cmath>
#include <ctype.h>
#include <climits>
#include <stdarg.h>
#include <string.h>

#include <string>
#include <iostream>
#include <typeinfo>

#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <map>
#include <unordered_map>
#include <algorithm>
#include <omp.h>



unsigned ARG_CHECK=1;

#define STRLEN  1023
#define MAX_STRLEN  STRLEN
char str[STRLEN+1];

template <class T> class wbvector;
template <class T> class wbMatrix;
template <class T> class wbarray;
template <class TD> class wbsparray;
template <class TQ, class TD> class QSpace;

class wbcomplex;
class wbstring;
class bitset;
class wbperm;
class wbindex;
class wbCMat;
class DATA;

class IndexLabels;
class ctrIdx;

using namespace std;


#if 1
   #define  INDEX_T size_t
   #define sINDEX_T long
   #define   PERM_T size_t
   #define  sPERM_T long
#else
   #define  INDEX_T unsigned
   #define sINDEX_T int
   #define   PERM_T unsigned
   #define  sPERM_T int
#endif

#if 1
   #define  SPIDX_T  size_t
   #define sSPIDX_T  long
#else
   #define  SPIDX_T  unsigned
   #define sSPIDX_T  int
#endif


#define WBINDEX  wbvector<INDEX_T>
#define WBPERM   wbvector<PERM_T>
#define WBIDXMAT wbMatrix<INDEX_T>


#define   UVEC        wbvector<unsigned>
#define C_UVEC  const wbvector<unsigned>
#define   DMAT        wbMatrix<double>
#define C_DMAT  const wbMatrix<double>
#define   UMAT        wbMatrix<unsigned>
#define C_UMAT  const wbMatrix<unsigned>
#define   IMAT        wbMatrix<int>
#define C_IMAT  const wbMatrix<int>
#define C_TMAT  const wbMatrix<T>
#define C_MEX   const mxArray
#define C_UINT  const unsigned
#define C_CHAR  const char

#define DMAT0   wbMatrix<double>()


#define __FL__ (file ? file : __FILE__), (line ? line : __LINE__)
#define F_L   (F ? F : __FILE__), (F ? L : __LINE__)
#define F_LF  (F ? F : __FILE__), (F ? L : __LINE__), __FUNCTION__
#define F_L_F (F ? F : __FILE__), (F ? L : __LINE__), (fct ? fct : __FUNCTION__)

#define FSTR  __FILE__ ": "
#define FLINE __FILE__, __LINE__
#define FL    __FILE__, __LINE__
#define FLF   __FILE__, __LINE__, __FUNCTION__
#define FCT   __FUNCTION__
#define FCL   __FUNCTION__, __LINE__


inline const char* shortPF(const char *P, const char *F, int L, char *s=NULL);

#ifndef PFL
#ifdef PROG
   #define PFL  shortPF(PROG,__FILE__,__LINE__), __LINE__
   #define PF_L shortPF(PROG,(F?F:__FILE__),(F?L:__LINE__)), (F?L:__LINE__)
#elif defined MATLAB_MEX_FILE
   #define PFL  shortPF(mexFunctionName(),__FILE__,__LINE__), __LINE__
   #define PF_L shortPF(mexFunctionName(),(F?F:__FILE__),(F?L:__LINE__)), (F?L:__LINE__)
#else
   #define PFL  shortPF("",__FILE__,__LINE__), __LINE__
   #define PF_L shortPF("",(F?F:__FILE__),(F?L:__LINE__)), (F?L:__LINE__)
#endif
#endif

#ifndef Inf
#define Inf FP_INFINITE
#endif

#if !defined(MAX)

   template <class T> inline const T& MAX(const T&a, const T&b) {
   if (a>b) return a; else return b; }

   template <class T> inline const T& MAX(const T&a, const T&b, const T&c) {
      if (a>b)
           { return (a>c ? a : c); }
      else { return (b>c ? b : c); }
   }

#endif
#if !defined(MIN)

   template <class T> inline
   const T& MIN(const T&a, const T&b) {
      if (a<b) return a; else return b;
   };

   template <class T> inline
   const T& MIN(const T&a, const T&b, const T&c) {
      if (a<b)
           { return (a<c ? a : c); }
      else { return (b<c ? b : c); }
   };

#endif
#if !defined(SWAP)

   template <class T>
   inline void SWAP(T &a, T &b) { T x=a; a=b; b=x; }

#endif

#if !defined(SSTR)
#define SSTR(a) a.sizeStr().data
#endif

#if !defined(SSTR_)
#define SSTR_(a) a->sizeStr().data
#endif

#if !defined(TOSTR)
#define TOSTR(a) a.toStr().data
#endif


#if !defined(POW)

   template <class T>
   inline T POW(const T &a, const unsigned &n) {
      T x = n ? a : 1;
      for (unsigned i=1; i<n; i++) x*=n;
      return x;
   }

#endif
#if !defined(SGN)

   template <class T>
   inline char SGN(const T &a) {
      if (a<0) return -1; else
      if (a>0) return +1;
      return 0;
   }

#endif
#if !defined(NUMCMP)

   template <class T>
   inline char NUMCMP(const T &a, const T &b) {
      if (a<b) return -1; else
      if (a>b) return +1;
      return 0;
   }

   template <class Ta, class Tb>
   inline char NUMCMP(const Ta &a, const Tb &b) {
      if (a<b) return -1; else
      if (a>b) return +1;
      return 0;
   }

#endif




#ifdef MAC
   #include "Include/nomex.hh"
   #include <Accelerate/Accelerate.h>
#elif defined(NOMEX)
   #include "nomex.hh"
#else


   #include <mex.h>
   #include <mat.h>



#endif

#ifdef MAIN
   #undef printf
#endif

#ifndef PROG
   #define PROG "(?program?)"
#endif


void dbstop(const char* file, int line);
void ExitMsg(const char* s="");
void doflush();

#if 0

#ifndef MWSIZE_MAX
   typedef int mwSize;
   typedef int mwIndex;
#else
   char* mwSize_check() {
      if (sizeof(mwSize)!=sizeof(int) || sizeof(mwIndex)!=sizeof(int)) {
         printf(
           "\n   WRN MatLab size data type `mwSize' does not match int!"
           "\n   WRN (%d,%d; %d)\n\n",(int)sizeof(mwSize),
         (int)sizeof(mwIndex), (int)sizeof(int) ); fflush(0);
      }
      return 0;
   };
   char *str_mws=mwSize_check();
#endif
#endif


#ifdef LOAD_CGC_QSPACE
   #define USE_WB_MPFR
#endif


#ifdef  __WBDEBUG__
#define       __WB_MEM_CHECK__
#undef NSAFEGUARDS
#elif defined __WB_MEM_CHECK__
#define __WBDEBUG__
#endif

#include "wblibx.h"

#include "wblog.h"
#include "wblog.c"

#include "memtrack.hh"

#include "util.hh"
#include "memlib.hh"
#include "mxlib.hh"

#ifdef USE_WB_MPFR
#include <gmp.h>
#include <mpfr.h>
#include "wbmpfr.hh"
#endif

#include "wbcomplex.hh"
#include "gathered.hh"
#include "wbsort.hh"

#include "wbvector.hh"
#include "wbstring.hh"
#include "wbbitset.hh"
#include "wbperm.hh"
#include "wbindex.hh"

#include "wbclock.hh"
#include "wbMatrix.hh"

#include "wbcmat.hh"

#include "mexlib.hh"

#include "wbMatIO.hh"

#include "wbarray.hh"
#include "wbsparray.hh"
#include "wbio.h"
#include "wbCMat.hh"

#include "wbblas.hh"
#include "wbMatrix_blas.hh"
#include "wbarray_blas.hh"

#include "wbopts.hh"
#include "wbidat.hh"
#include "data.hh"

#ifdef LOAD_CGC_QSPACE
#include "clebsch.hh"
#include "QSpace_aux.hh"
#include "QSpace.hh"
#endif


#ifdef USE_WB_MPFR
#include "wbmpfr.cc"
#endif

#ifdef __WB_MEM_CHECK__
#include "memtrack.cc"
#endif

#include "wbsort.cc"
#include "wbstring.cc"
#include "wbbitset.cc"
#include "memlib.cc"
#include "wbcomplex.cc"
#include "gathered.cc"
#include "util.cc"

#ifdef LOAD_CGC_QSPACE
#include "QSpace_aux.cc"
#include "QSpace.cc"
#endif

#include "wbvector.cc"
#include "wbperm.cc"
#include "wbindex.cc"
#include "wbMatrix.cc"
#include "wbarray.cc"
#include "wbsparray.cc"

#ifdef LOAD_CGC_QSPACE
#include "clebsch.cc"
#include "clebsch_aux.cc"
#include "clebsch_old.cc"
#endif

#include "wbclock.cc"
#include "mexlib.cc"
#include "wbCMat.cc"
#include "wbMatrix_blas.cc"
#include "wbarray_blas.cc"

#include "wbio.c"

#include "wbmsg.hh"


   template <class TD, class T>
   inline void DSET(TD &a, T* r, T* i, unsigned k) {
      if (i) if (i[k]!=0) wblog(FL,
      "ERR Missing out on imaginary data (%g) !?", double(i[k]));
      a=(TD)r[k];
   };


   template <class T>
   inline void DSET(wbcomplex &a, T* r, T* i, unsigned k) {
      a.set((double)r[k], i ? (double)i[k] : 0.);
   };


   class wb_init__ {
     public:
       wb_init__() : x(0) {
          srand((unsigned int)time((time_t*)NULL));

         #ifdef LOAD_CGC_QSPACE
          { unsigned i,k;
            if ((i=Wb::GetEnv(0,0,"CG_VERBOSE",k))==0) {
               if (k>9) wblog(FL,
                  "WRN %s() invalid env CG_VERBOSE=%d (0..9)",FCT,k);
               else {
                  CG_VERBOSE=k;
               }
            }
            else if (int(i)!=-1) wblog(FL,
               "WRN %s() invalid env CG_VERBOSE (e=%d)",FCT,i
            );
          }
         #endif
       };

     private:
       char x;
   };

   wb_init__ gwb_init;


#ifdef MAIN

#ifdef MATLAB_MEX_FILE
#warning "              \
Don't define MAIN with MEX files! \
this kills the whole matlab seesion!!        \
           see ExitMsg() in wblib.h"
#endif

   void ExitMsg(const char* s) {
      printf("%s ??? %s ???\n",shortFL(FL), s ? s:"");
   #ifdef DBSTOP
        dbstop(FL);
   #else
        exit(-1);
   #endif
   }

   void doflush() { fflush(0); };

#else

   void ExitMsg(const char* s) {
   #ifdef DBSTOP
        if (s) printf("%s >%s<\n",shortFL(FL),s ? s:"");
        dbstop(FL);
   #else
        printf("\n"); mexErrMsgTxt(s);
   #endif
      printf("\n%s dbstop - and on we go ...\n\n", shortFL(FL));
      mexErrMsgTxt(s);
   }

   void doflush() {
       fflush(0);
       if (!omp_in_parallel()) {
          mexEvalString("pause(0);");
       }
   };

#endif



#endif

