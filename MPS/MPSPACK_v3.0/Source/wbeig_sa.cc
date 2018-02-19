
/**********************************************************************/
/* stand-alone version of c_wbeigens.cc for MatLab customre support

# compiled to mex-type standalone

  gcc  -O2 -fopenmp -pthread -D_AMD_64 -o c_wbeigens_sa c_wbeigens_sa.cc \
     -I/usr/local/dist/DIR/matlab-R2011a/extern/include \
     -L/usr/local/dist/DIR/matlab-R2011a/bin/glnxa64 \
     -lstdc++ -lmx -lmex -lmat -lm -lmwlapack -lmwblas

# compiled to standalone using GotoBLAS2
# (shows similar CPU usage behavior as above)

  gcc  -Wall -o c_wbeigens_sa c_wbeigens_sa.cc \
     -L /home/a/Andreas.Weichselbaum/Source/GotoBLAS/  \
     -lstdc++ -lm -lgfortran -fopenmp -lpthread -lgoto2

# A. Weichselbaum ; Apr 03,2013
# *********************************************************************/


#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <ctype.h>
#include <climits>
#include <string.h>
#include <unistd.h>
#include <time.h>

#ifdef MATLAB_MEX_FILE
   #include <mex.h>
   #include <mat.h>
#endif

#define pINT ptrdiff_t

extern "C" {

#define dsyev dsyev_

    void dsyev (

     const char   &,
     const char   &,
     const pINT   &,
           double [],
     const pINT   &,
           double [],
           double [],
     const pINT   &,
           pINT   &
    );

#define ilaenv ilaenv_ 

   int ilaenv(
       const pINT   &,
       const char *,
       const char *,
       const pINT   &,
       const pINT   &,
       const pINT   &,
       const pINT   &
   );
 
};


template <class T>
class wbMatrix {

   public:

      wbMatrix() : data(NULL), dim1(0), dim2(0) {};

      wbMatrix(unsigned r, unsigned c=-1) : data(NULL), dim1(r), dim2(c){
         if (int(c)==-1) { dim2=c=dim1; }
         if (r && c) { new_data(r*c); }
      };

      ~wbMatrix() { if (data) { delete [] data; data=NULL; }};

      wbMatrix& init() { dim1=dim2=0;
         if (data) { delete [] data; data=NULL; }
         return *this;
      };

      wbMatrix& init(unsigned r, unsigned c=-1) {
         if (int(c)==-1) { c=dim1; }
         unsigned s=r*c; 
         if (s != dim1*dim2) { new_data(s); } else
         if (r && c) memset(data,0,s*sizeof(T));
         dim1=r; dim2=c;
         return *this;
      };

      unsigned numel() const { return dim1*dim2; };
      void info(const char *s) const;
      void print(const char *s=NULL) const;

      wbMatrix& operator= (const wbMatrix &M);

      wbMatrix& setRand(char sym=0);
      void eig(wbMatrix<T> &U, wbMatrix<T> &E);
       
#ifdef MATLAB_MEX_FILE
      mxArray* toMx() const;
#endif

      T *data;
      unsigned dim1, dim2;

   protected:
   private:

     void new_data(unsigned s, T* d=NULL) {
        if (data) { delete [] data; data=NULL; }
        if (s) { 
           data = new T[s];
           if (!data) {
              fprintf(stderr,"failed to allocate array [%d]",s);
              exit(-1);
           }
           if (d==NULL)
                memset(data,0,s*sizeof(T));
           else memcpy(data,d,s*sizeof(T));
        }
     }
};

template <class T>
wbMatrix<T>& wbMatrix<T>::operator= (const wbMatrix &M) {
    if (this!=&M) {
       unsigned s=M.numel();
       if (s!=dim1*dim2) 
            new_data(s,M.data);
       else memcpy(data,M.data,s*sizeof(T));

       dim1=M.dim1;
       dim2=M.dim2;
    }
    return *this;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::setRand(char sym) {
   static char first_call=1;
   double fac=2./(double)RAND_MAX;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   if (sym) {
      unsigned i,j;

      if (dim1!=dim2) { fprintf(stderr,
        "%s() symmetrization requires square matrix (%dx%d)",
       __FUNCTION__,dim1,dim2); exit(-1);
      }

      for (j=0; j<dim2; ++j) { data[j+j*dim1]=(T)(fac*rand()-1.);
      for (i=j+1; i<dim1; ++i) {
          data[i+j*dim1]=data[j+i*dim1]=(T)(fac*rand()-1.);
      }}
   }
   else {
      unsigned i,s=dim1*dim2;
      for (i=0; i<s; i++) data[i]=(T)(fac*rand()-1.);
   }

   return *this;
};

template <class T>
void wbMatrix<T>::info(const char *s) const {
   printf("  %-12s  %dx%d\n", s && s[0] ? s : "",dim1,dim2);
};

template <class T>
void wbMatrix<T>::print(const char *s) const {
   if (s && s[0]) printf("\n%s = [",s);

   unsigned i,j;
   for (i=0; i<dim1; ++i) { printf("\n");
   for (j=0; j<dim2; ++j) {
      printf(" %10.6g",data[i+dim1*j]);
   }}

   if (s && s[0]) printf("\n];\n");
   else printf("\n\n");
};

template <class T>
void wbMatrix<T>::eig(wbMatrix<T> &U, wbMatrix<T> &E) {

   wbMatrix<double> aux;
   unsigned n=dim1, lwork=3*n-1;
   pINT q=0;

   if (dim1!=dim2) { fprintf(stderr,
     "%s() symmetrization requires square matrix (%dx%d)",
    __FUNCTION__,dim1,dim2); exit(-1);
   }
   if (!dim1) { U.init(); E.init(); return; }

   if (n>255) {
      unsigned l=ilaenv(1,"dsytrd","U",(pINT)n,(pINT)n,(pINT)n,(pINT)n);
      printf("--> ispec=1: lwork=(%d+2)*%d  (minimal: 3*%d-1)\n",l,n,n);
      lwork=(l+2)*n;
   }

   U=*this; E.init(n,1);
   aux.init(lwork,1);

   dsyev (
      'V',
      'U',
       n,
       U.data,
       n,
       E.data,
       aux.data,
       aux.dim1,
       q
   );

   if (q) fprintf(stderr,"%s() DSYEV returned %d",__FUNCTION__,int(q));
};

#ifdef MATLAB_MEX_FILE

template <class T>
mxArray* wbMatrix<T>::toMx() const {

   unsigned i, n=dim1*dim2;
   mxArray *a=NULL;
   double *dr;

   mwSize s[2]; s[0]=dim1; s[1]=dim2;

   a=mxCreateNumericArray(2, s, mxDOUBLE_CLASS, mxREAL);
   if (a==NULL) {
      if (data!=NULL) fprintf(stderr,
         "\nERR %s() invalid usage (%lX, %lX; %d)",
       __FUNCTION__,(long)data,(long)a,n);
      return a;
   }
   if (!n) return a;

   dr=mxGetPr(a);
   if (dr==NULL) fprintf(stderr,
      "ERR failed to allocate mxArray (%dx%d)])",(int)s[0],(int)s[1]);
   if (data==NULL && n) fprintf(stderr,"ERR invalid input array");

   for (i=0; i<n; ++i) dr[i]=double(data[i]);
   return a;
};

#endif


void doflush() {
   fflush(0);
#ifdef MATLAB_MEX_FILE
   mexEvalString("pause(0);");
#endif
};


#ifdef MATLAB_MEX_FILE
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]){
#else
int main (int nargin, char **argin) {
#endif

    unsigned D=2046, nrepeat=4;
    unsigned j;

#ifdef MATLAB_MEX_FILE
    if (nargin==1) {
       const mxArray *a=argin[0];
       unsigned Din=0;
       if (!mxIsDouble(a) || mxIsSparse(a) || mxGetNumberOfElements(a)!=1) {
          fprintf(stderr,"\nERR invalid input (dimension) !??");
          exit(-1);
       }
       Din=(unsigned)mxGetPr(a)[0];
       if (int(Din)<=0 || Din>9999) {
          fprintf(stderr,"\nERR invalid input: dimension D=%d !?",Din);
          exit(-1);
       }
       D=Din;
    }
#else
    if (nargin==2) {
       unsigned Din=0;
       if (sscanf(argin[1],"%d",&Din)<=0 || Din>9999) {
          fprintf(stderr,"\nERR invalid input: dimension D=%d !??",Din);
          exit(-1);
       }
       D=Din;
    }
#endif

    wbMatrix<double> H(D,D), U;
    wbMatrix<double> E;

    printf("\n  diagonalization of %d matrices of dimension %d:\n\n",nrepeat,D);
    doflush();

    for (j=0; j<nrepeat; ++j) { H.setRand('s');
       printf("\r%4d/%d ...\r",j+1,nrepeat); doflush();
       H.eig(U,E);

#ifdef MATLAB_MEX_FILE
       if (nargout) {
          if (nargout==1) { argout[0]=E.toMx(); } else
          if (nargout>=2) {
             argout[0]=U.toMx(); argout[1]=E.toMx();
             if (nargout==3) argout[2]=H.toMx(); else
             if (nargout> 3) fprintf(stderr,"\nERR invalid usage\n");
          }
          break;
       }
#endif
       if (D==7) {
          H.print("H"); U.print("U"); E.print("E");
       }
    }
    printf("\r%40s\r","");
};


