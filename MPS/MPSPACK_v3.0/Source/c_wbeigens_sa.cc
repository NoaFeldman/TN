
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
   unsigned n=dim1;
   pINT q=0;

   if (dim1!=dim2) { fprintf(stderr,
     "%s() symmetrization requires square matrix (%dx%d)",
    __FUNCTION__,dim1,dim2); exit(-1);
   }
   if (!dim1) { U.init(); E.init(); return; }

   U=*this; E.init(n,1);
   aux.init(3*n-1,1);

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


int main (int nargin, char **argin) {

    unsigned D=2046, nrepeat=4;
    unsigned j;

    if (nargin==2) {
       unsigned Din=0;
       if (sscanf(argin[1],"%d",&Din)<=0 || Din>9999) {
          fprintf(stderr,"\nERR invalid input: dimension D=%d !?",Din);
          exit(-1);
       }
       D=Din;
    }

    wbMatrix<double> H(D,D), U;
    wbMatrix<double> E;

    printf("\n  diagonalization of %d matrices of dimension %d:\n\n",nrepeat,D);

    for (j=0; j<nrepeat; ++j) { H.setRand('s');
       printf("\r%4d/%d ...\r",j+1,nrepeat); fflush(0);
       H.eig(U,E);
       if (D==7) {
          H.print("H"); U.print("U"); E.print("E");
       }
    }
    printf("\r%40s\r","");
};


