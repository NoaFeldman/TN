// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

char USAGE[]="";

#include "wblib.h"


void tst_memcmp(double d1, double d2) {
    unsigned s=sizeof(double);

    wblog(FL, "memcmp(%4g, %4g) = %2d (%2d)", d1, d2,
      memcmp(&d1,&d2,s),
      memcmp(&d2,&d1,s)
    );
}


inline
int inl_function() { static int i=0; i++; return i; }

int reg_function() { static int i=0; i++; return i; }

void call_dgemm(unsigned n) { 
   wbMatrix<double> A(n,n), B(n,n), C;

   wbMatProd(A,B,C);

}


int main (int argc, char **argv) {


    WbClock tic("TestClock");
    WbClock toc("AuxClock");

    unsigned i,n=1<<16;
    double d;

    printf("\nTST call of clock itself (n=%d)\n\n",n);

    tic.resume(); for (i=0; i<n; i++) { toc.resume(); toc.stop(); }
    tic.stop();

    tic.info(); tic.reset();
    toc.info(); toc.reset();


    tic.resume(); for (i=0; i<n; i++) { inl_function(); }
    tic.stop();

    printf("\nTST call of inline function: %.3g sec/call\n\n",tic.getASM()/n);
    tic.info(); tic.reset();


    tic.resume(); for (i=0; i<n; i++) { reg_function(); }
    tic.stop();

    printf("\nTST call of regular function: %.3g sec/call\n\n",tic.getASM()/n);
    tic.info(); tic.reset();


    tic.resume(); for (i=0; i<n; i++) { d=sin(double(i)); }
    tic.stop();

    printf("\nTST call of sine() function: %.3g sec/call\n\n",tic.getASM()/n);
    tic.info(); tic.reset();


    tic.resume();
    for (unsigned D=0; D<10; D++) {
       toc.resume(); for (i=0; i<n; i++) { call_dgemm(D); }
       toc.stop();
       printf("  call of call_dgemm(%2d): %.3g sec/call\n",D,toc.getASM()/n);
       toc.reset();
    }

    printf("\nTST call of dgemm() function\n\n");
    tic.info(); tic.reset();

    return 0;
}

