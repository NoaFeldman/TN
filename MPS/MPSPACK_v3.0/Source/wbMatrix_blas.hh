#ifndef __WB_MATRIX_BLAS_HH__
#define __WB_MATRIX_BLAS_HH__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

void WbEigenSymmetric (const wbMatrix<double>&M,
     wbMatrix<double>&U, wbvector<double> &E);

template<class T>
void wbMatProd(
    const wbMatrix<T> &A,
    const wbMatrix<T> &B, wbMatrix<T> &C,
    char aflag='N', char bflag='N',
    const T afac=1., const T cfac=0.,
    const char i00flag=0
);

void wbCMatProd(
    const wbCMat &A0,
    const wbCMat &B0, wbCMat &C,
    char ia='N', char ib='N',
    const wbcomplex &afac=1., const wbcomplex &cfac=0.,
    const char i00flag=0
);

template<class T>
void wbMatProd(
    const wbMatrix<T> &A,
    const wbvector<T> &B, wbvector<T> &C,
    char aflag='N', const T afac=1., const T cfac=0.,
    const char i00flag=0
);

template<class T>
void wbSVD(
    const wbMatrix<T> &A,
    wbMatrix<T> &U,
    wbMatrix<double> &S,
    wbMatrix<T> &V
);


#endif

