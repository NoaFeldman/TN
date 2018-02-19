/*=================================================================
 *
 * YPRIME.C    Sample .MEX file corresponding to YPRIME.M
 *            Solves simple 3 body orbit problem 
 *
 * The calling syntax is:
 *
 *        [yp] = yprime(t, y)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2004 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.2 $ */
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define    T_IN    prhs[0]
#define    Y_IN    prhs[1]



#define    YP_OUT    plhs[0]

#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

static    double    mu = 1/82.45;
static    double    mus = 1 - 1/82.45;


static void testJIT_speed();
static void yprime(double yp[], double *t, double y[]);

void mexFunction( int nlhs, mxArray *plhs[], 
          int nrhs, const mxArray*prhs[] )
     
{ 
    double *yp; 
    double *t,*y; 
    unsigned int m,n; 
    
    
    if (nrhs != 3)
    mexErrMsgTxt("Three input arguments required."); 
    else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments."); 
    
    
    m = mxGetM(Y_IN); 
    n = mxGetN(Y_IN);
    if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) || 
    (MAX(m,n) != 4) || (MIN(m,n) != 1)) { 
    mexErrMsgTxt("YPRIME requires that Y be a 4 x 1 vector."); 
    } 
    
    YP_OUT = mxCreateDoubleMatrix(m, n, mxREAL); 
    
    yp = mxGetPr(YP_OUT);
    
    t = mxGetPr(T_IN); 
    y = mxGetPr(Y_IN);
        
    yprime(yp,t,y); 
    return;
    
};

static void yprime(double yp[], double *t, double y[]) {

    double    r1,r2;

    testJIT_speed();
    
    r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[2]*y[2]); 
    r2 = sqrt((y[0]-mus)*(y[0]-mus) + y[2]*y[2]);

    if (r1 == 0.0 || r2 == 0.0 ){
    mexWarnMsgTxt("Division by zero!\n");
    }
    
    yp[0] = y[1];
    yp[1] = 2*y[3]+y[0]-mus*(y[0]+mu)/(r1*r1*r1)-mu*(y[0]-mus)/(r2*r2*r2);
    yp[2] = y[3];
    yp[3] = -2*y[1] + y[2] - mus*y[2]/(r1*r1*r1) - mu*y[2]/(r2*r2*r2);

    return;
}


static void testJIT_speed() {

   double data[] = {
     0.17579619419810,
     0.98795081500675,
     0.95239245666016,
     0.04073345236204,
     0.46613719535510,
     0.87110782595894,
     0.73949946690158,
     0.88174190361038,
     0.66912314776612,
     0.81155247243690
  };

  double X1=data[0], X2=data[1], X3=data[2], X4=data[3];
  double T=data[4];

  unsigned i,N=(int)1E6;
  int it; double dt;

  mexPrintf("\nN=%7d: %9.4g, %9.4g, %9.4g, %9.4g, %9.4g",
  i, X1,X2,X3,X4, T);

  for (i=0; i<N; i++) {
      X1=(X1+X2+X3-X4)*T+i;
      X2=(X1+X2-X3+X4)*T+i;
      X3=(X1-X2+X3+X4)*T+i;
      X4=(-X1+X2+X3+X4)*T+i;
  }

  mexPrintf("\nN=%7d: %9.4g, %9.4g, %9.4g, %9.4g, %9.4g\n",
  i, X1,X2,X3,X4, T);

  return;
}


