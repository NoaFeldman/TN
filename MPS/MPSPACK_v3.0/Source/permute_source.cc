// See ~/Source/GSL/gsl-1.8/permutation/permute_source.c
// Wb,Apr19,06

/* permutation/permute_source.c
   In-place Permutations 

   permute:    OUT[i]       = IN[perm[i]]     i = 0 .. N-1
   invpermute: OUT[perm[i]] = IN[i]           i = 0 .. N-1

   PERM is an index map, i.e. a vector which contains a permutation of
   the integers 0 .. N-1.

   From Knuth "Sorting and Searching", Volume 3 (3rd ed), Section 5.2
   Exercise 10 (answers), p 617

 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough */


char USAGE[]="";

#include "wblib.h"

double checkClock(char restart=1) {

    double t;
    mxArray *a;

#ifdef MATLAB_MEX_FILE
    mexCallMATLAB(1, &a, 0, NULL, "toc");
    t=mxGetPr(a)[0]; mxDestroyArray(a);

    if (restart)
    mexCallMATLAB(0, NULL, 0, NULL, "tic");
#else
    wblog(FL,"ERR %s() not available outsite MatLab",FCT);
#endif
    return t;
};

template <class T>
void gsl_permute (const unsigned *p, T *data, const unsigned n) {
   unsigned i, k, pk;
   T t;

   for (i=0; i<n; i++) {
       k=p[i]; while (k>i) k=p[k];
       
       if (k<i) continue;
       
       pk=p[k];
       
       if (pk==i) continue;
       

       t=data[i];
     
       while (pk!=i) { data[k]=data[pk]; k=pk; pk=p[k]; };
       
       data[k]=t;
   }
}

template <class T>
void gsl_permute (wbArray<T> &A, const wbvector<unsigned> &P) {
   
   unsigned i,j,k,l, r=A.SIZE.len, s=A.SIZE.prod(0);

   wbvector<unsigned> iP, S2(A.SIZE), pvec(s), Ivec(r);
   unsigned *p=pvec.data, *I=Ivec.data;

   if (r!=P.len) wblog(FL, "ERR Size mismatch.");
   if (r==0) return;

   A.SIZE.select(P,S2); l=r-1;
   getIPerm(P,iP);

   for (i=0; i<s; i++) {
       for (j=I[iP[0]],k=1; k<r; k++) j = j*A.SIZE[k] + I[iP[k]];

       p[i]=j;

       I[l]++; k=l;
       while(I[k]>=S2[k] && k>0) { I[k]=0; ++I[--k]; }
   }

   A.SIZE=S2;

printf("t=%g\n", checkClock(0));
   gsl_permute (p, A.DATA, s);
}

template <class T>
void meschach_permute (unsigned *p, T *data, const unsigned n) {

   unsigned i, old_i, start=0;
   T t;

   while (start<n) {

       old_i=start; i=p[old_i];

       if (i>=n) { start++; continue; }

       t = data[start];
       while (1) {
          data[old_i] = data[i];
          p[old_i]=i+n; old_i=i; i=p[old_i];

          if (i>=n) break;
          if (i==start) {
              data[old_i]=t;
              p[old_i] = i+n;
              break;
          }
       }
       start++;
   }
};

template <class T>
void meschach_permute (wbArray<T> &A, const wbvector<unsigned> &P) {
   
   unsigned i,j,k,l, r=A.SIZE.len, s=A.SIZE.prod(0);

   wbvector<unsigned> iP, S2(A.SIZE), pvec(s), Ivec(r);
   unsigned *p=pvec.data, *I=Ivec.data;

   if (r!=P.len) wblog(FL, "ERR Size mismatch.");
   if (r==0) return;

   A.SIZE.select(P,S2); l=r-1;
   getIPerm(P,iP);

   for (i=0; i<s; i++) {
       for (j=I[iP[0]],k=1; k<r; k++) j = j*A.SIZE[k] + I[iP[k]];

       p[i]=j;

       I[l]++; k=l;
       while(I[k]>=S2[k] && k>0) { I[k]=0; ++I[--k]; }
   }

   A.SIZE=S2;

printf("t=%g\n", checkClock(0));
   meschach_permute (p, A.DATA, s);
}



      
      
      
      
      
      
      

        
      
        



          

      


      
      















void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
){
    unsigned i,q0;
    wbvector<double> t(12);

    wbArray<double> A("1,2,3,4,5,6,7,6,5,4"),B;
    wbvector<double> d(8);
    wbvector<unsigned> P,iP;
    
    for (i=0; i<d.len; i++) d[i]=i+1;        d.print("d");
    P=Index(0,d.len-1).flip();               P.print("P");

    gsl_permute     (P.data, d.data, d.len); d.print("d");"P");
    meschach_permute(P.data, d.data, d.len); d.print("d");"P");

    A.setRand(); A.put("A0"); A.info("A0");
    P=Index(0,A.SIZE.len-1).flip(); getIPerm(P,iP);

    P.put("P"); iP.put("iP");
    
    checkClock();

    gsl_permute(A, P);      A.put("A1"); A.info("A1"); t[0]=checkClock();
    gsl_permute(A,iP);      A.put("A2"); A.info("A2"); t[1]=checkClock();

    meschach_permute(A, P); A.put("A3"); A.info("A3"); t[2]=checkClock();
    meschach_permute(A,iP); A.put("A4"); A.info("A4"); t[3]=checkClock();

    A.Permute( P);          A.put("A5"); A.info("A5"); t[4]=checkClock();
    A.Permute(iP);          A.put("A6"); A.info("A6"); t[5]=checkClock();

    A.RPermute( P);         A.put("A7"); A.info("A7"); t[6]=checkClock();
    A.RPermute(iP);         A.put("A8"); A.info("A8"); t[7]=checkClock();

    mxPermute(A, P);        A.put("A9"); A.info("A9"); t[8]=checkClock();
    mxPermute(A,iP);        A.put("AA"); A.info("AA"); t[9]=checkClock();

    q0=10;



    mxArray *ap, *p1=P.toMx(1), *ip=iP.toMx(1), *ai[2] = { A.toMx(), p1 };
    char astr[]="A_";
    checkClock();
    for (unsigned q=q0; q<t.len; q++) {
        ai[1] = q%2==0 ? p1 : ip;
        astr[1]='a'+q-q0;

#ifdef MATLAB_MEX_FILE
        mexCallMATLAB(1,&ap, 2,ai, "permute");
        i=mexPutVariable("base",astr,ap); t[q]=checkClock();
        if (i) wblog(FL,"ERR e=%d",i);
#else
        wblog(FL,"ERR %s() not available outsite MatLab",FCT);
#endif
        mxDestroyArray(ai[0]); ai[0]=ap;
    }

    mxDestroyArray(ai[0]);
    mxDestroyArray(ai[1]);


    t.print("Elapsed time");
};


