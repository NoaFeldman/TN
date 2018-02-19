#ifndef __WB_CMAT_CC__
#define __WB_CMAT_CC__

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* init0 directly dumps MatLab data -> leads to transpose! */

void wbCMat::init0(const mxArray *a) {
   unsigned m,n;
   double *d;

   if (!mxIsDblMat(0,0,a,'C'))
   wberror(FLINE,"Need double array on input.");

   isdiag=0;

   m=mxGetM(a); n=mxGetN(a); d=mxGetPr(a);
   R.init(n,m,d); 

   d=mxGetPi(a);
   if (d==NULL)
        I.init();  
   else I.init(n,m,d); 
};

void wbCMat::init(const mxArray *a) {

   unsigned i,j,l,m,n,nnz,d=0,k=0;
   double *dr,*di=NULL;
   mwIndex *Ir,*Jc;
   int cflag;

   if (!mxIsDouble(a) || mxGetNumberOfDimensions(a)!=2)
   wberror(FLINE,"Need double 2D array on input.");

   m=mxGetM(a); n=mxGetN(a); R.init(m,n); dr=mxGetPr(a);

   if (mxIsSparse(a)) {
      nnz=mxGetNzmax(a);
      Ir=mxGetIr(a); Jc=mxGetJc(a);

      cflag=mxIsComplex(a);
      if (cflag) { I.init(m,n); di=mxGetPi(a); }

      for (isdiag=1, l=j=0; j<n; j++) {
          d=unsigned(Jc[j+1]-Jc[j]); if (d==0) continue;

          for (k=0; k<d; k++,l++) {
             i=(unsigned)Ir[l]; if (isdiag) if (i!=j) isdiag=0;
             if (dr[l]!=0.) R.data[i*n+j]=dr[l]; if (cflag)
             if (di[l]!=0.) I.data[i*n+j]=di[l];
          }
      }
      if (l!=nnz && l>0) wblog(FL,
      "WRN Severe data inconsistency w/ sparse array (%d,%d)", l, nnz);
   }
   else {
      isdiag=0;

      for (i=0; i<m; i++) 
      for (j=0; j<n; j++) R.data[k++]=dr[i+m*j];

      di=mxGetPi(a); if (di==NULL) { I.init(); return; } 
      I.init(m,n); k=0;

      for (i=0; i<m; i++) 
      for (j=0; j<n; j++) I.data[k++]=di[i+m*j];
   }
};


char wbCMat::isDiag() const {
   unsigned i,j,l=0, cflag=(I.data!=NULL);

   for (i=0; i<R.dim1; i++)
   for (j=0; j<R.dim2; j++,l++) 
   if (i!=j) {
      if (R.data[l]!=0.) return 0;
      if (cflag)
      if (I.data[l]!=0.) return 0;
   }

   return 1;
}


void wbCMat::operator+=(const wbCMat &Z) {
   if (isdiag) if (!Z.isdiag) isdiag=0;

   R+=Z.R;
   if (Z.I.data) {
      if (I.data) I+=Z.I; else I=Z.I;
   }
}

void wbCMat::operator*=(const wbcomplex &z) {

   if (isdiag) {
       unsigned i,l, s=MIN(R.dim1,R.dim2), cflag;
       wbcomplex zz;

       if (I.data==NULL && z.i!=0) I.init(R.dim1,R.dim2);

       cflag=(I.data!=NULL);

       for (i=0; i<s; i++) {
           l = i*R.dim2 + i;
           zz = z * wbcomplex(R.data[l], cflag ? I.data[l] : 0.);
           R.data[l]=zz.r;
           if (cflag) I.data[l]=zz.i;
       }

       return;
   }

   if (z.i==0.) {
       R*=z.r; I*=z.r;
       return;
   }

   wbCMat Z(*this);
   init(R.dim1,R.dim2);
   (*this).plus(Z,z);
}



wbCMat& wbCMat::plus(
   const wbCMat &Z,
   const wbcomplex &z
){
   if (R.dim1!=Z.R.dim1 || R.dim2!=Z.R.dim2) wblog(FL,
     "ERR Dimension mismatch (%dx%d + %dx%d)",
      R.dim1, R.dim2, Z.R.dim1, Z.R.dim2);
   if (this==&Z) wblog(FL, "ERR IO space must not overlap.");

   if (isdiag) if (!Z.isdiag) isdiag=0;

   if (z.r!=0)             R.add(Z.R,  z.r);
   if (z.i!=0 && Z.I.data) R.add(Z.I, -z.i);

   if (z.i!=0)             I.add(Z.R, z.i, I.data==NULL);
   if (z.r!=0 && Z.I.data) I.add(Z.I, z.r, I.data==NULL);

   return *this;
}

wbCMat& wbCMat::minus(const wbCMat &Z, const wbcomplex &z) {

   if (R.dim1!=Z.R.dim1 || R.dim2!=Z.R.dim2) wblog(FL,
     "ERR Dimension mismatch (%dx%d + %dx%d)",
      R.dim1, R.dim2, Z.R.dim1, Z.R.dim2);

   if (isdiag) if (!Z.isdiag) isdiag=0;

   if (z.r!=0)             R.minus(Z.R,  z.r);
   if (z.i!=0 && Z.I.data) R.minus(Z.I, -z.i);

   if (z.i!=0)             I.minus(Z.R, z.i, I.data==NULL);
   if (z.r!=0 && Z.I.data) I.minus(Z.I, z.r, I.data==NULL);

   return *this;
}


void wbCMat::set(unsigned i, unsigned j, const double &r) {
    if (isdiag) if (i!=j) isdiag=0;
    wbcomplex z(r,0); set(i,j,z);
}

void wbCMat::set(unsigned i, unsigned j, const wbcomplex &z) {
    if (isdiag) if (i!=j) isdiag=0;

   if (i>=R.dim1 || j>=R.dim2)
   wblog(FL,"ERR Index out of bounds (%d,%d; %dx%d).",i,j,R.dim1,R.dim2);

   R(i,j)=z.r;
   
   if (z.i==0) return;

   if (!I.data) I.init(R.dim1,R.dim2);
   I(i,j)=z.i;
}


double wbCMat::maxRelDiff(const wbCMat &Z) const {

   unsigned i, s=R.dim1*R.dim2;
   double *r1=R.data, *r2=Z.R.data, *i1, *i2, d, maxdiff=0., maxval=0.;
   wbcomplex z1, z2;

   if (!sameSize(Z))
   wblog(FL, "ERR Incompatible data objects.");

   if (I.data==NULL && Z.I.data==NULL) {
      for (i=0; i<s; i++) {
          d=r2[i]-r1[i]; maxdiff=MAX(maxdiff, fabs(d));
          maxval=MAX(maxval, MAX(fabs(r1[i]), fabs(r2[i])) );
      }
   }
   else {
      i1=I.data, i2=Z.I.data;
      for (i=0; i<s; i++) {
          z1.set(r1[i], i1 ? i1[i] : 0.);
          z2.set(r2[i], i2 ? i2[i] : 0.);
          d=(z2-z1).abs(); maxdiff=MAX(maxdiff, fabs(d));
          maxval=MAX(maxval, MAX(z1.abs(),z2.abs()) );
      }
   }
   return ((maxval==0.) ? 0. : maxdiff/maxval);
}


void wbCMat::toWbComplex(wbMatrix<wbcomplex> &Z) const {
   unsigned i, s=R.dim1*R.dim2;

   Z.init(R.dim1, R.dim2);
   for (i=0; i<s; i++) Z.data[i].r=R.data[i];

   if (I.data) {
      if (R.dim1!=I.dim1 || R.dim2!=I.dim2) 
      wblog(FL,"ERR Dimension mismatch.");

      for (i=0; i<s; i++) Z.data[i].i=I.data[i];
   }
}


void wbCMat::Contract(
   char idx,         
   const wbCMat &M,  
   wbCMat &Aout,     
   char mflag        

) const {

   unsigned i;
   const char flags[]="NTCc", i00=(Aout.R.data==NULL),
              tflag[]="TNcC"; 

   for (i=0; i<4; i++) if (flags[i]==mflag) break;
   if (i==4) wblog(FL,"ERR Invalid flag %d<%d>.",mflag,mflag);

   if (idx==1) {
       mflag=tflag[i]; 
       wbCMatProd(M, *this, Aout, mflag, 'N', 1., 1., i00);
   }
   else
   if (idx==2) {
       wbCMatProd(*this, M, Aout, 'N', mflag, 1., 1., i00);
   }
   else
   wblog(FL,"ERR Invalid contraction index %d.", idx);

   return;
}


wbcomplex wbCMat::scalarProd(const wbCMat &Z) const {

    unsigned i,s=R.dim1*R.dim2;
    double *d1, *d2, rr=0., ri=0., ir=0., ii=0.;

    if (R.dim1!=Z.R.dim1 || R.dim2!=Z.R.dim2)
    wblog(FL,"ERR Dimensions mismatch.");

    d1=  R.data;
    d2=Z.R.data;     for (i=0; i<s; i++) rr+=d1[i]*d2[i];

    if (I.data) {
        d1=  I.data; for (i=0; i<s; i++) ir+=d1[i]*d2[i];
    }
    if (Z.I.data) {
        d1=  R.data;
        d2=Z.I.data; for (i=0; i<s; i++) ri+=d1[i]*d2[i];
    }
    if (I.data && Z.I.data) {
        d1=  I.data;
        d2=Z.I.data; for (i=0; i<s; i++) ii+=d1[i]*d2[i];
    }

    return wbcomplex(rr+ii, ir-ri); 
}


void wbCMat::mat2mxc(mxArray* C, unsigned idx) const {

    unsigned i,j,k;
    double *d;
    mxArray *a;

    if (!mxIsCell(C))
    wberror(FLINE, "Need cell structure on input.");

    if (idx>=(unsigned)mxGetNumberOfElements(C))
    wberror(FLINE, "Index out of bounds.");

    a=mxCreateDoubleMatrix(R.dim1, R.dim2, !I.data ? mxREAL : mxCOMPLEX);
    mxSetCell(C,idx,a);

    d=mxGetPr(a); k=0;

    for (j=0; j<R.dim2; j++) 
    for (i=0; i<R.dim1; i++) d[k++]=R.data[i*R.dim2+j];

    if (!I.data) return;

    if (I.dim1!=R.dim1 || I.dim2!=R.dim2) 
        wblog(FL,"ERR Severe dimension inconsistency (%dx%d; %dx%d)",
        I.dim1, I.dim2, R.dim1, R.dim2);

    d=mxGetPi(a); k=0;

    for (j=0; j<I.dim2; j++) 
    for (i=0; i<I.dim1; i++) d[k++]=I.data[i*I.dim2+j];
}


void wbCMat::put(const char *vname, const char *ws) const {

    unsigned i,j,k;
    double *d1;

    mxArray *a=mxCreateDoubleMatrix(R.dim1, R.dim2,
        !I.data ? mxREAL : mxCOMPLEX);

    d1=mxGetPr(a); k=0;

    for (j=0; j<R.dim2; j++) 
    for (i=0; i<R.dim1; i++) d1[k++] = R.data[i*R.dim2+j];

    if (I.data!=NULL) {
        if (I.dim1!=R.dim1 || I.dim2!=R.dim2)
        wblog(FL,"ERR Severe dimension inconsistency (%dx%d; %dx%d)",
        I.dim1, I.dim2, R.dim1, R.dim2);

        d1=mxGetPi(a); k=0;

        for (j=0; j<I.dim2; j++) 
        for (i=0; i<I.dim1; i++) d1[k++] = I.data[i*I.dim2+j];
    }

    i=mexPutVariable(ws, vname, a);
    if (i==1)
    wblog(FL,"ERR Could not put variable %s into workspace %s.",vname,ws);

    mxDestroyArray(a);
}


void wbCMat::info(const char *vname) const {

    char istr[32];

    sprintf(istr, "%ldx%ld", R.dim1, R.dim2);

    printf("    %-16s 0x%7lX %8s %s", vname, (unsigned long)R.data, istr,
    I.data ? "complex" : "double" );

    if (I.data)
         printf(" (I::0x%7lX)\n", (unsigned long)I.data);
    else printf("\n");

    if ((I.data==NULL && (I.dim1!=0 || I.dim2!=0)) ||
        (I.data!=NULL && (I.dim1!=R.dim1 || I.dim2!=R.dim2)))
    wblog(FL,"ERR Severe data inconsistency (I: 0x%7lX, %dx%d)",
    I.data, I.dim1, I.dim2);
}


#endif

