#ifndef __WB_DATA_HCC__
#define __WB_DATA_HCC__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* data - here without symmtries
 *    D[i].R=real(A(:,:,i))
 *    D[i].I=imag(A(:,:,i))
 */

class DATA {

  public:

    DATA () {};

    explicit DATA (const mxArray *a) { init(a); };

    DATA (unsigned d1, unsigned d2, unsigned len) {
       D.init(len);
       for (unsigned i=0; i<len; i++) D[i].R.init(d1,d2);
    };

    DATA (const DATA &A) { *this = A; };

   ~DATA () {};

    void initSize(const DATA &A) {   


       if (D.len!=A.D.len) D.init(A.D.len); 

       for (unsigned i=0; i<D.len; i++) {
           const wbCMat &Z=A.D[i];
           D[i].R.init(Z.R.dim1, Z.R.dim2);
           
           if (Z.I.data==NULL) {
              if (D[i].I.data) D[i].I.init();
              continue;
           }

           D[i].I.init(Z.I.dim1, Z.I.dim2);
       }
    };

    wbMatrix<double>& R(unsigned i);
    wbMatrix<double>& I(unsigned i);

    char isReal() const;
    char sameSize(const DATA&) const;

    DATA& operator=(const DATA&);
    DATA& operator*=(const wbcomplex &);

    void init(const mxArray *a);

    void init() { D.init(); };
    void init(unsigned d1, unsigned d2, unsigned len) {
       unsigned i;

       if (D.len!=len) D.init(len);

       for (i=0; i<len; i++) {
           D[i].R.init(d1,d2);
           if (D[i].I.data) D[i].I.init();
       }
    };

    void setR(unsigned i, double r);
    void setI(unsigned i, double r);

    void plus (const DATA &A, const wbcomplex &z, DATA &AOUT) const;
    void minus(const DATA &A, const wbcomplex &z, DATA &AOUT) const;

    void plus (const DATA &A, DATA &AOUT) const; 
    void minus(const DATA &A, DATA &AOUT) const; 

    DATA& plus (const DATA &A, const wbcomplex &z=1.);
    DATA& minus(const DATA &A, const wbcomplex &z=1.);

    wbcomplex scalarProd(const DATA &A) const;

    void sum(wbCMat &s) const;
    void sum(wbCMat &s, wbvector<wbcomplex> z) const;

    void puti(const char *vname="ans", unsigned i=0, const char *ws="base") const;
    void put  (const char *vname="ans", const char *ws="base") const;

    void put3D(const char *vname="ans", const char *ws="base") const;
    mxArray* dat2mx3D() const;

    void info(const char* file, int line, const char *istr="ans") const;

    void setRand(double fac=1., double shift=0.);

    wbvector<wbCMat> D;   

  protected:
  private:

};


inline DATA& DATA::operator=(const DATA& A) {
   unsigned i;

   if (D.len!=A.D.len) D.init(A.D.len);

   for (i=0; i<D.len; i++) D[i]=A.D[i];
   return *this;
}

void DATA::init(const mxArray *a) {

    unsigned i,j,k,l,m,n,s, d=mxGetNumberOfDimensions(a);
    const mwSize *ip=mxGetDimensions(a);
    double *dd,*dr,*di;

    char isreal;

    if (!mxIsDouble(a) || d!=3)
    wblog(FL,"ERR A and B must be 3D double arrays.");
    if (mxIsSparse(a))
    wblog(FL,"ERR A and B must not be sparse.");

    m=(unsigned)ip[0]; n=(unsigned)ip[1]; d=(unsigned)ip[2];

    dr=mxGetPr(a);
    di=mxGetPi(a); isreal=(di==NULL); s=m*n;

    if (D.len!=d) D.init(d);

    for (k=0; k<d; k++, dr+=s, di+=s) {
        D[k].R.init(m,n); dd=D[k].R.data; l=0;
        for (i=0; i<m; i++)
        for (j=0; j<n; j++) dd[l++]=dr[i+m*j]; 

        if (isreal) {
           if (D[k].I.data) D[k].I.init();
           continue;
        }

        D[k].I.init(m,n); dd=D[k].I.data; l=0;
        for (i=0; i<m; i++)
        for (j=0; j<n; j++) dd[l++]=di[i+m*j]; 
    }
}

void DATA::put3D(const char *vname, const char *ws) const {
    mxArray *a=dat2mx3D();
    mxPutAndDestroy(FL,a,vname,ws);
}

mxArray* DATA::dat2mx3D() const {

    mwSize dim[3] = {0,0,0};
    unsigned i,j,k,l,m,n,s, cflag=0, e=0;
    double *dd, *dr, *di=NULL;
    mxArray *a;

    if (D.len==0) {
       a=mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
       return a;
    }

    m=dim[0]=D[0].R.dim1, n=dim[1]=D[0].R.dim2, dim[2]=D.len;

    for (i=0; i<D.len; i++) {
       if (!D[i].checkSize() || D[i].R.dim1!=m || D[i].R.dim2!=n) e++;
       if (D[i].I.data) cflag++;
    }
    if (e) wblog(FL,
    "ERR Cannot reshape data into 3D array (%d,%d,%d).", m,n,D.len);

    a=mxCreateNumericArray(3, dim, mxDOUBLE_CLASS,
      cflag ? mxCOMPLEX : mxREAL);

    s=m*n*D.len*sizeof(double);

    dr=mxGetPr(a); memset(dr,0,s); if (cflag) {
    di=mxGetPi(a); memset(di,0,s); }

    s=m*n;

    for (k=0; k<D.len; k++, dr+=s, di+=s) {
        dd=D[k].R.data; l=0;
        for (j=0; j<n; j++)
        for (i=0; i<m; i++) dr[l++]=dd[i*n+j]; 

        if (!D[k].I.data) continue;

        dd=D[k].I.data; l=0;
        for (j=0; j<n; j++)
        for (i=0; i<m; i++) di[l++]=dd[i*n+j]; 
    }
    return a;
}


DATA& DATA::operator*=(const wbcomplex &z) {
   for (unsigned i=0; i<D.len; i++) D[i]*=z;
   return *this;
}


inline char DATA::isReal() const {

   for (unsigned i=0; i<D.len; i++) if (D[i].I.data) return 0;

   return 1;
}

inline char DATA::sameSize(const DATA &A) const {
   if (D.len!=A.D.len) return 0;
   for (unsigned i=0; i<D.len; i++) if (!D[i].sameSize(A.D[i])) return 0;
   return 1;
}


inline wbMatrix<double>& DATA::R(unsigned i) {
   if (i>D.len)
   wblog(FL,"ERR Index out of bounds (%d,%d).", i, D.len);
   return D[i].R;
}

inline wbMatrix<double>& DATA::I(unsigned i) {
   if (i>D.len)
   wblog(FL,"ERR Index out of bounds (%d,%d).", i, D.len);

   if (D[i].R.data && !D[i].I.data)
   D[i].I.init(D[i].R.dim1, D[i].R.dim2);

   return D[i].I;
}


void DATA::setR(unsigned i, double r) {
   if (i>D.len)
   wblog(FL,"ERR Index out of bounds (%d,%d).", i, D.len);
   D[i].R.set(r);
}

void DATA::setI(unsigned i, double r) {
   if (i>D.len)
   wblog(FL,"ERR Index out of bounds (%d,%d).", i, D.len);

   if (r!=0. && !D[i].I.data)
   D[i].I.init(D[i].R.dim1, D[i].R.dim2);

   D[i].I.set(r);
}


void DATA::plus(const DATA &A, const wbcomplex &z, DATA &AOUT) const {
    if (&A==&AOUT) wblog(FL,"ERR A and AOUT must not overlap.");
    AOUT=*this; AOUT.plus(A,z);
}

void DATA::plus(const DATA &A, DATA &AOUT) const { 
    if (&A==&AOUT) wblog(FL,"ERR A and AOUT must not overlap.");
    AOUT=*this; AOUT.plus(A);
}

DATA& DATA::plus(const DATA &A, const wbcomplex &z) {

    if (D.len!=A.D.len)
    wblog(FL,"ERR Dimensions mismatch.");

    for (unsigned i=0; i<D.len; i++)
    D[i].plus(A.D[i], z);

    return *this;
}


void DATA::minus(const DATA &A, const wbcomplex &z, DATA &AOUT) const {
    if (&A==&AOUT) wblog(FL,"ERR A and AOUT must not overlap.");
    AOUT=*this; AOUT.minus(A,z);
}

void DATA::minus(const DATA &A, DATA &AOUT) const { 
    if (&A==&AOUT) wblog(FL,"ERR A and AOUT must not overlap.");

    AOUT=*this;
    AOUT.minus(A);
}

DATA& DATA::minus(const DATA &A, const wbcomplex &z) {
    if (D.len!=A.D.len)
    wblog(FL,"ERR Dimensions mismatch.");

    for (unsigned i=0; i<D.len; i++)
    D[i].minus(A.D[i], z);

    return *this;
}


wbcomplex DATA::scalarProd(const DATA &A) const {
    wbcomplex z=0;

    if (D.len!=A.D.len)
    wblog(FL,"ERR Dimensions mismatch.");

    for (unsigned i=0; i<D.len; i++)
    z+=D[i].scalarProd(A.D[i]);

    return z;
}

void DATA::sum(wbCMat &S) const {

   if (D.len==0) wblog(FL,
   "ERR DATA::sum - data array of length 0.");

   S=D[0];
   for (unsigned i=1; i<D.len; i++) S+=D[i];
}

void DATA::sum(wbCMat &S, wbvector<wbcomplex> z) const {
   char first=1;

   if (D.len==0) wblog(FL,
   "ERR DATA::sum - data array of length 0.");

   if (D.len!=z.len) wblog(FL,
   "ERR DATA::sum - dimension mismatch (%d,%d).", z.len, D.len);

   for (unsigned i=0; i<D.len; i++) if (z[i]!=double(0)) {
      if (first)
         { S=D[i]; S*=z[i]; first=0; }
      else S.plus(D[i],z[i]);
   }

   if (first) S.init( D[0].R.dim1, D[0].R.dim2 );
}

void DATA::put(const char *vname, const char *ws) const {
   unsigned i,j,k,l,dim1,dim2;
   double *d1,*d2;
   mxArray *a, *C=mxCreateCellMatrix(1,D.len);

   for (l=0; l<D.len; l++) {
       dim1=D[l].R.dim1;
       dim2=D[l].R.dim2;

       a=mxCreateDoubleMatrix(dim1, dim2, D[l].I.data ? mxCOMPLEX : mxREAL);
       mxSetCell(C,l,a);

       d1=mxGetPr(a); d2=D[l].R.data; k=0;
       for (j=0; j<dim2; j++) 
       for (i=0; i<dim1; i++) d1[k++] = d2[i*dim2+j];

       if (D[l].I.data==NULL) continue;

       if (D[l].I.dim1!=D[l].R.dim1 || D[l].I.dim2!=D[l].R.dim2)
       wblog(FL,"ERR Dimension mismatch.");

       d1=mxGetPi(a); d2=D[l].I.data; k=0;
       for (j=0; j<dim2; j++) 
       for (i=0; i<dim1; i++) d1[k++] = d2[i*dim2+j];
   }

#ifdef MATLAB_MEX_FILE
   i=mexPutVariable(ws, vname, C);
   if (i==1) wblog(FL,
      "ERR Could not put variable %s into workspace %s.",vname,ws);
#else
   wblog(FL,"ERR %s(%s,%s,%lX) not available outsite MatLab",
   FCT,vname?vname:"",ws?ws:"",&C);
#endif

   mxDestroyArray(C);
}

void DATA::puti(const char *vname, unsigned i, const char *ws) const {
   const size_t n=16; size_t l; char istr[n];

   if (i>D.len) wblog(FL,
      "ERR index out of bounds (%d,%d)", i, D.len);

   l=snprintf(istr,n,"%s_R%0d", vname, i); if (l>=n) wblog(FL,
     "ERR %s() string out of bounds (%s; %d; %d/%d)",FCT,istr,i,l,n); 
   D[i].R.put(istr,ws);

   l=snprintf(istr,n,"%s_I%0d", vname, i); if (l>=n) wblog(FL,
     "ERR %s() string out of bounds (%s; %d; %d/%d)",FCT,istr,i,l,n); 
   D[i].I.put(istr,ws);
};

void DATA::setRand(double fac, double shift) {
   for (unsigned i=0; i<D.len; i++) {
       R(i).setRand(fac,shift);
       I(i).setRand(fac,shift);
   }
}


void DATA::info(const char* file, int line, const char *vname) const {
   size_t i=0; const size_t n=64; int l=0;
   char istr[n];

   wblog(file,line, "DATA structure %s", vname);
   for (; i<D.len; ++i) {
       l=snprintf(istr,n,"%s[%d]",vname,int(i));
       if (l<0 || unsigned(l)>=n) wblog(FL,
          "ERR %s() string out of bounds\n(%s; %d; %d/%d)",
          FCT,istr,i,l,n);
       D[i].info(istr);
   }
};


#endif

