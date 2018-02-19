#ifndef __WB_MEXLIB_CC__
#define __WB_MEXLIB_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

// NB! mexCallMATLAB(), mexDisp, mexPutVariable, mexGetVariablePtr, ...
// require matlab environment
#ifdef MATLAB_MEX_FILE

double wbtoc(
   const char *F, int L, const char *istr, char restart
){
   double t=-1.; mxArray *a;

   if ((F==NULL) ? 1 : (F[0]==0 ? 1 : 0)) {
      mexCallMATLAB(0, NULL, 0, NULL, "tic");
      return t;
   }

   mexCallMATLAB(1, &a, 0, NULL, "toc");
   t=mxGetPr(a)[0]; mxDestroyArray(a);

   if (istr) {
      if (!istr[0])
           wblog(F_L,"Elapsed time is %.6f seconds.", t);
      else wblog(F_L,"Elapsed time is %.6f seconds (%s).", t, istr);
   }

   if (restart)
   mexCallMATLAB(0, NULL, 0, NULL, "tic");

   return t;
};


template <class T>
void mxPut(const T& x, const char *vname, const char* ws) {
    mxArray *a=numtoMx(x);
    mxPutAndDestroy(FL, a, vname, ws);
}

template <>
void mxPut(
    const wbvector<wbcomplex> &zz,
    const char *vname,
    const char* ws
){
    unsigned i;
    double *dd;

    mxArray *a=mxCreateDoubleMatrix(zz.len, 1, mxCOMPLEX);

    dd=mxGetPr(a); for (i=0; i<zz.len; i++) dd[i]=zz.data[i].r;
    dd=mxGetPi(a); for (i=0; i<zz.len; i++) dd[i]=zz.data[i].i;

    i=mexPutVariable(ws, vname, a);
    if (i==1){
        printf("Could not put variable %s into workspace %s.",vname,ws);
        wberror(FLINE,"");
    }

    mxDestroyArray(a);
    return;
}

void mxAnsAndDestroy(
    const char *F, int L, mxArray *a, const char* ws, char disp
){
    wblog(F_L,"ERR %s() simply assign to argout[0]!",FCT);

    unsigned i=mexPutVariable(ws, "ans", a);

    if (i) wblog(F_L,
       "ERR %s() mexPutVariable(%s,..) returned (%d)",FCT,ws,i);

    if (disp) {
       mexDisp(a);
    }

    mxDestroyArray(a);
}

void mxPutAndDestroy(
    const char *F, int L,
    mxArray *a, const char *vname, const char* ws
){
    if (F) wblog(F_L,"I/O putting variable '%s' to %s",vname,ws);

#ifdef MATLAB_MEX_FILE
    unsigned i=mexPutVariable(ws,vname,a);
    if (i) wblog(F_L,"ERR %s() "
       "mexPutVariable(%s,'%s',..) returned (%d)",FCT,ws,vname,i);
#else
    wblog(F_L,"WRN cannot put variable %s to workspace %s",vname,ws);
#endif
    mxDestroyArray(a);
}

void mxPutArray(const char *F, int L,
    mxArray *a, const char *vname, const char* ws, const char* dstr) {

    int i = mexPutVariable(ws, vname, a);

    if (i) wblog(F_L,
    "ERR Could not put variable %s into workspace %s (%d).",vname,ws,i);

    if (dstr!=NULL) {
       if (!strcmp(dstr,"destroy")) mxDestroyArray(a);
       else if (dstr[0]!=0)
       wblog(F_L,"ERR Invalid dstr >%s<", dstr);
    }
}

template<class T>
int getNumGlobal(const char *vname, T &x, const char* ws){

    const mxArray *a=mexGetVariablePtr(ws, vname);
    double dbl;

    if (!a || !mxIsNumber(0,0,a)) return 1;
    dbl=mxGetPr(a)[0]; x=(T)dbl;
    
    if ((double)x!=dbl) wblog(FL,
      "WRN ()%s invalid value for type '%s' (%g)",
    FCT, getName(typeid(T)).data, dbl);

    return 0;
};


template<>
int getNumGlobal(const char *vname, wbcomplex &x, const char* ws){
    const mxArray *a=mexGetVariablePtr(ws, vname);

    if (!a || !mxIsNumber(0,0,a,'C')) return 1;
    x.r=mxGetPr(a)[0];

    if (mxIsComplex(a)) x.i=mxGetPi(a)[0];
    else x.i=0.;

    return 0;
};

template<>
int getNumGlobal(const char *vname, wbstring &s, const char *ws) {
    const mxArray *a=mexGetVariablePtr(ws, vname);

    char istr[64];

    if (!a || !mxIsChar(a)) return 1;
    if (mxGetString(a,istr,63)) return 1;
    s=istr;

    return 0;
}

void mexDisp(const mxArray *a, const char *vname) {
    if (vname && vname[0]) printf("\n%s = \n\n", vname);
    mexCallMATLAB(0, NULL, 1, (mxArray**) &a,"display");
}

#else

double wbtoc(
   const char *F, int L, const char *istr,
   char restart __attribute__ ((unused))
){
   wblog(FL,"WRN %s(%s,%s) not available outsite MatLab",
   FCT,shortFL(F_L),istr ? istr:""); return 0;
}

template <class T>
void mxPut(const T& x, const char *vname, const char* ws) {
   wblog(FL,"ERR %s(%s,%s,%lX) not available outsite MatLab",
   FCT,vname?vname:"",ws?ws:"",&x);
};

void mxAnsAndDestroy(
    const char *F, int L, mxArray *a, const char* ws, char disp
){
   wblog(F_L,"ERR %s(%s,%lX;%d) not available outsite MatLab",
   FCT,ws?ws:"",&a,disp);
};

void mxPutAndDestroy(
    const char *F, int L,
    mxArray *a, const char *vname, const char* ws
){
   wblog(F_L,"ERR %s(%s,%s,%lX) not available outsite MatLab",
   FCT,vname?vname:"",ws?ws:"",&a);
};

void mxPutArray(const char *F, int L,
    mxArray *a, const char *vname, const char* ws,
    const char* dstr
){
   wblog(F_L,"ERR %s(%s,%s,%lX,%s) not available outsite MatLab",
   FCT,vname?vname:"",ws?ws:"",&a,dstr?dstr:"");
};

template<class T>
int getNumGlobal(
   const char *vname __attribute__ ((unused)),
   T &x              __attribute__ ((unused)),
   const char* ws    __attribute__ ((unused))
){
   wblog(FL,"ERR %s(%s,%s) not available outsite MatLab",
   FCT,vname?vname:"",ws?ws:""); return 0;
};

void mexDisp(const mxArray *a, const char *vname) {
   wblog(FL,"ERR %s(%s,%lX) not available outsite MatLab",
   FCT,vname ? vname:"",a);
};

#endif


int getFlag(const mxArray *a) {

    int n=mxGetNumberOfElements(a);

    if (n>1 || mxIsSparse(a)) usage(FLINE, "Invalid flag (scalar).");

    if (mxIsDouble(a)) {

        double *d=mxGetPr(a);

        if (d==NULL)
             return 0;
        else return ((int)d[0]);
    }

    else if (mxIsChar(a)) {

        char istr[2];
        if (mxGetString(a,istr,2)) usage(FLINE, "Invalid flag (char).");
        return istr[0];
    }

    usage(FLINE, "Invalid flag (number).");
    return -1;
}


template<class T> inline 
int mxGetNumber(const mxArray *a, T &d, const char qflag) {
    double dbl; int e=0;

    if (!a) {
       quietErrLog(FL,"mxArray is NULL",qflag);
       d=T(0); return 1;
    }
    if (mxIsComplex(a)) { e=1; } else
    if (mxIsSparse(a)) { e=2; } else
    if (mxGetNumberOfElements(a)!=1) { e=3; } else
    if (mxIsDouble(a)) { dbl=mxGetPr(a)[0]; } else
    if (mxIsLogical(a)) { dbl=mxGetLogicals(a)[0]; }
    else { e=-1; }

    if (e) {
       size_t l=0, n=64; char msg[n];
       l=snprintf(msg,n,"input must be scalar numerical value (%s; %ld)",
          mxGetClassName(a), mxGetNumberOfElements(a));
       if (l>=n) wblog(FL,
          "ERR %s() string out of bounds (%d/%d)",FCT,l,n); 
       quietErrLog(FL,msg,qflag); return 1;
    }

    d=T(dbl);
    if (double(d)!=dbl) {
       size_t l=0, n=128; char msg[n];
       l=snprintf(msg,n,"input is of type (%s <> %s)",
          mxGetClassName(a), getName(typeid(T)).data);
       if (l>=n) wblog(FL,
          "ERR %s() string out of bounds (%d/%d)",FCT,l,n); 
       quietErrLog(FL,msg,qflag); return 1;
    }
    
    return 0;
};


template<>
inline int mxGetNumber(const mxArray *a, wbcomplex &z, const char qflag) {
    double *d;

    if (!a) {
       quietErrLog(FL,"mxArray is NULL",qflag); return 1;
    }
    if (mxGetNumberOfElements(a)!=1 || !mxIsDouble(a) || mxIsSparse(a)) {
       size_t l=0, n=64; char msg[n];
       l=snprintf(msg,n,"input must be (complex) scalar (%s; %ld)",
          mxGetClassName(a), mxGetNumberOfElements(a));
       if (l>=n) wblog(FL,
          "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
       quietErrLog(FL,msg,qflag); return 1;
    }

    d=mxGetPr(a); z.r=d[0];
    d=mxGetPi(a); z.i=(d!=NULL) ? d[0] : 0;

    return 0;
}

template<>
inline int mxGetNumber(const mxArray *a, char &c, const char qflag) {
    if (!a) {
       quietErrLog(FL,"mxArray is NULL",qflag); return 1;
    }
    if (mxIsChar(a)) {
       char istr[2];
       if (mxGetString(a,istr,2)) {
          size_t l=0, n=64; char msg[n];
          l=snprintf(msg,n,"input must be single char (%ld)",
             mxGetNumberOfElements(a));
          if (l>=n) wblog(FL,
             "ERR %s() string out of bounds (%d/%d)",FCT,l,n); 
          quietErrLog(FL,msg,qflag); return 1;
       }
       c=istr[0]; return 0;
    }
    else
    if (mxIsDouble(a) || mxIsLogical(a)) {
       double d; if (mxGetNumber(a,d,qflag)) return 1;
       if (double(char(d))!=d) {
          quietErrLog(FL,"input must be single value",qflag);
          return 1;
       }
       c=char(d);
    }
    else {
       quietErrLog(FL,"input must be single char",qflag);
       return 1;
    }

    return 0;
}


inline int mxGetString(const mxArray *a, char *s) {
   if (mxGetString(a,s,127)) {
      sprintf(str,"%s:%d ERR Could not read string (%s) !?", FL, s);
      return 1;
   }
   return 0;
}

inline int mxGetString(const mxArray *a, wbstring &s) {
   char istr[128]; istr[0]=0;
   if (mxGetString(a,istr,127)) {
      sprintf(str,"%s:%d ERR Could not read string (%s) !?", FL, istr);
      return 1;
   }
   s=istr; return 0;
}


wbstring mxSize2Str(const mxArray *a) {
   if (!a) { return "(null)"; }

   size_t l=0, n=32; char s[n];
   int i=0, r=mxGetNumberOfDimensions(a);
   const mwSize *S=mxGetDimensions(a);

   for (; i<r; ++i) {
       l+=snprintf(s+l,n-l,"%s%d", i ? "x":"", unsigned(S[i]));
       if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   }
   return s;
};

wbstring mxTypeSize2Str(const mxArray *a) {
   size_t n=64; char s[n];
   if (!a) { strcpy(s,"(null)"); return s; }

   size_t l=sprintf(s,"%s: ",mxGetClassName(a));
   int r=mxGetNumberOfDimensions(a);
   const mwSize *S=mxGetDimensions(a);

   for (int i=0; i<r; ++i) {
      l+=snprintf(s+l,n-l,"%s%d",i?"x":"", (unsigned)S[i]);
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);
   }
   return s;
};


template<class T>
inline void mxGetVector(
   const char* F, int L,
   const mxArray *S, const char *vname, wbvector<T> &v
){
   if (!mxIsStruct(S) || !mxIsScalar(S)) wblog(F_L,
   "ERR Input must be scalar structure.");

   mxArray *a=mxGetField(S,0,vname);

   if (!a) wblog(F_L,
   "ERR Field %s does not exist.", vname);

   mxGetVector(F_L,a,v);
}

template<class T>
inline void mxGetVector(
   const char* F, int L,
   const mxArray *a, wbvector<T> &v
){
   size_t i,n=mxGetNumberOfElements(a);
   double *dr, *di;

   if (!mxIsVector(a) || !mxIsDouble(a) || mxIsSparse(a))
   wberror(F_L, "Input must be regular vector.");

   v.init(n); dr=mxGetPr(a); di=mxGetPi(a);

   if (di && typeid(T)!=typeid(wbcomplex)) {
      wblog(F_L, "WRN Skip imaginary part.");
      di=NULL;
   }

   for (i=0; i<n; i++) DSET(v[i], dr, di, i);
}

double mxGetNumber(
   const char *F, int L,
   const mxArray *S, const char *vname, const char qflag
){
   mxArray *a;
   wbcomplex z;

   if (!mxIsStruct(S)) wberror(F_L,
   "Need scalar structure on input.");
   if (mxGetNumberOfElements(S)>1) wberror(F_L,
   "Need SINGLE structure on input.");

   a=mxGetField(S,0,vname); 
   
   if (!a) wblog(F_L,
   "ERR Field %s does not exist in structure.", vname);

   if (mxGetNumber(a,z,qflag)) wberror(FL,str);
   if (z.i!=0.) wblog(FL,"WRN Input number is complex! (%g)", z.i);

   return z.r;
}


void getComplexData0( 
    const mxArray *a,
    wbMatrix<double> &R,
    wbMatrix<double> &I
){
    unsigned s,m,n=mxGetNumberOfDimensions(a);
    double *d;

    if (n!=2) wberror(FL, "Array dimension exceeds two.");

    m=mxGetM(a);
    n=mxGetN(a); R.init(m,n); I.init(m,n); s=m*n*sizeof(double);

    d=mxGetPr(a); memcpy(R.data,d,s);

    d=mxGetPi(a); if (d!=NULL) memcpy(I.data,d,s);
}


template<class T>
inline mxArray* cpyRange2Mx(const T* d, const wbvector<unsigned> &S0){

   unsigned i, n=S0.prod(0), r=S0.len, len=MAX(2U,r);
   mxArray *a=NULL; double *dr;

   mwSize s[len];

   for (i=0; i<r; ++i) s[i]=S0[i];
   if (i<len) {
      if (n) for (; i<len; i++) s[i]=1;
      else   for (; i<len; i++) s[i]=0;
   }

   a=mxCreateNumericArray(len,s, mxDOUBLE_CLASS, mxREAL);
   if (a==NULL) {
      if (d!=NULL) wblog(FL,
         "ERR %s() invalid usage (%lX, %lX; %d)",FCT,d,a,len);
      return a;
   }
   if (!n) return a;

   dr=mxGetPr(a);

   if (dr==NULL) wblog(FL,
      "ERR failed to allocate mxArray ([%s ; 0x%lX])",S0.toStr().data,d);
   if (d==NULL && n) wblog(FL,
      "ERR invalid input array ([%s ; 0x%lX])",S0.toStr().data,d);

   for (i=0; i<n; i++) dr[i]=double(d[i]);
   return a;
}


template<>
inline mxArray* cpyRange2Mx(const char* d0, const wbvector<unsigned> &S0){

   unsigned i, r=S0.len, n=S0.prod(0), len=MAX(2U,r);
   mxArray *a;

   mwSize s[len];
   for (i=0; i<r; i++) s[i]=S0[i];
   if (r<2 && n) s[1]=1;

   a=mxCreateCharArray(len,s);

   if (a==NULL) {
      if (!d0) wblog(FL,
         "ERR %s() invalid usage (%lX, %lX)",__FUNCTION__,d0,a);
      else return a;
   }

   mxChar *d = (mxChar*)mxGetData(a);
   if (!d) wblog(FL,"ERR invalid mxCharArray");

   for (i=0; i<n; ++i) { d[i]=mxChar(d0[i]);
      if (char(d[i])!=d0[i]) wblog(FL,
         "ERR %s() got data conversion error (%d/%d)",
         FCT,int(d[i]),int(d0[i])
      );
   }

   return a;
};


template<>
inline mxArray* cpyRange2Mx(const wbcomplex* d, const wbvector<unsigned> &S0){

   unsigned i, r=S0.len, n=S0.prod(0), len=MAX(2U,r), isz=0;
   mxArray *a; double *dr, *di;

   mwSize s[len];
   for (i=0; i<r; i++) s[i]=S0[i];
   if (r<2 && n) s[1]=1;

   for (i=0; i<n; i++) { if (d[i].i!=0) { isz=1; break; }}

   a=mxCreateNumericArray(
   len, s, mxDOUBLE_CLASS, isz ? mxCOMPLEX : mxREAL);

   if (a==NULL) { if (d!=NULL) wblog(FL,
   "ERR %s() invalid usage (%lX, %lX)",__FUNCTION__,d,a); else return a; }

   dr=mxGetPr(a); if (dr==NULL) wblog(FL,"ERR invalid mxArray"); if (isz) {
   di=mxGetPi(a); if (di==NULL) wblog(FL,"ERR invalid mxArray"); }

   if (isz) {
          for (i=0; i<n; i++) { dr[i]=d[i].r; di[i]=d[i].i; } }
   else { for (i=0; i<n; i++) { dr[i]=d[i].r; } }

   return a;
}


mxArray* matGetVariable(
    const char *F, int L,
    const char *fname,
    const char *vname,
    char force
){
    mxArray *a; int i; Wb::matFile f;

    i=f.open(F,L,fname,"r", force); if (!i) return NULL;
    a=matGetVariable(f.mfp,vname);
    return a;
}

mxArray* matGetVariableInfo(
    const char *F, int L,
    const char *fname,
    const char *vname,
    char force
){
    mxArray *a; int i; Wb::matFile f;
     
    i=f.open(F,L,fname,"r", force); if (!i) return NULL;
    a=matGetVariableInfo(f.mfp,vname);
    return a;
}

    
inline double* mexSetCell2Matrix(
    mxArray *a,
    int fid,
    unsigned d1,
    unsigned d2
){
    mxSetCell(a,fid,mxCreateDoubleMatrix(d1,d2,mxREAL));
    return mxGetPr(mxGetCell(a,fid));
}


template <class T>
void matvec2cell(
    const wbvector< wbMatrix<T> > &IM,
    mxArray *C,
    int fid
){
    unsigned i,j,k,l=0,d1,d2;
    mxArray *A, *ai;
    wbMatrix<T> *im;

    double *d; T *data;

    mxSetCell(C, fid, mxCreateCellMatrix(IM.len,1));
    A=mxGetCell(C, fid);

    for (k=0; k<IM.len; k++) {
        im=IM.data+k; d1=im->dim1; d2=im->dim2; data=im->data;

        ai=mxCreateDoubleMatrix(d1,d2,mxREAL);
        mxSetCell(A,k,ai);
        d=mxGetPr(ai);

        for (j=0; j<d2; j++) 
        for (i=0; i<d1; i++) d[l++]=(double)(data[i*d2+j]+1); 
    }
}


void mxCellVec2Data(const mxArray *C, DATA &A) {

    unsigned int i,j,k,l,m,n,nnz, cflag,d,isdiag, N=mxGetNumberOfElements(C);
    double *dr,*di=0,*dr0,*di0=0;
    mxArray *a;
    mwIndex *Ir,*Jc;

    if (!mxIsCell(C) || !mxIsVector(C))
    wblog(FL,"ERR Input must be cell vector of double matrizes.");

    A.D.init(N);

    for (k=0; k<N; k++) {
        a=mxGetCell(C,k); cflag=mxIsComplex(a);

        if (!mxIsDouble(a))
        wberror(FL,"Input must be cell vector of double matrizes.");

        m=mxGetM(a);
        n=mxGetN(a);

        A.D[k].R.init(m,n); dr=A.D[k].R.data; dr0=mxGetPr(a); if (cflag) {
        A.D[k].I.init(m,n); di=A.D[k].I.data; di0=mxGetPi(a); }

        if (!mxIsSparse(a)) { 

            for (l=i=0; i<m; i++) 
            for (  j=0; j<n; j++) dr[l++]=dr0[i+m*j];

            if (cflag)
            for (l=i=0; i<m; i++) 
            for (  j=0; j<n; j++) di[l++]=di0[i+m*j];
        }
        else {
            nnz=mxGetNzmax(a);
            Ir=mxGetIr(a); Jc=mxGetJc(a);

            for (isdiag=1, l=j=0; j<n; j++) { 
               d=unsigned(Jc[j+1]-Jc[j]); if (d==0) continue;

               for (k=0; k<d; k++,l++) {
                  i=(unsigned)Ir[l]; if (isdiag) if (i!=j) isdiag=0;
                  if (dr0[l]!=0.) dr[i*n+j]=dr0[l]; if (cflag)
                  if (di0[l]!=0.) di[i*n+j]=di0[l];
               }
            }
            if (l!=nnz && l>0) wblog(FL,
            "WRN Severe data inconsistency w/ sparse arrays (%d,%d)",l,nnz);

            A.D[k].isdiag=isdiag;
        }
    }
}


namespace Wb {

template <class TQ, class TD>
mxArray* mxCreateSparse(
   const char *F, int L, unsigned d1, unsigned d2,
   const wbMatrix<TQ> &IJ,
   const wbvector<TD> &D,
   const PERM_T *p,
   int *info
){
   const char lex=-1;
   INDEX_T n=D.len;
   bool gotp=(p!=NULL);

   if (D.data && ISCOMPLEX(D.data[0])) wblog(F,L,
      "ERR %s() got complex sparse data set (%dx%d; %d)",FCT,d1,d2,D.len); 
   if (IJ.dim1!=n || (n && IJ.dim2!=2)) wblog(F,L,
      "ERR %s() size mismatch (%dx%d : ij=%dx%d @ %d)",
       FCT,d1,d2,IJ.dim1,IJ.dim2,D.len); 

   if (!gotp && !IJ.recsSorted(+1,lex)) {
      wbMatrix<TQ> IX(IJ);
      wbperm P;

      IX.sortRecs(P,+1,lex);

      return mxCreateSparse(F_L,d1,d2,IX,D,P.data);
   }

   mxArray *a=mxCreateSparse(d1,d2,D.len,mxREAL);
   if (a==NULL) wblog(F,L,
      "ERR %s() failed to allocate sparse matrix (%dx%d; %d)",
       FCT,d1,d2,D.len);
   if (!n) { return a; }

   double *sr=mxGetPr(a);
   mwIndex *is=mxGetIr(a), *js=mxGetJc(a),k=0;
   unsigned i,j=0, nz=0, m=0; TQ *ij=IJ.data;

   if (!WbUtil<TQ>::isInt()) {
      wblog(FL,"ERR %s() invalid non-integer data type '%s'",FCT,
      getName(typeid(TQ)).data);
   }

   for (i=0; i<n; ++i, ij+=2) {
      const TD& x=(gotp ? D.data[p[i]] : D.data[i]);
      if (int(ij[0])<0 || ij[0]>=d1 || int(ij[1])<0 || ij[1]>=d2) wblog(FL,
         "ERR %s() index out of bounds (%d,%d) [%dx%d]",
         FCT,ij[0],ij[1],d1,d2); 
      if (x!=0) {
         if (i) {
            if (ij[-2]==ij[0] && ij[-1]==ij[1]) {
               sr[k]+=x; ++m; continue;
            }
            else if (ij[-1]>ij[1] || (ij[-1]==ij[1] && ij[-2]>ij[0])) {
               wblog(FL,"ERR %s() data not sorted (%d/%d) !?",FCT,i+1,n);
            }
            else ++k;
         }
         sr[k]=x;
         is[k]=ij[0];

         for (; j<=ij[1]; ++j) js[j]=k;
      }
      else if (ij[0]!=ij[1]) ++nz;
   }

   for (++k; j<=d2; ++j) js[j]=k;

   if (F) {
      if (nz) { wblog(F,L,
         "WRN %s() got off-diag zero data (%d/%d) !?",FCT,nz,D.len); 
      }
      if (m) { wblog(F,L,
         "WRN %s() got non-unique sparse data (%d/%d)",FCT,m,D.len);
      }
   }
   if (info) (*info)=nz+m;


   return a;
};

}


int mxAddField2Scalar(
   const char *F, int L,
   mxArray* S, const char *vname, mxArray *a
){
   int fid;

   if (!mxIsStruct(S))
      wberror(F_L,"scalar input structure required");
   if (mxGetNumberOfElements(S)>1)
      wberror(F_L,"single input structure required");

   fid=mxAddField(S,vname); 
   
   if (fid<0) wblog(F_L,
      "ERR failed to add field '%s' !? (%d)", vname,fid);

   if (a) mxSetFieldByNumber(S,0,fid,a);
   return fid;
};


int mxAddField(
    const char *F, int L,
    mxArray* S, const char *name
){
    int fid = mxAddField(S, name);

    if (fid<0) wblog(F,L,
    "ERR Failed to add field %s to structure (%d) ???", name, fid);

    return fid;
}

    
void mxAppendStructToStruct(
   const char *F, int L,
   mxArray *&S0, mxArray *S
){
   unsigned i,n;
   const char* name;

   if (!mxIsStruct(S0) || mxGetNumberOfElements(S0)!=1 ||
       !mxIsStruct(S ) || mxGetNumberOfElements(S )!=1) wblog(F,L,
   "ERR %s requires to scalar structures on input.",__FUNCTION__);

   n=mxGetNumberOfFields(S0);

   for (i=0; i<n; i++) {
      name=mxGetFieldNameByNumber(S0,i);

      if (mxGetFieldNumber(S,name)>=0) wblog(F,L,
      "WRN overwriting field `%s' in structure.", name);

      mxAddField2Scalar(F,L,
        S, name, mxGetFieldByNumber(S0,0,i)
      );

      mxSetFieldByNumber(S0,0,i,NULL);
   }

   mxDestroyArray(S0); S0=NULL;
}


int mxUpdateField(
    const char *F, int L,
    mxArray* S, const char *name, unsigned k, mxArray *a
){
    int fid=mxGetFieldNumber(S,name);
    unsigned n;

    if (fid>=0) {
       mxArray *x=mxGetFieldByNumber(S,k,fid);
       if (x) mxDestroyArray(x);
    }
    else fid=mxAddField(F,L, S, name);

    n=(unsigned)mxGetNumberOfElements(S);

    if (k>=n) wblog(FL,
    "ERR %s: index out of bounds (%d/%d)",FCT,k,n);

    mxSetFieldByNumber(S,k,fid,a);


    return fid;
}


void matPutVariable(const char *F, int L,
   const char *file, const char *name, mxArray *a, unsigned keep
){
   Wb::matFile f;

   f.Open(F,L,file,"u");

   if (matPutVariable(f.mfp, name, a)) wblog(F,L,
      "ERR Failed to write parameter`%s' to file `%s'", name, file);

   f.close();
   if (!keep) mxDestroyArray(a);
}


void mxReplaceField(
    const char *F, int L,
    mxArray* S, unsigned k, int fid, mxArray *a
){
    if (!S || !mxIsStruct(S)) wblog(F_L,
       "ERR %s() got %s",FCT, S ? mxGetClassName(S):"(empty)");

    mxArray *x=mxGetFieldByNumber(S,k,fid);
    if (x) mxDestroyArray(x);

    mxSetFieldByNumber(S,k,fid,a);
};

    
void mxSetFieldToNumber(
    mxArray *S,
    unsigned i,
    unsigned fid,
    wbcomplex z,
    char dblcheck
){
    mxArray *a;
    double *dd;

    if (dblcheck) {
        if (!mxIsStruct(S))
            wblog(FL,"ERR Invalid input structure.");
        if (i>=(unsigned)mxGetNumberOfElements(S))
            wblog(FL,"ERR Element index out of bounds.");
        if (fid>=(unsigned)mxGetNumberOfFields(S))
            wblog(FL,"ERR Field index out of bounds.");
    }

    a=mxCreateDoubleMatrix(1,1, z.i==0 ? mxREAL : mxCOMPLEX);
    mxSetFieldByNumber(S,i,fid,a);

    dd=mxGetPr(a); dd[0]=z.r; if (z.i==0) return;
    dd=mxGetPi(a); dd[0]=z.i;
}

    
int isFlagGlobal(const char *vname) {
    double flag=0;
    if (getNumGlobal(vname,flag)) return 0;
    return (flag!=0);
}

    
template<class T>
mxArray* vecVec2Mx(
    const wbvector< wbvector<T> > &V,
    const char tflag
){
    mxArray *C=mxCreateCellMatrix(1,V.len);

    for (unsigned i=0; i<V.len; i++)
    mxSetCell(C, i, V[i].toMx(tflag));

    return C;
}

template<class T>
mxArray* vecMat2Mx(
    const wbvector< wbMatrix<T> > &M,
    const char rawflag
){
    mxArray *C=mxCreateCellMatrix(1,M.len);

    for (unsigned i=0; i<M.len; i++)
    mxSetCell(C, i, M[i].toMx(rawflag));

    return C;
}

template<class T>
mxArray* vecVec2Mx(
    const wbvector< wbvector<T>* > &V,
    const char tflag
){
    mxArray *C=mxCreateCellMatrix(1,V.len);
    wbvector<T> emptyVec;

    for (unsigned i=0; i<V.len; i++)
    if (V[i]) mxSetCell(C,i,V[i]->toMx(tflag));
    else mxSetCell(C, i, emptyVec.toMx(tflag));

    return C;
}


#endif

