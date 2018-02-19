#ifndef __WB_MEXLIB_HH__
#define __WB_MEXLIB_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

   void mexDisp(const mxArray *a, const char *vname=0);

   double wbtoc(
      const char *F=NULL, int L=0,
      const char *istr="", char restart=0
   );

   void wbtic() { wbtoc(0,0); }

   template<class T>
   int mxGetNumber(const mxArray *a, T &d, const char qflag=0);

   double mxGetNumber(
      const char *F, int L,
      const mxArray *S, const char *vname, const char qflag=0
   );

   int getFlag(const mxArray *a);
   int mxGetString(const mxArray *a, char *s);
   int mxGetString(const mxArray *a, wbstring &s);

   wbstring mxTypeSize2Str(const mxArray *a);
   wbstring mxSize2Str(const mxArray *a);

   template<class T>
   void mxGetVector( const char* F, int L,
   const mxArray *S, const char *vname, wbvector<T> &v);

   template<class T>
   void mxGetVector(const char* F, int L,
   const mxArray *a, wbvector<T> &v);

   void getComplexData0(
   const mxArray *a, wbMatrix<double> &R, wbMatrix<double> &I);

   mxArray* matGetVariable(const char *F, int L,
      const char *fname, const char *vname, char force=0
   );
   mxArray* matGetVariableInfo(const char *F, int L,
      const char *fname, const char *vname, char force=0
   );

   double* mexSetCell2Matrix(mxArray *a, int fid, unsigned d1, unsigned d2);
   void mxSetFieldToNumber(mxArray *S, unsigned i, unsigned fid,
        wbcomplex z, char dblcheck=1);

   template<class T>
   int getNumGlobal(const char *vname, T &x, const char *ws="global");
   int isFlagGlobal(const char *vname);

   template<class T>
   mxArray* vecVec2Mx(
   const wbvector< wbvector<T> > &V, const char tflag=0);

   template<class T>
   mxArray* vecVec2Mx(
   const wbvector< wbvector<T>* > &V, const char tflag=0);

   template<class T>
   mxArray* vecMat2Mx(
   const wbvector< wbMatrix<T> > &M, const char rawflag=0);

   template<class T>
   inline mxArray* cpyRange2Mx(const T* d, const wbvector<unsigned> &S0);

   template <class T> inline
   size_t cpyRangeMx(const mxArray* a, T* b, size_t n=-1, char tcheck=1) {

      if (!a) return -1;

      mxClassID id=mxGetClassID(a);

      if (long(n)<0) n=mxGetNumberOfElements(a);

      if (id==mxDOUBLE_CLASS) {
         cpyRange((const double*)mxGetPr(a),b,n,tcheck);
      }
      else {
         const void *d=mxGetData(a);
         if (!d && n) wblog(FL,"ERR %s() got NULL pointer (%d)",FCT,n); 

         switch (id) {
            case mxCHAR_CLASS  : cpyRange(( mxChar *)d,b,n,tcheck); break;
            case mxINT8_CLASS  : cpyRange((  int8_T*)d,b,n,tcheck); break;
            case mxUINT8_CLASS : cpyRange(( uint8_T*)d,b,n,tcheck); break;
            case mxINT16_CLASS : cpyRange(( int16_T*)d,b,n,tcheck); break;
            case mxUINT16_CLASS: cpyRange((uint16_T*)d,b,n,tcheck); break;
            case mxINT32_CLASS : cpyRange(( int32_T*)d,b,n,tcheck); break;
            case mxUINT32_CLASS: cpyRange((uint32_T*)d,b,n,tcheck); break;
            case mxINT64_CLASS : cpyRange(( int64_T*)d,b,n,tcheck); break;
            case mxUINT64_CLASS: cpyRange((uint64_T*)d,b,n,tcheck); break;
            default: return -2;
         }
      }
      return n;
   };



   template <class T>
   void mxPut(const T& x, const char *vname="ans", const char* ws="caller");

   void mxPut(const wbvector<wbcomplex> &zz,
   const char* ="ans", const char* ="caller");

   template <class T>
   void matvec2cell(const wbvector< wbMatrix<T> > &IM, mxArray *C, int fid);

   void mxCellVec2Data(const mxArray *C, DATA &A);

   int mxAddField2Scalar(
       const char *F, int L,
       mxArray* S, const char *vname, mxArray *a=NULL
   );

   void mxAppendStructToStruct(
      const char *F, int L,
      mxArray *&S0, mxArray *S
   );

   int mxAddField(
       const char *F, int L,
       mxArray* S, const char *vname
   );

   void mxReplaceField(
       const char *F, int L,
       mxArray* S, unsigned k, int fid, mxArray *a
   );

   int mxUpdateField(
       const char *F, int L,
       mxArray* S, const char *name, unsigned k, mxArray *a
   );

   void matPutVariable(const char *F, int L,
      const char *file, const char *name, mxArray *a, unsigned keep=0
   );

    
namespace Wb {

   template <class TQ, class TD>
   mxArray* mxCreateSparse(
      const char *F, int L, unsigned d1, unsigned d2,
      const wbMatrix<TQ> &IJ, const wbvector<TD> &D,
      const PERM_T *p=NULL, int *info=NULL
   );

}

    
template<class T>
mxArray* numtoMx(const T &x) {
    mxArray *a=mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetPr(a)[0]=double(x); return a;
};


template<>
mxArray* numtoMx(const wbcomplex &z) {
    mxArray *a;
    if (z.i!=0.) {
       a=mxCreateDoubleMatrix(1,1,mxCOMPLEX);
       mxGetPr(a)[0]=z.r;
       mxGetPi(a)[0]=z.i;
    }
    else a=numtoMx(z.r);
    return a;
};

template<>
mxArray* numtoMx(const char &c) {
    if (isalpha(c)) {
       char s[2]={c,0};
       return mxCreateString(s);
    }
    else return numtoMx(double(c));
};


class mxStruct {
 public:

    mxStruct(const char *file, int line) : S(NULL) { init(file,line); };
   ~mxStruct() { 
       if (S) {
          wblog(FL,"WRN mxStruct() destroyed without further ado.");
          mxDestroyArray(S);
       }
   };

    void init(const char *file, int line) {
       F=file; L=line;
       if (S) mxDestroyArray(S);
       S=mxCreateStructMatrix(1,1,0,NULL);
    };

    void put(const char *name) {
       if (S) { mxPutAndDestroy(F.data,L,S,name); S=NULL; } else
       wblog(FL,"WRN mxStruct::put() not initialized yet - return");
       F.init(); L=0;
    };

    void save2Struct(mxArray *S0, const char *name) {
       if (S) {
          mxAddField2Scalar(F.data, L, S0, name, S);
          S=NULL;
       } else
       wblog(FL,"WRN mxStruct::save2Struct() not initialized yet");
       F.init(); L=0;
    };

    void checkNewField(const char *f=NULL, int l=0);

    mxStruct& operator()(const char *fn) { fld=fn; return *this; };

    template<class T>
    mxStruct& operator=(const T&x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, numtoMx(x));
       return *this;
    };

    mxStruct& operator=(mxArray *a) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, a);
       return *this;
    };

    template<class T>
    mxStruct& operator=(const wbvector<T>& x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, x.toMx());
       return *this;
    };

    template<class T>
    mxStruct& operator=(const wbMatrix<T>& x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, x.toMx());
       return *this;
    };

    mxArray *S;
    wbstring F, fld;
    unsigned L;

 protected:
 private:

};


template<>
mxStruct& mxStruct::operator=(const wbperm &P) {
   checkNewField(FL);
   mxAddField2Scalar(F.data, L, S, fld.data, P.toMx());
   return *this;
}

template<>
mxStruct& mxStruct::operator=(const wbindex &I) {
   checkNewField(FL);
   mxAddField2Scalar(F.data, L, S, fld.data, I.toMx());
   return *this;
}

inline void mxStruct::checkNewField(const char *f, int l) {
   int i; if (!f) { f=F.data; l=L; }

   if (!S) wblog(f,l,"ERR mxStruct() S not set yet.");
   if (fld.isEmpty()) wblog(f,l,"ERR mxStruct() field not specified.");

   i=mxGetFieldNumber(S,fld.data);
   if (i>=0) wblog(f,l,"ERR Field `%s' already exists!");
};



class MXPut {
 public:

    MXPut(const char *vn="ans")
     : S(0), vname(0),F(0), L(0) { init(0,0,vn); };

    MXPut(const char *file, int line, const char *vn="ans")
     : S(0), vname(0),F(0), L(0) {
       init(file,line,vn); if (file!=NULL) addFL();
    };

   ~MXPut() { put(); }

    void init() {
       if (F) { WB_DELETE(F); }
       if (vname) { WB_DELETE(vname); }
       L=0; if (S) { wblog(F_L,"WRN %s() S not reset yet (0x%lX) !??",S);
       S=0; }
    };

    void init(const char *file, int line, const char *vn=NULL) {
       if (F) { WB_DELETE(F); }
       if (file && file[0]) {
          WB_NEW(F,strlen(file)+1); strcpy(F,file); }

       if (vn && vn[0])
            { WB_NEW(vname,strlen(vn)+1); strcpy(vname,vn); }
       else { WB_NEW(vname,4); strcpy(vname,"ans"); }

       S=mxCreateStructMatrix(1,1,0,NULL); L=line;
    };

    void put(
       const char *ws __attribute__ ((unused)) = "caller",
       const char *vn=NULL
   ){
       if (S) {
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name !?",FCT,vname);
          const char s0[]="ans", *s=((vn && vn[0]) ? vn : vname);
          if (!s || !s[0]) s=s0;

#ifdef MATLAB_MEX_FILE
          #pragma omp critical
          if (F) wblog(F,L,"I/O putting structure '%s' to %s  ",s,ws);
          mxPutAndDestroy(NULL,0,S,s,ws); S=NULL;
#else
          wblog(F_L,"WRN cannot put '%s' to caller. ", s);
          mxDestroyArray(S); S=NULL;
#endif
       }
       init();
    };

    void save2(mxArray* &S0) { S0=S; S=NULL; init(); };

    mxArray* toMx() { mxArray* S0=S; S=NULL; init(); return S0; };

    void add2Struct(mxArray *S0) {
       if (S) { 
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name !?",FCT,vname);
          mxAddField2Scalar(F_L,S0,vname,S); S=NULL; 
       }
       init();
    };

    MXPut& addFL() {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       else {
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name !?",FCT,vname);
          if (mxGetFieldNumber(S,"FL")<0)
          mxAddField2Scalar(F_L,S,"FL",wbstring(shortFL(F_L)).toMx());
          else wblog(F_L,"WRN %s() already applied prior",FCT);
       }
       return *this;
    };

    template<class T>
    MXPut& add(const T&x, const char *fn) {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       else mxAddField2Scalar(F_L,S,fn,x.toMx());
       return *this;
    };

    template<class T>
    MXPut& addP(const T* x, const char *fn) {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       else wblog(F_L,
          "ERR %s() not defined yet for pointer to %s",
           FCT,getName(typeid(T)).data);
       return *this;
    };

    mxArray *S;
    char *vname, *F;
    unsigned L;

 protected:
 private:

};


template<>
MXPut& MXPut::addP(const mxArray* a, const char *fn) {
   if (!S) wblog(F_L,"ERR MXPut::%s() S not set yet",FCT); else
   if (!vname) wblog(F_L,
      "ERR %s() got empty variable name !?",FCT,vname);
   else mxAddField2Scalar(F_L,S, fn, (mxArray*)a);
   return *this;
};

#ifndef NOMEX
template<>
MXPut& MXPut::addP(const char* s, const char *fn)
{ return addP(mxCreateString(s),fn); };
#endif


template<>
MXPut& MXPut::add(const double &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const float &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const long double &x, const char *fn)
{ return addP(numtoMx(x),fn); };

#ifdef __WB_MPFR_HH__
template<>
MXPut& MXPut::add(const Wb::quad &x, const char *fn)
{ return addP(numtoMx(x),fn); };
#endif

template<>
MXPut& MXPut::add(const wbcomplex &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const char &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const unsigned &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const int &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const long &x, const char *fn)
{ return addP(numtoMx(x),fn); };

template<>
MXPut& MXPut::add(const unsigned long &x, const char *fn)
{ return addP(numtoMx(x),fn); };


#endif

