#ifndef __WB_MXLIB_HH__
#define __WB_MXLIB_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int isHelpIndicator(const mxArray *a) {
 
   char istr[4];

   if (mxGetString(a,istr,4)) return 0;

   if (!strcmp(istr, "-?") || !strcmp(istr, "-h")) return 1;

   return 0;
}

#ifndef NOMEX
int isHelpIndicator(const char *istr) {
 
   if (istr && (!strcmp(istr, "-?") || !strcmp(istr, "-h"))) return 1;
   return 0;
}
#endif

int mxIsEqual(const mxArray *a, const char *s) {
 
   int n=strlen(s)+1;
   char istr[n+1]; istr[0]=0;

   if (mxGetString(a,istr,n)) return 0;
   if (!strcmp(istr,s)) return 1;

   return 0;
}


   void mxPutArray(const char *file, int line,
        mxArray *a, const char *vname="ans", const char* ws="base",
        const char* dstr=""); 

   void mxAnsAndDestroy(
   const char *file, int line, mxArray *a, const char* ="caller", char disp=1);

   void mxPutAndDestroy(const char *file, int line,
   mxArray *a, const char *vname="ans", const char* ws="base");


inline int mxMatSafeGuard(
   const char *F, int L,
   const mxArray* a, unsigned r=2, char cflag=0, char type='d'
){
   if (!a) {
      if (L) {
         strcpy(str,"got a=null"); 
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return 1;
   }

   if (type=='d') {
      if (!mxIsDouble(a)) {
         if (L) {
            sprintf(str,"got data of type %s",mxGetClassName(a));
            if (F) wblog(F,L,"ERR %s() %s",FCT,str);
         }
         return 2;
      }
   }
   else if (type=='*') {
      if (!mxIsNumeric(a)) {
         if (L) {
            sprintf(str,"got data of type %s",mxGetClassName(a));
            if (F) wblog(F,L,"ERR %s() %s",FCT,str);
         }
         return 3;
      }
   }
   else wblog(F_L,"ERR %s() invalid type=%d<%c>",FCT,type,type);

   if (mxIsSparse(a)) {
      if (L) {
         strcpy(str,"got sparse input");
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return 4;
   }

   if (!cflag && mxIsComplex(a)) {
      if (L) {
         strcpy(str,"got complex input");
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return 5;
   }

   if (int(r)>0) {
      unsigned s=mxGetNumberOfDimensions(a);
      if (s>r) {
         if (s==2) {
            const mwSize *D=mxGetDimensions(a);
            if (D[0]==1 || D[1]==1) s=1;
         }
         if (s>r) {
            if (L) {
               sprintf(str,"got invalid rank (r=%d/%d)",s,r);
               if (F) wblog(F,L,"ERR %s() %s",FCT,str);
            }
            return 6;
         }
      }
   }

   return 0;
};


inline bool mxIsDblMat(const char *F, int L,
   const mxArray* a, char cflag=0, char dflag='d')
 { return !mxMatSafeGuard(F,L,a,2,cflag,dflag); }


inline bool mxIsDblArr(const char *F, int L,
   const mxArray* a, char cflag=0, char dflag='d')
 { return !mxMatSafeGuard(F,L,a,-1,cflag,dflag); }


inline bool mxIsDblVector(const char *F, int L,
   const mxArray *a, char cflag=0, char dflag='d')
 { return !mxMatSafeGuard(F,L,a,1,cflag,dflag); }

inline bool mxIsDblScalar(const char *F, int L,
   const mxArray *a, char cflag=0, char dflag='d') {
   if (mxMatSafeGuard(F,L,a,1,cflag,dflag)) return 0;
   return mxGetNumberOfElements(a)==1;
}

inline bool mxIsNumber(const char *F, int L,
   const mxArray *a, char cflag=0, char dflag='*') {
   if (mxMatSafeGuard(F,L,a,1,cflag,dflag)) return 0;
   return (mxGetNumberOfElements(a)==1);
}


inline bool mxIsScalar(const mxArray *a) {
   return (mxGetNumberOfDimensions(a)==2 && mxGetNumberOfElements(a)==1);
}

inline bool mxIsVector(const mxArray *a) {
   const mwSize *s=mxGetDimensions(a);
   const int n=mxGetNumberOfDimensions(a);
   return (n==2 && (s[0]==1 || s[1]==1));
}

inline bool mxIsIndex(const mxArray *a, int base=0) {
   const mwSize *s=mxGetDimensions(a);
   int const r=mxGetNumberOfDimensions(a);
   if (r==2 && (s[0]==1 || s[1]==1) && !mxIsComplex(a)) {
      const double *d=mxGetPr(a); if (!a) return 0;
      for (unsigned n=s[0]*s[1], i=0; i<n; i++) {
         if (d[i]<base || d[i]!=round(d[i])) return 0;
      }
      return 1;
   };
   return 0;
}


#ifndef MATLAB_MEX_FILE

#define mexPutVariable mexPutVariable_wbx
int mexPutVariable_wbx(const char *ws, const char *vname, const mxArray *a){
    wblog(FL,"WRN mexPutVariable(%s,%lX,%s) not available outside MatLab",
    vname?vname:"",a,ws?ws:""); return 1;
};

#define mexCallMATLAB mexCallMATLAB_wbx
int mexCallMATLAB_wbx(
    int nargout, mxArray **argin, int nargin, mxArray **argout,
    const char *cmd
){  wblog(FL,"ERR mexCallMATLAB(%d,%lX,%d,%lX,%s) not available outside MatLab",
    nargin,argin,nargout,argout,cmd?cmd:""); return 1;
};

#define mexGetVariable mexGetVariable_wbx
int mexGetVariable_wbx(const char *ws, const char *vname) {
    wblog(FL,"WRN mexGetVariable(%s,%s) not available outside MatLab",
    vname?vname:"",ws?ws:""); return 1;
};

#define mexGetVariablePtr mexGetVariablePtr_wbx
mxArray* mexGetVariablePtr_wbx(const char *ws, const char *vname) {
    wblog(FL,"WRN mexGetVariablePtr(%s,%s) not available outside MatLab",
    vname?vname:"",ws?ws:""); return NULL;
};

#define mexPrintf mexPrintf_wbx
int mexPrintf_wbx(const char *fmt, ...){
    wblog(FL,"WRN mexPrintf(%s,...) not available outside MatLab",
    fmt?fmt:""); return 1;
};

#define mexFunctionName mexFunctionName_wbx
const char* mex_function_name="(program)";
const char* mexFunctionName_wbx(const char *fmt, ...){
    wblog(FL,"WRN mexFunctionNamef(%s,...) not available outside MatLab",
    fmt?fmt:""); return mex_function_name;
};

#endif


#endif

