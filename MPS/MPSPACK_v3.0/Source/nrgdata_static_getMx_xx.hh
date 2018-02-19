#ifndef __WB_NRGDATA_HCC__
#define __WB_NRGDATA_HCC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

class NRGData {
  public:
    NRGData() : mx(NULL), fid(-1) {};
    NRGData(const char *s) : name(s), mx(NULL), fid(-1) {};
    NRGData(const NRGData &D) : name(D.name), mx(NULL), fid(D.fid) {};

   ~NRGData() {
       if (mx) {
          mxDestroyArray(mx);
          mx=NULL;
       }
       name.init(); fid=-1;
    };

    void init(const char *_name);

    bool isEmpty() const { return name.isEmpty(); };
    bool isQSpaceVec(int k=-1, char info=0);

    bool checkVec(unsigned r);
    bool checkVecVec(unsigned r);

    void update(const char *F, int L, mxArray *a, int k=-1);
    void print(const char *vstr="ans") const;

    const mxArray* getMx(int k=-1, char ionly=0);


    static const mxArray* getMx(const char *n);

    template<class T>
    static int getPara(const char *name, T &x, unsigned k, char i=0);

    static void updatePara(const char *F, int L,
    const char*, mxArray *, int k=-1);

    static void updateInfo(const char *F, int L,
    const char*, mxArray *, char keep=-1);

    static bool opExists(char kflag,
    NRGData &KK, NRGData &KT, NRGData &TK, NRGData &TT);

    static void SetupIO(const char *F, int L, const mxArray *a);
    static void SetupIO(const char *F, int L, const char *s);

    static void init() {
       if (mxNRG && !NAME.isEmpty())
       mexPutVariable("caller", NAME.data, mxNRG);

       NAME.init(); mxNRG=NULL; N=0; iter=0;
       getMx("");
    };

    wbstring name;
    mxArray *mx;
    int fid;


    static wbstring NAME;
    static mxArray* mxNRG;
    static unsigned N, iter;


  private:

    static const char* getFileName(unsigned k, char lenient=0) {
       if (!lenient && k>=N) wblog(FL,
       "ERR %s - index out of bounds (%d/%d)", __FUNCTION__, k, N);
       return getFileName_aux( k);
    };

    static const char* getFileNameI() {
       return getFileName_aux(-1);
    };

    static const char* getFileName_aux(int k) {
       static char f[128];
       if (k>=0)
            sprintf(f,"%s_%02d.mat", NAME.data, k);
       else sprintf(f,"%s_info.mat", NAME.data);
       return f;
    }
};

    wbstring NRGData::NAME="";
    mxArray* NRGData::mxNRG=NULL;
    unsigned NRGData::N=0;
    unsigned NRGData::iter=0;


void NRGData::SetupIO(const char *F, int L, const mxArray *a) {
   NAME.init(); N=0; iter=0;
   mxNRG=(mxArray*)a;
}

void NRGData::SetupIO(const char *F, int L, const char *s) {
   Wb::matFile f; unsigned n=strlen(s);
   mxArray *a;

   if (!n || n>118) wblog(FL,
   "ERR Invalid file/variable name `%s' (%d)", s, n);

   N=0; iter=0; mxNRG=NULL;
   NAME=s;

   a=mexGetVariable("caller",s);
   if (a) {
      wblog(F,L,"<i> I/O structure variable: %s", s);
      mxNRG=a;
      return;
   }

   if (!f.open(F,L, getFileName(0,'l'),"r"))
        wblog(F,L,"ERR Cannot open files of type %s_##.mat", s);
   else f.close();

   wblog(F,L,"<i> I/O files: %s_##.mat", s);
}


void NRGData::init(const char *s) {
   if (name.data || mx || fid>=0) { wblog(FL,
      "WRN Re-initializing variable `%s' (%s, 0x%lX, %d)",
       name.data ? name.data : "(N/A)", mx, fid);
      mx=NULL; fid=-1;
   }
   name=s;
} 


template<class T>
int NRGData::getPara(const char *name, T &x, unsigned k, char i) {

   int fid; mxArray *a;

   if (k>=N) wblog(FL,"ERR Index out of bounds (%d/%d)", k, N);

   if (mxNRG) {
      fid=mxGetFieldNumber(mxNRG, name);
      if (fid<0) {
         if (i) sprintf(str, "%s:%d Field `%s' could not be found.",
         FL, name); return 1;
      }

      a=mxGetFieldByNumber(mxNRG,k,fid); if (!a) {
         if (i) sprintf(str, "%s:%d Field `%s' is empty.",
         FL, name); return 1;
      }

      return mxGetNumber(a,x);

   } else {

      Wb::matFile F; int r;

      F.Open(FL, getFileName(k),"r");

      a=matGetVariable(F.mfp,name);

      if (!a || mxIsEmpty(a)) {
         if (i) sprintf(str, "%s:%d Field `%s' is empty.",
         FL, name); if (a); mxDestroyArray(a);
         return 1;
      }

      r=mxGetNumber(a,x); mxDestroyArray(a);
      return r;
   }
}


void NRGData::update(const char *F, int L,
   mxArray *a, int k
){
   if (k<0) k=iter; else
   if (k>=(int)N) wblog(FL,"ERR Index out of bounds (%d/%d)",k,N);

   updatePara(F,L, name.data, a, k);
}


void NRGData::updatePara(const char *F, int L,
   const char *name, mxArray *a, int k
){
   if (!name || !name[0]) wblog(FL,"ERR Invalid name (empty)");

   if (k<0) k=iter;
   if (mxNRG) {
      if (!NAME.isEmpty())
      mxUpdateField(F,L, mxNRG,name,k,a);
   }
   else {
      matPutVariable(F,L,getFileName(k),name,a);
   }
}


void NRGData::updateInfo(const char *F, int L,
   const char *name, mxArray *a, char keep
){
   if (!a) wblog(F,L,
   "WRN update info `%s' with NULL mxArray", name);

   if (mxNRG) {
      mexPutVariable("caller", name, a);
      if (!keep) mxDestroyArray(a);
   }
   else {
      matPutVariable(F,L,getFileNameI(),name,a,keep);
   }
}


bool NRGData::opExists(char tflag,
   NRGData &KK, NRGData &KT, NRGData &TK, NRGData &TT
){
   if (tflag) { 
      unsigned k=(N>2 ? N-2 : 0);
      return (TT.isQSpaceVec(k) &&
              KT.isQSpaceVec(k) && TK.isQSpaceVec(k));
   }
   else {
      return (KK.isQSpaceVec(0));
   }
}


bool NRGData::isQSpaceVec(int k, char info) {
    if (k<0) k=iter;

    if (k>=(int)N) wblog(FL,
    "ERR Index out of bounds (%d/%d)", k, N);

    if (isEmpty() || mxNRG && fid<0) {
       if (info) sprintf(str, "%s:%d NRG data `%s' does not exist.",
       FL, name.data); return 0;
    }

    const mxArray *a=getMx(k,'i');
    if (a && mxIsQSpaceVec(a)) return 1;

    return 0;
}


bool NRGData::checkVec(unsigned r) {

    if (mxNRG) {
       fid=mxGetFieldNumber(mxNRG, name.data);
       if (fid<0) {
          sprintf(str, "%s:%d Field `%s' could not be found.",
          FL, name.data); return 1;
       }

       if (!N) {
          const int *s=mxGetDimensions(mxNRG);
          const int n=mxGetNumberOfDimensions(mxNRG);
          if (n!=2 || (s[0]!=1 && s[1]!=1)) { sprintf(str,
             "%s:%d NRG data must be row structure vector (%d,%d)",
             FL, s[0], s[1]); return 1;
          }
          N=s[0]*s[1];
       }
       return !mxIsQSpaceVec(mxNRG, fid, r);
    }
    else {

       int i=mxIsQSpaceVec(NAME.data, name.data, r);
       if (i<0) return 1;
       if (!N) N=i; else
       if (i!=(int)N) {
          sprintf(str,"%s:%d Size mismatch (%d,%d)",
          FL, i, N); return 1;
       }
    }

    return 0;
};


bool NRGData::checkVecVec(unsigned r) {

    if (mxNRG) {
       fid=mxGetFieldNumber(mxNRG, name.data);
       if (fid<0) {
          sprintf(str, "%s:%d Field `%s' could not be found.",
          FL, name.data); return 1;
       }

       if (!N) {
          const int *s=mxGetDimensions(mxNRG);
          const int n=mxGetNumberOfDimensions(mxNRG);
          if (n!=2 || (s[0]!=1 && s[1]!=1)) { sprintf(str,
             "%s:%d NRG data must be row structure vector (%d,%d)",
             FL, s[0], s[1]); return 1;
          }
          N=s[0]*s[1];
       }
       return !mxsIsQSpaceVEC(mxNRG, fid, r);
    }
    else {

       int i=mxIsQSpaceVEC(NAME.data, name.data, r);
       if (i<0) return 1;
       if (!N) N=i;
       if (i!=(int)N) {
          sprintf(str,"%s:%d Size mismatch (%d,%d)",
          FL, i, N); return 1;
       }
    }

    return 0;
};


const mxArray* NRGData::getMx(int k, char ionly) {
    
if (k<0 && k!=-1) wblog(FL,"WRN k=%d", k);

    if (k<0) k=iter;
    if (k>=N) wblog(FL,"ERR getMx() - index out of bound (%d/%d)", k, N);

    if (mxNRG) {
         if (fid<0) return NULL;
         return mxGetFieldByNumber(mxNRG,k,fid);
    }
    else { 
       Wb::matFile F(FL,getFileName(k),"r");

       if (mx) mxDestroyArray(mx);
       if (ionly)
            mx = matGetVariableInfo(F.mfp, name.data);
       else mx = matGetVariable    (F.mfp, name.data);

       F.close(); return mx;
    }
};


const mxArray* NRGData::getMx(const char *n) {

    static mxArray* a=NULL;
    if (n==NULL) wblog(FL,"ERR calling getMx with NULL.");

    if (a) {
wblog(FL,"TST getMx(%s) - a=0x%lX", n ? n : "<>", a);
wblog(FL,"TST %s", mxGetClassName(a));
wblog(FL,"TST");
    }

    if (!n[0]) return NULL;

    if (mxNRG) {
       a=mexGetVariable("caller",n);
wblog(FL,"TST getMx(%s 0x%lX)", n, a);
    }
    else { 
       Wb::matFile F(FL,getFileNameI(),"r");
       a=matGetVariable(F.mfp, n);
wblog(FL,"TST getMx(%s 0x%lX)", n, a);
    }

    return a;
};


void NRGData::print(const char *vstr) const {
    
    printf("\nNRGData %s", vstr);

    printf("   name  :  `%s'", name.data);
    if (NAME.data)
         printf(" (%s: %s)\n", mxNRG ? "QSpaceRM" : "file", NAME.data);
    else printf("\n");

    printf("   mxNRG :  0x%lX (fid=%d", mxNRG, fid);
    if (mxNRG && fid>=0)
         printf(", %s)\n", mxGetFieldNameByNumber(mxNRG,fid));
    else printf("\n");
};


#endif

