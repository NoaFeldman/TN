#ifndef __WB_OPTS_HCC__
#define __WB_OPTS_HCC__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* class to manage optional input arguments */

class OPTS {

  public:

    OPTS () : mto(0) {};

    OPTS (const mxArray** ain, const unsigned n)
    : mto(0) { init(ain,n); };

    OPTS (const char *fname)
    : mto(0) { init(fname); };

    void init (const mxArray** ain, const unsigned n);
    void init (const char *fname);

   ~OPTS () { init(); }
    void init() { 
       name.init();

       if (MAT.mfp && aa.len) 
       for (unsigned i=0; i<aa.len; ++i) mxDestroyArray(aa[i]);

       aa.init();
    };

    int len() { return ( MAT.mfp ? -1 : (int)aa.len ); };

    template <class T>
    bool getOpt(const char*, int,
         const char* s, T &x, char force=0, const char *istr=0);

    template <class T>
    bool getOpt(const char*, int,
         const char* s, wbvector<T> &x, char force=0, const char *istr=0);

    template <class T>
    bool getOpt(const char*, int,
        const char* s, wbMatrix<T> &x, char force=0, const char *istr=0);

    bool getOpt(const char*, int, const char* s, const char *istr=0);

    template <class T>
    bool getOpt(const char* s, T &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0, s,x,force,istr);
    };

    template <class T>
    bool getOpt(const char* s,
       wbvector<T> &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0,x,force,istr);
    };

    template <class T>
    bool getOpt(const char* s,
       wbMatrix<T> &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0,x,force,istr);
    };

    bool getOpt(const char* s, const char *istr=0) {
       return getOpt(NULL,0,s,istr);
    };


    void checkAnyLeft(const char *MAT=NULL, int L=0);

  protected:
  private:

    wbvector<mxArray*> aa;
    wbvector<wbstring> name;

    unsigned mto;

    Wb::matFile MAT;

    int findOpt(const char* s);
};


int OPTS::findOpt(const char* s) {

   unsigned i, n=name.len;
   int k=-1;

   for (i=0; i<n; ++i)
   if (name[i].len) if (!strcmp(name[i].data,s)) {
      if (k<0) k=i; else ++mto;
   };

   return k;
}


void OPTS::init (const mxArray** ain, const unsigned n) {

   if (MAT.mfp) wblog(FL,
   "ERR OPTS already set by file >%s<", MAT.fname.data);

   aa.init(n, (mxArray**)ain);
   name.init(n);

   for (unsigned i=0; i<aa.len; ++i) {
       if (!mxIsChar(aa[i])) continue; else
       if (mxGetString(aa[i],str,128)) {
          wblog(FL,"WRN Error reading string in OPTS!?");
          continue;
       }

       name[i].init(str);
       aa[i]=NULL;"used" as name
   }
}


void OPTS::init (const char *fname) {

   if (aa.len) wblog(FL,
   "ERR OPTS already set by mxArray[] (%d)", aa.len);

   MAT.Open(FL, fname, "r");
}


bool OPTS::getOpt(
   const char *F, int L, const char* s, const char *istr
){ 
   unsigned r;

   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, s); r=(a!=NULL);
      if (a) {
         double dbl=0.;
         if (mxGetNumber(a,dbl)) wblog(FL,"ERR %s",str);
         if (dbl==0) r=0;
      }
      mxDestroyArray(a);
   }
   else {
      int i=findOpt(s); r=(i>=0);
      if (r) name[i].init();
   }

   if (r && F) {
      if (istr && istr[0])
           wblog(F,L," *  %s",istr);
      else wblog(F,L," *  %s",s);
   }

   return r;
}


template <class T>
bool OPTS::getOpt(
   const char *F, int L, const char* s, T &x0,
   char force, const char *istr
){ 
   unsigned k,e=0;
   mxArray *a; T x=0;

   if (MAT.mfp) {
      a = matGetVariable(MAT.mfp, s);
      if (!a) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }
   }
   else {
      int i=findOpt(s);
      if (i<0) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }

      k=(unsigned)(i+1);

      if (k>=name.len) ++e; else
      if (aa[k]==NULL) ++e; if (e) wblog(F_L,
      "ERR Input double missing in OPTS (%s).", s);

      a=aa[k];

      name[k-1].init(); aa[k]=NULL;
   }

   if (mxGetNumber(a,x)) wblog(F_L,
      "ERR reading OPT %s (%s)", s, str);

   if (F && x0!=x) {
      if (istr && istr[0])
           wblog(F,L," *  %-8s: %g [%s]", s,(double)x,istr);
      else wblog(F,L," *  %-8s: %g", s, (double)x);
   }

   x0=x;

   if (MAT.mfp) mxDestroyArray(a);

   return 1;
}


template<>
bool OPTS::getOpt(const char *F, int L,
   const char* s, wbstring &x, char force, const char *istr
){ 
   int lflag=0;
   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, s);
      if (!a) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }

      if (!mxIsChar(a)) wblog(F_L,
      "ERR Invalid type `%s' in OPTS (%s)", mxGetClassName(a), s);
      if (mxGetString(a,str,128)) wblog(F_L,
      "ERR Error reading string `%s' for OPTS (max.len=128)", s);

      if (x.data && !strcmp(x.data,str)) ++lflag;
      x=str; mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(s);
      if (i<0) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) ++e; else
      if (aa[k]!=NULL) ++e; if (e) wblog(F_L,
      "ERR Missing input string in OPTS (%s).", s);

      if (x.data && !strcmp(x.data,name[k].data)) ++lflag;
      x=name[k];

      name[k-1].init(); name[k].init();
   }

   if (F && lflag) {
      if (istr && istr[0])
           wblog(F,L," *  %-8s: %s [%s]",s,x.data,istr);
      else wblog(F,L," *  %-8s: %s", s, x.data);
   }

   return 1;
}

template<>
bool OPTS::getOpt(const char *F, int L,
   const char* s, mxArray* &a, char force, const char *istr
){ 
   if (MAT.mfp) {
      a=matGetVariable(MAT.mfp, s);
      if (a) aa.Append((mxArray*)a);
   }
   else {
      int i=findOpt(s);

      if (i<0) a=NULL;
      else {
         unsigned e=0, k=unsigned(i+1);

         if (k>=name.len) ++e; else
         if (aa[k]==NULL) ++e; if (e) wblog(F_L,
         "ERR Input double missing in OPTS (%s).", s);

         a=aa[k];

         name[k-1].init(); aa[k]=NULL;
      }
   }

   if (!a) {
      if (force) wblog(F_L,"ERR Missing input %s.", s);
      return 0;
   }

   if (F) {
      if (istr && istr[0])
           wblog(F,L," *  %-8s: (%s) [%s]",s,mxGetClassName(a),istr);
      else wblog(F,L," *  %-8s: (%s)",s,mxGetClassName(a));
   }

   return 1;
}


template<class T>
bool OPTS::getOpt(const char *F, int L,
   const char* s, wbvector<T>&x_, char force, const char *istr
){ 
   int lflag=0;
   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, s);
      if (!a) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }

      if (!mxIsDouble(a)) wblog(F_L,
      "ERR Invalid type `%s' in OPTS (%s)", mxGetClassName(a), s);

      wbvector<T> x; x.init(FL,a);
      if (x!=x_) { x.save2(x_); ++lflag; }
      mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(s);
      if (i<0) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) ++e; else
      if (aa[k]==NULL) ++e; if (e) wblog(F_L,
         "ERR Missing input data (vector) in OPTS (%s).", s);

      wbvector<T> x; x.init(FL,aa[k]);
      if (x!=x_) { x.save2(x_); ++lflag; }

      name[k-1].init(); aa[k]=NULL;
   }

   if (F && lflag) {
      if (istr && istr[0])
           wblog(F,L," *  %-8s: [%s] (%s)", s, x_.toStr().data,istr);
      else wblog(F,L," *  %-8s: [%s]", s, x_.toStr().data);
   }

   return 1;
}

template<class T>
bool OPTS::getOpt(const char *F, int L,
   const char* s, wbMatrix<T>&x_, char force, const char *istr
){ 
   int lflag=0;
   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, s);
      if (!a) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }

      if (!mxIsDouble(a)) wblog(F_L,
      "ERR Invalid type `%s' in OPTS (%s)", mxGetClassName(a), s);

      wbMatrix<T> x; x.init(FL,a);
      if (x!=x_) { x.save2(x_); ++lflag; }
      mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(s);
      if (i<0) {
         if (force) wblog(F_L,"ERR Missing input %s.", s);
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) ++e; else
      if (aa[k]==NULL) ++e; if (e) wblog(F_L,
         "ERR Missing input data (matrix) in OPTS (%s).", s);

      wbMatrix<T> x; x.init(FL,aa[k]);
      if (x!=x_) { x.save2(x_); ++lflag; }

      name[k-1].init(); aa[k]=NULL;
   }

   if (F && lflag) {
      if (istr && istr[0])
           wblog(F,L," *  %-8s: [%s] (%s)", s, x_.toStr().data,istr);
      else wblog(F,L," *  %-8s: [%s]", s, x_.toStr().data);
   }

   return 1;
}


void OPTS::checkAnyLeft(const char *F, int L) { 

   unsigned i,e=0;

   if (MAT.mfp) {
      MAT.init();
      for (i=0; i<aa.len; ++i) mxDestroyArray(aa[i]);
      return;
   }

   for (i=0; i<name.len; ++i) { if (name[i].len || aa[i]) {
       if ((++e)==1) printf("\n");

       if (name[i].len) {
            printf("%5d: '%s'\n", i+1, name[i].data);
       }
       else
       if (aa[i]) printf("%5d: ... %dx%d (%s)\n", i+1,
          (int)mxGetM(aa[i]),
          (int)mxGetN(aa[i]), mxGetClassName(aa[i])
       );
   }}

   if (e) {
      if (mto)
      wblog(F_L,"WRN Options appear more than once !?");

      printf("\n");
      wblog(F_L,"ERR unused/unrecognized options");
   }
}


#endif

