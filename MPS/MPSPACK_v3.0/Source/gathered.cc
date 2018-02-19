#ifndef __WB_GATHERED_CC__
#define __WB_GATHERED_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
    
void initSigHandler(const char i) {

   static sighandler_t prev_SIGINT;
   static sighandler_t prev_SIGABRT;
   static unsigned iflag=0;

   if (i=='i') {
      if (!iflag) {
         prev_SIGINT  = signal(SIGINT,  WBSigHandler);
         prev_SIGABRT = signal(SIGABRT, WBSigHandler); iflag=1;
      }
   }
   else
   if (i=='r') {
      if (!iflag) { wblog(gsh_F, gsh_L,
      "WRN SigHandler not initialized yet or already restored.");
      return; }

      signal(SIGINT,  prev_SIGINT );
      signal(SIGABRT, prev_SIGABRT); iflag=0;
   }
   else
   wblog(gsh_F, gsh_L,"ERR Invalid flag=%c<%d>",i,i);
}

void WBSigHandler(int i) {

   static unsigned NC0=3, icount=NC0, acount=NC0;
   fflush(0);

   if (i==99) {
      if (icount!=NC0 || acount!=NC0) {
         unsigned ic=icount, ac=acount;
         icount=acount=NC0;

         wblog(gsh_F, gsh_L, "ERR RECEIVED %s (%d) - TERMINATE",
         ic<ac ? "SIGINT":"SIGABRT", NC0-MIN(ic,ac));
      }
      return;
   }


   if (i==SIGINT ) {
      char istr[]="SIGINT";
      if (icount<=0 || icount>=NC0) printf("\n");
      if (icount<=0) {
         wblog(gsh_F,gsh_L,"ERR RECEIVED MULTIPLE %s - EXIT!%N",istr);
         wblog(FL,"TST %s()",FCT);
         fflush(0); exit(-1);
      }
      wblog(gsh_F, gsh_L, "*** received %s (%d) ***",istr,icount--);
      return;
   }

   if (i==SIGABRT) {
      char istr[]="SIGABRT";
      if (acount<=0 || acount>=NC0) printf("\n");
      if (acount<=0) { printf("\n"); 
         wblog(gsh_F,gsh_L,"ERR RECEIVED MULTIPLE %s - EXIT!%N",istr);
         wblog(FL,"TST %s()",FCT);
         fflush(0); exit(-1);
      }
      wblog(gsh_F, gsh_L, "*** received %s (%d) ***",istr,acount--);
      return;
   }

   wblog(gsh_F,gsh_L,"%N*** RECEIVED SIGNAL=%d ***\n*** Exit!", i);
   fflush(0); exit(1);
}





void dbstop(const char* F, int L) {

   wblog(F,L,"%NSIG Raising sigtrap ...%N");

   raise(SIGTRAP);
}

    
template<class T> inline 
wbstring num2Str(const T &x, const char *f) {
    if (f && f[0]) {
       wbstring s(32); snprintf(s.data,s.len,f,x);
       return s;
    }
    else return wbstring()<<x;
}


template<> inline 
wbstring num2Str(const float &x, const char *f) {
    if (f && f[0]) {
       wbstring s(32); snprintf(s.data,s.len,f,x);
       return s;
    }
    else return wbstring()<<x;
}

template<> inline 
wbstring num2Str(const double &x, const char *f) {
    if (f && f[0]) {
       wbstring s(32); snprintf(s.data,s.len,f,x);
       return s;
    }
    else return wbstring()<<x;
}

template<> inline 
wbstring num2Str(const wbcomplex &x, const char *f) {
    return x.toStr(f && f[0] ? f : "%.4g");
}


template<class T> inline
int num2int(const char *F, int L, const T &x){
    int i=int(x+0.5);
    if (F) {
       double e=fabs(x-i);
       if (e>1E-14) { wblog(F_L,"ERR %s() got %g !??",FCT,double(x)); }
    }
    if (T(i)<0) i=-1;
    return i;
};

template<> inline 
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const char &x){ return x; }

template<> inline 
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned char &x){ return x; }

template<> inline 
int num2int( const char *F __attribute__ ((unused)),
   int L  __attribute__ ((unused)), const int &x){ return x; }

template<> inline 
int num2int( const char *F __attribute__ ((unused)),
   int L  __attribute__ ((unused)), const long &x){ return x; }

template<> inline 
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned &x){ return x; }

template<> inline 
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned long &x){ return x; }


template<class T> inline
unsigned checkInt(const char *F, int L,
   const T* x, size_t n, double eps
){
    for (size_t i=0; i<n; ++i) {
        if (fabs(x[i]-round(x[i]))>eps) {
           if (F) wblog(F_L,
              "ERR %s() got data[%d]=%g (having %s) !??",FCT,
               i+1, double(x[i]), getName(typeid(T)).data);
           return (i+1);
        }
    }
    return 0;
};

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const int* x      __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const unsigned* x __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const char* x     __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }


inline const char* Wb::basename(const char *data, char c) {

    int i, n=strlen(data);
    for (i=n-1; i>=0; i--) if (data[i]==c) { i++; break; }

    return (data + (i>0 ? i : 0));
};


inline const char* Wb::strchri(const char *s, char c) {

   if (!s) wblog(FL,"ERR %s() got invalid null string",FCT);

   for (c=tolower(c); *s; ++s) {
      if (tolower(*s)==c) return s;
   }

   return NULL;
};


inline const char* Wb::strstri(const char *s1, const char *s2) { 

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%lX, %lX)",FCT,s1,s2);

   const char *a, *b;

   for(; *s1; ++s1) {
      for(a=s1, b=s2; *a && tolower(*a)==tolower(*b); ++a, ++b) {};
      if(!*b) return s1;
   }

   return NULL;    
};


inline const char* Wb::strstrw(const char *s1, const char *s2) { 

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%lX, %lX)",FCT,s1,s2);

   unsigned i=0, j=0, k=0;

   for(; s1[i]; ++i) { 
      for (j=i, k=0; s1[j] && s1[j]==s2[k]; ++j, ++k) {};
      if (!s2[k]) {
         if ((!i || !Wb::isword(s1[i-1])) && !Wb::isword(s1[j]))
         return s1+i;
      }
   }

   return NULL;    
};


inline const char* Wb::strstrwi(const char *s1, const char *s2) { 

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%lX, %lX)",FCT,s1,s2);

   unsigned i=0, j=0, k=0;

   for(; s1[i]; ++i) { 
      for (j=i, k=0; s1[j] && tolower(s1[j])==tolower(s2[k]); ++j, ++k);
      if (!s2[k]) {
         if ((!i || !Wb::isword(s1[i-1])) && !Wb::isword(s1[j]))
         return s1+i;
      }
   }

   return NULL;    
};


template <class T>
int charGetNumber(const char *F, int L, const char *s, T &x) {

    float fl=0;
    int i;

    if (!s) {
       sprintf(str,"cannot read number from empty string");
       return 1;
    }

    i=sscanf (s,"%f",&fl);
    x=(T)fl;

    if (!i || fl!=(float)x) {
       sprintf(str,"Invalid input for type `%s' (%g, %d)",
       getName(typeid(T)).data, fl, i); return 1;
    }


    return 0;
}


namespace Wb {

template<>
int GetEnv (const char *F, int L, const char *name, double &val) {

    char *s, *s2;
    double dbl;

    s=getenv(name);    if (!s ) return -1;
    dbl=strtod(s,&s2); if (!s2) return -2;
    val=dbl;

    if (F) wblog(F,L," *  getenv %-6s = %g", name, val);

    if (*s2 && isspace(*s2)) {
       wblog(FL,"WRN %s() got trailing white space (%s) !?",FCT,s);
       ++s2; while (*s2 && isspace(*s2)) ++s2;
    }


    return (*s2 ? -3 : 0);
};

template<>
int GetEnv (const char *F, int L, const char *name, int &val) {

    double dbl=0;
    int i=GetEnv(F,L,name,dbl); if (i) return i;

    i=int(dbl); if (double(i)!=dbl) return -2;
    val=i;

    if (F) wblog(F,L," *  getenv %-6s = %d", name, val);
    return 0;
}

template<>
int GetEnv (const char *F, int L, const char *name, unsigned &val) {

    double dbl=0;
    int i=GetEnv(F,L,name,dbl); if (i) return i;

    i=unsigned(dbl); if (double(i)!=dbl) return -2;
    val=i;

    if (F) wblog(F,L," *  getenv %-6s = %d", name, val);
    return 0;
};

template<>
int GetEnv (const char *F, int L, const char *name, char &val) {

    char *s=getenv(name); if (!s || strlen(s)!=1) return -1;
    val=s[0];

    if (F) wblog(F,L," *  getenv %-6s = '%c'(%d)", name, val,val);
    return 0;
}

template<>
int GetEnv (const char *F, int L, const char *name, wbstring &val) {

    char *s;

    s=getenv(name); if (!s) return -1;
    val=s;

    if (F) wblog(F,L," *  getenv %-6s = '%s'", name, val.data);

    return 0;
}

int strrep(
   const char *S0,
   const char *t0,
   const char *t2,
   char *Sout,
   size_t N,
   char gflag='g'
);

int strrep(
   const char *S0,
   const char *t0,
   const char *t2,
   char *Sout,
   size_t N,
   char gflag
){
   if (!S0 || !S0[0] || !t0 || !t0[0]) return 0;

   int nrep=0; 
   size_t k, l=0, m0=strlen(t0), m=(t2 ? strlen(t2):0);
   const char *s, *s0=S0; char *s2=Sout;

   while ((s=strstr(s0,t0))) { k=s-s0; ++nrep;
      if (k) {
         if (l+k>N) wblog(FL,
           "ERR string out of bounds (%d+%d / %d)",l,k,N); 
         strncpy(s2,s0,k); l+=k; s2+=k;
      }
      if (m) {
         if (l+m>N) wblog(FL,
           "ERR string out of bounds (%d+%d / %d)",l,m,N); 
         strncpy(s2,t2,m); l+=m; s2+=m;
      }
      s0+=(k+m0); if (!gflag) break;
   }

   k=strlen(s0);
   if (k) {
      if (l+k>N) wblog(FL,
        "ERR string out of bounds (%d+%d / %d)",l,k,N); 
      strncpy(s2,s0,k); l+=k; s2+=k;
   }
   s2[0]=0;

   return nrep;
};

wbstring repHome(const char *file) {
   
   size_t n=2*strlen(file);
   const char *s; char F0[n+1], F[n+1];
   char *f0=F0, *f=F; strcpy(f,file);

   if ((s=getenv("MEX" ))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$MEX", f,n)) return f; }
   if ((s=getenv("HOME"))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$HOME",f,n)) return f; }
   if ((s=getenv("LMA" ))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$LMA", f,n)) return f; }
   if ((s=getenv("USER"))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$USER",f,n)) return f; }

   return f;
};

};


wbstring getName(const std::type_info &type_id) {

   wbstring s;

   if (type_id==typeid(double       )) s="double";       else
   if (type_id==typeid(float        )) s="float";        else
   if (type_id==typeid(int          )) s="int";          else
   if (type_id==typeid(long         )) s="long";         else
   if (type_id==typeid(unsigned     )) s="unsigned";     else
   if (type_id==typeid(unsigned long)) s="ulong";        else
   if (type_id==typeid(char         )) s="char";         else
   if (type_id==typeid(size_t       )) s="size_t";       else
                                                     
   if (type_id==typeid(double*      )) s="double*";      else
   if (type_id==typeid(float*       )) s="float*";       else
   if (type_id==typeid(int*         )) s="int*";         else
   if (type_id==typeid(long*        )) s="int*";         else
   if (type_id==typeid(unsigned*    )) s="unsigned*";    else
   if (type_id==typeid(unsigned long*)) s="ulong*";      else
   if (type_id==typeid(char*        )) s="char*";        else
   if (type_id==typeid(size_t*      )) s="size_t*";      else
                                                       
   if (type_id==typeid(wbcomplex    )) s="wbcomplex";    else
   if (type_id==typeid(wbcomplex*   )) s="wbcomplex*";   else

   if (type_id==typeid(long double  )) s="long double";  else
   if (type_id==typeid(long double* )) s="long double*"; else

#ifdef __WB_MPFR_HH__
   if (type_id==typeid(Wb::quad )) s="Wb::quad" ; else
   if (type_id==typeid(Wb::qquad)) s="Wb::qquad"; else
   if (type_id==typeid(Wb::quad2)) s="Wb::quad2"; else
   if (type_id==typeid(Wb::quad3)) s="Wb::quad3"; else
#endif

   s=(char*)(type_id.name());

   return s;
}


bool isBaseType(const std::type_info &type_id) {

   if ( type_id==typeid(unsigned)
     || type_id==typeid(double)
     || type_id==typeid(float)
     || type_id==typeid(int)
     || type_id==typeid(char)
     || type_id==typeid(unsigned char)
     || type_id==typeid(long long)
     || type_id==typeid(size_t)
     || type_id==typeid(wbcomplex)
   ) return 1;

   return 0;
};

bool isComplexType(const std::type_info &type_id) {

   if ( type_id==typeid(wbcomplex)
   ) return 1;

   return 0;
};


template<class T>
void num2Fmt<T>::get_fmt(const char *F, int L, char *s) {

   if (!s || check_paras()) { wblog(F_L,
      "%N%NERR num2Fmt<> invalid type %d.%d%s (having `%s', %lX)%N%N",
      m,p,t,getName(typeid(T)).data,s);
   }

   s[0]='%'; ++s;
   if (m>0) {
      if (p>=0)
           sprintf(s,"%d.%d%s",m,p,t);
      else sprintf(s,"%d%s",m,t);
   }
   else {
     if (p>=0)
          sprintf(s,".%d%s",p,t);
     else sprintf(s,"%s",t);
   }
};


template<>
void num2Fmt<wbcomplex>::get_fmt(const char *F, int L, char *s) {

   if (m<0) m=6;
   if (p<0) p=4;
   if (t[0]<0) t[0]='g';

   if (!s || !strchr("gGeE",t[0])) { wblog(F_L,
      "%N%NERR num2Fmt<> invalid type %d.%d%s (having `wbcomplex', %lX)%N%N",
      m,p,t,s);
   }

   for (unsigned i=0; i<2; ++i) {
      s[0]='%'; ++s; if (i) { s[0]='+'; ++s; }
      if (m>0) {
         if (p>=0)
              sprintf(s,"%d.%d%s",m,p,t);
         else sprintf(s,"%d%s",m,t);
      }
      else {
        if (p>=0)
             sprintf(s,".%d%s",p,t);
        else sprintf(s,"%s",t);
      }
      s+=strlen(s);
   }
   s[0]='i';
};


template<class T>
int num2Fmt<T>::check_paras() {
   wblog(FL,"ERR num2Fmt<> invalid data type `%s'",
   getName(typeid(T)).data); return 1;
};


template<> inline
int num2Fmt<double>::check_paras() { 
   if (m<-1) m=8;
   if (p<-1) p=5;
   if (t[0]<0) t[0]='g'; else
   if (!Wb::strchri("ge",t[0])) return 't';
   return 0;
};

template<> inline
int num2Fmt<float>::check_paras() { 
   return ((num2Fmt<double>*)this)->check_paras();
};

template<> inline
int num2Fmt<int>::check_paras() {
   if (m<-1) m=4;
   if (p>=0) return 'p';
   if (t[0]<0) t[0]='d'; else
   if (!strchr("ducxX",t[0])) return 't';
   return 0;
};

template<> inline
int num2Fmt<unsigned>::check_paras() {
   return ((num2Fmt<int>*)this)->check_paras();
};

template<> inline
int num2Fmt<long>::check_paras() {
   if (m<-1) m=4;
   if (p>=0) return 'p';
   if (t[0]<0) { t[0]='l'; t[1]='d'; } else
   if (!strchr("lducxX",t[0])) return 't';
   return 0;
};

template<> inline
int num2Fmt<unsigned long>::check_paras() {
   return ((num2Fmt<long>*)this)->check_paras();
};

template<> inline
int num2Fmt<char>::check_paras() {
   if (m<-1) m=4;
   if (p>=0) return 'p';

   if (t[0]<0) t[0]='d';
   else if (!strchr("ducxX",t[0])) return 't';
   return 0;
};

template<> inline
int num2Fmt<char*>::check_paras() {
   if (t[0]<0) t[0]='s'; else
   if (t[0]!='s') return 't';
   return 0;
};


template <class T> inline
char* defaultFmt(char *fmt, 
const T &x __attribute__ ((unused)), int m, int p, char t) {

   const size_t n=16; size_t l=0;
   char dd[n]; dd[0]=0; fmt[0]=0;

   if (m>0) {
        if (p>=0)
             l=snprintf(dd,n,"%d.%d",m,p);
        else l=snprintf(dd,n,"%d",m);
   }
   else if (p>=0) l=snprintf(dd,n,".%d",p);

   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   if (typeid(T)==typeid(double) || typeid(T)==typeid(float)) {
      if (!t) t='g';
      sprintf(fmt,"%%%s%c",dd,t);
   }
   else if (typeid(T)==typeid(wbcomplex)) {
      if (!t) t='g';
      sprintf(fmt,"%%%s%c%%+%s%ci",dd,t,dd,t);
   }
   else if (
      typeid(T)==typeid(int) || typeid(T)==typeid(char) ||
      typeid(T)==typeid(unsigned)
    ){
      if (p>=0) wblog(FL,"WRN %s() got format *.%d* for type '%s'",
          FCT, p, getName(typeid(T)).data);
      if (!t) t='d';
      sprintf(fmt,"%%%s%c",dd,t);
   }
   else if (typeid(T)==typeid(char*)) {
       if (t && t!='s') wblog(FL,
          "ERR %s() got format type '%c'<%d> for type '%s'",
           FCT,t,t,getName(typeid(T)).data
       );
       sprintf(fmt,"%%%ss",dd);
   }
   else wblog(FL,
     "ERR %s() unsupported type '%s'",FCT,getName(typeid(T)).data);

   return fmt;
};



wbstring wbTimeStamp() {

   time_t     t  = time(NULL);
   struct tm *tb = localtime(&t);
   wbstring   s  = asctime(tb);

   unsigned n=strlen(s.data);
   if (n && s.data[n-1]<32) s.data[n-1]=0;

   return s;
}


template <class T1, class T2> inline 
void safeConvert(const char *F, int L, const T1 &x1, T2 &x2) {
   x2=T2(x1); if (T1(x2)!=x1) wblog(F_L,
      "ERR converting %s -> %s changes value\n%g, %g",
       getName(typeid(T1)).data, getName(typeid(T2)).data,
       double(x1), double(x2)
   );
}

template <> inline 
void safeConvert(const char *F, int L, const wbcomplex &x1, double &x2) {
   if (x1.i!=0.) wblog(F_L,
      "ERR cannot cast complex number to real (%.4g%+.4gi)",x1.r,x1.i);
   x2=x1.r;
}

template <> inline
void safeConvert(
   const char *F __attribute__ ((unused)), int L __attribute__ ((unused)),
   const double &x1, wbcomplex &x2
){ x2=wbcomplex(x1,0); }

   namespace Wb {

template <class T> inline 
bool isLower(const T* a, const T* b, const size_t n) {
   for (size_t i=0; i<n; i++) {
       if (a[i]<b[i]) return 1; else
       if (a[i]>b[i]) return 0;
   }
   return 0;
}

template <class T> inline 
bool isEqual(const T* a, const T* b, const size_t n) {

   for (size_t i=0; i<n; i++) { if (a[i]!=b[i]) return 0; }
   return 1;
}


template <class T> inline 
bool allEqual(const T* a, const size_t n, const T &x) {

   for (size_t i=0; i<n; i++) { if (a[i]!=x) return 0; }
   return 1;
}

template <class T> inline 
bool anyUnequal(const T* a, const size_t n, const T &x) {
   for (size_t i=0; i<n; i++) { if (a[i]!=x) return 1; }
   return 0;
};

template <class T> inline 
bool anyEqual(const T* a, const size_t n, const T &x) {

   for (size_t i=0; i<n; i++) { if (a[i]==x) return 1; }
   return 0;
};


template <class T> inline
void scale_eps(T &eps, const T* d, size_t n) {
   if (double(eps)>0) {
      if (double(eps)>1E-8) wblog(FL,
         "WRN %s() got eps=%g (ignore)",FCT,double(eps));
      else {
         T x=getdscale(d,n);
         if (x>T(1E-8)) eps*=x;"real" data is present
      }
   }
};

template <class T>
T getdscale(const T* data, size_t n) {

   if (n) { T d,x=0;
      for (size_t i=0; i<n; ++i) { d=ABS(data[i]); if (x<d) x=d; }
      return x;
   }
   return 0;
};

 }

template <class T> inline
T addRange2(const T* d, const size_t n) {
   T s=0; for (size_t i=0; i<n; ++i) s+=ABS2(d[i]);
   return s;
};

template <class T> inline
T addRange(const T* d, const size_t n) {
   T s=0; for (size_t i=0; i<n; ++i) s+=d[i];
   return s;
};

template <class T> inline
void addRange(const T* a, T* c, const size_t n) {
   for (size_t i=0; i<n; ++i) { c[i]+=a[i]; }
};

template <class T> inline
void addRange(
   const T* a, const T* b, T* c, size_t n, size_t stride
){
   size_t i=0;
   if (stride==1) { for (; i<n; ++i) c[i]=a[i]+b[i]; }
   else {
      if (!stride && n) wblog(FL,
         "ERR %s() got stride %d/%d",FCT,stride,n);
      for (n*=stride; i<n; i+=stride) c[i]=a[i]+b[i];
   }
};

template <class T> inline
void diffRange(
   const T* a, const T* b, T* c, size_t n, size_t stride
){
   size_t i=0;
   if (stride==1) { for (; i<n; ++i) c[i]=a[i]-b[i]; }
   else {
      if (!stride && n) wblog(FL,
         "ERR %s() got stride %d/%d",FCT,stride,n);
      for (n*=stride; i<n; i+=stride) c[i]=a[i]-b[i];
   }
};

template <class T> inline
void setRange2avg(T* a, size_t n) {
   T x=0; size_t i;
   for (i=0; i<n; i++) x+=a[i]; x*=(T(1)/n);
   for (i=0; i<n; i++) a[i]=x;
}

template <class T> inline
size_t nnzRange(const T* a, size_t n) {
   size_t i=0, m=0;
   for (; i<n; i++) if (a[i]!=0) m++;
   return m;
}


template <class T> inline
T rangeNormDiff2(
   const T* d1, const T* d2, size_t n
){
   if (!n || d1==d2) return 0;

   T x2, x2sum=0;

   for (size_t i=0; i<n; ++i) {
      x2=d1[i]-d2[i]; x2*=CONJ(x2);
      x2sum+=x2;
   }

   return x2sum;
};


template <class T> inline
T rangeNormDiff2(
   const T* d1, const T* d2, size_t n, const T &fac, size_t *k
){
   if (fac==1 && !k) { return rangeNormDiff2(d1,d2,n); }
   if (!n || (d1==d2 && fac==1)) { if (k) { (*k)=0; }; return 0; }

   T x2, x2sum=0;

   if (k) { T x2max=0;
      for (size_t i=0; i<n; ++i) {
         x2=d1[i]-fac*d2[i]; x2*=CONJ(x2);
         x2sum+=x2; if (x2>x2max) { x2max=x2; (*k)=i; }
      }
   }
   else {
      for (size_t i=0; i<n; ++i) {
         x2=d1[i]-fac*d2[i]; x2*=CONJ(x2);
         x2sum+=x2;
      }
   }

   return x2sum;
};


template <class T> inline
T rangeNorm2(const T* d, size_t n, size_t *k
){
   T x2sum=0;
   if (!n) { if (k) { (*k)=0; }; return x2sum; }

   if (k) {
      T x2, x2max=0; (*k)=0;

      for (size_t i=0; i<n; ++i) {
         x2=CONJ(d[i])*d[i]; if (x2>x2max) { x2max=x2; (*k)=i; }
         x2sum+=x2;
      }
   }
   else {
      for (size_t i=0; i<n; ++i) { x2sum+=CONJ(d[i])*d[i]; }
   }

   return x2sum;
};


template <class T> inline
T rangeMaxDiff(const T* d1, const T* d2, size_t n, size_t *k) {
   size_t imax=0, i=0; T x,xmax=0;
   
   if (!n || d1==d2) { if (k) (*k)=0; 
      if (!n) wblog(FL,"WRN %s() got empty input",FCT);
      return 0;
   }

   for (; i<n; i++) {
      x=ABS(d1[i]-d2[i]); if (x>xmax) { xmax=x; imax=i; }
   }

   if (k) (*k)=imax;
   return xmax;
};


template <class T, class T2> inline
void Wb::cpyStride(
   T2* d,
   const T* d0,
   unsigned n,
   const unsigned *idx, unsigned m,
   unsigned D2,
   unsigned D0,
   char add_flag
){
   unsigned i,k; if (!n || !m) return;

   if (int(D0)<0) D0=n;
   if (int(D2)<0) D2=n;
   
   if (D0<n || D2<n) wblog(FL,
      "ERR %s() stride too small (%d,%d/%d)",FCT,D0,D2,n);

   if (add_flag) {
      if (idx) { const T *d0_=d0;
         for (k=0; k<m; ++k, d+=D2) { d0 = d0_+ idx[k] * D0;
         for (i=0; i<n; ++i) d[i]+=d0[i];  }
      }
      else {
         for (k=0; k<m; ++k, d+=D2, d0+=D0)
         for (i=0; i<n; ++i) d[i]+=d0[i];
      }
   }
   else {
      if (idx) { const T *d0_=d0;
         for (k=0; k<m; ++k, d+=D2) { d0 = d0_+ idx[k] * D0;
         for (i=0; i<n; ++i) d[i]=d0[i]; }
      }
      else {
         for (k=0; k<m; ++k, d+=D2, d0+=D0)
         for (i=0; i<n; ++i) d[i]=d0[i];
      }
   }
};


template <class T> inline
void Wb::cpyStride(
   T* d, const T* d0, unsigned n,
   const unsigned *idx, unsigned m,
   unsigned D2, unsigned D0,
   char add_flag
){
   unsigned i,k; if (!n || !m) return;

   if (int(D0)<0) D0=n;
   if (int(D2)<0) D2=n;

   if (D0<n || D2<n) wblog(FL,
      "ERR %s() stride too small (%d,%d/%d)",FCT,D0,D2,n);
   if (gotMemOverlap(d0, (m-1)*D0+n, d, (m-1)*D2+n)) wblog(FL,
      "ERR %s() must not copy onto itself",FCT);

   if (add_flag) {
      if (idx) { const T *d0_=d0;
         for (k=0; k<m; k++, d+=D2) { d0 = d0_+ idx[k] * D0;
         for (i=0; i<n; i++) d[i]+=d0[i]; }
      }
      else {
         for (k=0; k<m; k++, d+=D2, d0+=D0)
         for (i=0; i<n; i++) d[i]+=d0[i];
      }
   }
   else {
      if (idx) { const T *d0_=d0;
         for (i=0; i<m; i++, d+=D2) { d0 = d0_+ idx[k] * D0;
         MEM_CPY<T>(d,n,d0); }
      }
      else {
         for (i=0; i<m; i++, d+=D2, d0+=D0)
         MEM_CPY<T>(d,n,d0);
      }
   }
};


template <class T, class T2> inline 
void cpyRange(const T* a, T2* c, size_t n, char check_type) {
   if (check_type) {
      for (size_t i=0; i<n; ++i) { c[i]=T2(a[i]);
         if (T(c[i])!=a[i]) wblog(FL,"ERR %s() got rounding error !?? "
            "(%g, %g)",FCT,double(a[i]),double(c[i])
         );
      }
   }
   else { for (size_t i=0; i<n; ++i) c[i]=T2(a[i]); }
};


template <class T> inline 
void cpyRange(const T* a, T* c, size_t n) {
   for (size_t i=0; i<n; ++i) c[i]=a[i];
}




template <class T> inline 
void cpyRangeR(const wbcomplex* a, T* c, size_t n) {
   for (size_t i=0; i<n; i++) c[i]=T(a[i].r);
}

template <class T> inline 
void cpyRangeI(const wbcomplex* a, T* c, size_t n) {
   for (size_t i=0; i<n; i++) c[i]=T(a[i].i);
}

template <class T> inline 
void cpyRangeA(const wbcomplex* a, T* c, size_t n) {
   for (size_t i=0; i<n; i++) c[i]=a[i].abs();
}


template <class T> inline 
bool uniformRange(const T* d, size_t n) {
   for (size_t i=1; i<n; ++i) if (d[i]!=d[0]) return 0;
   return 1;
};


template <class T> inline 
T maxRange(const T* d, size_t n) {

   if (!n) wblog(FL,"ERR %s() got empty range",FCT);

   T x=d[0];
   for (size_t i=1; i<n; ++i) { if (x<d[i]) { x=d[i]; }}
   return x;
};

template <class T> inline 
T minRange(const T* d, size_t n) {

   if (!n) wblog(FL,"ERR %s() got empty range",FCT);

   T x=d[0];
   for (size_t i=1; i<n; ++i) { if (x>d[i]) { x=d[i]; }}
   return x;
};


template <class T> inline
T prodRange(const T* d, size_t n) {
   T x=0;

   if (n) { size_t i=1;
      if (!d) wblog(FL,"ERR %s() got null space (n=%d)",FCT,n);
      for (x=d[0]; i<n; ++i) x*=d[i];
   }
   else {
      wblog(FL,"WRN %s() got empty range (returning 0)",FCT);
   }

   return x;
};

template <class T> inline 
T prodRange(const T* d, size_t n, T x) {
   if (n) { size_t i=1;
      if (!d) wblog(FL,"ERR %s() got null space (n=%d)",FCT,n);
      for (x=d[0]; i<n; ++i) x*=d[i];
   }
   return x;
};

template <class T> inline 
void timesRange(T* d, T x, size_t n, size_t stride) {
   if (!n) return;
   if (!d) wblog(FL,"ERR %s() got null space (n=%d)",FCT,n);

   if (stride==1) {
      for (size_t i=0; i<n; ++i) d[i]*=x;
   }
   else {
      if (!stride && n) wblog(FL,
         "ERR %s() got stride %d/%d",FCT,stride,n);
      for (size_t l=0, i=0; i<n; ++i, l+=stride) d[l]*=x;
   }
};


template <class T> inline
T prodRange(const T* d1, const T* d2, size_t n, T x) {
   if (n) {
      if (!d1 || !d2) wblog(FL,
         "ERR %s() got null space (%lX, %lX)",FCT,d1,d2);
      for (size_t i=0; i<n; ++i) x+=(d1[i]*d2[i]);
   }
   return x;
};

template <class T> inline
T prodRange(const T* d1, const T* d2, size_t n) {
   T x=0;

   if (n) {
      if (!d1 || !d2) wblog(FL,
         "ERR %s() got null space (%lX, %lX)",FCT,d1,d2);
      for (size_t i=0; i<n; ++i) x+=(d1[i]*d2[i]);
   }
   else wblog(FL,"WRN %s() got empty range",FCT);

   return x;
};


template<class T> inline
T overlap(
   const T* v1, const T* v2, size_t n,
   size_t stride,
   char tnorm __attribute__ ((unused))
){
   T x2=0;
   if (n) {
      if (!v1 || !v2) wblog(FL,
         "ERR %s() got null space (%lX, %lX)",FCT,v1,v2);

      if (stride==1) {
         for (size_t i=0; i<n; ++i) { x2+=v1[i]*v2[i]; }
      }
      else {
         if (!stride) wblog(FL,
            "ERR %s() got stride %d/%d",FCT,stride,n);
         for (size_t l=0, i=0; i<n; ++i, l+=stride) {
            x2+=CONJ(v1[l])*v2[l];
         }
      }
   }
   return x2;
};


template<> inline
wbcomplex overlap(
   const wbcomplex* v1, const wbcomplex* v2, size_t n,
   size_t stride, char tnorm
){
   if (tnorm) {
      wbcomplex x2=0;
      if (stride==1) {
         for (size_t i=0; i<n; ++i) { x2+=v1[i]*v2[i]; }
      }
      else {
         if (!stride && n) wblog(FL,
            "ERR %s() got stride %d/%d",FCT,stride,n);
         for (size_t l=0, i=0; i<n; ++i, l+=stride) { x2+=v1[l]*v2[l]; }
      }
      return x2;
   }
   else if (v1!=v2) {
      wbcomplex z2=0;
      if (stride==1) {
         for (size_t i=0; i<n; ++i) z2+=v1[i].conj()*v2[i];
      }
      else {
         if (!stride && n) wblog(FL,
            "ERR %s() got stride %d/%d",FCT,stride,n);
         for (size_t l=0, i=0; i<n; ++i, l+=stride) z2+=v1[l].conj()*v2[l];
      }
      return z2;
   }
   else {
      double nrm2=0;
      if (stride==1) {
         for (size_t i=0; i<n; ++i) { nrm2+=v1[i].abs2(); }
      }
      else {
         if (!stride && n) wblog(FL,
            "ERR %s() got stride %d/%d",FCT,stride,n);
         for (size_t l=0, i=0; i<n; ++i, l+=stride) { nrm2+=v1[l].abs2(); }
      }
      return nrm2;
   }
};



template<class T> 
T gs_project_range(
   T *v1, const T *v2,
   size_t n, size_t stride, char isnorm, char tnorm
){
   T x=overlap(v2,v1,n,stride,tnorm);

   if (!isnorm) {
      T nrm2=overlap(v2,v2,n,stride,tnorm);
      if (fabs(double(nrm2))<1E-12) wblog(FL,
         "WRN %s() got norm=%.3g !??",FCT,double(nrm2));
      x/=nrm2;
   }

   if (stride==1)
      for (size_t i=0; i<n; ++i) v1[i]-=x*v2[i];
   else {
      if (!stride && n) wblog(FL,
         "ERR %s() got stride %d/%d",FCT,stride,n);
      for (size_t l=0, i=0; i<n; ++i, l+=stride) v1[l]-=x*v2[l];
   }
   return x;
};


template <class T> inline
size_t replRange(T* a, size_t n, T x, T v) {

   size_t i,m=0;

   if (!isnan(x))
        for (i=0; i<n; i++) if (a[i]==x)     { a[i]=v; m++; }
   else for (i=0; i<n; i++) if (isnan(a[i])) { a[i]=v; m++; }

   return m;
}

template <class T> inline
size_t replRange(
   const char *F, int L, T *a, size_t n, T x, T v
){
   size_t m=replRange(a,n,x,v);
   if (m) wblog(F,L,"Skipped %d NaN's",m);
   return m;
}

inline size_t countNaN(double* a, size_t m) {
   size_t i,n=0;
   for (i=0; i<m; i++) if (std::isnan(a[i])) ++n;
   return n;
}


   template <class T>
   void cpyZRange(
      const double *R, const double *I, T *Z, size_t n
   ){
      wblog(FL,"ERR %s not applicable for type '%s'\n(%lX,%lX,%lX,%d)",
      FCT, getName(typeid(T)).data,R,I,Z,n);
   };


   template <>
   void cpyZRange(
      const double *R, const double *I, wbcomplex *Z, size_t n
   ){
      for (size_t i=0; i<n; i++)
      Z[i].set(R[i], I ? I[i] : 0);
   }


   template <class T>
   void splitZRange(
      const T *Z, double *R, double *I, size_t n
   ){
      wblog(FL,"ERR %s not applicable for type `%s'",
      FCT, getName(typeid(T)).data);
   }


   template <>
   void splitZRange(
      const wbcomplex *Z, double *R, double *I, size_t n
   ){
      for (size_t i=0; i<n; i++) {
      R[i]=Z[i].r; if (I) I[i]=Z[i].i; }
   }


namespace Wb {

template<class T>
void chopTiny_float(T *d, size_t n, T dref) { return; };

template<>
void chopTiny_float(double *d, size_t n, double ref) {

   size_t i; if (!n) return;

   if (ref<0) {
      for (ref=0, i=0; i<n; i++) ref=MAX(ref,::fabs(d[i]));
      if (ref==0) return;
   }

   if (ref!=0) {
      int k=int(std::log(ref)/std::log(2));
      if (k>=0) ref=(1<<k); else ref=1/(1<<-k);

      for (i=0; i<n; i++) {
         if (d[i]>0.) d[i]=float(d[i]+ref)-ref; else
         if (d[i]<0.) d[i]=float(d[i]-ref)+ref;
      }
   }
   else {
      for (i=0; i<n; i++)
      if (d[i]!=0.) d[i]=double(float(d[i]));
   }
};

};


void invertIndex(
    const wbvector<unsigned> &i1, unsigned N,
    wbvector<unsigned> &i2,
    char lflag
){
   wbvector<char> flag(N);
   unsigned i,k,n;

   for (i=0; i<i1.len; i++) {
       k=i1[i]; if (k<N) flag[k]=1; else
       wblog(FL, "ERR index out of bounds (%d/%d)", k, N);
   }
   n=flag.isum();

   if (n!=i1.len && lflag)
   wblog(FL, "WRN index set not unique!? (%d/%d)", n, i1.len);

   i2.init(N-n);
   for (k=i=0; i<N; i++) if (flag[i]==0) i2[k++]=i;
}


void setRand(wbvector<wbcomplex> &zz, double fac, double shift) {
   size_t i=0; fac/=(double)RAND_MAX;

   if (shift==0.)
   for (; i<zz.len; ++i)
        zz.data[i].set(fac*rand(), fac*rand());
   else
   for (; i<zz.len; ++i)
        zz.data[i].set(fac*rand()+shift,fac*rand()+shift);
}


template<class T>
void set2avg(T *dd, const size_t n) {
   size_t i;  T dbl=0; if (!n) return;
   for (i=0; i<n; i++) dbl+=dd[i]; dbl/= n;
   for (i=0; i<n; i++) dd[i]=dbl;
}

template<class T>
void getDiff(const wbvector<T> &xx, wbvector<T> &dx) {
   size_t i, n=xx.len; dx.init(n);
   for (i=1; i<n; i++) {
       dx[i-1] += (dx[i] = xx[i]-xx[i-1]);
       if (i>1) dx[i-1] *= 0.5;
   }
}


template<class T>
T wbIntTrapez(const T *xx, const T *yy, const size_t n) {
   T s=0;

   for (size_t i=1; i<n; i++)
   s += 0.5*(xx[i]-xx[i-1])*(yy[i]+yy[i-1]);

   return s;
}



template<class T>
void markSet(
   wbvector< wbvector<T>* > &E0,
   wbvector< wbvector<char> > &mark,
   int &Nkeep,
   double Etrunc,
   wbvector<T> &E,
   double eps,
   double b,
   int dmax,
   const char *dir"asc" (lowest Nkeep), "desc" (largest Nkeep)
){
   unsigned i,d,l,N,n=E0.len;
   wbvector<char> mm;
   wbperm P,iP;

   for (N=i=0; i<n; i++) N+=(E0[i]->len);

   E.init(N); mark.initDef(n);
   if (!N) return;

   for (l=i=0; i<n; i++, l+=d) {
      d=(E0[i]->len);  mark[i].init(d);
      memcpy(E.data+l, E0[i]->data, d*sizeof(T));
   }

   E.Sort(P); P.getIPerm(iP);

   if (eps>0) {
      unsigned i0=0;

      eps *= ABS(E.end()-E[0]);
      for (i=1; i<=N; i++) {
         if (i==N || ABS(E[i]-E[i-1])>eps) {
            d=i-i0; if (d>1)
            setRange2avg(E.data+i0,d);
            i0=i;
         }
      }

      wbvector<T> x; E.permute(x,iP);
      for (l=i=0; i<n; i++, l+=d) {
         d=E0[i]->len;
         memcpy(E0[i]->data, x.data+l, d*sizeof(T));
      }
   }

   if (Nkeep<0) Nkeep=N;
   if (Nkeep==0) return;
   if (Nkeep>=(int)N && Etrunc<=0) {
      for (i=0; i<n; i++) mark[i].set(1);
      Nkeep=N; return;
   }

   markSet(FL,E,mm,Nkeep,Etrunc,b,dmax,dir);

   mm.Permute(iP);

   for (l=i=0; i<n; i++, l+=d) {
      d=mark[i].len;
      memcpy(mark[i].data, mm.data+l, d*sizeof(char));
   }
};


template<>
void markSet(
   wbvector< wbvector<wbcomplex>* > &E0,
   wbvector< wbvector<char> > &mark,
   int &Nkeep,
   double Etrunc,
   wbvector<wbcomplex> &E,
   double eps,
   double b,
   int dmax,
   const char *dir_"asc" (lowest Nkeep), "desc" (largest Nkeep)
){
   unsigned i,d,l,N, n=E0.len;
   const char *s=NULL;

   wbvector<char> mm;
   wbvector<double> ee;
   wbperm P,iP;


   if (dir_ && dir_[0]) {
      if (!strncmp(dir_,"asc", 3)) { s=dir_+3; } else
      if (!strncmp(dir_,"desc",4)) { s=dir_+4; } else
      wblog(FL,"ERR %s() invalid direction >%s<",FCT,dir_);
   }

   for (N=i=0; i<n; i++) N+=(E0[i]->len);
   if (!N) return;

   E.init(N); ee.init(N); mark.initDef(n);

   for (l=i=0; i<n; i++, l+=d) {
      d=(E0[i]->len);  mark[i].init(d);
      memcpy(E.data+l, E0[i]->data, d*sizeof(wbcomplex));
   }

   if (s==NULL || !strcmp(s,"R")) {
      cpyRangeR(E.data,ee.data,N);
   }
   else if (!strcmp(s,"I")) {
      cpyRangeI(E.data,ee.data,N);
   }
   else if (!strcmp(s,"A")) {
      wbcomplex z=E.min();
      wbvector<wbcomplex> x(E); z.i=0; x-=z;
      cpyRangeA(x.data,ee.data,N);
   }
   else wblog(FL,"ERR markSet() invalid direction >%s<", dir_);

   ee.Sort(P); E.Select(P); P.getIPerm(iP);

   if (Nkeep<0) Nkeep=N;
   if (Nkeep==0) return;

   if (Nkeep>=(int)N && Etrunc<=0) {
      for (i=0; i<n; i++) mark[i].set(1);
      Nkeep=N; return;
   }

   markSet(FL,ee,mm,Nkeep,Etrunc,MAX(b,eps),dmax,dir_);


   mm.Permute(iP);

   for (l=i=0; i<n; i++, l+=d) {
      d=mark[i].len;
      memcpy(mark[i].data, mm.data+l, d*sizeof(char));
   }
}


template<class T>
void markSet(const char *F, int L,
   const wbvector<T> &E,
   wbvector<char> &mark,
   int &Nkeep,
   double Etrunc,
   double b,
   int dmax,
   const char *dir_"asc" (lowest Nkeep), "desc" (largest Nkeep)
){


   unsigned i,imax,l;
   unsigned N=E.len, nx;
   double x, xmax=0;
   char dir=+1;

   if (dir_ && dir_[0]) {
      if (!strncmp(dir_,"asc", 3)) dir=+1; else
      if (!strncmp(dir_,"desc",4)) dir=-1; else
      wblog(F,L,"ERR %s() invalid direction >%s<",FCT,dir_);
   }

   mark.init(N);

   if (Nkeep<0 || Nkeep>int(N)) Nkeep=N;
   if (Nkeep==0) return;

   if (Nkeep>=(int)N && Etrunc<=0) { mark.set(1); Nkeep=N; return; }
   if (N<1) wblog(FL,"ERR %s() N=%d !??",FCT,N);

   if (dir>0) {
      if (Etrunc>0) {
         x=Etrunc+E[0];
         for (l=1; l<N; ++l) {
            if (E[l-1]>E[l]) wblog(FL,"ERR %s() data not sorted!",FCT);
            if (E[l]>=x) break;
         }
         if (int(l)>Nkeep) l=Nkeep;
      }
      else { l=(unsigned)Nkeep; }

      if (l==N) { mark.set(1); Nkeep=N; return; }
      if (l==0) { Nkeep=0; return; }

      nx=(dmax>=0 ? (unsigned)dmax : MAX(16U,l/5));

      if (b>0) {
         for (imax=i=l; i<N && i<l+nx; ++i) {
            x=fabs(E[i]-E[i-1]);
            if (x<=b) { if (xmax<x) { xmax=x; imax=i; }}
            else break;
         }
         if (i>=l+nx) wblog(F,L,
            "WRN %s() reached Nkeep+dmax=%d+%d states\n"
            "(db=%.3g too small !? => %d @ %.3g)", FCT,l,nx,b,imax,xmax);

         Nkeep=l=imax;
      }
      else if (b<0 && nx && l<N) {

         double q=10/(-b*nx);
         unsigned NX=5/q, i2=MIN(N,l+NX);
         i=l-NX; if (int(i)<1) i=1;

         for (imax=l; i<i2; ++i) { x=q*int(i-l);
            x=exp(-x*x)*fabs(E[i]-E[i-1]);
            if (x>xmax) { xmax=x; imax=i; }
         }
         if (i==N) { x=q*int(i-l);
            x=0.1*exp(-x*x)*fabs(E[N-1]-E[0]);
            if (x>xmax) { xmax=x; imax=i; }
         }

         Nkeep=l=imax;
      }

      for (i=0; i<l; ++i) mark[i]=1;
   }
   else {
      if (Etrunc>0) {
         double x=Etrunc;
         for (l=N-1; l<N; --l) { 
            if (l && E[l-1]>E[l]) wblog(FL,"ERR %s() data not sorted!",FCT);
            if (E[l]<=x) break;
         }
         if (int(N-1-l)>Nkeep) { l = N-1-(unsigned)Nkeep; }
      }
      else { l=N-1-(unsigned)Nkeep; }

      if (l>N ) { mark.set(1); Nkeep=N; return; }
      if (l==N-1) { Nkeep=0; return; }

      nx=(dmax>=0 ? (unsigned)dmax : MAX(16U,(N-l)/5));

      if (b>0) {
         for (imax=i=l; i<N && (l-i)<nx; --i) {
            x=fabs(E[i+1]-E[i]);
            if (x<=b) { if (xmax<x) { xmax=x; imax=i; }}
            else break;
         }
         if ((l-i)>=nx) wblog(F,L,
            "WRN %s() reached Nkeep+dmax=%d+%d states (%d/%d)\n"
            "(db=%.3g too small !? => %d @ %.3g)",FCT,
            Nkeep,nx, Nkeep+nx,N, b,N-imax,xmax);
         l=imax; Nkeep=N-1-l;
      }
      else if (b<0 && nx && l<N ) {

         double q=10/(-b*nx);
         unsigned NX=5/q, i1=l-NX; if (int(i1)<0) i1=0;
         i=MIN(N-1,l+NX);

         for (imax=l; i>i1; --i) { x=q*int(i-l);
            x=exp(-x*x)*fabs(E[i+1]-E[i]);
            if (xmax<x) { xmax=x; imax=i; }
         }
         if (i==0) { x=q*int(i-l);
            x=0.1*exp(-x*x)*fabs(E[N-1]-E[0]);
            if (x>xmax) { xmax=x; imax=i-1; }
         }
         l=imax; Nkeep=N-1-l;
      }

      for (i=(l>N  ? 0 : l+1); i<N; ++i) mark[i]=1;
   }
};


#endif

