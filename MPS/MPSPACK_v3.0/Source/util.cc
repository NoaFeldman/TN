#ifndef __WB_UTIL_CC__
#define __WB_UTIL_CC__

/* -------------------------------------------------------------------- */
namespace Wb {
/* -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- *\
   using continued fractions to determine rational approximation
   of given number; see also $MLIB/cfract.m -- Wb,Jul07,10

   NB! if function fails to obtain rational representation within
   given number of iterations and tolerance, 
   then function returns 1. and P and Q remain unchanged
\* -------------------------------------------------------------------- */

inline int Rational(
   double &x0,
   long &P, long &Q,
   double *dx,
   unsigned *nmax,
   wbvector<double>* aa,
   unsigned niter,
   long pqmax,

   double eps,

   double eps2,

   char vflag
){
   unsigned i; int err=0;
   long p[niter+2], q[niter+2];
   double e,xi,r=0, x=x0, a[niter+1];

   if (eps<eps2 || eps2<=0) wblog(FL,
      "ERR %s() invalid eps=%g, eps2=%g",FCT,eps,eps2);

   if (pqmax<=0) { pqmax=long(MIN(1E12, 1/eps)); }
   if (pqmax<=0 || pqmax>(1L<<50)) {
      wblog(FL,"ERR %s() pqmax=%ld !??",FCT,pqmax);
   }

   if (std::fabs(x0)>pqmax || (x!=0 && std::fabs(1/x0)>pqmax)) {
      if (vflag) wblog(FL,"WRN %s() double out of bounds long int",FCT);
      if (dx) { (*dx)=0; }
      if (nmax) { (*nmax)=0; }; if (aa) { aa->init(); }
      return -1;
   }

   p[0]=1; p[1]=a[0]=xi=::round(x);
   q[0]=0; q[1]=1;

   if (vflag) printf(
      "\n  %s() finding rational approximation for %.16g ...\n"
      "\n%6d: %4g   | %8ld / %8ld\n",FCT,x, 0,a[0],p[1],q[1]);

   if (std::fabs(xi-x)<eps) {
      P=p[1]; Q=q[1];
      if (dx) { (*dx)=std::fabs(xi-x); }
      if (nmax) { (*nmax)=0; }; if (aa) { aa->init(1,a); }
      x0=xi; return 0;
   }

   if (std::fabs(x)>1) {
      eps *=std::fabs(x);
   }

   r=std::fabs(xi-x0);

   for (i=0; i<niter; ++i) {
      if ((e=std::fabs(x-a[i]))<eps || r<eps2) break;

      x=1./(x-a[i]); a[i+1]=::round(x);
      p[i+2]=a[i+1]*p[i+1]+p[i];
      q[i+2]=a[i+1]*q[i+1]+q[i];
      
      xi=double(p[i+2])/q[i+2];
      r=std::fabs(xi-x0);

      if (vflag) {
         printf("%6d: %4g   | %8ld / %8ld = %20.16g  @  %8.3g (%.3g)\n",
         i+1, a[i+1], p[i+2],q[i+2], xi,r, x-a[i+1]);
      }
      if (labs(p[i+2])>pqmax || labs(q[i+2])>pqmax) {
         ++i; err=1; break;
      }
   }
   P=p[i+1]; Q=q[i+1]; if (Q<0) { P=-P; Q=-Q; }
   if (i>=niter) err=2;

   if (vflag) {
      if (i<niter) printf("\n  "
         "Rational representation: %ld/%ld\n\n",P,Q);
      else printf("\n  "
         "Failed to find rational representation (%d @ %g)\n\n",
          niter,eps
      );
   }

   if (dx) (*dx)=r;
   if (aa) aa->init(i,a);
   if (nmax) (*nmax)=i;

   if (!err) x0=xi;
   return err;
};


template<class T>
double FixRational(
   const char *F  __attribute__ ((unused)),
   int L          __attribute__ ((unused)), 
   T *d           __attribute__ ((unused)), 
   unsigned n     __attribute__ ((unused)), 
   char rflag     __attribute__ ((unused)),
   unsigned niter __attribute__ ((unused)), 
   long pqmax     __attribute__ ((unused)),
   double eps     __attribute__ ((unused)),
   double eps2    __attribute__ ((unused)),
   double *rz     __attribute__ ((unused)),
   double *ra     __attribute__ ((unused)),
   char vflag     __attribute__ ((unused))
) { 
   return 0;
};

template<>
double FixRational(
   const char *F, int L, double *d, unsigned n,
   char rflag, unsigned niter,
   long pqmax,
   double eps,
   double eps2,
   double *rz_,
   double *ra_,
   char vflag
){
   unsigned i,m, nmax=0, nmissed=0; int e; long p,q;
   double x,a=0,r=0,x2,rz=0,ra=0,r2=0, epsi=eps;

   for (i=0; i<n; ++i) {
      x=std::fabs(d[i]); if (a<x) a=x;

      x=std::fabs(d[i]-::round(d[i])); if (r<x) r=x;
      if (x<eps) {
         d[i]=::round(d[i]); x2=(x*x);
         if (d[i]) r2+=x2; else rz+=x2;
         if (ra<x2) ra=x2;
      }
   }
   if (rz_) { (*rz_)=rz; } else { r2+=rz; }

   if (r<eps || int(niter)<1) {
      if (ra_) { (*ra_)=ra; }
      return r2;
   }

   if (a>1) epsi*=a;

   for (i=0; i<n; i++) {
       e=Rational(d[i],p,q,&r,&m,NULL,niter,pqmax,eps,eps2,vflag>1);
       if (!e) {
          if (nmax<m) nmax=m;
          if (r>epsi) {
             wblog(F_L,"WRN fixing number %.16g = %ld/%ld "
             "(%d; %.3g)",d[i],p,q,m,r);
          }
          x2=r*r; r2+=x2; if (ra<x2) ra=x2;
          continue;
       }
       if (!rflag) { 
          nmissed++; continue;
       }

       x=d[i]*d[i];

       if (x<eps) { nmissed++; continue; }

       e=Rational(x,p,q,&r,&m,NULL,niter,pqmax,
          x<0.01 ? eps *std::fabs(d[i]) : eps,
          x<0.01 ? eps2*std::fabs(d[i]) : eps2,
          vflag>1
       );
       if (!e) {
          if (nmax<m) nmax=m;
          x=std::sqrt(x); if (d[i]<0) x=-x;
          r=std::fabs(x-d[i]);

          if (r>epsi) {
             wblog(F_L,"WRN fixing %.16g = %ssqrt(%ld/%ld) (%.2g)",
               d[i], d[i]<0 ? "-":"", labs(p),labs(q),r);
             wblog(FL,"TST %s() pqmax=%ld, eps=%g, eps2=%g,niter=%d",
               FCT,pqmax,eps,eps2,niter); 
             double d2=d[i]*d[i];
             e=Rational(d2,p,q,&r,&m,NULL,niter,pqmax,eps,eps2,'v');
          }
          x2=r*r; r2+=x2; if (ra<x2) ra=x2;
          d[i]=x; continue;
       }

       nmissed++;
   }

   if (ra_) { (*ra_)=ra; }

   if (nmissed && vflag) wblog(F_L,
      "WRN %s() %d/%d values unaltered (%d)\r\\",FCT,nmissed,n,nmax);

   return r2;
};


template<>
double SkipZeros(double *d, unsigned n, double eps){
  
   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (std::fabs(d[i])<eps) { s2+=(d[i]*d[i]); d[i]=0; }
   }

   return s2;
};

template<>
double SkipZeros(long double *d, unsigned n, double eps){
  
   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (std::fabs(d[i])<eps) { s2+=(d[i]*d[i]); d[i]=0; }
   }

   return s2;
};

template<>
double SkipZeros(wbcomplex *z, unsigned n, double eps){
  
   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (std::fabs(z[i].r)<eps) { s2+=(z[i].r*z[i].r); z[i].r=0; }
       if (std::fabs(z[i].i)<eps) { s2+=(z[i].i*z[i].i); z[i].i=0; }
   }

   return s2;
};


inline wbstring size2Str(double n) {
   unsigned long x; unsigned l=0; const size_t m=16; char s[m];
   x=1<<30; if (n>x) { l=snprintf(s,m,"%.4g GB",n/x); } else {
   x=1<<20; if (n>x) { l=snprintf(s,m,"%.4g MB",n/x); } else {
   x=1<<10; if (n>x) { l=snprintf(s,m,"%.4g kB",n/x); } else {
                       l=snprintf(s,m,"%.4g bytes",n);
   }}}
   if (l>=m) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,m);
   return s;
};

inline const char* filename(const char *s) {

   if (s==NULL) { wblog(FL,"WRN %s() got NULL string",FCT); }
   else {
      for (int i=strlen(s)-1; i>=0; --i) {
      if (s[i]=='/' || s[i]=='\\') return (s+i+1); }
   }
   return s;
};


size_t getFilesize(const char* f) {

   struct stat S;
   if (stat(f, &S)!=0) { return 0; }
   return S.st_size;   
};


bool fexist(const char *f, char type) {

   if (!f || !f[0]) wblog(FL,"ERR %s() invalid (empty) file name",FCT);

   struct stat S;
   if (stat(f, &S)!=0) { return 0; }

   if (type=='f') return (S.st_mode & S_IFREG); else
   if (type=='d') return (S.st_mode & S_IFDIR); else
   return 1;
};


unsigned long getSize(const char *f, const char *tag) {

   char *s, *line=NULL, u[8];
   unsigned long val; unsigned i=0; size_t len=0; ssize_t n;

   FILE *fid=fopen(f,"r");
   if (fid==NULL) {
      wblog(FL,"WRN %s() failed to open process status file\n%s",FCT,f);
      return 0;
   }

   while ((n=getline(&line,&len,fid))>=0) { ++i;
      if (n) line[n-1]=0;
      if ((s=strstr(line,tag))) { u[0]=0;
         sscanf(s,"%*s %lu %s",&val,u);
         if (strlen(u)==2 && u[1]=='B') {
             if (u[0]=='k' || u[0]=='K') val*=(1<<10); else
             if (u[0]=='M') val*=(1<<20); else
             if (u[0]=='G') val*=(1<<30);
             else wblog(FL,"WRN %s() invalid unit >%s<",FCT,u);
         }
         else wblog(FL,
            "WRN %s() invalid unit >%s< (%d)",FCT,u,strlen(u));
         break;
      }
   }

   if (line) { free(line); line=NULL; }
   fclose(fid); fid=NULL;

   if (n<0) wblog(FL,"WRN %s() '%s' not found in %s",FCT,tag,f);

   return val;
};


template <class TI> inline 
TI sub2ind(const TI *s, const TI *I, unsigned n){
   if (n) {
      unsigned j=n-1; TI k=I[j];
      for (--j; j<n; --j) k=k*s[j]+I[j];
      return k;
   }
   else { wblog(FL,"ERR %s() got null set !??", FCT); return -1; }
};

template <class T1, class T2> inline 
void ind2sub(T1 k, const T1 *s, T2 *I, unsigned n) {

   unsigned x,i=0;
   for (; i<n; ++i) {
      if (s[i]) { x=k/s[i]; I[i]=k-x*s[i]; k=x; }
      else {
         wblog(FL,"ERR %s() got null-size %s",FCT,
         (wbvector<T1>(n,s,'r')).toStrD().data);
      }
   }
   if (k) wblog(FL,
      "ERR %s() index out of bounds\n(%d => [%s; %s])",FCT,
      (wbvector<T2>(n,I,'r')+1).toStr().data,
      (wbvector<T1>(n,s,'r')).toStrD().data
   );
};


int atoi(const char *s, int &iout,
   unsigned n,
   char white
){
   if (!s || !s[0]) { return -1; }
   unsigned i=0; if (int(n)<=0) n=9;

   if (white) {
      for (; s[i]; ++i) { if (s[i]!=' ') break; }
      if (i) { 
         if (!s[i]) { return -2; }
         s+=i; i=0;
      }
   }

   if (white<0) {
      if (s[0]=='+' || s[0]=='-') { if (!s[++i]) return -3; }
   }

   for (; s[i] && i<n; ++i) {
      if (!isdigit(s[i])) return -4;
   }; if (i==n) return -5;

   iout=std::atoi(s);
   return i;
};

};

wbstring cpu_time::toStr(char flag) {

   size_t l=0, n=32; char tstr[n];

   cpu_time dt=since(); 
   double t=(flag=='c' ? dt.tcpu : dt.tsys);

   if (t<0) {
      wblog(FL,"WRN %s() got negative time (%ld)",FCT,long(t)); t=-t; }

   if (t<60) {
      l=snprintf(tstr,n,"%.4g sec",t); }
   else if (t<1E14) {
      int s,m,h,d; long x=t,
      fac=24*3600; d=t/fac; x-=d*fac;
      fac=3600;    h=x/fac; x-=h*fac;
      fac=60;      m=x/fac; x-=m*fac; s=x; x-=s;

      if (d) 
           l=snprintf(tstr,n,"%d-%02d:%02d:%02d",d,h,m,s);
      else if (h || m>9)
           l=snprintf(tstr,n,   "%02d:%02d:%02d",  h,m,s);
      else l=snprintf(tstr,n,        "%02d:%02d",    m,s);

      if (t>2E9) l=-l;
   }
   else {
      l=snprintf(tstr,n,"(WRN %.3g yrs (%ld sec))",
      t/(3600*24*365), long(t));
      l=-l;
   }

   if (l>=n) {
     long tl=(flag=='c' ? dt.tcpu : dt.tsys);
     double yrs=double(tl)/(3600*24*365);
     wblog(FL,
        "WRN cpu_time() %s `%s'\n'%c' dt=%ld => %.12g yrs (%d) !??\n"
        "%15.12g %15.12g %% now [CPU/SYS]\n"
        "%15.12g %15.12g %% then",
        int(l)<0 ? "possibly invalid CPU time":"string out of bounds\n",
        tstr, flag,tl,yrs,l, dt.time('c'), dt.time('s'), tcpu, tsys
     ); 
   }

   return tstr;
};


void Wb::CleanUp(const char *F, int L, const char *istr) {

   static int do_cleanup=1;
   if (F==NULL && L) { do_cleanup=L; return; }

   if (do_cleanup) { do_cleanup=0; wbtop PS;
       wbstring s1(PS.VmSize2Str()), s2(PS.VmPeak2Str());

       wblog(F_L," *  %s() MEM size: %s (%s)", istr ? istr : FCT,
       s1.data, s2.data);

       wblog(F_L," *  %s() time: %s (%s)", istr ? istr : FCT,
       gCPUTime.toStr('s').data, gCPUTime.toStr('c').data);
    }
};


inline wbstring wbsys::MemTot2Str() const {
   return Wb::size2Str(getMemTot()); };

inline wbstring wbsys::MemFree2Str() const {
   return Wb::size2Str(getMemFree()); };

inline wbstring wbsys::SwapTot2Str() const {
   return Wb::size2Str(getSwapTot()); };

inline wbstring wbsys::SwapFree2Str() const {
   return Wb::size2Str(getSwapFree()); };

char wbsys::checkSwapSpace(const char *F, int L) const {
   static double xref=0.25;
   double x=getSwapTot(); if (x<=0) return 0;
   x=getSwapFree()/x; if (x>xref) { return 0; }

   xref-=0.05; wblog(F_L,"WRN free swap space: %s / %s",
   SwapFree2Str().data, SwapTot2Str().data);
   return 1;
};

inline wbstring wbtop::VmSize2Str() const {
   return Wb::size2Str(getVmSize()); };
inline wbstring wbtop::VmPeak2Str() const {
   return Wb::size2Str(getVmPeak()); };

char wbtop::runningLarge(const char *F, int L,
   double th0,
   double fac
){
   static double thr=th0; 
   wbsys S; unsigned long s0=S.getMemTot(), M=getVmPeak();
   if (M<thr*s0) { S.checkSwapSpace(F_L); return 0; }
  
   wblog(F_L,"%s() process is using %s",FCT,VmPeak2Str().data);
   thr*=fac; return 1;
};


template <class T>
size_t Wb::findfirst_sorted(const char *F, int L,
   const T* d0, const T* dd, size_t m, size_t N, size_t M, char lex
){
   if (!N || !M || m>=M) { if (F) wblog(F,L,
      "ERR %s() got empty/invalid input (%dx%d; %d",FCT,N,M,m); 
      return -1;
   }

   int q=Wb::cmpRange(dd,d0,m,lex);
   if (q>=0) {
      if (q>0) { if (F) wblog(F,L,
         "ERR %s() empty match (<first)",FCT);
         return -1;
      }
      return 0;
   }

   size_t i1=0, i2=N-1, i=(i2+i1)/2;
   int q2=Wb::cmpRange(dd+i2*M,d0,m,lex);

   if (q2<0) { if (F) wblog(FL,
      "ERR %s() empty match (>last)",FCT);
      return N;
   }

   while (i>i1) {
      q=Wb::cmpRange(dd+i*M,d0,m,lex);
      if (q<0) { i1=i; }
      else { 
         if (q>0 && q!=q2) wblog(F_L,"ERR %s() "
            "got unsorted data (upper i=%d: %d/%d) !?",FCT,i,q,q2);
         q2=q; i2=i;
      }
      i=(i1+i2)/2;
   }

   if (q2) { if (F) wblog(FL,
      "ERR %s() empty match (i=%d/%d @ %d,%d)",FCT,i,N,q,q2);
      return -1;
   }
   return i2;
};


template <class T>
size_t Wb::findlast_sorted(const char *F, int L,
   const T* d0, const T* dd, size_t m, size_t N, size_t M, char lex
){
   if (!N || !M || m>=M) { if (F) wblog(F,L,
      "ERR %s() got empty/invalid input (%dx%d; %d",FCT,N,M,m); 
      return -1;
   }

   int q1=Wb::cmpRange(dd,d0,m,lex);

   if (q1>0) { if (F) wblog(F,L,
      "ERR %s() empty match (<first)",FCT);
      return -1;
   }

   size_t i1=0, i2=N-1, i=(i2+i1)/2;
   int q=Wb::cmpRange(dd+i2*M,d0,m,lex);

   if (q<=0) {
      if (q<0) { if (F) wblog(F,L,
         "ERR %s() empty match (>last)",FCT);
         return N;
      }
      return i2;
   }

   while (i>i1) {
      q=Wb::cmpRange(dd+i*M,d0,m,lex);
      if (q>0) { i2=i; }
      else { 
         if (q<0 && q!=q1) wblog(F_L,"ERR %s() "
            "got unsorted data (lower i=%d: %d/%d) !?",FCT,i,q,q1);
         q1=q; i1=i;
      }
      i=(i1+i2)/2;
   }

   if (q1) { if (F) wblog(FL,
      "ERR %s() empty match (i=%d/%d @ %d,%d)",FCT,i,N,q1,q);
      return -1;
   }
   return i1;
};


#endif

