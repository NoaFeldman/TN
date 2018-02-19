#ifndef __WB_CLOCK_HH__
#define __WB_CLOCK_HH__

#include <set>

   class WbClock;
   set<WbClock*> gWbClocks;


class WbClock {

  public:

    WbClock(const char *s=NULL) {
       init(s); gWbClocks.insert(this); 
       TSC2sec(0);
    };

   ~WbClock() { info(); name.init(); gWbClocks.erase(this); };

    void Init() { info(); name.init(); }

    void init(const char *s=NULL) {
       tsec=tcpu=lcpu=lref=tref=nz=flag=used=0;
#ifdef PROG_TAG
       if (s) { sprintf(str,"%s:%s",PROG_TAG,s); name=str; }
       else { name = PROG_TAG; }
#elif defined(PROG)
       if (s) { sprintf(str,"%s:%s",PROG,s); name=str; }
       else { name = PROG; }
#else
       name = s ? s : "";
#endif
    };

    void reset() {
       tsec=tcpu=lcpu=lref=tref=nz=flag=0;
    };

    void start() {
       tsec=time((time_t *)NULL);
       tcpu=lcpu=lref=nz=0; flag=1;
       tref=clock();
       lref=rdtsc();
    };

    void resume() {
       if (!flag) { start(); return; }
       if (flag%2) return;

       tsec = time((time_t *)NULL)-tsec;
       tref = clock();
       lref = rdtsc();

       flag++;
    };

    int stop(const char *F=NULL, int L=0) {
       if (!flag) wblog(F_L,
          "ERR clock `%s' not yet started", name.data ? name.data:"");

       if (flag%2) { used++;
          clock_t t=clock();
          unsigned long long l2=rdtsc();
          flag++;

          lcpu += (l2-lref);
          tsec = time((time_t *)NULL)-tsec;
          if (t>tref) tcpu += t-tref;
          else if (t==tref) {
             if (l2<=lref) nz++; else
             tcpu += unsigned(TSC2sec(l2-lref)*CLOCKS_PER_SEC);
             return 1;
          }

       }
       return 0;
    };

    double getASM() const { return TSC2sec(lcpu); }
    wbstring gettime(const char *tag) const;

    mxArray* toMx(const char *F=NULL, int L=0);

    void Stop(const char *F, int L, const char *istr=NULL) {
       stop(F_L);
       Info(F_L, istr ? istr : name.data);
       reset();
    };

    void info(const char *istr=NULL) const;
    void Info(const char *F, int L, const char *istr) const;

    double getTickFreq() const;

    static unsigned long long rdtsc_2 () {
       unsigned long long l;
       __asm__ __volatile__("rdtsc" : "=A" (l) : );
       return l;
    };

    static unsigned long long rdtsc () {
       unsigned long low, high;
       unsigned long long l;
       __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high) : );
       l=high; l=l<<32 | low;
       return l;
    };

    void set_TSC_TO_SEC() { TSC_TO_SEC=1.0/getTickFreq(); };

    clock_t tsec;
    time_t  tref;
    unsigned nz, flag;
    unsigned long tcpu;
    unsigned long long lcpu, lref, used;
    wbstring name;

    static double TREF;
    
  protected:
  private:

    static double TSC_TO_SEC;

    template<class T>
    void disptime(T dt, const char *istr="  ") const;

    double TSC2sec(const unsigned long long &tsc) const {
       if (TSC_TO_SEC==0) {
          double tf=getTickFreq();
          wbstring h(32); gethostname(h.data,h.len);
          wblog(FL,"<i> CPU frequency rdtsc: %.3g GHz (%s)",tf*1E-9,h.data);
          TSC_TO_SEC=1.0/tf;
       }
       return (tsc*TSC_TO_SEC);
    };
};

   double WbClock::TREF=0;
   double WbClock::TSC_TO_SEC=0;


unsigned showAllClocks() {
   set<WbClock*>::iterator it=gWbClocks.begin();
   for (; it!=gWbClocks.end(); ++it) { (*it)->info(); }
   return gWbClocks.size();
};

unsigned initAllClocks() {
   set<WbClock*>::iterator it=gWbClocks.begin();
   for (; it!=gWbClocks.end(); ++it) { (*it)->info(); (*it)->reset(); }
   return gWbClocks.size();
};


namespace Wb {

class UseClock {
  public:

    UseClock(WbClock *c) : C(c) { if (C) C->resume(); };
   ~UseClock() { if (C) C->stop(); };

    void done() { if (C) { C->stop(); C=0; }};
    void stop() { if (C) C->stop(); };
    void resume() { if (C) C->resume(); };

  protected:
  private:

    WbClock *C;
};

};


inline wbstring WbClock::gettime(const char *tag) const {

   double dt=0; int it;

   if (!strcmp(tag,"int")) {
      if (flag%2)
           dt=time((time_t *)NULL)-tsec;
      else dt=tsec;
   }
   else if (!strcmp(tag,"cpu")) {
      if (flag%2)
           dt=double(clock()-tref)/CLOCKS_PER_SEC;
      else dt=double(tcpu)/CLOCKS_PER_SEC;
   }
   else if (!strcmp(tag,"asm")) {
      dt=TSC2sec(lcpu);
   }
   else wblog(FL,"ERR invalid time tag `%s'",tag);

   it=(int)dt;

   if (dt>10) {
      sprintf(str,"%02d:%02d:%02d",
      it/3600, (it/60)%60, it%60);
   }
   else {
      sprintf(str,"%02d:%02d:%02d (%.3g sec)", 
      it/3600, (it/60)%60, it%60, dt);
   }

   return str;
};


template<class T>
inline void WbClock::disptime(T it, const char *istr) const {
   printf("%s%02d:%02d:%02d", istr, it/3600, (it/60)%60, it%60);
};

template<>
inline void WbClock::disptime(double dt, const char *istr) const {
   unsigned it=(int)dt;
   if (it<10)
        printf("%s%8.3g",istr,dt);
   else disptime(it,istr);
};


void WbClock::info(const char *istr) const {
   int it; double dc,da;
   static int first_call=1;

   if (!flag || !used) return;

   if (flag%2) {
      it=time((time_t *)NULL)-tsec;
      dc=double(clock()-tref)/CLOCKS_PER_SEC;
   }
   else {
      it=tsec;
      dc=double(tcpu)/CLOCKS_PER_SEC;
   }

   if (nz) {
      double p0=1-nz/(double)(flag/2);
      sprintf(str,"/%5.1f%%", 100*p0);
   }  else str[0]=0;

   da=TSC2sec(lcpu);

   if (first_call) { wblog(FL,"Clock statistics%N"); printf(
   "  CLOCK NAME                     COUNT  "
   "  i_TIME()  d_CLOCK()    ASM()     /CALL\n\n"); first_call=0; }

   printf("%c %-24s%12d %s ",
   flag%2 ? '*':' ', istr ? istr : name.data, flag/2, str);

   disptime(it);
   disptime(dc);
   disptime(da);

   printf(" %9.3g ", da/(flag/2));

   if (TREF)
        printf("%6.1f%%\n", 100*dc/TREF);
   else printf("\n");
};

void WbClock::Info(const char *F, int L, const char *istr) const {

   size_t l=0, n=128, n0; char s[n];
   int it,ct,at; double dt;

   strcpy(s,"--------------------------------------------");
   n0=strlen(s);

   l=snprintf(s,n,"%s - time statistics (%d)",istr,flag/2);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   if (l<n0) s[l]=' ';

   it=tsec;
   dt=double(tcpu)/CLOCKS_PER_SEC; ct=(int)dt;
   at=(int)TSC2sec(lcpu);

   wblog(F_L,
       "%s\n"
       "    CPU time        : %02d:%02d:%02d (%g)\n"
       "    elapsed time    : %02d:%02d:%02d (%d)\n"
       "    ASM_TSC time    : %02d:%02d:%02d (%g)\n"
       "--------------------------------------------", s,
       ct/3600, (ct/60)%60, ct%60, dt,
       it/3600, (it/60)%60, it%60, it,
       at/3600, (at/60)%60, at%60, TSC2sec(lcpu)
   );

   if (nz) wblog(F_L,
   "NB! %d/%d time increments too short (0.)", nz, flag/2);
};


double WbClock::getTickFreq() const {

   unsigned i,n=32;
   clock_t t1,t2,nc=0;

   unsigned long long l1,l2,dl=0;

   for (i=0; i<n; i++) {
      l1=rdtsc(); t1=clock(); while (t1==clock()) {};
      l2=rdtsc(); t2=clock();

      if (i) {
         nc+=(t2-t1);
         dl+=(l2-l1);
      }
   }


   return (dl*(CLOCKS_PER_SEC/(double)nc));
}


#endif

