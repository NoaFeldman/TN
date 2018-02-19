#ifndef __WB_MEM_TRACK_CC__
#define __WB_MEM_TRACK_CC__

// NB! included in wblib.h only if __WBDEBUG__ is set
// see also memlib.hh (always included)


namespace Wb {

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

class MREC {
   public:

#ifdef EXTENDED_MEM_CHECK
      MREC() : id(0), len(0), unit(0), fline(0) {};
     ~MREC() { if (fline) { delete [] fline; fline=0; }};
#else
      MREC() : id(0), len(0), unit(0) {};
     ~MREC() {};
#endif

#ifdef EXTENDED_MEM_CHECK
      template <class T>
      const char* init(
         const char *F, int L, const T* p __attribute__ ((unused))
      ){
         unit=sizeof(T);
         if (F) {
            const std::type_info &tid = typeid(T);
            const char *s=shortFL(F,L), *t=tid.name();
            unsigned n=strlen(s) + strlen(t) + 5;
            if (fline) { delete [] fline; }
            fline = new char[n]; 
            snprintf(fline,n,"%s (%s)",s,t);
         }
         return fline;
      };
#else
      template <class T>
      const char* init(
         const char *F __attribute__ ((unused)),
         int L         __attribute__ ((unused)),
         const T* p    __attribute__ ((unused))
      ){
         unit=sizeof(T); return 0;
      };
#endif

      unsigned long long id;
      size_t len;
      unsigned unit;

#ifdef EXTENDED_MEM_CHECK
      char *fline;
#endif

   protected:
   private:
};

class WbListMTRACK {
   public:
     WbListMTRACK() : idx(0), totsize(0), maxsize(0) {
#ifndef __WB_MEM_QUIET__
        wblog(FL, "MTR ************ MEM TRACK : ON ********************");
#endif
     };

     map<void*, MREC> M;

     unsigned long long idx;
     unsigned long long totsize;
     unsigned long long maxsize;

     template<class T>
     void rm_ptr(const char* F, const int &L, T* &p) { if (p) {
        #pragma omp critical
        {  map<void*,MREC>::iterator im = M.find((void*)p);
           if (im == M.end()) wblog(F_L,"WRN memory to be freed "
              "not in list (0x%lx; #%d) !!?",p,omp_get_thread_num());
           else {
              totsize -= (im->second.len) * (im->second.unit);
              M.erase(im);
           }
        }
     }};

     template<class T>
     void add_ptr(const char* F, int L, T* &p) {

        if (p==NULL) wblog(F_L,"ERR got NULL pointer (out of memory !!?)");

        #pragma omp critical
        {  MREC &r=M[(void*)p];

           if (r.id!=0) wblog(F_L,"ERR entry with pointer "
              "%lx (%ld,%d,%d) already exists!!",p,r.id,r.len,r.unit);

           r.id=(++idx);
           r.len=1;
           r.init(F_L,p);

           totsize += sizeof(T);
           checkTotSize();
        }
     };

     template<class T>
     bool exists(T* &p) {
        map<void*,MREC>::iterator im = M.find(p);
        return (im!=M.end() ? 1 : 0);
     };

     template<class T>
     void tst_exists(const char* F, int L, T* &p) {
        wblog(FL,"TST 0x%lx %sin list", p, exists(p) ? "":"not ");
     };

     void checkTotSize() { if (totsize>maxsize) {
        size_t stot=(totsize>>30);
        if (stot) {
           if (!(stot>>2)) {
              size_t smax=(maxsize>>30);
              while (smax) { smax>>=1; stot>>=1; }
              if (stot) wblog(FL,
              " *  %s[] reached %s  ",FCT,totStr());
           }
           else {
              size_t smax=(maxsize>>32); stot>>=2;
              if (stot>smax) wblog(FL,
              " *  %s[] reached %s  ",FCT,totStr());
           }
        }
        maxsize = totsize;
     }};

     char* totStr(char mflag=0) { size_t n=16, l=0;

        if (mflag) { char s1[n], s2[n];
           memsize2Str(totsize,s1);
           memsize2Str(maxsize,s2); l=snprintf(sbuf,32,"%s (%s)",s1,s2);
        }
        else {
           memsize2Str(totsize,sbuf); l=strlen(sbuf);
        }

        if (l>=32) wblog(FL,
           "ERR %s() string out of bounds (%d/32)",FCT,l); 
        return sbuf;
     };

     void printSize(const char *F, int L, char mflag=0) {
        if (mflag) { char s1[16], s2[16]; wblog(F_L,
         " *  MEMTrack got %s (max. %s)  ",
              memsize2Str(totsize,s1), memsize2Str(maxsize,s2));
        }
        else { char s1[16]; wblog(F_L,
         " *  MEMTrack got %s   ", memsize2Str(totsize,s1));
        }
     };

   protected:
   private:

     char sbuf[32];
};

WbListMTRACK gML;


template<class T>
inline void NEW(const char* F, int L, T* &p, size_t n) {

   if (n) {
      p = new T[n];
      gML.add_ptr(F_L,p);
   }
   else {
      wblog(FL,"WRN allocate empty space (%s, %d)",shortFL(F,L),n);
      p=NULL;
   }
};



int MemCheck(const char *F="", int L=0,
    const char *task="",
    long dl=0, long db=0
){
    static unsigned long long nbytes=0, iref=0;
    static size_t len=0;

    if (task!=0 ? task[0]==0 : task==NULL) {
       wblog(F,L, "MTR info: %d entries (%s)",
       gML.M.size(), gML.totStr());
       return 0;
    }

    if (!strcmp(task,"start")) {
       len    = gML.M.size();
       nbytes = gML.totsize;
       iref   = gML.idx;
    }
    else if (!strcasecmp(task,"stop")) {
       size_t l=gML.M.size(); 
       unsigned long long b=gML.totsize;

       int m=0, e = long(l-len) != dl || long(b-nbytes) != db;

       if (e) {
          map<void*,MREC>::iterator im;
          for (im=gML.M.begin(); im!=gML.M.end(); ++im) {
              if (im->second.id>iref) { if (!(m++)) printf("\n"); if (m<20)
                 printf("%6d: %9ld %10lx %8d /%d\n",m,
                 (unsigned long)im->second.id,
                 (unsigned long)im->first, im->second.len, im->second.unit);
              }
          }
          if (m>20) wblog(FL,
             "%NTST %s() %d/%d entries skipped%N",FCT,m-20,l);
          else if (m) printf("\n");
       }

       wblog(F,L,
         "%s  %d entries (%+d) using %s (%+ld)",
          !strcmp(task,"STOP") ? (e ? "ERR" : "SUC") : "MTR",
          l, l-len, gML.totStr(), b-nbytes
       );
    }
    else if (!strcasecmp(task,"info")) {
       size_t l=gML.M.size(); long s=gML.totsize-nbytes;
       char ss[16];

       wblog(F,L,
         "MTR %d entries (%+d) using %s (%s%s; %lld)",
          l, l-len, gML.totStr(), s>=0 ? "+":"",
          memsize2Str(s,ss), gML.idx
       );
    }
    else if (!strcasecmp(task,"LIST")) {
       size_t i; map<void*, MREC>::iterator I; char s[32];
       MemCheck(F,L,"info",dl,db);
       for (i=0, I=gML.M.begin(); I!=gML.M.end() && i<501; ++I) {
           const MREC &r=I->second;
           snprintf(s,32,"%d * %d", r.len, r.unit);
           printf("%6d: 0x%7lx, id=%12ld, %12s   %s\n",++i,
           (unsigned long)I->first, (unsigned long)r.id, s,
#ifdef EXTENDED_MEM_CHECK
           r.fline ? r.fline : ""
#else
           ""
#endif
           );
       }
       if (I!=gML.M.end())
       printf("\n  ... and %ld others\n\n",long(gML.M.size()-i));
    }
    else wblog(F,L, "ERR MTRACK invalid flag >%s<", task);

    return 0;
};


}

#ifndef __WB_MEM_QUIET__

class wbdebug_dummy {
  public:

    wbdebug_dummy() {
       wblog(FL,"MTR ************ DEBUG MODE: ON ********************");
       Wb::MemCheck(FL,"info"); printf("\n");
    };

   ~wbdebug_dummy() {
       Wb::MemCheck(FL,"LIST");
    };

  protected:
  private: unsigned whatever;
};

   wbdebug_dummy wbdebug_dummy_var;

#endif

#endif

