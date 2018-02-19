#ifndef __WB_LOG_H__
#define __WB_LOG_H__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

// unsigned WBLOG_ERR_COUNT=0;
// unsigned WBLOG_ERR_CONT=0;

   void wbSetLogLevel(unsigned l);
   int wblog(const char* file, int line, const char *fmt, ...);

   unsigned wblog_checktag(const char *istr, const char *tag);
   char wblog_findtoken(const char *istr, const char *tok, int maxoffset=-1);

   void banner(
      unsigned n, const char *s,
      const char *istr=NULL, const char *fstr=NULL
   );


inline char* shortFL(const char *F, int L=-1, char *s=NULL) {
    static char fstr[512], ic=0;

    unsigned i=0, j=0;
    if (s==NULL) { s=fstr+128*ic; if ((++ic)>=4) ic=0; }

    if (F) {
       for (; F[i]; ++i) if (F[i]=='/') j=i+1;

       if (L>0)
             snprintf(s,128,"%s:%d",F+j,L);
       else  snprintf(s,128,"%s:",  F+j);
    }
    else s[0]=0;

    return s;
};


const char* shortPF(const char *P, const char *F, int L, char *s) {
    static char pstr[512], ic=0;

    unsigned i=0, j=0;
    if (s==NULL) { s=pstr+128*ic; if ((++ic)>=4) ic=0; }
    s[0]=0;

    if (F) {
       for (; F[i]; ++i) { if (F[i]=='/') j=i+1; }; F+=j; i=j=0;
       if (!P || !P[0]) return F;
       for (; F[i]; ++i) { if (F[i]=='.') j=i+1; }; if (!j) j=i;
    }
    else { return (P && P[0] ? P : s); }

    unsigned l1=strlen(P), l2=strlen(F);

    unsigned ltot=l1+l2+1, lmax=19-(L>0 ? log10(double(L)) : 0);

    if (ltot<lmax) {
       snprintf(s,128,"%s:%s",P,F);
    }
    else if (l1<9) { char fmt[16];
       unsigned n2x=l2-j, n2=(lmax-l1-n2x-3);
       snprintf(fmt,16,"%%.%ds'%%.%ds(%%s)", n2,n2x);
       snprintf(s,128,fmt, F, F+j, P);
    }
    else if (l2<9) { char fmt[16];
       unsigned n1=(lmax-l2-2);
       snprintf(fmt,16,"%%s(%%.%ds)", n1);
       snprintf(s,128,fmt, F, P);
    }
    else { char fmt[32];
       unsigned n2x=l2-j, n2=1+((l2-n2x-1)*lmax)/ltot, n1=lmax-n2-n2x-3;
       snprintf(fmt,32,"%%.%ds'%%.%ds(%%.%ds)",n2,n2x,n1);
       snprintf(s,128,fmt, F, F+j, P);
    }

    return s;
};


#endif

