#ifndef __WB_COMPLEX_CC__
#define __WB_COMPLEX_CC__

// tags: SET, set_imag, set_val, set_num
template <class T> inline T DSET_IMAG(T &d, double i) {
   if (i!=0) wblog(FL,
      "ERR %s() got imaginary data %.3g for %s",
       FCT,i,getName(typeid(T)).data);
   return d;
};

template <> inline wbcomplex DSET_IMAG(wbcomplex &z, double i) {
   z.i=i; return z;
};


wbstring wbcomplex::toStr(const char *fmt) const {
    wbstring zstr(32), zfmt(16);
    unsigned k,n;

    if (r==0) {
       if (i==0) {
           zstr.printf(FL,fmt,r);
           return zstr;
       }

       zfmt.printf(FL,"%si",fmt);
       zstr.printf(FL,zfmt.data,i);
       return zstr;
    }
 
    if (i==0) {
       zstr.printf(FL,fmt,r);
       return zstr;
    }

   
    for (n=strlen(fmt), k=0; k<n; ++k) { if (fmt[k]=='+') break; }
    if (k<n) {
       zfmt.printf(FL,"%s%si",fmt,fmt);
    }
    else {
       zfmt.printf(FL,"%s+%si",fmt,fmt);
       zfmt[n]='%'; zfmt[n+1]='+';
    }

    zstr.printf(FL,zfmt.data,r,i);
    return zstr;
};


inline wbstring toStr(const wbcomplex &z, const char* fmt) {
   return z.toStr(fmt);
}

inline wbstring toStr(const double &d, const char* fmt) {

   size_t l=0, n=32; char s[n];
   l=snprintf(s,n,fmt,d);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};


#endif

