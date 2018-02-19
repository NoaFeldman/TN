#ifndef __WB_STRING_CC__
#define __WB_STRING_CC__

// ----------------------------------------------------------------- //
// tags: filesize fsize kB MB GB

char* memsize2Str(double x, char *s, unsigned l) {
    if (int(l)<0) {
       if (x<5E3) sprintf(s,"%g bytes",x);        else
       if (x<1E6) sprintf(s,"%.3g kB",x/(1<<10)); else
       if (x<1E9) sprintf(s,"%.3gM",x/(1<<20));   else
                  sprintf(s,"%.3gG",x/(1<<30));
    }
    else {
       unsigned i=0; if (l) {
       if (x<5E3) i=snprintf(s,l,"%g bytes",x);        else
       if (x<1E6) i=snprintf(s,l,"%.3g kB",x/(1<<10)); else
       if (x<1E9) i=snprintf(s,l,"%.3gM",x/(1<<20));   else
                  i=snprintf(s,l,"%.3gG",x/(1<<30));
       }
       if (i>=l) wblog(FL,
         "WRN %s() insufficient length of string (%s; %d/%d)",
          FCT,s,strlen(s),l
       );
    }
    return s;
};


inline wbstring char2Str(char c) {
   size_t l=0, n=8; char s[n];
   if (isprint(c))
        l=snprintf(s,n,"'%c'",c);
   else l=snprintf(s,n, "%d", c);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   return s;
};

  wbstring repHome(const wbstring &file) {
     return Wb::repHome(file.data);
  };


template<class T>
wbstring wbstring::operator<< (const T& x __attribute__ ((unused)) ) {

   wblog(FL,"ERR %s() unsupported%s() type '%s'",
   FCT,getName(typeid(T)).data);

   return wbstring();
};

template<class T>
wbstring wbstring::operator<< (const T* x) {
   size_t n=16; char s[n];
   size_t l=snprintf(s,n,"%lX",x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   (*this)+=s; return *this;
};


template<>
wbstring wbstring::operator<< (const char* s) {
   (*this)+=s; return *this;
};


template<>
wbstring wbstring::operator<< (const double& x) {
   size_t n=16; char s[n];
   size_t l=snprintf(s,n,"%.4g",x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   (*this)+=s; return *this;
};

template<>
wbstring wbstring::operator<< (const float& x) {
   size_t n=16; char s[n];
   size_t l=snprintf(s,n,"%.4g",x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   
   (*this)+=s; return *this;
};

template<>
wbstring wbstring::operator<< (const wbcomplex& x) {
   return x.toStr();
};


template<>
wbstring wbstring::operator<< (const unsigned& x) {
   size_t n=16; char s[n];
   size_t l=snprintf(s,n,"%d",x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   (*this)+=s; return *this;
};

template<>
wbstring wbstring::operator<< (const int& x) {
   size_t n=16; char s[n];
   size_t l=snprintf(s,n,"%d",x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   (*this)+=s; return *this;
};

template<>
wbstring wbstring::operator<< (const char& x) {
   size_t n=8; char s[n];
   size_t l=snprintf(s,n,"%c<%d>",x,x);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   (*this)+=s; return *this;
};


inline wbstring& wbstring::push(const char *F, int L,
    const char* s,
    char extend  __attribute__ ((unused))
 ){
    if (s && s[0]) {
       unsigned n1=strlen(data), n2=strlen(s);
       if (n1+n2>=len) wblog(F_L,
          "ERR %s() string out of bounds (%d+%d/%d)",FCT,n1,n2,len
       );
       strcat(data+n1,s);
    }
    return *this;
};

inline wbstring& wbstring::pushf(const char *F, int L,
    const char *fmt, ...
 ){
    if (fmt && fmt[0]) {
       size_t l, n1=strlen(data), n2=len-n1;
       if (len<=1) wblog(F_L,"ERR %s() got empty string (%d)",FCT,len);
       if (n1>=len) wblog(F_L,
          "ERR %s() got string without null-terminator (%d)!",FCT,len);

       va_list args;
       va_start(args,fmt);
       l=vsnprintf(data+n1,n2,fmt,args);
       va_end(args);

       if (l>=n2) { data[len-1]=0; wblog(F_L,
          "ERR %s() string out of bounds (%d->%d/%d)%N%N> %s%N",
           FCT,n1,n1+l,len,data);
       }
    }
    return *this;
};


inline wbstring& wbstring::printf(const char *F, int L,
    const char *fmt, ...
 ){
    if (fmt && fmt[0]) {
       if (len<=1) wblog(F_L,"ERR %s() got empty string (%d)",FCT,len);

       va_list args;
       va_start(args,fmt);
       size_t l=vsnprintf(data,len,fmt,args);
       va_end(args);

       if (l>=len) { data[len-1]=0; wblog(F_L,
          "ERR %s() string out of bounds (%d/%d)%N%N> %s%N",
           FCT,l,len,data);
       }
    }
    return *this;
};

inline wbstring& wbstring::cpy(const char *F, int L,
    const char *s
 ){
    if (s) {
       size_t l=strlen(s);
       if (l+1>=len) wblog(F_L,
          "ERR %s() string out of bounds (%s; %d/%d)",FCT,s,l,len);
       if (s[0]) strcat(data,s);
       else data[0]=0;
    }
    else wblog(F_L,"ERR %s() got NULL string",FCT); 

    return *this;
};


#endif

