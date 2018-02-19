#ifndef __WB_STRING_HH__
#define __WB_STRING_HH__

// NB! convert (possible) long long on input to double
// NB! see header in wblibx.h 
   char* memsize2Str(double x, char *s, unsigned l);

   wbstring char2Str(char c);
   wbstring repHome(const wbstring &file);


class wbstring : public wbvector<char> {

  public:

    wbstring (unsigned n=0) : wbvector<char>(n) {};

    wbstring (unsigned n, const char c)
     : wbvector<char>(n+1) { set(c); data[len-1]=0; };

    wbstring (const char* i,unsigned n)
     : wbvector<char>(n+1,i) { data[len-1]=0; };

    wbstring (const wbstring &s)
     : wbvector<char>(strlen(s.data)+1, s.data) { data[len-1]=0; };

    wbstring (const char* s1)
     : wbvector<char>(strlen(s1)+1,s1) { data[len-1]=0; };

    wbstring (
       const char* s1, const char* s2,
       const char* s3=NULL, const char* s4=NULL
    ) : wbvector<char>() {
       unsigned 
          l1 = (s1 ? strlen(s1) : 0),
          l2 = (s2 ? strlen(s2) : 0),
          l3 = (s3 ? strlen(s3) : 0),
          l4 = (s4 ? strlen(s4) : 0); NEW(l1+l2+l3+l4+1);
       if (l1) strcpy(data,         s1);
       if (l2) strcpy(data+l1,      s2);
       if (l3) strcpy(data+l1+l2,   s3);
       if (l4) strcpy(data+l1+l2+l3,s4); data[len-1]=0;
    };


#ifndef NOMEX
    wbstring(const mxArray *a) : wbvector<char>() { init(0,0,a); };
    wbstring& init(const mxArray *a) { return init(0,0,a); }
#endif

    wbstring (const char *F, int L, const mxArray *a)
     : wbvector<char>() { init(F,L,a); };


    wbstring& init(unsigned  n=0) {  RENEW(n); return *this; };

    wbstring& init(const char *s) {
       RENEW(strlen(s)+1, s); data[len-1]=0;
       return *this;
    };

    wbstring& init(const char *F, int L, const mxArray *a) {
       int i,n;
       if (!a || !(n=mxGetNumberOfElements(a))) {
          RENEW(0); return *this;
       };

       if (!mxIsChar(a)) wblog(F_L,
       "ERR initializing string with type '%s'",mxGetClassName(a));

       RENEW(++n);
       i=mxGetString(a,data,n);

       if (i) wblog(F_L,"ERR failed to read string >%s<",data);
       return *this;
    };

    void set(char c)  {
       if (len) {
       memset(data, c, len-1); data[len-1]=0; }
    };

    wbstring& operator= (const char* i) {
        RENEW(strlen(i)+1, i); data[len-1]=0; return *this;
    };

    bool gotTerm0() {
        for (unsigned i=0; i<len; ++i) { if (data[i]==0) return 1; }
        return 0;
    };

    wbstring& push(const char *F, int L, const char* s, char extend=0);
    wbstring& push(const char* s, char extend=0) {
        return push(0,0,s,extend); };

    wbstring& pushf(const char *F, int L, const char *fmt, ...);
    wbstring& printf(const char *F, int L, const char *fmt, ...);
    wbstring& cpy(const char *F, int L, const char *s);

    wbstring& time() { return time(::time(NULL)); };
    wbstring& time(const time_t &t) { NEW(32);
       ctime_r(&t,data);
       return *this;
    };

    wbstring& time_sys() { return time_sys(::time(NULL)); };
    wbstring& time_sys(const time_t &t) { NEW(32);
       strftime(data,31,"%a %b %d %H:%M:%S %Z %Y",localtime(&t));
       return *this;
    };

    wbstring& Upper() {
        for (unsigned n=strlen(data), i=0; i<n; i++)
        data[i]=std::toupper(data[i]);
        return *this;
    };

    wbstring& Lower() {
        for (unsigned n=strlen(data), i=0; i<n; i++)
        data[i]=std::tolower(data[i]);
        return *this;
    };

    wbstring toupper() const {
    wbstring sout(*this); sout.Upper(); return sout; };

    wbstring tolower() const {
    wbstring sout(*this); sout.Lower(); return sout; };

    wbstring& Chomp() {
       if (data && data[0]) { int i=strlen(data);
          for (; i>=0; --i) {
             if (data[i]>20 && data[i]!=' ') { data[i+1]=0; break; }
          }
          if (i<0) data[0]=0; 
       }
       return *this;
    };

    char getopt(char c) {
        for (unsigned i=0; data[i]; i++) {
           if (data[i]==c) {
              data[i]=-1;
              return 1;
           }
        }
        return 0;
    }
    char getopt() {
        for (unsigned i=0; data[i]; i++) {
           if (data[i]>0) { return data[i]; }
        }
        return 0;
    }

    const char* basename(const char c='/') const {
        if (!data) return data;
        return Wb::basename(data,c);
    };

    mxArray* toMx() const { return mxCreateString(data); };

    void add2MxStruct(mxArray *S, unsigned i) const {
       mxSetCell(S,i,toMx());
    };

    template <class T>
    char* init2Fmt(const T& x __attribute__ ((unused)), int n=-1, int p=-1){
       char s[32]; num2Fmt<T>(FL,s,n,p);
       init(s); if (len>28) wblog(FL,
         "WRN %s() check format bounds (%d/32) !??",FCT,len);
       return data;
    };

    template <class T>
    wbstring operator<< (const T* x);
    template <class T>
    wbstring operator<< (const T& x);

    wbstring operator+ (const wbstring &s) const { return (*this)+s.data; };
    wbstring operator+ (const char *s) const {
        wbstring sout(strlen(data) + strlen(s));
        strcpy(sout.data, data); strcat(sout.data, s);
        return sout;
    };

    wbstring& operator+=(const char *s) {
       if (data && s) {
          unsigned n1=strlen(data), n2=strlen(s);
          if (n2) { Resize(n1+n2); strcat(data, s); }
       }
       else if (s) { (*this)=s; }
       return *this;
    };

    wbstring& operator+=(const wbstring &s) {
       (*this)+=s.data; return *this;
    };

    unsigned cat(unsigned &k, const char *s) {
       if (s) { unsigned i=0;
          for (; s[i]!=0 && k<len; ++i) { data[k++]=s[i]; }
          if (s[i]) wblog(FL,"ERR wbstring::cat out of bounds (%d/%d)",k,len);
       }
       data[k]=0; return k;
    };


    char& operator[] (unsigned i) const { return data[i]; };

    unsigned length()  const { return strlen(data); };
    char     isEmpty() const { return (data ? data[0]==0 : 1); };

    unsigned isName(int L=-1) const {
        unsigned i, n=strlen(data);
        if (L>=0 && L<(int)n) return 0;
        for (i=0; i<n; i++) if (!isalnum(data[i])) return 0;
        return n;
    };

    int operator== (const char *s) const {
        return (strcmp(data, s)==0);
    };
    int operator!= (const char *s) const {
        return (strcmp(data, s)!=0);
    };

    int operator== (const wbstring &s) const {
        return strcmp(data, s.data)==0;
    };
    int operator!= (const wbstring &s) const {
        return strcmp(data, s.data)!=0;
    };
    int operator<  (const wbstring &s) const {
        return strcmp(data, s.data)<0;
    };
    int operator<= (const wbstring &s) const {
        return strcmp(data, s.data)<=0;
    };
    int operator>  (const wbstring &s) const {
        return strcmp(data, s.data)>0;
    };
    int operator>= (const wbstring &s) const {
        return strcmp(data, s.data)>=0;
    };

    int firstLocation (const wbstring &, unsigned k=0);
    int firstLocation (const char*, unsigned offset=0);

    std::istream & getline (std::istream &istr, char c='\n');


  protected:
  private:

};


std::istream & wbstring::getline (std::istream &istr, char c) {
    const unsigned maxlen=256;
    char tmp[maxlen];
    istr.getline(tmp, maxlen-1, c);

    init(strlen(tmp)+1); strcpy (data, tmp);

    return istr;
}


int wbstring::firstLocation (const wbstring &s, unsigned k) {
    return firstLocation(s.data,k);
}

int wbstring::firstLocation (const char *s, unsigned offset) {
    int i, j, n=strlen(data), m=strlen(s);

    for (j=0, i=offset; i<n && j<m; i++, j++) {
        if (data[i] != s[j]) {
            i-=j;
            j=-1;
        }
    }

    if (j==m) return i-m;
    else return -1;
}


#include <queue>

class wblogbuf_struct {

 public:

   wblogbuf_struct() : file(NULL), line(0), text(NULL) {};

   wblogbuf_struct (const char *s) : file(NULL), line(0), text(NULL) {
      if (s) {
         unsigned n=strlen(s);
         WB_NEW(text,n+1); strcpy(text,s);
      }
   };

   wblogbuf_struct (const char *F, int L, const char *s)
    : file(NULL), line(L), text(NULL) { unsigned i=0,k=0,n;
      if (F) {
         for (; F[i]; i++) { if (F[i]=='/' || F[i]=='\\') k=i+1; }
         n=strlen(F+k); WB_NEW(file,n+1); strcpy(file,F+k);
      }
      if (s) { n=strlen(s); WB_NEW(text,n+1); strcpy(text,s); }
   };

   wblogbuf_struct(const wblogbuf_struct &S)
    : file(S.file), line(S.line), text(S.text
   ){
   };

  ~wblogbuf_struct (){
      if (file) { WB_DELETE(file); }
      if (text) { WB_DELETE(text); }
      line=0;
   };

   void unset() { file=NULL; line=0; text=NULL; };

   char *file; int line; char *text;

 protected:
 private:
};


class wblogbuffer {

  public:



   void push(const char *F, int L, const char *s) {
      wblogbuf_struct S(F,L,s);
      buf.push(S); S.unset();
   };

   void flush() {
      while (!buf.empty()) {
         wblogbuf_struct S=buf.front();
         wblog(
            S.file ? S.file : __FILE__,
            S.line ? S.line : __LINE__, S.text
         );
         S.unset();
         buf.pop();
      }
      doflush();
   };

 protected:
 private:

   queue<wblogbuf_struct> buf;
};

  wblogbuffer wblogBuf;


#endif

