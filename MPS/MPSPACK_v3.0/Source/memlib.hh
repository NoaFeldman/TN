#ifndef __WB_MEM_LIB_HH__
#define __WB_MEM_LIB_HH__

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //

template <class T>
bool gotMemOverlap(
   const T* src, size_t len0,
   const T* dest, size_t len2=-1 // default -1 => len0
);




template <class T>
inline void MEM_CPY_BASE(T* data, size_t n, T const* const d=NULL) {
#ifdef MLB_VFLAG
    wblog(FL,"TST base init"); 
#endif
    if (d==NULL)
         memset(data,0,n*sizeof(T));
    else memcpy(data,d,n*sizeof(T));
};

template <class T>
inline void MEM_CPY_BASE(
    T* data, size_t n, size_t n1, T const* const d1, T const* const d2
){
#ifdef MLB_VFLAG
    wblog(FL,"TST base init %d->%d",n1,n); 
#endif
    if (n1) {
       if (d1==NULL) wblog(FL,"ERR got (null) pointer"); else
       if (n1>n) { wblog(FL,"WRN length out of bounds (%d/%d)",n1,n); n1=n; }
       memcpy(data,d1,n1*sizeof(T));
    }
    if (n1<n) {
       if (d2) memcpy(data+n1,d2, (n-n1)*sizeof(T));
       else    memset(data+n1, 0, (n-n1)*sizeof(T));
    }
};



struct sptr_flags__ {

   sptr_flags__() : len(0), nrefs(0), isref(0) {};

  ~sptr_flags__() {
if (isref && nrefs) wblog(FL,"ERR %s() got (%d,%d) !??",FCT,isref,nrefs);
   };

   bool any() const { return (len || nrefs || isref); };
   bool isRef() const { 
      if (isref && nrefs) wblog(FL,
         "ERR %s() got (%d,%d) !??",FCT,isref,nrefs);
      return (nrefs || isref);
   };

   char check() const {
      if (len)
           { if (isref && nrefs) return 1; }
      else { if (isref || nrefs) return 2; }
      return 0;
   };

   wbstring toStr() const;

   size_t len;
   unsigned nrefs;
   char isref;
};

template <class T>
class wb_sptr {

  public:

   wb_sptr() : sptr_flags(NULL), data(NULL) {};

  ~wb_sptr() { 
if (check()) wblog(FL,"ERR %s() %d",FCT,check()); 
   delete_dref(); };

   void new_dref_bare(size_t n, char exflag=1);
   void new_dref(size_t n, const T *d0=NULL);

   void init2ref(const T* d0);

   unsigned add_dref(wb_sptr &b) const;
   void delete_dref();

   char check() const {
      if (data) {
         if (!sptr_flags) return 1;
         if (char e=sptr_flags->check()) return (10+e);
      }
      else if (sptr_flags) return 2;
      return 0;
   };

   unsigned nrefs() const {
      return (sptr_flags ? sptr_flags->nrefs : 0);
   };

   bool isRef() const { 
      return (sptr_flags && sptr_flags->isRef());
   };

   void Instantiate(const char *F=NULL, int L=0, size_t n=-1);

   wbstring toStr() const;


   sptr_flags__ *sptr_flags;
   T* data;

  protected:
  private:

};




template <class T>
class MEM_CPY {

  public:

     MEM_CPY(T* data, size_t n, const T* d=NULL) {

        if (WbUtil<T>::isPOD()) { MEM_CPY_BASE(data,n,d); }
        else { size_t i=0;
         #ifdef MLB_VFLAG
            wblog(FL,"TST class init"); 
         #endif
            if (d) { for (; i<n; ++i) data[i]=d[i]; }
            else {
               const T z=T();
               for (; i<n; ++i) data[i]=z;
            }
        }
     }

     MEM_CPY(T* data, size_t n, size_t n1, const T* d1, const T* d2) {

        if (WbUtil<T>::isPOD()) { MEM_CPY_BASE(data,n,n1,d1,d2); }
        else { size_t i=0;
         #ifdef MLB_VFLAG
            wblog(FL,"TST class init %d->%d",n1,n); 
         #endif
            if (n1) {
               if (d1==NULL) wblog(FL,"ERR got (null) pointer"); else
               if (n1>n) { wblog(FL,
                  "WRN length out of bounds (%d/%d)",n1,n); n1=n; }
               for (; i<n1; ++i) data[i]=d1[i];
            }
            if (d2) { for (; i<n; ++i) data[i]=d2[i-n1]; }
            else { 
               const T z=T();
               for (; i<n; ++i) data[i]=z;
            }
        }
     };

  private:

};


template <class T>
class MEM_CPY<T*> {

  public:

     MEM_CPY(T** data, size_t n, T* const* const d=NULL) {
#ifdef MLB_VFLAG
        wblog(FL,"TST got pointer space"); 
#endif
        MEM_CPY_BASE(data,n,d);
     }

     MEM_CPY(T** data, size_t n, size_t n1,
        T* const* const d1, T* const* const d2
     ){
#ifdef MLB_VFLAG
        wblog(FL,"TST got pointer space"); 
#endif
        MEM_CPY_BASE(data,n,n1,d1,d2);
     };

  private:

};


template <class T> inline 
void MEM_SET_BASE(T* data, size_t len, const T *a=NULL) {
#ifdef MLB_VFLAG
   wblog(FL,"TST base memset"); 
#endif
   if (a==NULL) memset(data,0,len*sizeof(T));
   else for (size_t i=0; i<len; i++) data[i]=(*a);
};

template <class T> inline
void MEM_SET_BASE(T** data, size_t len, T* const *a=NULL) {
#ifdef MLB_VFLAG
   wblog(FL,"TST base memset (pointer space)"); 
#endif
   if (a==NULL) memset(data,0,len*sizeof(T*));
   else for (size_t i=0; i<len; i++) data[i]=(*a);
};


template <class T>
class MEM_SET {

  public:

     MEM_SET(
        T* data      __attribute__ ((unused)),
        size_t len __attribute__ ((unused)),
        const T *a   __attribute__ ((unused)) = NULL
     ){
     };

  private:

};

template <class T>
class MEM_SET<T*> {

  public:

     MEM_SET(T** data, size_t len, T* const *a=NULL) {
#ifdef MLB_VFLAG
        wblog(FL,"TST pointer memset"); 
#endif
        MEM_SET_BASE(data,len,a);
     }

  private:

};


template <> inline MEM_SET<int>::MEM_SET(
   int* data, size_t len, const int *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<unsigned>::MEM_SET(
   unsigned* data, size_t len, const unsigned *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<unsigned long>::MEM_SET(
   unsigned long* data, size_t len, const unsigned long *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<long>::MEM_SET(
   long* data, size_t len, const long *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<char>::MEM_SET(
   char* data, size_t len, const char *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<double>::MEM_SET(
   double* data, size_t len, const double *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<float>::MEM_SET(
   float* data, size_t len, const float *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<wbcomplex>::MEM_SET(
   wbcomplex* data, size_t len, const wbcomplex *a
){ MEM_SET_BASE(data,len,a); }


#endif

