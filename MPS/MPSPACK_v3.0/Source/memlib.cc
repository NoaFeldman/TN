#ifndef __WB_MEM_LIB_CC__
#define __WB_MEM_LIB_CC__

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //

wbstring sptr_flags__::toStr() const {

   unsigned l=0, n=16; char s[n]; s[0]=0;

   if (l<n && len  ) l =snprintf(s,  n,  "l=%ld",len);
   if (l<n && nrefs) l+=snprintf(s+l,n-l,"%snrefs=%d", l?", ":"", nrefs);
   if (l<n && isref) l+=snprintf(s+l,n-l,"%s*ref*", l?", ":"");

   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   return s;
};


template <class T>
wbstring wb_sptr<T>::toStr() const {

   if (!data) {
      if (sptr_flags) {
         unsigned l,n=32; char s[n];
         l=snprintf(s,n,"null data (%s)",sptr_flags->toStr().data);
         if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
         return s;
      }
      else return "(null)";
   }
   else {
      unsigned l,n=64; char s[n];
      l=snprintf(s,n,"data=0x%px (%s)", data, sptr_flags->toStr().data);
      if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
      return s;
   }
};


template <class T> inline
void wb_sptr<T>::new_dref_bare(size_t n, char exflag) {

   char alloc_new=(n>0);
   if (data) {
      if (!sptr_flags) wblog(FL,"ERR %s() %s",FCT,toStr().data);
      if (n!=sptr_flags->len || sptr_flags->nrefs) delete_dref();
      else alloc_new=0;
   }

#if 0
   printf("TST %s() thread #%d: 0x%lX | 0x%lX 0x%lX | %d (%d)\n",
   FCT, omp_get_thread_num(), this, sptr_flags, data, n, exflag);
   fflush(0);
#endif

   if (alloc_new) {
      WB_NEW(data,n);
      WB_NEW(sptr_flags,1);
      sptr_flags->len=n;
if (!sptr_flags->len) wblog(FL,"ERR %s() len=%d",FCT,sptr_flags->len);
   }
   else if (n && exflag) { T z=T();
      for (size_t i=0; i<n; ++i) data[i]=z;
   }

   if (n && (!sptr_flags || !data))
   wblog(FL,"ERR %s() 0x%lX 0x%lX",FCT,sptr_flags,data);
};


template <class T> inline
void wb_sptr<T>::new_dref(size_t n, const T *d0) {
   new_dref_bare(n,0);
   if (n) {
      if (d0) 
           MEM_CPY<T>(data,n,d0);
      else MEM_SET<T>(data,n);
   }
};


template <class T> inline
void wb_sptr<T>::init2ref(const T* d0) {

   if (data==d0) return;
   if (!d0) { delete_dref(); return; }

   if (sptr_flags && sptr_flags->isref) {
      data=(T*)d0;
   }
   else { delete_dref();
      data=(T*)d0;
      WB_NEW(sptr_flags,1);
      sptr_flags->len=size_t(-1);
      sptr_flags->isref=1;
   }
};


template <class T> inline
unsigned wb_sptr<T>::add_dref(wb_sptr &b) const {

   b.delete_dref();

   if (!data) return 0;

   b.data=(T*)data;
   if (!sptr_flags) wblog(FL,
      "ERR %s() got uninitialized sptr flags !??",FCT);


   if (sptr_flags->isref) {
      WB_NEW(b.sptr_flags,1);
      b.sptr_flags->len   = sptr_flags->len;
      b.sptr_flags->isref = 1;
      return 1;
   }


   b.sptr_flags=(sptr_flags__*)sptr_flags;

   if (!b.sptr_flags->len) wblog(FL,
      "ERR %s() len=%d",FCT,b.sptr_flags->len);



#pragma omp critical
 { ++(b.sptr_flags->nrefs); }

   if (b.sptr_flags->nrefs>10) wblog(FL,
      "TST %s() 0x%lx nrefs=%d",FCT, b.data, b.sptr_flags->nrefs);

   return b.sptr_flags->nrefs;
};


template <class T> inline
void wb_sptr<T>::delete_dref() {

   if (data) {
      if (char e=check()) wblog(FL,
         "ERR %s() got %s e=%d !??", FCT,toStr().data,e);

      #pragma omp critical
      { if (sptr_flags->isref) {
           WB_DELETE(sptr_flags);
           data=NULL;
        }
        else if (sptr_flags->nrefs) {
           --(sptr_flags->nrefs);
           sptr_flags=NULL; data=NULL;
        }
        else {
           WB_DELETE(data);
           WB_DELETE(sptr_flags);
        }
      }
   }
   else if (sptr_flags && sptr_flags->any())
   wblog(FL,"ERR %s() got %s !??",FCT,toStr().data);
};


template <class T> inline
void wb_sptr<T>::Instantiate(const char *F, int L, size_t n_) {

   if (!isRef() || !data) { wblog(F_L,
      "WRN %s() no need to instantiate (%s at 0x%lx)",
      FCT, isRef()? "ref":"noref", data); return; }
   if (!sptr_flags) wblog(F_L,
      "ERR %s() got missing sptr flags !??",FCT);

   T* d0=data; size_t n;

   if (sptr_flags->isref) {
      if (long(n_)<0) wblog(F_L,
         "ERR %s() isref requires size specification (%ld)",FCT,n_);
      n=n_;
   }
   else {
      if (int(sptr_flags->len)<=0) wblog(F_L,"ERR %s() "
         "invalid sptr flag (len=%d !?)",FCT,sptr_flags->len);
      n=sptr_flags->len;
   }
   new_dref(n,d0);
};


template <class T>
bool gotMemOverlap(
   const T* src, size_t len0,
   const T* dest, size_t len2
){
   const char
     *L=(char*)src,  *R=(char*)(src+len0)-1,
     *a=(char*)dest, *b=(char*)(dest+(long(len2)<0 ? len0 : len2))-1;

   if (!len0 || !len2) return 0;
   return ((a<=L && b>=L) || (a<=R && b>=R) || (a>=L && b<=R) );
};


#endif

