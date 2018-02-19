#ifndef __WB_SIMPLE_VECTOR_CC__
#define __WB_SIMPLE_VECTOR_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

template <class T>
wbvector<T>& wbvector<T>::set(const wbvector<size_t> &I0, const T &x) {
   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      data[I[i]]=x;
   }
   return *this;
};

template <class T>
wbvector<T>& wbvector<T>::add(const wbvector<size_t> &I0, const T &x) {
   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      data[I[i]]+=x;
   }
   return *this;
};

template <class T>
bool wbvector<T>::anyUnequal(
   const wbvector<size_t> &I0, const T &x, size_t *k
 ) const {

   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      if (data[I[i]]!=x) { if (k) (*k)=I[i]; return 1; }
   }
   return 0;
};

template <class T>
bool wbvector<T>::anyEqual(
   const wbvector<size_t> &I0, const T &x, size_t *k
 ) const {

   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      if (data[I[i]]==x) { if (k) (*k)=I[i]; return 1; }
   }
   return 0;
};


template <class T>
wbvector<T>& wbvector<T>::Set(
   const wbvector<size_t> &I, const wbvector<T> &v
){
   if (!I.isUnique()) wblog(FL,"ERR index not unique");
   if (I.len!=v.len) wblog(FL,
      "ERR length mismatch (%d/%d)",I.len,v.len);

   for (size_t j,i=0; i<I.len; i++) { j=I[i];
      if (j>=len) wblog(FL,
      "ERR index out of bounds (%d: %d/%d)",i,j,len);

      data[j]=v[i];
   }
   return *this;
};


template <class T>
bool wbvector<T>::isNormal() const {
   wblog(FL,"ERR %s() not yet defined for type %s",FCT,
   getName(typeid(T)).data); return 0;
};


template <>
bool wbvector<double>::isNormal() const {
   for (size_t i=0; i<len; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <>
bool wbvector<float>::isNormal() const {
   for (size_t i=0; i<len; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <> bool wbvector<unsigned>::isNormal() const { return 1; };
template <> bool wbvector<size_t  >::isNormal() const { return 1; };
template <> bool wbvector<int     >::isNormal() const { return 1; };
template <> bool wbvector<char    >::isNormal() const { return 1; };
template <> bool wbvector<long    >::isNormal() const { return 1; };



template <class T>
wbvector<T>& wbvector<T>::RevertSigns() {

   size_t i=0;
   for (; i<len; ++i) { if (data[i]) break; }
   if (i<len) {
      if (data[i]>0 && (-data[i])>0) wblog(FL,"ERR %s() "
         "got unsigned data type %s",FCT,getName(typeid(T)).data);
      for (; i<len; ++i) data[i]=-data[i];
   }
   return *this;
};


template <class T>
void wbvector<T>::selectSU(const WBINDEX &I) {
    size_t i,n;
    WBINDEX mark(len);
    for (i=0; i<I.len; i++) {
       if (I[i]>=len)
       wblog(FL,"ERR Index out of bounds (%d/%d).", i,len);
       mark[I[i]]++;
    }
    for (n=i=0; i<mark.len; i++)
    if (mark[i]) mark[n++]=i; 

    mark.Resize(n);
    Select(mark);
};

template <class T>
wbvector<T>& wbvector<T>::select(
    const WBINDEX &I, wbvector<T> &a
) const {

    a.init(I.len); if (I.len==0) return a;

    for (size_t i=0; i<a.len; i++) {
       if (I[i]>=len) wblog(FL,
       "ERR Index out of bounds (%d/%d).", I[i],len);

       a.data[i]=data[I[i]];
    }
    return a;
};

template <class T>
void wbvector<T>::select(const WBINDEX &I, T* r) const {
    if (I.len==0) return;

    for (size_t i=0; i<I.len; i++) {
       if (I[i]>=len) wblog(FL,
       "ERR Index out of bounds (%d/%d).", I[i],len);

       r[i]=data[I[i]];
    }
};

template <class T>
wbvector<T>& wbvector<T>::Select(const WBINDEX &I) {

    if (isref) wblog(FL,"ERR must not resize vector reference!");
    size_t i; T* d0=data;

    if (I.len==0) { init(); return *this; }
    WB_NEW(data,I.len);

    for (i=0; i<I.len; i++) {
       if (I[i]>=len) wblog(FL,
          "ERR index out of bounds (%d/%d)", I[i],len);
       data[i]=d0[I[i]];
    }

    len=I.len;
    WB_DELETE(d0); return *this;
};

template <class T>
wbvector<T>& wbvector<T>::BlockSelect(
    const WBINDEX &I, const WBINDEX &D
){
    WBINDEX dc;
    size_t i,j,l,n; T *d0=data, *d;

    if (isref) wblog(FL,"ERR shall not resize vector reference!");
    if (I.len==0) { init(); return *this; }

    for (n=i=0; i<I.len; i++) {
       j=I[i]; if (j>=D.len) wblog(FL,
         "ERR index out of bounds (%d/%d).",j,D.len);
       n+=D[j];
    }

    i=D.cumsum0(dc);
    if (i!=len) wblog(FL,"ERR severe size mismatch (%d/%d; %d)",i,len,n);

    WB_NEW(data,n); len=n;

    for (n=l=i=0; i<I.len; ++i, l+=n) { j=I[i];
       d=d0+dc[j]; n=D[j];
       for (j=0; j<n; j++) data[l+j]=d[j];
    }

    WB_DELETE(d0); return *this;
};


template <class T>
wbvector<T>& wbvector<T>::Cat(const wbvector<wbvector<T> > &vv) {
   size_t i,n=0; T *d;

   for (i=0; i<vv.len; i++) n+=vv[i].len;
   init(n); d=data;

   for (i=0; i<vv.len; i++) {
      memcpy(d,vv[i].data,sizeof(T)*vv[i].len); d+=vv[i].len;
   }
   return *this;
};

template <class T>
wbvector<T>& wbvector<T>::Cat(const wbvector<wbvector<T>* > &vv) {
   size_t i,n=0; T *d;

   for (i=0; i<vv.len; i++) n+=(vv[i] ? vv[i]->len : 0);
   init(n); d=data;

   for (i=0; i<vv.len; i++) { if (vv[i]) {
       const wbvector<T> &v=*vv[i];
       memcpy(d,v.data,sizeof(T)*v.len); d+=v.len;
   }}
   return *this;
};


template <class T>
void wbvector<T>::info(const char *istr) const {

   size_t l=0, n=32; char sstr[n];
   size_t b=len*sizeof(T);
   wbstring tstr;

   if (typeid(T)==typeid(double)) tstr="double";
   else tstr.init(getName(typeid(T)).data);

   if (b<(1<<10)) l=snprintf(sstr,n,"%ld ",   b); else
   if (b<(1<<20)) l=snprintf(sstr,n,"%.3g kB",b/double(1<<10)); else
                  l=snprintf(sstr,n,"%.3g MB",b/double(1<<20));
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   printf("  %-12s %10d %12s  @ %p  %s vector%s\n",
   istr, len, sstr, (void*)data, tstr.data, isref ? "  ***ISREF***" : "");
}


template <class T>
void wbvector<T>::print(const char *istr, char mflag) const {

    mxArray *a;
    double *dd;

    if (mflag==1) {
        char fmt[8];

        if (typeid(T)==typeid(double))
             strcpy(fmt," %8g");
        else strcpy(fmt," %4d");

        if (istr[0]) printf("%s:",istr);

        for (size_t i=0; i<len; i++) printf(fmt,data[i]);
        printf("%s\n", isref ? " ***ISREF***" : "");

        return;
    }

    a=mxCreateDoubleMatrix(1,len,mxREAL);
    dd=mxGetPr(a);

    for (size_t i=0; i<len; i++)
    dd[i]=(double)data[i];

    if (!mflag) {
        if (istr[0])
        mexPrintf("\n%s = [1x%d double]\n\n", istr, len);
        else mexPrintf("\n");
    }
    else
    mexPrintf("\n%s = [%s\n", istr[0] ? istr : "ans",
      isref ? "  ***ISREF***" : "");

    mexCallMATLAB(0,NULL,1,&a, "disp");

    if (mflag>1) mexPrintf("];\n");
    mxDestroyArray(a);
    return;
}


template <class T> inline
wbindex& wbvector<T>::find(const T &x, wbindex &I, char iflag) const {

   if (!len) { I.init(); }
   else {
      size_t i=0, l=0;
      for (; i<len; ++i) { if (data[i]==x) ++l; }

      if (!iflag) { I.init(l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]==x) I[l++]=i; }}
      }
      else { I.init(len-l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]!=x) I[l++]=i; }}
      }
   }
   return I;
};


template <class T> inline
wbindex& wbvector<T>::find(const T &x, wbindex &I, wbindex &Ix) const {

   if (!len) { I.init(); Ix.init(); }
   else {
      size_t i=0, l=0, k=0;
      for (; i<len; ++i) { if (data[i]==x) ++l; }
      I.init(l); Ix.init(len-l); l=i=0;
      for (; i<len; ++i) { if (data[i]==x) I[l++]=i; else Ix[k++]=i; }
   }
   return I;
};


template <class T> inline
wbindex& wbvector<T>::findGT(const T &x, wbindex &I, char iflag) const {

   if (!len) { I.init(); }
   else {
      size_t i=0, l=0;
      for (; i<len; ++i) { if (data[i]>x) ++l; }

      if (!iflag) { I.init(l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]>x) I[l++]=i; }}
      }
      else { I.init(len-l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (!(data[i]>x)) I[l++]=i; }}
      }
   }
   return I;
};


template <class T> inline
wbindex& wbvector<T>::findGT(const T &x, wbindex &I, wbindex &Ix) const {

   if (!len) { I.init(); Ix.init(); }
   else {
      size_t i=0, l=0, k=0;
      for (; i<len; ++i) { if (data[i]>x) ++l; }; 
      I.init(l); Ix.init(len-l); l=i=0;
      for (; i<len; ++i) { if (data[i]>x) I[l++]=i; else Ix[k++]=i; }
   }
   return I;
};


template <class T> inline
size_t wbvector<T>::numelFind(const T &x) const {
   size_t n=0, i=0;
   for (; i<len; ++i) { if (data[i]==x) ++n; }
   return n;
};

template <class T> inline
size_t wbvector<T>::numelFindGT(const T &x) const {
   size_t n=0, i=0;
   for (; i<len; ++i) { if (data[i]>x) ++n; }
   return n;
};


template <class T> inline
wbindex& wbvector<T>::findRange(
   const T &x1, const T &x2, wbindex &I, const char *w
) const {

   if (len==0) { I.init(); return I; }

   bool i1=1,i2=1; size_t k=0;
   if (w) {
      if ((w[0]!='[' && w[0]!=']') || (w[1]!='[' && w[1]!=']') || w[2])
      wblog(FL,"ERR invalid interval specification '%s'",w);
      i1=(w[0]=='[');
      i2=(w[1]==']');
   }

   I.init(len);
   for (size_t i=0; i<len; i++) { const T &x=data[i];
      if (x<x1 || x>x2 || (x==x1 && !i1) || (x==x2 && !i2)) continue;
      I[k++]=i;
   }
   if (k) I.len=k; else I.init();

   return I;
};


template <class T>
size_t wbvector<T>::findValues(const wbvector<T> &B0, wbindex &Ia) const {

   size_t l=0,i=0,j=0, m=len, n=B0.len;
   char c;

   wbvector<T> A(*this), B(B0);
   wbperm P1,P2; wbindex I;
   T *a, *b;

   if (A.isEmpty() || B.isEmpty()) { Ia.init(); return 0; }

   A.Sort(P1);
   B.Sort(P2); I.init(m); a=A.data; b=B.data;

   while (i<m && j<n) { c=NUMCMP(a[i],b[j]);
      if (c<0) i++; else
      if (c>0) j++;
      else {
         do { I[l++]=i++; } while (i<m && a[i-1]==a[i]);
         do {        j++; } while (j<n && b[j-1]==b[j]);
      }
   }

   if (l) I.len=l; else I.init();

   P1.get(I,Ia);

   return Ia.len;
}


template <class T> inline
wbvector<T>& wbvector<T>::getI(size_t k, wbvector<T> &v) const {

    if (k>=len) wblog(FL,
       "ERR Index out of bounds (%d/%d) ???",k,len);

    v.init(len-1);
    for (size_t l=0, i=0; i<len; i++) if (i!=k) v.data[l++]=data[i];
    return v;
};


template <class T> inline
wbvector<T>& wbvector<T>::getI(
  const wbindex &I, wbvector<T> &v, char uflag
) const {

    wbindex J; I.invert(len,J,uflag);
    return select(J,v);
};


template <class T> inline
wbvector<T>& wbvector<T>::initI(
    const wbvector<T> &a, const wbindex &ia,
    const wbvector<T> &b, const wbindex &ib, char uflag
){
    if (!a.len && !b.len) {
       if (ia.len || ib.len) wblog(FL,
          "ERR %s() index out of bounds (%d+%d)",FCT,ia.len,ib.len);
       return init();
    }

    size_t i=0, j=0, na=a.len, nb=b.len, n=na+nb;
    char ma[n], *mb=ma+na; memset(ma,0,n*sizeof(char));

    if (uflag) {
       for (; i<ia.len; ++i) { j=ia[i];
          if (j>=na) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,a.len);
          if ((++ma[j])>1) wblog(FL,
             "ERR %s() index not unique (%d/%d)",FCT,j+1,a.len
          );
       }
       for (i=0; i<ib.len; ++i) { j=ib[i];
          if (j>=nb) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,b.len);
          if ((++mb[j])>1) wblog(FL,
             "ERR %s() index not unique (%d/%d)",FCT,j+1,b.len
          );
       }
       na-=ia.len; nb-=ib.len;
    }
    else {
       for (; i<ia.len; ++i) { j=ia[i];
          if (j>=a.len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,a.len);
          if ((++ma[j])==1) --na; else ma[j]=1;
       }
       for (i=0; i<ib.len; ++i) { j=ib[i];
          if (j>=nb) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,b.len);
          if ((++mb[j])==1) --nb; else mb[j]=1;
       }
    }

    init(na+nb);

    if (len) { j=0;
       if (na) {
          for (i=0; i<a.len; ++i) { if (!ma[i]) { data[j]=a[i]; ++j; } }
       }
       if (nb) {
          for (i=0; i<b.len; ++i) { if (!mb[i]) { data[j]=b[i]; ++j; } }
       }
       if (j!=len) wblog(FL,"ERR %s() %d/%d !??",FCT,j,len);
    }

    return *this;
};


template <class T> 
template <class T2>
wbvector<T>& wbvector<T>::initT(
   const char *F, int L, const wbvector<T2> &v) {

   RENEW(v.len);
   for (size_t i=0; i<len; ++i) {
      data[i]=T(v.data[i]);
      if (T2(data[i])!=v.data[i]) wblog(F_L,
         "ERR %s() got rounded value: %g (%s) -> %g (%s)",
         FCT, double(data[i]), getName(typeid(T2)).data,
         double(v.data[i]), getName(typeid(T)).data
      );
   }
   return *this;
};


template <class T> 
template <class TD> inline
wbvector<T>& wbvector<T>::initT(const wbsparray<TD> &S) {

   if (!S.isVector(FL)) wblog(FL,
      "ERR %s() invalid sparse vector (%s)",FCT,S.sizeStr().data);

   size_t i=0, j=0, n=S.IDX.dim1, m=S.IDX.dim2;
   const SPIDX_T *idx=S.IDX.data;

   for (; i<S.SIZE.len; ++i) { if (S.SIZE[i]>1) { idx+=i; break; }}

   RENEW(S.numel());

   for (i=0; i<n; ++i, idx+=m) { j=(*idx);
      if (j>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,len);
      data[j]=S.D.data[i];
   }

   return *this;
};


inline bool mxIsWbvector(
    const char *F, int L, const mxArray *a,
    size_t *n_, const char *istr, char dflag
){
    if (!a) return 0;
    size_t n=mxGetNumberOfElements(a); int i=0;
    if (!n) return 0;

    if (istr && istr[0])
         i=mxIsDblVector(0,0,a,0,dflag);
    else i=mxIsDblVector(F,L,a,0,dflag);

    if (!i) {
       if (F) wblog(F,L,
         "ERR %s() %s%snumeric vector required\n%s   (%s; %d)", FCT,
          istr?istr:"", istr?" ":"", str, mxGetClassName(a), n);
       return 0;
    }

    if (n_) (*n_)=n;
    return 1;
};


template <class T> inline
int wbvector<T>::init(
    const char *F, int L,
    const mxArray *a, const char *istr, char check_type
){
    size_t n=0;

    if (!a || !mxIsWbvector(F_L,a,&n,istr,'*')) {
       init(); return 1;
    }

    init(n);

    int q=cpyRangeMx(a,data,n,check_type);
    if (q<0) {
       sprintf(str,"%s%stype mismatch (%s; %d)", istr? istr:"",
          istr? " ":"", a ? mxGetClassName(a) : "null",q);
       if (F) wblog(F_L,"ERR %s() got %s",FCT,str);
       return 1;
    }

    return 0;
};

#ifdef USE_WB_MPFR




template <> inline
int wbvector<Wb::quad>::init(
    const char *F, int L,
    const mxArray *a, const char *istr, char check_type
){
    if (mxIsWbvector(0,0,a)) {
       wbvector<double> X; X.init(F_L,a,istr,check_type);
       this->initT(X); return 0;
    }

    wbarray<char> X; X.init(F_L,a, 0);
    if (X.isEmpty()) { init(); return 0; }
    if (X.SIZE.len!=2) wblog(FL,
       "ERR %s() expecting rank-2 char array (%s)",FCT,X.sizeStr().data);

    size_t i=0, m=X.SIZE[0], n=X.SIZE[1];
    unsigned base=(check_type<2 ? 30 : check_type);
    const char *sd=X.data;

    if (m==1) { i=m; m=n; n=i; i=0; }
    init(n);

    for (; i<n; ++i, sd+=m) {
       if (Wb::strnlen(sd,m)>=m) wblog(FL,
          "WRN %s() missing terminating null for string !??\n"
          "(%d: %d/%d)",FCT,i,strlen(sd),m
       );
       
       data[i].init_s(F_L,sd,base);
    }

    return 0;
};


template <> inline 
int wbvector< Wb::quad >::initMPFR(
   const char *F, int L, const mxArray *a, char base
){
   if (!a) {
      if (F) wblog(FL,"ERR %s() got null mxArray",FCT);
      return -1;
   };
   if (mxIsEmpty(a)) {
      if (mxIsChar(a) || mxIsDouble(a)) { init(); return 0; }
      if (F) wblog(FL,
         "ERR %s() got invalid type `%s'",FCT,mxGetClassName(a));
      init(); return -1;
   }
   if (!mxIsChar(a) && !mxIsNumeric(a)) {
      if (F) wblog(FL,
         "ERR %s() got invalid type `%s'",FCT,mxGetClassName(a));
      init(); return -1;
   }

   if (base==1 || base>62) wblog(FL,
      "ERR %s() invalid base=%d<%c>\n"
      "(interpreted as base for Wb::quad)",FCT,base,base); 

   size_t i,j,l, n=0; unsigned b=(base<=0 ? 30 : base);

   const mwSize r=mxGetNumberOfDimensions(a), *sp=mxGetDimensions(a);
   wbvector<size_t> S(r);

   for (i=0; i<r; ++i) S[i]=sp[i];

   if (r!=2 || S[0]>62 || S[0]<1) {
      if (F) wblog(FL,
         "ERR invalid input array (%s) !??",S.toStrf("","x").data);
      init(); return 1;
   }

   if (S[0]!=1)
        { n=S[0]; init(S[1]); }
   else { n=S[1]; init(S[0]); }

   wbarray<char> ss_(S[0],S[1]); char *ss=ss_.data;
   cpyRangeMx(a,ss,S[0]*S[1]);

   char s[n+1];

   for (i=0; i<len; ++i) {
      for (j=0; j<n; ++j) { if ((s[j]=ss[j])!=' ') break; }; l=j;
      for (   ; j<n; ++j) { s[j]=ss[j];
         if (!s[j] || s[j]==' ') { break; }
      }; s[j]=0; ss+=n;

      if (l>=j) wblog(FL,"ERR %s() got empty string (%d/%d)",FCT,l,j);
      
      try {
         data[i].init_s(F_L,s+l,b);
      }
      catch (...) {
         wblog(FL,"ERR %s() `%s'",FCT,s);
      }
   }

   return 0;
};

#endif


template <class T>
inline wbvector<T>& wbvector<T>::Permute(const wbperm &P, char iflag) {
   const PERM_T *p=P.data;

   if (P.isEmpty()) return *this;

   if (len!=P.len || !validPerm(P)) { if (P.len<12) wblog(FL,
      "ERR invalid permutation (%d/%d)\n[%s]",
       len, P.len, P.toStr().data); else wblog(FL,
      "ERR invalid permutation (%d/%d)", len, P.len); }

   if (!P.isIdentityPerm()) {
      wbvector<T> v(*this);
      if (iflag==0) {
             for (size_t i=0; i<len; ++i) data[i]=v.data[p[i]]; }
      else { for (size_t i=0; i<len; ++i) data[p[i]]=v.data[i]; }
   }

   return *this;
};


template <class T>
inline wbvector<T>& wbvector<T>::permute(
   wbvector<T> &v, const wbperm &P, char iflag
 ) const {

   if (P.isEmpty()) { v=*this; return v; }
   const PERM_T *p=P.data;

   if (len!=P.len || !validPerm(P)) {
      if (P.len<12) wblog(FL,
      "ERR invalid permutation (%d/%d)\n[%s]",len,P.len,P.toStr().data);
      else wblog(FL,
      "ERR invalid permutation (%d/%d)", len, P.len);
   }

   if (!P.isIdentityPerm()) { v.init(len);
      if (iflag==0) {
             for (size_t i=0; i<len; i++) v.data[i]=data[p[i]]; }
      else { for (size_t i=0; i<len; i++) v.data[p[i]]=data[i]; }
   }
   else v=(*this);
   
   return v;
};


template <class T> inline
wbvector<T>& wbvector<T>::blockPermute(
   const wbperm &P, wbvector<T> &v, char iflag
) const {

   if (isEmpty() || P.isEmpty() || P.isIdentityPerm()) {
       if (P.len && len%P.len) wblog(FL,
          "ERR %s() data mismatch (%d mod %d !?)",FCT,len,P.len); 
       v=*this; return v;
   }
   if (&v==this) {
      wbvector<T> x(*this);
      return x.blockPermute(P,v,iflag);
   }

   if (len%P.len || !validPerm(P)) wblog(FL,
      "ERR invalid permutation (%d/%d)", len, P.len);
   v.init(len);

   size_t i,j, m=len/P.len;
   const PERM_T *p=P.data; const T *d0=data; T *d=v.data;

   if (iflag==0) {
      for (i=0; i<P.len; ++i, d+=m) { d0=data+m*p[i];
      for (j=0; j<m; ++j) d[j]=d0[j]; }
   }
   else {
      for (i=0; i<P.len; ++i, d0+=m) { d=v.data+m*p[i];
      for (j=0; j<m; ++j) d[j]=d0[j]; }
   }

   return v;
};


template <class T> inline
int wbvector<T>::blockCompare(
   const wbperm &Pa, const wbvector &B, const wbperm &Pb
 ) const {

   if (len!=B.len ||
      (Pa.len &&   len%Pa.len) || (Pb.len && B.len%Pb.len))
      wblog(FL,"ERR %s() got length mistmatch (%d@%d, %d@%d)",
      FCT,len,Pa.len,B.len,Pb.len
   );

   if ((!Pa.len && !Pb.len) || Pa.len==len || Pb.len==B.len) {
      if (Pa.len) {
         if (Pb.len) {
            if (Pa.len!=Pb.len) wblog(FL,
               "ERR %s() got length mistmatch (%d@%d, %d@%d)",
               FCT,len,Pa.len,B.len,Pb.len
            );

            const PERM_T*pa=Pa.data, *pb=Pb.data;
            for (size_t i=0; i<len; ++i) {
               if (data[pa[i]]!=B.data[pb[i]]) {
                  return (data[pa[i]]<B.data[pb[i]] ? -1 : +1);
               }
            }
         }
         else {
            const PERM_T *pa=Pa.data;
            for (size_t i=0; i<len; ++i) {
               if (data[pa[i]]!=B.data[i]) {
                  return (data[pa[i]]<B.data[i] ? -1 : +1);
               }
            }
         }
      }
      else {
         if (Pb.len) {
            const PERM_T *pb=Pb.data;
            for (size_t i=0; i<len; ++i) {
               if (data[i]!=B.data[pb[i]]) {
                  return (data[i]<B.data[pb[i]] ? -1 : +1);
               }
            }
         }
         else {
            for (size_t i=0; i<len; ++i) {
               if (data[i]!=B.data[i]) {
                  return (data[i]<B.data[i] ? -1 : +1);
               }
            }
         }
      }
   }
   else {
      if (Pa.len) {
         if (Pb.len) {
            const PERM_T *pa=Pa.data, *pb=Pb.data;
            size_t j, i=0, m=len/Pa.len;"rank"
            const T *a, *b;

            if (Pa.len!=Pb.len) wblog(FL,
               "ERR %s() got invalid permutations (%d@%d, %d@%d)",
               FCT,len,Pa.len,B.len,Pb.len
            );

            for (; i<Pa.len; ++i) { a=data+m*pa[i]; b=B.data+m*pb[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) return (a[j]<b[j] ? -1 : +1);
               }
            }
         }
         else {
            const PERM_T *pa=Pa.data;
            size_t j, i=0, m=len/Pa.len;
            const T *a, *b=B.data;

            for (; i<Pa.len; ++i, b+=m) { a=data+m*pa[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) { return (a[j]<b[j] ? -1 : +1); }
               }
            }
         }
      }
      else {
         if (Pb.len) {
            const PERM_T *pb=Pb.data;
            size_t j, i=0, m=B.len/Pb.len;
            const T *a=data, *b;

            for (; i<Pb.len; ++i, a+=m) { b=B.data+m*pb[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) { return (a[j]<b[j] ? -1 : +1); }
               }
            }
         }
      }
   }

   return 0;
};


template <class T>
inline int wbvector<T>::set2Group(
   const wbperm &P, const WBINDEX &D,
   const T* S0,
   char iflag
){
   size_t i,j,l,d,n=D.len; int e=0;
   wbvector<T> X;

   if (!S0 && len!=P.len) wblog(FL,
   "ERR %s() severe size mismatch (%d/%d)", FCT, len, P.len);

   if (n==0) {
      if (P.len) wblog(FL,
      "ERR empty object inconsistency (%d,%d)", P.len, D.len);
      return 0;
   }

   if (S0==NULL || S0==data) { save2(X); S0=X.data; }

   init(n);

   for (l=i=0; i<n; i++, l+=d) { d=D[i];
      data[i]=S0[P[l]];

      for (j=1; j<d; j++) if (S0[P[l+j]]!=data[i]) {
         e++; if (iflag)
         sprintf(str,"%g/%g", (double)data[i], (double)S0[P[l+j]]);
      }
   }

   return e;
}


template<class T>
mxArray* wbvector<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbvector<T>::add2MxStruct(mxArray *S, unsigned i, char tflag) const {
   mxSetCell(S,i,toMx(tflag));
}


template<class T>
mxArray* wbvector<T>::mxCreateCell(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbvector<T>::add2MxCell(mxArray *S, unsigned i, char tflag) const {
   mxSetCell(S,i,toMx(tflag));
}


template<class T> inline
void wbvector<T>::add2MxStruct(
    mxArray *S, const char *vname, char tflag
) const {

    double *dd; int fid;
    mxArray *a;

    fid = mxAddField2Scalar(FL,S,vname);

    if (!tflag)
         a=mxCreateDoubleMatrix(1,len,mxREAL);
    else a=mxCreateDoubleMatrix(len,1,mxREAL);

    dd=mxGetPr(a);

    for (size_t i=0; i<len; i++) dd[i]=(double)data[i];

    mxSetFieldByNumber(S,0,fid, a);
};


template<class T> inline
wbstring wbvector<T>::toStrD() const { return toStrf("","x"); }

template<class T> inline 
wbstring wbvector<T>::toStr(int n, const char *sep) const {

    wbstring s(MAX((size_t)128, len*MAX(16,int(n)+6))), fmt;
    fmt.init2Fmt((T)0,n);

    for (size_t i=0; i<len; ++i) {
        if (i>0) { s.push(FL,sep); }
        s.pushf(FL,fmt.data,data[i]);
    }

    return s;
};

#ifdef USE_WB_MPFR

template<> inline 
wbstring wbvector<Wb::quad>::toStr(int n, const char *sep) const {

    wbstring s(MAX((size_t)128, len*MAX(16,int(n)+6))), fmt;
    fmt.init2Fmt((double)0,n);

    for (size_t i=0; i<len; ++i) {
        if (i>0) { s.push(FL,sep); }
        s.pushf(FL,fmt.data,double(data[i]));
    }

    return s;
};

#endif

template<class T> inline 
wbstring wbvector<T>::toStrf (
    const char *fmt0,
    const char *sep,
    const char stride,
    const char *sep2
) const {

    wbstring s(MAX((size_t)128, len*16)), fmt;
    char flag=1;

    if (fmt0==NULL ? 1 : (fmt0[0] ? 0 : 1)) {
        fmt.init2Fmt(T(0));
    } else fmt=fmt0;

    for (size_t i=0; i<len; ++i) {
       if (i>0) { if (flag) { s.push(FL,sep); } else { flag=1; }}
       s.pushf(FL,fmt.data,data[i]);
       if (stride && ((i+1)%stride)==0 && i+1<len) {
          s.push(FL,sep2); flag=0;
       }
    }

    return s;
};

#ifdef USE_WB_MPFR

template<> inline 
wbstring wbvector<Wb::quad>::toStrf (const char *fmt0,
   const char *sep, const char stride, const char *sep2
 ) const {

   wbvector<double> x; x.init(*this);
   return x.toStrf(fmt0,sep,stride,sep2);
};

#endif


#ifdef __WB_MPFR_HH__

template <> inline 
mxArray* wbvector< Wb::quad >::toMx(const char base) const {

   if (base==1 || base>62) wblog(FL,
      "ERR %s() invalid base=%d<%c>\n"
      "(interpreted as base for Wb::quad)",FCT,base,base); 

   size_t i,j, n=0, l=0; unsigned b=(base<=0 ? 30 : base);

   mwSize S[2];
   mxArray *a;

   wbstring s_; char *s=NULL;
   if (len) { data[0].toStr(s_,b); n=s_.len; s=s_.data; }

   if (len!=1)
        { S[0]=n; S[1]=len; }
   else { S[0]=len; S[1]=n; }

   if (len==1) {
      int k=double(data[0]);
      if (data[0]==Wb::quad(k)) {
         for (j=0; j<n && s[j]; ++j) { if (s[j]=='.') break; }
         if (j<n && s[j]) { l=j;
            for (++j; j<n && s[j]; ++j) { if (s[j]!='0')
               wblog(FL,"ERR %s() expecting integer (%s) !?",FCT,s);
            }; S[1]=l+1;

            a=mxCreateCharArray(2,S);
            if (!a) wblog(FL,"ERR %s() "
               "failed to allocate mxArray (%dx%d)",FCT,S[0],S[1]);
            mxChar* d=(mxChar*)mxGetData(a);
            for (j=0; j<l; ++j) { d[j]=s[j]; }; d[j]=0;

            return a;
         }
      }
   }

   wbarray<char> ss_(n,len); char *ss=ss_.data;
   size_t nmax=0;

   for (i=0; i<len; ++i, ss+=n) {
      if (i) { data[i].toStr(s_,b); s=s_.data;
         if (s_.len!=n) wblog(FL,
            "ERR %s() length of string changed (%d/%d)",FCT,s_.len,n
         );
      }

      for (j=0; j<n && s[j]; ++j) { if (s[j]!=' ') break; }

      if (s[j]=='-')
           { ss[0]=s[j]; }
      else { ss[0]=' '; --s; }

      for (s+=j, j=1; j<n && s[j]; ++j) { ss[j]=s[j]; }

      if (s[j]) wblog(FL,
         "ERR %s() string out of bounds (%d) !?",FCT,n); 
      if (nmax<(++j)) nmax=j;
   }

   if (nmax>n) wblog(FL,"ERR %s() got nmax=%d/%d",FCT,nmax,n);
   if (len!=1) S[0]=nmax; else S[1]=nmax;

   a=mxCreateNumericMatrix(S[0],S[1], mxINT8_CLASS, mxREAL);

   if (!a) wblog(FL,
      "ERR %s() failed to allocate mxArray (%dx%d)",FCT,S[0],S[1]);
   if (!len) return a;

   int8_T* dd = (int8_T*)mxGetData(a); ss=ss_.data;

   for (i=0; i<len; ++i, ss+=n, dd+=nmax) {
      for (j=0; j<nmax && ss[j]; ++j) { dd[j]=ss[j]; }
      for (l=j; j<nmax; ++j) { dd[j]=0;  }

      if (l>=nmax) wblog(FL,
         "ERR %s() no terminating null (Wb::quad string: %d/%d)\n"
         "'%s' (%d)",FCT,l,nmax,ss,strlen(ss)
      );
   }

   return a;
};

#endif


template <class T>
inline wbvector<T>& wbvector<T>::Move2_FE(
    const WBINDEX &I,
    char iflag
){
    size_t i,k;
    wbvector<char> mark(len);
    wbvector<T> x;

    if (I.isEmpty()) return *this;
    if (!isUniqueIdxSet(I,len)) {
        wblog(FL, "ERR Move2_FE Index out of bounds [%s; %d]",
        I.toStr().data, len); return *this;
    }

    for (i=0; i<I.len; i++) mark[I[i]]++;

    select(I,x);

    if (iflag=='E') {
        for (k=i=0; i<len; i++)
        if (mark[i]==0) {
            if (k!=i) data[k]=data[i];
            k++;
        }
        memcpy(data+k, x.data, x.len*sizeof(T));
    }
    else
    if (iflag=='F') {
        for (k=i=len-1; int(i)>=0; i--)
        if (mark[i]==0) {
            if (k!=i) data[k]=data[i];
            k--;
        }
        memcpy(data, x.data, x.len*sizeof(T));
    }
    else wblog(FL,"ERR Move2_FE - Invalid flag %c<%d>", iflag, iflag);

    return *this;
} 


template <class T>
inline wbvector<T>& wbvector<T>::Skip(size_t i1) {

    if (i1>=len) {
        wblog(FL, "ERR Skip - index out of range (%d,%d)\n[%s; %d]",
        i1, len, toStr().data, i1); return *this;
    }

    for (size_t i=i1+1; i<len; i++) data[i-1]=data[i];
    Resize(len-1);

    return *this;
}

template <class T>
inline wbvector<T>& wbvector<T>::Skip(
    const WBINDEX &I,
    char rflag
){
    size_t i,k;
    wbvector<char> mark(len); char *m=mark.data;

    for (k=0; k<I.len; k++) {
        if (I[k]<len) m[I[k]]++;
        else {
            wblog(FL, "ERR Skip - index out of range (%d,%d)\n[%s; %d]",
            I[k]+1, len, toStr().data, k+1); return *this;
        }
    }

    for (k=i=0; i<len; i++) if (m[i]==0) {
        if (k!=i) data[k]=data[i]; k++;
    }

    if (!rflag) Resize(k);
    else {
       if (k) len=k; else init();
    }

    return *this;
}


#endif

