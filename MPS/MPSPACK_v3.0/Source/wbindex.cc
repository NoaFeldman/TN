#ifndef __WB_INDEX_CC__
#define __WB_INDEX_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
// tags: group collect start front end

int wbindex::extend2Perm(INDEX_T N, wbperm &P, const char *pstr) const {

    INDEX_T i,k,pflag=0;
    wbvector<char> mark(N);

    if (len>=N) {
       if (len>N) wblog(FL,
       "ERR %s - severe size mismatch (%d/%d)",len,N);
       return 0;
    }

    for (i=0; i<len; i++) {
       if (data[i]>=N) { sprintf(str,
          "%s:%d %s - index out of bounds (%ld/%ld)",FLF,data[i],N);
          return 1;
       }
       if (mark[data[i]]++) { sprintf(str,
          "%s:%d %s - index not unique",FLF);
          return 1;
       }
    }

    if (pstr && pstr[0]) {
       if (!strcmp(pstr,"end")) pflag=1;
       else wblog(FL,"ERR Invalid pstr=`%s'",pstr);
    }

    P.init(N);

    if (pflag==0) {
       memcpy(P.data,data,len*sizeof(INDEX_T));

       for (k=len, i=0; i<N; i++)
       if (!mark[i]) P[k++]=i;
    }
    else {
       memcpy(P.data+N-len,data,len*sizeof(INDEX_T));

       for (k=i=0; i<N; i++)
       if (!mark[i]) P[k++]=i;
    }

    return 0;
}


wbindex& wbindex::invert(INDEX_T n, wbindex &I, char uflag) const {

    wbvector<char> M(n); char *m=M.data;
    INDEX_T i,l;

    for (i=0; i<len; i++) { l=data[i];
        if (l<n) m[l]=1; else wblog(FL,
        "ERR %s() index out of bounds (%d: %d/%d)",FCT,i,l,n);
    }

    l=M.nnz();
    if (uflag && l!=len) wblog(FL,
       "ERR %s() index not unique (%d/%d)",FCT,l,len);

    I.init(n-l);
    if (I.len) {
       for (l=i=0; i<n; ++i) if (!m[i]) I.data[l++]=i;
    }

    return I;
};


wbindex wbindex::flipIdx(INDEX_T ndim) const {
    wbindex a(len);

    if (ndim==0) {
       if (len==0) return *this; else { dbstop(FL);
       wblog(FL,"ERR wbindex - invalid ndim=%d (%d)", ndim, len); }
    }

    for (INDEX_T m=ndim-1, i=0; i<len; i++) {
        if (i>m) wblog(FL,"ERR wbindex - invalid ndim=%d (%d)",ndim,i);
        a.data[i]=m-data[i];
    }
    return a;
};


inline char wbindex::isUnique(INDEX_T imax, INDEX_T offset) const {
    wbvector<char> mark(imax);

    if (offset==0) {
       for (INDEX_T i=0; i<len; i++) {
          if (data[i]<imax) { if ((++mark[data[i]])>1) return 0; }
          else return 0;
       }
    }
    else {
       INDEX_T i=0, i0=offset; imax+=offset;
       for (; i<len; i++) {
          if (data[i]>=i0 && data[i]<imax) {
             if ((++mark[data[i]-i0])>1) return 0;
          }
          else return 0;
       }
    }

    return 1;
};


inline wbindex& wbindex::init(const ctrIdx &ic) {
   initT(ic.len,ic.data);
   return *this;
};

inline wbindex& wbindex::init(
    const char *F, int L, const mxArray *a, INDEX_T offset
){
    if (a==NULL) { init(); return *this; }

    INDEX_T i, m=mxGetM(a), n=mxGetN(a);
    const double *d0;

    if (!mxIsDblMat(0,0,a) || (m>1 && n>1)) wblog(F_L,
       "ERR %s() numeric vector required (%s; %dx%d)",
        FCT,mxGetClassName(a),m,n
    );

    if (!m || !n) { init(); return *this; }

    n*=m; init(n); d0=mxGetPr(a);

    for (i=0; i<n; i++) { data[i]=INDEX_T(d0[i]);

       if (double(data[i])!=d0[i]) wblog(F_L,
          "ERR %s() type mismatch (%d: %d %g)",i+1,data[i],d0[i]);

       if (offset) {
          if (data[i]<offset) {
             if (offset==1) wblog(F_L,
                "ERR %s() 1-based index required (%d: %d)",
                 FCT,i+1,data[i]);
             else wblog(F_L,
                "ERR %s() index out of bounds (%d: %d/%d)",
                 FCT,i+1,data[i],offset
             );
          }
          data[i]-=offset;
       }
    }

    return *this;
};


template <class T>
template <class T2>
void groupIndex<T>::getRecs0(
   const wbMatrix<T2> &AB, wbMatrix<T2> &X, INDEX_T m
 ) const {

   INDEX_T i0=0; if (sINDEX_T(m)<0) { m=-m; i0=AB.dim2-m;  }
   
   if (sINDEX_T(i0)<0 || i0+m>AB.dim2) wblog(FL,
      "ERR %s() size out of bounds (%d+%d/%d)",FCT,i0+m,AB.dim2);

   X.init(D.len,m); if (!D.len || !m) return;

   if (I0.len!=D.len+1) wblog(FL,
      "ERR %s() unexpected size %d/%d",FCT,I0.len-1,D.len);
   if (AB.dim1!=P.len) wblog(FL,"ERR %s() size mismatch "
      "(%dx%d <> %dx%d)",FCT,P.len, m, AB.dim1, AB.dim2);

   for (INDEX_T i=0; i<D.len; ++i) {
      X.recSetP(i,AB.ref(I0.data[i],i0));
   }
};


int ctrIdx::init(const char *F, int L, const char *s) {

   unsigned i=0, r=0, nsep=0, l=32;
   char sflag=1;

   if (!s) { if (F) {
      wblog(F,L,"ERR %s() got null string",FCT); }
      return -1;
   }

   for (; s[i]; ++i) {
      if (CTR_IS_NUM(s[i])) {
         if (sflag) { ++r; sflag=0; }
      }
      else if ( CTR_IS_SEP_ANY(s[i])) { ++nsep; sflag=1; }
      else if (!CTR_IS_CONJ   (s[i])) { break; }
   }

   if (s[i] || i>l) { if (F) wblog(F,L,"ERR %s() got %s "
      "contraction set '%s'", FCT, i>l ? "unexpected":"invalid", s);
      return -2;
   }
   if (!r) { if (F) wblog(F,L,
      "ERR %s() got empty contraction set '%s' !?",FCT,s);
      return r;
   }

   conj=0; l=i;

   if (nsep==0) {
      if (CTR_IS_CONJ(s[l-1])) { conj=1; --l; }
      wbvector<unsigned>::init(l);
      for (i=0; i<l; ++i) {
         if (!CTR_IS_NUM(s[i])) { if (F) wblog(F,L,
            "ERR %s() unexpected contraction set '%s'",FCT,s);
            return -10;
         }
         data[i]=s[i]-'1';
      }
   }
   else {
      int e=Str2Idx(s, (wbvector<unsigned>&)(*this), 1);

      if (r!=this->len) { if (F) wblog(F,L,
         "ERR %s() got rank mismatch %d/%d !?\n%s",FCT,r,len,s);
         return -3;
      }

      if (e<=0) { i=-e-1;
         if (!e) { if (F) {
            wblog(F,L,"ERR %s() got null string !?",FCT); }
            return -4;
         }
         if (s[i++]!=';') { if (F) wblog(F,L,
            "ERR %s() invalid conj in '%s' (%d: %d/;)",FCT,s,i,s[i-1]);
            return -5;
         }

         while (s[i]==' ' && i<l) { ++i; }
         if (CTR_IS_CONJ(s[i])) { conj=1; ++i; }
         else { if (F) wblog(F_L,































































































