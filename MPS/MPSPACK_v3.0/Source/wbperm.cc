#ifndef __WB_PERM_CC__
#define __WB_PERM_CC__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

inline char
wbperm::isValidPerm(const char *F, int L, PERM_T r) const {

   if (sPERM_T(r)>=0 && r!=len) { if (F) wblog(F,L,
      "ERR invalid permutation [%s; %d]",toStr().data,r);
      return 0; 
   }

   if (len) { wbvector<char> mark(len);
      for (PERM_T i=0; i<len; ++i) {
         if (data[i]>=len || mark[data[i]]++) { if (F) wblog(F,L,
            "ERR invalid permutation [%s]",toStr().data);
            return 0;
         }
      }
   }

   return 1;
};


wbperm& wbperm::init(const wbperm &p1_, const wbperm &p2_) {

   if (!p1_.len) {
      if (p2_.len) return init(p2_);
      else return init();
   }
   if (!p2_.len) { return init(p1_,'i'); }
   
   if (this==&p1_ || this==&p2_) {
      wbperm X; X.init(p1_,p2_);
      X.save2(*this); return *this;
   }

   const PERM_T *p1=p1_.data, *p2=p2_.data;
   PERM_T i=0, n=p1_.len; NEW(n);

   if (p1_.len!=p2_.len) wblog(FL,
      "ERR %s() length mismatch (%d/%d)",FCT,p1_.len,p2_.len);

   for (; i<n; ++i) {
      if (p1[i]>=n || p2[i]>=n) wblog(FL,"ERR %s() invalid "
         "permutations (%d: %d,%d/%d)",FCT,i,p1[i],p2[i],n);
      if (data[p1[i]]) wblog(FL,"ERR %s() invalid "
         "permutations (%d -> %d: %d)",FCT,i,p1[i],data[p1[i]]);
      data[p1[i]]=p2[i];
   }

   return *this;
};


wbperm& wbperm::init2Front(const wbvector<PERM_T> &I, PERM_T N) {

   PERM_T i,j;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   init(N);
   for (i=0; i<I.len; i++) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N); 
      data[i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (i<N)
   for (j=0; j<N; j++) { if (!m[j]) data[i++]=j; }

   return *this;
};

wbperm& wbperm::init2End(const wbvector<PERM_T> &I, PERM_T N) {
   PERM_T i,j, l=N-I.len;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   init(N);
   for (i=0; i<I.len; i++) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N); 
      data[l+i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (l)
   for (i=j=0; j<N; j++) { if (!m[j]) data[i++]=j; }

   return *this;
};


wbperm& wbperm::init2FrontB(PERM_T m, PERM_T N) {

   if (!m || m>N) wblog(FL,
      "ERR %s() m out of bounds (%d/%d)",FCT,m,N); 

   PERM_T i,l=N-m; init(N);

   for (i=0; i<m; ++i) { data[i]=i+l; }
   for (   ; i<N; ++i) { data[i]=i-m; }

   return *this;
};

wbperm& wbperm::init2EndB(PERM_T m, PERM_T N) {

   if (!m || m>N) wblog(FL,
      "ERR %s() m out of bounds (%d/%d)",FCT,m,N); 

   PERM_T i,l=N-m; init(N);

   for (i=0; i<l; ++i) { data[i]=i+m; }
   for (   ; i<N; ++i) { data[i]=i-l; }

   return *this;
};

wbperm& wbperm::init2Front(PERM_T k, PERM_T N) {
   if (k>=N) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,N);
   PERM_T i=0; RENEW(N); data[0]=k;
   for (; i<k; ++i) { data[i+1]=i; }; ++i;
   for (; i<N; ++i) { data[i  ]=i; };
   return *this;
};

wbperm& wbperm::init2End(PERM_T k, PERM_T N) {
   if (k>=N) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,N);
   PERM_T i=0; RENEW(N); data[N-1]=k;
   for (; i<k; ++i) { data[i  ]=i; }; ++i;
   for (; i<N; ++i) { data[i-1]=i; };
   return *this;
};


wbperm& wbperm::Cycle(PERM_T k1, PERM_T k2) {
   if (k1>=len || k2>=len) wblog(FL,
      "ERR %s() index out of bounds [%d %d]/%d !?",FCT,k1,k2,len);
   if (k1==k2) return *this;

   PERM_T x1=data[k1], k=k1;
   if (k1<k2) 
        { for (; k<k2; ++k) { data[k]=data[k+1]; }}
   else { for (; k>k2; --k) { data[k]=data[k-1]; }}
   data[k]=x1;

   return *this;
};


wbperm& wbperm::initCycle(PERM_T r, PERM_T k1, PERM_T k2) {

   if (k1>=r || k2>=r) wblog(FL,
      "ERR %s() index out of bounds [%d %d]/%d !?",FCT,k1,k2,len);
   if (k1==k2) { return Index(r); }

   RENEW(r);

   if (k1<k2) { PERM_T k=0;
      for (; k<k1; ++k) { data[k]=k; }
      for (; k<k2; ++k) { data[k]=k+1; }; data[k++]=k1;
      for (; k<r;  ++k) { data[k]=k; }
   }
   else { PERM_T k=r-1;
      for (; k>k1; --k) { data[k]=k; }
      for (; k>k2; --k) { data[k]=k-1; }; data[k--]=k1;
      for (; k<r;  --k) { data[k]=k; }
   }

   return *this;
};


wbperm& wbperm::Rotate(sPERM_T k) {

    if (!len) wblog(FL,
      "ERR %s() got empty permutation %d/%d !?",FCT,k,len);
    if ((k%=len)<0) k+=len;
    if (k) {
       wbperm P0(*this);
       for (PERM_T i=0; i<len; ++i) { data[(i+k)%len]=P0[i]; }
    }
    return *this;
};


wbperm& wbperm::initRotate(PERM_T r, sPERM_T k) {

    if (!r) {
       if (k) wblog(FL,
          "ERR %s() got empty permutation %d/%d !?",FCT,k,r);
       return Index(r);
    }

    if ((k%=r)<0) k+=r;
    if (!k) return Index(r);

    RENEW(r);
    for (PERM_T i=0; i<len; ++i) { data[(i+k)%len]=i; }
    return *this;
};


wbperm& wbperm::Permute(const wbperm &P, char iflag, PERM_T r){

   if (P.isEmpty() || P.isIdentityPerm(FL,r)) {} else
   if (  isEmpty() ||   isIdentityPerm(FL,r)) {
      if (iflag) P.getIPerm(*this); else init(P);
   }
   else {
      wbperm P_(r);
      PERM_T i=0, *&p=P_.data;
      const PERM_T *p1=P.data, *p0=data;

      if (iflag)
           { for (; i<r; ++i) p[p1[i]]=p0[i]; }
      else { for (; i<r; ++i) p[i]=p0[p1[i]]; }

      P_.save2(*this);
   }
   return *this;
};


#endif


