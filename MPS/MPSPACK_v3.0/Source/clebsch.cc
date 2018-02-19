#ifndef __WB_CLEBSCH_GORDAN_CC__
#define __WB_CLEBSCH_GORDAN_CC__

// ------------------------------------------------------------------ //
// See also header file clebsch.hh.
// Wb,Sep20,09
// ------------------------------------------------------------------ //

//====================================================================//
// QType :: most basic symmetry unit                                  //
//====================================================================//

inline bool QType::validType(const char *F, int L) const {

   if (type<=0 || type==QTYPE_ASEP) {
      if (F) wblog(F,L,
         "ERR %s() unspecified symmetry (%d/%d; %s)",
         FCT, type, QTYPE_NUM_TYPES, toStr().data);
      return 0;
   }
   if (type>=QTYPE_NUM_TYPES) {
      if (F) wblog(F,L,
         "ERR %s() type out of bounds (%d/%d)",FCT,type,QTYPE_NUM_TYPES);
      return 0;
   }

   if (type<QTYPE_ASEP) {
      if (type==QTYPE_ZN) { if (sub>=2) return 1;
         if (F) wblog(F,L,
            "ERR %s() invalid sub=%d for %s [>=2]",FCT,sub,QTYPE_STR[type]);
         return 0;
      }
      if (!sub) return 1; 
      if (F) wblog(F,L,
         "ERR %s() invalid sub=%d for %s [0]",FCT,sub,QTYPE_STR[type]);
      return 0;
   }

   if (type>QTYPE_ASEP && sub>=1) return 1; 

   if (F) wblog(F,L,
      "ERR %s() invalid sub=%d for %s [>=2]",FCT,sub,QTYPE_STR[type]);

   return 0;
};



int QType::init_s(const char *s) {

   int k=1, e=0;

   if (!s || !s[0]) { type=QTYPE_UNKNOWN; sub=0; return 1; }

   for (; k<QTYPE_NUM_TYPES; ++k) {
      if (!strcmp(s,QTYPE_STR[k])) {
         type=QTYPE_SET(k); sub=0;
         validType(FL); return 0;
      }
   }

   if (s[0]=='Z') {
      if ((e=this->atoi(s+1,k,3))>0 && k>1) {
         type=QTYPE_ZN; sub=k;
      } 
      else return (-60+e);
   }
   else if (s[0]=='S' && s[1]=='U') {
      if ((e=this->atoi(s+2,k,2))>0 && k>1) {
         type=QTYPE_SUN; sub=k-1;
      } 
      else return (-70+e);
   }
   else if (s[0]=='S' && s[1]=='p') {
      if ((e=this->atoi(s+2,k,2))>0 && k>=4 && !(k%2)) {
         type=QTYPE_SpN; sub=k/2;
      } 
      else return (-80+e);
   }
   else return -90;

   validType(FL); return 0;
};


template <class TQ> inline
unsigned QType::qdim(const TQ* q) const {

   if (isAbelian()) { return 1; }

   unsigned d=0;
   QSet<TQ> Q; { Q.init2(*this,q); }

   CDATA_TQ *C = gCG.find_cdata(Q);
   if (C && !C->isEmpty()) { d=C->cgd.cgsparray::dim(); }

   if (!d) {
      genRG_base<TQ,RTD> *R=gRG.find_RSet(*this,q);
      if (R && R->Z.dim2) d=R->Z.dim1;

   if (!d) {
      C=&gCG.getBUF(0,0,Q,'l');
      if (C->cgd.SIZE.len) wblog(FL,"ERR %s() "
         "got unexpected scalar CData\n%s",FCT,C->toStr().data);
      d=C->cgd.D.len;

      if (!d) {
         const genRG_base<TQ,RTD> &R=gRG.getR(FL,*this,q);
         if (R.Z.dim2) d=R.Z.dim1;
      }

   if (!d) {
      if (isSU2()) {
         int n=q[0];
         if (n<0 || double(n)!=q[0]) wblog(FL,
            "ERR invalid %s symmetry q=%g",toStr().data,double(q[0]));
         d=(n+1);
      }
      else wblog(FL,
        "ERR %s() failed to determine multiplet dimension for\n%s",
         FCT,Q.toStr().data
      );
   }}}

#ifdef NSAFEGUARDS
   if (isSU2()) {
      if (d!=unsigned(q[0]+1)) wblog(FL,
         "ERR %s() unexpected multiplet dimension (%d/%d)",
         FCT,d,q[0]+1,Q.toStr().data
      );
   }
#endif

   return d;
};


template <class TQ> inline
wbvector<unsigned>& QType::QDim(
   const TQ* q, unsigned r, unsigned stride, wbvector<unsigned> &S
 ) const {

   S.init(r);

   if (isAbelian()) { S.set(1); } else
   for (unsigned i=0; i<r; ++i, q+=stride) { S[i]=this->QDim(q); }

   return S; 
};


inline mxArray* QType::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
};


void QType::add2MxStruct(mxArray *S, unsigned i) const {
   if (!S) wblog(FL,"ERR");
   mxSetCell(S,i,toStr('t').toMx());
};



QVec& QVec::init(
   const char *F, int L, const char *s,
   unsigned D0
){

   if (!s || !s[0] || !strcmp(s,CGC_ALL_ABELIAN)) {
      init(); return *this;
   }

   unsigned l,i, n=strlen(s), i1=0, i2=0;
   unsigned D=(D0 ? 2*D0 : 8);
   char sep=0, sx[n+1];
   QType q;

   init(D);

   for (l=i=0; i<=n; i++) {
      if (s[i]==' ' || s[i]=='\t' || s[i]=='\n' ||
          s[i]==',' || s[i]==';') { sep++; continue; }

      if (sep || s[i]==0) {
         if (l>=D) { l++; break; }

         sx[i2]=0; q.init(F,L,sx+i1);
         data[l++]=q; i1=i2; sep=0;
      }
      sx[i2++]=s[i];
   }

   if (l>D) {
      wblog(F,L,"ERR invalid qtype specs: %s\n"
     "(string contains too many entries; %d/%d)",s,l,D);
   }

   if (l) len=l; else init();

   if ((D0 && D0!=l) || l>Qlen()) wblog(F,L,
   "ERR qtype length inconsistency (%d/%d; %d)",D0,l,Qlen());

   return *this;
};


unsigned QVec::Qlen() const {
   unsigned n=0, i=0;
   for (; i<len; ++i) { n+=data[i].qlen(); }
   return n;
};

unsigned QVec::Qlen(wbvector<unsigned> &dd) const {
   unsigned i=0, n=0;
   dd.init(len); for (; i<len; ++i) { n+=(dd[i]=data[i].qlen()); }
   return n;
};

unsigned QVec::Qlen(wbvector<unsigned> &dd, wbvector<unsigned> &dz) const {
   unsigned i=0, n=0;
   dd.init(len); for (   ; i<len; ++i) { n+=(dd[i]=data[i].qlen()); }
   dz.init(len); for (i=0; i<len; ++i) { dz[i]=data[i].qrank(); }
   return n;
};

unsigned QVec::Qpos(wbvector<unsigned> &dc) const {

   dc.init(len); if (len) { unsigned i=1;
      for (; i<len; ++i) {
         dc[i] = dc[i-1] + data[i-1].qlen(); 
      } 
      return (dc[i-1] + data[i-1].qlen());
   }
   return 0;
};


unsigned QVec::Qrank() const {
   unsigned n=0, i=0;
   for (; i<len; ++i) { n+=data[i].qrank(); }
   return n;
};


template <class TQ>
wbvector<unsigned>& QVec::QDim(
   const TQ *q, wbvector<unsigned> &S
 ) const {
   
   unsigned d0, i=0; S.init(len);

   for (; i<len; ++i, q+=d0) {
      d0=data[i].qlen();
      S[i]=data[i].QDim(q);
   }

   return S;
};


template <class TQ>
wbvector<unsigned>& QVec::QDim(
   const TQ *q, unsigned k, unsigned r, wbvector<unsigned> &S
 ) const {
   
   unsigned i, d0, n=0;

   if (k>=len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,len);

   if (!r) { S.init(); return S; }
   if (data[k].isAbelian()) { S.init(1).set(1); return S; }

   for (i=0; i<len; ++i, n+=d0) {
      d0=data[i].qlen();
      if (i==k) q+=n;
   }

   if (S.len!=r) S.init(r);
   for (i=0; i<r; i++, q+=n) { S[i]=data[k].QDim(q); }

   return S;
};


template <class TQ>
unsigned QVec::QDim(const TQ *q) const {

   unsigned i,d0, n=(len ? 1 : 0);

   for (i=0; i<len; ++i) {
      d0=data[i].qlen();
      n*=data[i].qdim(q); q+=d0;
   }

   return n;
};


template <class TQ>
unsigned QVec::QDim(const TQ *q, unsigned k) const {

   for (unsigned d0=0, i=0; i<len; ++i, q+=d0) {
      d0=data[i].qlen();
      if (i==k) return data[i].QDim(q); 
   }

   wblog(FL,"ERR index out of bounds (%d/%d)",k,len); 
   return 0;
};


template <class TQ>
wbMatrix<TQ>& QVec::getQsub(
   const wbMatrix<TQ>& Q0, const wbindex &I,
   wbMatrix<TQ>& QI,
   wbMatrix<TQ> *Qx
 ) const {

   if (Q0.isEmpty()) { QI.init(); if (Qx) Qx->init(); return QI; }

   unsigned i,j,k,l,n;
   wbvector<unsigned> s0,s,ic0,ic;
   const TQ *d0=Q0.data; TQ *d;
   
   Qlen(s0); n=s0.sum();

   if (Q0.dim2!=n) wblog(FL,
      "ERR %s() Q-set size mismatch (%dx%d/%d)",FCT,Q0.dim1,Q0.dim2,n);

   s0.cumsum0(ic0); ic0.select(I,ic);
   s0.get(I,s); QI.init(Q0.dim1,s.sum()); d=QI.data;

   for (i=0; i<Q0.dim1; ++i, d0+=Q0.dim2)
   for (j=0; j<I.len; ++j, d+=n) {
      for (n=s[j], k=ic[j], l=0; l<n; ++l) d[l]=d0[k+l];
   }

   if (Qx) {
      wbindex Ix; I.invert(len,Ix);
      s0.get(Ix,s); Qx->init(Q0.dim1,s.sum()); d=Qx->data;
      ic0.select(Ix,ic); d0=Q0.data;

      for (i=0; i<Q0.dim1; ++i, d0+=Q0.dim2)
      for (j=0; j<Ix.len; ++j, d+=n) {
         for (n=s[j], k=ic[j], l=0; l<n; ++l) d[l]=d0[k+l];
      }
   }

   return QI;
};


wbstring QVec::toStr(const char vflag) const {

   wbstring S;

   if (!len) { 
      if (vflag!='V')
           S=(vflag ? CGC_ALL_ABELIAN : "");
      else S="(all abelian)";
      return S;
   }

   char tflag=(vflag ? 0:1);
   S.init(8*len+1); S.data[0]=0;

   for (unsigned k=0, i=0; i<len; i++) {
      if (i) { S.cat(k,tflag? ",":"*"); }
      S.cat(k, data[i].toStr(tflag).data);
   }
   return S;
};



QDir& QDir::init(const IndexLabels& it) {

   init2val(it.len,+1);
   for (unsigned i=0; i<len; ++i) {
      if (it[i].isConj()) data[i]=-1;
   }
   return *this;
};


QDir& QDir::init(const char *F, int L,const char *s) {

   if (!s || !s[0]) { init(); }
   else {
      unsigned i=0, l=(s? strlen(s): 0);
      if (l>32) wblog(F_L,
         "ERR %s() unexpected QDir string '%s' (len=%d)",FCT,s?s:"",l);
      for (init(l); i<l; ++i) {
         if (s[i]=='+') { data[i]=+1; } else
         if (s[i]=='-') { data[i]=-1; } else
         if (s[i]=='0') { data[i]= 0; } else
         wblog(F_L,"ERR QDir::%s() invalid string `%s'",FCT,s);
      }
   }
   return *this;
};


QDir& QDir::init_iout(const char *F, int L,
   unsigned r, unsigned iout
){
   if (!r) {
      if (iout) wblog(F_L,
         "ERR %s() invalid iout=%d (having r=%d)",FCT,iout,r); 
      return init();
   }
   if (iout) {
      if ((--iout)>=r) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,iout+1,r); 
      init2val(r,+1); data[iout]=-1;
   }
   else { init(r); }

   return *this;
};


bool QDir::operator==(const char *s) const {

   if (!s) {
      if (!len) return 1;
      wblog(FL,"ERR %s() got null string",FCT); 
   }

   unsigned i=0; char q=1;
   for (; i<len && s[i]; ++i) {
      if (s[i]=='+') { if (data[i]<=0) { q=0; ++i; break; }} else
      if (s[i]=='-') { if (data[i]>=0) { q=0; ++i; break; }} else
      { q=-1; break; }
   }
   if (q>=0) { for (; i<len && s[i]; ++i) {
      if (s[i]!='+' && s[i]!='-') { q=-1; break; }
   }}

   if (q<0) wblog(FL,
      "ERR QDir::%s() invalid qdir='%s' (%d/%d)",FCT,s,i,len);
   return (i<len || s[i] ? 0 : q);
};


wbstring QDir::toStr() const {
   wbstring s(16); unsigned i=0, l=0, n=s.len;
   for (; i<len && l<len; ++i) {
      if (data[i]==+1) { s[i]='+'; ++l; } else
      if (data[i]==-1) { s[i]='-'; ++l; } else
      if (data[i]== 0) { s[i]='0'; ++l; } else
      l+=snprintf(s.data+l,n-l,"<%+d>",data[i]);
   }
   if (l>=n) { s.data[n-1]=0; wblog(FL,"ERR %s() "
      "string out of bounds (%d/%d)\n`%s...'",FCT,l,n,s.data);
   }
   return s;
};


wbstring QDir::toTag() const { wbstring s(len+1);
   for (unsigned i=0; i<len; ++i) {
      if (data[i]==+1) s[i]='p'; else
      if (data[i]==-1) s[i]='m'; else
      wblog(FL,"ERR %s() invalid qdir (%d: %d)",FCT,i,data[i]);
   }; s[len]=0;

   return s;
};


template <class TQ>
qset<TQ>::qset(const qset<TQ> &q1, char op, const qset<TQ> &q2) {

   if (q1.len!=q2.len) wblog(FL,"ERR size mismatch (%d/%d)",q1.len,q2.len);
   this->data=NULL;

   switch (op) {

      case '+':
         INIT(q1.len,q1.data);
         for (unsigned i=0; i<q1.len; i++) q1.data[i]+=q2.data[i];
         break;

      case '-':
         INIT(q1.len,q1.data);
         for (unsigned i=0; i<q1.len; i++) q1.data[i]-=q2.data[i];
         break;

      case '.':
         INIT(q1.len,q1.data,q2.len,q2.data);
         break;

      default: wblog(FL,"ERR invalid operator=%c<%d> ???",op,op);
   }
};


template <class TQ>
qset<TQ>& qset<TQ>::init(
   unsigned l1, const TQ *d1,
   unsigned l2, const TQ *d2
){
   unsigned i=0;
   wbvector<TQ>::RENEW(l1+l2);
   TQ *d=(this->data);

   if (l1!=l2) wblog(FL,
      "WRN cat() qlabel length mismatch (%d+%d)",l1,l2);
   if (l1 && (!d1 || !d2)) wblog(FL,
      "ERR cat(%d+%d) got NULL data (0x%Xl, 0x%Xl)",l1,l2,d1,d2);

   for (   ; i<l1; ++i) { d[i]=d1[i]; }; d+=l1;
   for (i=0; i<l2; ++i) { d[i]=d2[i]; }

   return *this;
};


template <class TQ>
qset<TQ>& qset<TQ>::init(
   unsigned l1, const TQ* d1,
   unsigned l2, const TQ *d2,
   unsigned l3, const TQ *d3
){
   unsigned i=0;
   wbvector<TQ>::RENEW(l1+l2+l3);
   TQ *d=(this->data);

   if (l1!=l2 || l1!=l3) wblog(FL,
      "WRN cat() qlabel length mismatch (%d+%d+%d)",l1,l2,l3);
   if (l1 && (!d1 || !d2 || !d3)) wblog(FL,"ERR cat(%d+%d+%d) "
      "got NULL data (0x%Xl, 0x%Xl, 0x%Xl)",l1,l2,l3,d1,d2,d3);

   for (   ; i<l1; ++i) { d[i]=d1[i]; }; d+=l1;
   for (i=0; i<l2; ++i) { d[i]=d2[i]; }; d+=l2;
   for (i=0; i<l3; ++i) { d[i]=d3[i]; }

   return *this;
};


template <class TQ>
int qset<TQ>::checkNQs(const QVec &qq, const char *F, int L) {

   unsigned i,d0, l=0;

   for (i=0; i<qq.len && l<this->len; i++) {
      d0=qq[i].qlen(); l+=d0;
   }

   if (i!=qq.len || l!=this->len) {
      if (F) wblog(F,L,"ERR qset/QType inconsistency (%d/%d, %d/%d)",
      i,qq.len,l,this->len); return -1;
   }

   return i;
};


template <class TQ>
bool qset<TQ>::isConsistent(
   const QVec &qq, int rank, const char *F, int L
) const {

   unsigned i,d0, l=0, m=(rank>0 ? rank : 1)*qq.len;
   if (qq.isEmpty()) wblog(F_L,"ERR got empty QVec"); 

   for (i=0; l<(this->len) && i<m; i++) { 
      d0=qq[i<qq.len ? i : i%qq.len].qlen();
      l+=d0;
   }

   if (i>qq.len && (rank<=0 || i%qq.len)) {
      if (F) wblog(F,L,
      "ERR incompatible QVec (%d/%d; %d)",i,qq.len,rank);
      return 0;
   }

   if (l!=this->len || i!=(rank>0 ? rank : 1)*qq.len) {
      if (F) wblog(F,L,"ERR qset out of bounds "
         "(%d/%d, %d/%d; %d)",l,this->len,i,qq.len,rank);
      return 0;
   }

   return 1;
};


template <class TQ>
unsigned qset<TQ>::QDim(const QVec &qq, int rank) const {

   unsigned i,j,d0, l=0, n=1;
   if (qq.isEmpty()) wblog(FL,"ERR got empty QVec"); 

   for (i=0; l<this->len && i<this->len; i++) { 
      if (rank<=0) { if (i>=qq.len) wblog(FL,
         "ERR QVec index out of bounds (%d/%d; %d)",i,qq.len,rank);
         d0=qq[i].qlen();
      }
      else d0=qq[i%qq.len].qlen();

      n*=qq[i].QDim(this->data+l);
      l+=d0;
   }

   if (i>qq.len) {
      if (rank<=0 || i%qq.len)
      wblog(FL,"ERR incompatible QVec (%d/%d; %d)",i,qq.len,rank);
   }

   if (l!=this->len || i>l || i%qq.len) wblog(FL,
   "ERR qset/QType mismatch (%d/%d, %d/%d)",i,qq.len,l,this->len);

   return n;
};



template <class TQ>
QSet<TQ>& QSet<TQ>::init(const char *F, int L,
   const CRef<TQ> &A_, const ctrIdx &ica,
   const CRef<TQ> &B_, const ctrIdx &icb,
   wbperm *Pcgd
){

   if (!A_.cgr || !B_.cgr) {
      if (!A_.isAbelian() || !B_.isAbelian()) wblog(FL,
         "ERR %s() invalid CRef data\n%s\n\%s",
         FCT,A_.toStr().data,B_.toStr().data);
      return init();
   }

   const CDATA_TQ &A=(*A_.cgr), &B=(*B_.cgr);

   unsigned na,nb, i,i_,j, l=0, nq=A.t.qlen(),
      ra=A.rank(F_L), rb=B.rank(F_L),
      nab=ra+rb;
   const wbperm &pa=A_.cgp, &pb=B_.cgp;
   char ma[nab], *mb=ma+ra; memset(ma,0,nab*sizeof(char));

   char sa=((A_.conj!=0) ^ (ica.conj!=0) ? -1 : +1),
        sb=((B_.conj!=0) ^ (icb.conj!=0) ? -1 : +1), xd=-sa*sb;

   if (A.t!=B.t || !nq) wblog(F_L,
      "ERR %s() got symmetry mismatch (%s, %s)",
      FCT, A.qStr().data, B.qStr().data
   );

   if (pa.len && !pa.isValidPerm(ra)) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pa.len,ra);
   if (pb.len && !pb.isValidPerm(rb)) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pb.len,rb);

   if (ica.len!=icb.len) wblog(FL,"ERR %s() "
      "invalid set of contraction indices (%d/%d)",FCT,ica.len,icb.len);
   if (A.qdir.len!=ra) wblog(FL,
      "ERR %s() invalid qdir set (A: %d/%d)",FCT,A.qdir.len,ra);
   if (B.qdir.len!=rb) wblog(FL,
      "ERR %s() invalid qdir set (B: %d/%d)",FCT,B.qdir.len,rb);

   for (i=0; i<ica.len; ++i) { j=pa.el1(ica[i]);
      if (j>=ra) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,ica.len);
      if ((++ma[j])>1) wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,ica.len
      );
   }
   for (i=0; i<icb.len; ++i) { j=pb.el1(icb[i]);
      if (j>=rb) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,icb.len);
      if ((++mb[j])>1) wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,icb.len
      );
      if (A.qdir[pa.el1(ica[i])] != xd*B.qdir[j]) wblog(FL,
         "ERR invalid CRef contraction [conj=-(%d,%d)(%d,%d)=%d]\n"
         "having %s @ [%s](%s) <> %s @ [%s](%s)\n"
         "hint: only pairs of {in/out} indices accepted",
         A_.conj, ica.conj, B_.conj, icb.conj, xd,
         A.qdir.toStr().data, pa.toStr().data, ica.toStr().data,
         B.qdir.toStr().data, pb.toStr().data, icb.toStr().data
      );
   }

   na=ra-ica.len;
   nb=rb-icb.len; nab=na+nb;

   t=A.t; qdir.init(nab); qs.init(nq*nab); 
   if (Pcgd) Pcgd->init(nab);

   if (!nab) return *this;
   if (nab>127) wblog(FL,
      "ERR %s() char index out of bounds (%d)",FCT,nab);

   const char *a=A.qdir.data, *b=B.qdir.data;
   const TQ *qa=A.qs.data, *qb=B.qs.data;
   char *c=qdir.data; TQ *qc=qs.data;

   if (Pcgd) {
      for (i=0; i<ra; ++i) if (!ma[i]) { ma[i]=-char(l); ++l; }
      for (i=0; i<rb; ++i) if (!mb[i]) { mb[i]=-char(l); ++l; }

      if (l!=nab) wblog(FL,
         "ERR %s() %d != %d + %d !??",FCT,l,na,nb);
      l=0;

      if (na) { for (i=0; i<ra; ++i) { 
         j=pa.el1(i); if (ma[j]<=0) { Pcgd->data[l++]=-ma[j]; }
      }}
      if (nb) { for (i=0; i<rb; ++i) {
         j=pb.el1(i); if (mb[j]<=0) { Pcgd->data[l++]=-mb[j]; }
      }}
      l=0;
   }

   if (nq==1) {
      if (ra) { for (i_=0; i_<ra; ++i_) { i=pa.el1(i_); if (ma[i]<=0) {
         c[l]=sa*a[i]; qc[l]=qa[i]; ++l;
      }}}
      if (rb) { for (i_=0; i_<rb; ++i_) { i=pb.el1(i_); if (mb[i]<=0) {
         c[l]=sb*b[i]; qc[l]=qb[i]; ++l;
      }}}
   }
   else {
      if (ra) { for (i_=0; i_<ra; ++i_) { i=pa.el1(i_); if (ma[i]<=0) {
         c[l]=sa*a[i];
         for (j=0; j<nq; ++j) { qc[j]=qa[i*nq+j]; }
         qc+=nq; ++l;
      }}}

      if (rb) { for (i_=0; i_<rb; ++i_) { i=pb.el1(i_); if (mb[i]<=0) {
         c[l]=sb*b[i];
         for (j=0; j<nq; ++j) { qc[j]=qb[i*nq+j]; }
         qc+=nq; ++l;
      }}}
   }

   if (l!=nab) wblog(FL,"ERR %s() %d/%d !??",FCT,l,nab);

   return *this;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType &t_, const qset<gTQ>&qs_,
   unsigned iout,
   char isref
){
   if (isref)
        { qs.init2ref(qs_); }
   else { qs=qs_; }

   unsigned n=t_.qlen(), r=qs.len/n;
   if (qs.len%n) wblog(FL,
      "ERR %s() invalid qset (%d @ %d)",FCT,qs.len,n);

   qdir.init_iout(FL, iout? r : 0, iout);

   t=t_; return *this;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType t_, const TQ* qs_, unsigned r, unsigned N,
   unsigned iout
){
   unsigned n=t_.qlen();
   t=t_; qs.init(r*n);

   if (qs_ && qs.len) Wb::cpyStride(qs.data,qs_,n,NULL,r,-1,N);

   qdir.init_iout(FL, iout? r : 0, iout);

   return *this;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType t_, const TQ* qs_, unsigned r, unsigned N,
   const IndexLabels &itags
){
   unsigned n=t_.qlen();

   t=t_; qs.init(r*n);
   if (qs_ && qs.len) Wb::cpyStride(qs.data,qs_,n,NULL,r,-1,N);

   if (itags.len!=r) {
      if (!itags.len) wblog(FL,
         "ERR %s() got empty itags (%d/%d)",FCT,itags.len,r); 
      else wblog(FL,"ERR %s() invalid itags "
         "[%d/%d: %s]",FCT,itags.len,r,itags.toStr().data
      );
   }
   qdir.init(itags);

   return *this;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::init_str(const char *F, int L, const char *s_) {

   unsigned i=0, nd=0, nc=0, r=0, nw=0, n=strlen(s_), e=0;
   wbstring S(s_);
   char *s=S.data, bflag=0;

   strncpy(s,s_,n);

   while (isspace(s[i])) { ++i; }
   if (i) { s+=i; n-=i; }

   for (i=0; i<n; ++i) { if (isspace(s[i])) { s[i]=0; break; }}
   if (i>=n || i<2 || i>6) wblog(FL,
      "ERR %s() invalid QSet specification (QType)\n%s",FCT,s);

   t.init(s); ++i; s+=i; n-=i;

   for (i=0; isspace(s[i]); ++i) {};

   if (s[i]=='(' || s[i]=='{' || s[i]=='[') {
      bflag=s[i]; for (++i; isspace(s[i]); ++i) {};
   }
   if (i) { s+=i; n-=i; }

   for (i=0; i<n; ++i) {
      if (isspace(s[i])) ++nw;
      else {
         if (isdigit(s[i])) ++nd; else
         if (s[i]==',' || s[i]==';') ++r; else
         if (s[i]=='*') ++nc; else
         if (bflag) {
            if (s[i]==')' && bflag=='(' ||
                s[i]=='}' && bflag=='{' ||
                s[i]==']' && bflag=='['
            ) { bflag=0; continue; } else e=1;
         }
         else e=2;
         if (e) wblog(FL,
            "ERR %s() invalid qset string\n%s (%d/%d)",FCT,s,i+1,n
         );
      }
   }

   ++r;

   if (nw==0 && nd%r==0) {
      unsigned l=0;
      qdir.init2val(r,+1); qs.init(nd);
      for (r=0, i=0; i<n; ++i) {
         if (isdigit(s[i])) qs[l++]=s[i]-'0'; else
         if (s[i]==',' || s[i]==';') ++r; else
         if (s[i]=='*') qdir[r]=-1;
      }
   }
   else wblog(FL,
      "ERR %s() string not yet interpreted (%d,%d,%d,%d)\n%s",
      FCT,nd,r,nc,nw, s_
   );
};


template <class TQ>
bool QSet<TQ>::isScalar() const {

   unsigned i=0, r=rank(FL), r2=r/2, m=t.qlen(), n2=r2*m;

   if (r<2) wblog(FL,"ERR %s() got r=%d QSet !?",FCT,r); 
   if (qdir.len!=r || qs.len!=r*m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,toStr().data);

   for (; i<r2; ++i) { if (qdir[i]<=0 || qdir[i+r2]>=0) return 0; }

   for (i=0; i<n2; ++i) { if (qs[i]!=qs[i+n2]) return 0; }

   if (t.isNonAbelian() && r%2) {
   for (i=2*n2; i<qs.len; ++i) { if (qs[i]!=0) return 0; }}

   return 1;
};


template <class TQ>
bool QSet<TQ>::isJ1Symbol() const {

   unsigned i=0, r=rank(FL), m=t.qlen();

   if (r<2) wblog(FL,"ERR %s() got r=%d QSet !?",FCT,r); 
   if (qdir.len!=r || qs.len!=r*m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,toStr().data);
   if (r!=2) return 0;

   for (; i<r; ++i) { if (qdir[i]<=0) return 0; }; i=0;

   for (unsigned l=2*m-1; i<m; ++i) { if (qs[i]!=qs[l-i]) return 0; }

   return 1;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::Sort(wbperm *cgp, char *conj, char iflag) {

   if (isEmpty() || isSorted()) {
      if (cgp) { cgp->init(); }; if (conj) { (*conj)=0; }
      return *this;
   }

   unsigned i=0, l=0, m=t.qlen(), r=(m ? qs.len/m : 0);
   wbMatrix<TQ> dQ(r,1+m);
   char cflag=0, *d=qdir.data;
   TQ *dq=dQ.data;
   wbperm P;

   if (!m || qs.len%m || qdir.len!=r || r>127) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,toStr().data);
   if (r<2) wblog(FL,"ERR %s() got rank-%d QSet !?",FCT,r);
   
   for (; i<r; ++i, dq+=dQ.dim2) {
      if (d[i]>0) ++l; else
      if (!d[i]) wblog(FL,"ERR %s() invalid qdir\n%s",FCT,toStr().data);
      dq[0]=(d[i]>0 ? 0 : 1);
      MEM_CPY<TQ>(dq+1, m, qs.data+i*m);
   }

   if (l<(r+1)/2) {
      cflag=1; dq=dQ.data;
      for (i=0; i<r; ++i, dq+=dQ.dim2) { dq[0]=(!dq[0]); }
   }

   dQ.sortRecs(P);

   if (2*l==r) { int q=0;
      dq=dQ.data+1;
      for (i=0; i<l; ++i, dq+=dQ.dim2) {
         if ((q=Wb::cmpRange(dq,dq+l*dQ.dim2,m))) {
            if (q>0) { P.Rotate(l); cflag=(2+!cflag); }
            else { i=r; }; break;
         }
      }
   }

   if (cflag%2) {
      if (conj) { (*conj)=1; }
      for (i=0; i<r; ++i) d[i]=-d[i];
   }
   else {
      if (conj) { (*conj)=0; }
   }

   qdir.Permute(P);

   if (P.isIdentityPerm()) {
      if (cgp) { cgp->init(); }
      return *this;
   }
   else if (cgp) {
      if (iflag)
           P.getIPerm(*cgp);
      else P.save2(*cgp);;
   }

   dq=dQ.data+1; if (cflag<2) l=0;
   for (i=0; i<r; ++i) {
      MEM_CPY<TQ>(qs.data+i*m, m, dq+((i+l)%r)*dQ.dim2);
   }

   return *this;
};


template <class TQ>
bool QSet<TQ>::isSorted() const {

   if (isEmpty() || !qdir.isSorted()) return 0;

   unsigned i=0, l=0, m=t.qlen(), r=(m ? qs.len/m : 0);
   if (!m || qs.len%m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,toStr().data);
   if (r<2) wblog(FL,
      "ERR %s() got rank-%d QSet !?",FCT,r);

   for (; i<r; ++i) { if (qdir[i]<=0) break; }
   if (i>1) {
      wbMatrix<TQ> X(i,m,qs.data,'r');
      if (!X.recsSorted(+1)) return 0;
   }

   for (l=i; i<r; ++i) {
      if (!qdir[i]) wblog(FL,
         "ERR %s() got invalid qdir=%s",FCT,qdir.toStr().data);
      if (qdir[i]>0) return 0;
   }

   if (l+1<r) {
      wbMatrix<TQ> X(r-l,m,qs.data+l*m,'r');
      if (!X.recsSorted(+1)) return 0;
   }

   if (l<(r+1)/2) return 0;
   if (2*l==r) { m*=l;
      if (Wb::cmpRange(qs.data, qs.data+m,m)>0) return 0;
   }

   return 1;
};


template <class TQ>
int QSet<TQ>::checkQ_U1(const char *F, int L) const {

   if (t.type!=QTYPE_U1) {
      if (F) wblog(F,L,"ERR %s() got %s",FCT,t.toStr().data); 
      return -1;
   }
   if (!qdir.len) { if (F) wblog(FL,
      "WRN %s() got empty QSet (%s)",FCT,toStr().data);
      return -2;
   }
   if (qdir.len!=qs.len) { if (F) wblog(FL,
      "ERR %s() qdir inconsistency (%s)",FCT,toStr().data);
      return -3;
   }

   unsigned i=0, nin=0, nout=0;
   TQ qin=0, qout=0;

   for (; i<qdir.len; ++i) {
      if (qdir[i]>0) { qin +=qs[i]; ++nin;  } else
      if (qdir[i]<0) { qout+=qs[i]; ++nout; } else wblog(FL,
     "ERR %s() got undetermined qdir=%s",FCT,qdir.toStr().data);
   }
   if (nin && nout) {
      if (qin!=qout) { if (F) wblog(F,L,
         "ERR %s() broken charge conservation (%g != %g !?!)",
         FCT,double(qin),double(qout));
         return 2;
      }
   }
   return 0;
};


template <class TQ>
int QSet<TQ>::checkQ_SU2(const char *F, int L) const {

   if (t.type!=QTYPE_SUN || t.sub!=1) { if (F) wblog(F,L,
      "ERR %s() got invalid type %s",FCT,t.toStr().data); 
      return -1;
   }
   if (qdir.len<2) { if (F) wblog(FL,
      "WRN %s() got rank-%d QSet (%s) !??",FCT,qdir.len,toStr().data);
      return -2;
   }
   if (qdir.len!=qs.len) { if (F) wblog(FL,
      "ERR %s() qdir inconsistency (%s)",FCT,toStr().data);
      return -3;
   }

   unsigned i=0, nin=0;

   for (; i<qdir.len; ++i) {
      if (qdir[i]>0) { ++nin; } else
      if (!qdir[i]) wblog(FL,
     "ERR %s() got undetermined qdir=%s",FCT,qdir.toStr().data);
   }

   if (!nin || nin==qdir.len) {
      wblog(F_L,"WRN %s() got all-in (out) %s",FCT,toStr().data);
      return 2;
   }

   if (qdir.len==2) {
      if (qs[0]!=qs[1]) { if (F) wblog(FL,
         "ERR %s() got scalar CGC data %s",FCT,toStr().data);
         return 3;
      }
   }
   else if (qdir.len==3) {
      unsigned i1,i2,i3=0;
      if (nin==1)
         for (i=0; i<qdir.len; ++i) { if (qdir[i]>0) { i3=i; break; }}
      else
         for (i=0; i<qdir.len; ++i) { if (qdir[i]<0) { i3=i; break; }}

      if (i3==0) { i1=1; i2=2; } else
      if (i3==1) { i1=0; i2=2; }
      else       { i1=0; i2=1; }

      if (qs[i1]<qs[i2]) { i=i1; i1=i2; i2=i; }

      if (qs[i3] > qs[i1]+qs[i2] || qs[i3] < qs[i1]-qs[i2]) {
         if (F) wblog(F,L,"ERR %s() "
            "invalid SU2 addition rule (%g + %g = %g !?!)\n%s",FCT,
            double(qs[i1]),double(qs[i2]),double(qs[i3]),toStr().data);
         return 4;
      }
   }

   return 0;
};


template <class TQ>
QSet<TQ>& QSet<TQ>::permute(
   QSet<TQ> &B, const wbperm &P, char iflag) const {

   unsigned r=rank(FL);
   if (!P.isValidPerm(FL,r)) wblog(FL,
      "ERR invalid permute on CData (%d/%d)",P.len,r);

   B.t=t;
   qdir.permute(B.qdir,P,iflag);
   qs.blockPermute(P,B.qs,iflag);

   return B;
};


template <class TQ>
wbstring QSet<TQ>::QStrS(const wbperm *cgp, char s_) const {
   wbstring s;
   unsigned k=0; char sep[2]=" ";

   if (t.qlen()<2 || qs.allIn(0,9)) sep[0]=0;

   if (cgp) {
      QDir qd; qset<TQ> qx;
      qdir.permute(qd,*cgp); k=qd.isSorted();
      qs.blockPermute(*cgp,qx);
      s=qx.toStrf("",sep,t.qlen(),";");
   }
   else {
      k=qdir.isSorted();
      s=qs.toStrf("",sep,t.qlen(),";");
   }

   if (k) { --k;
      unsigned i=0, r=0, n=s.len; char *c=s.data;
      for (; i<n; ++i) if (c[i]==';') {
         if ((++r)!=k) { c[i]=','; } else { c[i]=s_; }
      }
   }
   return s;
};


template <class TQ>
wbstring QSet<TQ>::QStr() const { 
   wbstring sout;

   unsigned i=0, j=0, l=0, sep=1, n, m=t.qlen(); 
   const TQ *q=qs.data;

   if (!t.isAbelian() && qs.wbvector<TQ>::allIn(0,9))
        { n=2*qs.len+qdir.len; sep=0; }
   else { n=4*qs.len+qdir.len; }
   sout.init(n+1);

   if (!n) return sout;

   char *s=sout.data;
   wbstring fmt; fmt.init2Fmt(qs[0]);

   for (; i<qdir.len && l<n; ++i, q+=m) { if (i) { s[l++]=','; }
      for (j=0; j<m && l<n; ++j) {
         if (j && sep) { s[l++]=' '; }
         l+=snprintf(s+l,n-l,fmt.data,q[j]);
      }
      if (qdir[i]<0 && l<n) { s[l++]='*'; }
   }

   if (l>=n) { s[n]=0; wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,s); }
   sout.len=l+1;

   return sout;
};


template <class TQ>
wbstring QSet<TQ>::toStr(const char *istr) const {
   wbstring sout;

   unsigned m=t.qlen(), n=127;

   if (!t.isAbelian() && qs.wbvector<TQ>::allIn(0,9)) { n=63; }
   if (istr && !istr[0]) istr=0;

   if (qs.len!=m*qdir.len) wblog(FL,"ERR %s() "
      "severe QSet inconsistency (%d!=%d*%d) !?",FCT,qs.len,m,qdir.len);


   sout.init(n);

   char *s=sout.data;
   unsigned l=snprintf(s,n,"%5s (",t.toStr().data);

   if (l<n) {
      l+=snprintf(s+l,n-l,"%s",QStr().data);
   }
   if (l<n) {
      if (istr)
           { l+=snprintf(s+l,n-l,") %s",istr); }
      else { s[l]=')'; s[++l]=0; }
   }

   if (l>=n) { s[n]=0; wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,s); }
   sout.len=l+1;

   return sout;
};


template <class TQ>
wbstring QSet<TQ>::toTag() const { wbstring s(64);

   unsigned l=0, n=s.len-1;

   char sep[2]=",";
   if (t.qlen()<=1 || qs.wbvector<TQ>::allIn(0,9)) sep[0]=0;

   if (!t.validType() || !qdir.len || !qs.len) wblog(FL,
      "ERR %s() got empty CData\n%s",FCT,toStr().data);

   l=snprintf(s.data,n,"%s[%s: %s]",
     t.toStr('t').data, qdir.toStr().data,
     qs.wbvector<TQ>::toStrf("",sep,t.qlen(),",").data
   );

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};


template <class TQ>
mxArray* QSet<TQ>::toMx() const {

   mxArray *S=mxCreateStruct(1,1);
   add2MxStruct(S,0); return S;

};


template <class TQ>
mxArray* QSet<TQ>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={"type","qset","qdir"};
   return mxCreateStructMatrix(m,n,3,fields);
};


template <class TQ>
void QSet<TQ>::add2MxStruct(mxArray *S, unsigned l) const {

   mxSetFieldByNumber(S,l,0, t.toMx());
   mxSetFieldByNumber(S,l,1, qs.toMx());
   mxSetFieldByNumber(S,l,2, qdir.toMx());
};


template <class TD>
template <class TQ>
cdata<TD>& cdata<TD>::init(const CRef<TQ> &R, char full) {

   if (!R.cgr) wblog(FL,"ERR %s() got NULL cgref !??",FCT);
   if (!R.cgw.len) wblog(FL,"ERR %s() got empty cgw !??",FCT);

   unsigned r=R.cgr->rank(FL), m=R.cgr->gotOM(FL);

   if (m<=1) {
      init(R.cgr->cgd) *= R.cgw[0];
      if (m) {
         if (r && this->SIZE[r]==1) {
            this->skipTrailingSingletons(FL,r-1);
            if (this->SIZE.len!=r-1) wblog(FL,
               "ERR %s() failed to remove OM singleton dimension",FCT
            );
         }
         else wblog(FL,"WRN %s() got "
            "invalid OM=%d CData (%d/%d)",FCT,m, R.cgr->cgd.SIZE.len, r
         );
      }
   }
   else if (R.cgw.len==m) {
      R.cgr->cgd.cgsparray::contract(FL,r,R.cgw, (cgsparray&)(*this));
   }
   else {
      unsigned i=0, l=R.cgw.len;
      wbvector<RTD> x(m); RTD z(0);

      for (; i<l; ++i) x[i]=R.cgw[i];
      for (; i<m; ++i) x[i]=z;

      R.cgr->cgd.cgsparray::contract(FL,r,x, (cgsparray&)(*this));
   }

   if (full) {
      if (R.cgp.len) Permute(R.cgp);
      if (R.conj) Conj();
   }

   return *this;
};


template <class TD>
cdata<TD>& cdata<TD>::init(
   const char *F, int L, const mxArray *C, unsigned k,
   const QType &q, unsigned r
){
   if (C==NULL) {
       if (k) wblog(FL,"ERR %s() got k=%d for NULL mxArray !??",FCT,k);
       cgsparray::init();
       return *this;
   }

   try { cgsparray::init(F_L,C,k); }
   catch (...) { wblog(F_L,
      "ERR %s() invalid (sparse) CG space info.cgs{} (%s)",
      FCT,mxGetClassName(C));
   }

   const bool isc=q.isAbelian();

   if (isc && !isScalar()) wblog(F_L,
      "ERR %s() CG scalar inconsistency (%d/%d %s)",
       FCT, isScalar(), isc, this->sizeStr().data
   );

   if (int(r)>=0) {
   this->addTrailingSingletons(FL,r); }

   return *this;
};


template <class TD>
cdata<TD>& cdata<TD>::initIdentity(
   const QType &t, SPIDX_T d, TD dval
){ 
   if (d!=1 || t.isNonAbelian())
        wbsparray<TD>::initIdentity(d,dval);
   else wbsparray<TD>::initScalar(dval);

   return *this;
};


template <class TD>
cdata<TD>& cdata<TD>::initIdentity(const wbvector<SPIDX_T> &S) {

   if (S.len==1 && S[0]==1) {
      wbsparray<TD>::initScalar(1.); }
   else if (S.len==2 && S[0]==S[1]) {
      wbsparray<TD>::initIdentity(S[0]);
   }
   else {
      wblog(FL,"ERR wbsparse::%s() non-square (%s) !??",
      FCT,S.toStrD().data);
   }
   return *this;
};


template <class TD>
unsigned cdata<TD>::getOM(
   const char *F, int L, const char *istr) const {

   double x2=this->norm2(), e=std::fabs(x2-std::round(x2));
   unsigned m=(x2+0.5), r=this->SIZE.len;

   if (e>1E-12) wblog(F_L,"ERR %s() "
      "got non-integer |cgc|^2 (%s%g @ %.3g)",FCT,istr?istr:"",x2,e);
   if (!r) {
      if (this->isDiag()) {
         if (m!=1) wblog(FL, "ERR %s() got |cgc|^2=%g\n"
            "for diagonal CData (D.len=%d) !?",FCT,x2,this->D.len); }
      else if (m) wblog(F_L,
         "ERR %s() got om=%g for empty CData !?",FCT,x2);
      return 0;
   }

   if (m<=1) return 0;

   const SPIDX_T *sd=this->SIZE.data;
   
   for (unsigned s=1, i=r-1; i<r; --i) {
      if ((s*=sd[i])==m) {
         if (r-i>1) wblog(FL,
            "WRN %s() got %d OM dimensions !?",FCT,r-i);
         return (r-i);
      }
      if (s>m) wblog(F_L,"ERR %s() OM cgn2=%g "
         "incompatbile with SIZE=[%s]",FCT,x2,this->sizeStr().data
      );
   }

   return 0;
};


template <class TD>
TD cdata<TD>::NormSignC(
   const char *F, int L,
   unsigned m,
   TD eps1, TD eps2
){
   TD nrm;

   if (this->D.len<=1) { if (!this->D.len) return 0; else {
      TD a=ABS(this->D[0]);
      if (m) wblog(FL,"ERR %s() got m=%d for scalars !?",FCT,m);
      if (a<eps1) {
         if (a>eps2 && F) wblog(F_L,
            "WRN %s() CGC noise %.3g (1: %g, %g)",
            FCT, double(this->D[0]), double(eps1), double(eps2));
         if (double(a)<CG_SKIP_EPS1) { this->D[0]=0; return 0; }
      }
      nrm=this->D[0]; this->D[0]=1;
      return nrm;
   }}

   this->Compress(F_L,CG_SKIP_EPS1);

   SPIDX_T i, n=this->D.len;
   unsigned r=this->rank(), l=this->SIZE.len;
   TD *d=this->D.data, x1=1E99, x2=0;

   if (m>2 || m>r || r<2 || (r!=l && l && r!=2)) wblog(FL,
      "ERR got OM=%d (r=%d/%d)",m,r,l);


   nrm=this->D.norm2();

   if (m) {
      nrm/=maxRange(this->SIZE.data+(r-m),m);
   }
   else if (r==2) {
   }
   else if (r>4) {
      wblog(FL,"WRN %s() got rank-%d cgd (OM=%d)",FCT,r,m);
   }

   nrm=SQRT(nrm);

   if (nrm<eps1) {
      if (nrm>eps2) wblog(FL,
         "WRN %s() CGC noise %.3g (2: %g; %g)",
         FCT, double(nrm), double(eps1), double(eps2));
      return 0;
   }

   for (i=0; i<n; ++i) if (ABS(d[i])>=eps1) {
      if (d[i]<0) { nrm=-nrm; }; break;
   }
   if (i==n) wblog(FL,
      "ERR %s() failed to determine sign (%g)",FCT,double(eps1));

   (this->D)/=nrm;

   CG::FixRational(FL, this->D.data, this->D.len, 4);
   CG::FixRational(FL, &nrm, 1, 4);



   return nrm;
};


template <class TD>
cdata<TD>& cdata<TD>::permute(
   cdata<TD> &B, const wbperm &P0, char iflag
) const {
   if (isScalar()) { B.init(*this); } else
   if (!this->SIZE.len && this->D.len>1) {
      if (!this->isDiag() || P0.len!=2) wblog(FL,"ERR %s() "
         "invalid diagonal (%s, D=%d)",FCT,P0.toStr().data,this->D.len);
      B.init(*this);
   }
   else { 
      unsigned r=this->SIZE.len;
      if (P0.len!=r && P0.len+1!=r && r!=3) wblog(FL,
         "ERR %s() invalid permutation [%s]",FCT,P0.toStr().data); 

      wbperm P(P0,this->SIZE.len);
      wbsparray<TD>::permute(P,B,iflag);
   }

   return B;
};



template <class TD>
TD cdata<TD>::contract(const char *F, int L,
   const wbindex &ica, const cdata<TD> &B,
   const wbindex &icb, cdata<TD> &Cin,
   cdata<TD> *Cx, wbperm *P
 ) const {

   TD nrm=1;
   cdata<TD> C;

   if (this->D.len<=1 && B.D.len<=1) {
      if (this->contract_scalar(F_L,wbindex(ica),B,wbindex(icb),C)) {
         if (Cx) wblog(FL,"ERR %s() got Cx for scalar C !?",FCT);
         if (!Cin.isEmpty()) wblog(FL,
            "ERR %s() got non-empty Cin (%s)",FCT,this->sizeStr().data);
         nrm=C.NormSignC(F,L);
         C.save2(Cin); return nrm;
      }
      else wblog(FL,"WRN %s() "
         "got non-scalar contraction for D.len=%d/%d !?",
         FCT,this->D.len, B.D.len
      );
   }

   unsigned ma=this->getOM(FL,"A: "), mb=B.getOM(FL,"B: ");

   this->cgsparray::contract(F,L,ica, B, icb, C, P ? *P : wbperm());

   nrm=C.NormSignC(F,L,ma+mb);


   if (Cx) Cx->init();

   if (Cin.isEmpty()) { C.save2(Cin); }
   else {
      if (C.SIZE!=Cin.SIZE) {
         if (!C  .SIZE.len) C  .diag2reg(); else
         if (!Cin.SIZE.len) Cin.diag2reg();
         if (C.SIZE!=Cin.SIZE) {
            MXPut(FL,"ix").add(C,"c").add(Cin,"c0").add(nrm,"nrm");
            wblog(FL,"ERR size inconsistency in C.cgs[] (%s; %s)",
            C.sizeStr().data, Cin.sizeStr().data);
         }
      }

      if (fabs(double(nrm))<1E-14) return nrm;

      double d=MAX(SPIDX_T(1),C.numel()), e,e2=std::sqrt(C.normDiff2(Cin)/d);
      if (e2>CG_EPS1) {
         if ((e=Cin.norm2())<CG_EPS2) C.save2(Cin);
         else if (Cx) C.save2(*Cx);
         else {
            MXPut(FL,"i_").add(C,"c").add(Cin,"c0").add(nrm,"nrm");
            wblog(FL,"ERR CGC difference: e=%g\n%s (%s) * %s (%s) => %s",
            e2, this->sizeStr().data, (ica+1).toStr().data,
            B.sizeStr().data, (icb+1).toStr().data, C.sizeStr().data);
         }
      }
      else if (e2>CG_EPS2) { wblog(FL,"WRN CGC difference: e=%g",e2); }
   }

   return nrm;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::init(const CRef<TQ> &R) {

   cgd.init(R);

   cstat.t=CGD_FROM_CGR;
   if (R.cgp.len)
        R.cgr->QSet<TQ>::permute((QSet<TQ>&)(*this),R.cgp);
   else (*this)=(QSet<TQ>&)(*R.cgr);

   if (R.conj) this->QSet<TQ>::Conj();

   return *this;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::RefInit(
   const char *F, int L, const CData<TQ,TD> &B, char bare) {

   (*this)=(QSet<TQ>&)B;
   cstat=B.cstat; cstat.t=CGD_REF_INIT;

   if (bare || !isSymmmetric()) {
      wbvector<SPIDX_T> S; B.getSize(S,bare);
      cgd.wbsparray<RTD>::init(S,0);
   }
   else {
      wbvector<SPIDX_T> S;
      wbvector<RTD> cgt; B.getSize(S); B.trace(FL,cgt);
      RefInit_auxtr(FL,S,cgt);
   }

   return *this;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::Reduce2Ref(const char *F, int L) {

   if (cstat.t==CGD_REF_INIT) { return *this; }

   wbvector<SPIDX_T> S; getSize(S);

   if (CG_VERBOSE>5) wblog(F_L," *  %s() %s",FCT,toStr().data); 

   if (!isSymmmetric()) {
      cgd.wbsparray<RTD>::init(S,0);
      cstat.t=CGD_REF_INIT;
   }
   else {
      wbvector<RTD> cgt; trace(FL,cgt);
      RefInit_auxtr(FL,S,cgt);
   }

   return *this;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::RefInit_auxtr(
   const char *F, int L,
   const wbvector<SPIDX_T> &S, const wbvector<RTD> &cgt
){
   unsigned r=this->qdir.len, m=0;

   if (S.len<r || S.len>r+1) wblog(FL,
      "ERR %s() unexpected size [%s] (r=%d)",FCT,S.toStr().data,r);

   if (!cgt.len) {
      if (!S.len) { cgd.init(); return *this; } else 
      if (isSymmmetric()) wblog(FL,
         "ERR %s() got empty cgt (S=[%s]) !?",FCT,S.toStr().data);
   }
   else {
      m=(S.len>r ? S[r] : 1);
      if (cgt.len!=m) wblog(FL,"ERR %s() "
         "got OM mismatch (cgt.len=%d/%d)",FCT,cgt.len,m);
      if (m>99) wblog(FL,"WRN %s() got unexpected OM=%d !?",FCT,m);
   }

   cgd.wbsparray<RTD>::init(S,m);
   if (m) {
      cgd.D=cgt;
      if (m>1) {
         unsigned i=1; SPIDX_T *idx=cgd.IDX.data+2*S.len-1;
         for (; i<m; ++i, idx+=S.len) { idx[0]=i; }
      }
   }

   cstat.t=CGD_REF_INIT;

   return *this;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::LoadRef(const char *F, int L)  {

   if (cstat!=CGD_REF_INIT) {
      return *this;
   }

   CDATA_TQ &C=gCG.getBUF(F_L, (QSet<TQ>&)*this, 'L');
   return *this;
};


template <class TQ, class TD>
unsigned CData<TQ,TD>::rank(const char *F, int L) const {

   unsigned n=this->t.qlen(), r=n? this->qs.len/n : -1;

   if (!n || this->qs.len%n) wblog(F_L,
      "ERR %s() invalid CData !?? (%d @ %d)",FCT,this->qs.len,n); 
   if (this->qdir.len!=r) wblog(F_L,
      "ERR %s() inconsistent qdir %s",FCT,this->toStr().data);
   if (!r && isEmpty()) {
      wblog(F_L,"WRN %s() got empty CData",FCT);
      return 0;
   }

   if (cgd.SIZE.len==r+1) {
      if (r<3 || !this->t.hasMP(r)) {
         if (cgd.SIZE[r]!=1) { MXPut(FL,"y").add(cgd,"c");
            wblog(F_L,"ERR %s() got OM for rank-%d CGS\n%s: %s",
            FCT, r, toStr().data, cgd.sizeStr().data);
         }
      }
   }
   else if (cgd.SIZE.len!=r) {
      if (r==2 && cgd.isDiag()) {
         return r;
      }
      if (!isAbelian()) wblog(F_L,
         "ERR %s() got CData rank mismatch (%s: %d/%d) !?",
         FCT, qStr().data, cgd.SIZE.len, r
      );
   }
   return r;
};


template <class TQ, class TD>
template <class T>
wbvector<T>& CData<TQ,TD>::getSize(
   wbvector<T> &S, char bare, const wbperm *cgp) const {

   unsigned i=0, l=(cgp ? cgp->len : 0);
   const SPIDX_T *s=cgd.SIZE.data;

   S.init(cgd.SIZE.len);

   if (!S.len) {
      cgd.wbsparray<TD>::getSize(S);
      if (S.len==2 && this->qdir.len==2) { 
         if (l && l!=2) wblog(FL,
            "ERR %s() invalid cgp=[%s] (2)",FCT,cgp->toStr().data);
         return S;
      }
      if (!S.len && isScalar()) {
         if (l && l!=this->qdir.len) wblog(FL,"ERR %s() "
            "invalid cgp=[%s] (%d)",FCT,cgp->toStr().data,this->qdir.len);
         return S.init2val(this->qdir.len,1);
      }
      if (S.len || l) wblog(FL,
         "ERR %s() got %s !?",FCT,toStr().data);
      return S;
   }

   if (l) {
      if (l+1<S.len || l>S.len) wblog(FL,
         "ERR %s() cgp/SIZE mismatch (%d/%d)",FCT,l,S.len);
      if (!cgp->isIdentityPerm()) {
         for (; i<l; ++i) { S.data[i]=s[cgp->data[i]]; }
      }
   }

   for (; i<S.len; ++i) { S.data[i]=s[i]; }

   if (bare) { i=0;
      if (S.len==this->qdir.len+1) {
         if (S.len>3) S.len=this->qdir.len; else i=1;
      }
      else if (S.len!=this->qdir.len) i=2;

      if (i) wblog(FL,
         "ERR %s() got invalid OM setting (%d/%d)\n%s",
         FCT,S.len,this->qdir.len, toStr().data
      );
   }

   return S;
};


template <class TQ, class TD>
int CData<TQ,TD>::checkConsistency(const char *F, int L) const {

   if (isEmpty()) return 1;

   unsigned r,d, i=0, l=0, n=this->t.qlen(); int e=0;

   if (!n || !this->qs.len || this->qs.len%n) { e=11; }
   else { r=this->qs.len/n;
      if (!cgd.SIZE.len) { if (!r) return 1;
         if (r==2 && cgd.isDiag()) {
            for (; i<r; ++i, l+=n) {
               d=gRG.qdim(FL,this->t,this->qs.data+l);
               if (d!=cgd.D.len) { e=31+i; break; }
            }
         }
         else return 30;
      }
      else if (cgd.SIZE.len<r || cgd.SIZE.len>r+1) { e=12; } else
      for (; i<r; ++i, l+=n) {
         d=gRG.qdim(FL,this->t,this->qs.data+l);
         if (d!=cgd.SIZE[i]) { e=21+i; break; }
      }
   }
   if (!e && r && cgd.SIZE.len>r) {
      if (cgd.SIZE[i]>=10)
         wblog(FL,"WRN %s() got OM=%d",FCT,cgd.SIZE[i]);
      if (cgd.SIZE[i]>=100) e=100;
   }

   if (e && F) {
      MXPut(FL).add(*this,"cg"); wblog(F_L,
     "ERR %s() invalid CData (e=%d)\n%s",FCT,e,toStr('v').data);
   }

#ifndef NSAFEGUARDS
   if (this->t.isAbelian()) {
      if (!cgd.isScalar()) { e=999;
          if (F) wblog(FL,"ERR %s() invalid %s CData (%s)",
          FCT, qStr().data, cgd.sizeStr().data);
      }
   }
#endif

   if (!e) {
      if (this->t.isAbelian()) return this->checkQ_abelian(F,L);
      if (this->t.isSU2()    ) return this->checkQ_SU2(F,L);
   }

   return e;
};


template <class TQ, class TD>
int CData<TQ,TD>::checkNormSign(
   const char *F, int L, char xflag) const {

   if (isRefInit() || isAbelian()) return 0;

   unsigned r=rank(), l=cgd.SIZE.len, m=getOM(F_L);

   if (m>1 && !this->t.hasMP(r)) {
       if (F) wblog(F_L,
          "ERR %s() got m=%d for rank-%d/%d CGC in %s !?",
          FCT,m,r,l,this->t.toStr().data);
      return 1;
   }
   if ((!cgd.isDiag() && (r+1<l || r>l)) || r<2) {
      if (F) wblog(F_L,
         "ERR %s() got m=%d for rank-%d/%d CGC !?",FCT,m,r,l);
      return 2;
   }
   if (!m || !cgd.D.len) {
      if (F) wblog(FL,
         "ERR %s() got empty CData !?\n%s",FCT,toStr().data);
      return 3;
   }
   if (r==2 && m>1) {
      if (F) wblog(FL,
         "ERR %s() got invalid rank-2 CGC (m=%d) !?",FCT,m);
      return 4;
   }
   if (CG::signFirstVal(F_L,cgd.D.data, cgd.D.len)<0) {
      if (F) wblog(F_L,
         "ERR %s() got unconventional sign of CData",FCT);
      return 5;
   }

   if (m<2 || !xflag) {
      RTD cn2=cgd.D.norm2(); if (m>1) { cn2/=m; }
      double e=std::fabs(double(cn2-1));

      if (e>1E-12) {
         if (F) wblog(F_L,"ERR %s() "
            "got |cgc|^2 = %g (r=%d; OM=%d)",FCT,double(cn2),r,m);
         return 11;
      }
   }
   else if (xflag!='X') {
      wbvector<RTD> cv2; cgd.norm2vec(r,cv2);
      RTD cn2=cv2.avg(); double e=std::fabs(double(cn2-1));

      if (cv2.len!=m) {
         if (F) wblog(FL,
            "ERR %s() got OM mismatch %d/%d !?",cv2.len,m);
         return 21;
      }
      if (!cv2.allEqual2(RTD(1),CG_SKIP_EPS2)) {
         if (!cv2.allEqual2(cn2,CG_SKIP_EPS2)) {
            if (F) wblog(FL,"ERR %s() got "
               "non-const normalization (%g) !?",FCT,double(cn2));
            return 22;
         }
         else {
            if (F) wblog(FL,"ERR %s() got "
               "non-normalized CData (%g) !?",FCT,double(cn2));
            return 23;
         }
      }
      if (e>1E-12) {
         if (F) wblog(FL,"ERR %s() got "
            "unexpected normalization %g (@ %g)",FCT,double(cn2),e);
         return 24;
      }
   }
   else {
      cdata__ E; wbindex ia(r,'i');
      cgd.cgsparray::contract(F_L,ia,cgd,ia, E);
      if (!E.isIdentityMatrix()) { if (F) {
         MXPut(F_L,"Ix").add(cgd,"cgd").add(E,"E");
         wblog(F_L,"ERR %s() CRef got non-orthonormal CData !?",FCT); }
         return 99;
      }
   }

   return 0;
};


template <class TQ, class TD>
int CData<TQ,TD>::checkQ(
   const char *F, int L, const QSet<TQ> &Q, const char *istr
 ) const {


   if ((*this)!=Q) {
      if (F || istr) wblog(F_L,
         "ERR %s() CRef inconsistency %s %s\n(%s <> %s)",
         FCT, istr? "in":"", istr? istr : "",
         this->toStr().data, Q.toStr().data
      );
      return 1;
   }

   return 0;
};


template <class TQ, class TD>
bool CData<TQ,TD>::checkAdditivityZ(const char *F, int L,
   const wbMatrix<double> &Z1,
   const wbMatrix<double> &Z2, const wbMatrix<double> &Z,
   TD eps
 ) const {

   unsigned nz=Z1.dim2; double e=0;
   wbvector<double> z1,z2,z3,zz;
   const SPIDX_T *I=cgd.IDX.data;

   checkConsistency(FL);
   if (cgd.numel()<=1) {
      if (!Z1.dim1 && !Z2.dim1 && !Z.dim1) return 1;
      if (!Z1.dim1 || !Z2.dim1 || !Z.dim1) wblog(FL,"ERR %s()",FCT);
   }

   if (cgd.SIZE.len<3 || cgd.SIZE.len>4) wblog(FL,
      "ERR CData::%s() inconsistency (cgd: %s)",
       FCT,cgd.sizeStr().data);

   for (SPIDX_T n=cgd.IDX.dim1, m=cgd.IDX.dim2, i=0; i<n; ++i, I+=m) {
       if (fabs(cgd.D[i])<eps) continue;

       z1.init2ref(nz,Z1.ref(I[0])); zz =z1;
       z2.init2ref(nz,Z2.ref(I[1])); zz+=z2;
       z3.init2ref(nz,Z .ref(I[2])); zz-=z3;

       if ((e=zz.norm2())>CG_EPS2) {
          MXPut(FL,"q").add(*this,"A").add(wbvector<SPIDX_T>(m,I)+1,"I")
           .add(z1,"z1").add(z2,"z2").add(z3,"z3").add(zz,"zz");
          wblog(F_L,
            "ERR z-labels not additive (%d @ %.3g => %d,%d,%d @ %.3g)",
             i+1, double(cgd.D[i]), I[0]+1, I[1]+1, I[2]+1, e
          );
       }
   }
   return 1;
};


template <class TQ, class TD>
bool CData<TQ,TD>::sameType(
   const CData<TQ,TD> &S, const char *F, int L
 ) const {

   if (this->t!=S.t) { 
      if (F) wblog(F_L,"ERR %s() got different type (%s; %s)",
          FCT, qStr().data, S.qStr().data);
      return 0;
   }

   if (this->qs.len!=S.qs.len) {
      if (F) wblog(F,L,"ERR %s() "
         "got different length in qs (%d/%d)",FCT,this->qs.len,S.qs.len);
      return 0;
   }

   if (cgd.SIZE!=S.cgd.SIZE) {
      if (F) wblog(F_L,"ERR %s() got different size in cgd\n"
         "(%s; %s)",FCT, cgd.sizeStr().data, S.cgd.sizeStr().data);
      return 0;
   }

   return 1;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::initOM( const char *F, int L) {

   unsigned r=this->rank(F_L);
   wbvector<SPIDX_T> &S=this->cgd.SIZE;

   this->cgd.wbsparray<TD>::checkSize(F_LF);
   if (!r) wblog(FL,
      "ERR %s() got empty CGDdata (%s)",FCT,sizeStr().data);

   if (S.len==r) { S.Append(1); } else
   if (S.len!=r+1) wblog(FL,
      "ERR %s() invalid size (len=%d/%d)",FCT,S.len,r);

   return *this;
};


template <class TQ, class TD>
double CData<TQ,TD>::normDiff(
  const char *F, int L, const CData<TQ,TD> &S, double eps
) const {

   double e, emax=0; sameType(S,F_L);

   e=this->qs.normDiff(S.qs);
      if (e>eps) return e;
      if (e>emax) emax=e;
   e=cgd.normDiff2(S.cgd); e=std::sqrt(std::fabs(e));

   if (e>eps) return e;
   if (e>emax) emax=e;

   return emax;
};


template <class TQ, class TD>
double CData<TQ,TD>::norm2(unsigned k) const {

   if (int(k)<0) { return this->cgd.D.norm2(); }

   unsigned m=getOM(FL);

   if (k>=m) wblog(FL,"ERR index out of bounds (%d/%d)",k,m);
   if (cgd.D.len!=cgd.IDX.dim1) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,cgd.IDX.dim1,cgd.D.len);

   if (m<2) { return cgd.D.norm2(); }
   else {
      double x2=0; const size_t *i4=cgd.IDX.data+3;
      for (SPIDX_T i=0; i<cgd.D.len; ++i, i4+=cgd.IDX.dim2) {
          if ((*i4)==k) x2+=(CONJ(cgd.D[i])*cgd.D[i]);
      }
      return x2;
   }
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::Sort(const char *F, int L) {

   checkConsistency(F_L); if (isEmpty(F_L)) {
      wblog(F_L,"WRN Sort() got empty CData !??");
      return *this;
   }

   if (this->t.isAbelian(F_L)) {
      if (cgd.SIZE.len!=1 || cgd.SIZE[0]!=1 || cgd.D.len!=1 || cgd.D[0]!=1)
      wblog(F_L,"ERR %s() invalid %s CData (%s)",
      FCT, qStr().data, cgd.sizeStr().data);
   }
   else {
      wbMatrix<TQ> Z; wbperm P;
      unsigned i,r, n=this->t.qlen(); const TQ *q=this->qs.data;

      if (!isrank(3,&r)) wblog(FL,
         "ERR %s() got rank-%d CData",FCT,r);

      for (i=0; i<3; ++i, q+=n) { Z=gRG.getZ(F_L,this->t,q);
         Z.sortRecs_float(P,-1);
         if (!P.isIdentityPerm()) wblog(FL,
            "WRN %s() rather permute gRG data [check source]",FCT); 
         cgd.Select0(P,i);
      }
      cgd.Compress();
   }

   return *this;
};


template <class TQ, class TD>
template <class DB>
bool CData<TQ,TD>::sameSizeR(
   const CData<TQ,DB> &B, unsigned *r_, const char *F, int L
 ) const {

   unsigned i=0, ra=this->rank(F_L), rb=B.rank(F_L);
   const SPIDX_T *sa=this->cgd.SIZE.data, *sb=B.cgd.SIZE.data;

   if (ra!=rb || ra<2 || this->t!=B.t) {
      if (F) wblog(F_L,
         "ERR %s() invalid input CG sets (%s <> %s; %s <> %s; %d/%d)",
          FCT, this->qStr().data, B.qStr().data,
          this->cgd.sizeStr().data, B.cgd.sizeStr().data, ra, rb);
      return 0;
   }

   for (i=0; i<ra; ++i) { if (sa[i]!=sb[i]) return 0; }

   if (r_) { (*r_)=ra; }

   return 1;
};


template <class TQ, class TD>
template <class T2>
int CData<TQ,TD>::sameSizeR(
   const char *F, int L, const wbvector<T2> &S
 ) const {

   if (cgd.isEmpty()) {
      if (isAbelian() || !S.len) return 1;
      return 0;
   }

   unsigned i=0, ra=this->rank(FL), m=gotOM(FL), m2=1;
   const SPIDX_T *sa=this->cgd.SIZE.data;

   if (S.len==ra+1) { m2=S[ra]; } else
   if (ra<2 || S.len<ra || S.len>ra+1) {
      if (F) wblog(F,L,
         "ERR %s() CData rank mismatch (%s: %s <> %s; %d/%d)",
          FCT, this->qStr().data,
          this->cgd.sizeStr().data, S.toStr().data, ra, S.len);
      return -1;
   }

   if (sa) {
      for (; i<ra; ++i) { if (S.data[i]!=sa[i]) return 0; }
   }
   else {
      SPIDX_T s=cgd.dim();
      if (ra!=2) wblog(FL,
         "ERR %s() invalid CData\n%s",FCT,toStr().data);
      for (; i<ra; ++i) { if (S.data[i]!=s) return 0; }
   }

   if (!m) m=1; else
   if (m2>m) return 2; else
   if (m2<m) return 3;

   return 1;
};


template <class TQ, class TD>
template <class DB>
bool CData<TQ,TD>::sameSizeR(
   const cdata<DB> &b, unsigned *r_, const char *F, int L
 ) const {

   unsigned i=0, ra=this->rank(F_L), rb=b.rank(F_L);
   const SPIDX_T *sa=this->cgd.SIZE.data, *sb=b.SIZE.data;

   if (ra!=rb || ra<2) {
      if (F) wblog(F_L,
         "ERR %s() invalid input CG sets (%s: %s <> %s; %d/%d)",
          FCT, this->qStr().data,
          this->cgd.sizeStr().data, b.sizeStr().data, ra, rb);
      return 0;
   }

   for (i=0; i<ra; ++i) { if (sa[i]!=sb[i]) return 0; }

   if (r_) { (*r_)=ra; }

   return 1;
};


template <class TQ, class TD>
bool CData<TQ,TD>::isSymmmetric() const {

   unsigned r=this->qdir.len;

   if (r%2) return 0;
   if (!r) {
      if (cgd.SIZE.len || cgd.IDX.dim1 || cgd.D.len>1)
         wblog(FL,"ERR %s() got invalid empty CGC !?",FCT);
      return 1;
   }

   unsigned i=0, r2=r/2;
   const char *qd=this->qdir.data;

   for (; i<r2; ++i) { if (qd[i]<=0 || qd[i+r2]>=0) return 0; }

   r2=this->qs.len/2;
   return (
      memcmp(this->qs.data, this->qs.data+r2, r2*sizeof(TQ)) ? 0 : 1
   );
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::AddMultiplicity(
   const char *F, int L, const cdata<TD> &c
){
   SPIDX_T i,l,m; unsigned n, r=0;
   wbsparray<TD> X;
   const cdata<TD> &a=this->cgd;

   a.checkSize(F_LF,"A");
   c.checkSize(F_LF,"C");

   if (!sameSizeR(c,&r)) wblog(F_L,
      "ERR %s() got incompatible CData\n%s: %s <> %s",FCT,
      this->toStr().data, a.sizeStr().data, c.sizeStr().data);
   if (r<3) wblog(FL,"ERR %s() unexpected rank-%d CData",FCT,r);

   wbvector<SPIDX_T> s(r+1);

   if (a.IDX.dim1 && CG_VERBOSE>2) wblog(PFL,
      "[+] %s() %s",FCT,toStr().data);

   for (i=0; i<a.SIZE.len; ++i) s[i]=a.SIZE[i];
   for (; i<s.len; ++i) s[i]=1;
   m=(s[r]++);

   X.init(s, a.IDX.dim1 + c.IDX.dim1);

   MEM_CPY<TD>(X.D.data,         a.D.len, a.D.data);
   MEM_CPY<TD>(X.D.data+a.D.len, c.D.len, c.D.data);

   for (n=a.IDX.dim2, i=0; i<a.IDX.dim1; ++i) {
      X.IDX.recSetP(i,a.IDX.rec(i),n);
   }
   for (n=c.IDX.dim2, l=i, i=0; i<c.IDX.dim1; ++i) {
      X.IDX.recSetP(l,c.IDX.rec(i),n);
      X.IDX(l++,r)=m;
   }

#ifndef NSAFEGUARDS
   { wbindex icm(r,'i');
     wbsparray<TD> E; X.contract(FL,icm,X,icm,E);
     TD nrm=TD(1);
     if (!E.isProptoId(nrm, CG_EPS1)) {
        MXPut(FL,"a").add(*this,"A")
        .add(c,"c").add(X,"X").add(E,"E").add(icm,"ic");
        wblog(FL,"ERR %s() got non-orthonormal OM space !?",FCT);
     }
   }
#endif

   X.save2(this->cgd);

   cstat.update_m();
   return *this;
};


template <class TQ, class TD>
bool CData<TQ,TD>::isBasicCG(const char *F, int L) const {


   if ((this->rank(F_L))!=3 || this->qdir!="++-") {
      if (F) wblog(F_L,
         "ERR invalid basic CG set (%s; %s)",
         this->sizeStr().data, this->qdir.toStr().data);
      return 0;
   }
   return 1;
};


template <class TQ, class TD>
SPIDX_T CData<TQ,TD>::getBasicCGSize(unsigned *m) const {

   unsigned r=rank(FL); SPIDX_T n=this->cgd.SIZE.len;
   const SPIDX_T *s=this->cgd.SIZE.data;

   if (n==r) { if (m) (*m)=1; }
   else {
      if (n!=r+1) wblog(FL,"ERR %s() got %d/%d !??",FCT,n,r+1);
      if (m) { (*m)=s[r]; }
   }

   if (r) { n=s[0]; for (unsigned i=1; i<r; ++i) n*=s[i]; }
   else { n=0; }

   return n;
};


template <class TQ, class TD>
wbsparray<TD>& CData<TQ,TD>::getBasicCG(
   unsigned k, wbsparray<TD> &a
 ) const {

   unsigned om=this->getOM(FL), r=this->rank(FL);

   if (k>=om) wblog(FL,"ERR OM out of bounds (k=%d/%d)",k+1,om);
   if (this->cgd.SIZE.len==r) {
      if (om!=1) wblog(FL,"ERR cdata inconsistency (%d)",om);
      a.init(this->cgd);
   }
   else {
      const cdata<TD> &c = this->cgd;
      wbvector<size_t> S(r,c.SIZE.data);

      const size_t *im=c.IDX.data+r;
      SPIDX_T i=0, l=0, N=c.D.len; unsigned m=c.IDX.dim2;

      for (; i<N; ++i, im+=m) { if ((*im)==k) { ++l; }}
      a.init(S,l);

      im=c.IDX.data+r;
      for (l=i=0; i<N; ++i, im+=m) { if ((*im)==k) {
         a.setRecP(l++, im-r, c.D[i]);
      }}
   }
   return a;
};


template <class TQ, class TD>
int CData<TQ,TD>::getCG_set(
   const char *F, int L, wbIndex &Idx, wbsparray<TD> &a
 ) const {

   const wbvector<SPIDX_T> &S2=this->cgd.SIZE;
   unsigned nq=this->t.qlen(), r=this->qdir.len, r2=S2.len, rm=r2-r;

   if (!r2) { r2=this->cgd.rank(); rm=r2-r; }

   if (!r || !r2 || r2<r || !nq || this->qs.len/nq!=r || this->qs.len%nq)
      wblog(F_L,"ERR %s() got invalid/inconsistent QSet\n{%d,%d,%d}: %s",
      FCT, r,r2,nq, ((QSet<TQ>*)this)->toStr().data
   );

   if (Idx.isEmpty()) {
      if (r2>r) { Idx.init(rm,S2.data+r); }
      else { a.init(cgd);
         if (cgd.D.len) {
            INDEX_T s=1; Idx.init(1,&s); ++Idx;
            return 1;
         }
         else return 0;
      }
   }
   else {
      if ((r2!=r && r2!=r+Idx.len) || (r2==r && Idx.len!=1)) wblog(F_L,
         "ERR %s() got rank mismatch in OM (%d+%d <> %d !?)",
         FCT,r,Idx.len,r2
      );
   }

   if (!(++Idx)) { a.init(); return 0; }
   Idx.checkValid(FL);

   const wbMatrix<SPIDX_T> &IDX=this->cgd.IDX;
   const wbvector<TD> &D=this->cgd.D;
   wbvector<SPIDX_T> S(r,S2.data);

   if (IDX.dim1!=D.len || IDX.dim2!=S2.len) wblog(FL,
      "ERR %s() cdata mismatch (%dx%d <> %d x %d)",
      FCT,IDX.dim1,IDX.dim2,D.len,S2.len);
   if (!IDX.dim1 || !IDX.dim2) wblog(FL,
      "ERR %s() got empty cdata (%dx%d <> %d x %d)",
      FCT,IDX.dim1,IDX.dim2,D.len,S2.len);

   INDEX_T
      i1=Wb::findfirst_sorted(NULL,0,
         Idx.data, IDX.data+r, rm, IDX.dim1,IDX.dim2,-1),
      i2=Wb::findlast_sorted(NULL,0,
         Idx.data, IDX.data+r, rm, IDX.dim1,IDX.dim2,-1),
      i=0, l=i2-i1+1;

   if (i1>i2) wblog(FL, 
      "ERR %s() got unsorted data (%d,%d) !?",FCT,i1,i2);

   if (i2>=IDX.dim1) { 
      if (i1<IDX.dim1) wblog(FL,"ERR %s() i1=%d, i2=%d",FCT,i1,i2);
      a.init(S,0); return 1;
   }


   a.init(S,l);
   for (; i<l; ++i, ++i1) {
      a.setRecP(i, IDX.rec(i1), D[i1]);
   }

   return 1;
};


template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::getJ1Symbol(
   const char *F, int L, CData<TQ,TD> &C) const {

   unsigned i=0, m=this->qs.len/3, l=2*m;
   const cdata__ &c=this->cgd;

   if (!this->isStd3()) wblog(F_L,"ERR %s() "
      "got invalid CData %s",FCT, ((QSet<TQ>*)this)->toStr().data);
   if (c.SIZE.len!=3 || c.SIZE[2]!=1) wblog(FL,"ERR %s() "
      "incompatible cgdata (%s)",FCT, c.sizeStr().data);
   if (c.IDX.dim2!=3) wblog(FL,"ERR %s() invalid cgdata (%s; %s)",
      FCT, c.IDX.sizeStr().data, c.sizeStr().data);
   for (; i<m; ++i) {
      if (this->qs.data[l+i]) {
         wblog(FL,"ERR %s() incompatible CData\n%s",
         FCT, ((QSet<TQ>*)this)->toStr().data);
      }
   }

   C=(const QSet<TQ>&)(*this); {
      C.cstat.init(CGD_J1SY_CGC);
      C.qs.len=l;
      C.qdir.len=2;
      cgd.skipTrailingSingletons(FL,C.cgd,2);
   }
   return C;
};


template <class TQ, class TD>
int CData<TQ,TD>::init(
   const char *F, int L, const mxArray *S, unsigned k,
   char bflag,
   char check
){
   if (!S) {
      init(); return 0;
   }
   if (!mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid CData (%s)",FCT,mxGetClassName(S));

   int itype = mxGetFieldNumber(S,"type"),
       iqset = mxGetFieldNumber(S,"qset"),
       iqdir = mxGetFieldNumber(S,"qdir"),
       icgd  = mxGetFieldNumber(S,"cgd" ),
       id    = mxGetFieldNumber(S,"cid" );

   if (itype<0 || iqset<0 || iqdir<0 || icgd<0 || id<0) wblog(F_L,
      "ERR invalid CData structure (%d,%d,%d,%d,%d)",
      itype,iqset,iqdir,icgd,id
   );

   unsigned n=mxGetNumberOfElements(S); int q=0;
   CDATA_TQ B;

   if (!n) wblog(FL,"ERR %s() got empty CData structure",FCT);
   if (int(k)>=0) {
      if (k>=n) wblog(FL,"ERR index out of bounds (%d/%d)",k+1,n); }
   else k=0;

   B.t    .init(FL,mxGetFieldByNumber(S,k,itype));
   B.qs   .init(FL,mxGetFieldByNumber(S,k,iqset),0,'!');
   B.qdir .init(FL,mxGetFieldByNumber(S,k,iqdir));
   B.cgd  .init(FL,mxGetFieldByNumber(S,k,icgd),0,B.t);
   B.cstat.init(FL,mxGetFieldByNumber(S,k,id));

   if (check) B.checkNormSign(FL);

   if (isEmpty()) {
      B.save2(*this); q=1;
   }
   else if (!bflag) {
      B.save2(*this); q=2;
   }
   else {
      if (B.cstat==CGD_REF_INIT) wblog(FL,
         "ERR %s() got REF_INIT for loaded CData !?\n%s",
         FCT,B.toStr().data);

      if (!sameAs(B,'l')) {
         MXPut(FL,"a").add(*this,"A").add(B,"B");
         wblog(FL,"ERR %s() "
           "got QSet mismatch%N%N   %s  |  %s%N<> %s  |  %s%N", FCT,
              toStr().data,   cstat.toStr('V').data,
            B.toStr().data, B.cstat.toStr('V').data
         );
      }

      if (cstat!=B.cstat) {
         if (this->olderThan(B)) {
            B.save2(*this); q=3;
         }
         else q=4;
      }
      else {
         if (cstat==CGD_REF_INIT) {
            B.save2(*this); q=5;
         }
         else q=6;
      }

      if (cstat==CGD_REF_INIT) wblog(FL,
      "ERR %s() got REF_INIT from RCStore !?\n%s",FCT,toStr().data);
   }

   return q;
};


template <class TQ, class TD>
mxArray* CData<TQ,TD>::toMx() const {

   mxArray *S=mxCreateStruct(1,1);
   add2MxStruct(S,0); return S;
};


template <class TQ, class TD>
mxArray* CData<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

#ifdef CG_CHECK_MW_PERM
   const char *fields[]={"type","qset","qdir","cid","cgd","P0"};
   return mxCreateStructMatrix(m,n,6,fields);
#else
   const char *fields[]={"type","qset","qdir","cid","cgd"};
   return mxCreateStructMatrix(m,n,5,fields);
#endif
};


template <class TQ, class TD>
void CData<TQ,TD>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      unsigned s=0; int k=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (k=mxGetFieldNumber(S,"qdir"))<0) wblog(FL,
      "ERR %s() missing mxCreateStruct() !??\n%lx, %d/%d, %d",
       FCT,S,i+1,s,k);
   }
   else if (!S) wblog(FL,"ERR %s() got non-initialized S",FCT);

   mxSetFieldByNumber(S,i,0, this->t.toMx());
   mxSetFieldByNumber(S,i,1, this->qs.toMx());
   mxSetFieldByNumber(S,i,2, this->qdir.toMx());
   mxSetFieldByNumber(S,i,3, cstat.toMx());
   mxSetFieldByNumber(S,i,4, cgd.toMx());

#ifdef CG_CHECK_MW_PERM
   if (gRG.buf.find(this->t)!=gRG.buf.end()) {
   if (this->qdir=="++-" && rank(FL)==3) {
      unsigned n=this->t.qlen();

      const genRG_base<TQ,RTD> &R =
         gRG.getR(FL,this->t,this->qs.data+2*n);
      mxSetFieldByNumber(S,i,5, R.P0.toMx());
   }}
#endif
};


template <class TQ, class TD>
void CData<TQ,TD>::info(const char *istr, const char *F, int L) const {

   unsigned e=0, n=this->t.qlen(), r=rank(F_L);

   checkConsistency(F_L);

   if (F) wblog(F,L,
      "CG Space %s :: %s",istr && istr[0] ? istr : "",qStr());
   else printf(
      "\nCG Space %s :: %s\n",istr && istr[0] ? istr : "",qStr());
   printf("\n");

   printf("  Q=[%s] (rank=%d)\n",this->QStr().data,r);

   if (cgd.SIZE.len)  cgd.info();

   if (e) wblog(FL,"ERR (%d)",e);
   printf("\n");
};


template <class TQ, class TD>
wbstring CData<TQ,TD>::toStr(char vflag) const {

   if (isEmpty()) return "(empty)";

   const unsigned n=127; char s[n+1];

   unsigned l=0, m=1;
   if (cgd.SIZE.len>this->qdir.len) m=cgd.SIZE[this->qdir.len];

   if (cstat.t>=CGD_NUM_TYPES) wblog(FL,"ERR %s() "
      "ctype out of range (%d/%d)",FCT,cstat.t,CGD_NUM_TYPES
   );


   l=snprintf(s,n,"%s",((const QSet<TQ>*)this)->toStr().data);
   if (!vflag && m>1 && l<n) { l+=snprintf(s+l,n-l,"[%d]",m); }
   if (cstat.t!=CGD_DEFAULT && l<n) {
      l+=snprintf(s+l,n-l," %s",CGD_TYPE_STR[cstat.t]);
   }

   if (vflag && l<n) {
      l+=snprintf(s+l,n-l," %s",cgd.sizeStr().data);
   }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};


template <class TQ, class TD>
wbstring CData<TQ,TD>::sizeStr() const {

   unsigned r=rank(FL);

   if (!r) {
      if (cgd.SIZE.len) wblog(FL,
         "ERR %s() got r=%d/%d !??",FCT,r,cgd.SIZE.len);
      return "(empty)";
   }
   else {
      unsigned l=0, n=128; char s[n];

      l+=snprintf(s+l,n-l,"%s", cgd.sizeStr().data);
      if (r<cgd.SIZE.len && l<n) {
         l+=snprintf(s+l,n-l," @ %d", unsigned(cgd.SIZE[r]));
      }
      if (l>n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
      return s;
   }
};



template <class TQ> inline
CRef<TQ>& CRef<TQ>::initBase(
   const CDATA_TQ* cgr_, unsigned i, unsigned m
){
   if (cgr_ && cgr_->t.isNonAbelian()) {
      cgr=(CDATA_TQ*)cgr_;
      if (m>(cgr->getOM(FL))) wblog(FL,"ERR %s() "
         "OM range out of bounds (%d/%d)",FCT,m,cgr->getOM(FL));
      if (i>=m) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,i+1,m);
      cgw.init(m); cgw[i]=1;
   }
   else {
      cgr=NULL; if (cgr_) cgr_->checkQ_abelian(FL);
      if (m!=1 || i) wblog(FL,
         "ERR %s() invalid abelian cgref (%d/%d)",FCT,i,m);
      cgw.init();
   }
   cgp.init(); conj=0; rtype=0;
   return *this;
};


template <class TQ> inline
unsigned CRef<TQ>::qdim(unsigned k, const QType *qt) const {

   if (isAbelian()) {
      if (qt && !qt->isAbelian()) wblog(FL,"ERR %s() "
         "abelian symmetry mismatch (%s)",FCT,qt->toStr().data);
      return 1;
   }
   if (!cgr) wblog(FL,"ERR %s() missing CData for %s",FCT);
   if (qt && (*qt)!=cgr->t) wblog(FL,"ERR %s() abelian symmetry "
      "mismatch (%s/%s)",FCT,cgr->t.toStr().data,qt->toStr().data);

   const wbvector<unsigned> &s=cgr->cgd.SIZE;
   if (k>=s.len) {
      if (!s.len || s.len==1 && s[0]==1) return 1;
      else wblog(FL,
         "ERR %s() dimension out bounds (%d/%d)",FCT,k,s.len
      );
   }

   return s[getP(k)];
};


template <class TQ>
CRef<TQ>& CRef<TQ>::init(
   const char *F, int L, const mxArray *a, unsigned k,
   char refC, const QSet<TQ> *QS
){
   if (!a || mxIsEmpty(a)) {
      if (k) wblog(FL,
         "ERR CRef::%s() got k=%d for null mxArray !?\n%s",FCT,k,
         QS ? QS->toStr().data : "(no QS specified)");
      return init();
   }

   unsigned n=mxGetNumberOfElements(a);

   int lflag,
      ix=-1,"cgt" ) // see below
      it=mxGetFieldNumber(a,"type"),
      iq=mxGetFieldNumber(a,"qset"),
      io=mxGetFieldNumber(a,"qdir"),
      id=mxGetFieldNumber(a,"cid" ),
      iw=mxGetFieldNumber(a,"cgw" ),
      is=mxGetFieldNumber(a,"size");

   wbvector<SPIDX_T> S;
   cgdStatus st;
   QSet<TQ> Q;

   if (it<0 || iq<0 || io<0 || id<0 || iw<0 || is<0)
      wblog(FL,"ERR %s() invalid CRef (%d %d %d; %d %d %d)",
      FCT, it,iq,io, id,iw,is);
   if (k>=n) wblog(F,L,
      "ERR CRef::%s() index out of bounds (%d/%d)",FCT,k,n);

   const mxArray
      *as=(is>=0 ? mxGetFieldByNumber(a,k,is) : NULL),
      *ad=(id>=0 ? mxGetFieldByNumber(a,k,id) : NULL),
      *aw=(iw>=0 ? mxGetFieldByNumber(a,k,iw) : NULL),
      *at=(it>=0 ? mxGetFieldByNumber(a,k,it) : NULL),
      *aq=(iq>=0 ? mxGetFieldByNumber(a,k,iq) : NULL),
      *ao=(io>=0 ? mxGetFieldByNumber(a,k,io) : NULL);

   init();
   cgw.init(F_L,aw);

   if (!as || !ad || !at || !aq || !ao) { int e=0;
      if (as && mxGetNumberOfElements(as)) e=1; else
      if (ad && mxGetNumberOfElements(ad)) e=2; else
      if (at && mxGetNumberOfElements(at)) e=3; else
      if (aq && mxGetNumberOfElements(aq)) e=4; else
      if (ao && mxGetNumberOfElements(ao)) e=5;
      if (e) wblog(FL,
         "ERR %s() got invalid CRef (CGR_ABELIAN; e=%d)",FCT,e);

      rtype=CGR_ABELIAN;
      return *this;
   }

   S     .init(F_L,as);
   st    .init(F_L,ad);
   Q.t   .init(F_L,at);
   Q.qs  .init(F_L,aq,0,'!');
   Q.qdir.init(F_L,ao);

   lflag=S.isEmpty();
   if (lflag ^ st.isEmpty()) wblog(FL,
      "ERR %s() got invalid minimal CRef data (%d,%d)\n"
      "(empty size and cid required to indicate initialization)!",FCT,
      lflag,st.toStr('v').data
   );

   if (QS) {
      if (Q!=(*QS)) wblog(FL,
         "ERR %s() got CGR QSet mismatch\n%s <> %s",
         FCT, QS->toStr().data, Q.toStr().data
      ); 
   }

   Q.Sort(&cgp,&conj,'i');
   CDATA_TQ &C=(
      refC && !lflag ? gCG.getBUF(0,0,Q,0) : gCG.getBUF(FL,Q,'L')
   );

   char gotC=(C.isEmpty() ? 0 : (C.isRefInit() ? -1 : 1));
   unsigned r=Q.rank();

   if (lflag) {
      if (!gotC) wblog(FL,"ERR %s() got minimal CRef "
         "for non-existing CData !?\n%s",FCT,Q.toStr().data);
      if (C!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%s <> %s)",FCT,C.toStr().data,Q.toStr().data);
      if (!cgw.len || cgw.len>C.getOM()) wblog(FL,
         "ERR %s() OM out of bounds (%d/%d)",FCT,cgw.len,C.getOM());

      cgr=(&C);

      if (CG_VERBOSE>2) wblog(FL,
         " *  got minimal CRef %s",toStr().data);
      return *this;
   }

   if (st==CGD_DEFAULT) wblog(FL,
      "ERR %s() got invalid CData status\n%s",FCT,st.toStr('v').data);
   if (Q.t.isAbelian()) {
      if (st.ctime) wblog(FL,
         "WRN %s() got cstat for abelian symmetry !?\n%s\n%s",
         FCT,Q.toStr().data,st.toStr('V').data
      );
   }
   else {
      if (!st.ctime) wblog(FL,
         "WRN %s() missing cstat for non-abelian symmetry !?\n%s\n%s",
         FCT,Q.toStr().data,st.toStr('V').data
      );
   }

   if (cgp.len && !cgp.isIdentityPerm()) {
      wbvector<SPIDX_T> Sx(S);
      if (cgp.len+1<S.len || cgp.len>S.len) wblog(FL,
         "ERR %s() invalid S.len=%d/%d",FCT,S.len,cgp.len);
      for (unsigned i=0; i<cgp.len; ++i) { Sx.el(cgp[i])=S[i]; }
      Sx.save2(S);
   }

   rtype=CGR_DEFAULT;

   if (!gotC) {
      if (!S.len) wblog(FL,"ERR %s() got empty cgr.size",FCT);

      C=Q; C.cstat=st;

      if (C.t.isAbelian()) {
         if (!S.allEqual(1)) wblog(FL,
            "ERR %s() invalid scalar size (%s)",FCT,C.cgd.sizeStr().data);
         C.cstat.init(CGD_ABELIAN);
         C.cgd.wbsparray<RTD>::init();
      }
      else if (refC) { 
         C.cstat.t=CGD_REF_INIT;
         ix=mxGetFieldNumber(a,"cgt");
         if (ix<0) C.cgd.init(S);
      }
      else wblog(PFL,
         "ERR %s() got CRef to undefined CData (%s)\n%s",
         FCT,PROG,C.toStr().data
      );
   }
   else if (gotC>0) {
      if (C.cstat.t==CGD_REF_INIT) wblog(FL,
         "ERR %s() cstat.t=%s !?",FCT,C.cstat.tstr());
      if (C!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%s <> %s)",FCT,C.toStr().data,Q.toStr().data);

      int q=C.sameSizeR(FL,S);
      if (q<=0) wblog(FL,
         "ERR %s() CRef size mismatch %s -> %s (%d)\n%s",FCT,
         C.sizeStr().data, S.toStr().data, cgw.len, C.toStr().data);
      if (q==2 || (r==2 && q!=1)) wblog(FL,
         "ERR %s() got input CRef with larger OM %s -> %s (%d)\n%s",
         FCT,C.sizeStr().data,S.toStr().data,cgw.len,C.toStr().data);

      if ((q=C.cstat.cmp(FL,st))) { int e=0;
         if (r<=2) e=1;
         if (q<0) e=2;
         if (e) wblog(PFL,
            "ERR %s() unexpected CData status mismatch (e=%d)"
            "\n   %s\n<> %s%N%N   %s%N-> %s%N", FCT, e, Q.toStr().data,
            C.toStr().data, st.toStr('V').data, C.cstat.toStr('V').data
        );
      }
   }
   else {
      if (C.cstat.t!=CGD_REF_INIT || !refC) wblog(FL,
         "ERR %s() got isref CData.cstat='%s' (%d) !?",
         FCT,C.cstat.toStr().data,refC);
      if (C!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%d <> %s)",FCT,C.toStr().data,Q.toStr().data);

      int q=C.sameSizeR(FL,S);
      if (q<=0) wblog(FL,
         "ERR %s() CRef size mismatch %s -> %s (%d)\n%s",FCT,
         C.sizeStr().data, S.toStr().data, cgw.len, C.toStr().data);

      if (!refC) {
      if ((q==1 && C!=st && st.t!=CGD_REF_INIT) ||
          (q!=1 && !C.cstat.sameAs(st,2))
        ) wblog(FL,
          "WRN %s() got CRef status mismatch\n   %s\n-> %s%N%N   "
          "%s-> %s", FCT, Q.toStr().data, ((QSet<TQ>&)C).toStr().data,
          C.cstat.toStr('V').data, st.toStr('V').data);
      }

      if (q==2) {
         ix=mxGetFieldNumber(a,"cgt");
      }
   }

   if (ix>=0) {
      wbvector<RTD> cgt(F_L,mxGetFieldByNumber(a,k,ix));

      if (C.cstat.t!=CGD_REF_INIT) wblog(FL,
         "ERR %s() got non-ref CData !?\n%s",FCT,C.toStr().data);
      if (r && S.len>r) {
         if (C.cgd.SIZE.len==S.len && C.cgd.SIZE[r]>S[r]) wblog(FL,
            "ERR %s() got reduced OM %d/%d !?",
            FCT,C.sizeStr().data,S.toStr().data
         );
      }

      C.RefInit_auxtr(FL,S,cgt);
   }

   cgr=(&C);

   if (cgw.len>cgr->getOM()) wblog(FL,
      "ERR %s() OM out of bounds (%d/%d)",FCT,cgw.len,cgr->getOM());

   return *this;
};


template <class TQ>
CRef<TQ>& CRef<TQ>::initDecompose(const char *F, int L,
   cdata__ c,
   CDATA_TQ *cgr_
){

   if (cgr_!=NULL) { cgr=cgr_; } else
   if (!cgr) wblog(FL,
      "ERR %s() CRef not yet initialized !??",FCT);
   if (cgr->t.isAbelian()) wblog(FL,
      "WRN %s() got Abelian symmetry !?\n%s",FCT,toStr().data);

   RTD x; double e=0;

   if (cgr->cgd.isEmpty()) {

      x=c.NormSignC(F,L);

      e=fabs(double(x)); if (e<1E-3) {
         MXPut(FL,"a").add(c,"c").add(double(x),"x");
         wblog(FL,"WRN %s() got small new coefficient !?? (%.3g)",FCT,e);
      }
      if (cgw.len) wblog(FL,
         "ERR %s() got cgw.len=%d for empty CData !?",FCT,cgw.len);

      cgw.init(1,&x);
      c.save2(((CDATA_TQ*)cgr)->cgd);

      cgp.init(); conj=0;

      gStore.save_CData(FL,*cgr);
      if (cgr->qdir.len==3 && cgr->isStd3()) wblog(FL,
         "WRN %s() decomposing std rank-3 CData !?\n%s",
         FCT,cgr->toStr().data
      );

      if (CG_VERBOSE>2 || (CG_VERBOSE>1 && cgr->rank(FL)<4)) wblog(PFL,
         "[+] CBUF[%d] new/dec: %s",gCG.BUF.size(),cgr->toStr().data);

      return *this;
   }

   unsigned r=cgr->rank(), m=cgr->getOM();

   if (r==2) {
      cgw.init(1); cgp.init(); conj=0;
      if (m!=1 || c.sameUptoFac(cgr->cgd, cgw.data)!=0) wblog(FL,
         "ERR %s() invalid scalar CData\n%s",FCT,cgr->toStr().data);
      return *this;
   }

   if (r<2) wblog(FL,
      "ERR %s() invalid rank-%d CData",r);
   if (c.SIZE.len!=r) wblog(FL,
      "ERR %s() got OM in input CData (%d/%d) !??",FCT,c.SIZE.len,r);

   cgr->checkNormSign(F_L,'x');

   if (m>1) { cdata__ X;
      wbindex
         ia(r,'i'),
         i1(1),
         im(1); im[0]=r;

      if (cgr->cgd.SIZE.len!=r+1) wblog(FL,
         "ERR %s() invalid CData (%d/%d)",FCT,cgr->cgd.SIZE.len,r);

      for (unsigned it=0; it<2; ++it) {

         c.cgsparray::contract(F_L,ia,cgr->cgd,ia, X);
         if (!X.isVector()) wblog(FL,
            "ERR %s() got rank-%d object (%s)",FCT,X.sizeStr().data);

         cgr->cgd.cgsparray::contract(F_L,im,X,i1,c,wbperm(),-1,1);

         if (it==0) {
            cgw.initT(X);
            e=c.norm(); if (e<CG_SKIP_EPS2) { break; }
         }
      }
   }
   else {
      for (unsigned it=0; it<2; ++it) {
         x=c.cgsparray::dotProd(F_L,cgr->cgd);

         c.Plus(FL,cgr->cgd,-x);
         if (it==0) {
            cgw.init(1,&x);
            e=c.norm(); if (e<CG_SKIP_EPS2) { break; }
         }
      }
   }

   if (e<1E-8) {
      if (e>CG_SKIP_EPS2) wblog(FL,
         "ERR %s() CGC ortho @ %.3g",FCT,e);
      return *this;
   }

   x=c.NormSignC(F,L);

   e=fabs(double(x)); if (e<1E-3) wblog(FL,
      "WRN %s() got small new coefficient !?? (%.3g)",FCT,e);

   cgw.Append(x);

   ((CDATA_TQ*)cgr)->AddMultiplicity(FL,c);
   ((CDATA_TQ*)cgr)->cstat.t=CGD_FROM_DEC;

   if (CG_VERBOSE>2 || (CG_VERBOSE>1 && r<4)) wblog(PFL,
      "[u] CBUF[%d] new/dec: %s",gCG.BUF.size(),cgr->toStr().data);
   cgr->checkConsistency(FL);

   if (r==2) {
      if (!isIdentityCGC()) wblog(FL,"ERR %s() "
         "got non-scalar rank-%d CGC (%s) !??",FCT,r,cgr->sizeStr().data);
      else wblog(FL,"ok. %s() got scalar rank-%d CGC\n"
         "%s @ %s",FCT, r, cgr->sizeStr().data, cgr->toStr().data
      );
   }

   gStore.save_CData(FL,*cgr);
   if (cgr->qdir.len==3 && cgr->isStd3()) wblog(FL,"WRN %s() "
      "decomposing rank-3 CData !?\n%s",FCT,cgr->toStr().data);

   return *this;
};


template <class TQ>
CRef<TQ>& CRef<TQ>::initIdentityR(
   const char *F, int L, const QType &t, const TQ *qs, unsigned dim,
   char xflag
){
   char isa=t.isAbelian(); cgp.init(); conj=0;

   if (isa && !xflag) {
      cgw.init(); cgr=NULL; rtype=CGR_ABELIAN;
   }
   else {
      cgr=&gCG.getIdentityC(F_L,t,qs,dim);
      cgw.init(1); 

      if (!cgr->cgd.D.len) {
         if (!isa) wblog(FL,
            "ERR %s() got invalid CData\n%s",FCT,cgr->toStr().data);
         cgw[0]=1;
      }
      else {
         cgw[0]=RTD(1)/RTD(cgr->cgd.D.data[0]);
      }
   }

   return *this;
};


template <class TQ>
CRef<TQ>& CRef<TQ>::initIdentityZ(
   const char *F, int L, const QType &t, const TQ *qs, unsigned dim,
   char xflag
){
   cgp.init(); conj=0;

   if (t.isAbelian() && !xflag) {
      cgw.init(); cgr=NULL; rtype=CGR_ABELIAN;
   }
   else {
      cgr=&gCG.getIdentityZ(F_L,t,qs,dim);
      if (int(dim)<=0) {
         dim=cgr->cgd.dim();
         if (!dim) wblog(FL,"ERR %s() got dim=%d !?",FCT,dim);
      }
      cgw.init(1); cgw[0]=SQRT(RTD(dim));
   }

   return *this;
};


template <class TQ>
CRef<TQ>& CRef<TQ>::Reduce2Identity(char xflag) {

   if (!cgr) {
      if (cgp.len || cgw.len) wblog(FL,
         "ERR %s() invalid CRef\n%s",FCT,toStr().data);
      return *this;
   }

   unsigned r=rank(FL);

   if (r<=2) {
      if (r<2)
           wblog(FL,"ERR %s() got rank-%d CRef !?",FCT,r);
      else wblog(FL,"WRN %s() got rank-%d CRef !?",FCT,r);
      return *this;
   }
   if (isScalar()) {
      if (!cgr->cgd.SIZE.len && cgr->cstat!=CGD_ABELIAN) wblog(FL,
         "ERR %s() invalid CData\n%s",FCT,cgr->toStr().data);
   }
  
   RTD one=1, cfac=one;
   unsigned d=1, n=cgr->t.qlen(), k=(cgp.len ? cgp[0] : 0);
   QSet<TQ> Q(*this);

   if (!Q.isScalar() || cgw.len!=1) {
      MXPut(FL,"a").add(*this,"R").add(*cgr,"C").add(cgr->cgd,"c");
      wblog(FL,"ERR %s() unexpected scalar CRef\n%s",FCT,toStr().data);
   }

   if (!cgr->cgd.SIZE.len) {
      cfac=cgr->getScalar();
   }
   else {
      d=cgr->cgd.SIZE.el(k);
      if (cgr->cgd.isIdentity(cgp,&cfac)!=0) { 
         MXPut(FL,"a").add(*cgr,"cgr").add(cgr->cgd,"cgd");
         wblog(FL,"ERR %s() failed to reduce to Id (%g)",
         FCT,cgr->cgd.D.len? double(cgr->cgd.D[0]) : -1.);
      }
   }

   cfac*=cgw[0];

   QType t=cgr->t;
   qset<TQ> qs(n, cgr->qs.data+k*n);

   initIdentityR(FL,t,qs.data, d, xflag);
   cgw*=cfac;

   return *this;
};


template <class TQ, class TD>
wbvector<RTD>& CData<TQ,TD>::trace(
   const char *F, int L, wbvector<RTD> &cgt) const {

   if (cgd.isEmpty() && cstat==CGD_ABELIAN) {
      TD x=1; return cgt.init(1,&x);
   }

   unsigned r=rank();

   if (cgd.SIZE.len==r || (!cgd.SIZE.len && r==2)) {
      TD x=cgd.trace(); return cgt.init(1,&x);
   }

   if (cgd.SIZE.len!=r+1) wblog(F_L,"ERR %s() "
      "invalid OM data (%d/%d) !?",FCT,cgd.SIZE.len, r);
   return cgd.trace(r,cgt);
};


template <class TQ>
bool CRef<TQ>::isAbelian(const char *F, int L) const {

   if (cgr) {
      if (!cgr->isAbelian()) return 0;
      if (cgp.len && cgp.len!=cgr->qdir.len) wblog(FL,
      "ERR %s() invalid ablian CRef (cgp=[%s])",FCT,cgp.toStr().data);
   }
   else if (rtype!=CGR_ABELIAN) return 0;

   if (cgw.len>1 || (cgw.len && cgw[0]!=1)) wblog(F_L,
      "ERR %s() invalid ablian CRef (cgw=[%s])",FCT,cgw.toStr().data);
   return 1;
};


template <class TQ>
bool CRef<TQ>::isScalar() const {

   if (!cgr) {
      if (cgp.len || conj || cgw.len>(rtype==CGR_CTR_SCALAR ? 1:0))
         wblog(FL,"ERR %s() got invalid scalar (cgw.len=%d; %g)",
         FCT,cgw.len,cgw.len?double(cgw[0]):NAN);
      if (rtype!=CGR_CTR_SCALAR && 
          rtype!=CGR_CTR_ZERO && rtype!=CGR_ABELIAN)
         wblog(FL,"ERR %s() got unexpected rtype=%d",FCT,rtype.t);
      return 1;
   }
   else {
      unsigned m=cgr->getOM(FL);
      if (cgw.len>m) wblog(FL,
         "ERR %s() cgw out of bounds (%d/%d)",FCT,cgw.len,m);
      if (m>1)
           return 0;
      else return cgr->cgd.isScalar();
   }
};


template <class TQ>
unsigned CRef<TQ>::numel() const {

   if (!cgr) {
      if (rtype==CGR_ABELIAN || rtype==CGR_CTR_SCALAR) return 1;
      else {
         wblog(FL,"WRN %s() got rtype=%d",FCT,rtype.t);
         return 0;
      }
   }
   return cgr->cgd.numel();
};


template <class TQ> inline
unsigned CRef<TQ>::Size(unsigned k) const {

   if (!cgr && rtype==CGR_ABELIAN) return 1;
   if (cgr && cgr->cgd.isScalar()) return 1;

   unsigned r=cgr->qdir.len;
   if (!r) wblog(FL,"ERR %s() unexpected CRef",FCT);
  
   return Size(k,r);
};


template <class TQ>
unsigned CRef<TQ>::Size(unsigned k, unsigned r) const {

   if (!cgr || cgr->cgd.isScalar()) return 1;

   const cdata__ &c=cgr->cgd;
   const wbvector<unsigned> &S=c.SIZE;

   if (!S.len && c.isDiag()) {
      if (r!=2 || k>=r) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k+1,r);
      if (cgp.len && cgp.len!=2) wblog(FL,
         "ERR %s() invalid cgp.len=%d/2",FCT,cgp.len);
      return c.D.len;
   }

   if (k>=S.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,S.len);

   if (S.len<r || S.len>r+1) {
      if (S.len || r!=2 || !c.D.len) wblog(FL,
         "ERR %s() rank mismatch (%s; %d/%d/%d; d=%d)",
          FCT,c.sizeStr().data, k+1,r,S.len, c.D.len
      );
      return c.D.len;
   }

   if (k<r && cgp.len) {
      if (cgp.len!=r) wblog(FL,"ERR %s() "
         "invalid cgp=[%s; %d]",FCT,(cgp+1).toStr().data,r);
      k=cgp.data[k];
   }

   return S.data[k];
};


template <class TQ>
unsigned CRef<TQ>::getOM(const char *F, int L) const {
   if (cgr) {
      unsigned r=cgr->rank(), m=cgr->getOM();
      if (!cgw.len || cgw.len>m) wblog(F_L,
         "ERR %s() invalid cgw (%d/%d) !?",FCT,cgw.len,m);
      if (cgp.len && cgp.len!=r && cgp.len!=(r-1)) wblog(F_L,
         "ERR %s() invalid cgp (len=%d/%d) !?",FCT,cgp.len,r);
      return m;
   }
   else if (isAbelian()) { return 1; }
   else {
      wblog(F_L,"ERR %s() unknown OM [%s]",FCT,toStr().data);
      return 0;
   }
};


template <class TQ>
unsigned CRef<TQ>::rankS(char lflag) const {
   if (cgr) {
      unsigned r=cgr->rankS();
      if (!lflag) { unsigned m=cgr->getOM();
         if (!cgw.len || cgw.len>m) wblog(FL,
            "ERR %s() invalid cgw (%d/%d) !??",FCT,cgw.len,m);
         if (cgp.len && cgp.len!=r) wblog(FL,
            "ERR %s() invalid cgp (len=%d/%d) !??",FCT,cgp.len,r
         );
      }
      return r;
   }
   else {
      if (!lflag) {
         if (cgw.len || cgp.len) wblog(FL,"ERR %s() "
         "cgw or cgp out of bounds (%d,%d/0)",FCT,cgw.len,cgp.len);
      }
      return 0;
   }
};


template <class TQ>
unsigned CRef<TQ>::rank(const char *F, int L, char lflag) const {
   if (cgr) {
      unsigned r=cgr->rank();
      if (!lflag) { unsigned m=cgr->getOM();
         if (!cgw.len || cgw.len>m) wblog(F_L,
            "ERR %s() invalid cgw (%d/%d) !??",FCT,cgw.len,m);
         if (cgp.len && cgp.len!=r && cgp.len!=(r-1)) wblog(F_L,
            "ERR %s() invalid cgp (len=%d/%d) !??",FCT,cgp.len,r
         );
      }
      return r;
   }
   else {
      wblog(F_L,"ERR %s() unknown rank (since cgref=NULL)",FCT);
      return 0;
   }
};


template <class TQ>
template <class T>
wbvector<T>& CRef<TQ>::getSize(wbvector<T> &S, char bare) const {

   if (cgr) {
      cgr->getSize(S, bare, cgp.len ? &cgp : NULL);
   }
   else {
      if (cgp.len || cgw.len) wblog(FL,"ERR %s() "
         "invalid scalar/empty CRef (%d,%d)",FCT,cgw.len,cgp.len);
      S.init();
   }

   return S;
};


template <class TQ>
int CRef<TQ>::HConjOpScalar() {

   if (!cgr || cgr->isEmpty()) {
      if (cgp.len) { wblog(FL,
         "WRN %s() got cgp.len=%d for scalar !?",FCT,cgp.len);
         cgp.init();
      }
      if (conj) { wblog(FL,
         "WRN %s() got conj=%d for scalar !?",FCT,conj);
         conj=0;
      }
      return 0;
   }

   unsigned i=0, r=cgr->qdir.len, n=(r ? cgr->qs.len/r : 0);
   const TQ *q=cgr->qs.data;

   if (cgp.len && cgp.len!=r) wblog(FL,
      "ERR %s() invalid permutation cgp.len=%d/%d !?",FCT,cgp.len,r);
   if (r<2 || cgr->qs.len%r || !n) wblog(FL,
      "ERR %s() got rank-%d QSpace (qlen=%d)",FCT,r,n);

   if (cgp.len || conj) {
      QDir qd(cgr->qdir);
         if (cgp.len) qd.Permute(cgp);
         if (conj) qd.Conj();
         if (r==3) qd.Conj(2);
      if (qd!=cgr->qdir) return 2;
      else { cgp.init(); conj=0; }
   }

   if (r!=3) {
      if (r!=2) {
         wblog(FL,"WRN %s() got rank-%d CGR data !?",FCT,r);
         return 1;
      }
      return 0;
   }

   if (memcmp(q,q+n,n)) {
      return 12;
   }

   for (q+=2*n; i<n; ++i) {
      if (q[i]!=0) return (30+i);
   }

   return 0;
};


template <class TQ>
int CRef<TQ>::checkAbelian(const char *F, int L) const {

   if (cgr) {
      if (!cgr->isAbelian()) {
         if (F) wblog(F,L,"ERR %s() "
            "got cgref with Abelian symmetry (%s)",FCT,cgr->toStr().data);
         return 1;
      }
      if (cgp.len && cgp.len!=cgr->qdir.len) wblog(FL,
         "ERR %s() got invalid abelian CRef\n%s",FCT,toStr().data);
   }

   int e=0;

   if (cgp.len && cgp.len!=cgr->qdir.len) e=2; else
   if (!cgw.len) { if (rtype!=CGR_ABELIAN) e=3; }
   else {
      if (cgw.len>1) e=4; else
      if (cgw[0]!=RTD(1)) e=5;
   }

   if (e && F) wblog(F,L,
      "ERR %s() invalid scalar (%d)\n%s",FCT,e,toStr().data);
   return e;
};


template <class TQ>
int CRef<TQ>::check(const char *F, int L) const {

   int e=0;

   if (!cgr) {
      if (rtype==CGR_ABELIAN) {
         if (cgp.len || cgw.len) e=1;
      }
      else if (rtype==CGR_CTR_SCALAR) {
         if (cgp.len || cgw.len!=1) e=2;
      }
      else if (cgp.len || cgw.len) e=3;

      if (e && F) wblog(F_L,"ERR %s() "
         "invalid scalar CRef (e=%d)\n%s",FCT,e,toStr().data);
      return e;
   }
   else {
      if (cgw.len<1) e=11; else
      if (cgw.len>cgr->getOM()) e=12;

      if (e && F) wblog(F_L,
         "ERR %s() invalid CRef data (e=%d, m=%d)\n%s",
         FCT,e,cgr->getOM(), toStr().data);
      return e;
   }

   return 0;
};


template <class TQ>
int CRef<TQ>::check(const char *F, int L, 
   const QType &q, const QDir &qdir
 ) const {

   int e=check(F,L); if (e) return e;

   if (q.isAbelian()) {
      if (cgr)  {
         if (sameQDir(FL,qdir)) return 1;
      }
      if (cgw.len) {
         if (cgw.len!=1 || cgw[0]!=1) { if (F) wblog(FL,"ERR %s() "
            "got cgw=[%s] for %s",FCT,cgw.toStr().data, q.toStr().data);
            return 3;
         }
      }
   }
   else {
      if (!cgr) { 
         if (rtype==CGR_CTR_SCALAR && cgw.len==1 && cgw[0]==1) {
            return 0;
         }
         if (qdir.len) {
            if (F) wblog(F,L,"ERR %s() got NULL cgref for %s: %s (%s)",
               FCT,q.toStr().data,qdir.toStr().data,cgw.toStr().data);
            return 12;
         }
         return 11;
      }
      if (sameQDir(FL,qdir)) {
         return 12;
      }
   }

   return 0;
};


template <class TQ>
bool CRef<TQ>::isDiagCGC(RTD eps) const {

   if (isScalar()) return 1;

   if (!isRefInit(FL)) {
      const cdata__ &c=cgr->cgd;
      if (c.numel()>1 && !c.isDiagMatrix(eps)) return 0;
   }
   else {
      unsigned m=getOM(FL);
      if (cgr->cgd.D.len!=m) wblog(FL,"ERR %s() "
         "got invalid CGD_REF data (%d/%d)",FCT,cgr->cgd.D.len,m);
      return (cgr->cgd.isSMatrix());
   }
   return 1;
};


template <class TQ>
bool CRef<TQ>::isIdentityCGC(RTD *nrm, RTD eps) const {

   if (isScalar()) {
      if (nrm) (*nrm)=RTD(1);
      return 1;
   }

   if (!cgr) wblog(FL,"ERR %s() got null cgr !?",FCT);
   if (!cgr->isSymmmetric()) return 0;

   unsigned d=1, r=rank(FL);

   if (!r || r!=cgr->qdir.len) wblog(FL,
      "ERR %s() got r=%d/%d !?",FCT,r,cgr->qdir.len);
   if (r%2 || (r<3 && cgw.len!=1) || !cgw.len) return 0;

   wbindex ip;
   cgr->qdir.findGT(0,ip);
   
   if (!ip.len || 2*ip.len!=cgr->qdir.len) return 0;
   if (ip.len==1 && cgr->cgd.isDiag()) d=cgr->cgd.dim();
   else {
      const SPIDX_T *s=cgr->cgd.SIZE.data;
      for (unsigned i=0; i<ip.len; ++i) d*=s[ip.data[i]];
   }

   if (cgr->cstat!=CGD_REF_INIT) {

      unsigned i=0; RTD x=0;
      wbvector<RTD> cgt; trace(FL,cgt);

      if (cgr->cstat==CGD_REF_INIT) wblog(FL,
         "ERR %s() got mixed ref_init\n%s",FCT,toStr().data);
      if (cgt.len<cgw.len || !cgt.len) wblog(FL,
         "ERR %s() got OM mismatch (%d/%d)",FCT,cgt.len,cgw.len);
      for (; i<cgw.len; ++i) {
          if (ABS(cgt[i])<eps) wblog(FL,
             "ERR %s() got small Id fac (%g)",FCT,double(cgt[i]));
          x+=cgw[i]/cgt[i];
      }
      if (nrm) { (*nrm)=x; }

      if (cgt.len==1) { x=RTD(1)/cgt[0];
         return cgr->cgd.isProptoId(x,eps);
      }
      else {
         return 1;
      }
   }
   else {
      if (!isRefInit() || !cgr->cgd.SIZE.len) wblog(FL,
         "ERR %s() got invalid CGD_REF data\n%s",FCT,toStr().data);
      const wbvector<RTD> &D=cgr->cgd.D;

      if (nrm) {
         if (D.len) {
            unsigned i=0; RTD x=0;
            if (D.len!=cgw.len || !D.len) wblog(FL,
               "ERR %s() got OM mismatch (%d/%d)",FCT,D.len,cgw.len);

            for (; i<D.len; ++i) {
               if (ABS(D[i])<eps) wblog(FL,
                  "ERR %s() got small Id fac (%g)",FCT,double(D[i]));
               x+=(cgw[i]/D[i]);
            }
            (*nrm)=x;
         }
         else (*nrm)=RTD(1);
      }
   }

   return 1;
};


template <class TQ>
RTD CRef<TQ>::NormStd(const char *F, int L, unsigned r)  {

   RTD nrm=1;
   if (isScalar()) return nrm;

   if (!cgr) wblog(F_L,"ERR %s() got null cgr !?",FCT);

   wbindex ip,in;
   unsigned i;

   if (int(r)>0 && r!=cgr->qdir.len) wblog(F_L,"ERR %s() "
      "got unexpected CRef of rank %d/%d",FCT,r,cgr->qdir.len);
   else {
      if (r!=cgr->qdir.len)
      wblog(F_L,"ERR %s() got r=%d/%d !?",FCT,r,cgr->qdir.len);
   }

   cgr->qdir.findGT(0,ip,in);

   if (in.len==1) i=in[0]; else
   if (ip.len==1) i=ip[0]; else {
      if (F) wblog(F_L,"ERR %s() "
         "got unexpected CRef %s",FCT,toStr().data);
      return nrm;
   }

   i=cgr->cgd.dim0(i);
   if (i!=1) { 
      if (int(i)<=0) wblog(F_L,"ERR %s() invalid dim=%d !?",FCT,i);
      nrm=SQRT(RTD(i)); cgw*=nrm;
      nrm=RTD(1)/nrm;
   }

   return nrm;
};


template <class TQ>
RTD CRef<TQ>::normExt(const char *F, int L) const {

   if (isScalar() ) return 1;
   if (rank(F_L)>2) return 1;


   if (!cgr) wblog(FL,"ERR %s() got null cgr !?",FCT);

   unsigned r=rank(F_L), d=cgr->cgd.maxDim(r);
   if (r!=cgr->qdir.len) wblog(FL,
      "ERR %s() got rank/OM mismatch (%d/%d)",FCT,r,cgr->qdir.len);

   if (int(d)<=0) wblog(FL,"ERR %s() invalid dim=%d !?",FCT,d);

   return SQRT(RTD(d));
};


template <class TQ>
RTD CRef<TQ>::norm2(char xflag) const {

   RTD w2=1;
   if (!cgr && isAbelian()) return w2;
   if (!cgw.len) wblog(FL,
      "ERR %s() got cgw.len=0 !?\n%s",FCT,toStr().data);
   w2=cgw.norm2();

   if (isEmpty()) {
      if (w2!=0) wblog(FL,
         "ERR %s() got |cgw|^2=%g for empty CRef !?",FCT,double(w2));
      return 0;
   }

   if (cgr) {
      if (xflag && !isRefInit()) { cgr->checkNormSign(FL,xflag); }
      else {
         unsigned r=cgr->rank(FL); check(FL);
         if (r<2) wblog(FL,"ERR %s() got r=%d",FCT,r);
      }
   }

   return w2;
};


template <class TQ> inline
bool CRef<TQ>::isSortedDQ(QSet<TQ> *QS, wbperm &pxt, char &cxt) const {

   if (!cgr) {
      if (pxt.len) pxt.init(); cxt=0;
      return 1;
   }

   QSet<TQ> Q(*cgr); 
      adapt(Q,'i'); if (QS) QS->init(Q);
      Q.Sort(&pxt, &cxt,'i');

   if (cgp.sameAs(pxt) && conj==cxt) {
      return 1;
   }
   return 0;
};


template <class TQ>
int CRef<TQ>::SortDegQ(const char *F, int L, QSet<TQ> *QS) {

   wbperm pxt; char cxt=0;
   if (isSortedDQ(QS,pxt,cxt)) return 0;

   if (cgr) {
      if (cgr->t.isAbelian()) {
         if (cgr->cgd.D.len && (cgr->cgd.D.len>1 || cgr->cgd.D[0]!=1))
            wblog(F_L,"ERR %s() invalid abelian CRef data\n%s",FCT,
            toStr().data
         );
      }
      else if (cgr->qdir.len<3) {
         if (cgr->qdir.len!=2) wblog(FL,"ERR %s()",FCT,toStr().data);
         if (isRefInit()) {
            const wbvector<unsigned> &S=cgr->cgd.SIZE;
            if (S.len!=2 || S[0]!=S[1] || cgr->cgd.D.len>1) wblog(F_L,
               "ERR %s() invalid ref CData (%d,%d)\n%s",
               FCT, cgr->cgd.SIZE.len, cgr->cgd.D.len, toStr().data
            );
         }
         else {
            if (cgr->cgd.SIZE.len || !cgr->cgd.D.len || cgw.len!=1)
               wblog(F_L,"ERR %s() invalid scalar CData (%d,%d)\n%s",
               FCT, cgr->cgd.SIZE.len, cgr->cgd.D.len, toStr().data
            );
         }
      }
      else {
         CRef<TQ> R;
         R=*this; pxt.save2(R.cgp); R.conj=cxt;

         gCX.contractDQ(FL,*this,R);

         R.save2(*this);
         return 2;
      }
   }

   pxt.save2(cgp); conj=cxt;
   return 1;
};


template <class TQ>
wbvector<RTD>& CRef<TQ>::trace(
   const char *F, int L, wbvector<RTD> &cgt) const {

   if (!cgr || !cgr->isSymmmetric()) { 
      if (F) wblog(F_L,"ERR %s() got %s",FCT,toStr().data);
      if (cgt.len) cgt.init();
      return cgt;
   }
   return cgr->trace(F_L,cgt);
};


template <class TQ>
RTD CRef<TQ>::trace() const {

   if (!cgr) { 
      if (cgw.len) wblog(FL,
         "ERR %s() invalid cgw.len=%d (cgr=0)",FCT,cgw.len);
      return 1;
   }

   if (!cgr->isSymmmetric()) wblog(FL,
      "ERR %s() got %s",FCT,toStr().data);

   wbvector<RTD> cgt;
   unsigned i=0, m=cgr->gotOM(FL);
   RTD x=0;

   if (m) wblog(FL,
      "WRN %s() in the presence of OM=%d !?",FCT,m);

   cgr->trace(FL,cgt);

   if (cgw.len>m || cgt.len!=m) wblog(FL,"ERR %s() "
      "got OM size mismatch (%d,%d/%d)",FCT,cgw.len,cgt.len,m);

   for (; i<cgw.len; ++i) { x+=(cgt[i]*cgw[i]); }

   return x;
};


template <class TQ>
CRef<TQ>& CRef<TQ>::Permute(const wbperm &P, char iflag, char isnew) {

   if (!P.len) return *this;

   if (!cgr) {
      if (cgw.len || cgp.len) wblog(FL,
         "ERR %s() invalid CRef\n%s",FCT,toStr().data);
      return *this;
   }

   unsigned r=cgr->rank(FL);

   if (isnew && cgp.len) wblog(FL,"ERR %s() got existing "
      "permutation cgp=[%s] !??",FCT,cgp.toStr().data);

   if (cgr->cgd.isScalar()) {
      if (cgw.len>1) wblog(FL,"ERR %s() "
         "got cgw.len=%d having %s !??",FCT,cgw.len,cgr->QStr().data);
   }

   cgp.Permute(P,iflag,r);

   if (cgp.isIdentityPerm()) {
      cgp.init();
   }
   else SortDegQ(FL);

   return *this;
};


template <class TQ>
char CRef<TQ>::sameSize(const CRef &B, char strict) {

   if ((cgr==NULL) ^ (B.cgr==NULL)) wblog(FL,
      "ERR CRef::%s() inconsistency in cgr",FCT);
   if (cgr) {
      if (!cgw.len || !B.cgw.len) wblog(FL,"ERR %s() "
         "got empty cgw data (%d/%d)",FCT,cgw.len,B.cgw.len);
      return cgr->sameSizeR(*B.cgr,strict);
   }
   else {
      if (cgw.len!=1 || B.cgw.len!=1) wblog(FL,"ERR %s() "
         "got invalid scalar cgw data (%d/%d)",FCT,cgw.len,B.cgw.len);
      return (cgw.len==B.cgw.len);
   }
};


template <class TQ>
RTD CRef<TQ>::NormSignR(
   const char *F, int L, RTD eps, RTD eps2
){
   if (!cgr) { checkAbelian(F_L); return 1; }

   unsigned i, m=cgr->getOM(F_L);
   RTD a,
      cfac=cgw.norm()/normExt(F_L);

   if (cgw.len>m) wblog(FL,
      "ERR %s() cgw out of bounds (%d/%d)",FCT,cgw.len,m);
   if (cfac<=eps) wblog(FL,
      "ERR %s() got small cfac (%.3g) !??",FCT,double(cfac));

   for (i=0; i<cgw.len; ++i) { if ((a=ABS(cgw[i]))>eps) {
      if (cgw[i]<0) { cfac=-cfac; }
      break;
   }}

   for (i=0; i<cgw.len; ++i) { a=ABS(cgw[i]);
      if (a<eps) {
         if (a>eps2) wblog(FL,
            "WRN %s() CGW noise %.3g (%g, %g; i=%d/%d)",
            FCT, double(cgw[i]), double(eps), double(eps2),i,m);
         cgw[i]=0;
      }
   }

   if (cfac!=1) { cgw *= (RTD(1)/cfac); }

   return cfac;
};


template <class TQ>
char CRef<TQ>::sameUptoFac(const CRef &B, RTD *fac_, RTD eps) const {

   if ((cgr==NULL) ^ (B.cgr==NULL)) wblog(FL,
      "ERR CRef::%s() inconsistency in cgr",FCT);
   if (cgr!=B.cgr) return 1;

   if (!cgr) {
      if (cgw.len || B.cgw.len) wblog(FL,"ERR %s() "
         "invalid CRef data (%d/%d; cgr=0)",FCT,cgw.len,B.cgw.len);
      if (fac_) { (*fac_)=1; }
      return 0;
   }

   if (!cgw.len && !B.cgw.len) {
      if (fac_) (*fac_)=1;
      return 0;
   }

   if (anyTrafo() || B.anyTrafo()) {
      if (!sameQSet(B)) return 2;
   }

   unsigned i=0, n1=cgw.len, n2=B.cgw.len,
      n=(n1<n2? n1:n2), m=cgr->getOM();
   const RTD *x1=cgw.data, *x2=B.cgw.data;
   RTD a=cgw.aMax(i);

   if (!cgw.len || cgw.len>m || !B.cgw.len || B.cgw.len>m) wblog(FL,
      "ERR %s() cgref out of bounds (%d,%d/%d)",FCT,cgw.len,B.cgw.len,m);
   if (a<=eps) wblog(FL,
      "ERR %s() got small CData (%.3g) !??",FCT,double(cgw.norm()));
   if (i>=n2 || ABS(x2[i])<=eps) return 3;

   RTD fac=x2[0]/x1[0]; if (fac_) (*fac_)=fac;

   for (i=1; i<n; ++i) {
      a=ABS(fac*x1[i]-x2[i]); if (a>eps) return 10;
   }

   if (i<n1) { for (; i<n1; ++i) { if (ABS(x1[i])>eps) return 11; }} else
   if (i<n2) { for (; i<n2; ++i) { if (ABS(x2[i])>eps) return 12; }}

   return 0;
};


template <class TQ>
char CRef<TQ>::sameAs(const CRef &B, RTD eps) const {

   if ((cgr==NULL) ^ (B.cgr==NULL)) wblog(FL,
      "ERR CRef::%s() inconsistency in cgr",FCT);
   if (cgr!=B.cgr) return 1;

   if (!cgw.len && !B.cgw.len) return 0;
   if (!cgr) {
      if (cgw.len || B.cgw.len || rtype!=CGR_ABELIAN || rtype!=B.rtype)
         wblog(FL,"ERR %s() invalid CRef data (%d/%d; %d/%d)",FCT,
         cgw.len,B.cgw.len, rtype.t, B.rtype.t);
      return 0;
   }

   unsigned i=0, n1=cgw.len, n2=B.cgw.len,
      n=(n1<n2? n1:n2), m=cgr->getOM();
   const RTD *x1=cgw.data, *x2=B.cgw.data;
   RTD a=cgw.aMax(i);

   if (!cgw.len || cgw.len>m || !B.cgw.len || B.cgw.len>m) wblog(FL,
      "ERR %s() cgref out of bounds (%d,%d/%d)",FCT,cgw.len,B.cgw.len,m);
   if (a<=eps) wblog(FL,
      "ERR %s() got small CData (%.3g) !??",FCT,double(cgw.norm()));
   if (i>=n2 || ABS(x2[i])<=eps) return 3;

   for (i=1; i<n; ++i) { if (ABS(x1[i]-x2[i])>eps) return 10; }

   if (i<n1) { for (; i<n1; ++i) { if (ABS(x1[i])>eps) return 11; }} else
   if (i<n2) { for (; i<n2; ++i) { if (ABS(x2[i])>eps) return 12; }}

   return 0;
};


template <class TQ>
char CRef<TQ>::got3(const char *F, int L,
   const QType &q, const qset<TQ> &J12, const qset<TQ> &J3
 ) const {

   if (!cgr) wblog(F_L,"ERR %s() got empty cgr",FCT);

   const CDATA_TQ &C = *cgr;

   if (q!=C.t) {
      if (F) wblog(FL,"ERR %s() symmetry mismatch (%s <> %s)",
         FCT,q.toStr().data,C.t.toStr().data);
      return 1;
   }
   if (C.qs.len!=J12.len + J3.len || J12.len!=2*J3.len) wblog(FL,
      "ERR %s() qset inconsistency (%d == 2*%d; %d)",FCT,
      J12.len,J3.len,C.qs.len);

   if (cgp.len) {
      if (cgp.len!=C.qdir.len) wblog(FL,"ERR %s() "
         "invalid cgp=[%s] %d",FCT,cgp.toStr().data,C.qdir.len);

      unsigned i=0, l=J12.len;
      qset<TQ> qx; C.qs.blockPermute(cgp,qx);
      const TQ *qs=qx.data;

      for (; i<J12.len; ++i) {
         if (qs[i]!=J12.data[i]) {
            if (F) wblog(FL,
               "ERR %s() qset mismatch: %s @ %s <> (%s | %s); i=%d",
               FCT, qx.toStr().data, (cgp+1).toStr().data,
               J12.toStr().data, J3.toStr().data, i+1);
            return 2;
         }
      }
      for (; i<J3.len; ++i) {
         if (qs[i+l]!=J3.data[i]) {
            if (F) wblog(FL,
               "ERR %s() qset mismatch (%s @ %s <> %s ; %s; %d)",
               FCT, qx.toStr().data, (cgp+1).toStr().data,
               J12.toStr().data, J3.toStr().data, i+l+1);
            return 3;
         }
      }
   }
   else {
      const TQ *qs=C.qs.data;
      if (Wb::cmpRange(qs, J12.data, J12.len)) return 2; qs+=J12.len;
      if (Wb::cmpRange(qs, J3 .data, J3 .len)) return 3;
   }

   return 0;
};


template <class TQ>
bool CRef<TQ>::sameQDir(const char *F, int L, const QDir &qd) const {

   if (!cgr) return 1;

   if (cgr->qdir.len!=qd.len) wblog(F_L,"ERR %s() QDir length mismatch "
      "(%s <> %s)",FCT, cgr->qdir.toStr().data, qd.toStr().data);
   if (cgp.len && cgp.len!=qd.len) wblog(F_L,
      "ERR %s() got length mismatch (%d/%d/%d)",
      FCT,cgr->qdir.len, qd.len, cgp.len
   );

   char e=0, cflag=gotConj(), pflag=gotPerm();
   if (!cflag && !pflag) return qd==cgr->qdir;

   unsigned i=0;
   const char *d=cgr->qdir.data, *d0=qd.data;

   if (pflag) {
      const PERM_T *p=cgp.data;

      if (cflag) {
         for (; i<qd.len; ++i) { if (d[p[i]]!=-d0[i]) { e=1; break; }}
      }
      else {
         for (; i<qd.len; ++i) { if (d[p[i]]!= d0[i]) { e=2; break; }}
      }
   }
   else {
      for (; i<qd.len; ++i) { if (d[i]!=-d0[i]) { e=3; break; }}
   }

   if (e) { sprintf(str,
      "%s <> %s (e=%d)", toStr().data, qd.toStr().data, e);
      if (F) wblog(F,L,"ERR %s() got\n%s",FCT,str);
      else   wblog(FL, "WRN %s() got\n%s",FCT,str); 
   }

   return (e? 0:1);
};


template <class TQ>
bool CRef<TQ>::sameQSet(const CRef &B) const {

   bool ipa=gotPerm(FL), ipb=B.gotPerm(FL);

   if (!cgr || cgr!=B.cgr) return 1;
   if (!ipa && !ipb) { return (conj!=B.conj); }

   const QSet<TQ> &Q=(QSet<TQ>&)(*cgr);
   unsigned i=0, r=Q.qdir.len;

   if ((cgp.len && cgp.len!=r) || (B.cgp.len && B.cgp.len!=r)) {
      wblog(FL,"ERR %s() got cgp.len=%d,%d/%d !?",
      FCT,cgp.len,B.cgp.len,r); 
   }

   if (Q.qs.blockCompare(cgp, Q.qs, B.cgp)) return 0;

   const PERM_T *pa=cgp.data, *pb=B.cgp.data;
   const char *d=Q.qdir.data;

   if (ipa) {
      if (ipb) {
         if (conj ^ B.conj) 
              for (; i<r; ++i) { if (d[pa[i]]!=-d[pb[i]]) return 0; }
         else for (; i<r; ++i) { if (d[pa[i]]!=+d[pb[i]]) return 0; }
      }
      else {
         if (conj ^ B.conj) 
              for (; i<r; ++i) { if (d[pa[i]]!=-d[i]) return 0; }
         else for (; i<r; ++i) { if (d[pa[i]]!=+d[i]) return 0; }
      }
   }
   else {
      if (ipb) {
         if (conj ^ B.conj) 
              for (; i<r; ++i) { if (d[i]!=-d[pb[i]]) return 0; }
         else for (; i<r; ++i) { if (d[i]!=+d[pb[i]]) return 0; }
      }
   }

   return 1;
}; 


template <class TQ>
int CRef<TQ>::gotSameCData(
   const char *F, int L, const CRef& B, char xflag
 ) const {

   if (cgr!=B.cgr) {
      if (F) wblog(F_L,"ERR %s() got CGC mismatch\n%s <> %s",FCT,
            cgr ?   cgr->toStr().data : "(null)",
          B.cgr ? B.cgr->toStr().data : "(null)"
      );
      return -1;
   }

   if (!cgr) {
      if (cgw.len || B.cgw.len || rtype!=CGR_ABELIAN || rtype!=B.rtype)
         wblog(FL,"ERR %s() invalid scalar cgref (%d/%d; %d/%d)",
         FCT, cgw.len,B.cgw.len, rtype.t, B.rtype.t);
      return 0;
   }

   unsigned r=cgr->rank(FL), n=MIN(cgw.len, B.cgw.len), m=cgr->getOM();

   if (r<2) { if (F) wblog(FL,
      "ERR %s() got rank-%d CData (%s)",FCT,r,cgr->sizeStr().data);
      return -2;
   }

   if (cgw.len>m || B.cgw.len>m) { if (F) wblog(FL,
      "ERR %s() CData out of multiplicity range (%d,%d/%d)",
      FCT,cgw.len,B.cgw.len,m);
      return -3;
   }

   if (!cgp.sameAs(B.cgp) || conj!=B.conj) {
      if (F) wblog(FL,
         "ERR %s() got cgp+conj mismatch (%s%s <> %s%s)", FCT,
           cgp.toStr(-1,"").data,   conj?"*":"",
         B.cgp.toStr(-1,"").data, B.conj?"*":"");
      return -4;
   }

   if (xflag) {
      char ira=isRefInit(), irb=B.isRefInit();
      if (ira || irb) {
         if (ira ^ irb) wblog(FL,
            "ERR %s() got ref mismatch (%d/%d)",FCT,ira,irb);
      }
      else { cgr->checkNormSign(F_L,xflag); }
   }

   return n;
};


template <class TQ>
bool CRef<TQ>::sameQSet(const QSet<TQ> &Q) const {

   if (!cgr || cgr->isEmpty()) {
      if (cgp.len || conj) wblog(FL,
         "WRN %s() invalid scalar (0x%lX)\ncgp.len=%d, conj=%d : %s",
         FCT, cgr, cgp.len, conj, Q.toStr().data);
      return Q.isEmpty();
   }

   if (anyTrafo()) {
      QSet<TQ> X((QSet<TQ>&)(*cgr));
         if (cgp.len) X.Permute(cgp);
         if (conj) X.Conj();
      return X==Q;
   }

   return ((QSet<TQ>&)(*cgr))==Q;
};


template <class TQ>
bool CRef<TQ>::affectsQDir() const {

   if (!cgr || !cgr->qdir.len || !cgp.len) return 0;
   if (cgp.len && cgp.len!=cgr->qdir.len) wblog(FL,
      "ERR %s() got length mismatch (%d/%d)",FCT,cgp.len,cgr->qdir.len
   );

   char cflag=gotConj(), pflag=gotPerm();
   if (!cflag && !pflag) return 0;

   if (pflag) {
      const char *d=cgr->qdir.data;
      const unsigned *p=cgp.data; unsigned i=0;

      if (cflag) {
         for (; i<cgp.len; ++i) { if (d[i]!=-d[p[i]]) return 1; }
      }
      else {
         for (; i<cgp.len; ++i) { if (d[i]!= d[p[i]]) return 1; }
      }
   }
   else {
      return 1;
   }

   return 0;
};


template <class TQ>
RTD CRef<TQ>::scalarProd(const char *F, int L, const CRef &B) const {

   RTD x2=0;

   unsigned i=0, n=gotSameCData(F_L,B,'!');
   if (int(n)<=0) {
      if (!n) return x2;
      else wblog(FL,"ERR %s() got e=%d",FCT,n);
   }

   for (; i<n; ++i) { x2 += cgw[i]*B.cgw[i]; }

   return x2;
};


template <class TQ>
RTD CRef<TQ>::normDiff2(const char *F, int L, const CRef &B) const {

   RTD x2=0;

   unsigned i=0, n=gotSameCData(F_L,B,'!');

   if (int(n)<=0) {
      if (!n) return x2;
      else wblog(FL,"ERR %s() got e=%d",FCT,n);
   }

   if (n) { RTD x;
      for (; i<n; ++i) { x=cgw[i]-B.cgw[i]; x2+=(CONJ(x)*x); }
   }

   if (cgw.len>n) {
      const RTD *x=cgw.data; n=cgw.len;
      for (; i<n; ++i) { x2+=(CONJ(x[i])*x[i]); }
   } else {
      const RTD *x=B.cgw.data; n=B.cgw.len;
      for (; i<n; ++i) { x2+=(CONJ(x)[i]*x[i]); }
   }

   return x2;
};


template <class TQ>
RTD CRef<TQ>::sumTimesEl(const char *F, int L, const CRef &B) const {

   if (!cgr && !B.cgr) {
      if (cgw.len || B.cgw.len || rtype!=CGR_ABELIAN || rtype!=B.rtype)
         wblog(FL,"ERR %s() invalid scalar cgref (%d/%d)",
         FCT,cgw.len,B.cgw.len, rtype.t, B.rtype.t);
      return RTD(1);
   }

   RTD x2=0;

   unsigned i=0, n=gotSameCData(F_L,B,'!');

   if (int(n)<=0) wblog(F_L,"ERR %s() got e=%d",FCT,n);

   for (; i<n; ++i) { x2+=cgw[i]*B.cgw[i]; }

   return x2;
};


template <class TQ>
RTD CRef<TQ>::safeCpy(
   const char *F, int L, const CRef<TQ> &B, CRef<TQ> &C
 ) const {

   RTD x=1; unsigned e=0;

   if (this->isEmpty() && B.isEmpty()) {
      if (!C.isEmpty()) wblog(F_L,"ERR CGdata inconsistency");
      return x;
   }

   if (!(e=sameUptoFac(B,&x))) { C=(*this); }
   else {
#ifdef MATLAB_MEX_FILE
      MXPut(FL).addP(this->toMX(),"a")
      .addP(B.toMX(),"b").add(C,"c").add(x,"x").add(e,"e");
#endif
      wblog(F_L,"ERR %s() CG data not the same (e=%d)",FCT,e);
   }

   return x;
};


template <class TQ>
wbstring CRef<TQ>::toStr() const {

   if (!cgr) {
      unsigned l=0, n=32; char s[n];

      if (cgw.len>1) wblog(FL,"ERR %s() got cgw.len=%d for %s (cgr=0)",
         FCT, cgw.len, rtype.tostr());
      if (!cgw.len) {
         if (rtype==CGR_ABELIAN)
              l=snprintf(s,n,"(abelian @ [1])");
         else l=snprintf(s,n,"(%s @ [])", rtype.tostr());
      }
      else {
         l=snprintf(s,n,"(%s @ [%s])",
         rtype==CGR_CTR_SCALAR ? "contracted": rtype.tostr(),
         cgw.toStr().data);
      }

      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
      return s;
   }
   else {
      wbstring s(118+10*cgw.len); unsigned l=0, n=s.len;

      if (rtype>=CGR_NUM_TYPES) wblog(FL,"ERR %s() "
         "cgrType out of range (%d/%d)",FCT,rtype.t,CGR_NUM_TYPES
      );

      if (!cgr) l=snprintf(s.data,n,"(null)");
      else {
         QSet<TQ> Q(*cgr); adapt(Q,'i');
         l=snprintf(s.data,n,"%s", Q.toStr(statStr().data).data);
      }

      if (l<n) {
         if (cgw.len==1)
              l+=snprintf(s.data+l,n-l,", w=%g",double(cgw[0]));
         else l+=snprintf(s.data+l,n-l,", w=[%s]",cgw.toStr().data);
      }

      if ((cgp.len || conj) && l<n) { wbperm P(cgp,'i');
         l+=snprintf(s.data+l,n-l,", p=%s%s",
         (P+1).toStrf("%d","").data, conj?"*":"");
      }

      if (l>=n) wblog(FL,"ERR %s() "
         "string out of bounds (%d/%d)%N%N%s%N%N",FCT,l,n,s.data);

      return s;
   }
};


template <class TQ>
wbstring CRef<TQ>::sizeStr() const {

   if (cgr) { return cgr->sizeStr(); }

   if (cgw.len) wblog(FL,
      "ERR %s() inconsistent cgw=[%s]",FCT,cgw.toStr().data);
   return (rtype==CGR_ABELIAN ? "(scalar)" : "[]");
};


template <class TQ>
wbstring CRef<TQ>::statStr(char vflag) const { 
   wbstring s(32);
   unsigned l=0, n=s.len;

   if (vflag) {
      l=snprintf(s.data,n,"%s",rtype.tostr());
      if ((cgp.len || conj) && l<n) l+=snprintf(s.data+l,n-l,
         " [%s]%s",(cgp+1).toStr().data, conj?"*":"");
      if (cgr && l<n) l+=snprintf(s.data+l,n-l," 0x%lX",(long)cgr);
   }

   if (l<n) {
   if (cgr && cgr->cstat.t) {
      if (rtype.t) {
         l+=snprintf(s.data+l,n-l,"%s;%s",
         rtype.tostr(), cgr->cstat.toStr().data);
      } else {
         l+=snprintf(s.data+l,n-l,"%s",cgr->cstat.toStr().data);
      }
   }
   else {
      if (rtype.t) {
         l+=snprintf(s.data+l,n-l,"%s",rtype.tostr());
      }
   }}

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};


template <class TQ>
mxArray* CRef<TQ>::toMx(char flag) const {

   mxArray* a=mxCreateStruct(1,1,flag);
   this->add2MxStruct(a,0,flag);
   return a;
};


template <class TQ>
mxArray* CRef<TQ>::mxCreateStruct(
   unsigned m, unsigned n, char flag) const {

   if (m && n) {
      if (flag==0) {
         const char *fields[]={
           "type","qset","qdir","cid","size","nnz","cgw","cgt"};"cgp","conj"
         return mxCreateStructMatrix(m,n,8,fields);
      }
      else {
         const char *fields[]={
           "type","qset","qdir","cgs","size","nnz","cgw"};"cgp","conj"
         return mxCreateStructMatrix(m,n,7,fields);
      }
   }
   else return NULL;
};


template <class TQ>
void CRef<TQ>::add2MxStruct(mxArray *S, unsigned k, char flag) const {

   if (cgr) {
      wbvector<size_t> S2; getSize(S2);
      if (!S2.len) wblog(FL,"WRN %s() got empty S2 !?",FCT);

      if (!cgr->QSet<TQ>::isSorted()) wblog(FL,
         "ERR %s() got non-sorted CData !?\n%s",FCT,cgr->toStr().data);

      if (anyTrafo()) {
         QSet<TQ> QS;
         if (!isSortedDQ(&QS)) wblog(FL,
            "ERR %s() got non-standard CRef !?\n%s",FCT,toStr().data);

         mxSetFieldByNumber(S,k,0, QS.t   .toMx('t'));
         mxSetFieldByNumber(S,k,1, QS.qs  .toMx());
         mxSetFieldByNumber(S,k,2, QS.qdir.toMx());
         mxSetFieldByNumber(S,k,4, S2     .toMx());
      }
      else {
         const CDATA_TQ &C=*cgr;
         mxSetFieldByNumber(S,k,0, C.t    .toMx('t'));
         mxSetFieldByNumber(S,k,1, C.qs   .toMx());
         mxSetFieldByNumber(S,k,2, C.qdir .toMx());
         mxSetFieldByNumber(S,k,4, S2     .toMx());
      }



      if (flag==0) {
         wbvector<RTD> cgt;
         mxSetFieldByNumber(S,k,3, cgr->cstat.toMx(rtype.t));

         mxSetFieldByNumber(S,k,7, this->trace(0,0,cgt).toMx());
      }
      else {
         mxSetFieldByNumber(S,k,3, cgr->cgd.toMx());
      }
   }

   mxSetFieldByNumber(S,k,5, numtoMx(
      !isRefInit() ? (
        double(cgr ? cgr->cgd.nnz() : (cgw.len ? 1 : 0))
      ) : -1.
   ));

   mxSetFieldByNumber(S,k,6, cgw.toMx(30));



};



template <class TQ, class TD>
x3map<TQ,TD>& x3map<TQ,TD>::initCtr(const char *F, int L,
   const CRef<TQ> &A_, const ctrIdx &ica_,
   const CRef<TQ> &B_, const ctrIdx &icb_, wbperm &cgp,
   char xflag
){
   if (!A_.cgr || !B_.cgr) wblog(FL,
      "ERR %s() got empty CRef data (0x%lX, 0x%lX)",FCT,A_.cgr,B_.cgr);

   const CDATA_TQ &A=(*A_.cgr), &B=(*B_.cgr);
   ctrIdx ica(ica_), icb(icb_); { A_.adapt(ica); B_.adapt(icb); }
   QSet<TQ> &Qc=(QSet<TQ>&)c;

   unsigned ka,kb, i,j, nq=A.t.qlen(),
      ra=A.rank(F_L), rb=B.rank(F_L),
      rc, l=ra+rb;
   const wbperm &pa=A_.cgp, &pb=B_.cgp;

   char sa=(ica.conj!=0 ? -1 : +1),
        sb=(icb.conj!=0 ? -1 : +1), xd=-sa*sb;

   char ma[l], *mb=ma+ra; memset(ma,0,l*sizeof(char));
   l=0;

#ifndef NSAFEGUARDS
   if (A.t!=B.t || !nq) wblog(F_L,"ERR %s() got symmetry mismatch "
      "(%s, %s)", FCT, A.qStr().data, B.qStr().data);
   if (ra>32 || rb>32) wblog(F_L,
      "ERR %s() unexpected QSpace tensors of rank %d, %d !?",FCT,ra,rb);

   if (pa.len && !pa.isValidPerm(ra)) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pa.len,ra);
   if (pb.len && !pb.isValidPerm(rb)) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pb.len,rb);

   if (ica.len!=icb.len || ica.len>ra || icb.len>rb) wblog(FL,
      "ERR %s() invalid set of contraction indices (%d/%d; %d,%d)",
      FCT,ica.len,icb.len,ra,rb);
   if (A.qdir.len!=ra) wblog(FL,
      "ERR %s() invalid qdir set (A: %d/%d)",FCT,A.qdir.len,ra);
   if (B.qdir.len!=rb) wblog(FL,
      "ERR %s() invalid qdir set (B: %d/%d)",FCT,B.qdir.len,rb);
#endif

   for (i=0; i<ica.len; ++i) { j=ica[i];
      if (j>=ra) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,ica.len);
      if (!ma[j]) { ma[j]=i+1; } else wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,ica.len
      );
   }
   for (i=0; i<icb.len; ++i) { j=icb[i];
      if (j>=rb) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,icb.len);
      if (!mb[j]) { mb[j]=i+1; } else wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,icb.len);
      if (A.qdir[ica[i]] != xd*B.qdir[j]) wblog(FL,
         "ERR %s() matching in/out pairs of indices required\n"
         "A: %s @ [%s]\nB: %s @ [%s] having xd=%d",FCT,
         A.toStr().data, ica.toStr().data,
         B.toStr().data, icb.toStr().data,xd
      );
   }

   ka=ra-ica.len;
   kb=rb-icb.len; rc=ka+kb;

   Qc.t=A.t; Qc.qdir.init(rc); Qc.qs.init(nq*rc); 

   cgp.init(); a.RefInit(0,0,A);
   pab.init(); b.RefInit(0,0,B);


   if (rc>32) wblog(FL,
      "ERR %s() unexpected QSpace rank (rc=%d)",FCT,rc);

   if (ka) { for (i=0; i<ra; ++i) { j=pa.el1(i);
      if (!ma[j]) { ma[j]=-(++l); }
   }}
   if (kb) { for (i=0; i<rb; ++i) { j=pb.el1(i);
      if (!mb[j]) { mb[j]=-(++l); }
   }}

   if (l!=rc) wblog(FL,"ERR %s() %d != %d + %d !??",FCT,l,ka,kb);

   wbperm P0(rc);
   PERM_T *p0=P0.data; l=0;


   const char *da=A.qdir.data, *db=B.qdir.data; char *dc=Qc.qdir.data;
   const TQ *qa=A.qs.data, *qb=B.qs.data; TQ *qc=Qc.qs.data; l=0;

   if (ka) { for (i=0; i<ra; ++i) { if (ma[i]<0) {
      dc[l]=sa*da[i]; p0[l++]=-ma[i]-1;
      MEM_CPY<TQ>(qc, nq, qa+i*nq); qc+=nq;
   }}}
   if (kb) { for (i=0; i<rb; ++i) { if (mb[i]<0) {
      dc[l]=sb*db[i]; p0[l++]=-mb[i]-1;
      MEM_CPY<TQ>(qc, nq, qb+i*nq); qc+=nq;
   }}}

   Qc.Sort(&pab,&conj);

   cgp.init(rc);
   if (pab.len)
        for (i=0; i<rc; ++i) { cgp[p0[pab[i]]]=i; }
   else for (i=0; i<rc; ++i) { cgp[p0[    i ]]=i; }
   if (cgp.isIdentityPerm()) cgp.init();


   char zflag=0, largeD=0;
   size_t DX=(1<<20);

   if (rc==2 
   #ifndef NSAFEGUARDS
      && (A.cgd.numel()>DX || B.cgd.numel()>DX)
   #endif
   ){
      const TQ *q=Qc.qs.data; unsigned m=Qc.qs.len/2;
      if (!Qc.qs.len || Qc.qs.len%2) wblog(FL,
         "ERR %s() %s !?",FCT,Qc.toStr().data);

      if (Qc.qdir=="+-") {
         for (i=0; i<m; ++i) { if (q[i]!=q[m+i]) { zflag=1; break; }}
      }
      else if (Qc.t.isSUN() && (Qc.qdir=="++" || Qc.qdir=="--")) {
         unsigned l=Qc.qs.len-1;
         for (i=0; i<m; ++i) { if (q[i]!=q[l-i]) { zflag=1; break; }}
      }
   }

   if (xflag && !zflag) {
      if (A.isRefInit()) A_.LoadRef(FL);
      if (B.isRefInit()) B_.LoadRef(FL);

      largeD=1*(A.cgd.D.len>DX) + 2*(B.cgd.D.len>DX);
      if (largeD && CG_VERBOSE>5) { largeD+=8; wblog(FL,
         "START %s()\n    %-35s @ %s (%s)\n    %-35s @ %s (%s)",FCT,
         A.toStr().data, ica.toStr().data, A.cgd.sizeStr().data,
         B.toStr().data, icb.toStr().data, B.cgd.sizeStr().data
      ); doflush(); }
   }

   RTD cfac=(zflag ? 0 : 1);

   if (!rc) {
      if (A_.cgr!=B_.cgr || ra!=ica.len || rb!=icb.len) wblog(FL,
         "ERR %s() invalid contraction to scalar\n%s\n%s",FCT,
         A_.cgr->toStr().data, B_.cgr->toStr().data); 

      if (xflag) {
         c.cgd.init();

         if (!zflag) {
         cfac=A.cgd.contract(
            F,L, wbindex(ica), B.cgd, wbindex(icb), c.cgd); }

         if (ABS(double(cfac))<CG_EPS1) wblog(FL,
            "ERR %s() got small x3 !?",FCT,double(cfac));

         c.cgd*=cfac; c.cgd.toFull(x3);
         c.cgd.init();
         x3.Reshape(A.getOM(),B.getOM(),1);
         x3.SkipTiny(CG_SKIP_REPS);

         rtype=CGR_CTR_SCALAR; cgr=NULL; 
      }
      else { x3.init(); rtype=CGR_DEFAULT; cgr=NULL; }

      if (largeD) {
         if (largeD & 1) gCG.BUF[(QSet<TQ>&)a].Reduce2Ref(FL);
         if (largeD & 2) gCG.BUF[(QSet<TQ>&)b].Reduce2Ref(FL);
         if (largeD & 8) wblog(FL,"FINISHED %s() contracted to "
            "scalar (@ %.3g)",FCT, x3.data?double(x3[0]):-1);
         doflush();
      }
      return *this;
   }

   if (xflag) {

      CDATA_TQ Cx(Qc);
      wbperm Pc;

      unsigned ma=A.gotOM(), mb=B.gotOM(),
         rcM = rc + (ma?1:0) + (mb?1:0);

      if (!ma) {
         if (!A.cgd.SIZE.len && !A_.isDiagCGC()) wblog(FL,
            "ERR %s() got empty A: %s !?",FCT, A.toStr().data);
         Pc.init(pab.len ? rcM : 0);
      }
      else {
         Pc.init(rcM).Cycle(ka,ka+kb);
      }

      if (!mb) {
         if (!B.cgd.SIZE.len && !B_.isDiagCGC()) wblog(FL,
            "ERR %s() got empty B: %s !?",FCT, B.toStr().data
         );
      }

      if (pab.len) {
         wbvector<PERM_T> px(pab.len,Pc.data);
         for (i=0; i<pab.len; ++i) Pc[i]=px[pab[i]];
      }

      if (!zflag) {
      cfac=A.cgd.contract(F,L,
         wbindex(ica), B.cgd, wbindex(icb), Cx.cgd, NULL, &Pc); }

      if (abs(double(cfac))<CG_EPS2) {
         if (zflag || CG_VERBOSE>5) { sprintf(str,"%s",Qc.toStr().data); }
         Qc.qdir.init(); Qc.qs.init();
         x3.init();
         cgr=NULL; rtype=CGR_CTR_ZERO;

         if (largeD) {
            if (largeD & 1) gCG.BUF[(QSet<TQ>&)a].Reduce2Ref(FL);
            if (largeD & 2) gCG.BUF[(QSet<TQ>&)b].Reduce2Ref(FL);
            if (largeD & 8) {
               if (zflag) wblog(FL,"FINISHED %s() "
                  "contracted to zero (by QSet!)\n%s",FCT,str);
               else wblog(FL,"FINISHED %s() "
                  "contracted to zero (@ %.4g)",FCT,double(cfac));
            }; doflush();
         }
         else if (zflag && CG_VERBOSE>5) wblog(FL,
         " *  contracted to zero (by QSet)\n    c: %s",FCT,str);

         return *this;
      }



      CDATA_TQ &Cb=gCG.getBUF(0,0,Qc,'l');
      unsigned mc=Cb.getOM();

      if (!mc || Cb.isEmpty()) {
         if (CG_VERBOSE>2) wblog(PFL,
         "[+] CBUF[%03d] new/dec: %s",gCG.BUF.size(), Qc.toStr().data);
         Cb=Qc; Cb.cstat.init(CGD_FROM_DEC);
      }
      else if (Cb!=Qc) wblog(FL,
         "ERR %s() got QSet mismatch\n%s <> %s",
         FCT,Cb.toStr().data,Qc.toStr().data
      );

      cgr=&Cb;


      if (ma<1) ma=1;
      if (mb<1) mb=1;

      wbvector< CRef<TQ> > cc(ma*mb);
      cdata<TD> cx;

      wbIndex Idx; i=0; l=0;

      while (Cx.getCG_set(FL,Idx,cx)) {
         if (cx.D.len) { l+=cx.D.len;
            cc[i].initDecompose(FL,cx,&Cb);
            if (mc<cc[i].cgw.len) mc=cc[i].cgw.len;
         }
         ++i;
      }

      if (l!=Cx.cgd.D.len) wblog(FL,
         "ERR %s() OM inconsistency (l=%d/%d)",FCT,l,Cx.cgd.D.len
      );

      c.RefInit(0,0,*cgr);

if (ma!=a.getOM(FL) || mb!=b.getOM(FL) || mc!=c.getOM(FL)) wblog(FL,
"ERR %s() %d/%d %d/%d %d/%d/%d", FCT, ma, a.getOM(FL), mb, b.getOM(FL), mc,
c.getOM(FL), cgr->getOM(FL));

      x3.init(cc.len,mc);

      for (i=0; i<cc.len; ++i) { l=cc[i].cgw.len;
      for (j=0; j<l     ; ++j) { x3(i,j)=cc[i].cgw[j]; }}

      x3.Reshape(ma,mb,mc);


      if (fabs(double(cfac))<1E-3) wblog(FL,
         "WRN %s() got small cfac=%.4g",FCT, double(cfac));
      x3*=cfac;

      double e=x3.SkipTiny(CG_SKIP_REPS);
      if (e && CG_VERBOSE>2) wblog(FL," *  %s() x3 @ %.3g",FCT,e);
   }
   else {
      x3.init(); c.cgd.init(); rtype=CGR_DEFAULT;
      cgr=NULL;
   }

   #ifndef NSAFEGUARDS
     if (cgr && !cgr->cstat.ctime)
     wblog(FL,"ERR %s() %s",FCT,cgr->cstat.toStr('V').data);
   #endif

   if (largeD) { if (cgr->cgd.D.len>DX) largeD+=4;
      if (largeD & 1) gCG.BUF[(QSet<TQ>&)a].Reduce2Ref(FL);
      if (largeD & 2) gCG.BUF[(QSet<TQ>&)b].Reduce2Ref(FL);
      if (largeD & 4) gCG.BUF[(QSet<TQ>&)c].Reduce2Ref(FL);
      if (largeD & 8) wblog(FL,
         "FINISHED %s() %s",FCT,x3Str().data);
      doflush();
   }
   return *this;
};


template <class TQ, class TD>
wbstring x3map<TQ,TD>::toStr(char vflag) const {

   wbstring s_(256);
   unsigned l=0, n=s_.len-1; char *s=s_.data;

   l=snprintf(s,n,"%s",((QSet<TQ>&)c).toStr().data);

   if (l<n && (pab.len || conj))
      l+=snprintf(s+l,n-l,", pab=[%s]%s",(pab+1).toStr().data,conj?"*":"");

   if (vflag && l<n) {
      if (vflag>2) {
         if (vflag=='v') vflag=1; else
         if (vflag=='V') vflag=2;
         else wblog(FL,"ERR %s() invalid vflag=%d !?",FCT,vflag);
      }
      l+=snprintf(s+l,n-l," => %s",x3Str(vflag-1).data);
   }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s_;
};


template <class TQ, class TD>
wbstring x3map<TQ,TD>::x3Str(char vflag) const {

   wbstring s_;

   unsigned l=0, n, n3=x3.numel(); char *s;
   wbstring dstr=wbvector<RTD>(n3,x3.data,'r').toStr();

   s_.init(strlen(dstr.data)+32);
   n=s_.len-1; s=s_.data;
  
   if (rtype==CGR_CTR_ZERO) {
      if (n3) wblog(FL,
         "ERR %s() got x3 data for ZERO CData (n3=%d) !?",FCT,n3);
      l=snprintf(s,n,"zero (%s)",rtype.tostr());
   }
   else if (rtype==CGR_CTR_SCALAR) {
      if (n3==1) 
           l=snprintf(s,n,"scalar [%s]",dstr.data); else
      if (n3<=4 || vflag) 
           l=snprintf(s,n,"scalar [%s] (%s)",dstr.data,x3.sizeStr().data);
      else l=snprintf(s,n,"scalar [%s]",dstr.data);
   }
   else {
      wbstring sstr=x3.sizeStr();
      if (n3==1)
         l=snprintf(s,n,"x3=%s (%s)",dstr.data,sstr.data);
      else if (n3<=4 || vflag)
         l=snprintf(s,n,"x3=[%s] (%s)",dstr.data,sstr.data);
      else l=snprintf(s,n,"x3 (%s)",sstr.data);
   }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d; %d/%d)%N%N%s...%N",
      FCT,l,n,x3.nnz(),n3,s);
   return s_;
};


template <class TQ, class TD>
mxArray* x3map<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={ "idc","idh"
      "a","ica","b","icb","c", "x3",
      "Pab","conj","rtype"
   };
   return mxCreateStructMatrix(m,n,10,fields);
};


template <class TQ, class TD>
void x3map<TQ,TD>::add2MxStruct(
   mxArray *S, unsigned l, const cgc_contract_id<MTI> *idc) const {

   if (!S) wblog(FL,"ERR %s() got null mxArray !?",FCT);
   if (l>=mxGetNumberOfElements(S)) wblog(FL,"ERR %s() "
      "index out of bounds (%d/%d)",FCT,l,mxGetNumberOfElements(S));

   ctrIdx ica,icb;

   wbvector<size_t> s3;
   if (cgr) cgr->getSize(s3,'b'); else s3.init();

#ifndef NSAFEGUARDS
   a.rank(FL); b.rank(FL); c.rank(FL);
#endif

   if (idc) {
      CData<TQ,RTD> a_,b_;
      const auto it=gCX.XBUF.find(*idc);
      if (it==gCX.XBUF.end()) wblog(FL,"ERR %s() key not found !?",FCT); 
      if (&(it->second)!=this) wblog(FL,"ERR %s() got invalid key !?",FCT); 

      idc->extract(a_,ica,b_,icb);

#ifndef NSAFEGUARDS
      if (a_!=(QSet<TQ>&)a || b_!=(QSet<TQ>&)b) wblog(FL,
         "ERR %s() QSet mismatch\n   %s | %s\n<> %s | %s",FCT,
         a_.toStr().data, b_.toStr().data,
         a.toStr().data, b.toStr().data
      );
#endif
   }


   mxSetFieldByNumber(S,l,1, a.toMx());
   mxSetFieldByNumber(S,l,3, b.toMx());

   if (idc) {
      mxSetFieldByNumber(S,l,0, idc->toMx());
      mxSetFieldByNumber(S,l,2, ica.toMx());
      mxSetFieldByNumber(S,l,4, icb.toMx());
   }

   if (rtype!=CGR_CTR_ZERO) {
      mxSetFieldByNumber(S,l,5, c .toMx());
   }

   mxSetFieldByNumber(S,l,6, x3.toMx());
   mxSetFieldByNumber(S,l,7, pab.toMx());
   mxSetFieldByNumber(S,l,8, numtoMx(conj));
   mxSetFieldByNumber(S,l,9, numtoMx(rtype.t));
};


template <class TM>
template <class TQ, class TD> inline
cgc_contract_id<TM>::cgc_contract_id(
   const CData<TQ,TD> &A, const ctrIdx &ica,
   const CData<TQ,TD> &B, const ctrIdx &icb
){
   unsigned n=256, ra=A.rank(FL), rb=B.rank(FL);
   char s0[n], *s=s0; unsigned i, l, m=sizeof(TM); wbstring s1;
   wbvector<unsigned> S;

   if (A.t!=B.t) wblog(FL,"ERR %s() qtype mismatch: "
      "%s <> %s",FCT, A.qStr().data, B.qStr().data);
   if (!A.qs.len || !B.qs.len || !ica.len || ica.len!=icb.len)
       wblog(FL,"ERR %s() got empty data set !??",FCT);
   if (A.qdir.len!=ra || B.qdir.len!=rb) wblog(FL,
      "ERR %s() qdir length mismatch (%d/%d; %d/%d",
      FCT, A.qdir.len, ra, B.qdir.len, rb);
   if (ica.len>ra || icb.len>rb) wblog(FL,
      "ERR %s() invalid contraction (%d/%d; %d/%d)",
      FCT, ica.len,ra, icb.len,rb);
   if (n<(l=4*(A.qdir.len+B.qdir.len)+16)) wblog(FL,
      "ERR %s() string size possibly too small (%d/%d)",FCT,n,l);

   set_CData(FL,s,(n-ica.len-1)/2, A);
   set_vec(FL,s, ica);
   set_val(FL,s, ica.conj); l=s-s0;

   set_CData(FL,s,n-l-icb.len-1, B, 0);
   set_vec(FL,s, icb);
   set_val(FL,s, icb.conj); l=s-s0;

   if (int(l)<2 || l+sizeof(TM)>=n) wblog(FL,"ERR %s() "
      "unexpected string length (%d/%d)",FCT,l,n);

   if (l%m) {
      for (l=m-(l%m), i=0; i<l; ++i) { s[i]=0; }
      s+=l;
   }

   l=s-s0; 
   if (l%m) wblog(FL,"ERR %s() l=%d @ %d",FCT,l,m);

   this->init(1+(l-1)/m,(TM*)s0);
};


template <class TM>
template <class TQ, class TD> inline
void cgc_contract_id<TM>::extract(
   CData<TQ,TD> &A, ctrIdx &ica,
   CData<TQ,TD> &B, ctrIdx &icb) const {

   const char *s = (const char*) this->data;

   get_CData(FL,s,A); A.rank(FL);
   get_vec(FL,s,ica);
   get_val(FL,s,ica.conj);

   get_CData(FL,s,B,0); B.t=A.t; B.rank(FL);
   get_vec(FL,s,icb);
   get_val(FL,s,icb.conj);

   if (ica.anyGT(20) || icb.anyGT(20)) wblog(FL,
      "ERR %s() unexpected contraction index (%s; %s)",
      FCT, ica.toStr().data, icb.toStr().data
   );

   const x3map<TQ,RTD> &M=gCX.XBUF[*this];
   if (!M.isEmpty()) {
      const wbvector<size_t> &Sa=M.a.cgd.SIZE, &Sb=M.b.cgd.SIZE;
      unsigned la=A.cgd.SIZE.len, lb=B.cgd.SIZE.len;

      A.cstat=M.a.cstat;
      if (la>Sa.len || A.cgd.D.len) wblog(FL,
         "ERR %s() unexpected OM setting (%s <> %s; %d)",
         FCT,M.a.sizeStr().data,A.sizeStr().data,A.cgd.D.len);
      else if (la<Sa.len) { A.cgd.SIZE=Sa; }

      B.cstat=M.b.cstat;
      if (lb>Sb.len || B.cgd.D.len) wblog(FL,
         "ERR %s() unexpected OM setting (%s <> %s; %d)",
         FCT,M.b.sizeStr().data,B.sizeStr().data,B.cgd.D.len);
      else if (lb<Sb.len) { B.cgd.SIZE=Sb; }
   }
};


template <class TQ>
CRef<TQ>& X3Map<TQ>::contract( const char *F, int L,
   const CRef<TQ> &A, const ctrIdx &ica_,
   const CRef<TQ> &B, const ctrIdx &icb_, CRef<TQ> &C
){


   if (!A.cgr || !B.cgr) wblog(F_L,
      "ERR %s() got empty cgref data (0x%lX, 0x%lX)",FCT,A.cgr,B.cgr);
   if (A.cgr->t!=B.cgr->t) wblog(F_L,
      "ERR %s() got mixed symmetry (%s, %s)",FCT,A.cgr->qStr().data,B.cgr);

   if (A.isAbelian()) {
      if (A.cgw.len!=1 || B.cgw.len!=1) wblog(F_L,
         "ERR %s() unexpected scalar CRef data (%d,%d)",
         FCT,A.cgw.len,B.cgw.len);
      C=A; C.cgw[0]*=B.cgw[0];

      QSet<TQ> Q(F_L, A, ica_, B, icb_);
      Q.checkQ_abelian(F_L);

      Q.Sort(&C.cgp, &C.conj,'i');

      CDATA_TQ &S=gCG.getBUF(FL,Q);

      if (S.isEmpty()) {
         S=Q; S.cstat.init(CGD_ABELIAN); S.cgd.init();
         if (CG_VERBOSE>2) 
         wblog(PFL,"[+] add2BUF() new/ctr: %s",S.toStr('v').data);
      }
      else { 
         if (S!=Q) wblog(F_L,
            "ERR %s() got CData inconsistency\n%s <> %s",
            FCT, S.toStr().data, Q.toStr().data);
         if (S.getScalar(F_L)!=1) wblog(F_L,
            "ERR %s() invalid scalar cgd\n%s",FCT,S.toStr().data
         );
      }

      C.cgr=&S;

      return C;
   }


   ctrIdx ica(ica_), icb(icb_); wbperm cgp;
   A.adapt(ica);
   B.adapt(icb);

   cgc_contract_id<MTI> idc(*A.cgr,ica, *B.cgr,icb);
   x3map<TQ,RTD> &Cm=XBUF[idc];





   char sq[3]={0,0,0};
   unsigned i,j,k, calc=0; int loaded=0,
      ma=A.cgw.len, Ma=A.getOM(),
      mb=B.cgw.len, Mb=B.getOM(), mc=0, Mc=0;

   if (!ma || ma>Ma || !mb || mb>Mb) wblog(FL,
      "ERR %s() invalid OM data (%dx%d; %dx%d)",FCT,ma,mb,Ma,Mb);

   if (Cm.isEmpty()) {
      if (gStore.load_XMap(0,0,idc)>0) loaded=1;
   };

   if (Cm.isEmpty()) { calc=1;
      Cm.initCtr(FL,A,ica_,B,icb_, cgp,'!');

      if (Cm.cgr) {
         mc=Cm.x3.SIZE[2]; Mc=Cm.cgr->getOM();
      }

      if (CG_VERBOSE>2 || (CG_VERBOSE>1 && (
           Cm.x3.numel()>1 || Cm.rtype==CGR_CTR_SCALAR ||
           (A.cgr->qdir.len + B.cgr->qdir.len - ica.len - icb.len) > 3
      ))) wblog(PFL,
         "[+] XBUF[%03d] %s\n    | %s @ %s->%s\n    | %s @ %s->%s\n"
         "    > %s", XBUF.size(), Cm.x3Str().data,
         A.toStr().data, ica_.toStr().data, ica.toStr().data,
         B.toStr().data, icb_.toStr().data, icb.toStr().data,
         Cm.toStr().data
      );

   }
   else {

      x3map<TQ,RTD> Cx;
      Cx.initCtr(FL,A,ica_,B,icb_, cgp, 0);

      if (loaded) {
         char e1=!Cx.a.sameAs(Cm.a,'L'), e2=!Cx.b.sameAs(Cm.b,'L');
         if (e1 || e2) { printf("\n");
            wblog(FL,
               "==> %s() x3map inconsistency (id=%ld)\n"
               "   %s @ %s -> %s\n   %s @ %s -> %s\n=> %s\n=! %s",
               FCT, MHash<MTI>()(idc),
               A.toStr().data, ica_.toStr().data, ica.toStr().data,
               B.toStr().data, icb_.toStr().data, icb.toStr().data
            ); printf("\n");

            if (e1) wblog(FL,
               "==> a: %s\n=! %s",Cx.a.toStr().data, Cm.a.toStr().data);
            if (e2) wblog(FL,
               "==> b: %s\n=! %s",Cx.b.toStr().data, Cm.b.toStr().data);
            wblog(FL,"ERR %s()",FCT);
         }
      }

      if ((Cm.rtype!=CGR_CTR_ZERO &&
          (QSet<TQ>&)(Cx.c)!=(QSet<TQ>&)(Cm.c)) || Cx.pab!=Cm.pab)
         wblog(FL,"ERR %s() QSet inconsistency (id=%ld)\n"
         "   %s @ %s -> %s\n   %s @ %s -> %s\n=> %s\n=! %s",
         FCT, MHash<MTI>()(idc),
         A.toStr().data, ica_.toStr().data, ica.toStr().data,
         B.toStr().data, icb_.toStr().data, icb.toStr().data,
         Cx.toStr().data, Cm.toStr().data
      );

      if (Cm.cgr && !Cm.x3.isEmpty()) {
         if (Cm.x3.SIZE.len!=3) wblog(FL,
            "ERR %s() invalid x3 (%s)",FCT,Cm.x3.sizeStr().data);
         mc=Cm.x3.SIZE[2];
         Mc=Cm.cgr->getOM();

         sq[0]=NUMCMP(Ma,Cm.x3.SIZE.data[0]);
         sq[1]=NUMCMP(Mb,Cm.x3.SIZE.data[1]);
         sq[2]=NUMCMP(Mc,mc);

         if (sq[2]<0) wblog(FL,
            "ERR %s() got x3 with larger OM !? %s (%dx%dx%d)%N%N   "
            "A: %-25s| a: %-25s @ %-3s%N   "
            "B: %-25s| b: %-25s @ %-3s%N   "
            " %-25s => c: %-25s%N   %s%N", FCT,Cm.x3.sizeStr().data,Ma,Mb,Mc,
            A.cgr->toStr().data, Cm.a.toStr().data, ica.toStr().data,
            B.cgr->toStr().data, Cm.b.toStr().data, icb.toStr().data,
            "",                  Cm.c.toStr().data,
            Cm.cgr->cstat.toStr('V').data
         );

         if (sq[0]>0 || sq[1]>0 || sq[2]>0) calc=2;
      }
      else {

         if (A.cgr->cstat!=Cm.a.cstat || B.cgr->cstat!=Cm.b.cstat) {
            if (CG_VERBOSE>1) wblog(FL,
               "TST %s() got altered cstat (%d,%d) -> recalc",FCT,
               A.cgr->cstat!=Cm.a.cstat, B.cgr->cstat!=Cm.b.cstat);
            calc=3;
         }

         if (Cm.x3.SIZE.len) {
            if (Cm.x3.SIZE[2]!=1) wblog(FL,"ERR %s() "
               "got unexpected x3 (%s)",FCT,Cm.x3.sizeStr().data);
            sq[0]=NUMCMP(Ma,Cm.x3.SIZE[0]);
            sq[1]=NUMCMP(Mb,Cm.x3.SIZE[1]);

            if (!calc && (sq[0] || sq[1])) wblog(FL,
               "ERR %s() got cstat inconsistency !?",FCT); 

            if (sq[0]<=0 && sq[1]<=0) calc=0;
         }
      }

      if (calc) {

         wbarray<RTD> x3_; Cm.x3.save2(x3_);
         if (x3_.isEmpty()) {
            if (Cm.rtype!=CGR_CTR_ZERO) wblog(FL,
               "ERR %s() got rtype=%s !?",FCT,Cm.rtype.tostr());
            x3_.init(size_t(0),size_t(0),size_t(0));
         }

         if (Mc<x3_.SIZE[2]) wblog(FL,
            "ERR %s() x3 got larger OM !? (%dx%dx%d <> %s)\n%s",
            FCT, Ma,Mb,Mc, x3_.sizeStr().data, Cm.cgr->toStr().data
         );

         Cm.initCtr(FL,A,ica_,B,icb_, cgp,'!');

         if (CG_VERBOSE>2) {
            wblog(PFL,"[u] XBUF[%03d] @ %s\n"
            "    %s @ %s -> %s\n    %s @ %s -> %s\n"
            "--> %s", XBUF.size(), Cm.x3Str().data,
            A.toStr().data, ica_.toStr().data, ica.toStr().data,
            B.toStr().data, icb_.toStr().data, icb.toStr().data,
            Cm.toStr().data);
         }

      #ifndef NSAFEGUARDS
         unsigned d1=x3_.SIZE[0], d2=x3_.SIZE[1], d3=x3_.SIZE[2];
         double e3;

         for (k=0; k<d3; ++k)
         for (j=0; j<d2; ++j)
         for (i=0; i<d1; ++i) x3_(i,j,k) -= Cm.x3(i,j,k);

         if ((e3=x3_.norm())>1E-12)
         wblog(FL,"ERR %s() got x3 data mismatch @ %.3g",FCT,e3); 
      #endif
      }
   }

   if (calc) gStore.save_XMap(FL,idc);

   if (Cm.rtype==CGR_CTR_ZERO) {
      if (Cm.cgr) wblog(FL,"ERR %s() got %s",FCT,Cm.cgr->toStr().data);
      C.init(); C.rtype=Cm.rtype;
      return C;
   }

   if (Cm.rtype==CGR_CTR_SCALAR) {
      const wbvector<unsigned> &S=Cm.x3.SIZE;

      if (Cm.cgr) wblog(FL,"ERR %s() got %s",FCT,Cm.cgr->toStr().data);
 
      if (S.len!=3 || 
          ma!=A.cgw.len || ma>S[0] || 
          mb!=B.cgw.len || mb>S[1] || S[2]!=1) wblog(FL,
         "ERR %s() unexpected OM setting (%dx%d; %dx%d; %s)",
         FCT, ma,mb, Ma,Mb, Cm.x3.sizeStr().data
      );

      C.init(); C.rtype=Cm.rtype;

      C.cgw.init(1); {
         RTD &ck=C.cgw.data[0];
         const RTD *a=A.cgw.data, *b=B.cgw.data;

         for (j=0; j<mb; ++j) {
         for (i=0; i<ma; ++i) { ck+=a[i]*b[j]*Cm.x3(i,j,0); }}
      }

      return C;
   }

   if (Cm.cgr) Mc=Cm.cgr->getOM();
   else wblog(FL,
      "ERR %s() got cgr==NULL !?\n%s",FCT,Cm.toStr().data);
   if (Cm.x3.SIZE.len!=3) wblog(FL,
      "ERR %s() invalid x3 (%s)",FCT,Cm.x3.sizeStr().data);
   if (!mc || mc>Mc) wblog(FL,
      "ERR %s() invalid mc=%d/%d !?",FCT,mc,Mc);

   if (ma>Cm.x3.SIZE[0] || mb>Cm.x3.SIZE[1] || mc>Cm.x3.SIZE[2])
      wblog(FL,"ERR x3map() size mismatch (m: %dx%dx%d <> %s)",
      ma,mb,mc, Cm.x3.sizeStr().data
   );

   C.cgr=Cm.cgr;
   C.cgp=cgp; C.conj=Cm.conj;

   if (calc) loaded=-loaded;
   if (loaded>0) {
      if (Cm.rtype!=CGR_DEFAULT) wblog(FL,
         "ERR %s() got rtype=%s !?",FCT,Cm.rtype.tostr());
      C.rtype=Cm.rtype;
   }
   else if (!calc) {
      C.rtype=Cm.rtype;
   }
   else C.rtype=CGR_DEFAULT;

   C.cgw.init(mc); {
      const RTD *a=A.cgw.data, *b=B.cgw.data;

      for (k=0; k<mc; ++k) { RTD &ck=C.cgw.data[k];
      for (j=0; j<mb; ++j) {
      for (i=0; i<ma; ++i) { ck+=a[i]*b[j]*Cm.x3(i,j,k); }}}
   }

   C.isRefInit();


   return C;
};



template <class TQ>
CRef<TQ>& X3Map<TQ>::contractDQ(const char *F, int L,
   const CRef<TQ> &A, CRef<TQ> &B
){
   if (B.cgr!=A.cgr || B.cgw.len!=A.cgw.len) wblog(FL,
      "ERR %s() invalid/unexpected B (0x%lx,0x%lx; %d/%d) !?",FCT,
      A.cgr, B.cgr, A.cgw.len, B.cgw.len);

   if (!A.cgr) {
      if (B.cgp.len) B.cgp.init(); B.conj=0;
      return B;
   }
   if (A.isAbelian()) {
      if (!A.cgr->isScalar() || A.cgw.len!=1 ||
         (A.cgr->cgd.D.len && A.cgr->cgd.D[0]!=1)) wblog(F_L,
         "ERR %s() unexpected scalar CRef data\n%s",
         FCT, A.toStr().data);
      return B;
   }

   unsigned r=A.rank(FL);
   if (r==2) {
      if (!A.isIdentityCGC()) wblog(FL,"ERR %s() got non-scalar "
         "rank-%d CGC (%s) !?",FCT,r,A.cgr->sizeStr().data);
      return B;
   }
   else if (r<2) wblog(FL,"ERR %s() "
      "got rank-%d CRef !?\n%s",FCT,r,A.cgr->sizeStr().data);

   unsigned i,j,k, m=A.cgw.len, M=A.getOM();
   char loaded=0, calc=0;

   wbperm pa(A.cgp,B.cgp);

   ctrIdx ica(pa,(A.conj!=0) ^ (B.conj!=0) ? 0 : 1);
   ctrIdx icb; { icb.Index(r); }
   ctrIdx ica_(ica), icb_(icb);
   wbperm cgp;

   A.adapt(ica_,'i');
   B.adapt(icb_,'i');

   cgc_contract_id<MTI> idc(*A.cgr,ica, *B.cgr,icb);
   x3map<TQ,RTD> &Cm=XBUF[idc];
   cdata__ a;

   if (!m || m>M) wblog(FL,
      "ERR %s() invalid OM data (%d/%d)",FCT,m,M);

   if (Cm.isEmpty()) {
      if (gStore.load_XMap(0,0,idc)>0) loaded=1;
   };

   if (Cm.isEmpty()) { calc=1;
      Cm.initCtr(FL,A,ica_,B,icb_, cgp,'!');

      if (Cm.cgr) wblog(FL,
         "ERR %s() got cgr data !?\n%s",FCT,Cm.cgr->toStr().data);
      if (Cm.c.qs.len || Cm.c.qdir.len) wblog(FL,
         "ERR %s() got non-empty c !?\n%s",FCT,Cm.c.toStr().data);
      if (Cm.x3.SIZE.len!=3 || Cm.x3.SIZE[2]!=1) wblog(FL,"ERR %s() "
         "got unexpected x3 data (%s)",FCT,Cm.x3.sizeStr().data);

      if (CG_VERBOSE>2 || (CG_VERBOSE>1 && Cm.x3.numel()>1)) wblog(PFL,
         "[+] XBUF[%03d]\n"
         "    | %s @ %s -> %s\n    | %s @ %s -> %s\n    > x3 is %s", 
         XBUF.size(),
         A.toStr().data, ica_.toStr().data, ica.toStr().data,
         B.toStr().data, icb_.toStr().data, icb.toStr().data,
         Cm.x3Str().data
      );
   }
   else {
      if (Cm.cgr || Cm.x3.SIZE.len!=3 || 
          Cm.x3.SIZE[0]!=Cm.x3.SIZE[1] || Cm.x3.SIZE[2]!=1)
      wblog(FL,"ERR %s() invalid x3 (%s)",FCT,Cm.x3.sizeStr().data);

      if (M!=Cm.x3.SIZE.data[0]){ calc=2;
         wbarray<RTD> x3_; Cm.x3.save2(x3_);
         if (M<x3_.SIZE.data[0]) wblog(FL,
            "ERR %s() x3map() size inconsistency (%dx%dx1 <> %s)",
            FCT,M,M, x3_.sizeStr().data
         );

         Cm.initCtr(FL,A,ica_,B,icb_, cgp,'!');

         if (CG_VERBOSE>2) wblog(PFL,
            "[u] XBUF[%03d] @ %s\n    %s @ %s -> %s\n--> %s",
            XBUF.size(), x3_.sizeStr().data, A.toStr().data,
            ica_.toStr().data, icb_.toStr().data, Cm.toStr().data
         );

      #ifndef NSAFEGUARDS
         unsigned d1=x3_.SIZE[0], d2=x3_.SIZE[1], d3=x3_.SIZE[2];
         double e3;

         for (k=0; k<d3; ++k)
         for (j=0; j<d2; ++j)
         for (i=0; i<d1; ++i) x3_(i,j,k) -= Cm.x3(i,j,k);

         if ((e3=x3_.norm())>1E-12)
         wblog(FL,"ERR %s() got x3 data mismatch @ %.3g",FCT,e3); 
      #endif
      }
   }

   if (calc) gStore.save_XMap(FL,idc);

   B.cgw.init(M);

   if (Cm.rtype==CGR_CTR_SCALAR) {
      const wbvector<unsigned> &S=Cm.x3.SIZE;

      if (Cm.cgr) wblog(FL,"ERR %s() got %s",FCT,Cm.cgr->toStr().data);

      if (m>M) wblog(FL,"ERR %s() invalid OM data (%d/%d)",FCT,m,M);
      if (S.len!=3 || S[0]!=M || S[1]!=M || S[2]!=1) wblog(FL,
         "ERR %s() OM out of bounds (%dx%dx1; %s)",
         FCT,M,M, Cm.x3.sizeStr().data
      );

      if (S.len==1) { B.cgw[0]=A.cgw[0]*Cm.x3[0]; }
      else {
         const RTD *a=A.cgw.data; RTD *b=B.cgw.data;

         for (j=0; j<M; ++j) {
         for (i=0; i<m; ++i) { b[j]+=a[i]*Cm.x3(i,j,0); }}
      }

      double x[2]={ double(A.cgw.norm2()), double(B.cgw.norm2()) };
      
      if (x[0]<CG_SKIP_EPS2 || x[1]<CG_SKIP_EPS2) wblog(FL,
         "ERR %s() got small cgw norm (%.3g, %.3g) !?",FCT,x[0],x[1]);
      if (ABS(x[1]-x[0])<CG_SKIP_EPS2) {
         return B;
      }

      a.init(A,0);
         if (!pa.isIdentityPerm()) a.Permute(pa);
         if (ica.conj) a.Conj();
      cdata__ b; b.init(B,0);

      a-=b;
   }
   else {
      if (Cm.rtype!=CGR_CTR_ZERO) wblog(FL,
         "ERR %s() unexpected rtype=%s !?",FCT,Cm.rtype.tostr()); 
      a.init(A,0);
   }

   RTD cfac=a.NormSignC(F,L);

   double x=double(cfac); if (ABS(x)<1E-3) wblog(FL,
      "WRN %s() got small new coefficient !?? (%.3g)",FCT,x);
   B.cgw.Append(cfac);

   ((CDATA_TQ*)(B.cgr))->AddMultiplicity(FL,a);
   ((CDATA_TQ*)(B.cgr))->cstat.t=CGD_FROM_DEC;

   if (CG_VERBOSE>2 || (CG_VERBOSE>1 && r<4)) wblog(PFL,
      "[+] CBUF[%d] new/dec: %s",gCG.BUF.size(),B.cgr->toStr().data);

   gStore.save_CData(FL,*B.cgr);
   if (B.cgr->qdir.len==3) wblog(FL,"WRN %s() "
      "decomposing rank-3 CData !?\n%s",FCT,B.cgr->toStr().data);

#ifndef NSAFEGUARDS
#endif

   return B;
};


template <class TQ>
mxArray* X3Map<TQ>::toMx() const {

   map<QType,
      map <wbvector<char>,
         map <wbvector<unsigned>,
            const cgc_contract_id<MTI>*
   > > > X1;

   wbvector<char> r3(3);
   wbvector<unsigned> S3; wbvector<size_t> s3;

   for (auto it=XBUF.begin(); it!=XBUF.end(); ++it) {
      const x3map<TQ,RTD> &M=it->second;

      r3[0]=M.a.rank(FL);
      r3[1]=M.b.rank(FL);
      r3[2]=M.c.rank(FL);

      if (M.cgr) M.cgr->getSize(s3,'b'); else s3.init();

      S3.Cat(&M.a.cgd.SIZE, &M.b.cgd.SIZE, &s3);

      X1[M.c.t][r3][S3]=&(it->first);
   }

   mxArray *S1=mxCreateStructMatrix(1,1,0,NULL);

   for (auto it1=X1.begin(); it1!=X1.end(); ++it1) {
      const auto &X2=it1->second;

      mxArray *S2=mxCreateStructMatrix(1,1,0,NULL);

   for (auto it2=X2.begin(); it2!=X2.end(); ++it2) {
      const auto &X3=it2->second;

      const wbvector<char> &r3=it2->first;
      mxArray *S3=x3map<TQ,RTD>().mxCreateStruct(X3.size(),1);
      unsigned l=0;

      for (auto it3=X3.begin(); it3!=X3.end(); ++it3, ++l) {

         const cgc_contract_id<MTI> &idc=*(it3->second);
         const auto itM=XBUF.find(idc);
         if (itM==XBUF.end()) wblog(FL,"ERR %s() key not found !?",FCT); 
         const x3map<TQ,RTD> &M=itM->second;

         CData<TQ,RTD> a_,b_;
         ctrIdx ica,icb;

         idc.extract(a_,ica,b_,icb);

         if (a_!=(QSet<TQ>&)M.a || b_!=(QSet<TQ>&)M.b) wblog(FL,
            "ERR %s() QSet mismatch\n   %s | %s\n<> %s | %s",FCT,
            a_.toStr().data, b_.toStr().data,
            M.a.toStr().data, M.b.toStr().data
         );

         M.add2MxStruct(S3,l,&idc);
      }

      sprintf(str,"rank%d%d%d",r3[0],r3[1],r3[2]);
      mxAddField2Scalar(FL,S2,str,S3);
   }
      mxAddField2Scalar(FL,S1,it1->first.toStr('t').data,S2);
   }

   wbvector<unsigned> nc(XBUF.bucket_count());
   for (unsigned i=0; i<nc.len; ++i) nc[i]=XBUF.bucket_size(i);

   mxAddField2Scalar(FL,S1,"MHash", MXPut(0,0)
      .add(XBUF.max_load_factor(),"max_load_factor")
      .add(XBUF.load_factor(),"load_factor")
      .add(nc,"filling_buckets")
      .add(XBUF.size(),"nr_entries")
      .add(nc.max(0),"max_collisions")
   .toMx());

   return S1;
};


template <class TQ>
unsigned X3Map<TQ>::add(
   const char *F, int L,
   const cgc_contract_id<MTI> &idc, const mxArray* S, unsigned k) {

   x3map<TQ,RTD> &Cm=gCX.XBUF[idc];

   CDATA_TQ a,b; ctrIdx ica,icb; int n;

   cgc_contract_id<MTI> id2;

   if (!Cm.a.isEmpty() || !Cm.b.isEmpty() || !Cm.x3.isEmpty())
      wblog(FL,"ERR %s() already got non-empty XBUF entry !?",FCT);

   if (!S) wblog(FL,"ERR %s() got null mxArray !?",FCT);
   if ((n=mxGetNumberOfElements(S))!=1) wblog(FL,
      "ERR %s() got invalid number of elements (%d)",FCT,n);

   id2.init(FL, mxGetFieldByNumber(S,k,0));
   if (idc!=id2) wblog(FL,"ERR %s() got contract ID mismatch !?",FCT); 

   idc.extract(a,ica,b,icb);

   Cm.a  .init(FL, mxGetFieldByNumber(S,k,1), 0,0,0);
      ica.init(FL, mxGetFieldByNumber(S,k,2));
   Cm.b  .init(FL, mxGetFieldByNumber(S,k,3), 0,0,0);
      icb.init(FL, mxGetFieldByNumber(S,k,4));
   Cm.c  .init(FL, mxGetFieldByNumber(S,k,5), 0,0,0);
   Cm.x3 .init(FL, mxGetFieldByNumber(S,k,6));
   Cm.pab.init(FL, mxGetFieldByNumber(S,k,7), 1);

   mxGetNumber(mxGetFieldByNumber(S,k,8), Cm.conj);
   mxGetNumber(mxGetFieldByNumber(S,k,9), Cm.rtype.t);

   if (Cm.rtype==CGR_CTR_ZERO || Cm.rtype==CGR_CTR_SCALAR) {
      if (!Cm.c.qdir.isEmpty() || !Cm.c.qs.isEmpty()) wblog(FL,
     "ERR %s() got %s (%s) !?",FCT,Cm.c.toStr().data,Cm.rtype.tostr());
      if (Cm.c.t.isUnknown()) Cm.c.t=Cm.a.t;
   }
   else {
      if (Cm.rtype!=CGR_DEFAULT) wblog(FL,
     "WRN %s() got rtype=%s",FCT,Cm.rtype.tostr());
   }

   n=0;
   if (Cm.rtype==CGR_CTR_ZERO) {
      if (!Cm.x3.isEmpty()) wblog(FL,"ERR %s()",FCT); }
   else {
      if (Cm.x3.SIZE.len!=3)
         wblog(FL,"ERR %s() [%s]",FCT,Cm.x3.sizeStr().data);
      if (Cm.rtype==CGR_DEFAULT) n=1; else
      if (Cm.rtype!=CGR_CTR_SCALAR) {
         wblog(FL,"ERR %s() [%s]",FCT,Cm.rtype.tostr());
      }
   }

   if (n) {
      CDATA_TQ &C = gCG.getBUF(0,0,(QSet<TQ>&)Cm.c,0);

      Cm.rtype=CGR_DEFAULT;
      Cm.cgr=&C;

      if (C.isEmpty()) {
         if (Cm.c.cstat!=CGD_REF_INIT) wblog(FL,
            "ERR %s() unexpected cstat (%s)",FCT,C.cstat.toStr().data);
         C=Cm.c;
      }
      else {
         int q=C.cstat.cmp(FL,Cm.c.cstat);
         if (q<0) {
            if (C.cstat==CGD_REF_INIT) {
               C=Cm.c;
            }
            else {
              wblog(FL,"ERR %s() got unexpected cstat sets !?"
                 "%N   %s%N<> %s%N%N   %s%N%N<> %s",FCT,
                 C.toStr().data, Cm.c.toStr().data,
                 C.cstat.toStr('V').data, Cm.c.cstat.toStr('V').data
              );
            }
         }
      }
   }

   if (CG_VERBOSE>5-n) wblog(PFL,
      "(+) XBUF[%03d] %s\n    | %s @ %s\n    | %s @ %s\n"
      "    > %s", gCX.XBUF.size(), Cm.x3Str().data,
      Cm.a.toStr().data, ica.toStr().data,
      Cm.b.toStr().data, icb.toStr().data, Cm.c.toStr().data
   );


   return 1;
};



template <class TQ>
void CStore<TQ>::getCData(
  const char *F, int L, const QType &q,
  const qset<TQ> &J1, const qset<TQ> &J2,
  wbvector < CDATA_TQ* > &S
){
  iMAP31 iJ12 = get_mp3_data(F,L, q,J1,J2);
  MAP32 &M2 = iJ12->second;
  
  S.init(M2.size()); if (S.len) {
     unsigned i=0;
     for (auto I2=M2.begin(); I2!=M2.end(); ++I2, ++i) {
        S[i]=&(I2->second);"S",FL);
     }
  }
};



template <class TQ>
const CDATA_TQ& CStore<TQ>::getCData(
  const char *F, int L,
  const QType &q, const qset<TQ> &J1, const qset<TQ> &J2,
  const qset<TQ> &J, char force
){
  if (J1.len!=J2.len || J1.len!=J.len) wblog(FL,
  "ERR qset length inconsistency (%d/%d/%d)",J1.len,J2.len,J.len);

  iMAP31 iJ12 = getCData(F,L,q,J1,J2,force);
  iMAP32 I2 = iJ12->second.find(J);

  if (I2!=iJ12->second.end()) {
     const CDATA_TQ &S=I2->second;
     if (S.q!=q) wblog(FL, "ERR QType inconsistency (%s/%s)",
         S.qStr().data,q.toStr().data);
     return S;
  }

  wblog(F,L,"WRN non-existing CData: %s [%s, %s] -> [%s]",
  q.toStr().data, J1.toStr(0).data, J2.toStr(0).data, J.toStr(0).data);

  return EMPTY_CGSTORE;
};


template <class TQ>
const CDATA_TQ& CStore<TQ>::getCData(const char *F, int L,
  const QVec &q, unsigned k,
  const qset<TQ> &J1, const qset<TQ> &J2,
  const qset<TQ> &J, char force
){
  unsigned i,l, d=q.Qlen();
  qset<TQ> j1,j2,j;

  if (k>=q.len) wblog(FL,
     "ERR %s() index out of bounds (%d/%d)",FCT,k+1,q.len); 
  if (J1.len!=J.len || J2.len!=J.len || J.len!=d) wblog(FL,
     "ERR %s() J1/J2/J inconsistency (%d,%d,%d/%d)",
      FCT,J1.len,J2.len,J.len,d);

  for (l=i=0; i<k; i++) {
     d=q[i].qlen(); l+=d;
  }; d=q[i].qlen();

  j1.init2ref(d,J1.data+l);
  j2.init2ref(d,J2.data+l);
  j .init2ref(d,J .data+l);

  return getCData(F,L,q[k],j1,j2,j,force);
};


template <class TQ>
const CDATA_TQ& CStore<TQ>::getIdentityC(
   const char *F, int L,
   const QType &t, const TQ *qs, unsigned dim, char loadRC
){
   if (t.isAbelian()) {
      if (dim!=1) wblog(FL,
         "ERR %s() abelian CGC with dim=%d !??",FCT,dim);
   }

   QSet<TQ> Q; { Q.init2(t,qs); }
   CDATA_TQ &C = getBUF(0,0,Q,loadRC);

   if (int(dim)<0) { dim=t.qdim(qs); } else 
   if (dim!=t.qdim(qs)) wblog(F_L,
      "ERR %s() got dim=%d having %s",FCT,dim,Q.toStr().data
   );

   if (C.cgd.isEmpty() && C.QSet<TQ>::isEmpty()) {
      RTD nrm=SQRT(RTD(1)/RTD(dim));

      C.t=t; C.qs=Q.qs; C.qdir.init_iout(FL,2,2);
      C.cgd.initIdentity(t,dim,nrm);
      C.cstat.init(CGD_IDENTITY);

      if (CG_VERBOSE>2) wblog(PFL,
         "[+] CBUF[%03d| Id] %s",BUF.size(),C.toStr('v').data);
      if (C!=Q) wblog(FL,
         "ERR %s() QSet inconsisteny !??\n%s",FCT,Q.toStr().data); 
      gStore.save_CData(FL,C);
   }
   else {
      if (C.t!=t || C.qs!=Q.qs || C.qdir!=Q.qdir) wblog(F_L,
         "ERR %s() CGC Id symmetry mismatch (%s <> %s)",
         FCT,C.toStr().data,Q.toStr().data);
      if (!C.isAbelian() && !C.cgd.isSMatrix(F_L,dim)) wblog(F_L,
         "ERR %s() CGC Id size mismatch (%s; %d)",FCT,C.sizeStr().data,dim);
      if (C.qdir!="+-") wblog(F_L,
         "ERR %s() invalid C.qdir=[%s]",FCT,C.qdir.toStr().data
      ); 
   }

   return C;
};


template <class TQ>
const CDATA_TQ& CStore<TQ>::getIdentityZ(
   const char *F, int L,
   const QType &t, const TQ *qs, unsigned dim, char loadRC
){
   if (t.isAbelian()) {
      if (dim!=1) wblog(FL,
         "ERR %s() abelian CGC with dim=%d !??",FCT,dim);
   }

   QSet<TQ> Q; { Q.initZ(t,qs); }
   CDATA_TQ &C = getBUF(0,0,Q,loadRC);

   if (int(dim)<0) { dim=t.qdim(qs); } else 
   if (dim!=t.qdim(qs)) wblog(F_L,
      "ERR %s() got dim=%d having %s",FCT,dim,Q.toStr().data
   );

   if (C.QSet<TQ>::isEmpty()) wblog(FL,
      "ERR %s() CData %s not yet defined",FCT,Q.toStr().data);

   if (C.t!=t || C.qs!=Q.qs || C.qdir!=Q.qdir) wblog(F_L,
      "ERR %s() CData symmetry mismatch (%s <> %s)",
      FCT,C.toStr().data,Q.toStr().data);
   if (!C.cgd.isSMatrix(F_L,dim)) wblog(F_L,
      "ERR %s() CGC Id size mismatch (%s; %d)",FCT,C.sizeStr().data,dim);
   if (C.qdir!="++") wblog(F_L,
      "ERR %s() invalid C.qdir=[%s]",FCT,C.qdir.toStr().data
   ); 

   return C;
};


template <class TQ>
void CStore<TQ>::getQfinal(
   const char *F, int L, const QType &t,
   const qset<TQ> &j1, const qset<TQ> &j2,
   wbMatrix<TQ> &jj,
   wbMatrix< const wbvector<CRef<TQ> >* > *ss
){

  iMAP31 iJ12 = get_mp3_data(F,L,t,j1,j2);
  MAP32 &M2 = iJ12->second;
  unsigned i=0, d=j1.len, j,m,m_, n=M2.size();

  if (d!=t.qlen()) wblog(FL,
     "ERR invalid qset record length (%d/%d)",d,t.qlen());
  if (!n) wblog(FL,"ERR %s() got empty map3 data !?",FCT);

  jj.init(n,d);
  if (ss) ss->init(n,1);

  for (iMAP32 I2=M2.begin(); I2!=M2.end(); ++I2, ++i) {
     const qset<TQ> &J=I2->first;
  

     if (J.len!=d) wblog(FL,
        "ERR qset record length inconsistency (%d/%d)",J.len,d);

     jj.recSetP(i,J.data);
     if (ss) {
        const wbvector< CRef<TQ> > &c3=(I2->second);
        ss->data[i]=&c3;

        if (!c3.len) wblog(FL,"ERR %s() got empty CG3s !?",FCT);

        for (j=0; j<c3.len; ++j) {
           if (!c3.data[j].cgr) wblog(FL,"ERR %s() got empty CRef "
              "(%d/%d; %s)",FCT, j+1, c3.len, t.toStr().data);

           const wbvector<unsigned> &S=c3.data[j].cgr->cgd.SIZE;
           if (S.len) {
              if (S.len<3 || S.len>4) wblog(FL,"ERR %s() "
              "unexpected size %s !?",FCT,c3.data[j].sizeStr().data);
           }

           if (j) { if (c3.data[j].cgr!=c3.data[0].cgr) wblog(FL,
              "ERR %s() got inconsistent CG3s\n"
              "(%d/%d: 0x%lX <> 0x%lX; %s)", FCT, j+1, c3.len,
              c3.data[j].cgr, c3.data[0].cgr, t.toStr().data);
           }
        }

        m=c3.data[0].cgr->getOM();

        if (m!=c3.len) wblog(FL,
           "ERR %s() OM size mismatch (%d/%d)",FCT, c3.len, m);
        if (m>1) {
           if (m!=(m_=c3[m-1].cgw.len)) wblog(FL,
              "ERR %s() OM size mismatch (%d/%d)",FCT,m_,m);
           if (c3[0].cgr->cgd.SIZE.len!=4 || m!=(m_=c3[0].Size(3)))
           wblog(FL,"ERR %s() OM size mismatch (%s; %d/%d)",
           FCT, c3[0].sizeStr().data, m_,m);
        }
     }
  }
};


template <class TQ>
void CStore<TQ>::getQfinal(
   const char *F, int L, const QVec &qvec,
   const qset<TQ> &q1, const qset<TQ> &q2, wbMatrix<TQ> &QQ,
   wbMatrix< const wbvector< CRef<TQ> >* > *SS
){
   unsigned i,l,d=0, D=qvec.Qlen();
   wbMatrix<TQ> X;
   qset<TQ> J1,J2;

   wbMatrix< const wbvector< CRef<TQ> >* > ss;

   if (D!=q1.len || D!=q2.len) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d/%d)", q1.len, q2.len, D);

   QQ.init(); if (SS) SS->init();

   for (l=i=0; i<qvec.len; ++i, l+=d) { d=qvec[i].qlen();
      J1.init(d,q1.data+l);
      J2.init(d,q2.data+l);

      getQfinal(F,L,qvec[i],J1,J2,X, SS ? &ss : NULL);
      QQ.ColKron(X);

      if (SS) SS->ColKron(ss);
   }
};


template <class TQ>
void CStore<TQ>::getQfinal(
  const char *F, int L, const QVec &qvec,
  const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, QMap<TQ> &M,
  char cgflag
){
  unsigned i,j;
  unsigned d1=Q1.dim1, d2=Q2.dim1, d=qvec.Qlen();

  wbvector<unsigned> D;
  wbMatrix< wbMatrix<TQ> > QQ;
  wbMatrix< wbMatrix< const wbvector< CRef<TQ> >* > > SS;
  wbperm P;

  qset<TQ> q1,q2;

  if (Q1.isEmpty() || Q2.isEmpty()) wblog(F,L,
     "ERR got empty input spaces (%d,%d)",Q1.isEmpty(),Q2.isEmpty());
  if (!Q1.isUnique() || !Q2.isUnique()) wblog(F,L,
      "ERR expect input to be sorted and unique");

   if (!d1 || !d2) {
      wblog(F,L,"WRN %s() got empty Q-Space (%d/%d)",d1,d2);
      M.init(); return;
   }
   if (d!=Q1.dim2 || d!=Q2.dim2) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d/%d; %s)",
       Q1.dim2, Q2.dim2, d, qvec.toStr().data);

   QQ.init(d1,d2); if (cgflag)
   SS.init(d1,d2);

   for (i=0; i<d1; i++)
   for (j=0; j<d2; j++) {
      q1.init2ref(d,Q1.rec(i));
      q2.init2ref(d,Q2.rec(j));
      getQfinal(F,L,qvec,q1,q2, QQ(i,j), cgflag ? &SS(i,j) : NULL);
   }

   M.init(FL,qvec,Q1,Q2,QQ, cgflag ? &SS : NULL);
};


template <class TQ>
int CStore<TQ>::getQfinal_zdim(
   const char *F, int L, const QVec &qvec,
   const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q,
   INDEX_T &s1, INDEX_T &s2, INDEX_T &s, unsigned *m
){
   wbMatrix< const CDATA_TQ* > SS;
   unsigned i,j,ss[4]={1,1,1,1};
   wbMatrix<TQ> QQ;

   getQfinal(F_L,qvec,q1,q2,QQ,&SS);

   if (q1.len!=q.len || QQ.dim2!=q.len) wblog(FL,
      "ERR %s() severe size inconsistency (%d,%d/%d)",FCT,q1.len,q.len);

   int id=QQ.findRec(q.data);
   if (id<0) {
      if (!F) return id;
      wblog(FL,"ERR %s() invalid total symmetry (%s): [%s; %s] => [%s]",
         FCT, qvec.toStr().data, q1.toStr().data, q2.toStr().data,
         q.toStr().data
      );
   }

   for (i=0; i<SS.dim2; ++i) {
      const wbvector<unsigned> &S=SS(id,i)->cgd.SIZE;
      if (S.len<3 || S.len>4) {
         if (S.len==1 && S[0]==1) continue; else wblog(FL,
         "ERR %s() got cgd size [%s]",FCT,S.toStrD().data);
      }
      for (j=0; j<S.len; ++j) ss[j]*=S.data[j];
   }

   for (j=0; j<4; ++j) {
      if (ss[j]<1) wblog(FL,
      "ERR %s() [%d %d %d %d] ???",FCT,ss[0],ss[1],ss[2],ss[3]); 
   }

   s1=ss[0]; s2=ss[1]; s=ss[2];
   if (m) (*m)=ss[3]; else s*=ss[3];

   return 0;
};


template <class TQ>
int CStore<TQ>::getQfinal_zdim(
   const char *F, int L, const QVec &qvec,
   const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q, 
   wbvector<INDEX_T> &s1, wbvector<INDEX_T> &s2, wbvector<INDEX_T> &s3,
   wbvector<unsigned> *m, wbMatrix< const wbvector< CRef<TQ> >* > *S3
){
   unsigned i,j,k,l,om;
   unsigned d1=Q1.dim1, d2=Q2.dim1, d3=Q.dim1, d=qvec.Qlen();

   wbvector<INDEX_T> D;
   wbvector<char> m3;
   wbindex i0,i2;
   wbMatrix< wbMatrix<TQ> > QQ;
   wbMatrix< wbMatrix< const wbvector< CRef<TQ> >* > > SS;
   wbperm P;

   qset<TQ> q1,q2;

   if (!d1 || !d2) {
      if (F) wblog(F,L,"WRN %s() got empty Q-Space (%d/%d)",d1,d2);
      s1.init(); s2.init(); s3.init();
         if (m) m->init();
         if (S3) S3->init();
      return -1;
   }
   if (!Q1.isUnique() || !Q2.isUnique() || !Q.isUnique()) wblog(F,L,
      "ERR expect input to be sorted and unique [%d,%d,%d]",
       Q1.isUnique(), Q2.isUnique(), Q.isUnique());
   if (d!=Q1.dim2 || d!=Q2.dim2 || d!=Q.dim2) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d,%d/%d; %s)",
       Q1.dim2, Q2.dim2, Q.dim2, d, qvec.toStr().data);

   QQ.initDef(d1,d2);
   SS.initDef(d1,d2);

   for (i=0; i<d1; ++i)
   for (j=0; j<d2; ++j) {
      q1.init2ref(d,Q1.rec(i));
      q2.init2ref(d,Q2.rec(j));
      getQfinal(F,L,qvec,q1,q2, QQ(i,j), &SS(i,j));
   }

   s1.init(d1).set(1);
   for (i=0; i<d1; ++i) {
       const wbMatrix< const wbvector< CRef<TQ> >* > &Si=SS(i,0);
       for (l=0; l<Si.dim2; ++l) {
          s1[i]*=Si(0,l)->el(0).Size(0);
       }
   }

   s2.init(d2).set(1);
   for (j=0; j<d2; ++j) {
       const wbMatrix< const wbvector< CRef<TQ> >* > &Sj=SS(0,j);
       for (l=0; l<Sj.dim2; ++l) {
          s2[j]*=Sj(0,l)->el(0).Size(1);
       }
   }

   if (Q.isEmpty()) return 0;

   s3.init(d3).set(1);
      if (m) m->init(d3).set(1);
      if (S3) S3->init(d3,SS(0,0).dim2);
   m3.init(d3);

   for (i=0; i<d1; ++i)
   for (j=0; j<d2; ++j) {
       const wbMatrix< const wbvector< CRef<TQ> >* > &Sij=SS(i,j);
       matchIndex(QQ(i,j),Q,i0,i2); if (!i2.len) continue;

       for (k=0; k<i0.len; ++k) { if (++m3[i2[k]]>1) continue;
          unsigned id=i2[k];
          if (S3) { S3->recSetP(id,Sij.ref(i0[k])); }

          if (Sij.dim2!=qvec.len) wblog(FL,"ERR %s() inconsistent CG3 "
             "(%d; %s)",FCT,Sij.dim2,qvec.toStr().data);

          for (l=0; l<Sij.dim2; ++l) { if (!qvec[l].isAbelian()) {
             const wbvector< CRef<TQ> > &r3=(*Sij(i0[k],l));
             om=r3.len;

             s3[id] *= r3[0].Size(2);
             if (om>1) {
                if (m)
                     { (*m)[id]*=om; }
                else { s3[id]*=om; }
             }
          }};
       }
       if (!m3.contains(0)) break;
   }

   if (m3.contains(0)) {
       matchIndex(Q1,Q,i0,i2);
       for (k=0; k<i0.len; ++k) { unsigned id=i2[k];
           if (m3[id]) continue; else m3[id]+=100;
           s3[id]=s1[i0[k]];
       }

       if (m3.contains(0)) {
       matchIndex(Q2,Q,i0,i2);
       for (k=0; k<i0.len; ++k) { unsigned id=i2[k];
           if (m3[id]) continue; else m3[id]+=101;
           s3[id]=s2[i0[k]];
       }}
   }

   if (int(i=m3.find1(0))>=0) {
      MXPut X(FL,"a");
      X.add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q").add(QQ,"QQ")
       .add(s1,"s1").add(s2,"s2").add(s3,"s3").add(m3,"m3").add(i,"i");
      if (m) X.add(*m,"m"); X.put("caller");

      if (F) wblog(F_L,
         "ERR %s() invalid output symmetry [%s]",FCT,Q.rec2Str(i).data);
      else return -1;
   }

   return 0;
};


template <class TQ> inline
void CStore<TQ>::add_CGData(
   const QType &q, const qset<TQ> &J1,const qset<TQ> &J2
){
   switch (q.type) {

     case QTYPE_U1  : add_CGData_U1  (q,J1,J2); break;
     case QTYPE_ZN  : add_CGData_ZN  (q,J1,J2); break;
     case QTYPE_P   : add_CGData_P   (q,J1,J2); break;

     case QTYPE_SUN :
     case QTYPE_SpN : {
        int i=gStore.load_Std3(FL,q,J1,J2);
        if (i<0) wblog(FL,"ERR %s() failed to load/generate\n"
           "CGC data for %s [%s; %s] (i=%d)",
           FCT,q.toStr().data, J1.toStr().data, J2.toStr().data,i);
        }
        break;

     default:
        wblog(FL,"ERR symmetry label %s (%d) not implemented yet",
        q.toStr().data,q.type);
   }
};


template <class TQ>
void CStore<TQ>::add_CGData_U1(
   const QType &q, const qset<TQ> &J1, const qset<TQ> &J2
){
   if (q.type!=QTYPE_U1 || J1.len!=J2.len || J1.len!=1) {
      wblog(FL,"ERR invalid %s qset record length (%d,%d)",
      q.toStr().data,J1.len,J2.len);
   }

   unsigned ix=0, mx=(J1[0]!=J2[0] ? 2 : 1);
   qset<TQ> J12(J1,J2), J21(J2,J1), J(1);

   if (J1[0]>J2[0]) {
      TQ *q=J12.data; J12.data=J21.data; J21.data=q;
   }
   J[0]=J1[0]+J2[0];


   if (CG_VERBOSE>1) wblog(PFL,"[+] map3[%02d] "
      "%-5s %s [%2g %2g; %2g ]", map3[q].size(), q.toStr().data,
      mx>1 ? "X":"=", double(J12[0]), double(J12[1]), double(J[0])
   );

   CDATA_TQ S; 


   S.init3(q,1U,1U,1U);
   S.cgd.init(); S.cstat.init(CGD_ABELIAN);

   S.qs[0]=J12[0];
   S.qs[1]=J12[1]; S.qs[2]=J[0];

   CDATA_TQ &Sb=getBUF(0,0,(QSet<TQ>&)S, 0);
   S.save2(Sb);

   for (; ix<mx; ++ix) {
      wbvector< CRef<TQ> > &cg3=map3[q][ix==0 ? J12 : J21][J];

      if (cg3.len!=1) {
         if (!cg3.len) cg3.init(1); else wblog(FL,
         "ERR %s() invalid abelian cg3 (len=%d)",FCT,cg3.len);
      }
      CRef<TQ> &R=cg3[0];

      if (R.stat() ^ (R.cgr!=NULL)) wblog(FL,
         "ERR %s() CRef status mismatch: %s", FCT, R.statStr().data);

      if (R.cgr) {
         if (!R.cgr->sameAs(Sb,'l')) wblog(FL,
            "ERR %s() inconsistent U1 CGC\n   %s\n<> %s",
            FCT,R.cgr->toStr().data, Sb.toStr().data);
         if (R.cgw.len!=1 || R.conj || (!ix && R.cgp.len) || (ix && 
            (R.cgp.len!=3 || R.cgp[0]!=1 || R.cgp[1]!=0))) wblog(FL,
            "ERR %s() invalid scalar CRef (%d)\n%s",
            FCT, ix+1, R.toStr().data
         );
      }
      else {
         if (R.cgw.len || R.cgp.len || R.conj) wblog(FL,"ERR %s() "
            "got unexpected empty CRef\n%s",FCT,R.toStr().data);
         R.cgr=&Sb;
         R.cgw.init(1); R.cgw[0]=1; if (ix) {
         R.cgp.initStr("2,1,3"); }

      }
   }
};


template <class TQ>
void CStore<TQ>::add_CGData_ZN(
   const QType &q, const qset<TQ> &J1, const qset<TQ> &J2
){
   if (q.type!=QTYPE_ZN || J1.len!=J2.len || J1.len!=1) {
      wblog(FL,"ERR invalid %s qset record length (%d,%d)",
      q.toStr().data,J1.len,J2.len);
   }

   qset<TQ> J12(J1,J2), J(1); J[0]=int(J1[0]+J2[0])%q.sub;

   if ( TQ(int(J1[0])%q.sub)!=J1[0] || TQ(int(J2[0])%q.sub)!=J2[0]) {
      for (unsigned i=0; i<2; ++i) J12[i]=(int(J12[i])%q.sub);
      wblog(FL,"WRN %s() invalid Z%d labels (%g,%g; %s => %g,%g)",
      FCT, q.sub, double(J1[0]), double(J2[0]), q.toStr().data,
      double(J12[0]), double(J12[1]));
   }


   if (CG_VERBOSE>1) wblog(PFL,
   "[+] map3[%02d]  %-7s  [ %3g %3g; %3g ]", map3[q].size(),
        q.toStr().data, double(J12[0]), double(J12[1]), double(J[0])
   );

   wbvector< CRef<TQ> > &cg3=map3[q][J12][J];

   if (cg3.len!=1) wblog(FL,
      "ERR %s() invalid abelian cg3 (len=%d)",FCT,cg3.len);
   CRef<TQ> &R=cg3[0];

   if (R.stat() ^ (R.cgr!=NULL)) wblog(FL,
      "ERR %s() CRef status mismatch: %s", FCT, R.statStr().data);

   if (R.cgr) {
      const CDATA_TQ &S=(*R.cgr);
      if (S.t!=q || S.cgd.D.len || S.cstat!=CGD_ABELIAN || S.qs.len!=3 ||
          S.qs[0]!=J1[0] || S.qs[1]!=J1[1] || S.qs[2]!=J1[2])
      wblog(FL,"ERR %s() inconsistent U1 CGC",FCT);
   }
   else {
      CDATA_TQ S; 
      S.init3(q,1U,1U,1U);
      S.cgd.init(); S.cstat.init(CGD_ABELIAN);
      S.qs[0]=J1[0];
      S.qs[1]=J2[0];
      S.qs[2]=J [0];

      CDATA_TQ &Sb=getBUF(FL,(QSet<TQ>&)S);
      S.save2(Sb);
   }
};


template <class TQ>
void CStore<TQ>::add_CGData_P(
   const QType &q, const qset<TQ> &J1, const qset<TQ> &J2
){
   if (q.type!=QTYPE_P || J1.len!=J2.len || J1.len!=1) {
      wblog(FL,"ERR invalid %s qset record length (%d,%d)",
      q.toStr().data,J1.len,J2.len);
   }
   if ( fabs(double(J1[0]))!=1 || fabs(double(J2[0]))!=1 ) {
      wblog(FL,"ERR %s() invalid %s labels (%g,%g)",
      FCT,q.toStr().data,double(J1[0]),double(J2[0])); 
   }

   qset<TQ> J12(J1,J2), J(1); J[0]=J1[0]*J2[0];


   if (CG_VERBOSE>1) wblog(PFL,
   "[+] map3[%02d]  %-7s  [ %3g %3g; %3g ]", map3[q].size(),
        q.toStr().data, double(J1[0]), double(J2[0]), double(J[0])
   );

   wbvector< CRef<TQ> > &cg3=map3[q][J12][J];

   if (cg3.len!=1) wblog(FL,
      "ERR %s() invalid abelian cg3 (len=%d)",FCT,cg3.len);
   CRef<TQ> &R=cg3[0];

   if (R.stat() ^ (R.cgr!=NULL)) wblog(FL,
      "ERR %s() CRef status mismatch: %s", FCT, R.statStr().data);

   if (R.cgr) {
      const CDATA_TQ &S=(*R.cgr);
      if (S.t!=q || S.cgd.D.len || S.cstat!=CGD_ABELIAN || S.qs.len!=3 || 
          S.qs[0]!=J1[0] || S.qs[1]!=J1[1] || S.qs[2]!=J1[2])
      wblog(FL,"ERR %s() inconsistent U1 CGC",FCT);
   }
   else {
      CDATA_TQ S; 
      S.init3(q,1U,1U,1U);
      S.cgd.init(); S.cstat.init(CGD_ABELIAN);
      S.qs[0]=J1[0];
      S.qs[1]=J2[0];
      S.qs[2]=J [0];

      CDATA_TQ &Sb=getBUF(FL,(QSet<TQ>&)S);
      S.save2(Sb);
   }
};


template <class TQ>
void CStore<TQ>::add2BUF(
   const char *F, int L, const CDATA_TQ &S,
   char nflag
){
   unsigned r=S.rank(F_L);
   if (!S.t.validType() || r!=3)
      wblog(FL,"ERR %s() got invalid rank-%d CData [%s]",
      FCT, r, S.qStr().data);
   if (S.cgd.SIZE.len<3 || S.cgd.SIZE.len>4)
      wblog(FL,"ERR %s() got invalid CData [%s: %s @ r=%d]",
      FCT, S.qStr().data, S.cgd.sizeStr().data, r
   );

   CDATA_TQ &B = getBUF(FL,QSet<TQ>(S));
   char gotB=(B.isEmpty()? 0 : B.isRefInit()? -1 : 1);

   if (gotB<=0) {
      if (gotB) {
         if (S.cgd.SIZE!=B.cgd.SIZE || (QSet<TQ>&)S!=(QSet<TQ>&)B)
         wblog(FL,"ERR %s() inconsistent size-ref data\n"
         "%s <> %s",FCT,S.toStr().data, B.toStr().data); 
      }
      if (CG_VERBOSE>2) wblog(PFL,
         "[+] CBUF[%02d] %s %s", BUF.size(), S.toStr().data,
         gotB==0 ? "": " [replaces Sref]"
      );
      B=S; 

      gStore.save_CData(FL,S);
   }
   else { str[0]=0;

      if (nflag)
         strcat(str,"already got existing entry"); else
      if (B.QSet<TQ>::operator!=(S))
         strcat(str,"inconsistent pre-existing data");
      else 
         strcat(str,"[TST double check this]");

      if (str[0]) {
         MXPut(FL,"q").add(S,"A").add(B,"A_").add(nflag,"nflag")
         .add(wbstring(str),"str");

         wblog(FL,"ERR %s() %s\nB: %s\nS: %s",
            FCT,str,B.toStr().data, S.toStr().data);
      }
   }
};


template <class TQ>
CDATA_TQ& CStore<TQ>::getBUF(
   const char *F, int L, const QSet<TQ> &Q, char loadRC 
){
   if (Q.isEmpty()) wblog(F_L,
      "ERR %s() got empty QSet !?",FCT);
   if (!Q.isSorted()) wblog(F_L,"ERR %s() "
      "CStore requires sorted QSet\n%s",FCT,Q.toStr().data
   );

   CDATA_TQ *C=NULL; char gotnew=0;

   auto it=BUF.find(Q);
   if (it!=BUF.end())
        { C=&(it->second); }
   else { C=&(BUF[Q]); gotnew=1; }

   if (Q.t.isAbelian()) { C->initAbelian(Q); return *C; }


   if (C->isEmpty() || C->cstat.t==CGD_REF_INIT) { if (loadRC) {

      int q=gStore.load_CData(0,0,Q,*C,'b');


      if (q>0) { if (CG_VERBOSE) { char tag[6]="";
         if (q==1) { if (CG_VERBOSE>4) strcpy(tag,"(+)"); }
         else if (q==3) {
            if ((C->qdir.len>2 && CG_VERBOSE>1) || CG_VERBOSE>4)
            strcpy(tag,"{+}"); }
         else if (q==4) {
            if ((C->qdir.len>2 && CG_VERBOSE>1) || CG_VERBOSE>4)
            strcpy(tag,"{-}"); }
         else if (q==5) {
            if (CG_VERBOSE>4) strcpy(tag,"(d)"); }
         else if (q==6) {
            if (CG_VERBOSE>5) strcpy(tag," ok"); }
         else wblog(FL,"ERR %s() load_CData() returned q=%d !?",FCT,q);

         if (tag[0]) wblog(PFL,
            "%3s CBUF[%03ld] %s",tag,BUF.size(),C->toStr().data
         );
      }}
      else { bool isJ1=Q.isJ1Symbol();
      if (Q.isStd3() || (isJ1 && Q.t.qlen()<2)) {
         unsigned m=Q.t.qlen();
         if (Q.qs.len != m*Q.qdir.len) wblog(FL,
            "ERR %s() invalid QSet %s",FCT,Q.toStr().data);
         qset<TQ> q1(m,Q.qs.data,'r'), q2(m,Q.qs.data+m,'r');
         gRG.getR(FL,Q.t,q1.data);
         gRG.getR(FL,Q.t,q2.data);

         gRG.buf[Q.t].getTensorProdReps_gen(q1,q2,'t');
      }
      else if (Q.qdir.len==2) {
         if (Q.isScalar()) {
            qset<TQ> J(Q.t.qlen(),Q.qs.data,'r');
            if (gStore.load_RSet(0,1,Q.t,J)>0) {
               gCG.getIdentityC(FL,Q.t,Q.qs.data,-1,0);
            }
         }
         else if (!isJ1) wblog(FL,"ERR %s() "
            "got rank-%d QSet !?\n%s",FCT,Q.qdir.len,Q.toStr().data
         );
      }}
   }}
   else if (loadRC=='L') { if (C->checkOM()) {
      wbstring file;
      int q=gStore.get_directory(F_L,file,Q,"cgd");
      if (q<=0) wblog(FL,
         "ERR missing CBUF %s (i=%d)\n%s",Q.toStr().data,q,file.data);

      Wb::matFile f(FL,file.data,"r");
      if (!f.mfp) wblog(FL,
         "ERR %s() failed to open file\n%s",FCT,file.data);

      mxArray *a=matGetVariable(f.mfp,"CRef");
      if (!a) wblog(FL,"ERR %s() "
         "missing variable 'cid' in file\n%s",FCT,file.data);
      
      CDATA_TQ Cr; Cr.init(FL,a,0,0,0);
      mxDestroyArray(a);

      if (Cr.cstat!=C->cstat) wblog(FL,
         "ERR %s() cstat mismatch with RCStore\n"
         "=> double check%N%N    %s%N -> %s",
         FCT,Cr.cstat.toStr('V').data,C->cstat.toStr('V').data);


   }}

   if (C->isEmpty()) {
      if (F || (CG_VERBOSE>1 && !gotnew)) { sprintf(str,
       " gCG.BUF[%03ld] %s (new/empty)",BUF.size(),Q.toStr().data);
         if (F) wblog(PF_L,"ERR %s",str);
         else wblog(PF_L,"TST %s",str);
      }
   }

   return *C;
};


template <class TQ>
void CStore<TQ>::testCG(const char *F, int L) {

   static TQ ncalls=0;

   QType q; qset<TQ> J12(2), J(1); TQ i;

   Info(F_L);

   for (i=0; i<3; i++) {
   q="A";   J12[0]=1; J12[1]=2+i+ncalls; J[0]=J12.sum();
   CDATA_TQ &S=gCG.map3[q][J12][J]; } ; ncalls+=i;

   Info(F_L);
   q="SU2"; J12[0]=2; J12[1]=1+ncalls; J[0]=J12.sum();
   { CDATA_TQ &S=gCG.map3[q][J12][J]; }

   Info(F_L);
   q="A";   J12[0]=3; J12[1]=0+ncalls; J[0]=J12.sum();
   { CDATA_TQ &S=gCG.map3[q][J12][J]; }

   Info(F_L);
   q="SU2"; J12[0]=4; J12[1]=1+ncalls; J[0]=J12.sum();
   { CDATA_TQ &S=gCG.map3[q][J12][J]; }

   Info(F_L);

   if (++ncalls>3) wblog(F_L,"ERR");
};

template <class TQ>
void CStore<TQ>::Info(const char *F, int L){
   unsigned m1=0,m2=0;

   for (iMAP30 I0=map3.begin(); I0!=map3.end(); ++I0) {
      const MAP31 &M1 = I0->second;

      for (ciMAP31 I1=M1.begin(); I1!=M1.end(); ++I1) {
         const MAP32 &M2 = I1->second;

         m1+=M1.size();
         m2+=M2.size();

      }
   }
   wblog(F,L,"<i> CGstore::%s has %d/%d/%d elements",FCT,map3.size(),m1,m2);
};


template <class TQ>
void CStore<TQ>::printStatus(const char *F, int L){

   unsigned i,k,i1,i2,n,n1,n2;

   wblog(F,L,"<i> CG store inventory");

   for (iMAP30 I0=map3.begin(); I0!=map3.end(); ++I0) {
      const QType &q = I0->first;
      const MAP31 &M1 = I0->second; i1=0;

      printf("\n"); banner(75,"="); i=M1.size();
      printf("  %s symmetries (%d %s):\n", q.toStr().data,
      i, i!=1 ? "entries" : "entry");

      if (q.isAbelian()) printf("   * "
         "(Q1,Q2) => Q (no z-labels), CG coefficients\n");
      else if (q.isSU2()) printf("   * "
         "(S1,S2) => S, [Sz1, Sz2, Sz=Sz1+Sz2], CG coefficients\n");

      banner(75,"=");

      for (ciMAP31 I1=M1.begin(); I1!=M1.end(); ++I1) {
         qset<double> J12; J12.initT(I1->first);
         const MAP32 &M2 = I1->second;

         n=0; n1=q.qdim(J12.data); n2=q.qdim(J12.data+J12.len/2);
         for (ciMAP32 I2=M2.begin(); I2!=M2.end(); ++I2) {
            n+=unsigned(round(I2->second.cgd.norm2()));
         }

         printf("\n%6d.  [ %s ]", ++i1, J12.toStr().data);
         if (M2.size()>1 || n1*n2!=n) {
            printf(" =>  %d symmetry sectors (%dx%d = %d states%s)\n",
            M2.size(),n1,n2,n, n1*n2!=n ? " ERR !??" : " - ok");
         }

         i2=1;
         for (ciMAP32 I2=M2.begin(); I2!=M2.end(); ++I2, ++i2) {
            const qset<TQ> &J = I2->first;
            const CDATA_TQ &S = I2->second;
            unsigned m=S.getOM(FL);

            printf("%s", M2.size()>1 ? "\n        " : " -> ");

            if (J==S.Q) { printf("%2d. ", i2); }
            else printf("WRN [%s] ?? ", J.toStr().data);

            printf("[%s](%d) ",S.QStrS().data, S.rank(F_L));

            for (k=0; k<m; k++){
               double dbl=S.cgd.norm2(k), e=std::fabs(dbl-round(dbl));
               if (e<=1E-12) dbl=round(dbl);
               printf("  # %g state%s",dbl,dbl!=1 ? "s" : "");
               if (e>1E-12) printf(" @ %.3g (WRN)",e);

               if (q.isSU2()) { dbl=round(dbl);
                  if (dbl!=J[0]+1) { printf("\n");
                     wblog(FL,"ERR %s dimensional inconsistency (%g/%g: %d %g)",
                     q.toStr().data, dbl, double(J[0]+1), J.len, double(J[0]));
                  }
               }
            }

         }
         if (M2.size()>1) printf("\n");
      }
      printf("\n");
   }
   printf("\n");
};


template <class TQ>
unsigned CStore<TQ>::add3(
   const QType &q, const mxArray* S, unsigned k
){
   if (!S || !mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid CStore.%s.std (%s <> cell)",
      FCT,q.toStr('t').data, S ? mxGetClassName(S) : "(null)");

   unsigned n=mxGetNumberOfElements(S), m=q.qlen();
   qset<TQ> J, J12;

   QDir qd;

   int itype = mxGetFieldNumber(S,"type"),
       iqset = mxGetFieldNumber(S,"qset"),
       iqdir = mxGetFieldNumber(S,"qdir");

   if (int(k)>=0) {
      if (k>=n) wblog(FL,"ERR CStore.%s.std "
         "index out of bounds (%d/%d)",q.toStr('t').data, k+1,n);
      n=k+1;
   } else k=0;

   if (itype<0 || iqset<0 || iqdir<0) wblog(FL,
      "ERR invalid CStore.%s.std (%d,%d,%d)",
      q.toStr('t').data, itype,iqset,iqdir);
   else {
      QSet<TQ> Q;
      Q.t.   init(FL,mxGetFieldByNumber(S,k,itype));
      Q.qs.  init(FL,mxGetFieldByNumber(S,k,iqset),0,'!');
      Q.qdir.init(FL,mxGetFieldByNumber(S,k,iqdir));

      if (Q.qs.len%3 || Q.qs.len/3!=m || Q.t!=q) wblog(FL,
         "ERR %s() invalid %s CRef data (%d/3*%d; %s)",
         FCT, q.toStr('t').data, Q.qs.len, m, Q.t.toStr('t').data);
      J12.init(2*m,Q.qs.data);
      qd=Q.qdir;
   }

   MAP32 &M3 = map3[q][J12];
   CRef<TQ> X;

   char xflag=(M3.size() ? 1 : 0);

   for (; k<n; ++k) {
      X.init(FL,S,k);
      const QSet<TQ> &Q=(const QSet<TQ>&)(*X.cgr);
      J.init(m,Q.qs.data+2*m);

      if (k) {
         if (Q.t!=q) wblog(FL,"ERR %s() qtype mismatch (%s <> %s)",
            FCT, Q.t.toStr('t').data, q.toStr('t').data); 
         if (Q.qdir!=qd) wblog(FL,"ERR %s() qdir mismatch (%s <> %s)",
            FCT, Q.qdir.toStr().data, qd.toStr().data
         );
      }

      wbvector< CRef<TQ> > &R=M3[J];
      unsigned l=0; double x;

      for (; l<R.len; ++l) {
         if (!X.cgr || X.cgr!=R[l].cgr) wblog(FL,
            "ERR %s() missing cgr reference space (%d: 0x%lX, 0x%lX)",
            FCT,l+1,X.cgr,R[l].cgr);
         if (X==R[l]) break;
         else {
            x=X.scalarProd(FL,R[l]);
            if (x>CG_EPS2) wblog(FL,
               "ERR %s() got non-orthogonal CG3 @ %.3g",FCT,x
            );
         }
      }

      if (xflag) {
         if (l==R.len) {
            if (!R.len) wblog(FL,"ERR %s() got empty CRef !??",FCT);
            else wblog(FL,
               "ERR %s() got CRef mismatch\n(%s <> %s)\n   %s\n-> %s)",
               FCT, X.toStr().data, R[0].toStr().data,
               X.cgr->cstat.toStr().data, R[0].cgr->cstat.toStr().data
            );
         }
      }
      else {
         if (l<R.len) wblog(FL,"ERR %s() "
            "got match with existing CRef (%d/%d) !?",FCT,l+1,R.len);
         R.Resize(l+1,&X);
      }
   }

   return (xflag ? 0 : n);
};


template <class TQ>
mxArray* CStore<TQ>::toMx(const QType &q_) const {

   const char *field0[]={"std"};

   mxArray *c,*S, *C=mxCreateStructMatrix(1,1,0,NULL);

   CDATA_TQ X_;

   unsigned i=0,j,j12;
   mxArray *a;

   for (auto I0=map3.begin(); I0!=map3.end(); ++I0) {
      const QType &q = I0->first;
      const MAP31 &M12 = I0->second;
   if (!q_.isKnown() || q_==q) {

      S=mxCreateStructMatrix(1,1,1,field0);
      c=map3_CreateStructMatrix(1,M12.size()); j12=0;

      for (auto I12=M12.begin(); I12!=M12.end(); ++I12, ++j12) {
         map3_add2MxStruct(c,j12, q, I12->first);
      }; mxSetFieldByNumber(S,0,0,c);

      map<unsigned,
         map <qset<TQ>,
            const CDATA_TQ* >
      > M;

      for (auto it=BUF.begin(); it!=BUF.end(); ++it) {
         if (it->first.t==I0->first) {
            const qset<TQ> &qs=it->first.qs;
            if (qs!=it->second.qs) wblog(FL,
               "ERR %s() got inconsistent key [%s] <> [%s]",
               FCT,qs.toStr().data,it->second.QStrS().data
            );
            M[it->second.rank(FL)][qs] = &(it->second);
         }
      }

      for (auto i1=M.begin(); i1!=M.end(); ++i1) { j=0;
         a=X_.mxCreateStruct(1,i1->second.size());
         for (auto i2=i1->second.begin(); i2!=i1->second.end(); ++i2) {
            i2->second->add2MxStruct(a,j++);
         }
         sprintf(str,"rank%d",i1->first);
         mxAddField2Scalar(FL,S,str,a);
      }
      mxAddField2Scalar(FL,C,q.toStr('t').data,S);
   }}

   wbvector<unsigned> nc(BUF.bucket_count());
   for (; i<nc.len; ++i) nc[i]=BUF.bucket_size(i);

   mxAddField2Scalar(FL,C,"QHash", MXPut(0,0)
      .add(BUF.max_load_factor(),"max_load_factor")
      .add(BUF.load_factor(),"load_factor")
      .add(nc,"filling_buckets")
      .add(BUF.size(),"nr_entries")
      .add(nc.max(0),"max_collisions")
   .toMx());

   return C;
};


mxArray* map3_CreateStructMatrix(unsigned m, unsigned n) {
   const char *fields[]={"J12","J","omult","cgr"};
   return mxCreateStructMatrix(m,n,4,fields);
};


template <class TQ>
void map3_add2MxStruct(
   mxArray *S, unsigned k, const QType &q, const qset<TQ> &J12) {

   const MAP32 &M3=gCG.map3[q][J12];
   wbMatrix<TQ> JJ(M3.size(),J12.len/2);
   wbvector<unsigned> OM(M3.size());

   unsigned i=0, j=0, l=0;
   CRef<TQ> R_;

   for (auto I3=M3.begin(); I3!=M3.end(); ++I3, ++j) {
      const qset<TQ> &J=I3->first;
      const wbvector< CRef<TQ> > &c3=I3->second;
      if (!c3.len) wblog(FL,"ERR %s() got empty cg3 !?",FCT);

      l+=(OM[j]=c3.len);
      JJ.recSet(j,(wbvector<TQ>&)J);

      if (c3[0].got3(FL,q,J12,J)!=0) wblog(FL,
         "ERR %s() qset mismatch (%s: %s %s <> %s: %s; %d)",
         FCT, q.toStr().data, J12.toStr().data, J.toStr().data,
         c3[0].qStr().data, c3[0].QStr().data,
         c3[0].got3(FL,q,J12,J)
      );

      for (i=1; i<c3.len; ++i) {
         if ((QSet<TQ>&)(*c3[i].cgr)!=(QSet<TQ>&)(*c3[0].cgr))
         wblog(FL,"ERR %s() severe cg3 inconsistency\n"
            "%d/%d: %s <> %s", FCT,i+1,c3.len,
            ((QSet<TQ>&)c3[i]).toStr().data,
            ((QSet<TQ>&)c3[0]).toStr().data
         );
      }
   }

   mxArray *a=R_.mxCreateStruct(l,1); l=j=0;

   for (auto I3=M3.begin(); I3!=M3.end(); ++I3, ++j) {
      const wbvector< CRef<TQ> > &c3=I3->second;

      for (i=0; i<c3.len; ++i, ++l) {
         c3[i].add2MxStruct(a,l);
      }
   }

   mxSetFieldByNumber(S,k,0, J12.toMx());
   mxSetFieldByNumber(S,k,1, JJ.toMx());
   mxSetFieldByNumber(S,k,2, OM.toMx());
   mxSetFieldByNumber(S,k,3, a);
};


template <class TQ, class TD>
mxArray* RStore<TQ,TD>::toMx(const QType &q_) const {

   mxArray *c, *C=mxCreateStructMatrix(1,1,0,NULL);

   for (ciMAP it=buf.begin(); it!=buf.end(); ++it) {
      if (!q_.isKnown() || q_==it->first) {
         const QType &q = it->first;
         const genRG_struct<TQ,TD> &B = it->second; 

         if (!B.q.isKnown()) wblog(FL,
            "ERR %s() got sym=%s / %s @ size=%d !?",
            FCT, B.q.toStr().data, it->first.toStr().data, B.RSet.size()
         ); 

         c=B.toMx();
         mxAddField2Scalar(FL,C,q.toStr('t').data,c);
      }
   }
   return C;
};


template <class TQ, class TD>
unsigned RStore<TQ,TD>::add(
   const QType &q, const mxArray* S, unsigned k
){
   static int itype=-1, iJ=-1, iZ=-1, iSp=-1, iSz=-1;
   const mxArray *ap, *az;

   if (!S || !mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid RStore.%s.store",FCT,q.toStr('t').data);
   if (k>=mxGetNumberOfElements(S)) wblog(FL,
      "ERR RStore.%s.store index out of bounds (%d/%d)",
      q.toStr('t').data, k+1, mxGetNumberOfElements(S));

   if (k==0) {
      int id[5]={ mxGetFieldNumber(S,"type"),
         mxGetFieldNumber(S,"J"),  mxGetFieldNumber(S,"Z"),
         mxGetFieldNumber(S,"Sp"), mxGetFieldNumber(S,"Sz")
      };
      if (itype<0) {
         if ((itype = id[0])<0 || (iJ = id[1])<0 || (iZ = id[2])<0 ||
             (iSp = id[3])<0 || (iSz = id[4])<0) wblog(FL,
            "ERR invalid RStore.%s.store (%d,%d,%d,%d,%d)",
             q.toStr('t').data, itype,iJ,iZ,iSp,iSz
         );
      }
      else {
         if (itype!=id[0] || iJ!=id[1] || iZ!=id[2] ||
             iSp!=id[3] || iSz!=id[4]) wblog(FL,
            "ERR invalid RStore.%s.store (%d/%d,%d/%d,%d/%d,%d/%d,%d/%d)",
             q.toStr('t').data, itype, id[0], iJ, id[1], iZ, id[2],
             iSp, id[3], iSz, id[4]
         );
      }
   }

   qset<TQ> qs(FL,mxGetFieldByNumber(S,k,iJ),0,'!');
   QType qk(FL,mxGetFieldByNumber(S,k,itype));

   if (qs.len!=q.qlen() || qk!=q) wblog(FL,
      "ERR %s() invalid %s RSet (%d/%d; %s)",
      FCT, q.toStr('t').data, qs.len, q.qlen(), qk.toStr('t').data
   );

   char xflag=0;
   genRG_base<TQ,TD> X, &R = buf[q].RSet[qs];

   if (CG_VERBOSE>4) {
      wblog(PFL,"(+) RBUF[%03d] %s [%s]",
      buf[q].RSet.size(), q.toStr().data, qs.toStr().data);
   }

   if (!R.isEmpty()) {
      xflag=1; R.save2(X);
   }

   R.q=q;
   R.J=qs;

   R.Z.init(FL, mxGetFieldByNumber(S,k,iZ));

   if (!(ap=mxGetFieldByNumber(S,k,iSp)) || mxGetNumberOfElements(ap)!=qs.len)
      wblog(FL,"ERR %s() invalid %s RSet(%d)->Sp",FCT,q.toStr('t').data,k);
   if (!(az=mxGetFieldByNumber(S,k,iSz)) || mxGetNumberOfElements(az)!=qs.len)
      wblog(FL,"ERR %s() invalid %s RSet(%d)->Sz",FCT,q.toStr('t').data,k);

   R.Sp.init(qs.len);
   R.Sz.init(qs.len);

   for (unsigned i=0; i<qs.len; ++i) {
      R.Sp[i].init(FL,ap,i);
      R.Sz[i].init(FL,az,i);
   }

   if (xflag) {
      if (R.q!=X.q || R.J!=X.J || R.Z!=X.Z || R.Sp!=X.Sp || R.Sz!=X.Sz) {
         MXPut(FL,"I").add(R,"R").add(X,"X");
         wblog(FL,"TST %s() q : %s",FCT,R.q ==X.q ? "same" : "different !??");
         wblog(FL,"TST %s() J : %s",FCT,R.J ==X.J ? "same" : "different !??");
         wblog(FL,"TST %s() Z : %s",FCT,R.Z ==X.Z ? "same" : "different !??");
         wblog(FL,"TST %s() Sp: %s",FCT,R.Sp==X.Sp? "same" : "different !??"); 
         wblog(FL,"TST %s() Sz: %s",FCT,R.Sz==X.Sz? "same" : "different !??"); 
      }
   }
   else {
      gCG.getIdentityC(FL,q,qs.data);
   }

   return (xflag ? 0 : 1);
};


template <class TQ, class TD>
void RStore<TQ,TD>::Info(const char *F, int L) const {

   wblog(F_L,"<i> RStore::%s() %s",
      FCT, buf.size() ? "...":"is empty"); 

   for (ciMAP it=buf.begin(); it!=buf.end(); ++it) {
      const QType &q = it->first;
      const genRG_struct<TQ,TD> &B = it->second; 
      unsigned i=0, l=0, n=B.RSet.size();

      printf("%12s ",q.toStr());

      for (crMAP I=B.RSet.begin(); I!=B.RSet.end(); ++I, ++i) {
         const qset<TQ> &q=I->first;
         l+=printf("%s[%s]",i?", ":"", q.toStr().data);
         if (l>50) {
            if (i+1<n) { printf(" ..."); }
            break;
         }
      }
      printf(" (%d)\n",n);
   }
};


template <class TQ, class TD>
genRG_base<TQ,TD>& RStore<TQ,TD>::getR(
   const char *F, int L, const QType &t, const TQ *qs
 ) const {

   qset<TQ> J;
   J.wbvector<TQ>::init2ref(t.qlen(),qs);

   genRG_base<TQ,TD> &R=gRG.buf[t].RSet[J];

   if (R.J.len && R.Z.dim2) return R;

   if (gRG.buf[t].RSet.size()==1) {
      gRG.buf[t].q=t; gRG.buf[t].initBase(FL);
      if (R.J.len && R.Z.dim2) return R;
   }

   if (!R.isEmpty()) wblog(FL,
      "ERR %s() got partially empty RSet !?",FCT);

   int q=gStore.load_RSet(0,1,t,J);
   if (q>0) {
      if (CG_VERBOSE>4 && L) wblog(PFL,
         "(+) RBUF[%03d] %s [%s]",gRG.buf[t].RSet.size(),
         t.toStr().data,J.toStr().data
      );
   }
   else wblog(FL,
      "ERR %s() %s irep (%s) not yet in RCStore (e=%d)",
      FCT, t.toStr().data, J.toStr().data, q
   );

   if (!R.J.len || !R.Z.dim2) wblog(F_L,
      "ERR %s() %s irep (%s) not yet generated (%dx%d)",
      FCT, t.toStr().data, J.toStr().data, R.Z.dim1, R.Z.dim2
   );

   return R;
};



int RCStore::get_directory(
   const char *F, int L, char *dir, unsigned n,
   const QType &t, const char *tag, const char *sub,
   char mflag
 ) const {

   unsigned l, m=0;

   if (!t.validType()) wblog(FL,
      "ERR %s() got invalid symmetry (%s)",FCT,t.toStr().data);
   if (!tag || !tag[0]) wblog(F_L,
      "ERR %s() got %s tag",FCT,tag?"empty":"null");
   if (!dir) wblog(FL,"ERR %s() got null string",FCT);

   l=snprintf(dir,n,
      "%s/%s/%s", root.data, t.toStr('t').data, tag);
   if (l+32>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n-32);

   if (!Wb::fexist(dir,'d')) { if (!mflag) { return -1; }

      int lflag=0;
      int e=0; mode_t p=(S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      m=strlen(tag); l-=(m+1); dir[l]=0;"/<tag>"

      if (!Wb::fexist(dir,'d')) {
         if ((e=mkdir(dir,p))==0) { ++lflag;
            printf("\n"
              "   NB! creating RCStore directories for symmetry %s\n"
              "   ==> %s\n",
            t.toStr().data, repHome(dir).data);
         }
         else wblog(F_L,
            "ERR failed to create RCStore directory (e=%d)\n%s",e,dir);
      }; dir[l++]='/';

      strcpy(dir+l,"RStore");
      if (!Wb::fexist(dir,'d')) {
         if ((e=mkdir(dir,p))==0) {
            if (++lflag==1) printf("\n");
            printf(  "    -> %s\n",repHome(dir).data);
         }
         else wblog(F_L,
            "ERR failed to create RCStore directory (e=%d)\n%s",e,dir);
      }

      strcpy(dir+l,"CStore");
      if (!Wb::fexist(dir,'d')) {
         if ((e=mkdir(dir,p))==0) {
            if (++lflag==1) printf("\n");
            printf(  "    -> %s\n",repHome(dir).data);
         }
         else wblog(F_L,
            "ERR failed to create RCStore directory (e=%d)\n%s",e,dir);
      }

      strcpy(dir+l,"XStore");
      if (!Wb::fexist(dir,'d')) {
         if ((e=mkdir(dir,p))==0) {
            if (++lflag==1) printf("\n");
            printf(  "    -> %s\n",repHome(dir).data);
         }
         else wblog(F_L,
            "ERR failed to create RCStore directory (e=%d)\n%s",e,dir);
      }

      if (lflag) printf("\n");


      strcpy(dir+l,tag);
      if (!Wb::fexist(dir,'d')) {
         wblog(FL,"ERR %s() invalid R+C store directory\n%s",FCT,dir);
      }; l+=m;
   }

   if (sub && (m=strlen(sub))>0) {
      if (l+m+1>n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l+m+1,n);
      dir[l++]='/'; strcpy(dir+l,sub);

      if (!Wb::fexist(dir,'d')) { if (!mflag) { return -2; }
         int e=0; mode_t p=(S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         e=mkdir(dir,p); if (e) wblog(FL,
            "ERR %s() failed to create directory (e=%d)\n%s",FCT,e,dir
         );
      }
      l+=m;
   }

   return l;
};


template <class TQ>
int RCStore::get_directory(const char *F, int L,
   wbstring &file, const QSet<TQ> &Q, const char *ext, char mflag) const {

   unsigned l,l2, m=28, n=100; { file.init(n+m); }
   char *f=file.data, *sub=file.data+n;

   if (Q.qdir.len<2 || Q.qdir.len>99) wblog(F_L,
      "ERR %s() got empty CData !?\n%s",FCT,Q.toStr().data);
   l2=snprintf(sub,m,"%s",Q.qdir.toStr().data);

   l=get_directory(F_L,f,n,Q.t,"CStore",sub,mflag);
   if (int(l)<0) { return l; }

   l2=snprintf(sub,m,"%s",Q.QStrS(0,';').data);

   if (ext && !strcmp(ext,"mp3")) {
      char *s=strstr(sub,";"); if (s) s[0]=0;";_"); }
   }

   l+=snprintf(f+l,n-l,"/(%s).%s",sub, ext?ext:"mat");

   if (l>=n || l2>=m) wblog(FL,
      "ERR %s() string out of bounds (%d/%d; %d/%d)\n%s",FCT,l,n,l2,m,f);
   else if (l) file.len=l+1;

   return Wb::fexist(f,'f');
};


template <class TQ, class TD>
int RCStore::get_directory(
   const char *F, int L, wbstring &file,
   const genRG_base<TQ,TD> &R, const char *ext, char mflag
 ) const {

   unsigned l, n=100; { file.init(n); }
   char *f=file.data;

   char sep[2]=" ";
   if (R.q.qlen()<2 || R.J.wbvector<TQ>::allIn(0,9)) sep[0]=0;

   if (!R.J.len || R.J.len>99) wblog(F_L,"ERR %s() "
      "got empty symmetry labels !? (J=[%s])",FCT,R.J.toStr().data);

   l=get_directory(F_L,f,n,R.q,"RStore",NULL,mflag);
   if (int(l)<0) { return 0; }

   l=snprintf(f+l,n-l,"/(%s).%s",R.J.toStrf("",sep).data, ext?ext:"mat");
   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,f);
   else if (l) file.len=l+1;

   return Wb::fexist(f,'f');
};


int RCStore::get_directory(
   const char *F, int L, wbstring &file,
   const cgc_contract_id<MTI> &idc, const char *ext, char mflag
 ) const {

   unsigned l,rc, n=256; { file.init(n); }
   char *f=file.data, sub[8];

   ctrIdx ica,icb;
   CData<gTQ,RTD> a,b;

   idc.extract(a,ica,b,icb);
   rc=a.qdir.len+b.qdir.len-ica.len-icb.len;

   if (!a.qs.len || a.qs.len>99 || !b.qs.len || b.qs.len>99) {
      wblog(F_L,"ERR %s() got invalid symmetry labels !?\n"
      "A: %s\nB: %s",FCT,a.toStr().data,b.toStr().data); }
   if (rc>99) wblog(FL,
      "ERR %s() got unexpected rc=%d !?\na: %s\nb: %s",
      FCT, rc, a.toStr().data, b.toStr().data
   );

   QDir cdir; cdir.init2val(rc,+1);
   l=rc - a.qdir.nconj(ica) - b.qdir.nconj(icb);
   for (; l<rc; ++l) cdir[l]=-1;

   if (rc)"rank%d",rc);
        { snprintf(sub,8,"%s",cdir.toStr().data); }
   else { strcpy(sub,"scalar"); }

   l=get_directory(F_L,f,n,a.t,"XStore",sub,mflag);
   if (int(l)<0) { return 0; }

   l+=snprintf(f+l,n-l,"/(%s)_%s%s(%s)_%s.%s",
      a.QStr().data, ica.toStr().data, ica.conj?"":" ",
      b.QStr().data, icb.toStr().data,"":" ",
      ext?ext:"mat"
   );

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,f);
   else if (l) file.len=l+1;

   return Wb::fexist(f,'f');
};


template <class TQ, class TD>
int RCStore::save_CData(const char *F, int L, const CData<TQ,TD> &A) {

   wbstring fs;
   int q=get_directory(F_L,fs,(const QSet<TQ>&)A,"cgd",'!');

   if (CG_VERBOSE>1 && F) wblog(PFL,
      "[%s] %s() %s",q?"W":"w",FCT,fs.data+root.len);

   if (!A.cstat.ctime) wblog(FL,
      "ERR %s() got empty time stamp / ID\n%s\n%s",
      FCT,A.toStr().data, A.cstat.toStr('V').data);
   if (A.cstat==CGD_REF_INIT) wblog(FL,
      "ERR %s() got REF_INIT CData !?\n%s\n%s",
      FCT,A.toStr().data, A.cstat.toStr('V').data
   );

   mxArray *a;
   Wb::matFile f(FL,fs.data,A.cgd.D.len<(1<<20) ? "w:std":"w");


   CData<TQ,TD> Ar; Ar.RefInit(FL,A);
   Ar.cstat.t=A.cstat.t;

   a=Ar.toMx();
   matPutVariable(f.mfp,"CRef",a);
   mxDestroyArray(a);

   a=A.cgd.toMx();
   matPutVariable(f.mfp,"cdata",a);
   mxDestroyArray(a);

   return q;
};


template <class TQ, class TD>
int RCStore::load_CData(
   const char *F, int L, const QSet<TQ> &Q, CData<TQ,TD> &A,
   char bflag
){
   int q; wbstring fs; mxArray *a;
   if ((q=get_directory(F_L,fs,Q,"cgd"))<=0) { if (F) wblog(FL,
      "ERR %s() missing CData file\n%s",FCT,fs.data);
      return 0;
   }

   if (CG_VERBOSE>4 && F) wblog(PFL,"[r] %s() %s",FCT,fs.data+root.len);

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   a=matGetVariable(f.mfp,"CRef");
   if (!a) wblog(FL,"ERR %s() "
      "failed to read variable 'CRef' from file\n%s",FCT,fs.data);
   q=A.init(FL,a,0,bflag,0);
   mxDestroyArray(a);

   a=matGetVariable(f.mfp,"cdata");
   if (!a) wblog(FL,"ERR %s() "
      "failed to read variable 'cdata' from file\n%s",FCT,fs.data);
   A.cgd.init(FL,a,0,A.t);
   mxDestroyArray(a);

   A.checkNormSign(FL);

   return q;
};


template <class TQ>
int RCStore::save_Std3(
   const char *F, int L, const QType &t, const qset<TQ> &J12) const {

   wbstring fs;

   QSet<TQ> Q; {
      Q.init3(t); if (J12.len) {
      memcpy(Q.qs.data,J12.data,J12.len*sizeof(TQ)); }
   }


   int q=get_directory(F_L,fs,Q,"mp3",'!');

   if (CG_VERBOSE>1 && F) wblog(PFL,
      "[%s] %s() %s",q?"W":"w",FCT,fs.data+root.len);

   Wb::matFile f(FL,fs.data,"w:std");

   mxArray *a=map3_CreateStructMatrix(1,1);
   map3_add2MxStruct(a,0,t,J12);

   matPutVariable(f.mfp,"map3",a);
   mxDestroyArray(a);

   return q;
};


template <class TQ>
int RCStore::load_Std3(const char *F, int L,
   const QType &t, const qset<TQ> &J1, const qset<TQ> &J2) const {

   if (!J1.len || J1.len!=t.qlen() || J1.len!=J2.len) wblog(FL,
      "ERR %s() invalid J12=[%s, %s] having %s",
      FCT,J1.toStr().data,J2.toStr().data,t.toStr().data);

   wbstring fs;

   QSet<TQ> Q; { Q.init3(t);
      memcpy(Q.qs.data,       J1.data, J1.len*sizeof(TQ));
      memcpy(Q.qs.data+J1.len,J2.data, J2.len*sizeof(TQ));
   }

   int q=get_directory(F_L,fs,Q,"mp3");
   if (q<=0) {
      if (!F) return -1;

      genRG_struct<TQ,RTD> &B=gRG.buf[t];
      if (!B.q.isKnown()) { B.checkInit(FL,t); }

      q=B.getTensorProdReps_gen(J1,J2,'!');
      if (q>0) {
         qset<TQ> J12(J1,J2);
         return gCG.map3[t][J12].size();
      }

      q=get_directory(F_L,fs,Q,"mp3");
   }

   if (q<=0) { if (F) wblog(FL,
      "ERR %s() missing map3 file\n%s",FCT,fs.data);
      return 0;
   }

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   mxArray *a, *S=matGetVariable(f.mfp,"map3");
   if (!S) wblog(FL,"ERR %s() "
      "failed to read variable 'map3' from file\n%s",FCT,fs.data);

   unsigned n=mxGetNumberOfElements(S); int m;
   int i3=(mxIsStruct(S) ? mxGetFieldNumber(S,"cgr") : -1);

   if (n!=1 || i3<0) wblog(FL,
      "ERR %s() invalid map3 data in RCStore (%d,%d)\n%s",
      FCT, n,i3, fs.data+root.len);

   a=mxGetFieldByNumber(S,0,i3);
   m=gCG.add3(t,a);

   if (CG_VERBOSE>4 && L) {
      Q.qs.len=2*J1.len; Q.qdir.len=2; wblog(PFL,
         "(+) map3[%03d] %s(%d)",gCG.map3[t].size(), Q.toStr().data, m
      );
   }

   if (m<=0) {
      if (m==0) wblog(FL,"%s() already got data internally",FCT);
      else wblog(FL,"ERR %s() m=%d",FCT,m);
   }



   mxDestroyArray(S);
   return m;
};


template <class TQ, class TD>
int RCStore::save_RSet(
   const char *F, int L, const genRG_base<TQ,TD> &R, char xflag) {

   wbstring fs;
   int q=get_directory(F_L,fs,R,"rep",'!');

   if (q>0 && xflag) {
      genRG_base<TQ,TD> X=R;
      gStore.load_RSet(0,0,X.q,X.J);

      double e=X.normDiff(FL,R);
      if (e>1E-12) wblog(FL,"ERR %s() "
         "got RSet[%s] inconsistency (%.3g)",FCT,X.J.toStr().data, e);
      return 2;
   }

   if (CG_VERBOSE>1 && F) wblog(PFL,
      "[%s] %s() %s",q?"W":"w",FCT,fs.data+root.len);

   mxArray *a=R.toMx();
   Wb::matFile f(FL,fs.data, R.Z.dim1<(1<<16)?"w:std":"w");

   matPutVariable(f.mfp,"RSet",a);
   mxDestroyArray(a);

   return q;
};


template <class TQ>
int RCStore::load_RSet(
   const char *F, int L, const QType &t, const qset<TQ> &J
){
   wbstring fs;
   genRG_base<TQ,RTD> R; R.q=t; R.J=J;

   if (get_directory(F_L,fs,R,"rep")<=0) { if (F) wblog(FL,
      "ERR %s() got non-existing file\n%s",FCT,fs.data);
      return 0;
   }

   if (CG_VERBOSE>4 && L) wblog(PFL,
      "(+) RBUF[%03d] %s [%s]",gRG.buf[t].RSet.size(),
      t.toStr().data,J.toStr().data);

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   unsigned n;
   mxArray *a=matGetVariable(f.mfp,"RSet");

   if (!a) wblog(FL,"ERR %s() "
      "failed to read variable 'RSet' from file\n%s",FCT,fs.data);
   if ((n=mxGetNumberOfElements(a))!=1) wblog(FL,
      "ERR %s() got invalid RSet (len=%d)",FCT,n);

   gRG.add(t,a,0);
   mxDestroyArray(a);

   return 1;
};


int RCStore::save_XMap(
   const char *F, int L, const cgc_contract_id<MTI> &idc
 ) const {

   const auto it=gCX.XBUF.find(idc);
   if (it==gCX.XBUF.end()) wblog(FL,"ERR %s() key not found !?",FCT); 
   const x3map<gTQ,RTD> &M=it->second;


   wbstring fs;
   int q=get_directory(F_L,fs,idc,"x3d",'!');

   if (CG_VERBOSE>1 && F) { 
      unsigned l=root.len, i=l;
      for (; i<fs.len && fs[i]; ++i) { if (fs[i]=='/') l=i+1; }
      wblog(PF_L,"[%s] %s",q?"W":"w",fs.data+l);
   }

   mxArray *a=M.toMx(&idc);

   Wb::matFile f(FL,fs.data,"w:std");

   matPutVariable(f.mfp,"x3m",a);
   mxDestroyArray(a);

   return q;
};


int RCStore::load_XMap(
   const char *F, int L, const cgc_contract_id<MTI> &idc
 ) const {

   wbstring fs;
   if (get_directory(F_L,fs,idc,"x3d")<=0) { if (F) wblog(FL,
      "ERR %s() got non-existing file\n%s",FCT,fs.data);
      return 0;
   }

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   unsigned n;
   mxArray *a=matGetVariable(f.mfp,"x3m");

   if (!a) wblog(FL,"ERR %s() "
      "failed to read variable 'RSet' from file\n%s",FCT,fs.data);
   if ((n=mxGetNumberOfElements(a))!=1) wblog(FL,
      "ERR %s() got invalid RSet (len=%d)",FCT,n);

   gCX.add(F_L,idc,a,0);
   mxDestroyArray(a);

   return 1;
};


template <class TQ, class TD>
double genRG_base<TQ,TD>::normDiff(const char *F, int L,
  const genRG_base<TQ,TD> &B, TD eps
) const {

   double e, emax=0;

   if (Sp.len!=B.Sp.len || Sz.len!=B.Sz.len || !Z.hasSameSize(B.Z))
       wblog(F_L,"ERR size mismatch (%d/%d; %d/%d; %s <> %s)",
       Sp.len, B.Sp.len, Sz.len, B.Sz.len,
       Z.sizeStr().data, B.Z.sizeStr().data);

   e=Z.normDiff(B.Z); if (e>eps) {
      if (F) {
          MXPut(F_L).add(Z,"Z").add(B.Z,"Z2").add(e,"e"); 
          wblog(F_L,"ERR %s() %.3g",FCT,e); 
      } else wblog(FL,"WRN %s() %.3g",FCT,e); 
      return e;
   }
   emax=MAX(emax,e);

   e=J.normDiff(B.J); if (e>eps) {
      if (F) {
         MXPut(F_L).add(Z,"Z").add(B.Z,"Z2")
         .add(J,"J").add(B.J,"J2").add(e,"e"); 
         wblog(F_L,"ERR %s() J (e=%.3g)",FCT,e);
      } else wblog(FL,"WRN %s() J (e=%.3g)",FCT,e);
      return e;
   }
   emax=MAX(emax,e);

   for (unsigned i=0; i<Sp.len; i++) {
      e=Sp[i].normDiff2(B.Sp[i])/Sp[i].numel();
      e=std::sqrt(ABS(e)); if (e>eps) {
         wblog(F_L,"TST %s() Sp[%d] (e=%.3g)",FCT,i+1,e);
         return e;
      }
      emax=MAX(e,emax);
   }
   if (emax && emax>1E-12) wblog(FL,"TST %s() e=%g",FCT,emax); 

   for (unsigned i=0; i<Sz.len; i++) {
      e=Sz[i].normDiff2(B.Sz[i]);
      e=std::sqrt(ABS(e)); if (e>eps) {
         wblog(F_L,"TST %s() Sz[%d] (e=%.3g)",FCT,i+1,e);
         return e;
      }
      emax=MAX(e,emax);
   }
   if (emax && emax>1E-12) wblog(FL,"TST %s() e=%g",FCT,emax); 

   return emax;
};


template <class TQ, class TD>
void genRG_base<TQ,TD>::checkOrthoCommRel(
   const char *F, int L, const genRG_base<TQ,TD> &B
 ) const {

   if (&B==this) {
      wblog(FL,"WRN %s() got same genRG_base set",FCT);
      return;
   }

   unsigned i=Sp.len+Sz.len, j=B.Sp.len+B.Sz.len, e=0;
   wbvector< const wbSparrayTD* > S1(i), S2(j);
   wbSparrayTD C; double x;

   for (j=i=0; i<  Sp.len; ++i, ++j) S1[j]=&(  Sp[i]);
   for (  i=0; i<  Sz.len; ++i, ++j) S1[j]=&(  Sz[i]);

   for (j=i=0; i<B.Sp.len; ++i, ++j) S2[j]=&(B.Sp[i]);
   for (  i=0; i<B.Sz.len; ++i, ++j) S2[j]=&(B.Sz[i]);

   for (i=0; i<S1.len; ++i) {
   for (j=0; j<S2.len; ++j) {
      S1[i]->comm(FL,*S2[j],C);
      if ((x=C.norm())>TD(1E-10)) { 
         sprintf(str,"[ S(%d), S(%d) ] @ %.3g",i+1,j+1,x);
         e=1; break;
      }
      S1[i]->comm(FL,*S2[j],C,'N','C');
      if ((x=C.norm())>TD(1E-10)) { 
         sprintf(str,"[ S(%d), S(%d)' ] @ %.3g",i+1,j+1,x);
         e=2; break;
      }
   }  if (e) break; }

   if (e) {
      MXPut(FL,"a").addP(str,"info").add(*this,"R1").add(B,"R2")
        .add(*S1[i],"a").add(*S2[j],"b").add(C,"c").add(x,"x");
      wblog(F_L,"ERR %s() got non-commuting symmetries\%s",FCT,str);
   }
};


template <class TQ, class TD>
genRG_base<TQ,TD>& genRG_base<TQ,TD>::Sort() {

    if (Z.isEmpty()) {
       MXPut(FL,"ans").add(*this,"R");
       wblog(FL,"ERR %s() got empty Z !?",FCT); 
    }

    unsigned i; wbperm P; wbMatrix<double> z2(Z);
    z2.FlipCols();
    z2.sortRecs_float(P,-1);
    
    Z.Set2Recs(P);

    for (i=0; i<Sp.len; ++i) { Sp[i].MatPermute(P); }
    for (i=0; i<Sz.len; ++i) { Sz[i].MatPermute(P); }

#ifdef CG_CHECK_MW_PERM
    if (P0.isEmpty()) P0=P;
    else P0.Select(P);
#endif

    return *this;
};


template <class TQ, class TD>
genRG_base<TQ,TD>& genRG_base<TQ,TD>::ApplyQFac(
   const wbvector<double> &qfac,
   const wbarray<double> *JM
){
   if (J.len!=Z.dim2) wblog(FL,
      "ERR %s() severe Q-label mismatch (%d/%d)",FCT,J.len,Z.dim2);
   if (qfac.len>1 && qfac.len!=Z.dim2) wblog(FL,
      "ERR %s() qfac: length mismatch (%d/%d)",FCT,qfac.len,Z.dim2);
   if (JM && (!JM->isSMatrix() || JM->SIZE[0]!=Z.dim2)) wblog(FL,
      "ERR %s() J-map size mismatch (%dx%d <> %s)",
       FCT,Z.dim1,Z.dim2, JM->sizeStr().data); 

   unsigned i,j;

   if (qfac.len) {
      if (qfac.len==1) { double x=qfac[0]; J*=x; Z*=x; }
      else {
         const double *x=qfac.data;
         for (j=0; j<Z.dim2; ++j) { J[j]*=x[j];
         for (i=0; i<Z.dim1; ++i) { Z(i,j)*=x[j]; }}
      }
   }

   if (JM) {
      wbarray<double> x2, x(1,J.len);
                           for (i=0; i<J.len; i++) x[i]= J[i];
      wbMatProd(x,*JM,x2); for (i=0; i<J.len; i++) J[i]=x2[i];
   }

   SkipTiny(); return *this;
};


template <class TQ, class TD>
void genRG_base<TQ,TD>::compareStdSU2(const char *F, int L) const {

   cgsparray xSp,xSz,xS2,xE; 
   wbvector<double> xsz;
   double e;

   if (!q.isSU2() || J.len!=1) wblog(F_L,
      "ERR invalid %s qset record [%s]",q.toStr().data, J.toStr().data);
   if (Sp.len!=1 || Sz.len!=1) wblog(F_L,
      "ERR got invalid SU(2) generator set (Sp[%d], Sz[%d]) !??",
      Sp.len,Sz.len);

   get_SU2mat(J[0],xsz,xSp,xSz,xS2,xE);
   xSz*=2; xsz*=2;

   if ((e=Sp[0].normDiff(xSp))>1E-12)
      wblog(F_L,"ERR got SU(2) inconsistency (Sp @ %.3g)",e); else
   if ((e=Sz[0].normDiff(xSz))>1E-12)
      wblog(F_L,"ERR got SU(2) inconsistency (Sz @ %.3g)",e); else
   if ((e=Z.normDiff(wbMatrix<double>().init2ref(xsz,'t')))>1E-12)
      wblog(F_L,"ERR got SU(2) inconsistency (Z @ %.3g)",e
   );
};


template <class TQ, class TD>
mxArray* genRG_base<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {
#ifdef CG_CHECK_MW_PERM
   const char *fields[]={"type","J","Z","Sp","Sz","P"};
   return mxCreateStructMatrix(m,n,6,fields);
#else
   const char *fields[]={"type","J","Z","Sp","Sz"};
   return mxCreateStructMatrix(m,n,5,fields);
#endif
};

template <class TQ, class TD>
void genRG_base<TQ,TD>::add2MxStruct(mxArray *S, unsigned i) const {

   if (q.type==QTYPE_UNKNOWN) {
      if (J.isEmpty())
         mxSetFieldByNumber(S,i,0, wbstring().toMx());
      else {
         mxSetFieldByNumber(S,i,0, q.toMx());
      }
   }
   else mxSetFieldByNumber(S,i,0, q.toMx());

   mxSetFieldByNumber(S,i,1, J.toMx());
   mxSetFieldByNumber(S,i,2, Z.toMx());
   mxSetFieldByNumber(S,i,3, Sp.toMx());
   mxSetFieldByNumber(S,i,4, Sz.toMx());

#ifdef CG_CHECK_MW_PERM
   mxSetFieldByNumber(S,i,5, P0.toMx());
#endif
};


template <class TQ, class TD>
void genRG_struct<TQ,TD>::initCommRel(
   const char *F, int L,
   const genRG_base<TQ,TD> &R, char vflag
){
   unsigned i,j,k,p,n, np=R.Sp.len, nz=R.Sz.len;
   wbsparray<TD> C; double x;
   TD dz;

   if (!np || !nz) wblog(F_L,"ERR %s() got empty RSet (%d,%d)",FCT,np,nz);

   CR.init(np);

   for (i=0; i<np; ++i) {
      if (R.Sp[i].isEmpty()) wblog(FL,
         "ERR %s() got empty R.Sp[%d/%d] !??",FCT,i+1,np);

      R.Sp[i].comm(FL,R.Sp[i],C,'N','C');
      if (C.norm2()<(TD)(1E-10)) {
         MXPut(FL,"ans").add(R,"R").add(i+1,"i");
         wblog(F_L,"ERR %s() [Sp,Sp'] has norm 0 !??",FCT);
      }

      CR[i].init(i,i,nz);

      wbvector<unsigned> &kk=CR[i].k;
      wbvector<double> &fac=CR[i].fac;

      for (p=k=0; k<nz; k++) {
         x=C.froNorm2(R.Sz[k])/R.Sz[k].froNorm2(R.Sz[k]);
         if (ABS(x)>1E-10) { kk[p]=k; fac[p++]=x; }
      }
      if (!p) {
         MXPut(FL,"i").add(R,"R").add(CR,"CR").add(C,"C")
          .add(p+1,"p").add(np,"np");
         wblog(F_L,"ERR %s() got no overlap of [Sp,Sp'] with Sz "
         "(%d/%d)",FCT,i+1,np);
      }
      kk.len=p; fac.len=p;
   }

   DZ.init(np,nz);

   for (i=0; i<np; ++i) {
   for (j=0; j<nz; ++j) {
      R.Sz[j].comm(FL,R.Sp[i],C);
      if (C.sameUptoFac(R.Sp[i],&dz)) { DZ(i,j)=dz;
         MXPut(FL,"cr").add(R,"R").add(DZ,"DZ").add(i+1,"i").add(j+1,"j");
         wblog(FL,"ERR %s() invalid Sp/Sz operators (%s):\n"
         "[Sz,Sp] ~ Sp not satisfied",FCT,q.toStr().data); }
      DZ(i,j)=dz;
   }}

   checkCommRel(F_L,R);

   if (vflag) {
      wblog(FL,
         "<i> summary of %s commutator relations",q.toStr().data);
      for (i=0; i<CR.len; i++) { n=CR[i].k.len;
      for (k=0; k<n; k++) {
         printf("%4d: %2d %2d %2d : %8g\n", i+1,
         CR[i].i+1, CR[i].j+1, CR[i].k[k]+1, CR[i].fac[k]);
      }}
   }

};


template <class TQ, class TD>
double genRG_struct<TQ,TD>::checkCommRel(
  const char *F, int L, const genRG_base<TQ,TD> &R
) const {

   unsigned i,j,p,l,m, n=CR.len, np=R.Sp.len, nz=R.Sz.len;
   wbsparray<TD> C; double e2, r2=0;
   QType qt=(q==QTYPE_UNKNOWN ? R.q : q);

   if (q!=R.q) {
      sprintf(str,"got QType mismatch: %s <> %s",
         q.toStr().data, R.q.toStr().data);

      if (q!=QTYPE_UNKNOWN)
           wblog(F_L,"ERR %s() %s",FCT,str);
   }

   if (qt!=QTYPE_UNKNOWN && (np!=nz || nz!=qt.sub)) wblog(F_L,
      "ERR %s() unexpected Sz[%d], Sp[%d], sub=%d\nhaving symmetry='%s'",
      FCT, nz, np, qt.sub, qt.toStr().data); 
   if (!R.Sp.len || !R.Sz.len) wblog(F_L,
      "ERR %s() CR not yet initialized",FCT);

   for (i=0; i<nz; i++) {
   for (j=i+1; j<nz; j++) {
      e2=R.Sz[i].froNorm2(R.Sz[j]);
      if (e2>1E-10) {
         if (F) {
           MXPut(FL,"I_").add(R,"R").add(i+1,"i").add(j+1,"j").add(C,"C");
           wblog(F_L,"ERR %s() CR inconsistency:\n"
           "z-ops not mutually orthogonal (%.3g)",FCT,e2);
         }
      }
      R.Sz[i].comm(FL,R.Sz[j],C); e2=C.norm2(); r2+=e2;
      if (e2>1E-10) {
         if (F) {
           MXPut(FL,"I_").add(R,"R").add(i+1,"i").add(j+1,"j").add(C,"C");
           wblog(F_L,"ERR %s() CR inconsistency ([Z,Z]: %.3g)",FCT,e2);
         }
      }
   }}

   for (l=0; l<n; l++) {
      const wbvector<double> &fac=CR[l].fac;
      const unsigned *k=CR[l].k.data; i=CR[l].i; j=CR[l].j;

      if (i>=np || j>=np || CR[l].k.max()>=nz) wblog(F_L,
         "ERR %s() index out of bounds (%d,%d,%s/%d)",
          FCT, i+1, j+1, (CR[l].k+1).toStr().data, n
      );

      R.Sp[i].comm(FL,R.Sp[j],C,'N','C'); m=fac.len;

      for (p=0; p<m; ++p) { C.Plus(FL,R.Sz[k[p]],-fac[p]); }

      e2=C.norm2(); r2+=e2;
      if (e2>1E-10) {
         if (F) {
             MXPut(FL,"q").add(R,"R").add(i+1,"i").add(j+1,"j")
             .add(CR[l].k+1,"k").add(fac,"fac").add(C,"C");
             wblog(F_L,"ERR %s() CR inconsistency (%d,%d: %.3g)",
             FCT,i+1,j+1,e2);
         }
      }
   }

   if (DZ.dim1!=np || DZ.dim2!=nz) wblog(FL,"ERR %s() "
      "got empty DZ (%s; %s)",FCT,q.toStr().data,DZ.sizeStr().data);

   for (i=0; i<np; ++i) {
   for (j=0; j<nz; ++j) {
      R.Sz[j].comm(FL,R.Sp[i],C); C.Plus(FL,R.Sp[i],-DZ(i,j));
      e2=C.norm2(); r2+=e2;

      if (e2>1E-10) {
         if (F) {
           MXPut(FL,"q").add(R,"R")
           .add(R.Sz[j],"Sz").add(R.Sp[i],"Sp").add(DZ(i,j),"sfac")
           .add(*this,"I").add(i+1,"i").add(j+1,"j").add(C,"C");
           wblog(F_L,"ERR %s() CR inconsistency ([Z,Sp]: %.3g)",FCT,e2);
         }
      }
   }}

   e2=sqrt(r2);
   if (e2>CG_EPS1) wblog(FL,"ERR %s() CR error = %.3g",FCT,e2); else
   if (e2>CG_EPS2) wblog(FL,"WRN %s() CR @ %.3g",FCT,e2);


   return r2;
};


template <class TQ, class TD>
unsigned genRG_struct<TQ,TD>::getTensorProdReps(const qset<TQ> &J) {

    unsigned i=0, l=0, m=0, n=0;
    wbvector< const qset<TQ>* > qq(RSet.size());

    for (crMAP I=RSet.begin(); I!=RSet.end(); ++I, ++i) {
       const genRG_base<TQ,TD> &R=I->second;
       if (!R.Sz.isEmpty() && R.Z.dim1 && R.Z.dim1<=32*q.qlen()) {
       qq[m++]=&(I->first); }
    }

    for (i=0; i<m; ++i) {
       l=getTensorProdReps_gen(J,*qq[i]);
       if (int(l)>0) n+=l;
    }

    return n;
};


template <class TQ, class TD>
unsigned genRG_struct<TQ,TD>::getTensorProdReps(
    unsigned dmax, char vflag) {

    unsigned i,j, l=0, n3=0;
    if (int(dmax)<0) { dmax=MAX(8U,2*q.sub); }
    unsigned dmax2=dmax*dmax;

    wbvector< const qset<TQ>* > qq(RSet.size());
    wbvector<unsigned> dd(RSet.size());

    for (crMAP I=RSet.begin(); I!=RSet.end(); ++I) {
       const genRG_base<TQ,TD> &R=I->second;
       if (!R.Sz.isEmpty() && R.Z.dim1 && R.Z.dim1<=dmax2) {
          qq[l]=&(I->first); 
          dd[l]=R.Z.dim1; ++l;
       }
    }

    if (!l) wblog(FL,"ERR %s() no RSets found (d<=%d) !?",FCT,dmax);

    if (vflag) printStatus(FL,"old");

    for (i=0; i<l; ++i)
    for (j=0; j<l; ++j) {
       if (dd[i]*dd[j]<=dmax2) {
       n3+=getTensorProdReps_gen(*qq[i],*qq[j]);
    }}

    if (vflag) printStatus(FL,"new");

    return n3;
};


template <class TQ, class TD>
void genRG_struct<TQ,TD>::printStatus(
   const char *F, int L, const char *istr
){
    unsigned n=0;

    wblog(F_L,"<i> %s status for %s: %d entries",
       istr && istr[0] ? istr:"current", q.toStr().data, RSet.size());
    wbvector<qset<TQ> > qq(RSet.size()); {
       for (crMAP I=RSet.begin(); I!=RSet.end(); ++I) {
          const genRG_base<TQ,TD> &R=I->second;
          printf("   %4d : [%s] (d=%d)%s\n",
          ++n, I->first.toStr().data, R.Z.dim1,
          R.Sz.isEmpty() || !R.Z.dim1 ? "  *** EMPTY ***":"");
       }
    }
};


template <class TQ, class TD>
void genRG_struct<TQ,TD>::initSU2(const qset<TQ> &J) {

   if (!q.isSU2() || J.len!=1) wblog(FL,"ERR %s() invalid setting "
      "(%s with J=[%s])",FCT,q.toStr().data,J.toStr().data
   ); 
      
   wbsparray<TD> Sp,Sz,S2,E; 
   wbvector<double> sz;

   genRG_base<TQ,TD> &R=RSet[J];

   get_SU2mat(J[0],sz,Sp,Sz,S2,E);

   if (R.Sp.len || R.Sz.len) {
      if (R.q!=q || R.J!=J) wblog(FL,
         "ERR %s() already got initialized (%s <> %s; [%s] <> [%s]) !??",
         FCT, R.q.toStr().data, q.toStr().data,
         R.J.toStr().data, J.toStr().data);
      if (R.Sp.len!=1 || R.Sz.len!=1) wblog(FL,
         "ERR %s() already got initialized (Sp[%d], Sz[%d]) !??",
         FCT, R.Sp.len, R.Sz.len);
      if (R.Sp[0].SIZE!=Sp.SIZE) wblog(FL,"ERR %s() "
         "already got Sp at %s (vs. %s) !??",
         FCT,R.Sp[0].sizeStr().data, Sp.sizeStr().data);
      if (R.Sz[0].SIZE!=Sz.SIZE) wblog(FL,"ERR %s() "
         "already got Sz at %s (vs. %s) !??",
         FCT,R.Sz[0].sizeStr().data, Sz.sizeStr().data
      );

      double e;
      if ((e=R.Sp[0].normDiff(Sp))>1E-12)
         wblog(FL,"ERR got SU(2) inconsistency (Sp @ %.3g)",e);
      if ((e=R.Sz[0].normDiff(Sz))>1E-12)
         wblog(FL,"ERR got SU(2) inconsistency (Sz @ %.3g)",e);
   }
   else {
      R.Sp.init(1); Sp.save2(R.Sp[0]);
      R.Sz.init(1); Sz.save2(R.Sz[0]); R.q=q; R.J=J;

      R.Z.init(sz.len,1);
      for (unsigned i=0; i<sz.len; ++i) { R.Z.data[i]=sz[i]; }
   }

#ifdef CG_CHECK_MW_PERM
   if (R.P0.len) R.P0.init();
#endif
};


template <class TQ, class TD>
int genRG_struct<TQ,TD>::getTensorProdReps_gen(
   const qset<TQ> &J1_, const qset<TQ> &J2_,
   char tflag
){
   unsigned np=0, n3=0, mp=(J2_<J1_);
   int n1, n2;

   const qset<TQ> &J1 (mp ? J2_ : J1_);
   const qset<TQ> &J2 (mp ? J1_ : J2_);

   qset<TQ> J12(J1,J2);
   qset<TQ> J21(J2,J1);
   
   mp=(J1!=J2 ? 2 : 1);

   n1=gCG.find_map3(q,J12);
   if (n1<=0 && !tflag) {
      gStore.load_Std3(FL,q,J1,J2);
      n1=gCG.find_map3(q,J12);
      if (n1>0) n3+=n1;
   }

   if (mp>1) {
      n2=gCG.find_map3(q,J21);
      if (n2<=0 && !tflag) {
         gStore.load_Std3(FL,q,J2,J1);
         n2=gCG.find_map3(q,J21);
         if (n2>0) n3+=n2;
      }

      if (n1!=n2) {
         wblog(FL,"WRN missing permuted standard CGC pair !?");

         if (n1) sprintf(str,"%d output multiplet%s",n1,n1!=1?"s":"");
         else strcpy(str,"empty");
         wblog(FL,"WRN (%s,%s) %s",J1.toStr().data,J2.toStr().data,str);

         if (n2) sprintf(str,"%d output multiplet%s",n2,n2!=1?"s":"");
         else strcpy(str,"empty");
         wblog(FL,"WRN (%s,%s) %s",J2.toStr().data,J1.toStr().data,str);
      }
      if (n1>0 && n2>0) return n3;
   }
   else if (n1>0) return n3;

   const genRG_base<TQ,TD> &G1=RSet[J1], &G2=RSet[J2];

   if (G1.isEmpty()) gStore.load_RSet(FL,q,J1);
   if (G2.isEmpty()) gStore.load_RSet(FL,q,J2);

   const wbvector< wbsparray<TD> >
         &Sp1=G1.Sp, &Sp2=G2.Sp,
         &Sz1=G1.Sz, &Sz2=G2.Sz;

   if (!Sp1.len && !tflag) { gRG.getR(FL,q,J1.data); }
   if (!Sp2.len && !tflag) { gRG.getR(FL,q,J2.data); }

   if (!Sp1.len || !Sp2.len || !Sz1.len || !Sz2.len) {
      if (tflag) return -1;
      if (q.isSU2()) {
         if (!Sp1.len || !Sz1.len) initSU2(J1);
         if (!Sp2.len || !Sz2.len) initSU2(J2);
      }
      else wblog(FL,
        "ERR %s() Sp/Sz not yet generated\n"
        "for %s ([%s]*[%s]: %d/%d; %d/%d)", FCT,
         q.toStr().data, J1.toStr().data, J2.toStr().data,
         Sp1.len, Sp2.len, Sz1.len, Sz2.len
      );
   }

   if (Sp1.len!=Sp2.len || Sz1.len!=Sz2.len) wblog(FL,
      "ERR %s()\ngot Sp/Sz inconsistency (Sp[%d/%d], Sz[%d/%d])",
      FCT, Sp1.len, Sp2.len, Sz1.len, Sz2.len);
   if (Sp1.len!=Sz1.len) wblog(FL,
      "ERR %s(): %s J1=[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, q.toStr().data, J1.toStr().data, Sp1.len, Sz1.len);
   if (Sp2.len!=Sz2.len) wblog(FL,
      "ERR %s(): %s J2=[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, q.toStr().data, J2.toStr().data, Sp2.len, Sz2.len
   );

   np=Sp1.len;

   unsigned i,it,im,m,d3, r=q.sub, nz=Sz1.len, nel=0,
         d1=Sp1[0].SIZE[0], d2=Sp2[0].SIZE[0], d12=d1*d2;
   int sr=0, isLargeD=d12>10000; double r2=0;
   wbvector< wbsparray<TD> > CM, Sp(np), Sz(nz);
   wbvector<genRG_base<TQ,TD> > RR;

   wbMatrix<unsigned> iOM;

   wbsparray<TD> X1,E1,E2;
   wbvector< wbsparray<TD> > U;
   wbMatrix<TQ> ssz,ss,sx,sz,Z;
   wbperm P, p213("2,1,3");
   wbindex I;
   TD nrm;

   const wbMatrix<double> &Z1=G1.Z, &Z2=G2.Z;
   wbvector<unsigned> dd;

   int ioC=-1, ioR=-1, io3=-1, saveC=0;

   if (G1.q!=G2.q) wblog(FL,"ERR %s() qtype inconsistency "
      "(%s; %s)",FCT,G1.q.toStr().data,G2.q.toStr().data);
   if (q.isNonAbelian() && (np!=r || nz!=r)) wblog(FL,
      "ERR %s()\ninvalid Sp/Sz set for %s (%d/%d)",
      FCT,q.toStr().data,np,nz
   );
   for (i=0; i<np; ++i) {
      if (!Sp1[i].isSMatrix(FL,d1) || !Sp2[i].isSMatrix(FL,d2))
         wblog(FL,"ERR %s()\nmatrix size mismatch (%s <> %s)",
         FCT,Sp1[i].sizeStr().data,Sp2[i].sizeStr().data
      );
   }
   for (i=0; i<nz; ++i) {
      if (!Sz1[i].isSMatrix(FL,d1) || !Sz2[i].isSMatrix(FL,d2))
         wblog(FL,"ERR %s()\nmatrix size mismatch (%s <> %s)",
         FCT,Sz1[i].sizeStr().data,Sz2[i].sizeStr().data
      );
   }
   if (Z1.dim1!=d1 || Z1.dim2!=nz) wblog(FL,
      "ERR %s() invalid Sz: %dx%d <> %dx%d (%s)",
      FCT, Z1.dim1, Z1.dim2, d1,nz, q.toStr().data);
   if (Z2.dim1!=d2 || Z2.dim2!=nz) wblog(FL,
      "ERR %s() invalid Sz: %dx%d <> %dx%d (%s)",
      FCT, Z2.dim1, Z2.dim2, d2,nz, q.toStr().data
   );

   E1.initIdentity(d1);
   E2.initIdentity(d2);

   for (i=0; i<np; ++i) { 
      Sp1[i].kron(FL,E2,Sp[i]); E1.kron(FL,Sp2[i],X1); Sp[i].Plus(FL,X1); }
   for (i=0; i<nz; ++i) { 
      Sz1[i].kron(FL,E2,Sz[i]); E1.kron(FL,Sz2[i],X1); Sz[i].Plus(FL,X1); }

   if (isLargeD) wblog(FL,
      "START getSymStates (%s) (%d*%d) \r\\",
      QSet<TQ>().init2(q,J1,J2).toStr().data, d1, d2
   ); doflush();

   r2=getSymmetryStates(FL,q,Sp,Sz,U,dd,RR,&iOM);

   if (isLargeD) wblog(FL,
      "FNSHD getSymStates (%s) (%dx%d)[%d] \r\\",
      QSet<TQ>().init2(q,J1,J2).toStr().data,d1,d2,RR.len
   ); doflush();

   for (it=0; it<RR.len; ++it) {
      const genRG_base<TQ,TD> &R=RR[it];
      const qset<TQ> J(R.J);

      CData<TQ,TD> X;
      genRG_base<TQ,TD> G, &G0=RSet[J]; d3=dd[it];

      X.init3(q,d1,d2,d3); X.cstat.init(CGD_NEW3);
      X.qs.init(J1,J2,J);

      if (!U[it].hasSize(d1*d2,d3)) wblog(FL,
         "ERR %s() size mismatch (%s <> %dx%dx%d)",
         FCT,U[it].sizeStr().data,d1,d2,d3
      ); 

      U[it].Reshape(d1,d2,d3);

      X.cgd.wbsparray<TD>::init(FL,d1,d2,d3, U[it]);

      if (CG::signFirstVal(FL,X.cgd.D.data, X.cgd.D.len)<0)
      wblog(FL,"ERR %s() *** WRONG SIGN ***",FCT);

      X.cgd.NormSignC();




      G.q=G1.q; G.J=J; G.Z=R.Z; G.Sp=R.Sp; G.Sz=R.Sz;
#ifdef CG_CHECK_MW_PERM
      G.P0=R.P0;
#endif

      G.SkipTiny();

      if (!G0.Sp.isEmpty() || !G0.Sz.isEmpty()) {
         double e=G.normDiff(0,0,G0,1E-10);
         if (e>1E-10) {
            MXPut(FL).add(J1,"J1").add(J2,"J2").add(J,"J")
              .add(G1,"G1").add(G2,"G2").add(G,"G").add(G0,"G0"); 
            wblog(FL,"ERR %s()\n"
              "inconsistency with earlier data (e=%.3g)",FCT,e);
         }
         r2+=checkCommRel(FL,R);
      }
      else {
         if (q.isSU2()) G.compareStdSU2(FL);
       
         try { r2+=checkCommRel(FL,R); }
         catch (...) {
            MXPut(FL).add(J1,"J1").add(J2,"J2").add(J,"J")
            .add(G1,"G1").add(G2,"G2").add(G,"G");
            wblog(FL,"ERR %s() inconsistency",FCT);
         };
         G.save2(G0);

         if (ioR<0) ioR=0;
         ioR+=gStore.save_RSet(0,0,G0);
      }

      if (G0.q.validType()) {
         gCG.getIdentityC(FL,q,J.data);
      }

      CData<TQ,TD> &Sb=gCG.getBUF(0,0,(QSet<TQ>&)X);

      if (!Sb.cgd.isEmpty()) {
         if (Sb.cstat==CGD_REF_INIT) {
            if (!Sb.cgd.SIZE.len) wblog(FL,
               "ERR %s() invalid CGD ref init (%d,%d)",FCT,
               Sb.cgd.SIZE.len, Sb.cgd.D.len);
            Sb.init();
         }
         else {
            if (Sb.t!=X.t || Sb.qs!=X.qs || Sb.qdir!=X.qdir) wblog(FL,
               "ERR %s() got CData inconsistency\n[%s] <> [%s]",
               FCT, Sb.QStr().data, X.QStr().data
            ); 
         }
      }

      wbvector< CRef<TQ> > &Sr=gCG.map3[q][J12][J];

      if (iOM.dim1)
           { im=iOM(it,0); m=iOM(it,1); }
      else { im=0; m=1; }

      if (!Sr.len) { Sr.init(m); } else
      if (Sr.len!=m) wblog(FL,
         "ERR %s() size mismatch map3.len=%d/%d !?",FCT,Sr.len,m);
      m=Sb.getOM(FL);


      if (mp>1) {
         CData<TQ,TD> X2;
         X.permute(X2,p213);
         sr=CG::signFirstVal(FL,X2.cgd.D.data, X2.cgd.D.len);
      }
        
      nel+=X.cgd.numel();

      if (im<m) {

         if (Sr.len>1) { unsigned m_;
wblog(FL,"TST %s() im=%d/%d : %s",FCT,im,m,Sb.toStr().data);
            Sr[im].initDecompose(FL,X.cgd, &Sb);

            if ((m_=Sb.getOM())>Sr.len || m_<m || m_>m+1) {
               MXPut(FL,"I_").add(X,"X").add(Sb,"Sb").add(Sr,"Sr");
               wblog(FL,"ERR %s() OM out of bounds (%d/%d)",FCT,m_,Sr.len);
            }
            if (m_>m) { m=m_; n3+=mp; saveC++; }
         }
         else {
            double e=X.normDiff(FL,Sb,1E-10);
            if (e<1E-10) { m=1;
               if (CG_VERBOSE>2) wblog(FL,
                  "TST %s() consistent (@ %.3g)",FCT,e); 
               Sr[im].initBase(&Sb,im,im+1); 
            }
            else {
               wblog(FL,"TST %s() 0x%lX  %d: %d %d",FCT,iOM.data,it+1,
                  Sb.getOM(FL), it<iOM.dim1 ? iOM(it,1):-1); 
               MXPut(FL).add(J1,"J1").add(J2,"J2")
                 .add(Sb,"Sb").addP(Sb.cgd.toMX(),"C0")
                 .add(X, "X" ).addP(X .cgd.toMX(),"C" )
                 .add(R,"R").add(Sp,"Sp").add(Sz,"Sz").add(U,"U")
                 .add(iOM,"iOM").add(dd,"dd")
                 .add(it+1,"it").add(d3,"d3").add(e,"e");
               wblog(FL,"ERR overwriting cgd: %s <> %dx%dx%d [%s] @ %.3g",
               X.cgd.sizeStr().data,d1,d2,d3,J.toStr().data,e);
            }
         }
      }
      else { n3+=mp; saveC++;
         if (im>m || m>=Sr.len || (im && Sb.cstat.t!=CGD_NEW3))
            wblog(FL,"ERR %s() im=%d/%d/%d (it=%d)\n%s",
            FCT, im,m,Sr.len, it, Sb.toStr().data
         );

         if (im) {
            Sb.AddMultiplicity(FL,X);
         }
         else {
            X.checkAdditivityZ(FL,G1.Z,G2.Z,R.Z);
            X.save2(Sb);
         }

         Sr[im].initBase(&Sb,im,im+1); 

         if ((m=Sb.getOM(FL))!=im+1) {
         wblog(FL,"ERR %s() m=%d/%d",FCT,m,im+1); }
      }

      Sr[im].NormStd(FL,3U);

      if (mp>1) {
         wbvector< CRef<TQ> > &Sx=gCG.map3[q][J21][J];

         if (!Sx.len) Sx.init(Sr.len); else
         if (Sx.len!=Sr.len) wblog(FL,
            "ERR %s() CRef inconsistency (%d/%d)",FCT,Sx.len,Sr.len);

         Sr[im].permute(p213,Sx[im]); if (sr<0) {
         Sx[im].cgw*=(-1); }

      }

      if (m!=Sr.len) continue;


      if (saveC){
         saveC=0; if (ioC<0) ioC=0;

         ioC+=gStore.save_CData(0,0,Sb);

         if (J.allZero()) {
            QSet<TQ> Q1(Sb);
            if (!J.len || Q1.qdir.len!=3) wblog(FL,"ERR %s()",FCT);
            if (m!=1) wblog(FL,"ERR %s() got J1 symbold @ om=%d !?",FCT,m);
            Q1.qdir.len=2;
            Q1.qs.len=2*J.len;

            CData<TQ,TD> &S1=gCG.getBUF(0,0,Q1, 0);
            if (S1.isEmpty()) {
               Sb.getJ1Symbol(FL,S1);
               int io1=gStore.save_CData(0,0,S1);

               if (CG_VERBOSE>1) {
                  wblog(PFL,"[+]  %c BUF[%03d] J1 symbol %s %8R",
                  io1? 'W':'w',gCG.BUF.size(), Q1.toStr().data,"=");
               }
            }
            else wblog(FL,"WRN %s() "
              "already got J1 symbol !?\n%s",FCT,Q1.toStr().data
            );
         }
      }

      if (CG_VERBOSE) {
         const unsigned l1=8, l2=12, l3=0;
         char s_[l1+l2+l3], *s1=s_, *s2=s_+l1;
         s1[0]=s2[0]=0;

         if (CG_VERBOSE>1) {
            if (it)
                 snprintf(s1,l1,"%s",it+1<RR.len?"":"___");
            else snprintf(s1,l1," %ld",RR.len);
         }
         else if (it==0) {
            snprintf(s1,l1,"%2ld%c",RR.len, mp>1?'X':'=');
         }

         if (iOM.dim1 && iOM(it,1)>1) 
            snprintf(s2,l2," @ OM=%d",iOM(it,1));

         wblog(PFL,"[+] %c%c map3[%02d] %-6s%-3s (%s) %4d%s",
            ioR<0? ' ': (ioR?'W':'w'), ioC<0? ' ': (ioC?'W':'w'),
            gCG.map3[q].size(), q.toStr().data,
            s1, Sb.QStr().data, d3, s2);
         doflush(); ioR=-1; ioC=-1;
      }
   }

   if (nel && nel!=d12*d12) wblog(FL,"ERR %s() "
      "state space inconsistency (%d/%dx%d)",FCT,nel,d12,d12);

   io3 =   gStore.save_Std3(0,0,q,J12); if (mp>1) {
   io3+=10*gStore.save_Std3(0,0,q,J21); }

   if (CG_VERBOSE || isLargeD || q.isLargeD(dd)) {
      if (dd.sum()!=d1*d2) wblog(FL,"ERR %s() "
         "got size inconsistency %d ~= %dx%d !?",FCT,dd.sum(),d1,d2);
      wblog(PFL,"%s  %c %s %dx%d = %d @ %.3g",
       io3%2 ? "WRN":"[+]", io3%2 ? 'W':'w',
         q.type==QTYPE_UNKNOWN ? "general" : q.toStr().data,
         d1, d2, d1*d2, sqrt(r2)); 
      if (!CG_VERBOSE) wblog(PFL,"==> [%s]", dd.toStr().data);
      fflush(0); doflush();
   }

   return n3;
};



template <class TQ, class TD>
double getSymmetryStates(const char *F, int L, const QType &q,
   const wbvector< wbSparrayTD > &Sp,
   const wbvector< wbSparrayTD > &Sz,
   wbvector< wbSparrayTD > &UK,
   wbvector<unsigned> &dd,
   wbvector<genRG_base<TQ,TD> > &RR,
   wbMatrix<unsigned> *iOM,
   char vflag
){
   unsigned i,j,ip,i0=0, D,d0=0, nu=0, m,found,
        it=0, nt, r=q.sub, np=Sp.len, nz=Sz.len;
   double r2=0;

   wbvector< wbSparrayTD* > uk;
   wbSparrayTD U,x1,x2,v0,vi;
   wbMatrix<TQ> z2,JJ;
   wbvector<INDEX_T> dJ;
   wbvector<TD> sz;
   wbperm P;

   TD x, eps=1E-8, eps2=1E-10;

#ifdef WB_SPARSE_CLOCK
   Wb::UseClock gss(&wbc_sparse_gss);
#endif


   if (!nz || np>nz) wblog(F_L,
      "ERR got invalid or empty Sp/Sz sets for %s (%d/%d)",
       q.toStr().data, np, nz);
   if (q.isNonAbelian(0,0,'l') && (np!=r || nz!=r)) wblog(FL,
      "ERR invalid number of Sz/Sp operators for %s (%d,%d/%d)",
       q.toStr().data, np,nz,r
   );

   D=Sz[0].dim();


   U.init(D,0);

   dd.init(D);
   uk.init(D);
   
   for (it=0; it<D; ++it) {
      v0.initz(D,1,1);
      for (; i0<D; ++i0) {
         v0.setRec(0, i0,0, 1.); if (!i0 && !it) break;

         wbMatProd(U,v0,x1,'C'); x=x1.norm();
         if (ABS(x-1)<eps) {
            if (ABS(x-1)>eps2) wblog(FL,
               "WRN %s() %x (%g,%g)",FCT,double(x),double(eps),double(eps2));
            continue;
         } else break;
      }

      if (i0==D) break;
      if (it) {
         wbMatProd(U,x1,x2);
         v0-=x2; v0.Normalize();
         wbMatProd(U,wbMatProd(U,v0,x1,'C'),x2);
         v0-=x2; v0.Normalize();
      }

      for (i=0; i<Sz.len; i++) {
          wbMatProd(Sz[i],v0,vi); if (vi.norm()<eps) continue;
          if (vi.sameUptoFac(v0)) wblog(FL,
             "ERR %s() failed to determine symmetry labels\n"
             "for starting vector (%d: %d,%d)",FCT,it+1,i0+1,i+1
          ); 
      }

      wbSparrayTD V(D,1);
      V.setCol(0,v0); found=1; m=0;

      while (found) { found=0;
         for (ip=0; ip<Sp.len; ++ip) {
            wbMatProd(Sp[ip],V,vi); x=vi.norm(); 
            if (x>eps) {
               vi*=(1/x); V.setCol(0,vi);
               ++found; ++m;
            }
         }
      }

      if (m && vflag && vflag!='v') wblog(FL,
         " *  applied %d Sp ops to get MW seed",m);
      found=1;


      while (found) { found=0;
      for (ip=0; ip<Sp.len; ++ip) {
          wbMatProd(Sp[ip],V,vi,'C');

          x=SQRT(vi.norm2()/vi.SIZE[1]);
          if (x<eps) {
             if (x>eps2) wblog(FL,"WRN %s() got %.3g [%g %g]",
                 FCT,double(x),double(eps),double(eps2));
             continue;
          }

          wbMatProd(V,vi,x1,'C');
          x=fabs(1-sqrt(x1.norm2()/vi.norm2()));
          if (x<eps) {
             if (x>eps2) wblog(FL,"WRN %s() got %.3g [%g %g]",
                 FCT,double(x),double(eps),double(eps2));
             continue;
          }


          wbMatProd(V,x1,x2); vi-=x2;

          wbMatProd(U,wbMatProd(U,vi,x1,'C'),x2); vi-=x2;
          x=x1.aMax(); if ((x*x)>eps) {
             MXPut(FL,"qU").add(U,"U").add(V,"V").add(vi,"vi")
               .add(Sp[ip],"Sp").add(x2,"x2");
             wblog(FL,"ERR %s() got overlap %.3g (%g)",
                FCT,double(x),double(eps)
             );
          }

          vi.OrthoNormalizeCols(FL,0,'x',eps);

          wbMatProd(U,wbMatProd(U,vi,x1,'C'),x2); vi-=x2;
          wbMatProd(V,wbMatProd(V,vi,x1,'C'),x2); vi-=x2;
          vi.OrthoNormalizeCols(FL,0,0,eps);

          V.Cat(2,vi,FL); ++found;

          if (V.SIZE[1]>D) wblog(FL,"ERR %s() "
             "too many vectors (%s; %d)",FCT,V.sizeStr().data,D); 
      }}




      U.Cat(2,V,FL); d0=V.SIZE[1]; dd[it]=d0;


      WB_NEW_1( uk[it], wbSparrayTD );
      V.save2(*uk[it]);

      if (U.SIZE.len && U.SIZE[1]==D) {
         ++it; break;
      }
      else if (U.SIZE[1]>D) {
         MXPut(FL).add(U,"U").add(V,"V").add(vi,"vi").add(x1,"x1")
         .add(it+1,"it").add(ip+1,"ip").add(dd,"dd");
         wblog(FL,"ERR %s() D=%d/%d !??",FCT,U.SIZE[1],D); 
      }
   }

   nt=it;
   if (U.SIZE[1]!=D || !nt) wblog(FL,
      "ERR %s() failed to obtain symmetry multiplets (%s/%d; %d)",
       FCT,U.sizeStr().data,D,nt);
   dd.len=nt;
   UK.init(nt);
   for (it=0; it<nt; ++it) {
      uk[it]->save2(UK[it]); WB_DELETE_1(uk[it]);
   }

#ifndef NSAFEGUARDS
   wbMatProd(U,U,x1,'C');
   if (!x1.isIdentityMatrix(eps2)) {
      MXPut(FL,"q").add(Sp,"SP").add(Sz,"SZ").add(U,"U").add(UK,"UK");
      wblog(FL,"ERR %s() new space not orthogonal (%d)",FCT,nt+1);
   }
#endif


   RR.init(nt); JJ.init(nt,nz);

   for (i=0; i<dd.len; ++i) { 
      genRG_base<TQ,TD> &R=RR[i];
      wbSparrayTD &V=UK[i]; d0=dd[i];

      R.Sp.init(np);
      for (j=0; j<np; j++) {
         wbMatProd(V,wbMatProd(Sp[j],V,x1),R.Sp[j],'C');
         R.Sp[j].Compress(0,0,CG_SKIP_REPS);
      }

      R.Sz.init(nz); R.Z.init(d0,nz);
      for (j=0; j<nz; j++) {
         wbMatProd(V,wbMatProd(Sz[j],V,x1),R.Sz[j],'C');
         if (!R.Sz[j].isDiagMatrix(eps2)) wblog(FL,
            "ERR %s() got non-diagonal z-operator !??",FCT);
         R.Sz[j].Compress(0,0,CG_SKIP_REPS);


         R.Sz[j].getDiag(FL,sz);
         R.Z.setCol(j,sz);
      }

      CG::FixRational(FL,R.Z.data,R.Z.numel(),4);

      findMaxWeight(q,R.Z,&R.J,&P); R.q=q;

      if (P.len && P.data[0]!=0 && vflag) {

          wblog(FL,"WRN state #1 in IREP-decomp is not MW "
            "(%s,%d; %d: %d/%d)", q==QTYPE_UNKNOWN ?
            "*":q.toStr('t').data, P.data[0]+1,i+1,d0,D);

          MXPut X(FL,"a"); X.add(Sp,"Sp").add(Sz,"Sz").add(UK,"UK")
           .add(dd,"dd"); if (iOM) X.add(*iOM,"M");
          X.add(D,"D").add(V,"V").add(R.Z,"Z").add(P,"P"); X.put();
      }

      if (!P.isIdentityPerm()) {
          R.Z.recPermute(P);
          V.ColPermute(P);

          for (j=0; j<np; ++j) R.Sp[j].MatPermute(P);
          for (j=0; j<nz; ++j) R.Sz[j].MatPermute(P);

          CG::rangeSignConvention(FL,V.D.data, V.D.len);
      }

      r2+=CG::FixRational(
          FL, V.D.data, V.D.len, 4, (TD)(CG_SKIP_EPS1), (TD)(CG_SKIP_EPS2)
      );
      nu+=V.D.len;

      V.Compress(FL,CG_SKIP_EPS1);

      JJ.recSetP(i,R.J.data);

#ifdef CG_CHECK_MW_PERM
      R.P0=P;
#endif
   }

   r2/=nu;
   if (sqrt(r2)>CG_EPS2) wblog(FL,"WRN %s() got r=%.3g",FCT,sqrt(r2));

   z2=JJ; z2.groupRecs(P,dJ);
   if (dJ.anyGT(1)) {
      if (iOM) {
         unsigned d, *iom; PERM_T *p=P.data;
         iOM->init(P.len,2);

         for (i=0; i<dJ.len; ++i) { d=dJ[i];
         for (j=0; j<d; ++j, ++p) { iom=iOM->rec(*p); iom[0]=j; iom[1]=d; }}
      }
      else {
         wblog(FL,"NB! %s() "
           "got outer multiplicity (OM<=%d)",FCT,dJ.max());
      }
#ifdef WB_SPARSE_CLOCK
      wbc_sparse_gss.stop();
#endif
      return r2;
   }
   else if (iOM) iOM->init();

#ifdef WB_SPARSE_CLOCK
   wbc_sparse_gss.stop();
#endif
   return r2;
};


template <class TQ, class TD>
void genRG_struct<TQ,TD>::put(const char *F, int L,
   const char *vname, const char *ws
) const {
   mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);
   
   if (i) wblog(F_L,
      "ERR failed to write variable `%s' (%d)",vname,i);
   if (F) wblog(F,L,"I/O putting variable '%s' to %s",vname,ws);

   mxDestroyArray(a);
};


template <class TQ, class TD>
mxArray* genRG_struct<TQ,TD>::toMx() const {

   unsigned i=0, n=RSet.size();
   const char* fields[]={"type","info","store"};
   mxArray *S, *a;

   if (!q.isKnown()) wblog(FL,"ERR %s() "
      "got invalid sym=%s (size=%d) !?",FCT,q.toStr().data,RSet.size());

   S=mxCreateStructMatrix(1,1,3,fields);
   a=MXPut().add(DZ,"DZ").add(CR,"CR").add(qs,"qdef").toMx();




   mxSetFieldByNumber(S,0, 0, q.toMx());
   mxSetFieldByNumber(S,0, 1, a);

   if (n) { i=0;
      for (crMAP I=RSet.begin(); I!=RSet.end(); ++I, ++i) {
         if (i==0) { a=I->second.mxCreateStruct(n,1); }
         I->second.add2MxStruct(a,i);
      }
      mxSetFieldByNumber(S,0, 2, a);
   }

   return S;
};



template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::initBase_SUN(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i,j, l, n=0, r=q.sub, N=r+1;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (N<2 || N>10) wblog(FL,"ERR %s() invalid SU(%d) !?",FCT,N);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (i=1; i<N; ++i) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];

       P.initz(N,N,1); P.setRec(0, i-1,i, 1.);

       Z.initDiag(N); Z[i]=-int(i); for (j=0; j<i; ++j) Z[j]=1; 

   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(N,r); X.initIdentity(N);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=N) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d) !??",FCT,q.toStr().data,R.Sz[i].isDiag(),sz.len,N);
      for (j=0; j<N; ++j) R.Z(j,i)=double(sz[j]);
   }

   findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() possibly already called earlier (%d)",
       FCT,R0.Sp.len,R0.Sz.len);
   R.save2(R0);

   int ioR=gStore.save_RSet(FL,R0,'q');

   
   if (CG_VERBOSE>4) {
      wblog(PF_L,
         "%s defining representation for %s: (%s)",
         ioR ? (ioR>1 ? "ok." : "[+]  W") : "[+]  w",
         R0.q.toStr().data, R0.J.toStr().data
      ); doflush();
   }



   unsigned m=0;
   n=MAX(3U,2*N); if (n>3) { --n; }; if (n>4) { --n; };
   for (i=0; i<n; ++i) {
      if (int(l=getTensorProdReps(R0.J))>0) m+=l;
   }; if (m) wblog(PFL,"generated %d initial multiplets (%d)",m,n);




   return *this;
};


template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::initBase_SpN(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i,j, r=q.sub, D=2*r;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (D<2 || D>10) wblog(FL,"ERR %s() invalid Sp(%d) !?",FCT,D);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (i=1; i<=r; i++) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];
       P.initz(D,D, i<r ? 2:1); Z.initDiag(D);

       if (i<r) {
          P.setRec(0, i-1, i, 1.);
          P.setRec(1, 2*r-i-1, 2*r-i, +1.);

          Z[i]=-int(i); for (j=0; j<i; ++j) Z[j]=1;
       }
       else {
          P.setRec(0, i-1, i, 1.);
          for (j=0; j<r; ++j) { Z[j]=1; }
       }

       for (j=0; j<r; ++j) Z[r+j]=-Z[r-1-j];
   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(D,r); X.initIdentity(D);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=D) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d) !??",FCT,q.toStr().data,R.Sz[i].isDiag(),sz.len,D);
      for (j=0; j<D; ++j) R.Z(j,i)=sz[j];
   }

   R.Sort();

   findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() possibly already called earlier (%d)",
       FCT,R0.Sp.len,R0.Sz.len);
   R.save2(R0);

   int ioR=gStore.save_RSet(FL,R0,'q');

   
   if (CG_VERBOSE>4) {
      wblog(PF_L,
         "%s defining representation for %s: (%s)",
         ioR ? (ioR>1 ? "ok." : "[+]  W") : "[+]  w",
         R0.q.toStr().data, R0.J.toStr().data
      ); doflush();
   }



   unsigned m=0, n=MAX(3U,r); int l;
   for (i=0; i<n; ++i) {
      if ((l=getTensorProdReps(R0.J))>0) m+=l;
   }; if (m) wblog(PFL,"generated %d initial multiplets (%d)",m,n);




   return *this;
};


#endif

