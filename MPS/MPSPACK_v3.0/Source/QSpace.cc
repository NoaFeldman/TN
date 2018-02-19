#ifndef __WB_QSPACE_COL_MAJOR_CC__
#define __WB_QSPACE_COL_MAJOR_CC__

// ----------------------------------------------------------------- //
// QSPACE
// ----------------------------------------------------------------- //
// tags: InitITags

template<class TQ, class TD>
void QSpace<TQ,TD>::init_itags(
   const char *F, int L, const char *s, unsigned k
){
   if (QIDX.dim2 && QDIM) {
      unsigned l,r=QIDX.dim2/QDIM;
      if (int(k)<0 && r) { k=-r; }

      l=itags.Init(F,L,s,k);

      if (QDIM) { if (r && r!=l) wblog(F_L,
         "ERR %s() mismatch in number of itags (%d/%d)",
          FCT,r,QIDX.dim2/QDIM); 
      }
      else if (r || QIDX.dim1 || QIDX.dim2) wblog(FL,
         "ERR %s() got invalid labels `%s' having (%dx%d;%d)",
          FCT,IT2STR__,QIDX.dim1,QIDX.dim2,QDIM
      ); 
   }
};



template<class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::SetTag(
   unsigned k, const iTag &b,
   char lflag
){

   if ((--k)>=itags.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,itags.len);
   iTag &a=itags[k]; bool conj=a.isConj();

   if (lflag) { a=b; if (conj) a.Conj(); }
   else {
      if (conj ^ b.isConj()) wblog(FL,
         "ERR %s() got conj-flag mismatch (%s -> %s)",
         FCT, a.toStr().data, b.toStr().data);
      a=b;
   }
   return *this;
};


template<class TQ, class TD>
void QSpace<TQ,TD>::initOType(const char *F, int L, const mxArray* a){

   if (!a) { otype=QS_NONE; return; }
   else {
      unsigned i=0; wbstring s(F,L,a);
      if (!s.data || !s.data[0]) { otype=QS_NONE; return; }

      for (i=1; i<QS_NUM_TYPES; ++i)
      if (!strcmp(s.data,QS_STR[i])) { otype=QS_TYPES(i); return; }

      if (i==QS_NUM_TYPES) wblog(FL,
      "ERR invalid object type '%s'",s.data);
   }
};


template<class TQ, class TD>
wbstring QSpace<TQ,TD>::otype2Str(const char *fmt) const {

   const char *s=0;

   switch (otype) {
     case QS_NONE     :
     case QS_AMATRIX  :
     case QS_OPERATOR : s=QS_STR[otype]; break;
     default : 
        sprintf(str,"(invalid: %d)",otype);
        return str;
   }
   if (fmt && fmt[0] && s && s[0]) {
      sprintf(str,fmt,s); return str;
   }
   return s;
};


template<class TQ, class TD>
template<class TDB>
unsigned QSpace<TQ,TD>::matchITags(
   const char *F, int L,
   const QSpace<TQ,TDB> &B, ctrIdx &ia, ctrIdx &ib
 ) const {

   unsigned n=0, ra=rank(FL), rb=B.rank(FL);

   if (itags.len!=ra || B.itags.len!=rb) wblog(F_L,
      "ERR %s() got empty or invalid set\n"
      "'%s' (%d/%d) <> '%s' (%d/%d)",FCT,IT2STR__,itags.len,ra,
      IT2STR(B), B.itags.len, rb
   );

   n=itags.match(B.itags,ia,ib);

   if (!n) { if (F) wblog(F,L,
     "ERR %s() got empty match\n"
     "'%s' <> '%s'", FCT,IT2STR__, IT2STR(B)
   ); }

   return ia.len;
};


template<class TQ, class TD>
void QSpace<TQ,TD>::PrependSingletons(unsigned R) {

    unsigned r; int ir=-1;

    if (!isConsistent(ir)) wberror(FL,str); r=(unsigned)ir;
    if (r>=R) { 
       if (r==R) return; else wblog(FL,
       "ERR PrependSingletons - index reduction !?? (%d->%d)", r,R);
    }
    else {
       unsigned i,s=QIDX.dim2*sizeof(TQ);
       wbMatrix<TQ> Q0(QIDX);
       TQ *q, *q0=Q0.data;

       if (itags.len) wblog(FL,
          "ERR %s() also need to adjust itags!  (%d/%d-%d)",
          FCT,itags.len,r,R);

       QIDX.init(QIDX.dim1, R*QDIM); q=QIDX.data+(R-r)*QDIM;
       for (i=0; i<QIDX.dim1; ++i) {
          memcpy(q,q0,s); q+=QIDX.dim2; q0+=Q0.dim2;
          DATA[i]->prependSingletons(R);
       }
    }
};


template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::setupCGR(const char *F, int L) {

    if (QIDX.dim1 && QIDX.dim2) {
    if (QDIM==0 || QIDX.dim2%QDIM) wblog(F_L,
       "ERR invalid QIDX (%dx%d; %d)",QIDX.dim1,QIDX.dim2,QDIM);
    if (qtype.len && qtype.Qlen()!=QDIM) wblog(F_L,
       "ERR qtype inconsistency (%d/%d)",QDIM,qtype.Qlen());
    }

    if (QIDX.dim1 && qtype.len) {
       CGR.init(QIDX.dim1,qtype.len);
    }
    else CGR.init();

    return *this;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::init2ref(const TD *v){
   unsigned i,s, n=DATA.len, N=0;
   
   for (i=0; i<n; ++i, v+=s) { s=DATA[i]->SIZE.prod(); N+=s;
      DATA[i]->init2ref(v);
   }
   return N;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::initIdentity(
   const QSpace<TQ,TD> **FF, unsigned n, char dflag
){
   unsigned i,d,l;
   wbvector<unsigned> D;
   wbperm P;

   wbvector< wbMatrix<TQ> > QQ;
   wbMatrix<TQ> Q1;

   wbvector< wbMatrix<unsigned> > SS;
   wbMatrix<unsigned> S1;

   clearQSpace(); if (n==0) return;

   QQ.initDef(n);
   SS.initDef(n);
   
   for (l=i=0; i<n; ++i) {
      const QSpace<TQ,TD> &F=*FF[i];
      if (!F.isEmpty() && !F.isConsistentR(2)) wberror(FL,str);
      if (i && !F.sameType(*FF[0])) wblog(FL,"ERR q-type mismatch");
      QQ[l]=F.QIDX;
      F.getSize(SS[l]); ++l;
   }

   if (l==0) return;

   QQ.len=l; Q1.CAT(1,QQ); d=Q1.dim2/2;
   SS.len=l; S1.CAT(1,SS); 

   if (S1.dim2==2) {
      Q1.Reshape(Q1.dim1*2,Q1.dim2/2);
      S1.Reshape(S1.dim1*2,S1.dim2/2);
   }
   else if (S1.dim2==1) {
      wbMatrix<TQ> QX; Q1.save2(QX);
      wbMatrix<TQ> Q2; QX.getBlock(0,d,Q1); QX.getBlock(1,d,Q2);
      if (Q1!=Q2) wblog(FL,
     "ERR diagonal operators not symmetric in Q");
   }

   if (Q1.dim1!=S1.dim1 || S1.dim2!=1) wblog(FL,
      "ERR size inconsistency (%d/%d,%d)",Q1.dim1,S1.dim1,S1.dim2);

   Q1.groupRecs(P,D);

   QDIM=Q1.dim2; QIDX.Cat(2,Q1,Q1);
   setupDATA();

   for (l=i=0; i<DATA.len; ++i, l+=d) { d=D[i];
      DATA[i]->initIdentity(S1[P[l]],dflag);
      if (!uniformRange(S1.data+P[l],d)) wblog(FL,
      "ERR size values inconsistent (not the same)");
   }

   qtype=FF[0]->qtype;

   initIdentityCGS();
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::initIdentity(
   const QSpace<TQ,TD> &A,
   char dflag
){
   wbMatrix<TQ> Q1;
   wbperm P;
   wbvector<INDEX_T> D;
   unsigned i,d,l;

   if (dflag && !strchr("rc",dflag))
   wblog(FL,"ERR %s() invalid dflag=%c<%d>",FCT,dflag,dflag);

   clearQSpace();
   
   if (A.isEmpty()) return;
   if (!A.isConsistentR(2)) wberror(FL,str);

   if (!A.isQSym(0,0,dflag)) { wblog(FL,
      "%s() requires symmetric QIDX\n%s",FCT,str);
       return;
   }

   A.QIDX.getBlock(0,A.QDIM,Q1);
   Q1.groupRecs(P,D);

   qtype=A.qtype;
   QDIM=Q1.dim2; QIDX.Cat(2,Q1,Q1); setupDATA();

   for (l=i=0; i<DATA.len; ++i, l+=D[i]) {
      d=A.DATA[P[l]]->SIZE.max();
      DATA[i]->initIdentity(d,dflag);
   }

   initIdentityCGS();
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::initIdentity(const wbMatrix<TQ> &Q) {
   clearQSpace();

   QDIM=Q.dim2; QIDX.Cat(2,Q,Q); setupDATA();
   for (unsigned i=0; i<DATA.len; ++i) {
      DATA[i]->init(1,1);
      DATA[i]->data[0]=1;
   }
}


template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initDiagonal(
   const wbMatrix<TQ> &Q,
   const wbvector<unsigned> &D,
   const wbvector<TD> &E,
   char rc
){
   unsigned i,j,d,l;

   if (D.len!=Q.dim1 || D.sum()!=E.len) wblog(FL,
      "ERR %s() severe size mismatch (%d==%d; %d==%d)",
       FCT, Q.dim1, D.len, D.sum(), E.len);
   if (rc && !strchr("rc",rc)) wblog(FL,
      "ERR %s() invalid rcflag=>%s<",FCT,rc
   );

   clearQSpace(); if (Q.isEmpty()) { return *this; }

   QIDX.Cat(2,Q,Q); QDIM=Q.dim2; setupDATA();

   for (l=i=0; i<D.len; ++i, l+=d) {
      wbarray<TD> &A=(*DATA[i]); d=D[i];

      if (rc) {
         if (rc=='r') A.init(1,d); else A.init(d,1);
         memcpy(A.data, E.data+l, d*sizeof(TD));
      }
      else {
         A.init(d,d);
         for (j=0; j<d; ++j) A.data[j+j*d]=E[l+j];
      }
   }
   return *this;
};

 
template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCGS(
   const char *F, int L
){
    unsigned i,j, r=2;
    wbvector<unsigned> S(r);

    if (!QIDX.dim1 || !QIDX.dim2) {
       CGR.init(0,qtype.len); return *this;
    }

    if (QDIM==0 || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);
    if (isref) wblog(FL,"ERR %s() got isref (bailing out)",FCT);
    if (QIDX.dim2/QDIM!=r) wblog(F_L,
       "ERR %s() for rank-%d objects only (%dx%d/%d; %s)",
        FCT,r,QIDX.dim1,QIDX.dim2,QDIM,qStr().data
    );

    if (qtype.len) {
       CGR.init(QIDX.dim1,qtype.len);
       wbvector<unsigned> qdc; qtype.Qpos(qdc);

       if (qtype.Qlen()!=QDIM)
          wblog(F_L,"ERR %s() qtype inconsistency (%s: %d/%d)",
          FCT, qStr().data,qtype.Qlen(),QDIM
       );

       for (i=0; i<CGR.dim1; ++i) { const TQ *qi=QIDX.rec(i);
       for (j=0; j<CGR.dim2; ++j) {
          CGR(i,j).initIdentityR(FL, qtype[j], qi+qdc[j]);
       }}
    }
    else CGR.init();

    if (!itags.len) itags.init_qdir("+-");
    else if (!itags.isOp(2)) wblog(FL,
       "ERR %s() invalid itags %s",FCT,itags.toStr().data
    );

    return *this;
};

 
template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCGS(
   const char *F, int L, const QSpace<TQ,TD> &H
){
    unsigned i,j, r=2;
    wbvector<unsigned> S(r);
    wbindex I1,I2;

    if (!QIDX.dim1 || !QIDX.dim2) {
       CGR.init(0,qtype.len); return *this;
    }

    if (QDIM==0 || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);
    if (qtype.isEmpty()) wblog(F_L,"ERR %s() got empty qtype",FCT);
    if (isref) wblog(FL,"ERR %s() got isref (bailing out)",FCT);
    if (qtype.Qlen()!=QDIM) wblog(F_L,
       "ERR %s() qtype inconsistency (%s: %d/%d)",
        FCT, qStr().data,qtype.Qlen(),QDIM); 
    if (QIDX.dim2/QDIM!=r) wblog(F_L,
       "ERR %s() for rank-%d objects only (%dx%d/%d; %s)",
        FCT,r,QIDX.dim1,QIDX.dim2,QDIM,qStr().data);
    if (qtype!=H.qtype) wblog(F_L,"ERR %s() mismatch in qtype (%s <> %s)",
        FCT, qStr().data, H.qStr().data); 

    i=(unsigned)matchIndex(QIDX, H.QIDX, I1, I2);
    if (i || QIDX.dim1!=I1.len) wblog(FL,
       "ERR %s() QIDX not fully contained in reference space (%d,%d)",
       FCT, I1.len, QIDX.dim1
    );

    CGR.init(QIDX.dim1,qtype.len);

    for (i=0; i<CGR.dim1; ++i) {
    for (j=0; j<CGR.dim2; ++j) { CGR(I1[i],j)=H.CGR(I2[i],j); }}

    return *this;
};


template <class TQ, class TD>
void QSpace<TQ,TD>::Diag2Vec(char tflag) {

   if (!isDiagMatrix()) wblog(FL,
   "ERR %s() called for non-diagonal object",FCT);

   for (unsigned i=0; i<DATA.len; ++i)
   DATA[i]->Diag2Vec(tflag);
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::ExpandDiagonal(
   unsigned i1,
   unsigned i2
){
   unsigned i,j1,j2, n=QDIM*sizeof(TQ), dflag=0;
   int r=-1;

   isConsistent(r,FL);
   if (!i1 || !i2 || int(i1)>r || int(i2)>r || i1==i2) wblog(FL,
      "ERR %s() invalid index (%d %d; %d)",FCT,i1,i2,r);

   i1--; i2--;

   j1=i1*QDIM;
   j2=i2*QDIM;

   for (i=0; i<DATA.len; ++i) {
      if (memcmp(QIDX.ref(i,j1), QIDX.ref(i,j2),n)) {
         continue;
      }

      const wbvector<unsigned> &S=DATA[i]->SIZE;


      if (S[i1]==S[i2]) continue;

      if (S[i1]!=1) {
         if (dflag) { if (dflag!='c') wblog(FL,
            "ERR %s() expecting `col'\n%s (%d,%d)",
             FCT, S.toStrD().data, i1+1, i2+1);
         } else dflag='c';
      }
      else
      if (S[i2]!=1) {
         if (dflag) { if (dflag!='r') wblog(FL,
            "ERR %s() expecting `row'\n%s (%d,%d)",
             FCT, S.toStrD().data, i1+1, i2+1);
         } else dflag='r';
      }

      DATA[i]->ExpandDiagonal(i1,i2);
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::Expand2Projector(const TD eps){

   unsigned i,dflag=0;
   wbvector<char> mark(DATA.len);

   wblog(FL,"TST check *** col-major ***");

   if (!isBlockDiagMatrix('D')) wblog(FL,
   "ERR Expand2Projector assumes block-diagonal structure !??");

   for (i=0; i<DATA.len; ++i) {
      const unsigned *S = DATA[i]->SIZE.data;

      if (S[0]!=1 && S[1]!=1) wblog(FL,
         "ERR ExpandDiagonal - expecting row/col\n%s: %dx%d",
          QIDX.rec2Str(i).data, S[0],S[1]
      );


      if (S[0]!=1) { if (dflag) { if (dflag!='c') wblog(FL,
         "ERR ExpandDiagonal - expecting `col' (%dx%d)",S[0],S[1]); }
          else dflag='c';
      } else
      if (S[1]!=1) { if (dflag) { if (dflag!='r') wblog(FL,
         "ERR ExpandDiagonal - expecting `row' (%dx%d)",S[0],S[1]); }
          else dflag='r';
      }

      DATA[i]->Expand2Projector(eps);
      if (DATA[i]->isEmpty()) mark[i]=1;
   }

   if (mark.anyUnequal(0)) {
      QSpace<TQ,TD> X; wbindex I;
      mark.find(I); getSub(I,X); X.save2(*this);
   }
}


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::init2DiffOp(
   const QSpace<TQ,TD> &A, unsigned ia,
   const QSpace<TQ,TD> &B, unsigned ib
){
   wbvector<INDEX_T> d1,d2,D;
   wbMatrix<TQ> Qa,Qb,Qx;
   wbindex Ia,I0;
   wbperm P1,P2;

   A.getQsub(ia,Qa).groupRecs(P1,d1); I0.initGroup0(P1,d1);
   B.getQsub(ib,Qb).groupRecs(P2,d2);

   Qa.getDiffSorted(Qb,Ia);
   if (!Ia.len) { init(); return *this; }

   unsigned i,j,d; wbvector<TD> E;

   Qa.getRecs(Ia,Qx); I0.Select(Ia);

   D.init(Ia.len);
   for (i=0; i<I0.len; ++i) { D[i]=A.DATA[I0[i]]->SIZE.elx(ia,1); }
   E.init(D.sum());

   initDiagonal(Qx,D,E,0);

   qtype=A.qtype; if (qtype.len) {
      wbvector<unsigned> qdc; qtype.Qpos(qdc);

      setupCGR();
      if (CGR.dim2!=A.CGR.dim2) wblog(FL,
         "ERR %s() size mismatch (%d/%d)",FCT,CGR.dim2,A.CGR.dim2);
      for (i=0; i<I0.len; ++i) {
         const CRef<TQ> *cg=A.CGR.rec(I0[i]);
         const TQ *qi = A.QIDX.ref(I0[i],ia*A.QDIM);
         for (j=0; j<CGR.dim2; ++j) {
            d=cg[j].Size(ia);
            CGR(i,j).initIdentityR(FL, qtype[j], qi+qdc[j], d);
         }
      }
   }
   else if (B.qtype.len) wblog(FL,
      "ERR %s() qtype mismatch (%s,%s)",FCT,
       A.qStr().data, B.qStr().data
   );

   itags.init_qdir("+-");

   return *this;
};


template <class TQ, class TD> inline
char QSpace<TQ,TD>::gotCGS(const char *F, int L) const {

   if (CGR.isEmpty()) {
      if (qtype.allAbelian() 
         || (!QIDX.dim1 && DATA.len==0)
         || (!QIDX.dim2 && DATA.len==1)
         ){ return 0; }
      wblog(F_L,"ERR %s() no CGR data for %s",FCT,qStr().data); 
   }
   if (qtype.isEmpty() || CGR.dim2!=qtype.len) wblog(F_L,
      "ERR %s() got invalid qtype (%d/%d)",FCT,qtype.len,CGR.dim2
   );
   if (CGR.dim1!=QIDX.dim1 ||
       CGR.dim2!=qtype.len || DATA.len!=QIDX.dim1) wblog(FL,
      "ERR %s() severe QSpace inconsistency\n(CGR: %dx%d <> %dx%d; %d)",
       FCT, CGR.dim1, CGR.dim2, QIDX.dim1, qtype.len, DATA.len
   );
   if (!QDIM || QIDX.dim2%QDIM || QDIM!=qtype.Qlen()) wblog(F_L,
      "ERR %s() invalid QDIM (QIDX.dim2: %d @ %d/%d)",
      FCT, QIDX.dim2, QDIM, qtype.Qlen()
   );
 
   unsigned i,j;
   char isa=isAbelian(), m=-1;
   QDir qdir(itags);

   if (!isa) {
      m=1+qtype.hasMP(rank(FL));
   }

   try {
      for (i=0; i<CGR.dim1; ++i)
      for (j=0; j<CGR.dim2; ++j) {
         CGR(i,j).check(FL,qtype[j],qdir);
      }
   } catch (...) {
      wblog(FL,"ERR %s() CGR(%d,%d)\nhaving %s [%dx%d]",
      FCT,i+1,j+1, qStr().data, CGR.dim1,CGR.dim2); 
   }

   return m;
};


template <class TQ, class TD> inline
char QSpace<TQ,TD>::gotCGX(const char *F, int L) const {

   unsigned i,j;
   for (j=0; j<CGR.dim2; ++j)
   for (i=0; i<CGR.dim1; ++i) {
      if (CGR(i,j).cgp.len || CGR(i,j).conj) {
         if (!CGR(i,j).cgr || CGR(i,j).cgr->cgd.isScalar()) wblog(F_L,
            "ERR %s() unexpected CRef(%d,%d) data (%s)",
            FCT,i+1,j+1,qStr().data);
         return 1;
      }
   }
   return 0;
};


template <class TQ, class TD>
template <class TDB>
bool QSpace<TQ,TD>::sameType(
   const QSpace<TQ,TDB> &B, int r, const char *istr
 ) const {

   if (istr) str[0]=0;
   if ((void*)this==(void*)&B) { if (r<0) return 1; }
   else {
      if (QIDX.dim2!=B.QIDX.dim2 || QDIM!=B.QDIM) {
         if (istr) sprintf(str,"%s: %s QIDX: %d/%d, QDIM: %d/%d",
         shortFL(FL), istr, QIDX.dim2, B.QIDX.dim2, QDIM, B.QDIM);
         return 0;
      }
      if (!qtype.sameType(B.qtype)) {
         if (istr) sprintf(str,"%s: %s qtype: %s/%s",
         shortFL(FL),istr, qStr().data, B.qStr().data);
         return 0;
      }
      if (QDIM!=qtype.Qlen()) {
         if (istr) { sprintf(str,"%s: %s QDIM/qtype: %d/%d",
             shortFL(FL),istr,QDIM,qtype.Qlen());
             return 0;
         }
         else wblog(FL,"ERR qtype not compatible with QDIM (%d/%d)",
         QDIM,qtype.Qlen());
      }
   };

   if (r<0) return 1;

   if (r==0) {
      if (QIDX.dim2) { if (istr) sprintf(str,
         "%s: %s rank = %d/0",shortFL(FL),istr,QIDX.dim2);
          return 0;
      }
      else return 1;
   }
   else if (QDIM==0 || QIDX.dim2/QDIM!=(unsigned)r) {
      if (istr) sprintf(str,"%s: %s rank = %d/%d",
          shortFL(FL),istr,QIDX.dim2/QDIM,r);
      return 0;
   }

   return 1;
};


template <class TQ, class TD>
bool QSpace<TQ,TD>::sameQDir(
   const char *F, int L, const QSpace<TQ,TD> &B
 ) const {

   if (qtype.len!=B.qtype.len || B.qtype.len!=B.CGR.dim2 || 
       qtype.len!=CGR.dim2) wblog(F_L,"ERR %s() qtype mismatch "
      "(%s <> %s)", FCT,qStr().data,B.qStr().data
   );

   if (!CGR.dim1 || !CGR.dim2) return 1;
   for (unsigned i=0, n=qtype.len; i<n; ++i) {
      if (!CGR(0,i).sameQDir(0,0,B.CGR(0,i))) {
         if (F) wblog(F_L,"ERR %s() got qdir mismatch\n"
            "(%s => %s <> %s => %s)", FCT,
              CGR(0,i).qdir2Str('V').data,   CGR(0,i).qdir2Str().data,
            B.CGR(0,i).qdir2Str('V').data, B.CGR(0,i).qdir2Str().data
         );
         return 0;
      }
   }

   return 1;
};


template <class TQ, class TD>
template <class TDB>
void QSpace<TQ,TD>::checkQ(
   const char *F, int L, const QSpace<TQ,TDB> &B, char aflag
 ) const {

   if (QDIM!=B.QDIM) wblog(F_L,
      "ERR %s() QDIM inconsistency (%d/%d)",FCT,QDIM,B.QDIM);

   if (qtype.len && QDIM!=qtype.Qlen()) wblog(F_L,
      "ERR %s() length inconsistency in qtype (%d/%d)",
       FCT, qtype.Qlen(), QDIM);

   if (qtype!=B.qtype) {
     if (aflag && isAbelian() && B.isAbelian() &&
        (!qtype.len || !B.qtype.len)) { aflag='A'; }
     if (aflag!='A') wblog(F_L,
        "ERR %s() qtype inconsistency ('%s' <> '%s')",
         FCT,qStr().data,B.qStr().data
     ); 
   }



   if (aflag!='A') {
      if (CGR.dim2!=B.CGR.dim2) wblog(F_L,
         "ERR %s() CGR inconsistency: %dx%d <> %dx%d",
         FCT, CGR.dim1, CGR.dim2, B.CGR.dim1, B.CGR.dim2
      );
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::ExpandQ(
   const wbMatrix<TQ> &QS,
   const wbvector<unsigned> &d,
   unsigned k,
   double ra,
   wbindex *I,
   wbvector<char> *mQ
){
   unsigned i,j,l,m,dj, r=rank(FL), s=QDIM*sizeof(TQ);
   wbMatrix<TQ> QQ,QT,QT2;
   wbindex I0,II;
   QSpace<TQ,TD> AX;
   wbvector<unsigned> S;
   wbvector<char> mark;
   TQ *q, *qref;
   double nrm1=0., nrm2=0.;

   if (k>=r) wblog(FL,
   "ERR %s - dimension out of range (%d/%d)",FCT,k,r);
   if (QS.dim1!=d.len || QS.dim2!=QDIM) wblog(FL,
   "ERR %s - dimension mismatch (%dx%d / %d,%d)", FCT,
    QS.dim1, QS.dim2, d.len, QDIM);

   getQsub(k,QQ); QQ.makeUnique();
   i=(unsigned)matchIndex(QQ, QS, I0, II);
   if (i || QQ.dim1!=I0.len) wblog(FL,
      "ERR %s - local Q[%d] not fully contained \nin input Q (%d,%d)",
      FCT, k+1, I0.len, QQ.dim1
   );

   if (mQ) { QIDX.blockSum(QDIM,QT); QT.makeUnique(); }

   QQ=QIDX; QQ.setBlock(k,QDIM,0);
   QQ.makeUnique(); m=QQ.dim1;

   QQ.Resize(QQ.dim1*QS.dim1, QQ.dim2);
   q=QQ.data+k*QDIM;

   for (i=0; i<QS.dim1; ++i) {
      qref=QS.rec(i);
      for (j=0; j<m; ++j, q+=QQ.dim2) {
         if (i) QQ.recSet(i*m+j,j);
         memcpy(q,qref,s);
      }
   }

   if (mQ) {
      QQ.blockSum(QDIM,QT2);
      matchIndex(QT, QT2, I0, II);
      mQ->init(QT2.dim1);
      for (i=0; i<II.len; ++i) mQ->data[II[i]]++;
   }

   i=(unsigned)matchIndex(QIDX, QQ, I0, II);

   if (i || QIDX.dim1!=I0.len) { QS.Print("Q0"); wblog(FL,
   "ERR %s - input QSpace not unique (%d; %d/%d) !??", FCT,
    i, QIDX.dim1, I0.len); }

   save2(AX); QQ.save2(QIDX); setupDATA();
   mark.init(QIDX.dim1);

   for (i=0; i<I0.len; ++i) {
      j=II[i]; mark[j]++;
      AX.DATA[I0[i]]->save2(*DATA[j]);
   }

   for (i=0; i<QIDX.dim1; ++i) { j=i/m; dj=d[j];
      if (mark[i]) {
         if (mark[i]!=1) wblog(FL,"WRN mark[%d]=%d !??", i, mark[i]);
         if (DATA[i]->SIZE[k]!=dj) wblog(FL,
         "ERR local dimension mismatch (%d/%d)",DATA[i]->SIZE[k],dj);
         nrm1+=DATA[i]->norm2();
      }
      else {
         for (r=i%m, l=0; l<QS.dim1; ++l, r+=m)
         if (mark[r]) { S=DATA[r]->SIZE; break; }
         if (l==QS.dim1) wblog(FL,"ERR Failed to find correct size !??");

         S[k]=dj;
         DATA[i]->init(S); mark[i]=-1;

         DATA[i]->setRand();
         nrm2+=DATA[i]->norm2();
      }
   }

   if (nrm2) {
      wblog(FL,"TST change rand norm %.3g -> %.3g", sqrt(nrm2), sqrt(nrm1)*ra);
      nrm2=ra*sqrt(nrm1/nrm2);

      for (i=0; i<QIDX.dim1; ++i) {
         if (mark[i]>0) continue;
         (*DATA[i]) *= nrm2;
      }
   }

   if (I) {
      wbperm P; I0.Sort(P);
      II.select(P,*I);
   }
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::recSetQ(unsigned k, unsigned i) {

   if (DATA.len!=QIDX.dim1 || (CGR.dim2 && CGR.dim1!=QIDX.dim1)) {
      wblog(FL,"ERR severe QSpace inconsistency (%d/%d/%d)",
      CGR.dim1, DATA.len, CGR.dim1);
   }

   if (i==k) return;
   if (i>=QIDX.dim1 || k>=QIDX.dim1) wblog(FL,
      "ERR index out of bounds ({%d,%d}/%d)",i,k,QIDX.dim1);

   QIDX.recSet(k,i);

   WB_DELETE_1(DATA[k]);
   DATA[k]=DATA[i]; DATA[i]=NULL;

   if (CGR.dim2) {
   for (unsigned j=0; j<CGR.dim2; ++j) {
      CGR(k,j)=CGR(i,j); CGR(i,j).init();
   }}
};


template <class TQ, class TD>
int QSpace<TQ,TD>::getQOverlap(
   const QSpace<TQ,TD> &B,
   wbMatrix<TQ> &QQ,
   wbMatrix<int> &IQ,
   char olonly
) const{

   unsigned i,i1,i2,d,k,l, m=QIDX.dim1, n=QIDX.dim2;
   wbvector<INDEX_T> D;
   wbperm P;

   if (B.QIDX.dim1==0) QQ=  QIDX; else
   if (  QIDX.dim1==0) QQ=B.QIDX;
   else {
      if (n!=B.QIDX.dim2 || rank(FL)!=B.rank(FL)) {
         sprintf(str,"%s incompatible QSpaces (QDIM=%d,%d; R=%d,%d).",
         shortFL(FL), QDIM, B.QDIM, rank(FL), B.rank(FL)); return 1;
      }
      QQ.Cat(1, QIDX, B.QIDX);
   }

   QQ.groupRecs(P,D);
   IQ.init(QQ.dim1,2);

   for (k=l=i=0; i<QQ.dim1; ++i, l+=d) {
      d=D[i];
      if (d==1) {
        if (!olonly) {
          if (k<i) QQ.recSet(k,i);
          i1=P[l]; if (i1<m)
               { IQ(k,0)=i1; IQ(k,1)=-1;   }
          else { IQ(k,0)=-1; IQ(k,1)=i1-m; }
          ++k;
        }
      }
      else if (d==2) {
          if (k<i) QQ.recSet(k,i);
          i1=P[l]; i2=P[l+1];

          if (i1<m && i2>=m) { i2-=m; } else
          if (i2<m && i1>=m) { SWAP(i1,i2); i2-=m; }
          else {
             sprintf(str,"%s:%d QSpaces have non-unique QIDX ???",
             FL); return 1;
          }

          IQ(k,0)=i1; IQ(k,1)=i2; 
          ++k;
      }
      else {
         sprintf(str,"%s:%d QSpaces have non-unique QIDX (%d) ???",
         FL, d); return 1;
      }
   }

   if (k<i) {
      QQ.Resize(k,QQ.dim2);
      IQ.Resize(k,IQ.dim2);
   }

   return 0;
}


template <class TQ, class TD>
bool QSpace<TQ,TD>::findDimQ(const TQ* q0, unsigned &d) const {
 
   unsigned i,j,r,n, s=QDIM*sizeof(TQ);
   const TQ *q=QIDX.data;
   int ir=-1;

   if (!isConsistent(ir)) wberror(FL,str);
   if (QIDX.isEmpty()) return 0;

   r=rank();
   n=QIDX.dim1*r;

   for (i=0; i<n; ++i, q+=QDIM) { if (!memcmp(q,q0,s)) break; }
   if (i==n) {
      wbvector<TQ> x(QDIM,q0);
      wblog(FL,"TST %s() %d/%d (%d)",FCT,i,n,QDIM); 
      MXPut(FL,"i").add(*this,"A").add(x,"q");
      return 0;
   }

   j=i%r; i=i/r; d=DATA[i]->SIZE[j];
   return 1;
};


template <class TQ, class TD>
bool QSpace<TQ,TD>::findDimQ(const TQ* q, unsigned k, unsigned &d) const {
 
   unsigned i,r, n=QIDX.dim1, s=QDIM*sizeof(TQ);
   TQ *q0=QIDX.data+k*QDIM;
   int ir=-1;

   if (!isConsistent(ir)) wberror(FL,str); r=(unsigned)ir;
   if (QIDX.isEmpty()) return 0;

   for (i=0; i<n; ++i, q0+=QIDX.dim2) {
      if (!memcmp(q0,q,s)) break;
   }  if (i==n) return 0;

   d=DATA[i]->SIZE[k];
   return 1;
};



template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isConsistent(
  int &r0, const char *F, int L, char level
) const {

    int dr=0, r=rank(F,L); str[0]=0;

    if (isEmpty()) {
       if (r0<0) r0=0;
       return 1;
    };
    if (r<0 && !F) return 0;

    if (r0>0) { if (r0!=r && (r0!=2 || r!=1)) { sprintf(str,
       "rank mismatch (is %d, should be %d)", r, r0);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }} else r0=r;

    if (itags.len && int(itags.len)!=r) { sprintf(str,
       "QSpace::itags inconsistency (%ldx%ld/%d, %ld/%d)",
        QIDX.dim1, QIDX.dim2, QDIM, itags.len,r); 
        if (F) wblog(F,L,"ERR %s",str); 
        return 0;
    }

    if (r==3 && otype==QS_OPERATOR) { dr=1; r=2; }

    if (QIDX.dim1!=DATA.len) wblog(F_L,
       "ERR %s() severe QSpace inconsistency (%dx%d/%d, %d)",
        QIDX.dim1, QIDX.dim2, QDIM, DATA.len
    ); 

    gotCGS(F_L);

    if (isref) return 1;

    if (level) { r=MAX(2,r+dr);
    for (unsigned i=0; i<DATA.len; ++i) {
       if (!DATA[i]) { sprintf(str,"QSpace got null space DATA[%d]",i+1);
           if (F) wblog(F,L,"ERR %s",str); 
           return 0;
       }

       if (DATA[i]->SIZE.len!=unsigned(r)) {
           if (DATA[i]->SIZE.len==1) continue;

           sprintf(str,"QSpace inconsistency (%s)\n\n"
             "DATA[%d] has rank %ld (size %s) expecting rank %d",
              shortFL(FL), i+1, DATA[i]->SIZE.len,
              DATA[i]->sizeStr().data, r0
           );
           if (F) wblog(F,L,"ERR %s",str); 
           return 0;
       }

       if (r0>=2) continue;

       if ((r0==1 && (DATA[i]->SIZE[0]!=1 && DATA[i]->SIZE[1]!=1)) ||
           (r0==0 && (DATA[i]->SIZE[0]!=1 || DATA[i]->SIZE[1]!=1)) ) {
           sprintf(str,"invalid size %s for rank-%d object",
           DATA[i]->sizeStr().data, r0);
           if (F) wblog(F,L,"ERR %s",str); 
           return 0;
       }

       if (r0==0 && i) { sprintf(str,
          "rank-%d object can only have max 1 element (%d)", r0, i+1);
           if (F) wblog(F,L,"ERR %s",str); 
           return 0;
       }
    }}

    return 1;
};


template <class TQ, class TD>
int QSpace<TQ,TD>::isOperator(unsigned *r2, char xflag) const {

   unsigned r=rank(FL), dr=0;
   if (otype==QS_OPERATOR) { char gotc=gotCGS(FL);
      if ((gotc>0 && r!=3) || (gotc<=0 && r!=2 && r!=3)) wblog(FL,
         "ERR %s() got rank-%d for %s",FCT,r,otype2Str().data);
      dr=r-2; r=2;
   }

   if (r2) {
      if (int(*r2)>0 && r!=(*r2)) { return 0; }
      (*r2)=r+dr;
   }

   if (int(r)<=0) return -1;
   if (otype!=QS_NONE && otype!=QS_OPERATOR) return -2;

   return (r==2 || (xflag && r%2==0)) ? r : 0;
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isDiagBlock(unsigned k) const {

    unsigned s=QDIM*sizeof(TQ);
    const TQ *q=QIDX.rec(k);

    if (!isConsistentR(2)) wberror(FL,str); 
    if (k>=QIDX.dim1) wblog(FL,
       "ERR index out of bounds (%d/%d)",k,QIDX.dim1);

    if (memcmp(q, q+QDIM, s)) return 0;
    if (DATA[k]->SIZE[0]!=DATA[k]->SIZE[1]) return 0;

    return 1;
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isBlockDiagMatrix(
    const char *F, int L, char dflag) const {

    unsigned r=rank(F_L); if (r%2) return 0;

    unsigned i, n=QIDX.dim1,
       r2=r/2, m=r2*QDIM, sq=m*sizeof(TQ), s2=r2*sizeof(TQ);
    TQ *q=QIDX.data;

    if (dflag && r>2) { wblog(FL,
       "WRN %s() dflag=%d will be ignored (r=%d)",FCT,dflag,r);
       dflag=0;
    }

    for (i=0; i<n; ++i, q+=QIDX.dim2) {
       if (memcmp(q, q+m, sq)) {
          if (F) wblog(F_L,"--> %s(): non-matching Q(%d,:)",FCT,i+1); 
          return 0;
       }

       const wbvector<unsigned> &S=DATA[i]->SIZE;

       if (dflag) {
          if (S.len!=2 || (S[0]!=1 && S[1]!=1)){
             if (F) wblog(F_L,"--> %s(): data{%d} got size %s [dflag]",
                FCT, i+1, DATA[i]->sizeStr().data,dflag); 
             return 0;
          }
       }
       else {
          if (S.len!=r) wblog(FL,"ERR %s() got "
             "rank mismatch (%d; %d)",FCT,DATA[i]->SIZE.toStr().data,r);
          if (memcmp(S.data, S.data+r2, s2)!=0) {
             if (F) wblog(F_L,"--> %s(): data{%d} got size %s",FCT, i+1,
                DATA[i]->sizeStr().data); 
             return 0;
          }
       }
    }

    return 1;
}


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isDiagMatrix(const TD eps) const {

    if (!isBlockDiagMatrix()) return 0;

    for (unsigned i=0; i<DATA.len; ++i) {
       if (!(DATA[i]->isDiagMatrix(eps))) return 0;
    }

    if (gotCGS(FL)>0) {
       for (unsigned  n=CGR.numel(), i=0; i<n; ++i) {
          if (!CGR.data[i].isDiagCGC(eps)) return 0;
       }
    }

    return 1;
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isIdentityMatrix(const TD eps) const {

   if (!isBlockDiagMatrix(FL)) return 0;

   if (CGR.isEmpty() || qtype.allAbelian()) {
      for (unsigned i=0; i<DATA.len; ++i)
      if (!(DATA[i]->isIdentityMatrix(eps))) return 0;
   }
   else {
      if (CGR.dim1!=DATA.len) wblog(FL,"ERR %s() "
         "size mismatch (%d / %dx%d)",FCT,DATA.len,CGR.dim1,CGR.dim2);
      RTD x; TD dval; double cfac; unsigned i=0, j=0;

      for (; i<DATA.len; ++i) {
         for (cfac=1, j=0; j<CGR.dim2; ++j) {
            if (!CGR(i,j).isIdentityCGC(&x)) { return 0; }
            cfac*=double(x);
         }
         if (!cfac) { return 0; }

         dval=TD(1)/cfac;
         if (!(DATA[i]->isProptoId(dval,eps))) { return 0; }
      }
   }

   return 1;
};


template <class TQ, class TD> inline
char QSpace<TQ,TD>::hasIdentityCGS(const RTD eps) const {

    RTD x;
    for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
       if (!CGR.data[i].isIdentityCGC(&x,eps)) {
          return 0;
       }
    }

    return 1;
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isHConj(
   const char *F, int L, double eps, char vflag) const {

   bool cgflag=(gotCGS(FL)>0);

   if (cgflag) {
      unsigned r=rank(F,L);

      if (r<2 || (r!=3 && !isQSym(F,L)) || (r==3 && getDIM(2)>1)) {
         if (vflag || F) {
            if (r<2 || r==3) { sprintf(str,"%s got rank-%d QSpace (%s)",
               shortFL(FL),r, sizeStr('v').data); }
            else { sprintf(str,
               "%s got non-symmetry Q data (r=%d)", shortFL(FL),r); }
            if (F) wblog(F_L,"ERR %s() %s",FCT,str);
         }
         return 0;
      }

      QSpace<TQ,TD> X; wbperm P;
      double e=0;

      try {
         if (r!=3) 
              P.initTranspose(r);
         else P.initStr("2,1,3");
         permute(P,X).Conj(); if (r==3) X.ConjOpScalar();
         X-=(*this); e=X.norm();
      }
      catch (...) {
        if (r!=3) wblog(FL,
           "ERR %s() failed to subtract block-transpose QSpace\n"
           "ERR %s (p=[%s]; r=%d)",FCT,sizeStr('v').data,P.toStr().data,r);
        else wblog(FL,
           "ERR %s() failed to subtract transpose IROP\n"
           "ERR %s (p=[%s]; r=%d)",FCT,sizeStr('v').data,P.toStr().data,r
        );
      }
      if (e>eps) { double q=norm(); if (q>1) eps*=q;
         if (e>eps) {
            if (vflag || F) {
               sprintf(str,"%s got non-symmetric data (%g)",shortFL(FL),e);
               if (F) wblog(F,L,"ERR %s() %s",FCT);
            }
            return 0;
         }
      }
      return 1;
   }
   else {
      return isSym_aux(F,L,FCT,eps,'s',vflag);
   }
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isAHerm(
   const char *F, int L, double eps, char vflag) const {

   bool cgflag=(gotCGS(FL)>0);

   if (cgflag) {
      if (!isQSym(F,L)) {
         if (vflag) wblog(FL,"ERR %s() got size mismatch",FCT);
         return 0;
      }
      unsigned r=rank(FL); QSpace<TQ,TD> X;
      wbperm P; P.initTranspose(r); permute(P,X); X+=(*this);
      double e=X.norm();
      if (e>eps) { if (vflag) wblog(FL,
         "ERR %s() got non-symmetric data (%g)",FCT,e);
         return 0;
      }
      else return 1;
   }
   else {
      return isSym_aux(F,L,FCT,eps,'s',vflag);
   }
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isSym_aux(
  const char *F, int L, const char *fct,
  RTD eps, char symflag, char vflag
) const {

    wbMatrix<TQ> Q1,Q2;
    wbperm P1,P2,iP2,P;

    unsigned i,j,k=QDIM,K; int e=0, r=-1;
    gotCGS(FL);

    if (isEmpty()) return 1;
    str[0]=0; isConsistent(r,FL);

    if (r%2) {
       sprintf(str,"%s() requires even-rank (%d)",fct,r);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    K=(unsigned)r/2;

    Q1=QIDX; Q1.sortRecs(P1);

    P.init((unsigned)r).Rotate(K);
    QIDX.blockPermute(P,Q2); Q2.sortRecs(P2);

    if (Q1!=Q2) { 
       sprintf(str,"QIDX[:,1] does not match QIDX[:,2]");
       if (F || vflag) wblog(F_L,"ERR %s",str);
       return 0;
    }
    if (!Q1.isUniqueSorted()) {
       sprintf(str,"WRN QIDX is not unique!");
       if (F || vflag) wblog(F_L,"ERR %s",str);
       return 0;
    }

    P2.getIPerm(iP2);

    for (i=0; i<DATA.len; ++i) {
       j=P1[iP2[i]]; if (j<i) continue;

       const wbarray<TD> &di=(*DATA[i]);
       const wbarray<TD> &dj=(*DATA[j]);

       if (symflag=='s') {
          if (i==j)
               { if (!di.isHConj(   eps)) e=1; }
          else { if (!di.isHConj(dj,eps)) e=2; }
       }
       else if (symflag=='a') {
          if (i==j)
               { if (!di.isAHerm(   eps)) e=3; }
          else { if (!di.isAHerm(dj,eps)) e=4; }
       }
       else wblog(F,L,"ERR %s() invalid flag %c<%d>",FCT,symflag,symflag);

       for (k=0; !e && k<CGR.dim2; ++k) {
          if (CGR(i,k).cgr && !CGR(i,k).cgr->cgd.isScalar()) {
          const cdata__ &ci=(CGR(i,k).cgr->cgd);
          const cdata__ &cj=(CGR(j,k).cgr->cgd);

          if (symflag=='s') {
             if (i==j)
                  { if (!ci.isHConj(   eps)) { e=11; break; }}
             else { if (!ci.isHConj(cj,eps)) { e=12; break; }}
          }
          else {
             if (i==j)
                  { if (!ci.isAHerm(   eps)) { e=13; break; }}
             else { if (!ci.isAHerm(cj,eps)) { e=14; break; }}
          }
       }}

       if (!e && !CGR.isEmpty()) {
          cdata__ ci,cj; char gotc=0;
          
          for (k=0; k<CGR.dim2; ++k) { if (CGR(i,k).cgr) {
             if (gotc) {
                ci.Kron(CGR(i,k).cgr->cgd); if (i!=j) {
                cj.Kron(CGR(j,k).cgr->cgd); }
             }
             else { gotc=1;
                ci=(CGR(i,k).cgr->cgd); if (i!=j) {
                cj=(CGR(j,k).cgr->cgd); };
             }
          }}

          if (symflag=='s') {
              if (i==j)
                   { if (!ci.isHConj(   eps)) { e=21; break; }}
              else { if (!ci.isHConj(cj,eps)) { e=22; break; }}
          }
          else {
              if (i==j)
                   { if (!ci.isAHerm(   eps)) { e=23; break; }}
              else { if (!ci.isAHerm(cj,eps)) { e=24; break; }}
          }
       }

       if (e) {
          size_t l=0, n=128; char istr[n];
          if (e<10) { l=snprintf(istr,n,
             "%s %s() non-matching data{%d} <> data{%d} (e=%d)",
              shortFL(FL),fct, i+1, j+1,e);
          }
          else { l=snprintf(istr,n,
             "%s %s() non-matching cgs{%d,%d} <> cgs{%d,%d} (e=%d)",
              shortFL(FL),fct, i+1, k+1, j+1, k+1,e);
          }
          if (l>=n) wblog(FL,
             "ERR %s() string out of bounds (%d/%d)",FCT,l,n);

          if (str[0]) strcat(str,"\n");
          strcat(str,istr);

          if (vflag) {
             const char *s=strstr(istr,"non-match");
             if (s) wblog(F_L,"%s ",s);
          }
          else if (vflag=='!') { MXPut(FL,"qs")
             .add(*this,"A").add(i+1,"i").add(j+1,"j")
             .add(e<10 ? -1 : k+1,"k").add(e,"e");
          }

          return 0;
       }
    }

    return 1;
};
 

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isQSym(const char *F, int L, char dflag) const {

    wbMatrix<TQ> Q1,Q2;
    wbperm P1,P2,iP2,P;

    unsigned i,j,l,k=QDIM,r2; int e=0, r=-1;
    gotCGS(FL);

    if (isEmpty()) return 1;
    str[0]=0; isConsistent(r,FL);

    if (r%2) {
       sprintf(str,"%s() requires even-rank (%d)",FCT,r);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    r2=(unsigned)r/2;

    Q1=QIDX; Q1.sortRecs(P1);

    P.init((unsigned)r).Rotate(r2);
    QIDX.blockPermute(P,Q2); Q2.sortRecs(P2);

    if (Q1!=Q2) { 
       sprintf(str,"QIDX[:,1] does not match QIDX[:,2]");
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }
    if (!Q1.isUniqueSorted()) {
       sprintf(str,"WRN QIDX is not unique!");
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    P2.getIPerm(iP2);

    if (dflag && !strchr("rc",dflag)) wblog(FL,
       "ERR %s() invalid dflag=%c<%d>",FCT,dflag,dflag);
    if (dflag && r!=2) wblog(FL,
       "ERR %s() expecting rank-2 for dflag=%c (%d)",FCT,dflag,r);

    for (i=0; i<DATA.len; ++i) {
       j=P1[iP2[i]]; if (j<i) continue;

       const wbvector<unsigned> &si=DATA[i]->SIZE, &sj=DATA[j]->SIZE;

       if (!dflag) {
          for (l=0; l<r2; ++l) 
          if (si[l]!=sj[l+r2] || si[l+r2]!=sj[l]) { sprintf(str,
             "data skew dimensional (%d,%d: %s vs. %s)",i+1,j+1,
              si.toStrD().data, sj.toStrD().data);
              return 0;
          }
       }
       else
       if (i!=j || (dflag=='r' && si[0]!=1) || (dflag!='r' && si[1]!=1)) {
           sprintf(str, "data expected to be diagonal (%c<%d>) "
           "\n(%d,%d: %s vs. %s)", dflag, dflag, i+1,j+1,
           si.toStrD().data, sj.toStrD().data);
           return 0;
       }

       for (k=0; !e && k<CGR.dim2; ++k) {

          if (!CGR(i,k).cgr || !CGR(j,k).cgr) {
             if (CGR(i,k).cgr || CGR(j,k).cgr || !qtype[k].isAbelian()) {
             wblog(FL,"ERR %s() (%d,%d; %d): 0x%lX 0x%lX (%s) !?",
                FCT,i,j,k, CGR(i,k).cgr, CGR(j,k).cgr,
                qtype[k].toStr().data
             ); }
             continue;
          }

          wbvector<unsigned> sa, sb;
             CGR(i,k).getSize(sa,'b');
             CGR(j,k).getSize(sb,'b');

          if (sa.len!=sb.len || (sa.len && int(sa.len)!=r)) wblog(FL,
             "ERR %s(%d,%d;%d) CGC size mismatch: %s <> %s (%d)\n%s\n%s",
             FCT,i+1,j+1,k+1, sa.toStrD().data, sb.toStrD().data, r,
             CGR(i,k).toStr().data, CGR(j,k).toStr().data
          );
          if (!sa.len) continue;

          for (l=0; l<r2; ++l) {
          if (sa[l]!=sb[l+r2] || sa[l+r2]!=sb[l]) {
             sprintf(str,
               "data skew dimensional (%d,%d (%d): %s vs. %s)",
                i+1,j+1,k+1, sa.toStrD().data, sb.toStrD().data
             );
             if (F) wblog(F,L,"ERR %s",FCT,str); 
             return 0;
          }}
       }
    }

    return 1;
};


template <class TQ, class TD>
template <class T2>
bool QSpace<TQ,TD>::hasSameQ(const QSpace<TQ,T2> &B) const {

    const QSpace<TQ,T2> &A=(*this);

    unsigned i,j,d,k,l,m; int r=-1;
    wbperm pA, pB;
    wbvector<unsigned> dA, dB;
    wbMatrix<TQ> QAk, QBk;
    wbindex Ia, Ib;

    if (isEmpty() || B.isEmpty()) return 1;
    if (QIDX.dim2!=B.QIDX.dim2 || QDIM!=B.QDIM) return 0;

    if (!A.isConsistent(r)) { A.info("A"); wberror(FL,str); }
    if (!B.isConsistent(r)) { B.info("B"); wberror(FL,str); }

    for (k=0; (int)k<r; ++k) {
       A.getQsub(k, QAk); QAk.groupRecs(pA,dA);
       B.getQsub(k, QBk); QBk.groupRecs(pB,dB);

       for (l=i=0; i<dA.len; ++i, l+=d) {
          d=dA[i]; m=dA[i]=DATA[pA[l]]->SIZE[k];

          for (j=1; j<d; ++j)  
          if (m!=(DATA[pA[l+j]]->SIZE[k])) {"AA");
          wblog(FL,"ERR severe QSpace inconsistency (%d,%d) ???",
          m,DATA[pA[l+j]]->SIZE[k]); }
       }

       for (l=i=0; i<dB.len; ++i, l+=d) {
          d=dB[i]; m=dB[i]=B.DATA[pB[l]]->SIZE[k];

          for (j=1; j<d; ++j)  
          if (m!=(B.DATA[pB[l+j]]->SIZE[k])) {"BB");
          wblog(FL,"ERR severe QSpace inconsistency (%d,%d) ???",
          m, B.DATA[pB[l+j]]->SIZE[k]); }
       }

       matchIndex(QAk, QBk, Ia, Ib);
       
       for (i=0; i<Ia.len; ++i)
       if (!QAk.recEqual(Ia[i], QBk.rec(Ib[i]))) return 0;
    }

    return 1;
}


template<class TQ, class TD>
template<class TDA, class TDB>
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCG(
   const wbvector< QSpace<TQ,TDA> > &A, const wbindex &Ia,
   const wbvector< QSpace<TQ,TDB> > &B, const wbindex &Ib,
   char vflag
){
   unsigned i,l,n1,n2, QDIM, dim2=A[0].QIDX.dim2;
   wbMatrix<TQ> Qa,Qb;
   wbvector<INDEX_T> Sa,Sb;
   QVec qtype=A[0].qtype;
   QMap<TQ> M;

   wbvector<char> m1(A.len), m2(B.len);

   if (!A.len) wblog(FL,"ERR got empty input space 1");
   if (!B.len) wblog(FL,"ERR got empty input space 2");

   for (i=0; i<A.len; ++i) if (A[i].isEmpty()) { wblog(FL,
       "WRN got empty space A[%d]",i+1); m1[i]++; }
   for (i=0; i<B.len; ++i) if (B[i].isEmpty()) { wblog(FL,
       "WRN got empty space B[%d]",i+1); m2[i]++; }

   n1=m1.nnz(); n2=m2.nnz();
   if (n1 || n2) {
      wbvector< QSpace<TQ,TDA> > X(A.len);
      wbvector< QSpace<TQ,TDB> > Y(B.len);
      for (l=i=0; i<A.len; ++i) {
         if (!A[i].isEmpty()) X[l++].init2ref(A[i]); }
      if (l) X.len=l; else wblog(FL,
         "ERR got %d empty input space 1",A.len);
      for (l=i=0; i<B.len; ++i) {
         if (!B[i].isEmpty()) Y[l++].init2ref(B[i]); }
      if (l) Y.len=l; else wblog(FL,
         "ERR got %d empty input space 1",B.len);
      wblog(FL,"WRN %s() => skip empty input spaces", FCT);

      return initIdentityCG(X,Ia,Y,Ib,vflag);
   }

   iTag itA, itB;
     itA.init(FL,A,Ia);
     itB.init(FL,B,Ib);

   init(); QDIM=A[0].QDIM;

   for (i=1; i<A.len; ++i) if (A[i].qtype!=qtype) wblog(FL,
      "ERR qtype mismatch of input space 1 (%d: %s; %s)",
       i+1, A[i].qStr().data, qStr().data);
   for (i=0; i<B.len; ++i) if (B[i].qtype!=A[0].qtype) wblog(FL,
      "ERR qtype mismatch of input space 2 (%d: %s; %s)",
       i+1, B[i].qStr().data, qStr().data);
   for (i=1; i<A.len; ++i)
   if (A[i].QDIM!=QDIM || A[i].QIDX.dim2!=dim2) wblog(FL,
      "ERR %s() size mismatch (%d: %d/%d %d/%d)",
       i+1, A[i].QDIM, QDIM, A[i].QIDX.dim2, dim2);
   for (i=1; i<B.len; ++i)
   if (B[i].QDIM!=QDIM || B[i].QIDX.dim2!=B[0].QIDX.dim2) wblog(FL,
      "ERR %s() size mismatch (%d: %d/%d %d/%d)",
       FCT,i+1, B[i].QDIM, QDIM, B[i].QIDX.dim2, B[0].QIDX.dim2);

   getQDimGen(A,Qa,Sa,Ia);
   getQDimGen(B,Qb,Sb,Ib);

   if (qtype.isEmpty()) {
      qtype.init(QDIM);
      for (i=0; i<QDIM; ++i) qtype[i]=QTYPE_U1;
   }

   gCG.getQfinal(FL,qtype,Qa,Qb,M,'c');
   M.getIdentityQ(FL,Sa,Sb, *this, vflag);

   if (itags.len!=3) wblog(FL,"ERR %s() itags.len=%d",FCT,itags.len);
   if (!itags[0].t) itags[0]=itA;
   if (!itags[1].t) itags[1]=itB;


   return *this;
};


template<class TQ, class TD>
template<class TDA>
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCG(
   const wbvector< QSpace<TQ,TDA> > &A, const wbindex &ia,
   char zflag
){
   unsigned i, dim2=A[0].QIDX.dim2;
   wbMatrix<TQ> Qa,Qb;
   wbvector<INDEX_T> Sa,Sb;
   wbMatrix<INDEX_T> Sc;

   if (!A.len) wblog(FL,"ERR got empty input space 1");

   init(); QDIM=A[0].QDIM; qtype=A[0].qtype;

   for (i=1; i<A.len; ++i) if (A[i].qtype!=qtype) wblog(FL,
      "ERR qtype mismatch of input space 1 (%d: %s; %s)",
       i+1, A[i].qStr().data, qStr().data);
   for (i=1; i<A.len; ++i)
   if (A[i].QDIM!=QDIM || A[i].QIDX.dim2!=dim2) wblog(FL,
      "ERR %s() size mismatch (%d: %d/%d %d/%d)",
       i+1, A[i].QDIM, QDIM, A[i].QIDX.dim2, dim2);

   getQDimGen(A,Qa,Sa,ia,&Sc);

   if (qtype.isEmpty()) {
      qtype.init(QDIM);
      for (i=0; i<QDIM; ++i) qtype[i]="A";
   }

   itags.init(2);
   itags[1]=itags[0].init(FL,A,ia);
   itags[1].Conj();


   QIDX.Cat(2,Qa,Qa);
   setupDATA();
   setupCGR();

   unsigned j=0, d[qtype.len]; TQ *qij;
   for (; j<qtype.len; ++j) d[j]=qtype[j].qlen();

   for (i=0; i<DATA.len; ++i) {
      DATA[i]->initIdentity(Sa[i]);

      for (qij=Qa.rec(i), j=0; j<CGR.dim2; ++j) {
         CGR(i,j).initIdentityR(
            FL, qtype[j], qij, Sc.data ? Sc(i,j) : -1);
         qij+=d[j];
      }
   }

   if (!zflag) return *this;

   itags[1].Conj();

   for (i=0; i<DATA.len; ++i) {
      for (qij=QIDX.rec(i), j=0; j<CGR.dim2; ++j) {
         qtype[j].getDual(qij,qij+QDIM);
         CGR(i,j).initIdentityZ(
            FL, qtype[j], qij, Sc.data ? Sc(i,j) : -1);
         qij+=d[j];
      }
   }

   if (zflag<0) Conj();

   return *this;
};


template <class TQ, class TD>
TD QSpace<TQ,TD>::norm2() const {

   if (gotCGS(FL)>0) {
      unsigned i=0, j=0; TD c2, x2=0;
      for (i=0; i<DATA.len; ++i) {
         for (c2=1,j=0; j<CGR.dim2; ++j) {
            if (!CGR(i,j).isEmpty()) c2*=CGR(i,j).norm2();
         }
         x2+=c2*(DATA[i]->norm2());
      }
      return x2;
   }
   else {
      TD x2=0;
      for (unsigned i=0; i<DATA.len; ++i) x2+=(DATA[i]->norm2());
      return x2;
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::checkNorm(
   const char *F, int L, double nrm, char vflag, double eps
 ) const {

   double dbl=norm();
   
   if (fabs(dbl-nrm)>1E-12) wblog(F,L,
      "ERR invalid norm %.4g/%.4g (%.3g)",dbl,nrm,(nrm-dbl)/nrm);
   else if (vflag) wblog(F,L,
      "TST norm = %.4g/%.4g (%.3g)",dbl,nrm,(nrm-dbl)/nrm
   );
};


template <class TQ, class TD>
unsigned QSpace<TQ,TD>::RemoveZLabels(
  const char *F, int L, wbMatrix<TQ> &QQ, wbMatrix<TQ> *Z
) const {

   if (QDIM<=0) wblog(FL,"ERR %s() missing QDIM (=%d)",FCT,QDIM);

   unsigned mz=qtype.Qrank(), m0=qtype.Qlen(), m=m0+mz;
   
   if (&QQ==&QIDX) { if (QDIM!=m) wblog(F_L,
      "ERR no zlabels included in *this (%d <> %d+%d)",QDIM,m0,mz); }
   else if (QDIM!=m0) wblog(F_L,
      "ERR severe QSpace inconsistency (%d/%d (+%d))",QDIM,m0,mz);
   if (QQ.dim2%m) wblog(F_L,
      "ERR incompatible QQ-label space (%dx%d/(%d+%d))",
       QQ.dim1,QQ.dim2,m0,mz);

   if (QQ.isEmpty()) {
      QQ.dim1=0; if (m && QQ.dim2)
      QQ.dim2=m0*(QQ.dim2/m);
      return m0;
   }

   int r=QQ.dim2/m;
   unsigned i,j,n, l=0, i0=0, iz=0, N=QQ.dim1*r;

   wbvector<unsigned> d0_,dz_;
   wbvector<char> isz_(m);
   char *isz=isz_.data;

   TQ *qz=NULL, *qq=QQ.data; 

   qtype.Qlen(d0_,dz_);
   if (m0!=d0_.sum() || mz!=dz_.sum()) wblog(FL,
      "ERR %s() severe qdim-inconsistency: %d/%d, %d/%d",
       FCT, m0, d0_.sum(), mz, dz_.sum());

   if (Z) { Z->init(QQ.dim1,mz*r); qz=(Z ? Z->data : NULL); }
   if (!mz) return m0;

   for (l=i=0; i<qtype.len; ++i)
   for (l+=d0_[i], n=dz_[i], j=0; j<n; ++j) isz[l++]=1;

   for (l=i=0; i<N; ++i) {
      for (j=0; j<m; ++j, ++l) {
         if (isz[j]==0) qq[i0++]=qq[l];
         else { if (qz) qz[iz]=qq[l]; ++iz; }
      }
   }

   if (!i0) wblog(F_L,"ERR removing z-labels leaves nothing !?");
   else if (i0 != QQ.dim1*m0*r || iz != QQ.dim1*mz*r) wblog(F_L,
     "ERR %s() severe size inconsistency\n%dx(%d*%d)=%d =>\n"
     "[ %dx(%d*%d)) = %d/%d ] + [ %dx(%d*%d)) = %d/%d ]",FCT,
       QQ.dim1,m, r, QQ.dim1*QQ.dim2,
       QQ.dim1,m0,r, QQ.dim1*m0*r, i0,
       QQ.dim1,mz,r, QQ.dim1*mz*r, iz
   ); 

   QQ.dim2=m0*r;
   return m0;
};


template <class TQ, class TD>
void QSpace<TQ,TD>::initOpZ_WET(const char *F, int L,
   const QVec &qvec,
   const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q,
   const wbarray<TD> &D3,
   const double eps1, const double eps2
){
   unsigned d,i,j,k,l,id,nd, k1,k2, is,i0, r=0, im, got_omult=0;
   unsigned i1,i2,i3, j1,j2,j3, l1,l2,l3, N1,N2,N3, n1,n2,n3;
   unsigned qdimz, qdimz3, dz;

   wbvector<INDEX_T> d1,d2,d3, I1,I2,I3, M1,M2,M3, s1,s2,s3, D,I,S;
   wbvector<unsigned> m;
   wbvector<TD> d0,dc;
   wbvector<double> cc;
   wbarray<char> mark;
   wbIndex m1;
   double x;
   TD w;

   wbMatrix<TQ> Q1r(Q1),Q2r(Q2),Q3r(Q), Q1u(Q1),Q2u(Q2),Q3u(Q);
   wbMatrix<TQ> QZ,QX,Z,q2,q3,qc; qset<TQ> q1;
   wbMatrix<double> x1,x2;
   const wbperm P3("3,1,2");
   wbperm p1,p2,p3, P;
   wbindex J1,J2,J3;
   QMap<TQ> M;


   clearQSpace(); qtype=qvec; QDIM=qvec.Qlen(); dz=qvec.Qrank();

   qdimz=QDIM+dz;
   qdimz3=3*qdimz;

   if (Q1.dim2!=qdimz || Q2.dim2!=qdimz || Q.dim2!=qdimz) wblog(FL,
      "ERR inconsistency in number of symmetry labels\n"
      "(%d,%d,%d/%d)", Q1.dim2, Q2.dim2, Q.dim2, qdimz);
   if (D3.SIZE.len!=3) wblog(FL,
      "ERR expecting rank-3 array (%s)",D3.sizeStr().data);

   N1=D3.SIZE[0];
   N2=D3.SIZE[1];
   N3=D3.SIZE[2];

   if (Q1.dim1!=N1 || Q2.dim1!=N2 || Q.dim1!=N3) wblog(FL,
      "ERR dimension mismatch: symmetry labels <> data\n"
      "(%dx%dx%d <> %s)",Q1.dim1,Q2.dim1,Q.dim1,D3.sizeStr().data
   );

   RemoveZLabels(F_L,Q1r); Q1r.groupRecs(p1,d1,Q1u,-1,1,&J1);;
   RemoveZLabels(F_L,Q2r); Q2r.groupRecs(p2,d2,Q2u,-1,1,&J2);;
   RemoveZLabels(F_L,Q3r); Q3r.groupRecs(p3,d3,Q3u,-1,1,&J3);;

   gCG.getQfinal_zdim(F,L,qvec,Q2u,Q3u,Q1u,s2,s3,s1,&m);
      n1=Q1u.dim1;
      n2=Q2u.dim1;
      n3=Q3u.dim1;


   for (i=0; i<n1; ++i) {
      if (d1[i]%s1[i]==0) d1[i]/=s1[i];
      else {
         MXPut(FL,"a").add(D3,"D3").add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q")
         .add(qvec,"qvec").add(Q1r,"Q1r").add(Q2r,"Q2r").add(Q3r,"Q3r")
         .add(Q1u,"Q1u").add(Q2u,"Q2u").add(Q3u,"Q3u")
         .add(J1,"J1").add(J2,"J2").add(J3,"J3")
         .add(d1,"d1").add(s1,"s1").add(s2,"s2").add(s3,"s3")
         .add(m,"m").add(i+1,"i");
         wblog(FL,"ERR incomplete WET multiplet space (1)\n"
         "%d/%d(m=%d); iq=%d/%d, D=%d)",d1[i],s1[i],m[i],i+1,n1,Q1.dim1);
      }
   }
   for (i=0; i<n2; ++i) {
      if (d2[i]%s2[i]==0) d2[i]/=s2[i];
      else wblog(FL,
        "ERR incomplete WET multiplet space (2: %d/%d; %d) !??",
         d2[i],s2[i],Q2.dim1
      );
   }
   for (i=0; i<n3; ++i) {
      if (d3[i]%s3[i]==0) d3[i]/=s3[i];
      else wblog(FL,
        "ERR incomplete WET multiplet space (3: %d/%d; %d) !??",
         d3[i],s3[i],Q.dim1
      );
   }

   M1.init(N1); M2.init(N2); M3.init(N3);
   I1.init(n1); I2.init(n2); I3.init(n3);

   for (i=0; i<N1; i+=d) { j=J1[i]; d=s1[j];
      for (l=1; l<d; ++l) {
         if (Q1r.recCompare(i,i+l)) wblog(FL,
          "ERR %s() multiplet not grouped in z-labels (1: %d/%d)",
           FCT,i,Q1.dim1); 
      }
      for (id=I1[j]++, l=0; l<d; ++l) M1[i+l]=id;
   }
   if (I1!=d1) wblog(FL,"ERR %s() missed multiplets !?",FCT);

   for (i=0; i<N2; i+=d) { j=J2[i]; d=s2[j];
      for (l=1; l<d; ++l) {
         if (Q2r.recCompare(i,i+l)) wblog(FL,
          "ERR WET: multiplet not grouped in z-labels (2: %d/%d)",
           i,Q2.dim1); 
      }
      for (id=I2[j]++, l=0; l<d; ++l) M2[i+l]=id;
   }
   if (I2!=d2) wblog(FL,"ERR %s() missed multiplets !?",FCT);

   for (i=0; i<N3; i+=d) { j=J3[i]; d=s3[j];
      for (l=1; l<d; ++l) {
         if (Q3r.recCompare(i,i+l)) wblog(FL,
          "ERR WET: multiplet not grouped in z-labels (3: %d/%d)",
           i,Q.dim1); 
      }
      for (id=I3[j]++, l=0; l<d; ++l) M3[i+l]=id;
   }
   if (I3!=d3) wblog(FL,"ERR %s() missed multiplets !?",FCT);

   n3=D3.numel(); mark.init(D3.SIZE);
   for (nd=i=0; i<n3; ++i) { x=fabs(D3.data[i]);
      if (x>eps1) { ++nd; mark.data[i]=1; }
      else if (x>eps2) wblog(F,L,
      " *  ignoring small value data[%d]=%.3g @ %g",i+1,D3.data[i],eps2);
   }
   if (!nd) { 
      wblog(FL,"WRN %s() got all-zero matrix elements !??",FCT);
      init(); return;
   }

   QZ.init(nd,qdimz3); d0.init(nd);
   QIDX.init(nd,3*QDIM);
   setupDATA(FL); setupCGR(FL);

   for (i3=0; i3<N3; i3+=n3) { j3=J3[i3]; n3=s3[j3];
   for (i2=0; i2<N2; i2+=n2) { j2=J2[i2]; n2=s2[j2];
   for (i1=0; i1<N1; i1+=n1) { j1=J1[i1]; n1=s1[j1];
      d0.len=QZ.dim1=nd; l=0;

      for (k=0; k<n3; ++k) { l3=i3+k;
      for (j=0; j<n2; ++j) { l2=i2+j;
      for (i=0; i<n1; ++i) { l1=i1+i; if (!mark(l1,l2,l3)) continue; 
         QZ.recSetB(l,0,qdimz,Q1.rec(l1));
         QZ.recSetB(l,1,qdimz,Q2.rec(l2));
         QZ.recSetB(l,2,qdimz,Q .rec(l3)); d0[l++]=D3(l1,l2,l3);
      }}}
      if (!l) continue;

      d0.len=QZ.dim1=l;
      QX=QZ; QX.SkipTiny_float();
      RemoveZLabels(F_L,QX,&Z);

      q1.init(  QDIM, QX.ref(0,     0));
      q2.init(1,QDIM, QX.ref(0,  QDIM));
      q3.init(1,QDIM, QX.ref(0,2*QDIM));

      gCG.getQfinal(F,L,qvec,q2,q3,M,'c');

      m1.init();
      M.getCGZlist(F,L,qc,dc,&I,NULL,&q1,&m1);

      if (!I.allEqual(I[0])) wblog(FL,"ERR %s() CData not unique",FCT);

      qc.BlockPermute(P3);  cc.init(1); cc[0]=1;
      QZ.isref=d0.isref=1; {
         QZ.SkipTiny_float(); QZ.sortRecs(P); d0.Permute(P);
         qc.SkipTiny_float(); qc.sortRecs(P); dc.Permute(P);
      }; QZ.isref=d0.isref=0;

      if (m1.numel()>1) {
         unsigned l, mq=m1.numel();
         wbvector< wbMatrix<TQ> > QC(mq+1);
         wbvector< wbvector<TD> > DC(mq+1);
         wbindex Ia,Ib;

         wblog(FL," *  %s() got outer multiplicity [%s]",
            FCT, m1.SIZE.toStr().data); 

         for (l=0; l<=mq; ++l) {
            if (l) { 
               if (!M.getCGZlist(F,L,qc,dc,&I,NULL,&q1,&m1)) break;
               qc.BlockPermute(P3);
               qc.SkipTiny_float(); qc.sortRecs(P); dc.Permute(P);
            }
            qc.findUnique(Ia);
            qc.getRecs(Ia,QC[l]);
            dc.select(Ia,DC[l]);
         }
         if (l!=mq) wblog(FL,"ERR %s() !?? (%d/%d)",FCT,l,mq); 
         
         QZ.findUnique(Ia);
         QZ.getRecs(Ia,QC[l]);
         d0.select(Ia,DC[l]); ++l;

         qc=QC[0];
         for (l=1; l<QC.len; ++l) {
            matchSortedIdx(qc,QC[l],Ia,Ib);
            QC[l]=QC[l].getRecs(Ib,qc); DC[l].Select(Ib);
            if (qc.isEmpty()) break;
         }
         if (qc.dim1<QC.len) wblog(FL,
            "ERR %s() insufficient number of unique recs (%d/%d)",
            FCT,qc.dim1,QC.len
         ); l=DC.len-2;
         if (DC[l].len!=DC[l+1].len) wblog(FL,
            "ERR %s() mismatch of unique recs (%d/%d)",
            FCT,DC[l].len,DC[l+1].len
         );
         for (l=0; l<QC.len; ++l) {
            matchSortedIdx(qc,QC[l],Ia,Ib);
            if (Ia.len!=qc.dim1) wblog(FL,
               "ERR %s() %d/%d",FCT,Ia.len,qc.dim1);
            QC[l].Set2Recs(Ib); DC[l].Select(Ib);
         }

         wbarray<double> b,c,X2,Xc,iX,X(qc.dim1, mq);
         for (l=0; l<mq; ++l) { X.setCol(l,DC[l].data); }

         wbMatProd(X,X,X2,'C');
         wbMatVProd(X,DC[l],b,'C');
         wbVMatProd(wbInverse(FL,X2,iX),b,c,'C');

         cc.init(c.numel(), c.data); m1.reset();
         wbVMatProd(X,c,Xc); x=Xc.normDiff(DC[l]);
         if (x>1E-12) wblog(FL,
            "ERR %s() failed to match outer mult. (%.3g)",FCT,x);
         else wblog(FL,
            " *  %s() matched outer multiplicity (@ %.2g)\n==> "
            "[%s] (%.5g)",FCT,x,cc.toStr().data,cc.norm2()); 

         if (!M.getCGZlist(F,L,qc,dc,&I,NULL,&q1,&m1,&cc)) wblog(FL,
            "ERR %s() failed to get CGZ list (%s)",FCT,m1.toStr().data);

         qc.BlockPermute(P3);
         qc.SkipTiny_float(); qc.sortRecs(P); dc.Permute(P);

      }

      if (d0.len!=dc.len) {

         MXPut X(FL,"a"); wbvector<double> i(3);
         X.add(D3,"D3").add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q")
          .add(Q1r,"Q1r").add(Q2r,"Q2r").add(Q3r,"Q3r")
          .add(Q1u,"Q1u").add(Q2u,"Q2u").add(Q3u,"Q3u")
          .add(QZ,"QZ").add(d0,"d0").add(qc,"qc").add(dc,"dc")
          .add(eps1,"eps1").add(eps2,"eps2")
          .add(q1,"q1").add(q2,"q2").add(q3,"q3").add(m,"m").add(M,"M");
         i[0]=i1; i[1]=i2; i[2]=i3; X.add(i+1,"i");
         i[0]=j1; i[1]=j2; i[2]=j3; X.add(i+1,"j");
         i[0]=n1; i[1]=n2; i[2]=n3; X.add(i,"n");
         i[0]=N1; i[1]=N2; i[2]=N3; X.add(i,"N");
         wblog(FL,"ERR WET data mismatch (len=%d/%d)",d0.len,dc.len);
      }

      x=SQRT(double(QZ.normDiff2(qc))); if (x>eps1) {
        MXPut(FL,"q").add(dc,"dc").add(qc,"qc").add(QZ,"QZ").add(d0,"d0"); 
        wblog(FL,"ERR %s() data mismatch (%.3g)",FCT,x);
      }

      k1=QZ.findUniqueRecSorted1(); if (int(k1)<0) {
         MXPut(FL,"a").add(QZ,"QZ").add(k1,"k1");
         wblog(FL,"ERR %s() failed to find unique Q3-record",FCT);
      }
      k2=qc.findRecSorted(QZ.rec(k1)); if (int(k2)<0)
         wblog(FL,"ERR %s() failed to find matching Q3-record",FCT);
      w=d0[k1]/dc[k2];

      x1.init(FL,QZ,d0).sortRecs();
      x2.init(FL,qc,dc);
      x2.scaleCol(x2.dim2-1,w).sortRecs();

      x=x1.normDiff(x2); if (x>eps1) {
         wbvector<double> i(3);
         MXPut X(FL,"a"); X.add(D3,"D3")
          .add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q").add(QZ,"QZ").add(M,"M")
          .add(d0,"d0").add(qc,"qc").add(dc,"dc").add(w,"w");
         i[0]=i1; i[1]=i2; i[2]=i3; X.add(i+1,"i");
         i[0]=n1; i[1]=n2; i[2]=n3; X.add(i,"n");
         i[0]=N1; i[1]=N2; i[2]=N3; X.add(i,"N");
         i[0]=j1; i[1]=j2; i[2]=j3; X.add(i+1,"j");
         X.add(x1,"x1").add(x2,"x2");
        wblog(FL,"ERR WET data mismatch (%.3g)",x);
      }
      x=fabs(w); if (x>1E3 || x<1E-3)
        wblog(FL,"WRN %s() got WET=%.3g !??",FCT,w);

      i0=I[0]; I.init(3); m1.reset();
      for (im=0; im<cc.len; ++im, ++r) {
         QIDX.recSetB(r,0,QDIM,q1.data);  I[0]=M1[i1];
         QIDX.recSetB(r,1,QDIM,q2.data);  I[1]=M2[i2];
         QIDX.recSetB(r,2,QDIM,q3.data);  I[2]=M3[i3];

         if (!(++m1)) wblog(FL,
            "ERR %s() counter on outer multiplicity !??",FCT);

         for (is=0; is<qvec.len; ++is) {
            if (M.cg3(i0,is)->len!=cc.len) { wblog(FL,
               "ERR %s() QMap mismatch (%d/%d)",
               FCT,M.cg3(i0,is)->len, cc.len);
            }

            CGR(r,is)=M.cg3(i0,is)->el(im);
            CGR(r,is).Permute(P3);
         }

         wbarray<TD> &a = DATA[r]->init(d1[j1], d2[j2], d3[j3]);
         a.element(I)=w*cc[im];

         if (fabs(w*cc[im])<1E-12) wblog(FL,"WRN %s() "
            "got w = %g * %g = %g (im=%d)",FCT,w,cc[im],w*cc[im],im);

         if (im) ++got_omult;
      }
   }}}

   if (!r || r>QIDX.dim1) {
      MXPut(FL,"a").add(s1,"s1").add(s2,"s2").add(s3,"s3")
      .add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q")
      .add(Q1u,"Q1u").add(Q2u,"Q2u").add(Q3u,"Qu").add(mark,"mark");
      wblog(FL,"ERR %s (r=%d)",FCT,r);
   }

   QIDX.dim1 = DATA.len = CGR.dim1 = r;

   otype=QS_OPERATOR;

   itags.init_qdir("-++");

   Conj();

   makeUnique(got_omult ? 'v':0);

};



template <class TQ, class TD>
void QSpace<TQ,TD>::reduceMatEl(
   const char *F, int L, 
   const wbMatrix<TQ> &QQ, const wbvector<TD> &dd,
   const wbMatrix<TQ> &Qc, const wbvector<TD> &dc,
   wbvector<unsigned> &J,
   wbvector<unsigned> &Jc,
   wbvector<TD> &dr
) const {

   unsigned i,j,k,m,n,r; char e=0;
   wbvector<unsigned> d1,d2,D1,D2;
   wbperm pp,pc;
   wbindex I1,I2,i1,i2;
   wbvector<TD> xx,x1,x2;
   double x=0;

   wbMatrix<TQ> Z1,Z2,q2;
   QSpace<TQ,TD> A1,A2;

   if (QQ.dim1!=dd.len || Qc.dim1!=dc.len || QQ.dim2!=Qc.dim2) wblog(F,L,
      "ERR dimension mismatch (%d/%dx%d; %d/%dx%d)",
       dd.len, QQ.dim1, QQ.dim2, dc.len, Qc.dim1, Qc.dim2
   );


   A1.QIDX=QQ; A1.qtype=qtype;
   A1.RemoveZLabels(F,L,&Z1); A1.QIDX.groupRecs(pp,d1,Z1);

   A2.QIDX=Qc; A2.qtype=qtype;
   A2.RemoveZLabels(F,L,&Z2); A2.QIDX.groupRecs(pc,d2,Z2);

   i=matchIndex(A1.QIDX,A2.QIDX,I1,I2);
   if (i) wblog(F,L,"ERR got non-unique Q-set (%d)",i);
   if (I1.len!=A1.QIDX.dim1) {
       wbindex Ix; I1.invert(A1.QIDX.dim1,Ix);
       for (i=0; i<Ix.len; ++i) A1.QIDX.recPrint(Ix[i],"Q");
       wblog(F,L,
         "ERR got %d additional symmetry sectors (%d/%d matching)",
          A1.QIDX.dim1-I1.len, I1.len, A1.QIDX.dim1
       ); 
   }

   pc.BlockSelect(I2,d2); Qc.getRecs(pc,q2);

   matchIndex(QQ,q2,i1,i2);

   if (i1.len!=QQ.dim1) { ++e;
       wbindex Ix; i1.invert(QQ.dim1,Ix);
       for (i=0; i<Ix.len; ++i) QQ.recPrint(Ix[i],"QZ");
       wblog(F,L,"==> ERR got additional z-labels (%d/%d)",
       QQ.dim1-i1.len,QQ.dim1); 
   }

   if (i2.len!=q2.dim1) { ++e;
       wbindex Ix; i1.invert(q2.dim1,Ix);
       for (unsigned i=0; i<Ix.len; ++i) q2.recPrint(Ix[i],"Qall");
       wblog(F,L,
          "==> ERR got missing z-labels (%d/%d) !??",
           q2.dim1-i2.len,q2.dim1
       ); 
   }

   if (e) ExitMsg("\n");"ERR"); does


   if (dd.len<pp.len || dd.len%pp.len || pc.max()>=dc.len) wblog(F,L,
      "ERR index out of bounds (%d/%d, %d/%d)",
       dd.len, pp.len, pc.max(), dc.len
   );

   m=dd.len/pp.len;

   d1.cumsum0(D1); d2.cumsum0(D2); dr.init(m*I1.len);
   J .init(I1.len);
   Jc.init(I1.len);


   for (i=0; i<I1.len; ++i) { n=d1[I1[i]];
      if (n!=d2[I2[i]]) wblog(F,L,
      "ERR size mismatch (%d/%d)", d1[I1[i]], d2[I2[i]]);

      x1.init(n); x2.init(n); xx.init(n);
      k=D2[i]; Jc[i]=pc[k]; for (j=0; j<n; ++j) x2[j]=dc[pc[k+j]];
      k=D1[i]; J[i]=pp[k];

      for (r=0; r<m; ++r) {
         for (j=0; j<n; ++j){  x1[j]=dd[pp[k+j]+r*pp.len];

             if (fabs(x1[j])<1E-12 || fabs(x2[j])<1E-12) wblog(F,L,
                "ERR got matrix element (%d,%d): %.4g, %.4g",
                 i+1,j+1,x1[j],x2[j]
             );

             xx[j]=x1[j]/x2[j];
             if (j) { if (fabs(x-xx[j])>1E-12) {
#ifdef MATLAB_MEX_FILE
                 MXPut(FL).add(QQ,"QQ").add(dd,"dd").add(Qc,"qc")
                 .add(dc,"dc").add(pp+1,"pp").add(pc+1,"pc")
                 .add(D1+1,"d1").add(D2+1,"d2");
#endif
                 wblog(F,L,
                   "ERR inconsistent reduced matrix element (%d,%d/%d; %d,%d/%d):\n"
                   "%.4g/%.4g = %.4g (%.4g) @ %.4g", 
                    pp[D1[i]]+r*pp.len+1, pp[D1[i]+j]+r*pp.len+1, dd.len,
                    pc[D2[i]]+1, pc[D2[i]+j]+1, dc.len,
                    x1[j], x2[j], xx[j], x, xx[j]-x
                 );
             }}
             else x=xx[j];
         }
         dr[i+r*I1.len]=xx.avg();
      }
   }
};


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Select(const wbindex& I, char data_only) {

   INDEX_T *ix; wbindex Ix;
   char gotcgr=(!CGR.isEmpty());

   if (data_only) {
      if (QIDX.dim1!=I.len)
         wblog(FL,"ERR %s() dimension mismatch (%d/%d; %d/%d)",
         FCT, QIDX.dim1, I.len, CGR.dim1, DATA.len
      ); 
   }

   if (gotcgr && (CGR.dim1!=DATA.len || CGR.dim2!=qtype.len)) wblog(FL,
      "ERR %s() CGR size mismatch (%dx%d ; %dx%d)",
      FCT, CGR.dim1, CGR.dim2, DATA.len, qtype.len
   );

   if (isref) wblog(FL,"ERR %s() got QSpace isref=%d",FCT,isref);
   if (!I.isUnique(DATA.len)) wblog(FL,
      "ERR %s() got invalid index (%d/%d)",FCT,I.max()+1,DATA.len); 

   I.invert(DATA.len,Ix); ix=Ix.data;
   if (!data_only) QIDX.Set2Recs(I);

   if (DATA.len) {
      wbvector< wbarray<TD>* > DD(DATA); DD.select(I,DATA);
      for (unsigned i=0; i<Ix.len; ++i) {
         wbarray<TD>* &d=DD[ix[i]]; WB_DELETE_1(d);
      }
   }
   
   if (gotcgr) {
      CGR.Set2Recs(I);
   }

   return *this;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::Append(
   const wbvector<TQ> &Q,
   const wbarray<TD> &D
){

   wbarray<TD> *d=0; WB_NEW_1(d,wbarray<TD>); d->init(D);

   QIDX.appendRow(Q);
   DATA.Append(d);
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::Append(unsigned k) {
   unsigned i, n=DATA.len;

   QIDX.Resize(QIDX.dim1+k, QIDX.dim2);
   DATA.Resize(n+k);

   for (i=0; i<k; ++i) {
   WB_NEW_1(DATA[n+i],wbarray<TD>); }
}


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Append(
   const char *F, int L, const QSpace<TQ,TD> &B, char uflag
){
   if (B.isEmpty()) return *this;
   if (this==&B) wblog(F_L,"ERR %s() requires distinct B",FCT);
   if (isEmpty()) { *this = B; return *this; }

   int r=-1;

   if (!  isConsistent(r)) wberror(F,L,str);
   if (!B.isConsistent(r)) wberror(F,L,str);

   if (QDIM!=B.QDIM || QIDX.dim2!=B.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
       FCT,QDIM, B.QDIM, QIDX.dim2, B.QIDX.dim2
   );
   if (B.QIDX.dim1!=B.DATA.len || QIDX.dim1!=DATA.len) wblog(FL,
      "ERR %s() severe QIDX/DATA size inconsistency\n%d/%d; %d/%d",
       FCT,B.QIDX.dim1,B.DATA.len,QIDX.dim1,DATA.len
   );
   if (qtype!=B.qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
       FCT,B.qStr().data,qStr().data
   );
   if (itags!=B.itags) wblog(F_L,
      "ERR %s() itag inconsistency ('%s' <> '%s')",
      FCT, itags.toStr().data, B.itags.toStr().data
   );
   if (!CGR.isEmpty() || !B.CGR.isEmpty()) {
      if (qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
          qStr().data,B.CGR.dim1,B.CGR.dim2,CGR.dim1,CGR.dim2
      );
      if (B.CGR.dim1!=B.QIDX.dim1 || CGR.dim1!=QIDX.dim1 ||
          B.CGR.dim2!=B.qtype.len || CGR.dim2!=B.CGR.dim2) wblog(FL,
          "ERR %s() severe QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
          FCT, B.CGR.dim1, B.CGR.dim2, CGR.dim1, CGR.dim2, qtype.len
      ); 
   }

   QIDX.appendRows(B.QIDX.dim1, B.QIDX.data);

   unsigned i=DATA.len, n;
   wbarray<TD> *a=NULL;

   DATA.Append(B.DATA.len, B.DATA.data);
   for (n=DATA.len; i<n; ++i) {
       WB_NEW_1(a,wbarray<TD>); *a = (*DATA[i]);
       DATA[i]=a;
   }

   if (!CGR.isEmpty())
   CGR.appendRows(B.CGR.dim1,B.CGR.data);

   if (uflag) makeUnique(uflag);

   return *this;
};


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Cat(
   const char *F, int L,
   const QSpace<TQ,TD> &A, const QSpace<TQ,TD> &B,
   TD afac, TD bfac, char uflag
){

   if (A.isEmpty()) {
      if (this!=&B   ) { (*this)=B;     }
      if (bfac!=TD(1)) { (*this)*=bfac; }; return *this;
   }
   if (B.isEmpty()) {
      if (this!=&A   ) { (*this)=A;     }
      if (afac!=TD(1)) { (*this)*=afac; }; return *this;
   }
   if (&A==this || &B==this) {
      QSpace<TQ,TD> X; X.Cat(F,L,A,B,afac,bfac,uflag);
      return X.save2(*this);
   }

   int r=-1;
   if (!A.isConsistent(r)) wberror(F,L,str);
   if (!B.isConsistent(r)) wberror(F,L,str);

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
      FCT, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2
   );
   if (A.QIDX.dim1!=A.DATA.len || B.QIDX.dim1!=B.DATA.len) wblog(FL,
      "ERR %s() severe QIDX/DATA size inconsistency\n%d/%d; %d/%d",
      FCT,A.QIDX.dim1,A.DATA.len,B.QIDX.dim1,B.DATA.len
   );
   if (A.qtype!=B.qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
      FCT,A.qStr().data,B.qStr().data
   );
   if (A.itags!=B.itags) wblog(F_L,
      "ERR %s() itag inconsistency ('%s' <> '%s')",
      FCT, A.itags.toStr().data, B.itags.toStr().data
   );
   if (!A.CGR.isEmpty() || !B.CGR.isEmpty()) {
      if (A.qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
         A.qStr().data,A.CGR.dim1,A.CGR.dim2,B.CGR.dim1,B.CGR.dim2
      );
      if (B.CGR.dim1!=B.QIDX.dim1 || A.CGR.dim1!=A.QIDX.dim1 ||
         B.CGR.dim2!=B.qtype.len || A.CGR.dim2!=B.CGR.dim2) wblog(FL,
         "ERR %s() severe QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
         FCT,A.CGR.dim1, A.CGR.dim2, A.CGR.dim1, A.CGR.dim2, A.qtype.len
      ); 
   }

   initQ(A);

   QIDX.Cat(1,A.QIDX,B.QIDX);
   CGR .Cat(1,A.CGR, B.CGR );

   unsigned i,n, na=A.DATA.len;
   wbarray<TD> *a;
   TD one=1;

   DATA.Cat(A.DATA, B.DATA);
   for (n=DATA.len, i=0; i<n; ++i) {
       WB_NEW_1(a,wbarray<TD>); *a = (*DATA[i]);
       DATA[i]=a;
   }

   if (afac!=one) { for (i=0; i<na; ++i) (*DATA[i])*=afac; }
   if (bfac!=one) { for (i=na; i<n; ++i) (*DATA[i])*=bfac; }

   if (uflag) makeUnique(uflag);

   return *this;
};



template <class TQ, class TD>
int QSpace<TQ,TD>::Append2AndDestroy(
   const char *F, int L, QSpace<TQ,TD> &A, char unique
){
   if (isEmpty()) return 0;
   if (this==&A) wblog(F,L,"ERR %s() requires distinct A",FCT);

   int r=-1,e=0;
 
   if (A.isEmpty()) { save2(A); return e; }

   if (!  isConsistent(r)) wberror(F,L,str);
   if (!A.isConsistent(r)) wberror(F,L,str);

   if (QDIM!=A.QDIM || QIDX.dim2!=A.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
       FCT,QDIM, A.QDIM, QIDX.dim2, A.QIDX.dim2
   );
   if (A.QIDX.dim1!=A.DATA.len || QIDX.dim1!=DATA.len) wblog(FL,
      "ERR %s() severe QIDX/DATA size inconsistency\n%d/%d; %d/%d",
       FCT,A.QIDX.dim1,A.DATA.len,QIDX.dim1,DATA.len
   );
   if (A.qtype!=qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
       FCT,A.qStr().data,qStr().data
   );
   if (!CGR.isEmpty() || !A.CGR.isEmpty()) {
      if (qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
          qStr().data,A.CGR.dim1,A.CGR.dim2,CGR.dim1,CGR.dim2
      );
      if (A.CGR.dim1!=A.QIDX.dim1 || CGR.dim1!=QIDX.dim1 ||
          A.CGR.dim2!=A.qtype.len || CGR.dim2!=A.CGR.dim2) wblog(FL,
          "ERR %s() severe QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
          FCT, A.CGR.dim1, A.CGR.dim2, CGR.dim1, CGR.dim2, qtype.len
      ); 
   }

   if (unique && A.hasQOverlap(*this)) { ++e;
      sprintf(str,"%s ERR objects have QIDX overlap", shortFL(F,L));
   }

   A.QIDX.appendRows(QIDX.dim1, QIDX.data);

   A.DATA.Append(DATA.len, DATA.data); if (!CGR.isEmpty()) {
   A.CGR.appendRows(CGR.dim1,CGR.data); }

   QIDX.init(0,QIDX.dim2); DATA.init(); CGR.init(); init();

   return e;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::SkipZeroData(const double eps, char bflag) {

   if (!QIDX.dim2 && QIDX.dim1==1 && DATA.len==1 && CGR.dim1==1) {
      double cfac=1;
      for (unsigned j=0; j<CGR.dim2; ++j) {
         const CRef<TQ> &cj=CGR(0,j);
         if (!cj.isScalar()) wblog(FL,
            "ERR %s() got invalid contraction to scalar\n%s",
            FCT,cj.toStr().data);
         if (cj.cgw.len==1) {
            cfac*=double(cj.cgw[0]); }
         else if (cj.rtype!=CGR_ABELIAN || cj.cgw.len) wblog(FL,
            "ERR %s() invalid scalar or empty\n%s",FCT,cj.toStr().data
         );
      }
      if (cfac<1E-12) wblog(FL,"WRN %s() got cfac=%g !?",FCT,cfac);
      (*DATA[0])*=cfac; CGR.init(1,0);
   }

   unsigned i,j, nz=0;
   bool isz, cgflag=(gotCGS(FL)>0);

   if (isref) wblog(FL,"ERR must not change QSpace reference!");
   if (cgflag && (CGR.dim1!=DATA.len || CGR.dim2!=qtype.len)) wblog(FL,
      "ERR %s() CGR size mismatch (%dx%d ; %dx%d)",
       FCT, CGR.dim1, CGR.dim2, DATA.len, qtype.len
   );

   for (i=0; i<DATA.len; ++i) { wbarray<TD> &a=(*DATA[i]);

      if (cgflag) {
         double cfac=1;
         for (j=0; j<CGR.dim2; ++j) { cfac *= CGR(i,j).norm2(); }

         cfac=sqrt(cfac);

         if (cfac<eps) { isz=1;
            double x=a.aMax(); if (x>1E-6) { wblog(FL,
              "WRN %s() skipping |cgc|=%.3g having max(|data|)=%.3g !?",
              FCT,cfac,x
            ); }
         }
         else {
            if (cfac<1) cfac=1;
            isz=(a.isZero(eps/cfac,bflag));
         }
      }
      else { isz=(a.isZero(eps,bflag)); }

      if (isz) { a.init(); ++nz; }
   } 

   if (nz)
        return SkipEmptyData(FL);
   else return nz;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::SkipEmptyData(const char *F, int L) {

   unsigned j, i=0, l=0, cgflag=(CGR.data ? 1:0);

   if (QIDX.dim1!=DATA.len || (cgflag && CGR.dim1!=DATA.len))
      wblog(FL,"ERR %s() unexpected QSpace (%d/%d/%d)",
      FCT,QIDX.dim1,DATA.len,CGR.dim1);
   if (isref) wblog(F_L,"ERR got QSpace reference");

   for (; i<DATA.len; ++i) {
      if (DATA[i]->isEmpty()) {
         WB_DELETE_1(DATA[i]); if (cgflag) {
         for (j=0; j<CGR.dim2; ++j) { CGR(i,j).init(); }}
      }
      else {
         if (l<i) {
            QIDX.recSet(l,i);
            DATA[l]=DATA[i]; if (cgflag) {
            for (j=0; j<CGR.dim2; ++j) CGR(i,j).save2(CGR(l,j)); }
         }; ++l;
      }
   }

   if (!l) {
      DATA.init();
      QIDX.init(); if (cgflag) CGR.init();
      clearQSpace();
   }
   else if ((i=(DATA.len-l))) {
      QIDX.dim1=l; DATA.len=l; if (cgflag) CGR.dim1=l;
   }

   return i;
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::gotZeroData(unsigned i, const double eps) const {

   double x2, n2=1;

   if (i>=DATA.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,i,DATA.len);

   if (CGR.data) {
      if (i>=CGR.dim1) wblog(FL,
         "ERR %s() index out of bounds (%d/%dx%d)",FCT,i,CGR.dim1,CGR.dim2);
      for (unsigned j=0; j<CGR.dim2; ++j) {
         const CRef<TQ>& c=CGR(i,j);
         if (!c.isRefInit()) {
            n2*=(x2=double(c.norm2()));
            if (x2<1E-12) wblog(FL,"WRN %s() got small CGC (%g)",FCT,x2);
            if (!n2) return 1;
         }
      }
   }
   if (!DATA[i]) {
      wblog(FL,"WRN %s() got NULL data[%d]",FCT,i);
      return 1;
   }
   n2 *= DATA[i]->norm2();

   return (std::sqrt(n2)<=eps);
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::skipZeroOffDiag(double eps, char bflag){

   unsigned i,k,s; int r=-1;
   wbindex I(DATA.len);
   TQ *qi;

   if (isEmpty()) return 0;
   isConsistent(r,FL);

   if (r%2)
   wblog(FL,"ERR %s() requires even-rank (%d)",FCT,r);

   s=QDIM*sizeof(TQ);

   for (k=i=0; i<DATA.len; ++i) { qi=QIDX.rec(i);
      if (memcmp(qi,qi+QDIM,s) && DATA[i]->isZero(eps,bflag)) {
         WB_DELETE_1(DATA[i]);
      }
      else { 
         if (k<i) DATA[k]=DATA[i];
         I[k++]=i;
      }
   }

   i=DATA.len-k;

   if (k==DATA.len) return 0;
   if (k==0) {
      DATA.init();
      clearQSpace(); return i;
   }

   I.len=k; DATA.len=k;

   QIDX.Set2Recs(I);

   return i;
}


template <class TQ, class TD>
void QSpace<TQ,TD>::TransferSpace(const char *F, int L,
   QSpace<TQ,TD> &B, unsigned dim,
   const QSpace<TQ,double> &W, double eps
){
   unsigned i=0; int e=0;
   char gotpartial=0;

   QSpace<TQ,TD> Aout, Bout, &A=*this;
   if ((void*)&A==(void*)&B) wblog(FL,"ERR got same space for A and B");

   if (this->isEmpty()) return;
   if (B.isEmpty()) { B.QDIM=A.QDIM; B.QIDX.init(0,A.QIDX.dim2); }

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) {
      wblog(F_L,"ERR dimension mismatch (%d,%d; %d,%d)",
      A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
   }

   if (!W.isEmpty()) {
      for (i=0; i<W.DATA.len && gotpartial<3; ++i) {
         const wbarray<double> &w=(*W.DATA[i]);
         if (!w.isVector()) wblog(FL, 
            "ERR %s() invalid weight space (%s)",FCT,w.sizeStr().data);
         for (unsigned n=w.numel(), j=0; j<n && gotpartial<3; ++j) {
            if (!(gotpartial & 1) && w.data[j]<eps) gotpartial|=1;
            if (!(gotpartial & 2) && w.data[j]>=eps) gotpartial|=2;
         }
      }
   }
   else gotpartial=1;

   if (gotpartial==2) {
      return;
   }
   if (gotpartial==1) {
      TransferSpace(F_L,B,dim); return;
   }

   wbindex Ia,Ib,Ja,Jw,I2; --dim;

   e=matchIndex(A.QIDX,B.QIDX,Ia,Ib);
   if (e) wblog(FL,"ERR %s() got non-unique QIDX !??",FCT); 
   Ia.invert(QIDX.dim1,I2);

   wbMatrix<TQ> QA,QW;
   A.getQsub(dim,QA);
   W.getQsub(1,QW);

   try {
      matchIndexU(FL,QA,QW,Ja,'f'); }
   catch (...) {
      wblog(F_L,"ERR %s() failed to match input weights uniquely",FCT);
   }

   if (QA.dim1!=Ja.len) wblog(FL,
      "ERR %s() failed to fully match QIDX with W",FCT);

   try {
      for (i=0; i<Ia.len; ++i) {
         DATA[Ia[i]]->Append2(FL,
          *B.DATA[Ib[i]], dim, *W.DATA[Ja[Ia[i]]], eps);
      }

      QSpace<TQ,TD> AX; getSubInit(I2,AX);
      for (i=0; i<I2.len; ++i) {
         DATA[I2[i]]->Append2(FL,
          *AX.DATA[i], dim, *W.DATA[Ja[I2[i]]], eps);

      }
      AX.Append2AndDestroy(FL,B);
   }
   catch (...) {
   #pragma omp critical
      MXPut(FL,"i").add(A,"A").add(B,"B").add(dim,"dim").add(W,"W")
        .add(QA,"QA").add(QW,"QW").add(Ia,"Ia").add(Ib,"Ib").add(Ja,"Ja");
      if (Ia.len && Ja.len)
         wblog(F_L,"ERR %s() %d/%d: %d,%d,%d",
         FCT, i+1,Ia.len, Ia[i]+1, Ib[i]+1, Ja[Ia[i]]+1);
      else wblog(F_L,"ERR %s() %d [%d %d]",FCT,i,Ia.len,Ja.len);
   }

   A.SkipEmptyData(F_L);
   B.SkipEmptyData(F_L);
};


template <class TQ, class TD>
void QSpace<TQ,TD>::TransferSpace(const char *F, int L,
   QSpace<TQ,TD> &B, unsigned dim
){
   unsigned i=0; int e=0, r=rank(FL);

   QSpace<TQ,TD> &A=*this;
   if ((void*)&A==(void*)&B) wblog(FL,"ERR got same space for A and B");

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) {
      wblog(F_L,"ERR dimension mismatch (%d,%d; %d,%d)",
      A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
   }
   if (dim<1 || (int)dim>r) wblog(FL,
      "ERR %s() dimension out of bounds",FCT,dim,r); 

   wbindex Ia,Ib,I2; --dim;

   e=matchIndex(A.QIDX,B.QIDX,Ia,Ib);
   if (e) wblog(FL,"ERR %s() got non-unique QIDX !??",FCT); 
  
   try {
      for (i=0; i<Ia.len; ++i) {
      DATA[Ia[i]]->Append2(*B.DATA[Ib[i]],dim); }
   }
   catch (...) { wblog(F_L,"ERR %s() %d/%d",FCT,i+1,Ia.len); }

   Ia.invert(QIDX.dim1,I2);
   if (I2.len) {
      QSpace<TQ,TD> AX;
      saveSub2(AX,I2); AX.Append2AndDestroy(F_L,B);
   }

   A.init();
};


template <class TQ, class TD>
void QSpace<TQ,TD>::TimesEl(const QSpace<TQ,TD> &B, const char qmatch) {

   unsigned i,j,k; int i1,i2;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   wbstring mark(DATA.len);

   bool cgflag=(!CGR.isEmpty());

   if (isEmpty()) return;
   if (B.isEmpty()) { clearQSpace(); return; }
   checkQ(FL,B);

   if (getQOverlap(B,QQ,IQ,'!'))
   wberror(FL,str);

   for (i=0; i<IQ.dim1; ++i) {
      i1=IQ(i,0); i2=IQ(i,1);

      if (i1<0 || i2<0) {
         if (qmatch) wblog(FL,"ERR QSpaces must match (%d,%d).",i1,i2);
         continue;
      }

      DATA[i1]->TimesEl(*B.DATA[i2]);
      if (cgflag) {
         RTD c2fac=1;
         for (j=0; j<CGR.dim2; ++j) {
            c2fac*=CGR(i1,j).sumTimesEl(FL,B.CGR(i2,j));
            CGR(i1,j).initCtrScalar();
         }
         (*DATA[i1]) *= c2fac;
      }

      if ((++mark[i1])>1) wblog(FL,
      "ERR %s() mark[%d]=%d ???", FCT, i1, mark[i1]);
   }

   for (k=i=0; i<DATA.len; ++i) {
      if (mark[i]) {
         if (k!=i) recSetQ(k,i);
         ++k;
      }
   }

   if (k) {
      QIDX.dim1=k; if (cgflag) { CGR.dim1=k; }
      DATA.len=k;
   }
   else clearQSpace();
};


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::plus_plain(
   const QSpace &B, QSpace &C, TD bfac, char vflag
 ) const {

   unsigned i,j; int i1, i2, cgflag=0;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;

   if (B.isEmpty()) { C=*this; return C; }
   if (  isEmpty()) { C=B; C*=bfac; return C; }

   if (qtype.len || B.qtype.len) {
      if (qtype==B.qtype) cgflag=1; else
      if ((qtype.len && B.qtype.len) || !isAbelian() || !B.isAbelian())
         wblog(FL,"ERR %s() qtype inconsistency ('%s' vs. '%s')",
         FCT,qStr().data,B.qStr().data
      ); 
   }

   if (getQOverlap(B, QQ, IQ)) wberror(FL,str);

   C.init(QQ, isEmpty() ? B.QDIM : QDIM);
   C.initQ(FL,*this,1,&B);

   if (itags.len || B.itags.len) {
      if (!B.itags.len) C.itags=  itags; else
      if (!  itags.len) C.itags=B.itags;
      else {
         if (itags!=B.itags) wblog(FL,
            "ERR %s() got different iLabels\n%s <> %s",
            FCT,IT2STR__,IT2STR(B));
         C.itags=itags;
      }
   }
   
   checkQ(FL,B, cgflag? 0 : 'a');
   C.setupCGR();

   if (CGR.dim2!=B.CGR.dim2 || CGR.dim2!=C.CGR.dim2) wblog(FL,
      "ERR cdata inconsistency (%d/%d/%d)",CGR.dim2,B.CGR.dim2,C.CGR.dim2);
   if (CGR.dim1 && CGR.dim1!=QIDX.dim1) wblog(FL,
      "ERR cdata inconsistency (%d/%d)", CGR.dim1, QIDX.dim1);
   if (B.CGR.dim1 && B.CGR.dim1!=B.QIDX.dim1) wblog(FL,
      "ERR cdata inconsistency (%d/%d)", B.CGR.dim1, B.QIDX.dim1);
   if (C.CGR.dim1 && C.CGR.dim1!=C.QIDX.dim1) wblog(FL,
      "ERR cdata inconsistency (%d/%d)", C.CGR.dim1, C.QIDX.dim1
   );

   unsigned n1=0, n2=0, n12=0;
   double bfac_c=1;

   for (i=0; i<IQ.dim1; ++i) { i1=IQ(i,0); i2=IQ(i,1);
       
       if (i1>=0 && i2>=0) { ++n12; bfac_c=1;
          if (cgflag) {
             for (j=0; j<CGR.dim2; ++j)
             bfac_c *= CGR(i1,j).safeCpy(FL, B.CGR(i2,j), C.CGR(i,j));
          }

          DATA[i1]->plus(*B.DATA[i2], *C.DATA[i], bfac*bfac_c);
       }
       else {
          if (i1>=0) { ++n1;
             C.DATA[i]->set(*DATA[i1]); if (cgflag) {
             for (j=0; j<CGR.dim2; ++j) C.CGR(i,j)=CGR(i1,j); }
          }
          else if (i2>=0) { ++n2;
             C.DATA[i]->set(*B.DATA[i2], bfac); if (cgflag) {
             for (j=0; j<CGR.dim2; ++j) C.CGR(i,j)=B.CGR(i2,j); }
          }
          else wblog(FL,"ERR i1=%d, i2=%d ???", i1, i2);
       }
   }

   if (vflag && n12) wblog(FL,
      "  * %s() w/out any outer multiplicity (%d+%d+%d)",
      FCT,n1,n2,n12); 

if (CGR.data && CGR.data[0].isEmpty()) wblog(FL,"ERR %s()",FCT);

   return C;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::makeUnique(char omFlag) {
   if (!hasMP())
        { return makeUnique_plain(omFlag); }
   else { return makeUnique_ortho(omFlag); }
};


template <class TQ, class TD>
unsigned QSpace<TQ,TD>::makeUnique_plain(char omFlag) {

   unsigned ib,j,j0,k,l, m=0, m1=0, m2=0, mx=0, mb=0, db,r;
   PERM_T const* p; bool cgflag=(!CGR.isEmpty());
   wbvector<INDEX_T> D,I;
   wbperm P;

   if (isEmpty()) { return mb; }

   if (QIDX.dim1!=DATA.len || !QDIM || QIDX.dim2%QDIM) wblog(FL,
      "ERR severe QSpace size mismatch (%d/%d,%d)",
       DATA.len, QIDX.dim1, QDIM
   );

   if (omFlag<0) { omFlag=(cgflag ? 'Q': 0); }
   else if (omFlag) { if (!cgflag) omFlag=0; else {
      if (CGR.dim1!=QIDX.dim1 || CGR.dim2!=qtype.len) wblog(FL,
         "ERR cdata inconsistency (%dx%d/%d <> %d/%d",
          CGR.dim1, CGR.dim2, qtype.len, QIDX.dim1, QDIM
      ); 
      NormCG();
   }}

   if (omFlag) {
      wbMatrix<TQ> X(QIDX); X.groupRecs(P,D);
   }
   else {
      QIDX.groupRecs(P,D);
   }

   p=P.data; I.init(DATA.len);

   for (db=k=l=ib=0; ib<D.len; ++ib, k+=db) { db=D[ib];
      wbvector<char> mark(db); mark.set(1);

      for (j0=db, j=0; j<db; ++j) {
         if (gotZeroData(p[k+j])) { mark[j]=0; ++mx; }
         else if (j0>j) { j0=j; }
      }
      if (j0==db) { j0=0; --mx; }

      j=mark.nnz(); if (j>1) m+=(j-1);


      for (m1=0; j0<db; ++m1) {
         wbarray<TD> &a0=(*DATA[p[k+j0]]);
         I[l++]=p[k+j0]; mark[j0]=0;


         for (j=j0+1; j<db; ++j) { if (!mark[j]) continue;
            wbarray<TD> &a2=(*DATA[p[k+j]]);

            if (j0==0 && a2.SIZE!=a0.SIZE) {
#ifdef MATLAB_MEX_FILE
               MXPut(FL).add(*this,"A");
#endif
               wblog(FL,"ERR severe size inconsistency "
               "(%d: %d,%d): %s + %s = ??", ib+1,k+j0+1,k+j+1,
               a2.sizeStr().data, a0.sizeStr().data);
            }

            if (!omFlag) {
               if (CGR.data) { double e2=0;
                  const CRef<TQ> *c0=CGR.ref(p[k+j0]);
                  const CRef<TQ> *c2=CGR.ref(p[k+j ]);
                  for (r=0; r<CGR.dim2; ++r) {
                     e2=sqrt(double(c0[r].normDiff2(FL,c2[r])));
                     if (e2>1E-12) { wblog(FL,
                        "ERR %s() got non-identical CGC(%d,:) <> CGC(%d,:)",
                        FCT,p[k+j0]+1,p[k+j]);
                     }
                  }
               }

               a0+=a2; a2.init(); mark[j]=0; ++m2;
               if (cgflag) {
                  CRef<TQ> *c2=CGR.ref(p[k+j]);
                  for (r=0; r<CGR.dim2; ++r) c2[r].init();
               }
               continue;
            }
            else {
               CRef<TQ> *c0=CGR.ref(p[k+j0]);
               CRef<TQ> *c2=CGR.ref(p[k+j ]); double e2=0;

               for (r=0; r<CGR.dim2; ++r) {
                  if (c0[r].cgr!=c2[r].cgr) wblog(FL,
                     "ERR CGR data mismatch (%d,%d)\n%s <> %s",
                     k+j0+1, k+j+1, c0[r].toStr().data, c2[r].toStr().data
                  );

                  if (e2<1E-12)
                  e2=sqrt(double(c0[r].normDiff2(FL,c2[r])));
               }
               if (e2<1E-12) {
                  a0+=a2; a2.init(); ++m2;
                  for (r=0; r<CGR.dim2; ++r) c2[r].init();
                  mark[j]=0;
               }
               else wblog(FL,"ERR %s() got e2=%g",FCT,e2);
            }
         }

         for (; j0<db; ++j0) { if (mark[j0]) break; }

         if (j0<db && !omFlag) wblog(FL,
            "ERR %s() failed to add up all blocks (CG: %d)\n"
            "ERR hint: got outer multiplicity?",FCT,mark.nnz()
         );
      }
      if (mb<m1) mb=m1;
   }

   if (m!=m2 && omFlag!='Q') {"\n"); 
       wblog(FL,"  * got %d/%d block%s with identical Q-labels "
         "(%d;%d/%d) \r\\", m-m2,DATA.len-m2-mx, m!=1 ? "s":"",m2,mx,
         DATA.len); if (omFlag=='v' || omFlag=='V') printf("\n");
   }

   if (!l) wblog(FL,"ERR l=%d !??",l);
   else I.len=l;

   this->Select((wbindex&)I, omFlag ? 0:'d');

   return mb;
};



template <class TQ, class TD>
unsigned QSpace<TQ,TD>::makeUnique_ortho(char vflag) {

   if (isEmpty() || !isConsistent(FL)) { return 0; }
   if (!hasMP()) { wblog(FL,"WRN %s() no OM expected "
        "(%s; %d; om=%d)\nusing makeUnique_plain() instead",
         FCT, qtype.toStr().data, rank(FL), gotCGS(FL));
      return makeUnique_plain();
   }

   QSpace<TQ,TD> C;
   wbvector<INDEX_T> D,Dc;
   wbperm P;

   C.initQ(*this); C.QIDX=QIDX;
   C.QIDX.groupRecs(P,D); D.cumsum0(Dc,'x');
   if (C.QIDX.dim1==QIDX.dim1) return 1;

   unsigned ic,N,nc=D.len, ep=0, nz=0;
   wbMatrix<unsigned> om(nc,qtype.len);
   wbvector<unsigned> OM(nc), OMc;
   wbMatrix< OrthoCGS<TQ> > M(nc,CGR.dim2);

   for (ic=0; ic<D.len; ++ic) {
      if (D[ic]==1) {
         OM[ic]=1; continue;
      }

      const wbvector<INDEX_T> ia(D[ic],P.data+Dc[ic],'r');

         for (unsigned j=0; j<qtype.len; ++j) {
            om(ic,j)=M(ic,j).getCGS_ortho(FL,*this,ia,j);
            if (!om(ic,j)) wblog(FL,
               "WRN %s() got om(%d,%d)=%d",FCT,ic,j,om(ic,j));
         }; OM[ic]=om.recProd(ic);
   }
   if (ep) wblog(FL,"ERR %s() ep=%d",FCT,ep); 

   N=OM.cumsum0(OMc);

   C.QIDX.Resize2Mult(OM);
   if (C.QIDX.dim1!=N) wblog(FL,
      "ERR %s() size inconsistency %d/%d",FCT,N,QIDX.dim1);
   C.setupDATA().setupCGR();

   for (ic=0; ic<nc; ++ic) {
      if (D[ic]==1) {
         INDEX_T j=0, ia=P[Dc[ic]], l=OMc[ic];
         DATA[ia]->save2(*C.DATA[l]); { for (; j<CGR.dim2; ++j)
         CGR(ia,j).save2(C.CGR(l,j)); }
      }
      else if (D[ic]>1) {
         const wbvector<INDEX_T> ia(D[ic],P.data+Dc[ic],'r');
         const wbvector< OrthoCGS<TQ> > Mc(M.dim2,M.rec(ic),'r');
          
            nz+=combineDATA_ortho(FL,*this,ia, Mc, C, OMc[ic]);
      }
   }
   if (ep) wblog(FL,"ERR %s() ep=%d",FCT,ep); 
   if (nz) C.SkipZeroData();

   C.save2(*this);

   N=OM.max();
   if (vflag && N>4) wblog(FL,"WRN %s() got OM=%d",FCT,N);

   return N;
};


template <class TQ, class TD> inline
TD QSpace<TQ,TD>::trace() const {

    unsigned i,j,K,s; int r=-1; gotCGS(FL);
    TQ *q=QIDX.data;
    TD tsum=0,t=0;

    if (isEmpty()) return 0;
    isConsistent(r,FL);

    if (r%2) wblog(FL,"ERR %s() requires even-rank (%d)",FCT,r);

    K=QDIM*(unsigned)r/2; s=K*sizeof(TQ);

    for (i=0; i<DATA.len; ++i, q+=QIDX.dim2)
    if (!memcmp(q, q+K, s)) {
       t=DATA[i]->trace();
       for (j=0; j<CGR.dim2; ++j) { t*=CGR(i,j).trace(); }
       tsum+=t;
    }

    return tsum;
};


template <class TQ, class TD>
double QSpace<TQ,TD>::maxDiff(const QSpace<TQ,TD> &B) const {

   unsigned i,m=0; int i1, i2;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   double dmax=0, d=0;

   if (QDIM!=B.QDIM || QIDX.dim2!=B.QIDX.dim2) {
      sprintf(str,"%s:%d QDIM mismatch (%d,%d; %d,%d)", FL,
      QDIM, B.QDIM, QIDX.dim2, B.QIDX.dim2); return NAN;
   }

   if (getQOverlap(B, QQ, IQ)) wberror(FL,str);

   for (i=0; i<IQ.dim1; ++i) {
       i1=IQ(i,0); i2=IQ(i,1);
       
       if (i1>=0 && i2>=0) {
          d=DATA[i1]->maxDiff(*B.DATA[i2]); ++m;
          if (ISNAN(d)) return d;
       }
       else {
          if (i1>=0) d=  DATA[i1]->aMax(); else
          if (i2>=0) d=B.DATA[i2]->aMax(); else
          wblog(FL,"ERR i1=%d, i2=%d ???", i1, i2);
       }
       dmax=MAX(dmax,d);
   }

   if (!m) wblog(FL, "WRN no Q-overlap");

   return dmax;
}


template <class TQ, class TD>
TD QSpace<TQ,TD>::normDiff2(const QSpace<TQ,TD> &B) const {

   unsigned i;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   TD x=0;

   if (isEmpty() || B.isEmpty()) {
      wblog(FL,"WRN %s() with empty object (%d,%d) !?",
      FCT, isEmpty(), B.isEmpty()); return x;
   }

   if (this==&B) return x;

   if (QIDX!=B.QIDX) {
      int i1,i2;

      if (getQOverlap(B, QQ, IQ)) wberror(FL,str);

      for (i=0; i<IQ.dim1; ++i) {
          i1=IQ(i,0); i2=IQ(i,1);

          if (i1>=0 && i2>=0) 
             x+=(DATA[i1]->normDiff2(*B.DATA[i2]));
          else {
             if (i1>=0) x+=(  DATA[i1]->norm2()); else
             if (i2>=0) x+=(B.DATA[i2]->norm2());
             else wblog(FL,"ERR i1=%d, i2=%d ???", i1, i2);
          }
      }
   }
   else {
      for (i=0; i<DATA.len; ++i)
      x+=DATA[i]->normDiff2(*B.DATA[i]);
   }

   return x;
};


template <class TQ, class TD>
TD QSpace<TQ,TD>::sumData() const {

    TD x=0;
    for (unsigned i=0; i<DATA.len; ++i) x+=(DATA[i]->sum());
    return x;
}


template <class TQ, class TD> inline 
void QSpace<TQ,TD>::contractMat(
   unsigned ia, const QSpace<TQ,TD> &B, unsigned ib,
   QSpace<TQ,TD> &C
) const {

   unsigned i, ra=rank(FL), rb=B.rank(FL);

   if (ia==0 || ib==0) wblog(FL,
      "ERR contractMat index is 1-based (%d,%d)!", ia, ib);
   if (ia>ra || ib>rb) wblog(FL,
      "ERR contractMat index out of range (%d/%d; %d/%d)!", ia,ra,ib,rb);
   if (rb!=2) wblog(FL,
      "ERR contractMat requires 2-D object! (%d)", rb);

   ctrIdx ica(1,&(--ia)), icb(1,&(--ib));

   wbperm P(ra);
   i=ia; P[i++]=ra; for (; i<ra; ++i) P[i]=i-1;
   
   contract(ica, B, icb, C, P);
};


#ifdef WB_CLOCK
   WbClock wbc_qs_cidx("QS::contract.getIdx");
#endif


template <class TQ, class TD>
template <class TB, class TC>
int QSpace<TQ,TD>::contract_getIdxSet(const char *F, int L,
   const ctrIdx &ica,
   const QSpace<TQ,TB> &B, const ctrIdx &icb,
   wbindex &Ia, wbindex &Ib, wbvector<INDEX_T> &Dc,
   QSpace<TQ,TC> &C
 ) const {

   unsigned i,l, ra=rank(F_L), rb=B.rank(F_L); int e;
   wbindex ika, ikb;
   wbMatrix<TQ> QAc, QAk, QBc, QBk;

   wbMatrix<TQ> &QC=C.QIDX;
   IndexLabels &idC=C.itags;

#ifdef WB_CLOCK
   wbc_qs_cidx.resume();
#endif

   Ia.init(); Ib.init();

   if (QDIM!=B.QDIM) wblog(F_L,
      "ERR contract() incompatible objects (QDIM=%d/%d)",
       F ? "QSpace::" : "", QDIM, B.QDIM); 
   if (ica.len!=icb.len || !ica.len) wblog(F_L,
      "ERR invalid contraction [%s] <> [%s]",
       ica.toStr().data, icb.toStr().data
   );

   if (!isUniqueIdxSet(ica,ra) || !isUniqueIdxSet(icb,rb) ||
       ica.len!=icb.len) wblog(F_L,"ERR invalid contraction indices\n"
      "[%s; %d], [%s; %d] (QDIM=%d/%d)",
       ica.toStr().data, ra, icb.toStr().data, rb, QDIM, B.QDIM
   );

   const wbindex ica_(ica), icb_(icb);

   ica_.invert(ra,ika);   getQsub(ica_,QAc,QAk);
   icb_.invert(rb,ikb); B.getQsub(icb_,QBc,QBk);

   e=matchIndex(QAc,QBc,Ia,Ib);

   QC.Cat(2,QAk,Ia,QBk,Ib);


   if (!ika.isEmpty() || !ikb.isEmpty()
     || Ia.isEmpty()
   ){
      wbperm P; QC.groupRecs(P,Dc);
      Ia.Permute(P); Ib.Permute(P);
   }
   else {
      Dc.init(1); Dc[0]=Ia.len;
      QC.init(1,0);
   }

   if (itags.len || B.itags.len) {
      if (itags.len && itags.len!=ra) wblog(FL,
         "ERR %s() length mismatch in index-labels (A: %d/%d)",
          FCT,itags.len, ra);
      if (B.itags.len && B.itags.len!=rb) wblog(FL,
         "ERR %s() length mismatch in index-labels (B: %d/%d)",
          FCT,B.itags.len, rb
      );

      idC.init(ika.len+ikb.len);

      if (itags.len) {
         for (i=0; i<ika.len; ++i) idC[i]=itags[ika[i]];
      }
      else {
         char s[8]; i=0;

         if (ika.len>=ica.len) {
            for (; i<icb.len; ++i) idC[i]=B.itags[icb[i]];
         }

         for (; i<ika.len; ++i) {
            if (ika[i]==2 && otype==QS_OPERATOR && rank()==3)
                 strcpy(s,"op");
            else sprintf(s,"a%03d",int(ika[i]+1));
            idC[i].init(FL,s);
         }
      }

      if (B.itags.len) {
         for (l=i, i=0; i<ikb.len; ++i) idC[l+i]=B.itags[ikb[i]];
      }
      else {
         char s[8]; l=i; i=0;

         if (ikb.len>=icb.len) {
            for (; i<ica.len; ++i) idC[l+i]=itags[ica[i]];
         }

         for (; i<ikb.len; ++i) {
            if (ikb[i]==2 && B.otype==QS_OPERATOR && B.rank()==3)
                 strcpy(s,"op");
            else sprintf(s,"b%03d",int(ikb[i]+1));
            idC[l+i].init(FL,s);
         }
      }

      if (itags.len && B.itags.len) { char zflag=0;
         for (i=0; i<ica.len; ++i)  {
            if (!itags.sameAs(ica[i], B.itags[icb[i]], &zflag))
            wblog(FL,"ERR contract() got itag mismatch\n"
              "%s [%s] <> %s [%s] @ i=%d %N",
               IT2STR__, (ica+1).toStr().data,
               IT2STR(B),(icb+1).toStr().data, i+1
            );
         }

         if (zflag<0) {
            if (!(ica.conj^icb.conj)) {
               wblog(FL,"WRN contract() applying implicit conjA flag");
               ((ctrIdx&)ica).Conj();
            }
         }
         if (ica.conj) { for (i=0;       i<ika.len; ++i) idC[i].Conj(); }
         if (icb.conj) { for (i=ika.len; i<idC.len; ++i) idC[i].Conj(); }
      }

   }

#ifdef WB_CLOCK
   wbc_qs_cidx.stop();
#endif

   return e;
};


#ifdef WB_CLOCK
 WbClock wbc_qs_cact("QS::contract.actual");
#endif

template <class TQ, class TD>
template <class TDB, class TDC>
QSpace<TQ,TDC>& QSpace<TQ,TD>::contract(const char *F, int L,
   const ctrIdx &ica, const QSpace<TQ,TDB> &B,
   const ctrIdx &icb,
   QSpace<TQ,TDC> &C, const wbperm &P,
   char cgflag
) const {

   if ((void*)this==(void*)&C || (void*)&B==(void*)&C) {
      QSpace<TQ,TDC> X;
      contract(F,L,ica,B,icb,X,P,cgflag); X.save2(C);
      return C;
   }

#ifdef WB_CLOCK
   wbc_qs_cact.resume();
#endif

   unsigned N,nc,ic,ra,rb,rc, nx=0, ep=0;
   wbvector<INDEX_T> D,Dc;
   wbindex Ia,Ib;

   C.clearQSpace(); if (isEmpty() || B.isEmpty()) return C;

   if (QDIM!=B.QDIM) wblog(F,L,
      "ERR QSpace() incompatible contraction (QDIM: %d/%d)",
       QDIM, B.QDIM);

   ra=rank(FL); rb=B.rank(FL); rc=(ra+rb)-ica.len-icb.len;

   if (!P.isEmpty() && (!validPerm(P) || P.len!=rc)) wblog(FL,
      "ERR %s() invalid permutation [%s; %d]",
      FCT,P.toStr().data,rc
   );

   contract_getIdxSet(FL,ica,
      B, icb, Ia,Ib,D,
      C
   );
   D.cumsum0(Dc,'x');

   if ((qtype.len || B.qtype.len) && qtype!=B.qtype) {
      if (isAbelian() && B.isAbelian() && (!qtype.len || !B.qtype.len)) {
         if (cgflag<0) cgflag=0; else
         if (cgflag) wblog(FL,
            "ERR %s() got cgflag=%d for all-abelian",FCT,cgflag
         );
      }
      else wblog(FL,"ERR %s() qtype inconsistency ('%s' vs. '%s')",
         FCT, qStr().data, B.qStr().data
      ); 
   }

   if (cgflag<0) cgflag=gotCGS(FL);
   
   C.initQ(F,L,*this,0, cgflag<=0? NULL : &B);

   nc=D.len;

   if (cgflag<=0) {
      C.qtype.init(); C.CGR.init();
      C.setupDATA();

      for (unsigned ic=0; ic<nc; ++ic) {
         INDEX_T i0=Dc[ic], d=D[ic];
         const wbvector<INDEX_T>
            ia(d,Ia.data+i0,'r'),
            ib(d,Ib.data+i0,'r');

         try {
            contractDATA_plain(FL,*this,ia,ica,B,ib,icb,C,ic);
         }
         catch (...) { ++ep; }
      }

      if (ep) wblog(FL,"ERR %s() ep=%d",FCT,ep); 
   }
   else {
      wbMatrix<unsigned> om(nc,qtype.len);
      wbvector<unsigned> OM(nc), OMc;
      wbMatrix< OrthoCGS<TQ> > M(nc,qtype.len);


      for (ic=0; ic<nc; ++ic) {
         INDEX_T i0=Dc[ic], d=D[ic];
         const wbvector<INDEX_T>
            ia(d,Ia.data+i0,'r'), ib(d,Ib.data+i0,'r');
          
            for (unsigned j=0; j<qtype.len; ++j) {
               om(ic,j)=M(ic,j).contractCGS_ortho(
                  FL,*this,ia,ica,B,ib,icb,j
               );
            }
            
            OM[ic]=om.recProd(ic);
      }
      if (ep) wblog(FL,"ERR %s() ep=%d",FCT,ep); 

      N=OM.cumsum0(OMc);

      C.QIDX.Resize2Mult(OM);
      if (C.QIDX.dim1!=N) wblog(FL,
         "ERR %s() size inconsistency %d/%d",FCT,N,C.QIDX.dim1);
      C.setupDATA().setupCGR();


      for (ic=0; ic<nc; ++ic) { if (OM[ic]) {

         INDEX_T i0=Dc[ic], d=D[ic];
         const wbvector<INDEX_T>
            ia(d,Ia.data+i0,'r'), ib(d,Ib.data+i0,'r');
         const wbvector< OrthoCGS<TQ> > Mc(M.dim2,M.rec(ic),'r');
          
            nx+=contractDATA_ortho(FL,
               *this,ia,ica, B,ib,icb, Mc, C, OMc[ic]
            );
      }}
      if (ep) wblog(FL,"ERR %s() ep=%d",FCT,ep); 
   }

   if (nx) {
      C.SkipEmptyData(FL);
   }

   C.checkScalarCGS(FL);
   C.Permute(P);

   if (C.rank(FL)<3) {
      C.NormCG();
   }


#ifdef WB_CLOCK
   wbc_qs_cact.stop();
#endif

   return C;
};


template <class TQ, class TD>
void QSpace<TQ,TD>::contract(
   unsigned i1, unsigned i2,
   QSpace<TQ,TD> &C
) const {

   unsigned i,k,m,n,N, r=rank(FL);
   wbindex ik, Iu, ic(2);
   wbMatrix<TQ> Qc, Qk;

   if (i1==i2 || !i1 || !i2 || i1>r || i2>r) wblog(FL,
      "ERR contract() got trace indices (%d,%d / %d)",i1,i2,r);

   ic[0]=i1-1; ic[1]=i2-1;
   ic.invert(r,ik);

   getQsub(ic,Qc,Qk);

   if (!ik.isEmpty()) {
      Qk.makeUnique_ig(Iu);
   }
   else {
      Qk.init(1,0); n=QIDX.dim1;
      Iu.init(n+1); Iu[0]=n; for (i=0; i<n; ++i) Iu[i+1]=i;
   }

   C.init(Qk,QDIM); N=Qk.dim1;

   for (i=n=0; n<N; ++n) {
       for (m=Iu[i++], k=0; k<m; ++k, ++i) {
           DATA[Iu[i]]->contract(i1,i2, *(C.DATA[n]));
       }
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::checkScalarCGS(const char *F, int L) {

   if (CGR.data && qtype.len!=CGR.dim2) wblog(F_L,"ERR %s() "
      "size mismatch (%dx%d/%d)",FCT,CGR.dim1,CGR.dim2,qtype.len);
   if (!CGR.dim2) return;
  
   unsigned i=0,j=0;
   char isa[qtype.len];
   
   for (; j<CGR.dim2; ++j) {
      if ((isa[j]=qtype[j].isAbelian())!=0) ++i;
   }

   if (i) {
      for (i=0; i<CGR.dim1; ++i)
      for (j=0; j<CGR.dim2; ++j) { if (isa[j]) CGR(i,j).checkAbelian(FL); }
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::disp_cgs(
   const char *F, int L, const char *istr
 ) const {

   unsigned i,j; wblog(F_L,"TST CGR data (%s)",istr);  
   for (i=0; i<CGR.dim1; ++i) {
   for (j=0; j<CGR.dim2; ++j) {
      printf(" %8s",CGR(i,j).sizeStr().data);
   }; printf("\n"); }

};


template <class TQ, class TD>
wbstring QSpace<TQ,TD>::sizeStrQ(const char *F, int L) const {
   size_t l=0, n=32; char s[n];
   isConsistent(F_L);
   if (QDIM)
        l=snprintf(s,n,"%dx(%dx%d)",QIDX.dim1,QIDX.dim2/QDIM,QDIM);
   else l=snprintf(s,n,"%dx%d/%d",QIDX.dim1,QIDX.dim2,QDIM);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};


template <class TQ, class TD>
wbstring QSpace<TQ,TD>::sizeStr(char vflag) const {
   if (isEmpty()) { return wbstring(); }
   if (isAbelian()) {
      wbvector<INDEX_T> D; getDim(D);
      return D.toStrD();
   }
   else {
      wbvector<INDEX_T> D,DD; getDim(D,&DD);
      if (vflag) { char s[64]; int i=
         snprintf(s,64,"%s (%s)",D.toStrD().data,DD.toStrD().data);
         if (i>=64) wblog(FL,
            "ERR %s() got s='%s...' !??  (%d)",FCT,s,i);
         return s;
      }
      else { return DD.toStrD(); }
   }
};

 
template <class TQ, class TD> inline
void QSpace<TQ,TD>::init_gRG_zdim_bare(const char *F, int L) const {

    if (!QIDX.dim1 || !QIDX.dim2) return;

    if (gotCGS(F,L)<=1) { return; }

    unsigned i,j,s, n, r=QIDX.dim2/QDIM;
    const TQ *q=QIDX.data;
    qset<TQ> J;

    for (i=0; i<QIDX.dim1; ++i)
    for (j=0; j<r; ++j) {
        for (s=0; s<qtype.len; ++s, q+=n) {
           n=qtype[s].qlen();
           if (n<=1) continue;
           if (!CGR(i,s).cgr) wblog(FL,
              "ERR %s() cgr required for %s",FCT,qtype[j].toStr().data);

           J.init(n,q);

           genRG_base<TQ,RTD> &R=gRG.buf[qtype[s]].RSet[J];
           wbvector<unsigned> S; CGR(i,s).getSize(S);

           if (R.J.isEmpty() && R.Z.isEmpty()) {
              if (j>=r || r!=S.len) wblog(FL,
                 "ERR %s() index out of bounds (%d/%d;%d)",FCT,j,r,S.len); 
              R.J=J; R.Z.init(S[j],n).set(99);
           }
           else if (R.J!=J || R.Z.dim1!=S[j] || R.Z.dim2!=n) { wblog(FL,
             "ERR %s() severe size mismatch (%d/%d) => %d [%s]\n"
             "%dx%d <> %dx%d", FCT, j,r,s,qtype[s].toStr().data,
              R.Z.dim1, R.Z.dim2, S[j], n);
           }
        }
    }
};








template <class TQ, class TD>
void QSpace<TQ,TD>::EigenSymmetric(
   QSpace<TQ,TD>      &AK, QSpace<TQ,TD>     &AT,
   QSpace<TQ,double>  &EK, QSpace<TQ,double> &ET,
   wbvector<double>   &Etot,
   wbMatrix<unsigned> &DD,
   int &Nkeep,
   double Etrunc,
   double *E0,
   const wbperm &PA,
   double deps,
   double db,
   int dmax,
   const char *sdir
 ) const {

   if ((void*)&AK==(void*)this || (void*)&AT==(void*)this || 
       (void*)&EK==(void*)this || (void*)&ET==(void*)this ){

      unsigned e=0; QSpace<TQ,TD> X;
      if ((void*)&AK==(void*)this) { AK.save2(X); ++e; }
      if ((void*)&AT==(void*)this) { AT.save2(X); ++e; }
      if ((void*)&EK==(void*)this) { EK.save2(X); ++e; }
      if ((void*)&ET==(void*)this) { ET.save2(X); ++e; }

      if (e!=1) wblog(FL,
         "ERR multiple overlap within output space (%d) !??",e);

      X.EigenSymmetric(
        AK,AT,EK,ET,Etot,DD, Nkeep,Etrunc,E0,PA,deps,db,dmax,sdir);
      return;
   }

   unsigned i,d,n,l,nk,nt,r2; int r=-1; bool cgflag=(gotCGS(FL)>0);
   wbMatrix<TQ> Q1,Q2,Qtot;
   wbarray<TD> MM,UU;
   wbMatrix<INDEX_T> S1(1,1), S2;
   wbvector<INDEX_T> D,D1(1), D2;
   wbvector<TD> Ek,Et;
   QSpace<TQ,TD> A1,E1;

   wbindex Ik,It;
   wbperm P;

   wbvector< EIG_QS<TQ,TD> > EIG;
   wbvector< wbvector<char> > mark;
   wbvector< wbvector<TD>* > ESET;

   if (!isConsistent(r)) wberror(FL,str);
   if (r%2) wblog(FL,"ERR %s() requires even-rank (%d)",FCT,r);
   r2=unsigned(r)/2;

   EK.init(); AK.init(); ET.init(); AT.init(); 
   if (isEmpty()) return;

   if (itags.len) {
      if (!itags.got_op_labels(FL) || unsigned(r)!=itags.len) wblog(FL,
         "ERR %s() got invalid operator labels '%s' (%d)",
         FCT,itags.toStr().data,r
      );
   }
   else {
      if (!CGR.isEmpty() || (qtype.len && !qtype.allAbelian()))
      wblog(FL,"ERR %s() got missing itags (empty)",FCT); 
   }

   if (cgflag) { RTD x;
      if (r2>1) wblog(FL,
         "ERR %s() not implemented yet for rank-%d (qtype '%s')",
          FCT,r,qStr().data
      ); 

      if (!isBlockDiagMatrix(FL)) {
         MXPut(FL).add(*this,"H");
         wblog(FL,"ERR %s() H is not block-diagonal (%d)",FCT,cgflag);
      }
     
      for (n=CGR.numel(), i=0; i<n; ++i) {
         if (CGR[i].isScalar()) continue;
         if (!CGR[i].isIdentityCGC(&x,1E-12)) {
            MXPut(FL,"q").add(*this,"H").add(Nkeep,"Nkeep")
             .add(Etrunc,"Etrunc").add(i+1,"i").add(n,"n");
            wblog(FL,
             "ERR expecting identity as CGR coefficients (%s; %d/%d)",
              CGR[i].cgr ? CGR[i].cgr->sizeStr().data : "???",i+1,n
           ); 
         }
      }
   }

   getQsum( wbindex(r2,'i'), Q1).groupRecs(P,D);

   ESET.init(D.len);
   EIG.init(D.len);

   if (cgflag) DD.init(D.len,2); else DD.init(D.len,1);

   for (l=i=0; i<D.len; ++i, l+=d) { d=D[i];

       getSub(wbindex(d,P.data+l),A1,'r')
       .toBlockMatrix(MM, r2,
          EIG[i].Q, Q2,
          EIG[i].S, S2,
          EIG[i].D, D2
       );

       if (EIG[i].Q!=Q2 || EIG[i].S!=S2) {
          MXPut(FL,"IH").add(*this,"this").add(P,"P").add(D,"D")
          .add(EIG[i].Q,"Q1").add(Q2,"Q2").add(EIG[i].S,"S1").add(S2,"S2")
          .add(EIG[i].D,"D1").add(D2,"D2").add(wbindex(d,P.data+l),"idx");

          if (!isHConj())
               wblog(FL,"ERR %s() requires symmetric QSpace!",FCT);
          else wblog(FL,"ERR %s() Q sectors couple to others\n"
          "=> quantum symmetries not preserved !?", FCT);
       }

       wbEigenS(MM, EIG[i].U, EIG[i].E);

       ESET[i]=&(EIG[i].E); DD(i,0)=(MM.SIZE.len>1 ? MM.SIZE[1]: 0);
   }

   if (cgflag) {
      wbvector<INDEX_T> S1,Sc; wbindex Ia,Ib; INDEX_T i1,i2;
      wbMatrix<INDEX_T> SC;
      getQDim(Q2,S1,&SC); SC.recProd(Sc);

      Q1.init(EIG.len,Q2.dim2);
      for (i=0; i<D.len; ++i) {
         if (EIG[i].Q.dim2!=Q1.dim2) wblog(FL,"ERR %s() "
            "size mismatch (%dx%d/%d",FCT,EIG[i].Q.dim2,Q1.dim2);
         Q1.recSetP(i,EIG[i].Q.data);
      }

      i=matchIndex(Q1,Q2,Ia,Ib);
      if (i || Ia.len!=Q1.dim1) wblog(FL,"ERR %s() "
         "failed to match q-labels (%d/%d)",FCT,Ia.len,Q1.dim1);

      for (i=0; i<D.len; ++i) { i1=Ia[i]; i2=Ib[i];
         DD(i1,1)=DD(i1,0)*Sc[i2];
      }
   }

   markSet(ESET,mark,Nkeep,Etrunc,Etot,deps,db,dmax,sdir);
   nk=nt=0;

   if (E0) { (*E0)=Etot[0]; Etot-=(*E0); }

   for (l=i=0; i<D.len; ++i, l+=d) { d=D[i];

       EIG[i].Q.blockSum(QDIM,Qtot);
       if (!Qtot.recAllEqual()) {
          EIG[i].put("EIG"); wblog(FL,"ERR %d/%d ???", i+1, D.len); }
       Qtot.set2Rec(0);

       mark[i].find(Ik,It);
       EIG[i].E.select(Ik,Ek); if (E0) Ek-=(*E0);
       EIG[i].E.select(It,Et); if (E0) Et-=(*E0);

       S1[0]=D1[0]=Ik.len; nk+=Ik.len;
       if (Ik.len) {
          EIG[i].U.select0(Ik,1,UU);

          A1.initFromBlockMatrix(UU, EIG[i].Q, Qtot, QDIM, EIG[i].S, S1);
          E1.initDiagonal(Qtot, D1, Ek, 'r');

          A1.Append2AndDestroy(FL,AK);
          E1.Append2AndDestroy(FL,EK);
       }

       S1[0]=D1[0]=It.len; nt+=It.len;
       if (It.len) {
          EIG[i].U.select0(It,1,UU);

          A1.initFromBlockMatrix(UU, EIG[i].Q, Qtot, QDIM, EIG[i].S, S1);
          E1.initDiagonal(Qtot, D1, Et, 'r');

          A1.Append2AndDestroy(FL,AT);
          E1.Append2AndDestroy(FL,ET);
       }

#ifndef NSAFEGUARDS
       if (EIG[i].D.sum()!=mark[i].len || Ik.len+It.len!=mark[i].len) {
          wblog(FL,"ERR %d: %d %d == %d+%d",
          i, EIG[i].D.sum(), mark[i].len, Ik.len, It.len);
       }
#endif
   }

   if (nk+nt!=Etot.len) wblog(FL, 
      "ERR inconsistency in state count (%d+%d==%d)", nk,nt,Etot.len);

   if (cgflag>0) {
      AK.qtype=qtype; AK.initIdentityCGS(FL,*this);
      AT.qtype=qtype; AT.initIdentityCGS(FL,*this);
      EK.qtype=qtype; EK.initIdentityCGS(FL,*this);
      ET.qtype=qtype; ET.initIdentityCGS(FL,*this);
   }

   EK.itags.init_qdir("+-"); ET.itags=EK.itags;

   AK.itags.wbvector<iTag>::init(r2+1,itags.data);
      AK.itags.end()=EK.itags[1];
   AT.itags=AK.itags;


   if (!PA.isEmpty() && !PA.isIdentityPerm()) {
      if (!AK.isEmpty()) AK.Permute(PA);
      if (!AT.isEmpty()) AT.Permute(PA);
   }

   AK.Sort();
   AT.Sort();
   
   if (!AK.QIDX.isUniqueSorted()) 
   wblog(FL,"ERR AK does not have unique QIDX ???");
   
   if (!AT.QIDX.isUniqueSorted()) 
   wblog(FL,"ERR AT does not have unique QIDX ???");
};


template <class TQ, class TD>
void QSpace<TQ,TD>::Eigen_CSymmetric(
   QSpace<TQ,TD> &AK, QSpace<TQ,TD> &AT,
   QSpace<TQ,TD> &EK, QSpace<TQ,TD> &ET, wbvector<TD> &Etot,
   wbvector<unsigned> &DD,
   int &Nkeep,
   double Etrunc,
   TD *E0,
   const wbperm &PA,
   double deps,
   double db,
   int dmax,
   const char *sort
) const {

   if (&AK==this || &AT==this || &EK==this || &ET==this){

      unsigned e=0; QSpace<TQ,TD> X;
      if (&AK==this) { AK.save2(X); ++e; }
      if (&AT==this) { AT.save2(X); ++e; }
      if (&EK==this) { EK.save2(X); ++e; }
      if (&ET==this) { ET.save2(X); ++e; }

      if (e!=1) wblog(FL,
      "ERR multiple overlap within output space (%d) !??",e);

      X.Eigen_CSymmetric(
        AK,AT,EK,ET,Etot,DD,Nkeep,Etrunc,E0,PA,deps,db,dmax,sort);
      return;
   }

   unsigned i,d,r2,l,nk,nt; int r=-1;
   wbMatrix<TQ> Q1,Q2,Qtot;
   wbarray<TD> MM,UU;
   wbMatrix<unsigned> S1(1,1),S2;
   wbvector<unsigned> D,D1(1),D2;
   wbvector<TD> Ek,Et;
   QSpace<TQ,TD> A1,E1;

   wbindex Ik,It;
   wbperm P;

   wbvector< EIG_QS<TQ,TD> > EIG;
   wbvector< wbvector<char> > mark;
   wbvector< wbvector<TD>* > ESET;

   if (!isConsistent(r)) wberror(FL,str);
   r2=unsigned(r)/2;

   if (r%2) wblog(FL,
      "ERR EigenSymmetric() requires even-rank (%d)", r);

   EK.init(); AK.init(); ET.init(); AT.init(); 
   if (isEmpty()) return;

   getQsum( wbindex(r2,'i'), Q1).groupRecs(P,D);

   ESET.init(D.len);
   EIG.init(D.len); DD.init(D.len);

   for (l=i=0; i<D.len; ++i, l+=d) {
       d=D[i];

       getSub( wbindex(d,P.data+l),A1,"ref")
      .toBlockMatrix(MM, r2,
          EIG[i].Q, Q2,
          EIG[i].S, S2,
          EIG[i].D, D2
       );

       if (EIG[i].Q!=Q2 || EIG[i].S!=S2) {
#ifdef MATLAB_MEX_FILE
          mxStruct S(FL);
          S("this")=toMx(); S("P")=P; S("D")=D;
          S("Q1")=EIG[i].Q; S("Q2")=Q2;
          S("S1")=EIG[i].S; S("S2")=S2;
          S("D1")=EIG[i].D; S("D2")=D2;
          S("idx")=wbindex(d,P.data+l);

          S.put("Hstru");
#endif
          if (!isHConj())
          wblog(FL,"ERR %s - requires symmetric QSpace!",FCT); else
          wblog(FL,"ERR %s - Q sectors couple to others\n"
          "Quantum symmetries not preserved ???", FCT);
       }

       wbEigen_CS(MM, EIG[i].U, EIG[i].E,'R','T');

       ESET[i]=&(EIG[i].E); DD[i]=(MM.SIZE.len>1 ? MM.SIZE[1]: 0);
   }

   markSet(ESET,mark,Nkeep,Etrunc,Etot,deps,db,dmax,sort);
   nk=nt=0;

   if (E0) { (*E0)=Etot[0]; Etot-=(*E0); }

   for (l=i=0; i<D.len; ++i, l+=d) {
       d=D[i];

       EIG[i].Q.blockSum(QDIM,Qtot);
       if (!Qtot.recAllEqual()) {
          EIG[i].put("EIG"); wblog(FL,"ERR %d/%d ???", i+1, D.len); }
       Qtot.set2Rec(0);

       mark[i].find(Ik,It);
       EIG[i].E.select(Ik,Ek); if (E0) Ek-=(*E0);
       EIG[i].E.select(It,Et); if (E0) Et-=(*E0);

       S1[0]=D1[0]=Ik.len; nk+=Ik.len;
       if (Ik.len) {
          EIG[i].U.select0(Ik,1,UU);

          A1.initFromBlockMatrix(UU, EIG[i].Q, Qtot, QDIM, EIG[i].S, S1);
          E1.initDiagonal(Qtot, D1, Ek, 'r');

          A1.Append2AndDestroy(FL,AK);
          E1.Append2AndDestroy(FL,EK);
       }

       S1[0]=D1[0]=It.len; nt+=It.len;
       if (It.len) {
          EIG[i].U.select0(It,1,UU);

          A1.initFromBlockMatrix(UU, EIG[i].Q, Qtot, QDIM, EIG[i].S, S1);
          E1.initDiagonal(Qtot, D1, Et, 'r');

          A1.Append2AndDestroy(FL,AT);
          E1.Append2AndDestroy(FL,ET);
       }

       if (EIG[i].D.sum()!=mark[i].len || Ik.len+It.len!=mark[i].len)
       wblog(FL,"ERR %d: %d %d == %d+%d", i, EIG[i].D.sum(), mark[i].len,
       Ik.len, It.len);
   }

   if (nk+nt!=Etot.len) wblog(FL, 
   "ERR inconsistency in state count (%d+%d==%d)", nk,nt,Etot.len);


   if (!PA.isEmpty() && !PA.isIdentityPerm()) {
      if (!AK.isEmpty()) AK.Permute(PA);
      if (!AT.isEmpty()) AT.Permute(PA);
   }

   AK.Sort();
   AT.Sort();
   
   if (!AK.QIDX.isUniqueSorted()) 
   wblog(FL,"ERR AK does not have unique QIDX ???");
   
   if (!AT.QIDX.isUniqueSorted()) 
   wblog(FL,"ERR AT does not have unique QIDX ???");
}


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::map2Vec(const char *F, int L,
   wbvector<TD> &V, const wbvector<unsigned> &S
){
   if (!S.isEmpty()) {
      unsigned i,s, n=DATA.len, e=0, N=0;

      if (n!=S.len) ++e; else
      for (i=0; i<n && !e; ++i) {
         s=DATA[i]->SIZE.prod(); N+=s;
         if (S[i]!=s) ++e;
      }
      if (N!=V.len) ++e;

      if (e) wblog(F,L,
      "ERR %s() size changed (%d->%d; %d)",FCT,V.len,S.sum(),__LINE__);
   }
   else {
      unsigned N=getDataSize();
      if (V.len!=N) V.init(N);
   }

   return map2Vec(V.data);
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::map2VecI(const char *F, int L,
   wbvector<TD> &V, const char ref,
   const wbvector<unsigned> &S
){
   if (!S.isEmpty()) {
      unsigned i, n=DATA.len, e=0;
      if (n!=S.len) ++e;

      for (i=0; i<n && !e; ++i) if (DATA[i]->SIZE.prod()!=S[i]) ++e;
      if (e) wblog(F,L,
      "ERR %s() size changed (%d->%d; %d)",FCT,V.len,S.sum(),__LINE__);
   }
   else {
      unsigned N=getDataSize();
      if (V.len!=N) wblog(F,L,
      "ERR %s() size changed (%d->%d; %d)",FCT,V.len,N,__LINE__);
   }

   if (ref)
        return init2ref(V.data);
   else return map2Vec(V.data,'i');
}


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::map2Vec(TD *v, char Iflag
){
   unsigned i,s, n=DATA.len, N=0;
   
   for (i=0; i<n; ++i, v+=s) { s=DATA[i]->SIZE.prod();
      if (Iflag)
           memcpy(DATA[i]->data, v, s*sizeof(TD));
      else memcpy(v, DATA[i]->data, s*sizeof(TD));
      N+=s;
   }
   return N;
};


#endif

