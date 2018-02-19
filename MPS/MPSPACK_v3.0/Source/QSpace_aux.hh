#ifndef __WB_QSpace_AUX_HH__
#define __WB_QSpace_AUX_HH__

template <class TQ, class TD> class QSpace;
template <class TQ, class TD> class CPAT;

// ================================================================= //
// ensure ortogonal CGC spaces through decomposition/orthogonalization 
// w.r.t. existing CGC coefficient spaces; required for:
//   * QSpace::CGC contraction in the presence of outer multiplicity
//     => see contractCGS_ortho(), followed by contractDATA_ortho().
//   * make CGC spaces of given QSpace unique
//     => see makeUnique_ortho()
// Wb,Jul10,14
// ================================================================= //

template <class TQ>
class OrthoCGS {

  public:

    OrthoCGS() : cgr(NULL), conj(0), rtype(CGR_DEFAULT) {
    };


    OrthoCGS& init() {
       Q.init(); cgr=NULL; conj=0; rtype=CGR_DEFAULT;
       cgp.init(); cgQ.init(); cgR.init(); 
       return *this;
    };

    OrthoCGS& init(const CRef<TQ> &R) {
       Q.init(R); cgr=R.cgr; cgp=R.cgp; conj=R.conj; rtype=R.rtype;
       cgQ.init(); cgR.init();
       return *this;
    };

    template <class TA, class TB>
    unsigned contractCGS_ortho(
       const char *F, int L,
       const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
       const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
       unsigned isym
    );

    template <class TD>
    unsigned getCGS_ortho(
       const char *F, int L,
       QSpace<TQ,TD> &A, C_UVEC &Ia, unsigned isym
    );

    unsigned getOM() const { return getOM(wbvector<unsigned>()); }


   QSet<TQ> Q;
   const CDATA_TQ *cgr;

   wbperm cgp;
   char conj;

   cgrType rtype;

   wbarray<RTD> cgQ;
   wbarray<RTD> cgR;

  protected:
  private:
};

template <class TQ, class TA, class TB, class TC>
int contractDATA_ortho(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   const wbvector< OrthoCGS<TQ> > &Mc,
   QSpace<TQ,TC> &C, unsigned ic
);

template <class TQ, class TA, class TB, class TC>
void contractDATA_plain(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic
);

template <class TQ,class TD>
unsigned combineDATA_ortho(
   const char *F, int L,
   const QSpace<TQ,TD> &A, C_UVEC &Ia,
   const wbvector< OrthoCGS<TQ> > &MC, QSpace<TQ,TD> &C, unsigned ic
);


template<class TQ, class TD>
void mxInitQSpaceVec(
   const char *F, int L, const mxArray* C,
   wbvector< QSpace<TQ,TD> > &Fk,
   const char ref=0,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL
);

template<class TQ, class TD>
void mxInitQSpaceVecVec(
   const char *F, int L, const mxArray* C,
   wbvector< wbvector< QSpace<TQ,TD> > > &Fk, const char ref=0
);

template<class TQ, class TD>
void mxcInitQSpaceVec(
   const char *F, int L, const mxArray* C,
   wbvector< QSpace<TQ,TD> > &Fk, const char ref=0
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &Fk
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &F, wbindex &I
);

template<class TQ, class TD>
mxArray* QSpaceVecVec2Mx(
   const wbvector< wbvector< QSpace<TQ,TD> > > &F,
   const char ref=0
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mxc(
   const wbvector< QSpace<TQ,TD> > &F,
   const char ref=0
);

template<class TQ, class TD>
void putQSpaceVec(
   const wbvector< QSpace<TQ,TD> > &Fk,
   const char *vname, const char *ws="caller"
){
   mxArray *a=Fk.toMx();
   mxPutAndDestroy(FL,a,vname);
};

template <class TQ, class TD>
void getDiffOp(
   const QSpace<TQ,TD> &A, unsigned ia,
   const QSpace<TQ,TD> &B, unsigned ib, const QSpace<TQ,TD> &C
);

template<class TQ, class TD>
void getQall(const wbvector< QSpace<TQ,TD> > &F, wbMatrix<TQ> &QA);


template <class TQ, class TD>
class EIG_QS {
 public:

    EIG_QS() {};

    wbMatrix<TQ> Q;
    wbMatrix<INDEX_T> S;
    wbvector<INDEX_T> D;

    wbarray<TD> U;
    wbvector<TD> E;

    mxArray* toMx() const;
    void put(const char *vname, const char *ws="caller") {
       mxArray *a=toMx();
       mxPutAndDestroy(FL,a,vname,ws);
    };

 protected:
 private:
};

template <class TQ, class TD>
mxArray* EIG_QS<TQ,TD>::toMx() const {

   mxArray *a=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,a, "Q", Q.toMx());
   mxAddField2Scalar(FL,a, "S", S.toMx());
   mxAddField2Scalar(FL,a, "D", D.toMx());
   mxAddField2Scalar(FL,a, "U", U.toMx());
   mxAddField2Scalar(FL,a, "E", E.toMx());

   return a;
};


template <class TQ, class TD>
class blockSpaceQS {

  public:


     void init(unsigned d, unsigned qdim, const TQ *q, const QVec &qv);
     void init();

     void initWeights();
     void setWeight(unsigned i, unsigned j, double w);

     void addCell(
     unsigned i, unsigned j, const QSpace<TQ,TD> &H4, unsigned k);

     void cell2Full();
     void Eig();

     void initA(
        wbvector<unsigned> &I,
        QSpace<TQ,TD> &A, QSpace<TQ,TD> &H, TD *E0
     );

     mxArray* toMx() const;
     mxArray* mxCreateStruct(unsigned m, unsigned n) const;
     void add2MxStruct(mxArray *S, unsigned i) const;

     void put(const char *vname, const char *ws="base") const
     {  put(0,0,vname,ws); }

     void put(const char *file, int line,
        const char *vname, const char *ws="base"
     ) const {
        mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);

        if (i) wblog(__FL__,
        "ERR failed to write variable `%s' (%d)",vname,i);
        if (file) wblog(file,line,
        "I/O putting variable `%s' to `%s'",vname,ws);

        mxDestroyArray(a);
     };

     QVec qvec;
     wbMatrix< wbarray<TD> > C;

     wbMatrix<TQ> QQ;
     wbMatrix<unsigned> SS;
     wbvector<unsigned> D;
     wbvector<TQ> Q;

     wbarray<TD> M;
     wbarray<TD> W;
     wbarray<TD> U;
     wbvector<TD> E;

   protected:
   private:
};


template <class TQ, class TD>
mxArray* blockSpaceQS<TQ,TD>::toMx() const {
   mxArray *S=mxCreateStruct(1,1);
   this->add2MxStruct(S,0);
   return S;
}

template <class TQ, class TD>
mxArray* blockSpaceQS<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={
      "qvec","C","M","Q","size","dim","D","Q","weight","U","E"
   };
   return mxCreateStructMatrix(m,n,11,fields);
}

template <class TQ, class TD>
void blockSpaceQS<TQ,TD>::add2MxStruct(mxArray *s, unsigned i) const {

   mxSetFieldByNumber(s,i, 0,qvec.toMx());
   mxSetFieldByNumber(s,i, 1,C .toMx());
   mxSetFieldByNumber(s,i, 2,M .toMx());
   mxSetFieldByNumber(s,i, 3,QQ.toMx());
   mxSetFieldByNumber(s,i, 4,SS.toMx());
   mxSetFieldByNumber(s,i, 5,D .toMx());
   mxSetFieldByNumber(s,i, 6,numtoMx(D.sum()));
   mxSetFieldByNumber(s,i, 7,Q .toMx());
   mxSetFieldByNumber(s,i, 8,W .toMx());
   mxSetFieldByNumber(s,i, 9,U .toMx());
   mxSetFieldByNumber(s,i,10,E .toMx());
}


template <class TQ, class TD> inline
void blockSpaceQS<TQ,TD>::init() {
   D.init(); SS.init(); W.init(); U.init(); E.init();
   C.init(); QQ.init(); Q.init(); M.init();
};

template <class TQ, class TD> inline
void blockSpaceQS<TQ,TD>::init(
   unsigned d, unsigned qdim, const TQ *q, const QVec &qv
){
   D.init(d);
   SS.init(d,2);
   QQ.init(d,qdim); Q.init(qdim/2,q); qvec=qv;
   W.init(d,d);
   C.initDef(d,d); M.init(); U.init(); E.init();
};

template <class TQ, class TD> inline
void blockSpaceQS<TQ,TD>::initWeights() {
   unsigned d=D.len;
   W.init(d,d);
};

template <class TQ, class TD> inline
void blockSpaceQS<TQ,TD>::setWeight(unsigned i, unsigned j, double w) {
   if (W.SIZE.len!=2) wblog(FL,"ERR invalid rank-%d",W.SIZE.len);
   if (i>=W.SIZE[0] || j>=W.SIZE[1]) wblog(FL,
   "ERR index out of bounds (%d,%d) %s",i,j,W.sizeStr().data); 
   if (W(i,j)!=0 && W(i,j)!=w) wblog(FL,
   "ERR W(%d,%d) already set (%.3g, %.3g)",i,j,W(i,j),w);

   W(i,j)=w;
};


template <class TQ, class TD>
void blockSpaceQS<TQ,TD>::addCell(
   unsigned i, unsigned j, const QSpace<TQ,TD> &H4, unsigned k
){
   if (k>=H4.DATA.len) wblog(FL,"ERR index out of bounds");

   wbarray<TD> &B=*(H4.DATA[k]);
   unsigned d1=0,d2=0, *s=B.SIZE.data,*s2=NULL, d=QQ.dim2;
   unsigned r=B.SIZE.len, r2=r/2;
   TQ *q=H4.QIDX.rec(k), *q2=NULL;

   if (B.isEmpty()) { wblog(FL,"WRN got empty block."); return; }

   if (i>=C.dim1 || j>=C.dim2 || i>=QQ.dim1 || j>=QQ.dim1) wblog(FL,
      "ERR index out of bounds (%d,%d/%d %d/%d)",
       i,C.dim1,QQ.dim1,j,C.dim2);
   if (H4.QIDX.dim2!=2*d || SS.dim2!=r2) wblog(FL,
      "ERR Q/Size mismatch (%d/2*%d; %d/%d)",
       H4.QIDX.dim2,QQ.dim2,SS.dim2,r2);

   if (r==4) {
      B.SIZE.prod2(d1,d2); s2=s+r2; q2=q+d;
   }
   else if (r==2) {
      d1=d2=B.SIZE.prod(); s2=s; q2=q;
   }
   else wblog(FL,"ERR invalid rank-%d object",r);

   if (D[i]) {
      if (D[i]!=d1 || Wb::cmpRange(s,SS.ref(i),r2)) wblog(FL,
         "ERR %s() size mismatch (%d i: %d,%d):\n%s <> %dx%d (%s)",
        __FUNCTION__,k+1,i+1,j+1,B.sizeStr().data,d1,d2,
          SS.rec2Str(i,"","x").data);
      if (Wb::cmpRange(q,QQ.ref(i),d)) wblog(FL,
         "ERR %s() Q mismatch (%d i: %d,%d):\n[%s] <> [%s]",
        __FUNCTION__,k+1,i+1,j+1,H4.QIDX.rec2Str(k).data,
          QQ.rec2Str(i).data
      );
   }
   else { D[i]=d1;
      cpyRange(s,SS.ref(i),r2);
      cpyRange(q,QQ.ref(i),d);
   }

   if (D[j]) {
      if (D[j]!=d2 || Wb::cmpRange(s2,SS.ref(j),r2)) wblog(FL,
         "ERR %s() size mismatch (%d j: %d,%d):\n%s <> %dx%d (%s)",
        __FUNCTION__,k+1,i+1,j+1,B.sizeStr().data,d1,d2,
          SS.rec2Str(j,"","x").data);
      if (Wb::cmpRange(q2,QQ.ref(j),d)) wblog(FL,
         "ERR %s() Q mismatch (%d j: %d,%d):\n[%s] <> [%s]",
        __FUNCTION__,k+1,i+1,j+1,H4.QIDX.rec2Str(k).data, 
          QQ.rec2Str(j).data
      );
   }
   else { D[j]=d2;
      cpyRange(s2,SS.ref(j),r2);
      cpyRange(q2,QQ.ref(j),d);
   }

   C(i,j).Plus(B,W(i,j),'i');
};


template <class TQ, class TD>
void blockSpaceQS<TQ,TD>::cell2Full() {

   unsigned i,j,r, d1=C.dim1, d2=C.dim2, n=d1*d2;
   wbvector<unsigned> cD,I0(2);
   wbarray<TD> B;

   if (n==0 || D.sum()==0) { init(); return; }

   D.cumsum0(cD); n=D.sum(); M.init(n,n);

   for (i=0; i<d1; i++)
   for (j=0; j<d2; j++) { r=C(i,j).SIZE.len;
      if (r==0) continue;
      if (r!=2 && r!=4) wblog(FL,"ERR invalid rank (r=%d)",r);

      I0[0]=cD[i]; I0[1]=cD[j];

      C(i,j).toMatrixRef(B,2);
      M.addBlock(I0, B, r==2 ? 'd' : 0);
   }

};


template <class TQ, class TD>
void blockSpaceQS<TQ,TD>::initA(
   wbvector<unsigned> &I,
   QSpace<TQ,TD> &A, QSpace<TQ,TD> &H, TD *E0
){
   unsigned i,i1,m,n=QQ.dim1, d=Q.len, d2=2*d, s=d*sizeof(TQ), s2=2*s;
   wbvector<unsigned> S(3);
   wbvector<TD> ee;
   wbarray<TD> X;
   TQ *q;

   if (U.SIZE.len!=2 || U.SIZE[0]!=D.sum()) wblog(FL,
      "ERR invalid eigenvector space (%s/%d)",U.sizeStr().data,D.sum());
   if (n!=D.len || QQ.dim2!=d2) wblog(FL,
      "ERR internal size mismatch (%d/%d,%d/%d)",
       SS.dim1,QQ.dim1,D.len,QQ.dim2,2*d);
   if (SS.dim1!=n || SS.dim2!=2) wblog(FL,
      "ERR invalid size specs %dx%d (%d/%d)",SS.dim1,SS.dim2,n,2);

   U.select0(I,1,X);
   A.init(n,3,d); q=A.QIDX.data; m=A.QIDX.dim2;

   for (i1=i=0; i<n; i++, q+=m) { d=D[i];
      memcpy(q,QQ.rec(i),s2);
      memcpy(q+d2,Q.data,s );
      X.select0(i1,i1+d-1,1,*A.DATA[i]); i1+=d;

      S[0]=SS(i,0); S[1]=SS(i,1); S[2]=I.len;
      A.DATA[i]->Reshape(S);
   }
   A.qvec=qvec;

   E.select(I,ee); if (E0) ee-=(*E0); d=Q.len;

   H.init(1,2,d); q=H.QIDX.data;
   H.DATA[0]->init2Vec(ee);

   memcpy(q  , Q.data,s);
   memcpy(q+d, Q.data,s);

   H.qop.init(qvec.len);
   H.qvec=qvec;
};


template <class TQ, class TD>
void blockSpaceQS<TQ,TD>::Eig() {

   if (!M.isEmpty()) {
      if (!M.isHConj()) wblog(FL,"WRN H4 block not hconj !?");
      wbEigenS(M,U,E);
   }
   else { U.init(); E.init(); }
};


int mxIsQSpace(
   const char *F, int L,
   const mxArray *a, unsigned &rank,
   char cflag=0,
   unsigned k=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   const char *istr=NULL
);

inline int mxIsQSpace(
   const char *F, int L, const mxArray *a, char cflag=0,
   const char *istr=NULL
){
   unsigned r=-1; return mxIsQSpace(F,L,a,r,cflag,-1,NULL,NULL,istr);
};

inline int mxIsQSpace(
   const mxArray *a, char cflag=0, const char *istr=NULL
){
   unsigned r=-1; return mxIsQSpace(NULL,0,a,r,cflag,-1,NULL,NULL,istr);
};

bool mxIsQSpaceArr(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   int arrdim=2,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0
);

bool mxIsQSpaceVecVec(const mxArray *C, unsigned r=-1, char cflag=0);

inline bool mxIsQSpaceVec(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,1,rmin,rmax,cflag); };

inline bool mxIsQSpaceVecR23(
   const char *F, int L,
   const mxArray *S,
   char cflag=0
){ 
   unsigned rank=-1, rmin=2, rmax=3;
   return mxIsQSpaceArr(F,L,S,rank,1,&rmin,&rmax,cflag);
};

inline bool mxIsQSpaceMat(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,2,rmin,rmax,cflag); };

inline bool mxIsEmptyQSpace(const mxArray *a, unsigned k=0);

inline
bool mxIsScalarQSpace(const mxArray *a, unsigned k=0);

inline bool mxIsQSpaceScalar(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,0,rmin,rmax,cflag); };

bool mxsIsQSpaceVec(const mxArray *S, int fid, unsigned r=-1, char cflag=0);
bool mxsIsQSpaceVEC(const mxArray *S, int fid, unsigned r=-1, char cflag=0);
   
int mxIsQSpaceVec(
   const char *F, int L,
   const char *fname, const char *vname,
   unsigned rank=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0,
   unsigned N=-1
);

int mxIsQSpaceVEC(
   const char *F, int L,
   const char *fname, const char *vname,
   unsigned rank=-1,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL,
   char cflag=0
);

bool mxIsQSpaceVecOrEmpty(const mxArray *S, int rank, char cflag=0) {
   if (mxIsEmpty(S)) return 1;
   return mxIsQSpaceVec(FL,S,rank,NULL,NULL,cflag);
};


inline int mxIsQSpaceOrEmpty(
   const mxArray *a, unsigned rank=-1, char cflag=0, unsigned k=-1
){ if (mxIsEmpty(a)) return 1;
   else return mxIsQSpace(NULL,0,a,rank,cflag,k);
};

inline int mxIsQSpaceOrEmpty(const char *F, int L,
   const mxArray *a, unsigned rank=-1, char cflag=0, unsigned k=-1
){ if (mxIsEmpty(a)) return 1;
   else return mxIsQSpace(F,L,a,rank,cflag,k);
};


#endif

