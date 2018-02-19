#ifndef __WB_QSPACE_COL_MAJOR_HH__
#define __WB_QSPACE_COL_MAJOR_HH__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //
// QSpace - set of objects (possibly complex matrices)
// organized by generalized index set storeed as row vectors (recs)
// in index matrix.
// Wb,Aug08,05  Apr07,06
//
// Introduction of non-abelian QSpace through variable 'qtype'.
// see also clebsch.hh
// Wb,Sep25,09
// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

#ifdef WB_CLOCK
#define WBC_QSPACE_IO
#endif

enum QS_TYPES {
   QS_NONE, QS_OPERATOR, QS_AMATRIX,
   QS_NUM_TYPES
};

const char* QS_STR[QS_NUM_TYPES] = {
  "", "operator", "A-matrix",
};


template <class TQ, class TD>
class QSpace {

  public:
    
    QSpace(unsigned n=0) : QDIM(n), otype(QS_NONE), isref(0) {};

    QSpace(unsigned m, unsigned n, unsigned nq)
     : QDIM(nq), otype(QS_NONE), isref(0) { init(m,n,nq); };

    QSpace(const QSpace &B)
     : QIDX(B.QIDX), CGR(B.CGR),
       qtype(B.qtype), QDIM(B.QDIM),
       otype(QS_NONE), itags(B.itags), isref(0)
    {
       setupDATA();
       for (unsigned i=0; i<B.DATA.len; ++i) {
           DATA[i]->init(*(B.DATA[i]));
       }
    };

    QSpace(const char *F, int L, const mxArray *a, char ref=0) 
     : otype(QS_NONE), isref(0) { init(F,L,a,ref); };
    QSpace(const mxArray *a, char ref=0)
     : otype(QS_NONE), isref(0) { init(FL,a,ref); };

   ~QSpace() { clearQSpace(); };

    QSpace& clearQSpace() {
        QIDX.init(); qtype.init(); itags.init();
        clearDATA(); CGR.init(); isref=0;
        return *this;
    };

    void clearDATA() {
        if (DATA.len) {
           if (!isref && !DATA.isref) {
              for (unsigned i=0; i<DATA.len; ++i) WB_DELETE_1(DATA[i]);
           }
           DATA.init();
        }

    };

    unsigned rank(const char *F=NULL, int L=0) const {
       if (!QDIM) {
          if (QIDX.dim2) wblog(F_L,"ERR QSpace inconsistency "
             "(QIDX: %dx%d/%d)",QIDX.dim1, QIDX.dim2, QDIM);
          return 0;
       }

       if (QIDX.dim2%QDIM) wblog(F_L,"ERR QSpace inconsistency "
          "(QIDX: %dx%d/%d)",QIDX.dim1, QIDX.dim2, QDIM);

       return (QIDX.dim2/QDIM);
    };

    unsigned len() const {
       if (QIDX.dim1!=DATA.len) wblog(FL, 
          "ERR severe qspace inconsistency (QIDX.dim1=%d, len=%d).",
           QIDX.dim1, DATA.len);
       return QIDX.dim1;
    };

    TQ* Qref(unsigned r, unsigned k) { unsigned l=k*QDIM;
        if (!QDIM || r>=QIDX.dim1 || (l+QDIM)>QIDX.dim2) wblog(FL,
           "ERR %s() index out of bounds (%d/%d; %d*%d/%d)",
           FCT,r,QIDX.dim1,k,QDIM,QIDX.dim2
        );
        return (QIDX.data + r*QIDX.dim2 + l);
    };

    const TQ* Qref(unsigned r, unsigned k) const {
        unsigned l=k*QDIM;
        if (!QDIM || r>=QIDX.dim1 || (l+QDIM)>QIDX.dim2) wblog(FL,
           "ERR %s() index out of bounds (%d/%d; %d*%d/%d)",
           FCT,r,QIDX.dim1,k,QDIM,QIDX.dim2
        );
        return (QIDX.data + r*QIDX.dim2 + l);
    };


    void init() { clearQSpace(); };

    void init(unsigned m, unsigned n, unsigned nq) {
       clearQSpace();
       QDIM=nq; QIDX.init(m,nq*n);
       setupDATA();
    };

    QSpace& init(const QSpace &B) {
        clearQSpace();
        QIDX=B.QIDX; QDIM=B.QDIM; qtype=B.qtype; otype=B.otype;
        setupDATA(); itags=B.itags; CGR=B.CGR;

        if (qtype.len!=CGR.dim2) wblog(FL,
           "ERR %s() got CGR size mismatch (%dx%d <> %dx%d)",
           FCT,CGR.dim1,CGR.dim2,QIDX.dim1,qtype.len
        );

        for (unsigned i=0; i<B.DATA.len; ++i)
        DATA[i]->init(*(B.DATA[i]));

        return *this;
    };

    void init(const char *file, int line,
    const mxArray *S, const char ref=0, const unsigned k=0);

    void init(const wbMatrix<TQ> &Q, const unsigned D,
       const QVec *qv=NULL
    ){
       clearQSpace();
       QIDX=Q; QDIM=D; if (qv) qtype=(*qv); else qtype.init();
       setupDATA(FL);
    };

    QSpace& setupDATA (const char *file=NULL, int line=0) {
        if (!QIDX.dim1) { clearDATA(); return *this; }

        if (QIDX.dim2 && (!QDIM || QIDX.dim2%QDIM)) wblog(__FL__,
           "ERR invalid QIDX (%dx%d; %d)",QIDX.dim1,QIDX.dim2,QDIM);
        clearDATA();

        if (isref>1) wblog(FL,"ERR %s() got isref=%d",FCT,isref);

        DATA.init(QIDX.dim1);

        if (!isref) {
           for (unsigned i=0; i<DATA.len; ++i) {
           WB_NEW_1(DATA[i],wbarray<TD>); }
        }

        return *this;
    };

    char gotCGS(const char *F=NULL, int L=0) const;

    char gotCGX(const char *F=NULL, int L=0) const;

    QSpace& Enlarge(unsigned m) {
       if (!isConsistent(FL) || !QIDX.dim1) wblog(FL,
          "ERR %s() can't enlarge empty QSpace (%d)",FCT,QIDX.dim1);
       if (!m) return *this;

       unsigned i,j, n0=QIDX.dim1, n=n0+m;
       gotCGS(FL);

       QIDX.Resize(n, QIDX.dim2);
       DATA.Resize(n);
       for (i=n0; i<n; ++i) { WB_NEW_1(DATA[i],wbarray<TD>); }

       CGR.Resize(n,CGR.dim2);
       return *this;
    };

    QSpace& Trim(unsigned n) {
       if (!isConsistent(FL,0) || n>QIDX.dim1) wblog(FL,
       "ERR %s() can't trim to LARGER QSpace (%d/%d)",FCT,n,QIDX.dim1);
       if (n==QIDX.dim1) { return *this; }

       unsigned i,j, n2=QIDX.dim1;
       gotCGS(FL);

       QIDX.Resize(n, QIDX.dim2);

       for (i=n; i<n2; ++i) { WB_DELETE_1(DATA[i]); }
       DATA.Resize(n);

       CGR.Resize(n,CGR.dim2);

       return *this;
    };

    void init2scalar(const TD &x) {
        clearQSpace();
        DATA.init(1);

        WB_NEW_1(DATA[0],wbarray<TD>(1,1));
        DATA[0]->data[0]=x;
    };

    void init2ref(const mxArray *S){
       init(FL,S,'r',0);
    };

    void init2ref(const QSpace &);
    unsigned init2ref(const TD *v);

    void initOpZ_WET(const char *F, int L, const QVec &qvec,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2,
       const wbMatrix<TQ> &Q, const wbarray<TD> &D3,
       const double eps1=1E-10, const double eps2=1E-14
    );

    void reduceMatEl(const char *F, int L,
       const wbMatrix<TQ> &Q1, const wbvector<TD> &dd,
       const wbMatrix<TQ> &Q2, const wbvector<TD> &dc,
       wbvector<unsigned> &I1, wbvector<unsigned> &I2,
       wbvector<TD> &dr
    ) const;

    void unRef() {

        if (!isref && !qtype.isref && !itags.isref && 
            !QIDX.isref && !DATA.isref && !CGR.isref
        ) return;

        wblog(FL,"WRN %s() [%d%d%d%d%d%d]",FCT,
        isref, QIDX.isref, DATA.isref, CGR.isref, qtype.isref, itags.isref);

        QIDX.unRef(); qtype.unRef(); itags.unRef();

        if (DATA.isref) DATA.unRef(); else
        for (unsigned i=0; i<DATA.len; ++i) DATA[i]->unRef();
        
        if (CGR.isref) CGR.unRef();
    };

    void initIdentity(const QSpace &A, char dflag=0);
    void initIdentity(const wbMatrix<TQ> &Q);

    void initIdentity(
    const QSpace **FF, unsigned n, char dflag=0);

    void initIdentity(
       const wbvector< QSpace > &FF, char dflag=0
    ){
       wbvector< const QSpace* > ff(FF.len);
       for (unsigned i=0; i<FF.len; ++i) ff[i]=&FF[i];
       initIdentity(ff.data, ff.len, dflag);
    };

    void initIdentity(
       const QSpace &H,
       const wbvector< QSpace > &FF, char dflag=0
    ){
       wbvector< const QSpace* > ff(FF.len+1);
       ff[0]=&H; for (unsigned i=0; i<FF.len; ++i) ff[i+1]=&FF[i];
       initIdentity(ff.data, ff.len, dflag);
    };

    QSpace& setupCGR(const char *F=NULL, int L=0);
    QSpace& initIdentityCGS(const char *F=NULL, int L=0);
    QSpace& initIdentityCGS(const char *F, int L, const QSpace<TQ,TD>& H);

    void init_gRG_zdim_bare(const char *F, int L) const;

    template<class TDA>
    QSpace<TQ,TD>& initIdentityCG(
       const wbvector< QSpace<TQ,TDA> > &A, const wbindex &ia,
       char zflag=0
    );

    template<class TDA, class TDB>
    QSpace<TQ,TD>& initIdentityCG(
       const wbvector< QSpace<TQ,TDA> > &A, const wbindex &ia,
       const wbvector< QSpace<TQ,TDB> > &B, const wbindex &ib,
       char vflag='v'
    );

    template<class TDA, class TDB>
    QSpace<TQ,TD>& initIdentityCG(
       const wbvector< QSpace<TQ,TDA> > &A,
       const wbvector< QSpace<TQ,TDB> > &B,
       const char *pstr=NULL
    ){ wbindex i;
       initIdentityCG(A,i,B,i);
       if (pstr && pstr[0]) Permute(wbperm(pstr));
       return *this;
    }

    QSpace& initDiagonal(
        const wbMatrix<TQ> &Q,
        const wbvector<unsigned> &S,
        const wbvector<TD> &D,
        char rcflag
    );

    QSpace& init2DiffOp(
       const QSpace &A, unsigned ia,
       const QSpace &B, unsigned ib
    );

    template <class TDA>
    QSpace& initQ(const QSpace<TQ,TDA> &A, char full=1){
       if (isref) wblog(FL,"WRN %s() got QSpace ref!",FCT);
       qtype=A.qtype; QDIM=A.QDIM; if (full) {
       otype=A.otype; itags=A.itags; }
       return *this;
    };

    template <class TDA, class TDB>
    QSpace& initQ(const char *F, int L,
       const QSpace<TQ,TDA> &A, char full=1, const QSpace<TQ,TDB> *B=NULL
    ){
       if (B) A.checkQ(F,L,*B);
       return initQ(A,full);
    };

    void initQ(unsigned n, unsigned r, const QSpace &A){
       clearQSpace();
       QDIM=A.QDIM; qtype=A.qtype; QIDX.init(n,r*QDIM);
       setupDATA(FL);
       if (qtype.len) setupCGR(FL);
    };

    void init_itags(const char *F, int L, const char *s, unsigned k);

    QSpace& SetTag(unsigned k, const iTag &b, char lflag=0);

    QSpace& SetConjTags() {
       for (unsigned k=0; k<itags.len; ++k) { itags.data[k].Conj(); }
       return *this;
    };

    QSpace& SetConjTags(unsigned k1, unsigned k2) {
       if (!k1 || !k2) wblog(FL,
          "ERR %s(%d) expects 1-based index",FCT,k1,k2);
       itags.SetConj(0,0,k1-1);
       itags.SetConj(0,0,k2-1); return *this;
    };

    QSpace& SetConjTag(unsigned k) {
       if (!k) wblog(FL,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetConj(0,0,k-1); return *this;
    };

    QSpace& SetConjTag(const char* F, int L, unsigned k) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetConj(F,L,k-1); return *this;
    };

    QSpace& SetConj(const char* F, int L, unsigned k) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetConj(F,L,--k);
       for (unsigned n=CGR.numel, i=0; i<n; ++i) {
           CGR[i].Conj(k);
       }
       return *this;
    };

    QSpace& SetFlag(const char *F, int L, unsigned k, unsigned l) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetFlag(F_L,k-1,l); return *this;
    };

    QSpace& UnsetFlag(const char *F, int L, unsigned k, unsigned l) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.UnsetFlag(F_L,k-1,l); return *this;
    };

    QSpace& UnsetFlags() { itags.UnsetFlags(); return *this; };

    template<class TDB>
    unsigned matchITags(
       const char *F, int L, const QSpace<TQ,TDB> &B,
       ctrIdx &ia, ctrIdx &ib) const;

    unsigned RemoveZLabels(
       const char *F, int L, wbMatrix<TQ> &QQ, wbMatrix<TQ> *Z=NULL) const;

    void RemoveZLabels(const char *F, int L, wbMatrix<TQ> *Z){
       if (QDIM==0 && qtype.len) { QDIM=qtype.Qlenz(); }
       QDIM=RemoveZLabels(F_L,QIDX,Z);
    };

    void Expand2Projector(const TD eps=0);
    void ExpandDiagonal(unsigned i1=1, unsigned i2=2);
    void Diag2Vec(char tflag=0);

    void ExpandQ(
       const wbMatrix<TQ> &QK, const wbvector<unsigned> &d,
       unsigned k, double ra, wbindex *I=NULL,
       wbvector<char> *mQ=NULL
    );

    unsigned long getDataSize(wbMatrix<unsigned> &S) const;
    unsigned long getDataSize(wbvector<unsigned> &S) const;
    unsigned long getDataSize() const;
    unsigned long getCGSSize() const;
    wbstring totSize2Str() const;

    unsigned long getSize() const;

    unsigned map2Vec (const char *F, int L, wbvector<TD> &V,
    const wbvector<unsigned> &S=wbvector<unsigned>());

    unsigned map2VecI(const char *F, int L, wbvector<TD> &V,
    const char ref=0, const wbvector<unsigned> &S=wbvector<unsigned>());

    unsigned map2Vec(TD *v, char Iflag=0);

    void ExpandQ(
    CPAT<TQ,TD> &CP, double nrm, wbvector<unsigned> &xflag, char disp=0);
    
    void ResetRec(unsigned k) {
       isConsistent(FL); if (k>=QIDX.dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,QIDX.dim1);

       QIDX.setRec(k,0); DATA[k]->init();
       for (unsigned j=0; j<CGR.dim2; ++j) CGR(k,j).init();
    };

    void setDATA(const TD c) {
        for (unsigned i=0; i<DATA.len; ++i)
        DATA[i]->set(c);
    };

    void setRand(TD nrm=0);

    void markUnequalZero(const TD m=1, const TD eps=0) {
        for (unsigned i=0; i<DATA.len; ++i)
        DATA[i]->markUnequalZero(m,eps);
    };

    void PrependSingletons(unsigned R);

    void swap(QSpace &B) {
        wblog(FL,"TST check this");
        if (isref!=B.isref) wblog(FL,
        "ERR cannot swap data with isref=%d,%d", isref, B.isref);
        SWAP(QDIM,B.QDIM);
        SWAP(QIDX,B.QIDX); SWAP(qtype,B.qtype); SWAP(CGR,B.CGR);
        SWAP(DATA,B.DATA);

        SWAP(otype,B.otype);
        SWAP(itags,B.itags);
    };

    QSpace& save2(QSpace &B) {
        if (isref || QIDX.isref || DATA.isref || CGR.isref) wblog(FL,
           "ERR must not save referenced data (isref=%d,%d,%d)",
           isref, QIDX.isref,DATA.isref);
        for (unsigned i=0; i<DATA.len; ++i) {
           if (DATA[i]->isRef()) wblog(FL,"%NERR using save2 with "
              "referenced data\n(data[%d] has isref=%d)",
              i,DATA[i]->isRef()
           );
        }

        B.clearQSpace();

        DATA.save2(B.DATA);
        QIDX.save2(B.QIDX); B.QDIM=QDIM;
        CGR.save2(B.CGR);

        itags.save2(B.itags);
        qtype.save2(B.qtype);
        B.otype=otype; otype=QS_NONE;

        B.isref=isref; isref=0;

        return B;
    };


    char isAbelian() const {
       return qtype.allAbelian();
    };

    bool hasMP() const {
       if (CGR.isEmpty() || qtype.allAbelian()) return 0;
       return qtype.hasMP(rank(FL));
    };


    bool isEmpty() const {
       if (QIDX.dim1!=DATA.len && (QIDX.dim1 || DATA.len>1)) {
          wblog(FL,"ERR severe QSpace inconsistency (%d,%d)",
          QIDX.dim1, DATA.len);
       }
       return (QIDX.data==NULL);
    };

    bool isConsistent(
       int &r, const char *F=NULL, int L=0, char level=2) const;

    bool isConsistent(
      const char *F=NULL, int L=0, char level=2
    ) const { int r=-1; return isConsistent(r,F,L,level); };

    bool isConsistentR(
      int r, const char *F=NULL, int L=0, char level=2
    ) const { return isConsistent(r,F,L,level); };

    bool isScalar() const {
        if (!isConsistentR(0)) { info("this"); wberror(FL,str); }

        wbvector<unsigned> &S = DATA[0]->SIZE;

        if (QIDX.dim1!=1 || S[0]!=1 || S[1]!=1) wblog(FL,
        "ERR severe QSpace inconsistency (scalar %d,%d)",
        QIDX.dim1, S[0], S[1]);

        return 1;
    };

    bool isComplex() const {
        isConsistent(FL);

        for (unsigned n=DATA.len, i=0; i<n; ++i)
        if (DATA[i]->isComplex()) return 1;

        return 0;
    };

    bool isDiagBlock(unsigned k) const;

    bool isBlockDiagMatrix(const char *F, int L, char dflag=0) const;
    bool isBlockDiagMatrix(char dflag=0) const {
         return isBlockDiagMatrix(NULL,0,dflag); }

    bool isDiagMatrix(const TD eps=1e-14) const;
    bool isIdentityMatrix(const TD eps=1e-14) const;
    char hasIdentityCGS(const RTD eps=1e-14) const;

    bool isR3Op() const {
       return (rank(FL)==3 && otype==QS_OPERATOR);
    };

    int isOperator(unsigned *dr, char xflag=0) const;

    bool isOperator() const { return (isOperator(NULL,0)>0 ? 1 : 0); };

    bool isQSym(const char *F=0, int L=0, char dflag=0) const;

    bool isHConj(
      const char *F=0, int L=0, double eps=1e-14, char vflag=0
    ) const;

    bool isAHerm(
      const char *F=0, int L=0, double eps=1e-14, char vflag=0
    ) const;

    int isSingleton(unsigned k) const;

    template <class T2>
    bool hasSameQ(const QSpace<TQ,T2> &B) const;

    bool hasSameQ(
       unsigned r1, unsigned k1,
       unsigned r2, unsigned k2
    ) const {
       if (memcmp( Qref(r1,k1), Qref(r2,k2), QDIM*sizeof(TQ) )) return 0;
       if (DATA[r1]->dim0(k1) != DATA[r2]->dim0(k2)) {
          wblog(FL,
             "WRN %s()=1, yet DATA size mismatch (%s <> %s; %d,%d)",
              FCT, DATA[r1]->sizeStr().data,
                   DATA[r2]->sizeStr().data, k1+1, k2+1); 
          return 0;
       }
       return 1;
    };

    template <class T2>
    bool hasSameRank(const QSpace<TQ,T2> &B) const {
        return (QIDX.dim2==B.QIDX.dim2 && QDIM==B.QDIM);
    };

    template <class TDB>
    void checkQ(
       const char *F, int L, const QSpace<TQ,TDB> &B,
       char aflag=0
     ) const;

    template <class TDB>
    bool sameType(
       const QSpace<TQ,TDB> &B, int r=-1, const char *istr=NULL
    ) const;


    bool sameQDir(const char *F, int L, const QSpace &B) const;

    bool hasQOverlap(const QSpace &A) const;

    template <class TDB>
    int hasQOverlap(
       unsigned ica, const QSpace<TQ,TDB> &B, unsigned icb,
       const char type=0
    ) const;

    int getQOverlap(
       const QSpace &B,
       wbMatrix<TQ> &QQ,
       wbMatrix<int> &IQ,
       char olonly=0
    ) const;

    bool findDimQ(const TQ* q, unsigned &d) const;
    bool findDimQ(const TQ* q, unsigned k, unsigned &d) const;


    bool operator==(const QSpace &B) const {
       if (this!=&B) {
          if (!(QIDX==B.QIDX)) return 0;
          if (!(qtype==B.qtype)) return 0;
          if (!(CGR==B.CGR)) return 0;
          if (!DATA.deepEqualP(B.DATA)) return 0;
       }
       return 1;
    };

    bool operator!=(const QSpace &B) const { return (!((*this)==B)); }

    template<class T>
    void operator*= (const T fac);

    template<class T>
    QSpace& times(const T fac, QSpace &B) const {
       B=*this; B*=fac; return B;
    };

    QSpace& operator= (const QSpace &A) {
      return init(A);
    };

    void operator+= (const QSpace &Q);
    void operator-= (const QSpace &Q);

    QSpace& Cat(const char *F, int L,
       const QSpace<TQ,TD> &A, const QSpace<TQ,TD> &B,
       TD afac=1, TD bfac=1, char uflag=1
    );

    QSpace& plus(
       const QSpace &B, QSpace &C, TD bfac=1, char vflag=0) const;

    QSpace& minus(
       const QSpace &B, QSpace &C, TD bfac=1, char vflag=0) const;

    QSpace& Plus (const QSpace &A, TD bfac=1);
    QSpace& Minus(const QSpace &A, TD bfac=1);

    QSpace& plus_plain(
       const QSpace &A, QSpace &B, TD bfac=1, char vflag=0) const;


    QSpace& NormCG(char skipzeros=1);

    TD norm2() const;

    TD norm() const { return sqrt(norm2()); }

    TD normDiff2(const QSpace &B) const;
    TD scalarProd(const QSpace &B) const;

    void checkNorm(
       const char *F, int L, double nrm, char vflag=1, double eps=1E-12
     ) const;

    void TimesEl(const QSpace &B, const char qmatch=0);

    TD sumData() const;
    TD trace() const;

    QSpace& TensorProd(const QSpace &B, QSpace &C,
    const char aflag='N', const char bflag='N') const;

    void TensorProdUnity(const wbMatrix<TQ> &Q, QSpace &C) const;

    void EigenSymmetric(wbvector<TD> &Etot) const {
       QSpace Ak,At;
       QSpace<TQ,double> Ek,Et;
       wbMatrix<unsigned> DB; int Nkeep=-1;

       EigenSymmetric(Ak,At,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric(
       QSpace &Ak, QSpace &Ek,
       wbvector<TD> &Etot
    ) const {
       QSpace At;
       QSpace<TQ,double> Et;
       wbMatrix<unsigned> DB; int Nkeep=-1;

       EigenSymmetric(Ak,At,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric(QSpace &Ak, QSpace &Ek) const {
       QSpace At;
       QSpace Et;
       wbMatrix<unsigned> DB; int Nkeep=-1;
       wbvector<TD> Etot;

       EigenSymmetric(Ak,At,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric(
       QSpace     &Ak, QSpace     &At,
       QSpace<TQ,double> &Ek, QSpace<TQ,double> &Et, wbvector<double> &Etot,
       wbMatrix<unsigned> &DD, int &Nkeep,
       double Etrunc=0, double *E0=NULL,
       const wbperm &P=wbperm(), double eps=0., double b=0., int dmax=-1,
       const char *sort=NULL
    ) const;

    void Eigen_CSymmetric(
       QSpace &Ak, QSpace &At,
       QSpace &Ek, QSpace &Et, wbvector<TD> &Etot,
       wbvector<unsigned> &DD, int &Nkeep,
       double Etrunc=0, TD *E0=NULL,
       const wbperm &P=wbperm(), double eps=0., double b=0., int dmax=-1,
       const char *sort=NULL
    ) const;

    void GroupIndizes(unsigned K);

    void Append(unsigned k);
    void Append(const wbvector<TQ> &, const wbarray<TD> &D);

    QSpace& Append(const char *F, int L, const QSpace &A, char uflag=1);

    int Append2AndDestroy(
    const char *F, int L, QSpace &A, char unique=0);

    void TransferSpace(const char *F, int L,
       QSpace &B, unsigned k,
       const QSpace<TQ,double> &W, double eps
    );

    void TransferSpace(const char *F, int L,
       QSpace &B, unsigned k
    );

    unsigned makeUnique(char omFlag=-1);
    unsigned makeUnique_plain(char omFlag=-1);
    unsigned makeUnique_ortho(char vflag=0);

    bool isZero(const double eps=0., char flag=0) const;
    bool hasZeroData (const double eps=0., char flag=0) const;
    unsigned SkipZeroData(const double eps=0., char flag=0);

    unsigned skipZeroOffDiag(double eps=0., char bflag=0);

    unsigned SkipEmptyData(const char *F=0, int L=0);

    bool gotZeroData(unsigned i, const double eps=0.) const;

    void Set(unsigned k, const wbvector<TQ> &Q, const wbarray<TD> &D);
    void SetQ(unsigned i, unsigned k, const TQ *d);

    QSpace& Select(const wbindex &I, char data_only=0);

    template <class TM>
    QSpace& SelectMarked(const wbvector<TM> &mark){
       if (mark.len!=QIDX.dim1) wblog(FL,
          "ERR %s() size inconsistency (%d/%d)",FCT,mark.le,QIDX.dim1);
       wbindex I; mark.find(I);
       return Select(I);
    };

    void Sort(const QSpace &A);
    void Sort() {
       wbperm P; QIDX.sortRecs(P); DATA.Select(P);
       if (qtype.len) CGR.recPermute(P);
       else if (!CGR.isEmpty()) wblog(FL,"ERR %s() "
          "invalid all-abelian CGR (%dx%d)",FCT,CGR.dim1,CGR.dim2
       );
    };

    unsigned cgsDimScalar(unsigned k) const;
    unsigned cgsDim(unsigned k, unsigned r) const;
    unsigned getDim(unsigned k, INDEX_T *D=NULL) const;

    unsigned getDIM(unsigned k) const {
       INDEX_T D=0; getDim(k,&D);
       return D;
    };

    void getDim(wbvector<INDEX_T> &D, wbvector<INDEX_T> *DD=NULL) const;

    wbvector<INDEX_T> getDim() const {
       wbvector<INDEX_T> D; getDim(D); return D;
    }

    void getQDim(unsigned k,
       wbMatrix<TQ> &Q,
       wbvector<INDEX_T>&S,
       wbMatrix<INDEX_T>*SC=NULL
    ) const;

    void getQDim(wbMatrix<TQ> &Q,
       wbvector<INDEX_T> &S,
       wbMatrix<INDEX_T>*SC=NULL
    ) const;

    void getDRange(TD &dmin, TD &dmax) const;

    void getDistQtot(wbMatrix<TQ> &Qtot, wbvector<double> &w) const;
    INDEX_T getDQtot() {
       wbMatrix<TQ> Q; QIDX.blockSum(QDIM,Q); Q.makeUnique();
       return Q.dim1;
    };

    wbMatrix<TQ>& getQsub(const wbindex &I, wbMatrix<TQ> &QI) const;
    wbMatrix<TQ>& getQsub(
    const unsigned k, wbMatrix<TQ> &Qk, char iflag=0) const;

    void getQsub(const wbindex &I, wbMatrix<TQ> &QI, wbMatrix<TQ> &Q2) const;

    wbMatrix<TQ> getQsub(const unsigned k) const {
       wbMatrix<TQ> Qk; getQsub(k,Qk); return Qk;
    };

    wbMatrix<TQ>& getQtot(wbMatrix<TQ> &Qtot, char uflag=0) const;
    wbMatrix<TQ>& getQsum(const wbindex &J, wbMatrix<TQ> &Q) const;

    wbMatrix<TQ>& getQfinal(const wbindex &J, wbMatrix<TQ> &Q) const;

    QSpace& getSub(const wbindex &I, QSpace &A, char ref=0) const;

    QSpace& getSubInit(const wbindex &I, QSpace &A) const;

    void saveSub2(QSpace &A, const wbindex &I);

    void initFromBlockMatrix(
       const wbarray<TD> &MM,
       const wbMatrix<TQ> &Q1,
       const wbMatrix<TQ> &Q2,
       const unsigned QDIM,
       const wbMatrix<INDEX_T> &S1,
       const wbMatrix<INDEX_T> &S2,
       const wbperm &P=wbperm()
    );

    void toBlockMatrix(
       wbarray<TD> &MM,
       unsigned K,
       wbMatrix<TQ> &Q1,
       wbMatrix<TQ> &Q2,
       wbMatrix<INDEX_T> &S1,
       wbMatrix<INDEX_T> &S2,
       wbvector<INDEX_T> &D1,
       wbvector<INDEX_T> &D2,
       const wbperm &P=wbperm()
    ) const;

    void toFull(
       wbarray<TD> &A,
       wbvector< wbMatrix<TQ> > &Q,
       wbvector< wbvector<unsigned> > &S
    ) const;

    void toFull2(
       wbarray<TD> &A,
       wbMatrix<TQ> &Q,
       wbMatrix<INDEX_T> &S
    ) const;

    void ExtendQ(const QSpace &B);

    void takeMainQS(char disp=1);

    void setQ(unsigned k, const QSpace &B, unsigned k0);
    void setQRec(unsigned r, unsigned j, const TQ *q0, const TQ *q2=NULL);

    bool hasEqualQ(unsigned i1, unsigned i2) const;

    double maxDiff(const QSpace &B) const;

    void PermuteQ(const wbperm &P);
    void permuteQ(const wbperm &P, wbMatrix<TQ>&) const;

    QSpace& Permute(const wbperm &P, char iflag=0);
    QSpace& permute(const wbperm &P, QSpace &B, char iflag=0) const;

    void permute (const char *pstr, QSpace &B, unsigned offset=1
    ) const {
       permute(wbperm(pstr,offset),B);
    };

    void Permute(const char* pstr, unsigned offset=1) {
       Permute(wbperm(pstr,offset));
    };

    void PermuteFirstTo (unsigned k) {
       unsigned r=rank(FL);
       if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
       wbperm P2(r); P2.mvFirstTo(k); Permute(P2);
    };
    void permuteFirstTo (unsigned k, QSpace &C) {
       if (this==&C) { PermuteFirstTo(k); return; }
       unsigned r=rank(FL);
       if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
       wbperm P2(r); P2.mvFirstTo(k); permute(P2,C);
    };

    void PermuteLastTo (unsigned k) {
       unsigned r=rank(FL);
       if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
       wbperm P2(r); P2.mvLastTo(k); Permute(P2);
    };
    void permuteLastTo (unsigned k, QSpace &C) {
       if (this==&C) { PermuteLastTo(k); return; }
       unsigned r=rank(FL);
       if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
       wbperm P2(r); P2.mvLastTo(k); permute(P2,C);
    };

    void checkScalarCGS(const char *F=0, int L=0);
    void FlipQ();

    void contractMat(
       unsigned ia, const QSpace &B, unsigned ib,
       QSpace &C
    ) const;

    template <class TB, class TC>
    int contract_getIdxSet(const char *F, int L,
       const ctrIdx &ica, const QSpace<TQ,TB> &B, const ctrIdx &icb,
       wbindex &Ia, wbindex &Ib, wbvector<INDEX_T> &Dc, 
       QSpace<TQ,TC> &C
    ) const;

    template<class TDB, class TDC>
    QSpace<TQ,TDC>& contract(const char *F, int L,
       const ctrIdx &ica, const QSpace<TQ,TDB> &B, const ctrIdx &icb,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
     ) const;

    template<class TDB, class TDC>
    QSpace<TQ,TDC>& contract(const char *F, int L,
       unsigned ia, const QSpace<TQ,TDB> &B, unsigned ib,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
     ) const {

       if (ia==0 ||ib==0) wblog(F,L,
          "ERR %s() index is 1-based (%d,%d)!",FCT,ia,ib);
       ctrIdx ica(1,&(--ia)), icb(1,&(--ib));
       return contract(F,L,ica,B,icb,C,P,cgflag);
    };

    template<class TDB, class TDC>
    QSpace<TQ,TDC>& contract(
       const char *sa, const QSpace<TQ,TDB> &B, const char *sb,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
    ) const { return contract(FL,sa,B,sb,C,P,cgflag); };

    template<class TDB, class TDC>
    QSpace<TQ,TDC>& contract(
       const ctrIdx &ica, const QSpace<TQ,TDB> &B, const ctrIdx &icb,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
    ) const { return contract(FL,ica,B,icb,C,P,cgflag); };

    template<class TDB, class TDC>
    QSpace<TQ,TDC>& contract(
       unsigned ia, const QSpace<TQ,TDB> &B, unsigned ib,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
    ) const { return contract(FL,ia,B,ib,C,P,cgflag); };

    template<class TDB, class TDC> inline
    QSpace<TQ,TDC>& contract(const char *F, int L,
       const char *sa, const QSpace<TQ,TDB> &B, const char *sb,
       QSpace<TQ,TDC> &C, const wbperm &P =wbperm(), char cgflag=-1
    ) const {

       ctrIdx ica(F_L,sa), icb(F_L,sb);
       return contract(F,L,ica,B,icb,C,P,cgflag);
    };

    void contract(unsigned ia, unsigned ib, QSpace &C) const;

    mxArray* QIDX_toMx() const;
    mxArray* INFO_toMx() const;
    mxArray* DATA_toMx() const;
    mxArray* DATA_save2Mx(char vflag=0);

    mxArray* toMx() const;
    mxArray* save2Mx(char vflag=0);

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    void save2MxStruct(mxArray *S, unsigned i, char tst=0, char vflag=0);

    mxArray* mxCreateCell(unsigned m, unsigned n) const {
       wblog(FL,"ERR %s()",FCT); return 0; };
    void add2MxCell(mxArray *S, unsigned i, char tst=0) const {
       wblog(FL,"ERR %s()",FCT); };

    void put (
       const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       wblog(F,L,"I/O putting '%s' to %s", vname, ws);
       put(vname,ws);
    };

    void put (const char *vname, const char *ws="caller") const {
       mxArray *S=toMx();
       mxPutAndDestroy(0,0,S,vname,ws);
    };

    wbstring sizeStrQ(const char *F=NULL, int L=0) const;

    wbstring sizeStr(char vflag=0) const;
    wbstring qStr() const { return qtype.toStr(); };

    wbstring qrec2Str(unsigned k) const {
        return QIDX.rec2Str(k,""," ",  QDIM," ;");
    };

    void disp_cgs(const char *F, int L, const char *istr) const;

    void info(const char *vname) const;
    void print(const char *vname, char dflag=0) const;

    void initOType(const char *F, int L, const mxArray* a);

    wbstring otype2Str(const char *fmt="%s") const;

    QSpace& ConjOpScalar() {
       unsigned r=rank(FL);
       if (r!=3) {
          if (r!=2) wblog(FL,
             "WRN %s() got rank-%d QSpace !?",FCT,r);
          return *this;
       }
       for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
          CGR[i].HConjOpScalar();
       }

       if (itags.len!=3) wblog(FL, 
          "ERR %s() got itags.len=%d/3 !?",FCT,itags.len);
       itags[2].Conj();

       return *this;
    };

    QSpace& SortDegQ(const char *F=0, int L=0) {
       for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
          CGR.data[i].SortDegQ(F_L); }
       return *this;
    };

    QSpace& Conj() {
       unsigned i,n; SetConjTags();
       if (typeid(TD)==typeid(wbcomplex)) {
       for (n=DATA.len,    i=0; i<n; ++i) DATA[i]->Conj(); }
       for (n=CGR.numel(), i=0; i<n; ++i) { CGR[i].Conj(); }
       return SortDegQ(FL);
    };

    QSpace& HConj() {
       QSpace B; this->save2(B);
       return B.hconj(*this);
    };

    QSpace& hconj(QSpace &B) const {
       if (!isEmpty()) {
          if (!QDIM || QIDX.dim2%QDIM) wblog(FL,
             "ERR %s() severe QSpace inconsistency (%dx%d/%d)",
             FCT, QIDX.dim1, QIDX.dim2, QDIM);
          unsigned r=QIDX.dim2/QDIM;
          if (r%2 && (r!=3 || otype!=QS_OPERATOR)) wblog(FL,"ERR %s() "
             "requires even-rank object or operator (%dx%d/%d; %s)",
              FCT, QIDX.dim1, QIDX.dim2, QDIM, QS_STR[otype]);
          transp(B); B.Conj(); B.otype=otype;
          return B;
       }
       return B.clearQSpace();
    };

    QSpace& transp(QSpace &B) const {

       if (!isEmpty()) { unsigned r=rank(FL);
          if (r%2==0) { wbperm P2;
             P2.initTranspose(r); permute(P2,B);
          }
          else if (r==3 && otype==QS_OPERATOR) {
             permute("2,1,3", B);
          }
          else wblog(FL,"ERR %s got rank-%d object",FCT,r);
       }
       else B.clearQSpace();

       return B;
    };

    const QSpace& opA(char tflag, QSpace& X) const {
       if (!strchr("NTC",tflag))
       wblog(FL,"ERR %s() invalid flag %c<%d>",FCT,tflag,tflag);

       if (tflag=='N' || isEmpty()) return (*this);
       else {
          wbperm P2("2 1"); unsigned r=rank(FL);

          if (r!=2) wblog(FL,
          "ERR %s() only applicable to rank-2 objects (%d)",FCT,r);

          permute(P2,X); if (tflag=='C') X.Conj();
          return X;
       }
    };


    wbMatrix<TQ> QIDX;

    wbvector<wbarray<TD>*> DATA;
    wbMatrix< CRef<TQ> > CGR;

    QVec qtype;

    unsigned QDIM;

    QS_TYPES otype;
    IndexLabels itags;

    char isref;



  protected:
  private:

    bool isSym_aux(
      const char *F, int L, const char *fct,
      RTD eps=1E-14, char symflag='s', char vflag=0
    ) const;

    void recSetQ(unsigned k, unsigned i);

};


template <class TQ, class TD>
template <class T> inline
void QSpace<TQ,TD>::operator*= (const T fac) {
   if (fac==T(1)) return;
   for (unsigned i=0; i<DATA.len; ++i)
   (*DATA[i]) *= fac;
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::operator+= (const QSpace &B) {
   QSpace<TQ,TD> A(*this);
   A.plus(B, *this);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::operator-=(const QSpace &B) {
   QSpace<TQ,TD> A(*this);
   A.plus(B, *this,-1);
}

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::Plus(const QSpace &B, TD bfac) {
   QSpace<TQ,TD> A(*this);
   A.plus(B, *this, bfac);
   return *this;
}

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::Minus(const QSpace &B, TD bfac) {
   QSpace<TQ,TD> A(*this);
   A.plus(B, *this,-bfac);
   return *this;
}

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::minus(
   const QSpace &B, QSpace &C, TD bfac, char vflag
) const {
   plus(B, C,-bfac,vflag); return C;
}

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::plus(
   const QSpace &B, QSpace &C, TD bfac, char vflag
 ) const {

   if (gotCGS(FL)<2) { return plus_plain(B,C,bfac,vflag); }

   C=B; if (bfac!=1.) { C*=bfac; }
   C.Append(FL,*this, vflag<1 ? 'q' : vflag);

   return C;
};


#ifdef WB_CLOCK
 WbClock wbc_qs_scprod("QS::scalarProd");
#endif

template <class TQ, class TD>
TD QSpace<TQ,TD>::scalarProd(const QSpace<TQ,TD> &B) const {

   unsigned i; int i1,i2;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   TD x=0;

   if (this==&B) return norm2();

#ifdef WB_CLOCK
 wbc_qs_scprod.resume();
#endif

   if (isEmpty() || B.isEmpty()) {
      wblog(FL,"WRN Overlap with empty object (%d,%d) !?",
      isEmpty(), B.isEmpty()); return x;
   }

   if (getQOverlap(B, QQ, IQ, 1 )) wberror(FL,str);

   for (i=0; i<IQ.dim1; ++i) {
       i1=IQ(i,0); i2=IQ(i,1);
       if (i1<0 || i2<0) wblog(FL, "ERR i1=%d, i2=%d ???", i1, i2);

       x+=(DATA[i1]->scalarProd(*B.DATA[i2]));
   }

#ifdef WB_CLOCK
 wbc_qs_scprod.stop();
#endif
   return x;
}


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::TensorProd(
   const QSpace<TQ,TD> &B,
   QSpace<TQ,TD> &C,
   const char aflag, const char bflag
) const {

   if (isEmpty() || B.isEmpty()) { C.clearQSpace(); return C; }

   if (gotCGS(FL)>0 || B.gotCGS(FL)>0) wblog(FL,
      "ERR %s() not implemented for non-abelian symmetries (got %s)",
      FCT,qtype.toStr().data);

   if (aflag!='N' || bflag!='N') {
      QSpace<TQ,TD> AX, BX;
      const QSpace<TQ,TD> *Ap=this, *Bp=&B;

      if (!strchr("NTC",aflag) || rank(FL)!=2) wblog(FL,
         "ERR invalid aflag=%c<%d> (rank-%d)",aflag,aflag,rank(FL));
      if (!strchr("NTC",bflag) || B.rank(FL)!=2) wblog(FL,
         "ERR invalid bflag=%c<%d> (rank-%d)",bflag,bflag,B.rank(FL));

      if (aflag!='N') {
          permute("2,1",AX); if (aflag=='C') AX.Conj();
          Ap=&AX;
      }
      if (bflag!='N') {
         B.permute("2,1",BX); if (bflag=='C') BX.Conj();
         Bp=&BX;
      }

      return Ap->TensorProd(*Bp,C,'N','N');
   }

   unsigned i,ia,ib,l,r;
   const unsigned n1=QIDX.dim1, n2=B.QIDX.dim1, R=rank(FL), s=QDIM*sizeof(TQ);
   const TQ *qa, *qb; TQ *q;

   if (R!=B.rank(FL) || QDIM!=B.QDIM) wblog(FL,
      "ERR %s() severe size mismatch (%d,%d)",FCT,R,B.rank(FL));

   C.clearQSpace();
   C.QIDX.init(n1*n2, 2*QIDX.dim2); C.QDIM=QDIM; C.setupDATA();

   C.initQ(FL,*this,'f',&B);

   if (itags.len || B.itags.len) { C.itags.init(2*R);
      if (!B.itags.len) {
         if (itags.len!=R) wblog(FL,
            "ERR %s() invalid itag %s",FCT, itags.toStr().data);
         for (i=0; i<R; ++i) {
            C.itags[2*i]=C.itags[2*i+1]=itags[i];
         }
      }
      else if (!itags.len) {
         if (B.itags.len!=R) wblog(FL,
            "ERR %s() invalid itag %s",FCT, B.itags.toStr().data);
         for (i=0; i<R; ++i) {
            C.itags[2*i]=C.itags[2*i+1]=B.itags[i];
         }
      }
      else {
         if (itags.len!=B.itags.len || itags.len!=R || !sameQDir(FL,B))
            wblog(FL,"ERR %s() itag/qdir mismatch (%s <> %s)",FCT,
            itags.toStr().data,B.itags.toStr().data);
         for (i=0; i<R; ++i) {
            C.itags[2*i  ]=  itags[i];
            C.itags[2*i+1]=B.itags[i];
         }
      }
   }

   for (i=ib=0; ib<n2; ++ib) {
      qb=B.QIDX.rec(ib);

      for (ia=0; ia<n1; ++ia, ++i) {
         qa=QIDX.rec(ia); q=C.QIDX.rec(i);

         for (l=r=0; r<R; ++r, l+=QDIM) {
            memcpy(q,qa+l,s); q+=QDIM;
            memcpy(q,qb+l,s); q+=QDIM;
         }
         DATA[ia]->tensorProd(*B.DATA[ib], *C.DATA[i]);
      }
   }

   return C;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::GroupIndizes(unsigned K) {
   unsigned i; int r=-1;

   if (!isConsistent(r)) wberror(FL,str);

   if (!K || unsigned(r)%K) wblog(FL,
   "ERR %s() cannot block QSpace (%d/%d)",FCT,r,K);

   QDIM=K*QDIM;

   for (i=0; i<DATA.len; ++i)
   DATA[i]->GroupIndizes(K);
}


template <class TQ, class TD>
unsigned long QSpace<TQ,TD>::getDataSize(wbMatrix<unsigned> &S) const {

   unsigned long s=0; unsigned i,m,n=DATA.len;

   if (n!=QIDX.dim1) wblog(FL,
   "ERR severe size mismatch (%d/%d)",QIDX.dim1,DATA.len);

   if (isEmpty()) { S.init(); return; }

   m=DATA[0]->SIZE.len; S.init(n,m);

   for (i=0; i<n; ++i) { wbvector<unsigned> &s=DATA[i]->SIZE;
      if (s.len!=m) wblog(FL,"ERR severe size mismatch (%d/%d)",s.len,m);
      S.recSetP(i,s.data);
      s+=s.prod(0);
   }
   return s;
};

template <class TQ, class TD> inline
unsigned long QSpace<TQ,TD>::getDataSize(wbvector<unsigned> &S) const {
   unsigned long s=0; S.init(DATA.len);
   for (unsigned i=0; i<DATA.len; ++i) { s+=(S[i]=DATA[i]->SIZE.prod()); }
   return s;
};

template <class TQ, class TD> inline
unsigned long QSpace<TQ,TD>::getDataSize() const {
   unsigned long s=0;
   for (unsigned i=0; i<DATA.len; ++i) s+=DATA[i]->SIZE.prod();
   return s;
};

template <class TQ, class TD> inline
unsigned long QSpace<TQ,TD>::getCGSSize() const {
   if (gotCGS(FL)>0) {
      unsigned i,j; unsigned long s=0;
      for (i=0; i<CGR.dim1; ++i) {
      for (j=0; j<CGR.dim2; ++j) s+=CGR(i,j).numel(); }
      return s;
   }
   return 0;
};

template <class TQ, class TD> inline
wbstring QSpace<TQ,TD>::totSize2Str() const {

   unsigned long ss[3]={
       sizeof(TD)*getDataSize(),
       sizeof(TQ)*(QIDX.dim1*QIDX.dim2),
       sizeof(double)*getCGSSize()
   };
   memsize2Str(ss[0]+ss[1]+ss[2],str);
   return wbstring(str);
};


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isZero(const double eps, char flag) const {

   for (unsigned i=0; i<DATA.len; ++i)
   if (!(DATA[i]->isZero(eps,flag))) return 0;

   return 1;
}

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::hasZeroData(const double eps, char flag) const {

   for (unsigned i=0; i<DATA.len; ++i)
   if (DATA[i]->isZero(eps,flag)) return 1;

   for (unsigned n=CGR.numel(), i=0; i<n; ++i)
   if (CGR[i]->isZero(eps,flag)) return 1;

   return 0;
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::Set(
    unsigned k, const wbvector<TQ> &Q, const wbarray<TD> &D
){

   QIDX.recSet(k,Q);

   if (k>=DATA.len)
   wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,DATA.len);

   DATA[k]->init(D);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::SetQ(unsigned i, unsigned k, const TQ *d){

   if (i>=QIDX.dim1 || k>=rank(FL)) wblog(FL,
   "ERR index out of bounds (%d/%d; %d;%d)",i,QIDX.dim1,k,rank(FL));

   memcpy(QIDX.data+QIDX.dim2*i+QDIM*k, d, QDIM*sizeof(TQ));
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::Sort(const QSpace<TQ,TD> &R) {
   int e;
   wbindex P,Ir; wbperm Pr;

   if (QIDX.dim1==0) return;

   e=matchIndex(QIDX,R.QIDX, P, Ir);

   if (e || P.len!=QIDX.dim1) wblog(FL,
   "ERR QSpace R not complete to be used as sorting reference (%d/%d).",
    P.len, QIDX.dim1);

   Ir.Sort(Pr); P.Select(Pr);

   QIDX.recPermute(P); if (!CGR.isEmpty())
   CGR .recPermute(P);
   DATA.Select(P);
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDim(
   wbvector<INDEX_T> &D, wbvector<INDEX_T> *DD
 ) const {

   unsigned k; int r=-1;

   if (!isConsistent(r)) {
      wblog(FL,"ERR %s() %s",FCT,str);
   }
   if (r<0) wblog(FL,"ERR invalid rank r=%d !??",r);

   if (otype==QS_OPERATOR) r=rank(FL);

   D.init(r);
   if (DD) { DD->init(r);
          for (k=0; (int)k<r; ++k) D[k]=getDim(k,DD->data+k); }
   else { for (k=0; (int)k<r; ++k) D[k]=getDim(k); }
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::getDim(unsigned k, INDEX_T *D_) const {

   INDEX_T i,s, l=0, Dm=0, Ds=0;
   wbvector<INDEX_T> dd;
   wbMatrix<TQ> Qk;
   wbperm pp;

   char cgflag=gotCGS(FL);

   if (QIDX.dim1==0) { if (D_) *D_=0; return 0; }

   getQsub(k,Qk).SkipTiny_float().groupRecs(pp,dd);


   if (cgflag==1 && qtype.hasMP()>1) { cgflag=2; }
   if (cgflag>0) {
      unsigned j,ip,e=0,dc,r=rank(FL);

      for (l=i=0; i<dd.len; ++i, l+=dd[i-1]) { ip=pp[l];
         const wbarray<TD> &a=*DATA[ip];
         if (a.SIZE.len!=r) wblog(FL,
            "ERR %s() rank mismatch (%s; %d)",FCT,a.sizeStr().data,r);
         s=a.SIZE[k]; Dm+=s;

         for (dc=1,j=0; j<CGR.dim2; ++j) {
            dc*=CGR(ip,j).Size(k);

         }
         if (cgflag<=1) {
            unsigned d0=qtype.QDim(QIDX.ref(ip,QDIM*k)); if (d0!=dc)
            wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,dc,d0);
         }
         Ds+=(s*dc);
      }

      if (e) wblog(FL,
      "WRN %s() CGR data not yet initialized !?? (%d)",FCT,e);
   }
   else {
      for (l=i=0; i<dd.len; ++i, l+=dd[i-1]) {
         Dm+=(s=DATA[pp[l]]->SIZE[k]);
         if (qtype.len) Ds += s * qtype.QDim(Qk.rec(i));
      }
   }

   if (D_) { D_[0] = (Ds ? Ds : Dm); }

   return Dm;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::cgsDim(unsigned k, unsigned r) const {
   
   if (k>=QIDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,QIDX.dim1); 

   unsigned d=1, R=rank(FL);
   if (r>=R) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,r,R);
   if (!CGR.dim2) return d;

   if (CGR.dim1!=QIDX.dim1) wblog(FL,
      "ERR %s() size mismatch (%dx%d <> %dx%d)",
       FCT,QIDX.dim1,QIDX.dim2,CGR.dim1,CGR.dim2); 

   for (unsigned j=0; j<CGR.dim2; ++j) {
      const cdata__ &c=(*CGR(k,j)); if (c.isScalar()) continue;
      if (c.SIZE.len!=R) wblog(FL,
         "ERR %s() rank mismatch (%s; %d)",FCT,c.sizeStr().data,R); 
      d*=c.SIZE.data[r];
   }

   return d;
};


template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::cgsDimScalar(unsigned k) const {
   
   if (k>=QIDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,QIDX.dim1); 
   unsigned d=1, r=rank(FL);

   if (r!=2) wblog(FL,"ERR %s() invalid operator rank-%d",FCT,r);
   if (CGR.dim2) {
      if (CGR.dim1!=QIDX.dim1) wblog(FL,"ERR %s() size mismatch "
         "(%dx%d <> %dx%d)",FCT,QIDX.dim1,qtype.len,CGR.dim1,CGR.dim2);
      for (unsigned j=0; j<CGR.dim2; ++j) {
         if (CGR(k,j).cgr) { d*=CGR(k,j).cgr->cgd.dim(); }
      }
   }

   if (!d) wblog(FL,"ERR %s() got dim(cgd)=%d !?",FCT,d);

   return d;
};


template <class TQ, class TD>
void getQDimGen(wbvector< const QSpace<TQ,TD>* > &A,
   wbMatrix<TQ> &Q, wbvector<INDEX_T> &S,
   wbindex I,
   wbMatrix<INDEX_T> *SC_=NULL
){
   if (!A.len) wblog(FL,"ERR got empty input space");

   unsigned i,j,k,l,ic,d, d2=0, m,nq,p,i1,i2, n=0; INDEX_T *ip;
   unsigned r=-1, rx, QDIM=A[0]->QDIM;
   const QS_TYPES &otype=A[0]->otype;
   const QVec &qtype=A[0]->qtype;
   wbMatrix<INDEX_T> S0,SC,II;
   wbvector<INDEX_T> D;
   wbperm P;
   
   char cgflag=(SC_!=0);

   nq=A[0]->CGR.dim2;

   for (k=0; k<A.len; ++k) { const QSpace<TQ,TD> &Ak=(*A[k]);
      if (!A[k]) wblog(FL,"ERR %s() got null object",FCT);
      n+=Ak.QIDX.dim1;

      if (I.len) { r=rx=Ak.rank(FL); }
      else { rx=r;
         if (int(r=Ak.isOperator(&rx,'x'))<=0) wblog(FL,
            "ERR %s() got invalid rank-%d operator (%d/%d; '%s')",
            FCT,rx,k+1,A.len,A[0]->otype2Str().data
         );
         if (Ak.otype!=otype) wblog(FL,
            "ERR %s() object type inconsistency (%s ; %s)",
            FCT, Ak.otype2Str().data, A[0]->otype2Str().data
         );
      }

      if (!k) d2=r*QDIM;

      if (Ak.qtype!=qtype) wblog(FL,
         "ERR %s() Q-type inconsistency (%s ; %s)",
          FCT, Ak.qStr().data, qtype.toStr().data);
      if (Ak.QDIM!=QDIM || Ak.QIDX.dim2<d2) wblog(FL,
         "ERR %s() rank inconsistency (%d: %d/%d %d/%d)",
          FCT, k+1, Ak.QDIM, QDIM, Ak.QIDX.dim2, d2);
      if (Ak.QIDX.dim1!=Ak.DATA.len)
          wblog(FL,"ERR %s() size mismatch (%dx%d <> %d)",
          FCT, Ak.QIDX.dim1, Ak.QIDX.dim2, Ak.DATA.len
      );

      if (Ak.gotCGS(FL)<=0) {
         if (k==0 && cgflag) cgflag=0;
         else if (cgflag) {
            wblog(FL,"ERR %s() CGR inconsistency (%d)",FCT,k+1);
         }
      }
      if (cgflag && k && Ak.CGR.dim2!=nq) wblog(FL,
      "ERR %s() CGR inconsistency (%d: %d/%d)",FCT,k+1,Ak.CGR.dim2,nq); 
   }

   if (!r) wblog(FL,"ERR got empty QSpaces !?");

   if (!I.len) I.Index(r);
   else if (I.anyGE(r)) wblog(FL,
     "ERR %s() index out of bounds (%s; %d)",FCT,(I+1).toStr().data,r);

   S0.init(n,I.len);
   II.init(n,I.len*2);
   Q .init(n,I.len*QDIM);

   if (cgflag) SC.init(n,nq*I.len).set(1);

   for (l=k=0; k<A.len; ++k) {
      const QSpace<TQ,TD> &Ak=(*A[k]);
      const wbMatrix<TQ> &qq=Ak.QIDX; m=Ak.DATA.len;

      for (i=0; i<m; ++i, ++l) {
         const wbvector<INDEX_T> &s=Ak.DATA[i]->SIZE;
         if (s.len!=rx) wblog(FL,
            "ERR %s() rank mismatch (data[%d]: %s/%d)",
             FCT,i+1,Ak.DATA[i]->sizeStr().data, rx);
         Q .recSetB(l,I,QDIM, qq.rec(i));
         S0.recSetB(l,I,1, s.data);

         if (cgflag) {
            for (j=0; j<nq; ++j) {
               wbvector<INDEX_T> s; Ak.CGR(i,j).getSize(s);
               if (!s.len || s.allEqual(1)) continue;
               if (s.len<r || s.len>rx) wblog(FL,
                  "ERR %s() rank mismatch (cgr(%d,%d): %s/%d)",
                   FCT,i+1,j+1, Ak.CGR(i,j).sizeStr().data, r
               );
               for (ic=j, i1=0; i1<I.len; ++i1, ic+=nq) {
                  SC(l,ic)=s[I.data[i1]];
               }
            }
         }

         ip=II.rec(l);
         for (p=j=0; j<II.dim2; j+=2, ++p) { ip[j]=k; ip[j+1]=I[p]; }
      }
   }


   Q .Reshape(n*I.len,QDIM);
   S0.Reshape(n*I.len,1);
   II.Reshape(n*I.len,2);

   if (cgflag)       
   SC.Reshape(n*I.len,nq);


   Q.groupRecs(P,D); S.init(D.len);
   if (cgflag) SC_->init(D.len,nq);


   for (l=i=0; i<D.len; ++i) { d=D[i]; i1=P[l++];
      S[i]=S0(i1,0);
      if (cgflag) SC_->recSetP(i,SC.ref(i1));
      for (j=1; j<d; ++j, ++l) { i2=P[l];
         if (S0.recCompare(i1,i2)) wblog(FL,
            "ERR %s() size inconsistency\nDATA: %d:%d <> %d:%d",
             FCT, II(i1,0)+1, II(i1,1)+1, II(i2,0)+1, II(i2,1)+1
         );
         if (cgflag && SC.recCompare(i1,i2)) wblog(FL,
            "ERR %s() size inconsistency\nCGS: %d:%d <> %d:%d",
             FCT, II(i1,0)+1, II(i1,1)+1, II(i2,0)+1, II(i2,1)+1
         );
      }
   }
};


template <class TQ, class TD> inline
void getQDimGen(const wbvector< QSpace<TQ,TD> > &A,
   wbMatrix<TQ> &Q, wbvector<INDEX_T> &S,
   const wbindex &I,
   wbMatrix<INDEX_T> *SC=NULL
){
   wbvector< const QSpace<TQ,TD>* > Ap(A.len);
   for (unsigned i=0; i<A.len; ++i) Ap[i]=&A[i];
   getQDimGen(Ap,Q,S,I,SC);
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQDim(wbMatrix<TQ> &Q,
   wbvector<INDEX_T> &S,
   wbMatrix<INDEX_T>*SC
) const {

   wbindex I;
   wbvector< const QSpace* > A(1); A[0]=this;

   getQDimGen(A,Q,S,I,SC);
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQDim(unsigned k,
   wbMatrix<TQ> &Q,
   wbvector<INDEX_T> &SD,
   wbMatrix<INDEX_T> *SC
) const {

   unsigned i,j,q, r=rank(FL),n=DATA.len;
   wbvector<INDEX_T> D;
   wbindex Ig;
   wbperm P;

   if (k>=r) wblog(FL,"ERR %s() dim out of bounds (%d/%d)",FCT,k,r);
   if (QIDX.dim1!=n) wblog(FL,
      "ERR %s() severe QSpace inconsistency (%d/%d)",FCT,QIDX.dim1,n);

   getQsub(k,Q).groupRecs(P,D,-1,1,&Ig);
   SD.init(D.len).set(-1);

   for (i=0; i<n; ++i) {
      const wbvector<INDEX_T> &si=DATA[i]->SIZE;
      INDEX_T &s=SD.el(Ig[i]);
      if (si.len>r || k>=si.len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d,%d)",FCT,k+1,si.len,r); 
      if (int(s)>=0) {
         if (s!=si.data[k]) wblog(FL,
            "ERR %s() size mismatch DATA(%d): %s [S(%d)=%d]",
            FCT,i+1,DATA[i]->sizeStr().data, k+1, s
         ); 
      }
      else s=si.data[k];
   }

   if (SC==NULL) return;
   if (gotCGS(FL)<=0) { SC->init(); return; }

   SC->init(D.len,CGR.dim2).set(-1);

   for (i=0; i<n; ++i) {
      for (j=0; j<CGR.dim2; ++j) {
         INDEX_T &c=SC->el(Ig[i],j);
         q=CGR(i,j).qdim(k,&qtype[j]);

         if (int(c)>=0) {
            if (c!=q) wblog(FL,
               "ERR %s() size mismatch CGR(%d,%d): %s [S(%d)=%d]",
               FCT,i+1,j+1, CGR(i,j).sizeStr().data, k+1, q
            );
         } else c=q;
      }
   }
};


template <class TQ, class TD> inline
int QSpace<TQ,TD>::isSingleton(unsigned k) const {
   unsigned i;

   if (k>=rank(FL)) wblog(FL,
   "ERR %s() index out of bounds (%d,%d)",FCT,k,rank(FL));

   for (i=0; i<DATA.len; ++i)
   if (DATA[i]->SIZE[k]!=1) break;

   return (i==DATA.len);
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDRange(TD &dmin, TD &dmax) const {
    unsigned i,n=DATA.len;

    if (!isConsistent()) { info("this"); wberror(FL,str); }

    for (i=0; i<n; ++i) if (DATA[i] && !DATA[i]->isEmpty()) {
       dmin=dmax=DATA[i]->data[0];
       break;
    }

    for (i=0; i<n; ++i) if (DATA[i] && !DATA[i]->isEmpty()) {
       dmin=MIN(dmin, DATA[i]->min());
       dmax=MAX(dmax, DATA[i]->max());
    }
}


template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::getSub(
   const wbindex &I, QSpace<TQ,TD> &A, char ref
 ) const {

   if ((void*)this==(void*)&A) {
      QSpace<TQ,TD> X; A.save2(X);
      X.getSub(I,A,ref); return A;
   }

   unsigned i,j;
   char gotcg=gotCGS(FL);

   A.init();
   QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype; A.itags=itags;

   A.isref=(ref ? 1 : 0);

   A.setupDATA();

   if (gotcg>0) {
      if (CGR.dim1!=DATA.len) wblog(FL,
         "ERR %s() severe size mismatch (CGR: %d/%dx%d)",
          FCT,DATA.len,CGR.dim1,CGR.dim2);
      A.setupCGR();
   }

   for (i=0; i<I.len; ++i) {
      if (I[i]>=DATA.len) wblog(FL,
         "ERR %s() index out of bounds (%d,%d)",FCT,I[i],DATA.len);

      if (A.isref) 
           A.DATA[i]=DATA[I[i]];
      else A.DATA[i]->init(*DATA[I[i]]);

      if (gotcg>0) {
         for (j=0; j<CGR.dim2; ++j) {
            A.CGR(i,j)=CGR(I[i],j);
         }
      }
   }

   return A;
};


template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::getSubInit(
    const wbindex &I,
    QSpace<TQ,TD> &A
) const {

    if (this==&A) {
       QSpace<TQ,TD> X; A.save2(X);
       return X.getSubInit(I,A);
    }

    unsigned i,j;
    char cgflag=gotCGS(FL);

    A.init();
       QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype;
    A.setupDATA();

    if (cgflag) { A.setupCGR();

       for (i=0; i<I.len; ++i) {
          if (I[i]>=QIDX.dim1) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,I[i]+1,QIDX.dim1);
          for (j=0; j<CGR.dim2; ++j) A.CGR(i,j)->init(*CGR(I[i],j));
       }
    }

    return A;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::saveSub2(QSpace<TQ,TD> &A, const wbindex &I) {

    if (this==&A) wblog(FL,
       "ERR %s() got overlapping input and output space!",FCT);
    if (isref) wblog(FL,"ERR %s() called for ISREF!",FCT); 

    if (!I.len) { A.init(); return; }

    unsigned i,j;
    char cgflag=gotCGS(FL);

    A.init();
    QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype;
    A.setupDATA();
    A.setupCGR();

    for (i=0; i<I.len; ++i) {
       if (I[i]>=DATA.len) wblog(FL,
          "ERR %s() index out of bounds (%d,%d)",FCT,I[i],DATA.len);

       DATA[I[i]]->save2(*A.DATA[i]);
       WB_DELETE_1(DATA[I[i]]);

       if (cgflag) {
          for (j=0; j<CGR.dim2; ++j) { cdata__* &c=CGR(I[i],j); 
             if (DATA[I[i]]->isref || c->isref) wblog(FL,
                "ERR %s() got ISREF space (%d/%d, %d)",FCT,i+1,I.len,j+1); 
             c->save2(*A.CGR(i,j));
             WB_DELETE_1(c);
          }
       }
    }

    wbindex I2; I.invert(QIDX.dim1,I2);

    QIDX.Set2Recs(I2);
    DATA.Select(I2); if (cgflag) { CGR.Set2Recs(I2); }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::toBlockMatrix(
   wbarray<TD> &MM,
   unsigned r2,
   wbMatrix<TQ> &Q1,
   wbMatrix<TQ> &Q2,
   wbMatrix<INDEX_T> &S1,
   wbMatrix<INDEX_T> &S2,
   wbvector<INDEX_T> &D1,
   wbvector<INDEX_T> &D2,
   const wbperm &P
 ) const {

   INDEX_T i,j,d,l,M,N, r,s, dim1,dim2; unsigned rk=-1;
   wbvector<size_t> s1,s2;
   wbvector<INDEX_T> d1,d2, cD1,cD2, g1, g2;
   wbMatrix<unsigned> mark;
   wbarray<TD> Mi;
   wbindex I1,I2;
   wbperm P1, P2, iP1, iP2;
   char cgflag=gotCGS(FL);

   if (!isConsistent((int&)rk)) wblog(FL,"%s",str);
   if (isEmpty()) {
      if (r2 || !P.isEmpty()) wblog(FL,
         "ERR data inconsistency (%d,%d)", r2, P.len);
      MM.init(); Q1.init(); Q2.init(); S1.init(); S2.init();
      return;
   }

   if (int(r2)<0) r2=rk/2;

   if (r2>rk || int(rk)<=0) wblog(FL,
      "ERR index out of bounds (%d/%d)",r2,rk);
   if (!r2 || r2==rk) wblog(FL,
      "ERR got empty group (%d/%d)",r2,rk);

   if (cgflag>0) {
      if (gotCGX() || rk!=2) wblog(FL,
      "ERR %s() got unexpected non-abelian CGR (%d,%d)",FCT,cgflag,rk);
   }

   if (!P.isEmpty()) {
      if (P.len!=unsigned(rk)) wblog(FL,
         "ERR Dimensions mismatch (%d,%d)",P.len,rk);
      I1.init(r2, P.data);
      I2.init(unsigned(rk)-r2, P.data+r2);
   }
   else {
      I1.Index(0, r2-1);
      I2.Index(r2,unsigned(rk)-1);
   }

   getQsub(I1,Q1).groupRecs(P1,d1); P1.getIPerm(iP1); M=d1.len;
   getQsub(I2,Q2).groupRecs(P2,d2); P2.getIPerm(iP2); N=d2.len;

   S1.init(M,I1.len); D1.init(M);
   S2.init(N,I2.len); D2.init(N);

   g1.init(P1.len); g2.init(P2.len);
   for (l=i=0; i<d1.len; ++i) for (d=d1[i], j=0; j<d; ++j) g1[l++]=i;
   for (l=i=0; i<d2.len; ++i) for (d=d2[i], j=0; j<d; ++j) g2[l++]=i;

   g1.Permute(iP1);
   g2.Permute(iP2); mark.init(M+1,N+1);

   for (i=0; i<P1.len; ++i) {
      r=g1[i]; s=g2[i];

      DATA[i]->SIZE.select(I1,s1);
      DATA[i]->SIZE.select(I2,s2);

      if (mark(r,s)++) wblog(FL,"ERR %s() block %d not unique",FCT,i+1);

      if (!mark(r,N)++) { S1.recSet(r,s1); D1[r]=s1.prod(0); }
      else if (S1.recCompare(r,s1)) wblog(FL,
      "ERR %s() block size mismatch (%d, r: %d/%d)\n[%s, %s]",
       FCT, i, r+1, M, S1.rec2Str(r,"","x").data, s1.toStrD().data);

      if (!mark(M,s)++) { S2.recSet(s,s2); D2[s]=s2.prod(0); }
      else
      if (S1.recCompare(r,s1)) wblog(FL,
      "ERR %s() block size mismatch (%d, c: %d/%d)\n[%s, %s]",
       FCT, i, s+1, N, S1.rec2Str(r,"","x").data, s1.toStrD().data);
   }

   dim1=D1.sum(); dim2=D2.sum();
   MM.init(dim1,dim2); D1.cumsum0(cD1); D2.cumsum0(cD2);

   for (i=0; i<DATA.len; ++i) {
      DATA[i]->toMatrixRef(Mi,r2,P);
      Mi.copyStride( MM.ref(cD1[g1[i]], cD2[g2[i]]), dim1);
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::initFromBlockMatrix(
   const wbarray<TD> &MM,
   const wbMatrix<TQ> &Q1,
   const wbMatrix<TQ> &Q2,
   const unsigned QD,
   const wbMatrix<INDEX_T> &S1,
   const wbMatrix<INDEX_T> &S2,
   const wbperm &P
){
   INDEX_T i,j,k,l,m,n,s, dim1, dim2;
   wbvector<INDEX_T> D1, D2, cD1, cD2, S;
   wbvector<TQ> Q;
   wbarray<TD> Ai;
   const TD *d0; TD *d;

   char cgflag=gotCGS(FL);
   if (cgflag>0) {
      if (Q1.dim1!=1 || Q2.dim1!=1) wblog(FL,
         "ERR %s() got unexpected non-abelian data (%d,%d; %d)",
         FCT,Q1.dim1,Q2.dim1,cgflag
      );
   }

   if (Q1.dim1!=S1.dim1 || Q2.dim1!=S2.dim1) wblog(FL,
      "ERR severe block index mismatch (%d,%d; %d,%d)",
   Q1.dim1, Q2.dim1, S1.dim1, S2.dim1);

   D1.init(S1.dim1);
   for (i=0; i<S1.dim1; ++i) D1[i]=S1.recProd(i);
   D1.cumsum0(cD1);

   D2.init(S2.dim1);
   for (i=0; i<S2.dim1; ++i) D2[i]=S2.recProd(i);
   D2.cumsum0(cD2);

   if (MM.rank()!=2) wblog(FL,
      "ERR %s() only applicable to rank-2 objects (%d).",FCT,MM.rank());
   dim1=MM.SIZE[0]; dim2=MM.SIZE[1];

   if (dim1!=D1.sum() || dim2!=D2.sum()) wblog(FL,
      "ERR severe block size mismatch (%dx%d; %dx%d) ???",
       dim1, dim2, D1.sum(), D2.sum());
   if (!QD || Q1.dim2%QD || Q2.dim2%QD) wblog(FL,
      "ERR QDIM does not match size (%d; %d,%d)",
       QD, Q1.dim2, Q2.dim2);

   init(D1.len*D2.len, S1.dim2+S2.dim2, QD);

   for (k=i=0; i<D1.len; ++i)
   for (  j=0; j<D2.len; ++j, ++k) {
       m=D1[i]; n=D2[j]; if (m==0 || n==0) continue;

       Ai.init(m,n); d=Ai.data; s=m*sizeof(TD);
       d0=MM.ref(cD1[i],cD2[j]);

       for (l=0; l<n; ++l, d+=m, d0+=dim1) 
       memcpy(d,d0,s);

       S.init(S1.dim2, S1.rec(i), S2.dim2, S2.rec(j));
       Q.init(Q1.dim2, Q1.rec(i), Q2.dim2, Q2.rec(j));

       QIDX.recSet(k,Q);
       Ai.Reshape(S); Ai.save2(*DATA[k]);
   }

   for (k=i=0; i<DATA.len; ++i) {
      if (DATA[i]->isEmpty()) continue;
      if (k<i) { QIDX.recSet(k,i); DATA[k]=DATA[i]; }
      ++k;
   }

   if (k<DATA.len) {
      QIDX.Resize(k,QIDX.dim2);
      DATA.Resize(k);
   }

   if (!P.isEmpty()) Permute(P);
};


template <class TQ, class TD>
void QSpace<TQ,TD>::toFull(
   wbarray<TD> &A,
   wbvector< wbMatrix<TQ> > &QB,
   wbvector< wbvector<unsigned> > &SB
) const {

   if (gotCGS(FL)>0 || gotCGX()) wblog(FL,
      "ERR %s() not implemented yet for %d %d",FCT,
      gotCGS(FL), gotCGX(FL)
   );

   INDEX_T i,j,k,r,M=QIDX.dim1; int ir=-1;
   wbvector<INDEX_T> N,SA,D;
   wbarray<unsigned> mark;
   wbMatrix<INDEX_T> G;
   wbindex I,J;
   wbperm P,iP;

   if (!isConsistent(ir)) wblog(FL,"%s",str);
   if (isEmpty()) {
      A.init(); QB.init(); SB.init();
      return;
   }

   r=(unsigned)ir;

   QB.initDef(r); SB.initDef(r); D.init(r);
   I.init(r); P.init(r); iP.init(r); G.init(M,r); N.init(r);

   for (k=0; k<r; ++k) {
      getQsub(k,QB[k]).groupRecs(P,D); P.getIPerm(iP);

      I.BlockIndex(D).Permute(iP);
      G.setCol(k,I.data);

      N[k]=D.len; SB[k].init(N[k]);
   }

   mark.init(N+1);

   for (i=0; i<M; ++i) {
      G.getRec(i,I);

      if (mark(I)++) wblog(FL,
         "ERR %s() block (%s) not unique (%d)",
          FCT,(I+1).toStrf("",",").data, mark(I)
      );

      J.wbvector<INDEX_T>::operator=(N);

      for (k=0; k<r; ++k) {
         wbvector<INDEX_T> &SBk=SB[k];

         J[k]=j=I[k]; if (k) J[k-1]=N[k-1];
         if (!mark(J)++) {
            SBk[j]=DATA[i]->SIZE[k];
         }
         else if (SBk[j]!=DATA[i]->SIZE[k]) wblog(FL,
           "ERR %s - block size mismatch (%d: %d,%d)",
            FCT,i+1, SBk[j], DATA[i]->SIZE[k]
         );
      }
   }

   SA.init(r); for (k=0; k<r; ++k) SA[k]=SB[k].sum();
   A.init(SA);

   for (i=0; i<DATA.len; ++i) {
      G.getRec(i,I);
      A.setBlock(SB, I, *DATA[i]);
   }
}


template <class TQ, class TD>
void QSpace<TQ,TD>::toFull2(
   wbarray<TD> &A,
   wbMatrix<TQ> &QB,
   wbMatrix<INDEX_T> &SB
) const {

   if (gotCGS(FL)>0 || gotCGX()) wblog(FL,
      "ERR %s() not implemented yet for %d %d",FCT,
      gotCGS(FL), gotCGX(FL)
   );

   unsigned i,j,k,r,K,N,M=QIDX.dim1; int ir=-1;
   wbMatrix<TQ> Q1,Q2;
   wbvector< wbvector<INDEX_T> > S2;
   wbvector<INDEX_T> s,D,SA;
   wbMatrix<unsigned> mark; wbMatrix<INDEX_T> IJ;
   wbindex I,I1,I2;
   wbarray<TD> a;
   wbperm P,iP;

   if (!isConsistent(ir)) wblog(FL,"%s",str);
   if (isEmpty()) {
      A.init(); QB.init(); SB.init();
      return;
   }

   if (ir%2) wblog(FL,"ERR %s requires even-rank object.", FCT);
   r=(unsigned)ir;
   K=r/2;

   I.Index(0,K-1); getQsub(I,Q1);
   I.Index(K,r-1); getQsub(I,Q2);

   QB.Cat(1,Q1,Q2);

   QB.groupRecs(P,D); P.getIPerm(iP);

   I.BlockIndex(D); I.Select(iP);

   I1.init(M, I.data);
   I2.init(M, I.data+M);

   N=D.len;
   mark.init(N,N+1);
   SB.init(N,K);
   IJ.init(M,2);

   for (k=0; k<M; ++k) {
      i=IJ(k,0)=I1[k]; j=IJ(k,1)=I2[k];

      if (mark(i,j)++) wblog(FL,
      "WRN %s - block (%d,%d) not unique (%d)", FCT, i,j, mark(i,j));

      s.init2ref(K, DATA[k]->SIZE.data);
      if (mark(i,N)++) {
         if (!SB.recEqual(i,s.data)) wblog(FL,
         "ERR %s - block size mismatch (%d: %s != %s)",FCT,k+1,
         SB.rec2Str(i,"","x").data, s.toStrD().data);
      }
      else SB.recSet(i,s);

      s.init2ref(K, DATA[k]->SIZE.data+K);
      if (mark(j,N)++) {
         if (!SB.recEqual(j,s.data)) wblog(FL,
         "ERR %s - block size mismatch (%d: %s != %s)",FCT,k+1,
         SB.rec2Str(j,"","x").data, s.toStrD().data);
      }
      else SB.recSet(j,s);
   }

   S2.init(2); SB.recProd(S2[0]); S2[1]=S2[0];
   SA.init(2); SA[0]=SA[1]=S2[0].sum();

   A.init(SA);

   for (i=0; i<DATA.len; ++i) {
      I.init2ref(2,IJ.rec(i));

      a.init2ref(*DATA[i]);
      a.GroupIndizes(K);

      A.setBlock(S2,I,a);
   }
}


template <class TQ, class TD>
void QSpace<TQ,TD>::ExtendQ(const QSpace<TQ,TD> &B) {

   unsigned i,k,m,n=B.QIDX.dim1;
   wbvector<char> mark;
   wbindex Ia,Ib;

   if (B.isEmpty()) return;
   if (isEmpty()) {
      (*this)=B;
      for (i=0; i<DATA.len; ++i) DATA[i]->set(TD(0));
      return;
   }

   if (QIDX.dim2!=B.QIDX.dim2) wblog(FL,
   "ERR %s - size mismatch (%d/%d)", FCT, QIDX.dim2, B.QIDX.dim2);
   
   matchIndex(QIDX,B.QIDX,Ia,Ib);

   mark.init(n);
   for (i=0; i<Ib.len; ++i) mark[Ib[i]]++;

   m=mark.count(0);
   if (m==0) return;

   k=QIDX.dim1;
   Append(m);

   for (i=0; i<n; ++i) {
      if (mark[i]) continue;
      QIDX.recSetP(k,B.QIDX.rec(i));
      DATA[k]->init(B.DATA[i]->SIZE);
      ++k;
   }
};


template <class TQ, class TD>
void QSpace<TQ,TD>::takeMainQS(char disp) {

   unsigned i,j,n,i0,imax=0;
   wbvector<INDEX_T> d, dc, I;
   wbperm P;

   wbMatrix<TQ> QT;
   wbvector<TD> Z;
   TD z, zmax=0, ztot=0, eps=1E-3;

   getQtot(QT);
   QT.groupRecs(P,d);

   d.cumsum0(dc); Z.init(d.len);

   for (i=0; i<d.len; ++i) {
      for (z=0, n=d[i], i0=dc[i], j=0; j<n; ++j)
      z+=DATA[P[i0+j]]->norm2();

      Z[i]=z;

      if (z>zmax) { zmax=z; imax=i; }
      ztot+=z;
   }


   if (disp) {
      z=fabs((ztot-zmax)/ztot);
      sprintf(str,"%s Skipping weight %.3g of total %.3g (%.3g%s)\n"
      "    Q=[%s]", (z>eps ? "WRN" : " · "), ztot-zmax, ztot, z*100.,
      z>0.1 ? "% !!?" : "%", QT.rec2Str(imax).data);

      wblog(FL,str);
   }

   I.init(d[imax], P.data+dc[imax]);

   QIDX.Set2Recs(I);
   DATA.Select(I);
}


template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsub(
    const unsigned k, wbMatrix<TQ> &Qk, char iflag
) const {

    unsigned i, r=rank(FL), s=QDIM*sizeof(TQ), offset=k*QDIM;

    if (iflag) {
       wbindex I; I.Index_ex(r,k); getQsub(I,Qk);
       return Qk;
    }

    if (k>=r) wblog(FL,
       "ERR %s() index out of bounds (%dx%d/%d; %d)",
        FCT, QIDX.dim1, QIDX.dim2, QDIM, k);

    Qk.init(QIDX.dim1, QDIM);

    for (i=0; i<QIDX.dim1; ++i)
    memcpy(Qk.rec(i), QIDX.rec(i)+offset, s);

    return Qk;
};


template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsub(
    const wbindex &I, 
    wbMatrix<TQ> &QI  
) const {
    unsigned i,k, r=rank(FL);
    TQ *q;

    if (&QI==&QIDX) wblog(FL,
    "ERR %s() QIN and QOUT space are the same!",FCT);
    for (i=0; i<I.len; ++i) if (I[i]>=r) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,I[i],r);

    QI.init(QIDX.dim1, I.len*QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        const TQ *q0=QIDX.rec(i); q=QI.rec(i);

        for (k=0; k<I.len; ++k, q+=QDIM)
        cpyRange(q0+I[k]*QDIM, q, QDIM);
    }

    return QI;
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQsub(
    const wbindex &I,
    wbMatrix<TQ> &QI,
    wbMatrix<TQ> &Q2
) const {

    unsigned i, k, r=rank(FL), s=QDIM*sizeof(TQ);
    wbindex I2;
    TQ *q1,*q2; const TQ *q0;


    I.invert(r,I2);

    if (&QI==&QIDX || &Q2==&QIDX)
    wblog(FL,"ERR getQsub requires distinct QOUT space!");

    QI.init(QIDX.dim1, I .len*QDIM);
    Q2.init(QIDX.dim1, I2.len*QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        q0=QIDX.rec(i); q1=QI.rec(i); q2=Q2.rec(i);

        for (k=0; k<I .len; ++k) memcpy(q1+k*QDIM, q0+I [k]*QDIM, s);
        for (k=0; k<I2.len; ++k) memcpy(q2+k*QDIM, q0+I2[k]*QDIM, s);
    }
}


template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQtot(
    wbMatrix<TQ> &Qtot,
    char uflag
) const {

    unsigned i,k,l,r, r0=rank(FL);
    TQ *q, *q0;

    if (&Qtot==&QIDX)
    wblog(FL,"ERR getQtot requires distinct QOUT space!");

    Qtot.init(QIDX.dim1, QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        q0=QIDX.rec(i); q=Qtot.rec(i);
        for (l=r=0; r<r0; ++r)
        for (k=0; k<QDIM; ++k, ++l) q[k]+=q0[l];
    }

    if (uflag) Qtot.makeUnique();

    return Qtot;
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDistQtot(
   wbMatrix<TQ> &Q, wbvector<double> &w
) const {
   unsigned i,j,l,d;
   wbperm P; wbvector<unsigned> D;
   double n2;

   QIDX.blockSum(QDIM,Q); Q.groupRecs(P,D);
   w.init(D.len);

   for (l=i=0; i<D.len; ++i, l+=d) { d=D[i];
      for (n2=0, j=0; j<d; ++j) n2 += (DATA[P[l+j]]->norm2());
      w[i]=sqrt((double)n2);
   }
}


template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsum(
  const wbindex &J, wbMatrix<TQ> &Q
) const {

    unsigned i,j, r=rank(FL);
    const TQ *q0; TQ *q;

    if (J.isEmpty()) { Q.init(); return Q; }

    for (i=0; i<J.len; ++i) if (J[i]>=r) 
    wblog(FL,"ERR Block index out of range (%d/%d)",J[i],r);

    Q.init(QIDX.dim1,QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
       q=Q.rec(i); q0=QIDX.rec(i);
       for (j=0; j<J.len; ++j) addRange(q0+J[j]*QDIM, q, QDIM);
    }

    return Q;
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::FlipQ() {
   wbperm P(rank(FL)); P.Flip();
   PermuteQ(P);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::PermuteQ(const wbperm &P) {

   unsigned i,k,r=rank(FL), s=QDIM*sizeof(TQ);
   wbMatrix<TQ> Q0(QIDX);
   TQ *q0, *qp;

   if (!QDIM) wblog(FL,"WRN QDIM=0!");
   if (P.isEmpty()) return;
   if (isEmpty() || QIDX.dim2==0) return;

   if (!P.isValidPerm(r)) wblog(FL,
   "ERR %s() invalid permutation [%s] (%d)",FCT,P.toStr().data,r);

   for (i=0; i<QIDX.dim1; ++i) {
       q0=Q0.rec(i); qp=QIDX.rec(i);
       for (k=0; k<P.len; ++k) memcpy(qp+k*QDIM, q0+P[k]*QDIM, s);
   }
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::permuteQ(const wbperm &P, wbMatrix<TQ> &Q) const {

   if (&Q==&QIDX) {
      wblog(FL,"ERR You shall not permuteQ onto itself");
   }
   else {
      unsigned i, k, r=rank(FL), s=QDIM*sizeof(TQ);
      const TQ *q0; TQ *qp;

      if (P.isEmpty()) { Q=QIDX; return; }

      if (!P.isValidPerm(r))
      wblog(FL,"ERR invalid permutation [%s].", P.toStrD().data);

      Q.init(QIDX.dim1, QIDX.dim2);

      for (i=0; i<QIDX.dim1; ++i) {
          q0=QIDX.rec(i); qp=Q.rec(i);
          for (k=0; k<P.len; ++k) memcpy(qp+k*QDIM, q0+P[k]*QDIM, s);
      }
   }
}


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Permute(
   const wbperm &P_, char iflag
){
   unsigned i,n,r=rank(FL);
   if ((!r && isEmpty()) || !P_.len || P_.isIdentityPerm(FL,r)) {
      return *this;
   }
   if (P_.len!=r) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P_.len,r);

   wbperm P(P_,iflag);
   PermuteQ(P);

   for (i=0; i<DATA.len; ++i) DATA[i]->Permute(P);
   if (itags.len) itags.Select(P);

   for (n=CGR.numel(), i=0; i<n; ++i) {
      CGR[i].Permute(P);
   }

   return *this;
};


template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::permute(
   const wbperm &P_, QSpace<TQ,TD> &B, char iflag
 ) const {

   unsigned r=rank(FL);
   if (P_.isEmpty() || P_.isIdentityPerm(FL,r)) {
      B.init(*this); return B;
   }

   unsigned i,n;
   wbperm P(P_,iflag); B.init();

   permuteQ(P,B.QIDX);
   B.QDIM=QDIM; B.qtype=qtype; B.itags=itags;

   if (B.itags.len) B.itags.Select(P);

   B.setupDATA();
   for (n=DATA.len, i=0; i<n; ++i) {
      DATA[i]->permute(*B.DATA[i],P);
   }

   if (CGR.isEmpty()) { B.CGR.init(); }
   else {
      B.setupCGR(); n=CGR.numel();
      if (n!=B.CGR.numel()) wblog(FL,"ERR %s() got CGR size mismatch "
         "(%s <> %s)",FCT,CGR.sizeStr().data,B.CGR.sizeStr().data);
      for (i=0; i<n; ++i) {
         CGR[i].permute(P,B.CGR[i]);
      }
   }

   return B;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::setQ(unsigned k,
    const QSpace<TQ,TD> &B, unsigned k0
) {

    unsigned i, s=QDIM*sizeof(TQ), dd0=B.QIDX.dim2, dd2=QIDX.dim2;
    TQ *d0 = B.QIDX.data+k0*QDIM, *d2 = QIDX.data+k*QDIM;

    if (QIDX.dim1!=B.QIDX.dim1 || QDIM!=B.QDIM) wblog(FL,
       "ERR %s() incompatible objects (%d,%d; %d,%d)",
        FCT, QIDX.dim1, B.QIDX.dim1, QDIM, B.QDIM);

    if (rank(FL)<=k || B.rank(FL)<=k0) wblog(FL,
       "ERR %s() index out of bounds (%d/%d; %d/%d)",
        FCT, k, rank(FL), k0, B.rank(FL));

    for (i=0; i<QIDX.dim1; ++i, d0+=dd0, d2+=dd2)
    memcpy(d2,d0,s);
}


template <class TQ, class TD> inline
void QSpace<TQ,TD>::setQRec(
    unsigned r, unsigned j, const TQ *q0, const TQ *q2
){
    if (r>=QIDX.dim1 || QDIM*(j+1)>QIDX.dim2) wblog(FL,
       "ERR %s() index out of bounds (%d,%d; %d,%d)",
        FCT, r, j, QIDX.dim1, QDIM);

    TQ *q = QIDX.rec(r)+j*QDIM;
    if (q2==NULL)
         for (unsigned i=0; i<QDIM; ++i) q[i]=q0[i];
    else for (unsigned i=0; i<QDIM; ++i) q[i]=q0[i]+q2[i];
}


template <class TQ, class TD> inline
bool QSpace<TQ,TD>::hasEqualQ(unsigned i1, unsigned i2) const {

    unsigned i, s=QDIM*sizeof(TQ), D=QIDX.dim2, r=rank(FL);
    TQ *d1 = QIDX.data+i1*QDIM, *d2 = QIDX.data+i2*QDIM;

    if (i1>=r || i2>=r) wblog(FL,
    "ERR %s() index out of bounds (%d,%d; %d)",FCT,i1,i2,r);

    for (i=0; i<QIDX.dim1; ++i, d1+=D, d2+=D)
    if (memcmp(d1,d2,s)) return 0;

    return 1;
};


template <class TQ, class TD>
bool QSpace<TQ,TD>::hasQOverlap(const QSpace<TQ,TD> &A) const {

    wbindex Ia, Ib;
    if (qtype!=A.qtype) {
       wblog(FL,"WRN got objects with different qtype");
       return 0;
    }

    matchIndex(QIDX, A.QIDX, Ia, Ib);
    return (Ia.len!=0);
};


template <class TQ, class TD>
template <class TDB>
int QSpace<TQ,TD>::hasQOverlap(
   unsigned ica, const QSpace<TQ,TDB> &B, unsigned icb,
   const char type
) const {

   unsigned ra=rank(FL), rb=B.rank(FL); int i=0;
   wbMatrix<TQ> QA,QB;
   wbindex Ia,Ib;

   if (!ica || !icb) wblog(FL,
      "ERR %s() index must be 1-based (%d,%d)",FCT,ica,icb);
   ica--; icb--;

   if (QDIM!=B.QDIM) wblog(FL,
      "ERR %s() incompatible objects\nQDIM=%d,%d",
       FCT, QDIM, B.QDIM);
   if (ica>=ra || icb>=rb) wblog(FL,
      "ERR %s() index out of bounds\n%d/%d, %d/%d (QDIM=%d,%d)",
       FCT, ica,ra, icb,rb, QDIM, B.QDIM);
 
     getQsub(ica,QA).makeUnique();
   B.getQsub(icb,QB).makeUnique(); matchIndex(QA,QB,Ia,Ib);

   if (type=='<') i=(Ia.len==QA.dim1); else
   if (type=='>') i=(Ib.len==QB.dim1); else
   if (type=='=') i=(Ib.len==QB.dim1 && Ib.len==QB.dim1);
   else i=(Ia.len!=0);

   return i;
}


template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::QIDX_toMx() const {

   unsigned k, r;
   mxArray *a;

   if (QIDX.isEmpty())
   return mxCreateCellMatrix(0,0);

   if (!QDIM || QIDX.dim2%QDIM) wblog(FL,
      "ERR %s() invalid QSpace (QIDX: %dx%d; %d)",
       FCT,QIDX.dim1,QIDX.dim2,QDIM
   );

   r=QIDX.dim2/QDIM; a=mxCreateCellMatrix(1,r);
   if (!a) wblog(FL,"ERR could not allocate cell array (%d)", r);

   for (k=0; k<r; ++k)
   mxSetCell(a, k, getQsub(k).toMx());

   return a;
}

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::INFO_toMx() const {

   if ((qtype.isEmpty() && CGR.isEmpty() &&
       itags.isEmpty() && otype==QS_NONE)
   ){ return mxCreateCellMatrix(0,0); }

   const char *fields[] = { "qtype","otype","itags","cgr" };

   mxArray *S=mxCreateStructMatrix(1,1,4,fields);

   mxSetFieldByNumber(S,0,0,qStr().toMx());
   mxSetFieldByNumber(S,0,1,otype2Str().toMx());
   mxSetFieldByNumber(S,0,2,itags.toMx());
   mxSetFieldByNumber(S,0,3,CGR.toMx());

   return S;
}

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::DATA_toMx() const {

   mxArray *a=mxCreateCellMatrix(DATA.len, DATA.len ? 1 : 0);
   if (!a) wblog(FL,"ERR failed to allocate cell array (%d)",DATA.len);

   for (unsigned k=0; k<DATA.len; ++k)
   mxSetCell(a, k, DATA[k]->toMx());

   return a;
}

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::DATA_save2Mx(char vflag) {

   mxArray *a=mxCreateCellMatrix(DATA.len, DATA.len ? 1 : 0);
   if (!a) wblog(FL,"ERR failed to allocate cell array (%d)",DATA.len);

   if (vflag && getDataSize()<(1<<30)) vflag=0;
   if (vflag) wblog(FL,"TST %s() %d entries ...",FCT,DATA.len); 

   for (unsigned k=0; k<DATA.len; ++k) {
      mxSetCell(a, k, DATA[k]->toMx());
      if (vflag && k%10==0) printf("  %6d/%d ...  \r",k+1,DATA.len);
      DATA[k]->init();
   }

   return a;
}


#ifdef WBC_QSPACE_IO
   WbClock wbc_qs_toMx("QS::toMx");
#endif

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[] = { "Q","data","info" };
   return  mxCreateStructMatrix(m,n,3,fields);
}

template <class TQ, class TD>
void QSpace<TQ,TD>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   mxArray *a;

#ifdef WBC_QSPACE_IO
   wbc_qs_toMx.resume();
#endif

   if (tst) {
      unsigned s=0; int q=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (q=mxGetFieldNumber(S,"qtype"))<0) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,q);
   }

   a=QIDX_toMx(); mxSetFieldByNumber(S,i,0,a);
   a=DATA_toMx(); mxSetFieldByNumber(S,i,1,a);
   a=INFO_toMx(); mxSetFieldByNumber(S,i,2,a);

#ifdef WBC_QSPACE_IO
   wbc_qs_toMx.stop();
#endif
};


template <class TQ, class TD>
void QSpace<TQ,TD>::save2MxStruct(
   mxArray *S, unsigned i, char tst, char vflag
){
   mxArray *a;

#ifdef WBC_QSPACE_IO
   wbc_qs_toMx.resume();
#endif

   if (tst) {
      unsigned s=0; int q=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (q=mxGetFieldNumber(S,"qtype"))<0) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,q);
   }

   a=QIDX_toMx();         mxSetFieldByNumber(S,i,0,a);
   a=DATA_save2Mx(vflag); mxSetFieldByNumber(S,i,1,a);
   a=INFO_toMx();         mxSetFieldByNumber(S,i,2,a);

   init();

#ifdef WBC_QSPACE_IO
   wbc_qs_toMx.stop();
#endif
}


template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::toMx() const {

   mxArray *S=mxCreateStruct(1,1);
   this->add2MxStruct(S,0);
   return S;
}

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::save2Mx(char vflag) {

   mxArray *S=mxCreateStruct(1,1);
   this->save2MxStruct(S,0,0,vflag);
   return S;
}


#ifdef WBC_QSPACE_IO
   WbClock wbc_qs_mxInit("QS::mxInit");
#endif

template <class TQ, class TD>
void QSpace<TQ,TD>::init(
   const char *F, int L,
   const mxArray *S,
   const char ref,
   unsigned k 
){
   unsigned i,j,l,m,n,r,dim1=0,dim2=0; unsigned rk=-1;
   mxArray *aq,*ad,*ai;
   wbvector< wbMatrix<TQ>  > MQ;
   wbvector< wbMatrix<TQ>* > mq;
   const char cflag = (typeid(TD)==typeid(wbcomplex) ? 'C' : 0);
   char refD=0, refC=1;

#ifdef WBC_QSPACE_IO
   wbc_qs_mxInit.resume();
#endif

   clearQSpace();

   if (!S || mxIsEmpty(S) || mxIsEmptyQSpace(S,k)) return;

   aq=mxGetField(S,k,"Q");
   ad=mxGetField(S,k,"data");
   ai=mxGetField(S,k,"info");

   if (!mxIsQSpace(F_L,S,rk,cflag,k) || !aq || !ad) wblog(F_L,
      "ERR got invalid QSpace (%d,%d)",aq==0,ad==0);

   if (ref) {
      if (ref=='r') { refD        = ref; } else
      wblog(FL,"ERR %s() got invalid ref=%c<%d>",FCT,ref,ref);

      if (WbUtil<TD>::isComplex()) refD=0;
   }

   m=mxGetM(aq); r=mxGetN(aq);
   
   if ((m!=1 && r!=1 && m && r) || mxGetNumberOfDimensions(aq)>2) {
       wblog(F,L,"ERR %s() invalid QSpace\n"
      "Q must be blocked into row cell vector (%dx%d; %d)",
       FCT,m,r,mxGetNumberOfDimensions(aq));
   }

   r*=m; MQ.init(r);

   for (i=0; i<r; ++i) {
      MQ[i].init(mxGetCell(aq,i));

      if (i==0) {
         dim1=MQ[i].dim1; dim2=MQ[i].dim2;
      }
      else if (MQ[i].dim1!=dim1 || MQ[i].dim2!=dim2) { wblog(F,L,
        "ERR %s() invalid QSpace\nsize mismatch in Q data (%dx%d; %dx%d)",
         FCT, MQ[i].dim1, MQ[i].dim2, dim1, dim2);
      }
   }

   mq.init(MQ.len); for (i=0; i<MQ.len; ++i) mq[i] = &MQ[i];

   QIDX.CAT(2, (const wbMatrix<TQ>**) mq.data, mq.len);
   QDIM=dim2;

   setupDATA();

   if (ai && !mxIsEmpty(ai)) {
      if (mxIsStruct(ai)) {
         mxArray *a; QVec qv; unsigned dq=0;

         a=mxGetField(ai,0,"qtype");
         if (a) qtype.init(F,L,a); else qtype.init();

         
         initOType(F_L,mxGetField(ai,0,"otype"));
         
         itags.init(F_L,mxGetField(ai,0,"itags"));

         if (itags.len && itags.len!=r) wblog(FL,
            "ERR %s() invalid number of itags (%s; %d/%d)",
             FCT,IT2STR__,itags.len,r);

         a=mxGetField(ai,0,"cgr");
         if (a && !mxIsEmpty(a)) {

             wbvector<unsigned> qdc; qtype.Qpos(qdc);
             QSet<TQ> Q;

             setupCGR();
             m=mxGetM(a); n=mxGetN(a);

             if (m==1 && !QIDX.dim1) {
                if (!mxIsScalarQSpace(S) || m>1 || r) wblog(FL,
                   "ERR %s() invalid scalar QSpace",FCT);
                QIDX.init(1,0); QDIM=qtype.Qlen(); r=2;
                setupDATA(); setupCGR();
             }
             else
             if (!QDIM || QIDX.dim2%QDIM || QIDX.dim2/QDIM!=r) wblog(FL,
                "ERR %s() empty QDIM while info.cgr is specified (%d/%d)",
                 FCT,m,n);
             else 
             if (mxGetNumberOfDimensions(ad)>2 || 
                 m!=QIDX.dim1 || m!=CGR.dim1 ||
                 n!=qtype.len || n!=CGR.dim2 || qtype.Qlen()!=QDIM)
             wblog(FL,"ERR %s() invalid dimensions for cell array "
                "info.cgr\n{%s} cgr: %dx%d / %dx%d; QIDX: %dx%d/%d",
                 FCT, qStr().data, m,n, CGR.dim1, CGR.dim2,
                 QIDX.dim1, QIDX.dim2, QDIM
             );

             if (!mxIsCell(a)&& !mxIsStruct(a)) wblog(FL,
                "ERR %s() cell array expected for info.cgr field (%s)",
                 FCT,mxGetClassName(a)
             );

             for (l=j=0; j<n; ++j)
             for (  i=0; i<m; ++i, ++l) {
                 Q.init(qtype[j], QIDX.rec(i)+qdc[j], r, QDIM, itags);
                 CGR(i,j).init(F_L,a,l, refC, &Q);
             }
         }

         if (qtype.len) { dq=qtype.Qlen();
            if (dq!=QDIM && (QIDX.dim2 || CGR.dim2)) wblog(FL,"ERR "
               "incompatible number of symmetry labels (%d/%d; %s; %s)",
               QDIM, dq, QIDX.sizeStr().data, CGR.sizeStr().data
            );
         }
      }
      else wblog(F,L,"ERR invalid qtype (structure required)");

      if (qtype.allAbelian()) {
         if (!CGR.isEmpty()) {
            for (i=0; i<CGR.dim1; ++i)
            for (j=0; j<CGR.dim2; ++j) { if (CGR(i,j).cgr) {
               const cdata__ &c=(CGR(i,j).cgr->cgd);
               if (!c.isScalar()) wblog(FL,
                  "ERR %s() invalid abelian CGC space (%s)",
                   FCT,c.sizeStr().data);
               if (c.D.data && c.D.data[0]!=1) wblog(FL,
                  "ERR %s() invalid abelian CGC coefficient (%.4g)",
                   FCT,double(c.D.data[0])
               );
            }}
            CGR.init();"TST %s()",FCT); 
         }

         qtype.init();
      }
   }
   else qtype.init();


   m=mxGetM(ad); n=mxGetN(ad);

   if ((m>1 && n>1) || mxGetNumberOfDimensions(ad)>2) wblog(F,L,
      "ERR QSpace() cell VECTOR expected for data (%d,%d)",m,n);

   m*=n;

   if (DATA.len!=m) {
      if (m==1 && DATA.len==0 && QIDX.isEmpty() && CGR.isEmpty()) {
         QIDX.init(1,0); setupDATA();
         if (r<2) r=2;
      }
      else wblog(F,L,
        "ERR QSpace() dimension mismatch (QIDX: %dx%d, data: %d)",
         QIDX.dim1, QIDX.dim2, m
      );
   }
   
   if (r==1) r=2;

   for (i=0; i<m; ++i) {
      DATA[i]->init(mxGetCell(ad,i), refD);
      DATA[i]->addSingletons(r);
   }


#ifdef WBC_QSPACE_IO
   wbc_qs_mxInit.stop();
#endif
};



template <class TQ, class TD> inline
void QSpace<TQ,TD>::init2ref(const QSpace<TQ,TD> &A) {

   clearQSpace();
    
   QIDX.init2ref(A.QIDX); QDIM=A.QDIM;
   DATA.init2ref(A.DATA);
   qtype.init2ref(A.qtype); otype=A.otype;
   itags.init2ref(A.itags);

   if (!A.CGR.isEmpty())  CGR.init2ref(A.CGR);

   isref=1;
};


template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::NormCG(char skipzeros) {
   
   if (CGR.isEmpty()) return *this;

   wbvector<unsigned> mark;
   unsigned i,j; RTD cfac, one=1;

   gotCGS(FL);
   unRef();

   for (i=0; i<CGR.dim1; ++i) {
      for (cfac=one, j=0; j<CGR.dim2; ++j) {
         cfac *= CGR(i,j).NormSignR();
      }
      if (cfac!=1.) {
         if (cfac!=0) { (*DATA[i])*=TD(cfac); }
         else if (skipzeros) {
            if (!mark.len) { mark.init(CGR.dim1); }
            mark[i]=1;
         }
      }
   }

   if (mark.len) {
      wbindex I;
      Select(mark.findzeros(I));
   }

   return *this;
};


template <class TQ, class TD> inline
void QSpace<TQ,TD>::setRand(TD nrm) {

   unsigned i;
   for (i=0; i<DATA.len; ++i) DATA[i]->setRand();

   if (nrm>0) {
      double fac=double(nrm)/sqrt(double(norm2()));
      for (i=0; i<DATA.len; ++i) (*DATA[i]) *= fac;
   }
}


template <class TQ, class TD>
void QSpace<TQ,TD>::info(const char *vname) const {

   wbvector<INDEX_T> D,DX;
   char cgflag=gotCGS();
   unsigned r=(QDIM!=0 ? QIDX.dim2/QDIM : 0);
   size_t l=0, n=64; char tstr[n];

   if (QDIM && QIDX.dim2 && QIDX.dim2%QDIM) wblog(FL,
      "ERR %s() %d/%d = ?",FCT, QIDX.dim2, QDIM);

   l=snprintf(tstr,n,"QSpace<%s,%s> `%s'",
      getName(typeid(TQ)).data, getName(typeid(TD)).data,
      qtype.len ? qStr().data : "");
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   getDim(D,&DX);

   printf("  %-8s %-22s %3dx{%dx%d}  %10s",
      vname,tstr, r,QIDX.dim1, QDIM, D.toStrD().data);
   if (cgflag>0) printf(" [%12s]", DX.toStr().data);
   printf("%s\n", isref? "  *REF*":"");
};


template <class TQ, class TD>
void QSpace<TQ,TD>::print(const char *vname, char dflag) const {

   size_t l=0, n=32;
   char tstr[n], fmt[]="%2g", cgflag=gotCGS(FL);

   if (typeid(TQ)!=typeid(double)) fmt[2]='d';

   l=snprintf(tstr,n,"<%s,%s>",
      getName(typeid(TQ)).data, getName(typeid(TD)).data);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   printf("\n  %-9s %g-D QSpace %-20s 0x%lX  %dx(%g*%d) %s\n\n",
      vname, QIDX.dim2 ? QIDX.dim2/double(QDIM) : 0, tstr,
      (unsigned long)DATA.data, QIDX.dim1, QIDX.dim2/double(QDIM),
      QDIM, isref ? " *REF*" : ""
   );


   for (unsigned j,s, n=QIDX.dim1, i=0; i<n; ++i) {
      printf("%4d: { %s } %8s ", i+1,
      i<QIDX.dim1 ? QIDX.rec2Str(i, fmt, " ", QDIM, " ;").data : "???",
      i<DATA.len  ? DATA[i]->sizeStr().data : "???");

      if (i>=DATA.len) {
         printf("   ???\n");
         continue;
      }

      if (dflag) {
         if (!(DATA[i]->printdata(vname))) printf("\n");
      }
      else {
         if (cgflag>0) { printf("|");
            for (j=0; j<CGR.dim2; ++j)
            printf("%8s",CGR(i,j)->sizeStr().data); printf(" |");
         }

         s=DATA[i]->SIZE.prod();

         if (s==0) printf("   []\n"); else
         if (s==1 && DATA[i]->SIZE.len==1)
            printf("   %s\n", num2Str(DATA[i]->data[0]).data);
         else {
            s *= sizeof(TD);
            if (s<(1<<10)) printf("   %5d B\n", s); else
            if (s<(1<<20)) printf("   %5.2f kB\n", s/double(1<<10));
            else           printf("   %5.2f MB\n", s/double(1<<20));
         }
      }
   }
   printf("\n");
};


#endif

