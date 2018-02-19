#ifndef __WB_NRGDATA_HCC__
#define __WB_NRGDATA_HCC__

template <class TQ, class TD>
class NRGData;

#define NRG_FLEN 128

unsigned NRG_N, NRG_ITER;
char NRG_FILE[NRG_FLEN];

unsigned STRICT_ITER0=1;

WbClock tic_nrgIO("NRGData I/O");



template <class TQ, class TD>
class NRGIndex {

  public:

    NRGIndex() : idx(-1), nrg(NULL) {};
    NRGIndex(const NRGIndex &I) : idx(I.idx), nrg(I.nrg) {};


    void init() { idx=-1; nrg=NULL; };
    void init(int i, NRGData<TQ,TD> *d) { idx=i; nrg=d; }
    void init(int i, const char *t)     { idx=i; tag=t; }

    mxArray* toMx() const;
    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    void put (const char *vname, const char *ws="base") const {
       mxArray *S=toMx();
       mxPutAndDestroy(FL,S,vname,ws);
    };

    int idx;
    NRGData<TQ,TD> *nrg;
    wbstring tag;

  private:

};


template <class TQ, class TD>
mxArray* NRGIndex<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[] = {"idx","nrg","tag"};
   return  mxCreateStructMatrix(m,n,3,fields);
};

template <class TQ, class TD>
void NRGIndex<TQ,TD>::add2MxStruct(
   mxArray *S, unsigned i, char tst __attribute__ ((unused))
 ) const {

   mxSetFieldByNumber(S,i,0, numtoMx(idx));
   mxSetFieldByNumber(S,i,1, numtoMx((unsigned long)(nrg)));
   mxSetFieldByNumber(S,i,2, mxCreateString(tag.isEmpty() ? tag.data:""));
};


namespace Wb {

   template <class TQ, class TD>
   int updateOp(const char *F, int L,
      const NRGData<TQ,TD> &NRG,
      wbvector< QSpace<TQ,TD> > &F12,
      const QSpace<TQ,TD> &A1,
      const QSpace<TQ,TD> &A2,
      const wbvector< QSpace<TQ,TD> > &FKK,
      unsigned iter 
   );

   template<class TT>
   void initLocalOp(const char *F, int L,
      const mxArray* a, TT &C,
      int k=0, char force=0, char *name=NULL
   );

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbvector< QSpace<TQ,TD> > &C
   ){ mxInitQSpaceVec(F,L,a,C); }

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbMatrix< QSpace<TQ,TD> > &C
   ){ mxInitQSpaceMat(F,L,a,C); }

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbvector< wbvector< QSpace<TQ,TD> > > &C
   ){ mxInitQSpaceVecVec(F,L,a,C); }
};


template <class TQ, class TD>
class NRGData {
  public:

    NRGData(const char *s=NULL)
     : name(s), calc(0), store(0), kloc(0), MXloc(0), MX(NULL), aux(NULL) {};

    NRGData(const NRGData &A)
     : name(A.name), nrgIdx(A.nrgIdx), calc(0), store(0), kloc(0),
       MXloc(0), MX(NULL), aux(NULL)
    {  wblog(FL,"WRN Rule #101 - you shall not copy NRGData!");
    };

   ~NRGData() {
       if (MX && !NAME.isEmpty() && !MXloc) wblog(FL,
          "WRN ** data return `%s' not implemented yet **\n(%lX; %d)",
           NAME.data, MX, MXloc);
       init();
    };

    void initName(const char *n) { name=n; };
    void init(const char *F, int L, const char *XX, int iter=-1);
    void init(const char *F, int L, const mxArray* a,
       wbvector< QSpace<TQ,TD> > &C0,
       int k=0, char force=0, const char *vname0=0
    );

    void init(char flag=0) {
       if (!flag) { name.init(); NAME.init(); }

       nrgIdx.init(); calc=0; store=0;
       if (aux) {
          mxDestroyArray(aux);
          aux=NULL;
       }
       if (MX && MXloc) {
          mxDestroyArray(MX);
          MXloc=0; MX=NULL;
       }

       KK.init(); TK.init(); K.clearQSpace();
       KT.init(); TT.init(); T.clearQSpace(); A.clearQSpace();
    };

    bool isQSpaceVec(const char *nm, int iter=-1);

    bool checkVec(const char *F, int L,
       const char *XX, unsigned R, const char *istr="");

    bool checkVecVec(const char *XX, unsigned r);

    void forceCalc();
    void getRelevantOps(wbindex &I) const;

    wbvector< QSpace<TQ,TD> >& getVec(const char *t) {
       if (t && strlen(t)==2)
       switch (t[0]) {
          case 'K': return (t[1]=='K' ? KK : KT);
          case 'T': return (t[1]=='K' ? TK : TT);
       }
       wblog(FL,"ERR invalid tag `%s'", t ? t : "(null)");
       return KK;
    };

    QSpace<TQ,TD>& getQSpace(char t1, char t2, unsigned i=0) {
       char t[3]={t1,t2,0};
       return getQSpace(t,i);
    };

    QSpace<TQ,TD>& getQSpace(const char *t, unsigned i=0);
    QSpace<TQ,TD>& getQSpace(char t) {
       switch (t) {
          case 'A': return A;
          case 'K': return K;
          case 'T': return T;
       }
       wblog(FL,"ERR invalid tag `%c'(%d)",t,t);
       return A;
    };

    void info(const char *vstr=NULL) const;
    void dispIdx(const char *F, int L) const;
    void info(const char *F, int L, const char *istr=" * ") const;
    void disp_otype(const char *F=0, int L=0) const;

    const mxArray* getMxInfo(const char *n);
    mxArray* getMxData(
    const char *nm, int iter=-1, char ionly=0);

    mxArray* toMx() const;

    void put(const char *vname, const char* ws="base") {
       mxPutAndDestroy(FL, toMx(), vname);
    };

    void put(const char *F, int L, const char *vname, const char* ws="base"
    ){ wblog(F,L,"I/O putting '%s' to %s", vname, ws);
       mxPutAndDestroy(FL, toMx(), vname);
    };

    template<class T>
    int getPara(const char *name, T &x, int iter=-1, char i=0) {
       mxArray *a=getMxData(iter);
       return mxGetNumber(a,x);
    }

    bool opExists(char tflag=1);

    void initOp(unsigned n=1, char _calc=1, char _store=0) {
       KK.initDef(n); TK.initDef(n); nrgIdx.initDef(n);
       KT.initDef(n); TT.initDef(n); calc=_calc; store=_store;
    };

    bool initOp(const char *F, int L,
       const wbvector< QSpace<TQ,TD> > &C0,
       char tflag=1,
       NRGData<TQ,TD>* CR=NULL,
       const char *istr=""
    );

    bool initOp(const char *F, int L,
       const wbvector< wbvector< QSpace<TQ,TD> > > &CC,
       NRGData<TQ,TD>* CR=NULL,
       const char *istr="",
       char tflag=1 
    );

    void initSym(const char *tag);

    void updatePara(
       const char *F, int L,
       const char *tag, mxArray *, int iter=-1
    ) const;

    void updatePara(const char *F, int L,
       const char *tag, int iter=-1,
       const char *tag2=NULL, char all=0
    ) const;

    void updateOp(const char *F, int L,
       const NRGData<TQ,TD> &A,
       const NRGData<TQ,TD> &B, unsigned nop, int iter=-1
    );

    void updateOp(const char *F, int L,
       const NRGData<TQ,TD> &A, unsigned nop, int iter=-1
    ) { updateOp(F,L,A,A,nop,iter); };

    void unsetAllRefs();
    unsigned skipOpZeros(double eps=DEPS, char bflag='b');

    void storeOp(const char*, int, int iter=-1, char flag=4) const;

    void times(TD fac, const char *flag);

    void setupIO(const char*, int,
    const char*, const mxArray*, char logf=0);

    void updateInfo(
       const char *F, int L,
       const char*, mxArray *, char keep=-1
    ) const;

    void updateInfO(
       const char *F, int L,
       const char *vtag, mxArray *, char keep=-1
    ) const;

    const char* getFileName(unsigned k, char lenient=0) const {
       if (!lenient && k>=NRG_N) wblog(FL,
       "ERR %s - index out of bounds (%d/%d)", __FUNCTION__, k, NRG_N);
       return getFileName_aux( k);
    };


    wbstring name;
    wbvector< NRGIndex<TQ,TD> > nrgIdx;
    char calc, store;
    int  kloc;
    char MXloc;

    wbvector< wbvector< QSpace<TQ,TD> > > CI;
    wbvector< QSpace<TQ,TD> > KK, KT, TK, TT;
    QSpace<TQ,TD> A, K, T;

    wbvector<char> issym;


    mxArray* MX;
    wbstring NAME;

    mxArray *aux;

  private:


    const char* getFileNameI() const {
       return getFileName_aux(-1);
    };

    const char* getFileName_aux(int k) const {
       int i;
       if (k>=0)
            i=snprintf(NRG_FILE,NRG_FLEN,"%s_%02d.mat", NAME.data, k);
       else i=snprintf(NRG_FILE,NRG_FLEN,"%s_info.mat", NAME.data);
       if (i>=NRG_FLEN) wblog(FL,"ERR File name exceeds length (%d/127)",i);
       return NRG_FILE;
    }
};



template<class TQ, class TD>
void NRGData<TQ,TD>::setupIO(
   const char *F, int L,
   const char *s, const mxArray *a, char logf
){
   
   if (s==NULL || s[0]==0) {
      if (!a || (!mxIsCell(a) && !mxIsStruct(a))) wblog(F,L,
      "ERR invalid input NRG data (%s)", a ? mxGetClassName(a):"NULL");

      NAME.init(); NRG_N=0; NRG_ITER=0;
      MX=(mxArray*)a;
   }
   else {
      Wb::matFile f; unsigned n=strlen(s);
      tic_nrgIO.resume();

      if (!n || n+16>NRG_FLEN) wblog(FL,
      "ERR invalid file/variable name `%s' (%d/%d)", s, n, NRG_FLEN);

      NRG_N=0; NRG_ITER=0; MX=NULL;
      NAME=s;

#ifdef MATLAB_MEX_FILE
      mxArray *a=mexGetVariable("caller",s);
      if (a && (mxIsStruct(a) || mxIsCell(a))) {
         if (logf) wblog(F,L,"<i> I/O structure variable: %s", s);

         MX=a;
         return;
      }
#endif

      if (!f.open(F,L, getFileName(0,'l'),"r"))
           wblog(F,L,"ERR cannot open files of type %s_##.mat", s);
      else f.close();

      tic_nrgIO.stop();

      if (logf) {
         if (strlen(s)<30) {
            wblog(F,L,"<i> using NRG data: %s_##.mat", s);
         }
         else {
            if (getcwd(str,256)==NULL) wblog(FL,"ERR failed read cwd");
            wblog(F,L,"<i> using NRG data ...");
            printf("\n   pwd: %s\n   I/O: %s_##.mat\n\n",
            Wb::repHome(str).data, Wb::repHome(s).data);
         }
      }
   }
};


template<class TQ, class TD>
void NRGData<TQ,TD>::init(const char *F, int L, const char *XX, int iter) {

   unsigned n=strlen(XX);
   mxArray *a;
   
   if (!n || n>2) wblog(FL,"ERR invalid operator tag `%s'", XX);
   if (iter<0) iter=NRG_ITER;

   strcpy(str,name.data); if (XX[0]!='A') strcat(str,XX);
   a=getMxData(str,iter);

   if (!strcmp(XX,"KK")) mxInitQSpaceVec(F,L,a,KK); else
   if (!strcmp(XX,"KT")) mxInitQSpaceVec(F,L,a,KT); else
   if (!strcmp(XX,"TK")) mxInitQSpaceVec(F,L,a,TK); else
   if (!strcmp(XX,"TT")) mxInitQSpaceVec(F,L,a,TT); else
   if (!strcmp(XX,"K" )) K.init(F,L,a); else
   if (!strcmp(XX,"T" )) T.init(F,L,a); else
   if (!strcmp(XX,"A" )) A.init(F,L,a); else
   wblog(FL,"ERR invalid operator tag '%s'", XX);
};



template<class TQ, class TD>
void NRGData<TQ,TD>::init(
   const char *F, int L, const mxArray* a,
   wbvector< QSpace<TQ,TD> > &C0,
   int k,
   char force,
   const char *vname0
){
   char nstr[32]; nstr[0]=0;

   Wb::initLocalOp(F,L,a,C0,k,force,nstr);

   if (vname0 && vname0[0]) { name=vname0; }
   else if (nstr[0]) { name=nstr; }

}


template<class TT>
void Wb::initLocalOp(const char *F, int L,
   const mxArray* a,
   TT &C,
   int k,
   char force,
   char *name
){
   size_t n=12; char istr[n];
   if (k)
        snprintf(istr,n,"arg #%d",k);
   else strcpy(istr,"operator");

   if (!a || mxIsEmpty(a)) {
      if (force) wblog(F,L,"ERR invalid %s%s",istr,a ? " (empty)":"");
      return;
   }

   if (mxIsChar(a)) {
      const unsigned n=32; wbstring aux;

      if (name==NULL) { aux.init(n+1); name=aux.data; }
      name[0]=0;

      if (mxGetString(a,name,n) || name[0]==0) wblog(F,L,
         "ERR invalid string of %s%s",istr,name[0] ? "":" (empty)");

      a=mexGetVariablePtr("caller",name);

      if (!a) wblog(F,L,
      "ERR Can not read QSpace `%s' (%s) from workspace",name,istr);
   }
    
   Wb::initLocalOp_aux(F,L,a,C);

   if (force && C.isEmpty()) wblog(F,L,
   "ERR invalid QSpace (%s is empty)",istr);
};


template<class TQ, class TD>
void NRGData<TQ,TD>::forceCalc() {

   wbvector< NRGIndex<TQ,TD> > &I=nrgIdx;
   if (I.len==0) return;

   for (unsigned i=0; i<I.len; ++i) {
      if (I[i].idx<0) continue;
      if (I[i].nrg==NULL) I[i].idx=-1;
      if (!I[i].tag.isEmpty()) {
         I[i].init(); I[i].idx=-1;
      }
   }

   calc=store=1;
}


template<class TQ, class TD>
void NRGData<TQ,TD>::updateOp(
   const char *F, int L,
   const NRGData<TQ,TD> &A,
   const NRGData<TQ,TD> &B,
   unsigned nop,
   int iter
){
   unsigned i,j,ndup, read, read2, m=nrgIdx.len;
   const char* tags[4] = { "KK","KT","TK","TT" };
   int e=0,k;

   wbvector< QSpace<TQ,TD> > X, X2;

   if (iter<0) iter=NRG_ITER;
   if (nrgIdx.len==0) return;

   for (ndup=i=0; i<m; ++i) if (nrgIdx[i].idx>=0) ++ndup;


   unsetAllRefs();

   if (ndup<nrgIdx.len) {
   if (nop==4) {
         Wb::updateOp(F,L,*this, TT, A.T, B.T, KK, iter);
         Wb::updateOp(F,L,*this, TK, A.T, B.K, KK, iter);
         Wb::updateOp(F,L,*this, KT, A.K, B.T, KK, iter);
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter);
   }
   else if (nop==3) {
         Wb::updateOp(F,L,*this, TT, A.T, B.T, KK, iter);
         Wb::updateOp(F,L,*this, TK, A.T, B.K, KK, iter);
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter);

      if (KT.len!=TK.len) KT.initDef(TK.len);
      for (i=0; i<TK.len; ++i) TK[i].transp(KT[i]);
   }
   else if (nop==1) {
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter);
   }
   else {
      wblog(FL,"ERR %s() invalid nop (%d)", FCT, nop);
   }}

   if (e) wblog(FL,
      "ERR %s() operator %s[%d] updates to empty QSpace!",
       FCT,name.data,iter
   );

   if (ndup) {
   for (i=0; i<nop; ++i) {

      wbvector< QSpace<TQ,TD> > &C = getVec(tags[i]);
      if (C.len!=m) C.initDef(m);

      for (read=read2=j=0; j<m; ++j) {
         k=nrgIdx[j].idx; if (k<0) continue;

         if (nrgIdx[j].nrg) {
            wbvector< QSpace<TQ,TD> > &R=nrgIdx[j].nrg->getVec(tags[i]);

            if (k>=(int)R.len) wblog(FL,
            "ERR Index out of bounds (%d/%d)",k,R.len);

            C[j].init2ref(R[k]);
         }
         else if (nrgIdx[j].tag.isEmpty()) {
            if (!read) { read=1;
               sprintf(str,"%s%s", name.data, tags[i]);
               mxInitQSpaceVecR23(F,L,getMxData(str,iter),X);
            }
            if (k>=(int)X.len) wblog(FL,
            "ERR Index out of bounds (%d/%d; %d)",k,X.len, nrgIdx.len);

            X[k].save2(C[j]);
         }
         else {
            if (!read2) { read2=1;
               sprintf(str,"%s%s", nrgIdx[j].tag.data, tags[i]);
               mxInitQSpaceVecR23(F,L, getMxData(str,iter),X2);
            }
            if (k>=(int)X2.len) wblog(FL,
            "ERR Index out of bounds (%d/%d)",k,X2.len);

            X2[k].save2(C[j]);
         }
      }
   }}
}


template<class TQ, class TD>
void NRGData<TQ,TD>::unsetAllRefs() {

   unsigned i,j,m=0;
   const char* tags[4] = { "KK","KT","TK","TT" };

   for (i=0; i<4; ++i) {
      wbvector< QSpace<TQ,TD> > &C = getVec(tags[i]);

      for (j=0; j<C.len; ++j) {
         if (nrgIdx[j].idx<0 || !C[j].isref) continue;
         C[j].init(); ++m;
      }
   }

}


template<class TQ, class TD>
unsigned NRGData<TQ,TD>::skipOpZeros(double eps, char bflag) {

   wbvector< QSpace<TQ,TD> > *XX[] = { &KK, &KT, &TK, &TT };
   unsigned i,j,n,r=0;

   for (i=0; i<4; ++i) {
      n=XX[i]->len;
      for (j=0; j<n; ++j)
      r += (*XX[i])[j].SkipZeroData(eps,bflag);
   }

   return r;
}


template<class TQ, class TD>
void NRGData<TQ,TD>::times(TD fac, const char *flag) {

   unsigned i,j,n=(flag ? strlen(flag) : 0);

   if (n==2) {
      if (!strcmp(flag,"XX")) {
         char* tags[4] = { "KK","KT","TK","TT" };
         for (i=0; i<4; ++i) {
            wbvector< QSpace<TQ,TD> > &X=getVec(tags[i]);
            for (j=0; j<X.len; ++j) X[j]*=fac;
         }
      }
      else {
         wbvector< QSpace<TQ,TD> > &X=getVec(flag);
         for (j=0; j<X.len; ++j) X[j]*=fac;
      }
   }
   else if (n==1) {
      if (flag[0]=='X') {
         getQSpace('K')*=fac;
         getQSpace('T')*=fac;
      }
      else {
         getQSpace(flag[0])*=fac;
      }
   }
   else wblog(FL, "ERR invalid flag `%s'", flag ? flag : "(null)");
}


template<class TQ, class TD>
void NRGData<TQ,TD>::storeOp(
  const char *F, int L, int iter, char flag
) const {

   unsigned i,ndup, m=nrgIdx.len;
   const char* tags[4] = { "KK","KT","TK","TT" };

   if (iter<0) iter=NRG_ITER;
   if (flag>4) wblog(FL,"ERR flag out of bounds (%d/%d)", flag, 4);

   for (ndup=i=0; i<m; ++i) if (nrgIdx[i].idx>=0) ++ndup;

   if (store && ndup<nrgIdx.len) {
      for (i=0; (char)i<flag; ++i)
      updatePara(F,L, tags[i], iter);
   }
} 


template<class TQ, class TD>
void NRGData<TQ,TD>::updatePara(const char *F, int L,
  const char *tag, mxArray *a, int iter
) const {
   sprintf(str,"%s%s", name.data, tag);
   if (!str[0]) wblog(F,L,"ERR invalid name (empty)");

   if (iter<0) iter=NRG_ITER;
   if (MX) {
      if (!NAME.isEmpty())
      mxUpdateField(F,L, MX,str,iter,a);
   }
   else {
      tic_nrgIO.resume();
      matPutVariable(F,L,getFileName(iter),str,a);
      tic_nrgIO.stop();
   }
}


template<class TQ, class TD>
void NRGData<TQ,TD>::updatePara(const char *F, int L,
   const char *tag,
   int iter,
   const char *tag2,
   char all
) const {

   const wbvector< QSpace<TQ,TD> > *XX=0;
   unsigned n=strlen(tag);
   const char es[]="";
   mxArray *a=0;

   if (n==1) {
      switch(tag[0]) {
         case 'A': a=A.toMx(); if (!tag2) tag2=es; break;
         case 'K': a=K.toMx(); break;
         case 'T': a=T.toMx(); break;
         default : wblog(FL,"ERR invalid tag `%s'", tag);
      }
   }
   else
   if (n==2) {
      switch(tag[0]) {
         case 'K': XX = (tag[1]=='K' ? &KK : &KT); break;
         case 'T': XX = (tag[1]=='K' ? &TK : &TT); break;
         default : wblog(FL,"ERR invalid tag `%s'", tag);
      }

      if (all) a=XX->toMx();
      else {
         wbindex I; getRelevantOps(I);
         a=QSpaceVec2Mx(*XX,I);
      }
   }
   else wblog(FL,"ERR invalid tag=`%s'", tag);

   sprintf(str,"%s%s", name.data, tag2 ? tag2 : tag);
   if (!str[0] || !name.isName()) wblog(F,L,"ERR invalid name `%s'",str);

   if (iter<0) iter=NRG_ITER;
   if (MX) {
      if (!NAME.isEmpty())
         mxUpdateField(F,L, MX,str,iter,a);
      else wblog(FL,
        "ERR nowhere to store RHO! Hint: specify NRG-data\n"
        "by name to workspace variable or set locRho."
      );
   }
   else {
      tic_nrgIO.resume();
      matPutVariable(F,L,getFileName(iter),str,a);
      tic_nrgIO.stop();
   }
}


template<class TQ, class TD>
inline void NRGData<TQ,TD>::getRelevantOps(wbindex &I) const {

   unsigned i,k;
   I.init(nrgIdx.len);

   for (k=i=0; i<nrgIdx.len; ++i) {
      if (nrgIdx[i].idx<0) I[k++]=i;
   }

   I.len=k;
}


template<class TQ, class TD>
QSpace<TQ,TD>& NRGData<TQ,TD>::getQSpace(
   const char *t, unsigned i
){
   int d=-1;

   if (!t || !index("KT",t[0]) || !index("KT",t[1]) || t[2])
   wblog(FL,"ERR invalid tag `%s'", t ? t : "(null)");

   if (t[0]=='K') {
      if (t[1]=='K')
           { if (i<KK.len) return KK[i]; else d=KK.len; }
      else { if (i<KT.len) return KT[i]; else d=KT.len; }
   }
   else {
      if (t[1]=='K')
           { if (i<TK.len) return TK[i]; else d=TK.len; }
      else { if (i<TT.len) return TT[i]; else d=TT.len; }
   }

   wblog(FL,"ERR index out of bounds (%s: %d/%d)",t,i,d);
   return A;
};


template<class TQ, class TD>
void NRGData<TQ,TD>::updateInfo(const char *F, int L,
  const char *vname, mxArray *a, char keep
) const {

   if (!a) wblog(F,L,
   "WRN update info `%s' with NULL mxArray", vname);

   if (MX) {
      mexPutVariable("caller", vname, a);
      if (!keep) mxDestroyArray(a);
   }
   else {
      tic_nrgIO.resume();
      matPutVariable(F,L,getFileNameI(),vname,a,keep);
      tic_nrgIO.stop();
   }
}


template<class TQ, class TD>
void NRGData<TQ,TD>::updateInfO(const char *F, int L,
  const char *vtag, mxArray *a, char keep
) const {

   if (!name.data || !name.data[0]) wblog(F,L,
   "WRN update info `%s' of unnamed operator", vtag);

   sprintf(str,"%s%s",
      name.data ? name.data : "",
      vtag ? vtag : ""
   );

   updateInfo(F,L,str,a,keep);
}


template<class TQ, class TD>
bool NRGData<TQ,TD>::opExists(char tflag) {

   if (tflag) { 
      unsigned i, k=(NRG_N>2 ? NRG_N-2 : 0);
      const char* tag[3] = { "KK","KT","TK" };

      for (i=0; i<3; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]); 
         if (!isQSpaceVec(str,k)) return 0;
      }
   }
   else {
      sprintf(str,"%sKK",name.data); 
      return isQSpaceVec(str,0);
   }

   return 1;
}


template<class TQ, class TD>
bool NRGData<TQ,TD>::initOp(const char *F, int L,
   const wbvector< QSpace<TQ,TD> > &C0,
   char tflag,
   NRGData<TQ,TD>* CR,
   const char *istr
){
   unsigned i,j,m,iter,e,ndup=0;
   const char* tag[3] = { "KK","KT","TK" };
   unsigned ntag=sizeof(tag)/sizeof(char*);
   wbvector< QSpace<TQ,TD> > NRG;
   mxArray *a; int r=-1;

   nrgIdx.initDef(C0.len); kloc=0;

   if (C0.len==0) {
      init(1); return 0;
   }

   if (istr && istr[0]) sprintf(str,"%s ",istr); else str[0]=0;

   if (!name.isName(99)) wblog(F,L,
      "ERR invalid name for %soperator (%d,`%s')",
       str, nrgIdx.len, name.data
   );

   for (i=0; i<C0.len; ++i) {
      r=-1;
      if (!C0[i].isConsistent(r) || r<2 || r>4) wblog(F,L,
         "ERR invalid operator (got rank-%d QSpace) ???\n"
         "%s (%s)",r,shortFL(FL),str);
      if (C0[i].QDIM!=C0[0].QDIM) wblog(FL,
         "ERR invalid operator (QDIM=%d, %d) ???\n%s:%d %s",
          C0[i].QDIM, C0[0].QDIM,FL,str);
   }

   KK.initDef(C0);

   for (i=0; i<C0.len; ++i) {
      if (nrgIdx[i].idx>=0) continue;
      for (j=i+1; j<C0.len; ++j)
      if (C0[j]==C0[i]) {
         nrgIdx[j].init(i,this);
         ++ndup;
      }
   }

   if (CR) {
      m=CR->KK.len;
      for (i=0; i<C0.len; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         for (j=0; j<m; ++j) if (C0[i]==CR->KK[j]) {
            nrgIdx[i].init(j,CR);
            ++ndup;
            break;
         }
      }
      sprintf(str,"%s%s",CR->name.data,"KK_"); 
      a=getMxData(str,0);
      if (a && mxIsQSpaceVec(FL,a)) {
         mxInitQSpaceVecR23(F,L,a,NRG); m=NRG.len;
         for (i=0; i<C0.len; ++i) {
            if (nrgIdx[i].idx>=0) continue;
            for (j=0; j<m; ++j) if (C0[i]==NRG[j]) {
               nrgIdx[i].init(j,CR->name.data);
               ++ndup;
               break;
            }
         }
      }

      if (ndup==nrgIdx.len) {
         return ndup;
      }
   }


   e=0;

   if (tflag) {
      iter=(NRG_N>2 ? NRG_N-2 : 0);
      for (i=0; i<ntag; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]); 
         if (!isQSpaceVec(str,iter)) { ++e; break; }
      }
   }

   iter=0;

   sprintf(str,"%s%s",name.data,"KK"); 
   if (!isQSpaceVec(str,iter)) ++e;

   sprintf(str,"%s%s",name.data,"KK_"); 
   a=getMxData(str,iter); if (!a || !mxIsQSpaceVecR23(FL,a)) ++e;

   if (e) {
      calc=store=1;
      return ndup;
   }

   mxInitQSpaceVecR23(F,L,a,NRG);

   for (i=0; i<C0.len; ++i) {
      if (nrgIdx[i].idx>=0 && nrgIdx[i].tag.isEmpty()) continue;
      for (j=0; j<NRG.len; ++j) if (C0[i]==NRG[j]) {
         nrgIdx[i].tag.init();
         nrgIdx[i].idx=j; ++ndup;
         break;
      }
   }

   calc  = (ndup==C0.len) ? 0 : 1;
   store = (ndup==C0.len) ? 0 : 1;

   return ndup;
};


template<class TQ, class TD>
bool NRGData<TQ,TD>::initOp(const char *F, int L,
   const wbvector< wbvector< QSpace<TQ,TD> > > &CC,
   NRGData<TQ,TD>* CR,
   const char *istr,
   char tflag
){
   unsigned i,j,l,m,iter,q=0,e, ndup=0;
   const char* tag[3] = { "KK","KT","TK" };
   unsigned ntag=sizeof(tag)/sizeof(char*);
   const mxArray *a; int r=2;

   wbvector< wbvector< QSpace<TQ,TD> > > NRG;

   nrgIdx.initDef(CC.len);

   if (CC.len==0) {
      init(1); return 0;
   }

   if (istr && istr[0])
   sprintf(str,"%s ",istr); else str[0]=0;

   if (!name.isName(99)) wblog(F,L,
      "ERR invalid name for %soperator (%d,`%s')",
       str, nrgIdx.len, name.data
   );


   kloc=-1;

   for (i=0; i<CC.len; ++i) {
   for (l=j=0; j<CC[i].len; ++j) { const QSpace<TQ,TD> &x = CC[i][j];
      if (x.isEmpty()) continue;
      l=j; if (!q) q=x.QDIM;

      if (!x.isConsistent(r) || r!=2) wblog(F,L,
         "ERR invalid operator (rank-%d QSpace) ???\n"
         "--> %s %s",r, shortFL(FL), str);
      if (x.QDIM!=q) wblog(FL,
         "ERR invalid operator (QDIM=%d, %d) ???\n--> %s %s",
          x.QDIM, q, shortFL(FL), str
      );
   } kloc=MAX(kloc,(int)l); }

   if (kloc<0) wblog(FL,"ERR all empty operators !??");
   kloc=MAX(0,kloc-1);

   KK.initDef(CC.len);
   CI.initDef(CC.len);
   for (i=0; i<CC.len; ++i) { 
      if (!CC[i].len) wblog(FL,"ERR empty operator set !??");
      CI[i].initDef(CC[i]);
      KK[i]=CC[i][0];
   }

   for (i=0; i<CC.len; ++i) { if (nrgIdx[i].idx>=0) continue;
      for (j=i+1; j<CC.len; ++j)
      if (CC[i].isEqual(CC[j])) { ++ndup;
         nrgIdx[j].init(i,this);
      }
   }

   if (CR && CR->CI.len) { m=CR->CI.len;
      for (i=0; i<CC.len; ++i) {
         if (nrgIdx[i].idx>=0) continue;

         for (j=0; j<m; ++j)
         if (CC[i].isEqual(CR->CI[j])) {
            nrgIdx[i].init(j,CR);
            ++ndup; break;
         }
      }

      sprintf(str,"%s%s",CR->name.data,"CI_"); 
      a=getMxInfo(str);
      if (a) { if (mxIsQSpaceVecVec(a,2)) {
         mxInitQSpaceVecVec(F,L,a,NRG); m=NRG.len;
         for (i=0; i<CC.len; ++i) {
            if (nrgIdx[i].idx>=0) continue;
            for (j=0; j<m; ++j) if (CC[i]==NRG[j]) {
               nrgIdx[i].init(j,CR->name.data);
               ++ndup; break;
            }
         }} else wblog(FL,
         "WRN %s::CI_ not a QSpace cell vector !??", CR->name.data);
      }

      if (ndup==nrgIdx.len) { return ndup; }
   }


   e=0;

   if (tflag && NRG_N>2) {
      iter=NRG_N-2;
      for (i=0; i<ntag; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]); 
         if (!isQSpaceVec(str,iter)) { ++e; break; }
      }
   }

   iter=0;
   sprintf(str,"%s%s",name.data,"KK"); 
   if (!isQSpaceVec(str,iter)) ++e;

   sprintf(str,"%s%s",name.data,"CI_"); 
   a=getMxData(str,iter); if (!a || !mxIsQSpaceVec(FL,a)) ++e;

   if (e) { calc=store=1;
      wblog(FL," *  calculate all operators %s (%d)",
      name.data, nrgIdx.len); return ndup;
   }

   mxInitQSpaceVecVec(F,L,a,NRG); m=NRG.len;

   for (i=0; i<CC.len; ++i) {
      if (nrgIdx[i].idx>=0 && nrgIdx[i].tag.isEmpty()) continue;
      for (j=0; j<m; ++j) if (CC[i].isEqual(NRG[j])) {
         nrgIdx[i].tag.init();
         nrgIdx[i].idx=j; ++ndup;
         break;
      }
   }

   calc  = (ndup==CC.len) ? 0 : 1;
   store = (ndup==CC.len) ? 0 : 1;

   return ndup;
}


template<class TQ, class TD>
void NRGData<TQ,TD>::initSym(const char *tag) {

  if (strlen(tag)==2) {
     const wbvector< QSpace<TQ,TD> > &K=getVec(tag);
     unsigned i,n=K.len;
     issym.init(n);

     for (i=0; i<n; ++i) {
        if (K[i].isHConj()) issym[i]=+1; else
        if (K[i].isAHerm()) issym[i]=-1;
     }
  }
  else if (strlen(tag)==1) {
     const QSpace<TQ,TD> &K=getQSpace(tag);
     issym.init(1);
     if (K.isHConj()) issym[0]=+1; else
     if (K.isAHerm()) issym[0]=-1;
  }
  else wblog(FL,"ERR invalid tag `%s'", tag ? tag : "(null)");
}


template<class TQ, class TD>
bool NRGData<TQ,TD>::isQSpaceVec(const char *nm, int iter) {

    if (iter<0) iter=NRG_ITER;
    if (iter>=(int)NRG_N) wblog(FL,
       "ERR Index out of bounds (%d/%d)", iter, NRG_N);

    const mxArray *a=getMxData(nm,iter,'i');
    if (a && mxIsQSpaceVecR23(FL,a)) return 1;

    return 0;

}


template<class TQ, class TD>
bool NRGData<TQ,TD>::checkVec(
   const char *F, int L,
   const char *XX, unsigned R, const char *istr
){
   unsigned l,n=128;
   char nm[n]; int fid;
   
   l=snprintf(nm,n,"%s%s", name.data, XX);
   if (l>=n) wblog(FL,"ERR strlen out of bounds (%d/%d)",l,n);

   if (MX) {
      const mwSize *s=mxGetDimensions(MX);
      const int n=mxGetNumberOfDimensions(MX);

      fid=mxGetFieldNumber(MX, nm);
      if (fid<0) {
         if (istr) { 
            wblog(F,L,"WRN invalid %s", istr);
            wblog(FL, "ERR field `%s' could not be found",nm);
         }
         return 1;
      }

      if (n!=2 || (s[0]!=1 && s[1]!=1)) {
         if (istr) {
            wblog(F,L,"WRN invalid %s (%s)",istr,nm);
            wblog(FL, "ERR NRG data must be row structure vector (%d,%d)",
            (int)s[0],(int)s[1]);
         }
         return 1;
      }

      if (!NRG_N) NRG_N=unsigned(s[0]*s[1]); else
      if (NRG_N!=unsigned(s[0]*s[1])) {
         if (istr) {
            wblog(F,L,"WRN invalid %s (%s)",istr,nm);
            wblog(FL, "ERR size mismatch (%d/%d)", NRG_N, int(s[0]*s[1]));
         }
         return 1;
      }

      return !mxsIsQSpaceVec(MX, fid, R);
   }
   else {

      int i=mxIsQSpaceVec(F_L,NAME.data, nm, R, 0,0,0, NRG_N);
      if (i<0) {
          wblog(F,L,"ERR invalid %s (%s)%N%N%s%N",istr?istr:"???", nm, str);
          return 1;
      }

      if (!NRG_N) NRG_N=i; else
      if (i!=(int)NRG_N) {
         wblog(F,L,"WRN invalid %s (%s)", istr?istr:"???", nm);
         wblog(FL,"ERR size mismatch (%d,%d)", i, NRG_N);
         return 1;
      }
   }

   return 0;
};


template<class TQ, class TD>
bool NRGData<TQ,TD>::checkVecVec(const char *XX, unsigned R) {

   unsigned l,n=128;
   char nm[n]; int fid;

   l=snprintf(nm,"%s%s",name.data,XX);
   if (l>=n) wblog(FL,"ERR strlen out of bounds (%d/%d)",l,n);
   
   if (MX) {
      fid=mxGetFieldNumber(MX, nm);
      if (fid<0) {
         sprintf(str, "%s:%d Field `%s' could not be found.",
         FL, nm); return 1;
      }

      if (!NRG_N) {
         const mwSize *s=mxGetDimensions(MX);
         const int n=mxGetNumberOfDimensions(MX);
         if (n!=2 || (s[0]!=1 && s[1]!=1)) { sprintf(str,
            "%s:%d NRG data must be row structure vector (%d,%d)",
             FL, (int)s[0], (int)s[1]); return 1;
         }
         NRG_N=unsigned(s[0]*s[1]);
      }
      return !mxsIsQSpaceVEC(MX, fid, R);
   }
   else {

      int i=mxIsQSpaceVEC(FL,NAME.data, nm, R);
      if (i<0) return 1;
      if (!NRG_N) NRG_N=i;
      if (i!=(int)NRG_N) {
         sprintf(str,"%s:%d Size mismatch (%d,%d)",
         FL, i, NRG_N); return 1;
      }
   }

   return 0;
}


template<class TQ, class TD>
mxArray* NRGData<TQ,TD>::getMxData(const char *nm, int iter, char ionly) {
    
    if (iter<0) iter=NRG_ITER;
    if (iter>=(int)NRG_N) wblog(FL,
    "ERR %s - index out of bound (%d/%d)", __FUNCTION__, iter, NRG_N);

    if (MX) {
       int fid=mxGetFieldNumber(MX, nm);
       if (fid<0) return NULL;
       return mxGetFieldByNumber(MX,iter,fid);
    }
    else { 
       Wb::matFile F(FL,getFileName(iter),"r");
       tic_nrgIO.resume();

       if (aux) mxDestroyArray(aux);
       if (ionly)
            aux = matGetVariableInfo(F.mfp, nm);
       else aux = matGetVariable    (F.mfp, nm);

       tic_nrgIO.stop();
       F.close(); return aux;
    }
};


template<class TQ, class TD>
const mxArray* NRGData<TQ,TD>::getMxInfo(const char *n) {

    if (!n || !n[0]) wblog(FL,
    "ERR calling %s with NULL", __FUNCTION__);

    if (aux) { mxDestroyArray(aux); aux=NULL; }

    if (MX) {
       return mexGetVariablePtr("caller",n);
    }
    else { 
       Wb::matFile F(FL,getFileNameI(),"r");

       tic_nrgIO.resume();
       aux=matGetVariable(F.mfp, n);
       tic_nrgIO.stop();

       return aux;
    }
}


template<class TQ, class TD>
mxArray* NRGData<TQ,TD>::toMx() const {

    mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

    mxAddField2Scalar(FL, S, "name",    name.toMx());
    mxAddField2Scalar(FL, S, "calc",    numtoMx(calc));
    mxAddField2Scalar(FL, S, "store",   numtoMx(store));
    mxAddField2Scalar(FL, S, "nrgIdx",  nrgIdx.toMx());

    mxAddField2Scalar(FL, S, "CI",      CI.toMx());
                                       
    mxAddField2Scalar(FL, S, "KK",      KK.toMx());
    mxAddField2Scalar(FL, S, "KT",      KT.toMx());
    mxAddField2Scalar(FL, S, "TK",      TK.toMx());
    mxAddField2Scalar(FL, S, "TT",      TT.toMx());
    mxAddField2Scalar(FL, S, "A",       A.toMx());
    mxAddField2Scalar(FL, S, "K",       K.toMx());
    mxAddField2Scalar(FL, S, "T",       T.toMx());

    mxAddField2Scalar(FL, S, "nrgNAME", NAME.toMx());
    mxAddField2Scalar(FL, S, "nrgMX",   num2Str(MX,"0x%lX").toMx());
    mxAddField2Scalar(FL, S, "nrgN",    numtoMx(NRG_N));
    mxAddField2Scalar(FL, S, "nrgITER", numtoMx(NRG_ITER));

    return S;
}


template<class TQ, class TD>
void NRGData<TQ,TD>::info(const char *vstr) const {
    
    unsigned i,nl=0;
    wbstring s(nrgIdx.len+1);

    for (i=0; i<nrgIdx.len; ++i) {
       if ( nrgIdx[i].idx<0) s[i]='C';
       if (!nrgIdx[i].nrg) s[i]='R';
       else if ( nrgIdx[i].nrg==this)
            s[i]='1';
       else s[i]='2';
    }

    sprintf(str,"%d/%d; calc=%d (%s; %d); store=%d",
    NRG_ITER+1, NRG_N, calc, s.data, nrgIdx.len, store);

    if (MX) {
       printf("\n%s Structure `%s' (0x%lX), %s\n",
       vstr && vstr[0] ? vstr : (name.data ? name.data : "ans"),
       NAME.data ? NAME.data : "", MX, str);
    }
    else {
       printf("\n%s File structure `%s', %s\n",
       vstr && vstr[0] ? vstr : (name.data ? name.data : "ans"),
       NAME.data ? NAME.data : "", str);
    }

    if (!KK.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<KK.len; ++i) {
          sprintf(str,"KK[%d]", i);
          KK[i].info(str);
       }
    }

    if (!KT.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<KT.len; ++i) {
          sprintf(str,"KT[%d]", i);
          KT[i].info(str);
       }
    }
    if (!TK.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<TK.len; ++i) {
          sprintf(str,"TK[%d]", i);
          TK[i].info(str);
       }
    }
    if (!TT.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<TT.len; ++i) {
          sprintf(str,"TT[%d]", i);
          TT[i].info(str);
       }
    }

    nl=0;

    if (!A.isEmpty()) { if (!nl) { printf("\n"); ++nl; }; A.info("A"); }
    if (!K.isEmpty()) { if (!nl) { printf("\n"); ++nl; }; K.info("K"); }
    if (!T.isEmpty()) { if (!nl) { printf("\n"); ++nl; }; T.info("T"); }

    if (nl) printf("\n");
};


template <class TQ, class TD>
void NRGData<TQ,TD>::dispIdx(const char *F, int L) const {

   wblog(F,L,"%NNRGIndex %s (calc=%d, store=%d)%N", name.data, calc, store);

   for (unsigned i=0; i<nrgIdx.len; ++i) {
      if (!nrgIdx[i].tag.isEmpty())
           strcpy(str, nrgIdx[i].tag.data);
      else if (nrgIdx[i].nrg)
           strcpy(str, nrgIdx[i].nrg==this ? "(this)" : "(that)");
      else if (nrgIdx[i].idx>=0)
           strcpy(str, "(file)");
      else strcpy(str, "(calc)");

      printf("%6d: %4d  %s:%lX\n", i, nrgIdx[i].idx, str, nrgIdx[i].nrg);
   }

   if (nrgIdx.len) printf("\n");
}

template <class TQ, class TD>
void NRGData<TQ,TD>::info(const char *F, int L, const char *istr) const {

   unsigned i, nt1=0, nt2=0, nx1=0, nx2=0, nc=0;
   wbstring mark(nrgIdx.len+1);
   wbstring S(128);

   for (i=0; i<nrgIdx.len; ++i) {

      if (nrgIdx[i].nrg && nrgIdx[i].nrg!=this) { ++nx1; mark[i]='x'; } else
      if (!nrgIdx[i].tag.isEmpty()) { ++nx2; mark[i]='X'; } else
      if (nrgIdx[i].nrg==this) { ++nt1; mark[i]='r'; } else
      if (nrgIdx[i].idx>=0) { ++nt2; mark[i]='R'; }
      else { ++nc; mark[i]='*'; }
   }


   if (nc) {
      S.pushf(FL,", %s%d/%d op%s",
      calc ? "calc ":"", nc, KK.len, KK.len!=1 ? "s":"");
   }
   else {
      S.pushf(FL,", %d op%s%s", KK.len, KK.len!=1 ? "s":"",
      calc ? ", calc!??":", all stored or referenced");
   }

   if (store)
        S.push(FL," and store");
   else S.push(FL," (not stored)");

   if (kloc      ) S.pushf(FL,", nloc=%d", kloc+1);

   if (strcmp(istr,name.data)) {
      wblog(F,L," *  %s (%3s): [%s] %s",
      name.data, istr, mark.data, S.data[0] ? S.data+2 : "");
   }
   else {
      wblog(F,L," *  %s [%s] %s",
      name.data, mark.data, S.data[0] ? S.data+2 : "");
   }
};


template <class TQ, class TD>
void NRGData<TQ,TD>::disp_otype(const char *F, int L) const {

   wbstring s(128); unsigned i=0;

   if (name.data) s.cpy(FL,name.data);

   if (KK.len) {
      s.pushf(FL,"%sKK: ",i?"; ":"");
      for (i=0; i<KK.len; ++i) s.pushf(FL,"%d",KK[i].otype);
   }
   if (TK.len) {
      s.pushf(FL,"%sTK: ",i?"; ":"");
      for (i=0; i<TK.len; ++i) s.pushf(FL,"%d",TK[i].otype);
   }
   if (KT.len) {
      s.pushf(FL,"%sKT: ",i?"; ":"");
      for (i=0; i<KT.len; ++i) s.pushf(FL,"%d",KT[i].otype);
   }
   if (TT.len) {
      s.pushf(FL,"%sTT: ",i?"; ":"");
      for (i=0; i<TT.len; ++i) s.pushf(FL,"%d",TT[i].otype);
   }

   wblog(F_L,"TST %s() %s ",FCT,s);
};



template <class TQ, class TD> inline
void updateOp_L(const char *F, int L,
   const QSpace<TQ,TD> &Xin,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout
){
   QSpace<TQ,TD> Xk; unsigned r=-1;


   if (Xin.isEmpty()) { wblog(F_L,
      "ERR %s() got empty L-operator",FCT); return; }
   if (Xin.isOperator(&r)<=0) { Xin.info("L-op");
      wblog(F_L,"ERR %s() got invalid L-operator (%d)",FCT,r);
   }


   Xin.contract(F,L,2,A2,1,Xk);
   A1.contract(F,L,"1,3;*",Xk, (r!=3 ? "1,3":"1,4"), Xout);

   if (r==3) Xout.Permute("1,3,2");

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};


template <class TQ, class TD> inline
void updateOp_s(const char *F, int L,
   const QSpace<TQ,TD> &x,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout
){
   QSpace<TQ,TD> Xk; unsigned r=-1;


   if (x.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if (x.isOperator(&r)<=0) { x.info("local-op");
      wblog(FL,"ERR %s() got invalid L-operator (%d)",FCT,r);
   }

   A2.contract(F,L,3,x,2, Xk);
   A1.contract(F,L,"1,3;*", Xk,"1,3", Xout);

   Xout.otype=x.otype;
   Xout.SkipZeroData();
};


template <class TQ, class TD>
inline void updateOp_Ls(
   const QSpace<TQ,TD> &Xin,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout
){
   QSpace<TQ,TD> Xk; unsigned r=-1;

   if (Xin.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if (Xin.isOperator(&r,'x')<=0 || r!=4) { Xin.info("Ls-op");
      wblog(FL,"ERR %s() got invalid Ls-operator (%d)",FCT,r);
   }

   Xin.contract("3,4",A2,"1,3", Xk);
   A1.contract("1,3;*",Xk,"1,2",Xout);

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD>
inline void updateOp_sL(
   const QSpace<TQ,TD> &Xin,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout
){
   QSpace<TQ,TD> Xk; unsigned r=-1;

   if (Xin.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if (Xin.isOperator(&r)<=0 || r!=4) { Xin.info("sL-op");
      wblog(FL,"ERR %s() got invalid sL-operator (%d)",FCT,r);
   }

   Xin.contract("4,3",A2,"1,3", Xk);
   A1.contract("1,3;*",Xk,"2,1",Xout);

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD> inline
void updateOp_loc(const char *F, int L,
   const NRGData<TQ,TD> &NRG, unsigned i,
   const QSpace<TQ,TD> &XL,
   const QSpace<TQ,TD> &x,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout,
   unsigned iter
){
   QSpace<TQ,TD> xA,XA;
   unsigned r1=-1, r2=-1;
   bool o1=0, o2=0;

   if (!x.isEmpty()) {
      if (x.isOperator(&r2)<=0) { x.info("local_op"); wblog(FL,
         "ERR %s() got invalid local operator (%d)",FCT,r2); }
      o2=(int(r2)==3);

      wblog(FL,"TST contracting %s[%d] onto A[%d].s",NRG.name.data,i,iter);
      A2.contract(F,L,3,x,2,xA);
   }
   else xA.init2ref(A2);

   if (iter || !XL.isEmpty()) {
      if (XL.isOperator(&r1)<=0) { x.info("L-op"); wblog(FL,
         "ERR %s() got invalid L-operator (%d)",FCT,r1); }
      o1=(int(r1)==3);

      if (iter<1) wblog(FL,"TST contracting onto A[%d].L",iter);
      if (o1 && o2) { o1=o2=0;
           XL.contract(F,L,"2,3",xA,"2,4",XA); }
      else XL.contract(F,L,2,xA,2,XA);
   }
   else XA.init2ref(xA);

   if (o1)
        A1.contract(F,L,"1,3;*",XA,"1,4", Xout,"1,3,2");
   else A1.contract(F,L,"1,3;*",XA,"1,3", Xout);

   if (XL.otype!=x.otype) wblog(FL,
      "WRN %s() got different operator types (%s; %s)",
       FCT, QS_STR[XL.otype],QS_STR[x.otype]
   ); 

   Xout.otype=XL.otype;
   Xout.SkipZeroData();
};


template <class TQ, class TD>
int Wb::updateOp(const char *F, int L,
   const NRGData<TQ,TD> &NRG,
   wbvector< QSpace<TQ,TD> > &F12,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   const wbvector< QSpace<TQ,TD> > &FKK,
   unsigned iter
){
   unsigned i, r=-1, rx, j=0, m=0, nop=FKK.len;
   int e=0;

   const wbvector< wbvector<QSpace<TQ,TD> > > &CI=NRG.CI;
   const wbvector< NRGIndex<TQ,TD> > &nrgIdx=NRG.nrgIdx;

   if (A1.isEmpty() || A2.isEmpty()) { F12.initDef(nop); return 0; }

   if (nrgIdx.len!=nop) wblog(FL,
      "ERR %s() severe size mismatch (nrgIdx: %d,%d)",FCT,nrgIdx.len,nop);

   if (iter>0) {
      for (i=0; i<nop; ++i) { rx=r;
         if (nrgIdx[i].idx>=0) continue;
         if (!FKK[i].isConsistent(FL)
         || int(r=FKK[i].isOperator(&rx))<=0) break;
      }
   }
   else if (CI.isEmpty()) {
      for (i=0; i<nop; ++i) { rx=r;
      if (!FKK[i].isConsistent(FL)
      || int(r=FKK[i].isOperator(&rx,'x'))<=0) break; }
   }
   else {
      if (nop!=CI.len) wblog(FL,
         "ERR size inconsistency (%d/%d)",nop,CI.len);

      for (i=0; i<nop; ++i) { m=CI[i].len;
      for (j=0; j<m; ++j) { rx=r;
         if (!CI[i][j].isEmpty() && (!CI[i][j].isConsistent(FL)
         || int(r=CI[i][j].isOperator(&rx))<=0)) break;
      }}
   }

   if (i<nop || j<m) wblog(F,L,
      "ERR invalid operator rank-%d (iter=%d, op %d/%d)",r,iter,i+1,nop);
   if (A1.rank(F,L)!=3 || A2.rank(F,L)!=3 || A1.QDIM!=A2.QDIM) wblog(F,L,
      "ERR %s() severe rank inconsistency (%s; %s)",
       FCT, A1.sizeStrQ().data, A2.sizeStrQ().data
   );

   if (typeid(TD)!=typeid(double)) wblog(F,L,
   "ERR A is expected to be real! (%s)",getName(typeid(TD)).data);


   if (F12.len!=nop)
   F12.initDef(nop);


   if ((int)iter>NRG.kloc) {
      for (i=0; i<nop; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         updateOp_L(F,L, FKK[i], A1, A2, F12[i] );
         wbtop().runningLarge(FL);

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }
   }
   else if (!CI.isEmpty()) {
      QSpace<TQ,TD> E;

      for (i=0; i<nop; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         updateOp_loc(F,L,NRG,i,
            FKK[i], iter+2>CI[i].len ? E : CI[i][iter+1],
            A1, A2, F12[i],iter
         );

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }
   }
   else {
      char Lflag = A1.getDim(1)>1 || A2.getDim(1)>1;
      unsigned l=0, tlen=16; char tag[tlen];

      wbstring mark(nop,'*');

      static time_t last_call_iter0=0;
      time_t tnow=time(0);

      for (i=0; i<nop; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         r=FKK[i].rank(FL);

         if (r==2 || (r==3 && FKK[i].otype==QS_OPERATOR)) {
            if (Lflag &&
                FKK[i].hasQOverlap(1,A1,1,'<') &&
                FKK[i].hasQOverlap(2,A1,1,'<')
            ){
               updateOp_L(F,L,FKK[i],A1,A2, F12[i]);
               mark[i]='L';
            }
            else
            if (FKK[i].hasQOverlap(1,A1,3,'<') &&
                FKK[i].hasQOverlap(2,A1,3,'<')
            ){
               updateOp_s(F,L,FKK[i], A1, A2, F12[i]);
               mark[i]='s';
            }
            else
            if (FKK[i].hasQOverlap(1,A1,2,'<') &&
                FKK[i].hasQOverlap(2,A1,2,'<')
            ){ 
               mark[i]='R';
            }
            else wblog(FL,
              "ERR %s() failed to deal input operator with A%d",FCT,iter);
         }
         else if (r==4) {

            wbvector<int> m1,m2;
            int mm[]= {
                FKK[i].hasQOverlap(1,A1,1,'<'),
                FKK[i].hasQOverlap(2,A2,3,'<'),
                FKK[i].hasQOverlap(3,A2,1,'<'),
                FKK[i].hasQOverlap(4,A2,3,'<'),

                FKK[i].hasQOverlap(1,A1,3,'<'),
                FKK[i].hasQOverlap(2,A2,1,'<'),
                FKK[i].hasQOverlap(3,A2,3,'<'),
                FKK[i].hasQOverlap(4,A2,1,'<')
            };
            m1.init(4,mm); m2.init(4,mm+4);

            if (m1.allUnequal(0)) {
               updateOp_Ls(FKK[i], A1, A2, F12[i]);
               mark[i]='T';
            }
            else 
            if (m2.allUnequal(0)) {
               updateOp_sL(FKK[i], A1, A2, F12[i]);
               mark[i]='t';
            }
            else {
               if (STRICT_ITER0) { mark[i]='?';
               if (mark.count('?')==1 && last_call_iter0+5<tnow){
               wblog(FL,
                 "WRN %NFailed to fully associate rank-4 tensor product\n"
                 "operator with first site. [ %s; %s ]\n"
                 "A1: [ %s ];  A2: [ %s ]%N%N"
                 "NB! Hint: A0 may not be complete due to truncation "
                 "in NRGWilsonQS already at iter=0 !?%N",
                  m1.toStr().data, m2.toStr().data,
                  A1.getDim().toStrf("%d","x").data,
                  A2.getDim().toStrf("%d","x").data);
               }}

               updateOp_Ls(FKK[i], A1, A2, F12[i]);
            }
         }
         else wblog(FL,
        "ERR %s() invalid rank-%d operator (iter=%d)",FCT,r,iter);

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }


      l=snprintf(tag,tlen,"FDM %2d/%d",iter,NRG_N); {
         if (NRG.name.data && NRG.name.data[0] && l<tlen)
            l+=snprintf(tag+l,tlen-l," %s:",NRG.name.data);
         if (l>=tlen) wblog(FL,
        "ERR %s() string out of bounds (%d/%d)",FCT,l,tlen);
      }

      if (!nop)
         wblog(FL,"%s got empty spectral ops",tag);
      else {
         unsigned i1=0, m=mark.count('*');
         for (i=0; i<nop; ++i) { if (mark[i]!='*') { i1=i; break; }}
         for (; i<nop; ++i) { if (mark[i]!='*') {
            if (mark[i]!=mark[i1]) break;
         }}

         if (m==nop) { wblog(FL,
            "%s no spectral ops to contract (all referenced)",
            tag,mark.data);
         }
         else if (i==nop && !m) {
            char ostr[32];
            if (nop>1) 
                 sprintf(ostr,"all %d spectral ops",nop);
            else sprintf(ostr,"spectral operator");

            if (mark[0]!='R') wblog(FL,
               "%s contracting %s onto %c",tag,ostr,mark[0]);
            else wblog(FL,
               "%s %s already assumed in %c",tag,ostr,mark[0]
            );
         }
         else if (i==nop && m) {
            if (mark[i1]!='R') wblog(FL,"%s contracting %d/%d "
               "spectral ops onto [%s]",tag,nop-m,nop,mark.data);
            else wblog(FL,"%s all %d/%d spectral "
               "ops already assumed in %s",tag,nop-m,nop,mark.data);
         }
         else {
            wblog(FL,"%s contracting "
            "%d/%d spectral ops onto [%s]",tag,nop-m,nop,mark.data);
         }
      }

      last_call_iter0=tnow;
   }

   if (iter==0) {
      for (i=0; i<F12.len; ++i) {
         if (F12[i].itags.len>2 && !F12[i].itags[2]) {
            F12[i].itags[2].init("op");
         }
      }
   }

   wbtop().runningLarge(FL);
   return e;
};


template <class TQ, class TD>
unsigned applyZ0(const QSpace<TQ,TD> &Z0, QSpace<TQ,TD> &CK){

   unsigned r=CK.rank(FL);
   const QS_TYPES otype=CK.otype;

   if (r==2 || r==3) {
      if (r==3 && CK.otype!=QS_OPERATOR) wblog(FL,
         "ERR %s() expecting otype==%s for rank-%d op (%s)",
          FCT,QS_STR[QS_OPERATOR],r,QS_STR[CK.otype]);
      if (Z0.rank(FL)!=2) wblog(FL,
         "ERR %s() Z0 of rank-2 required for rank-%d op (%d)",
          FCT,r,Z0.rank(FL)
      );

      if (!Z0.hasQOverlap(2, CK,1,'>')) {
         wblog(FL,"ERR %s() Z0 does not match operator space",FCT);
         return 1;
      }
      Z0.contract(2,CK,1,CK);
      CK.otype=otype;
   }
   else if (r==4) {
      if (Z0.rank(FL)!=4) wblog(FL,
         "ERR tensor Z0 required for tensor operator (rank=%d/%d)",
          Z0.rank(FL), r
      );

      if (!Z0.hasQOverlap(3, CK,1,'>') ||
          !Z0.hasQOverlap(4, CK,2,'>')) {
         wblog(FL,"ERR %s() Z0 does not match operator space",FCT);
         return 2;
      }
      Z0.contract("3 4",CK,"1 2",CK);
      CK.otype=otype;
   }
   else wblog(FL,"ERR %s() invalid rank=%d.",FCT,r);

   return 0;
};


template <class TQ, class TD>
unsigned getDLoc(
   NRGData<TQ,TD> &A, wbSigHandler *sig=NULL, unsigned *Dk=NULL,
   char fflag=0
){
   unsigned iter, N=NRG_N; INDEX_T d0=0, d=0;
   char cgflag=-1;

   if (Dk) { (*Dk)=0; }; printf("\r");

   for (iter=0; iter<N; ++iter) { if (sig) sig->call99();
      printf("  %s: %2d/%d loading AK ...   \r",FCT,iter,N);
      fflush(0);

      A.init(FL,"K",iter);

      if (iter==0) cgflag=A.K.gotCGS(FL);

      if (A.K.isEmpty()) {
         if (iter==N-1) continue; else wblog(FL,"ERR %s() "
        "empty AK at intermediate iteration %d/%d",FCT,iter,N);
      }
      else {
         A.K.isConsistentR(3,FL);
         if (cgflag<=1 && !A.K.QIDX.isUnique()) wblog(FL,
            "ERR invalid AK at iteration %d/%d",iter,N
         );
      }

    #ifdef NSAFEGUARDS
      if (fflag) { QSpace<TQ,TD> E;
         A.K.contract("1,3;*",A.K,"1,3",E);
         if (!E.isIdentityMatrix(1E-12)) {
            MXPut(FL).add(A,"A").add(E,"E");
            wblog(FL,"ERR NRG[%d].AK not in LRs index order",iter,N);
         }
      }
    #endif

      A.K.getDim(2,&d);

      if (iter>2 || d0) {
         if (d!=d0) {
            MXPut(FL,"a").add(A.K,"K").add(d,"d").add(d0,"d0")
              .add(iter,"iter").add(N,"N");
            wblog(FL,"ERR %s() inconsistent local dimension "
              "(d=%d/%d @ k=%d/%d)",FCT,d,d0,iter,N
            );
         }
      }
      else if (iter==2 || iter+2==N) { d0=d; }
      else {
          A.init(FL,"T",iter);
          if (!A.T.isEmpty()) { d0=d; }
      }

      if (Dk) {
         d=A.K.getDim(1);
         if (*Dk<d) { *Dk=d; continue; }
      }
      if (!fflag && iter>2) break;
   }

   if (!d0) wblog(FL,"ERR %s() "
      "failed to get local dimension (d=%d; %d/%d)",FCT,d0,iter,N);
   return d0;
};



#endif

