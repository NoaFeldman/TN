#ifndef __WB_ORTHO_HCC__
#define __WB_ORTHO_HCC__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// library routine to orthonormalize PSI
// by grouping its indizes into two groups

template <class TQ, class TD>
mxArray* mpsOrthoQS(
   const QSpace<TQ,TD> &PSI, // main input, i.e. Psi1 contracted with Psi2
   QSpace<TQ,TD> &A1,
   QSpace<TQ,TD> &A2,
   unsigned K, unsigned Nkeep,
   double stol,        // relative tolerance of singular values to keep
   const char *info
// =d: display;  =[sS]: (extended) store into mx structure
// examples: "ds", "dS"
);

char mxGetLRDir(
   const char *F, int L, const mxArray *argin, unsigned k
);

template <class TQ, class T1, class T2, class T3>
void twoSiteInit(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &PSI, unsigned &K,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
);

template <class TQ, class T1, class T2, class T3>
void twoSiteFinal(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &A1,
   QSpace<TQ,T3> &A2,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
);


template <class TQ, class TD>
class QBlock {
  public:
    QBlock() {};

    QSpace<TQ,TD> AS;
    wbindex I0;

    wbMatrix<TQ> Q1,Q2;
    wbMatrix<unsigned> S1,S2;
    wbarray<TD> U,V;
    wbvector<double> S;

    wbindex Ik;

  protected:
  private:
};


template <class TQ, class TD>
class SVDData {
  public:
    SVDData() {};

    void blockSVD(const QSpace<TQ,TD> &PSI, const unsigned K);

    unsigned dmrgTruncate(
       unsigned Nkeep0, double stol,
       wbvector<double> &SV, QSpace<TQ,TD> &A1, QSpace<TQ,TD> &A2,
       mxArray *S=NULL, char disp=0
    );

    void getQSpaceUSV(
       QSpace<TQ,TD> &U, QSpace<TQ,double> &S, QSpace<TQ,TD> &V,
       mxArray *I
    ) const;

    mxArray* toMx() const;

    void put(const char *vname, const char *ws="base") const {
       put(0,0,vname,ws);
    };
    void put(const char *F, int L,
       const char *vname, const char *ws="base"
    ) const {
       mxArray *a=toMx();
       int i=mexPutVariable(ws,vname,a);
       if (i) wblog(F,L,"ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F,L,"I/O putting variable `%s' to `%s'",vname,ws);
       mxDestroyArray(a);
    };

    wbvector< QBlock<TQ,TD> > QB;

  protected:
  private:
};


template <class TQ, class TD>
void SVDData<TQ,TD>::blockSVD(
   const QSpace<TQ,TD> &PSI,
   const unsigned K
){
   unsigned s,l,n, d, r=PSI.rank();
   wbvector<unsigned> D,D1,D2;
   wbMatrix<TQ> Q1;
   wbarray<TD> MM;
   wbindex I1(K,'i');
   wbperm P;

   if (PSI.isEmpty()) wblog(FL,"ERR PSI is empty ???");
   if (K>=r) wblog(FL,"ERR Group size exceeds rank (%d,%d)", K, r);

   PSI.getQsum(I1,Q1).groupRecs(P,D); n=D.len;
   QB.init(n);

   for (l=s=0; s<n; s++, l+=d) {
      QBlock<TQ,TD> &b=QB[s]; d=D[s];

      b.I0.init(d,P.data+l);
      PSI.getSub(b.I0, b.AS);

      b.AS.toBlockMatrix(MM, K, b.Q1, b.Q2, b.S1, b.S2, D1, D2);


      wbSVD(MM, b.U, b.S, b.V);
   }
}



template <class TQ, class TD>
unsigned SVDData<TQ,TD>::dmrgTruncate(
   unsigned Nkeep0, double stol,
   wbvector <double> &SV,
   QSpace<TQ,TD> &A1,
   QSpace<TQ,TD> &A2,
   mxArray *S,
   char disp
){

   unsigned i,j,i0,d,n, nk, nt, nq=QB.len, QDIM, flag=0, Nkeep=Nkeep0;
   wbvector<char> mark;
   wbperm P;
   int ix;

   wbMatrix<TQ> QS;
   wbMatrix<unsigned> DS; 
   wbarray<TD> VS;
   QSpace<TQ,TD> A;

   double snorm,sfac=1,sref,eps, s2t=0, stol2=stol*stol;

   QDIM=QB[0].AS.QDIM; if (!QDIM) wblog(FL,"ERR QDIM=%d ???", QDIM);
   str[0]=0;

   for (n=i=0; i<nq; i++) n+=QB[i].S.len;
   SV.init(n);

   for (i0=i=0; i<nq; i++, i0+=d) {
      const wbvector<double> &s=QB[i].S; d=s.len;
      memcpy(SV.data+i0, s.data, d*sizeof(double));
   }

   SV.Sort(P,-1);

   if (Nkeep>SV.len) Nkeep=SV.len;

   snorm=SV.norm();
   stol*=snorm;

   eps=1E-6*stol;

   if (SV.last()>stol) { if (SV.len>Nkeep) {
       if (stol!=0) {
          sprintf(str,"WRN all SVD>stol (Smin=%.3g (%g), %d/%d)",
          SV.last(), stol, Nkeep, SV.len); flag++;
       }
       else {
          sprintf(str,"keep all (stol=%g; Smin=%.3g, %d/%d)",
          stol, SV.last(), Nkeep, SV.len);
       }
   }}
   else {
      ix=SV.findClosestSorted(stol);
      if (ix<0) wblog(FL,"ERR find closest to stol (%d) failed !??",ix);

      i=unsigned(ix+1);
      if (i+1<Nkeep) { sprintf(str,
        "ok  Nkeep: %2d->%2d (Smin=%.3g < %.3g < stol=%.3g)",
         Nkeep, i+1, SV.last(), SV[i], stol);

         Nkeep=i+1;
      }
      else if (Nkeep<SV.len) { sprintf(str,
        "WRN Nkeep too small (%d/%d: Smin=%.3g < %.3g <= stol=%.3g)",
         Nkeep, SV.len, SV.last(), SV[MAX(0,int(Nkeep)-1)], stol);
         flag++;
      }
   }

   sref = MAX(eps, SV[MIN(Nkeep,SV.len)-1] - eps);

   for (nk=nt=i=0; i<nq; i++) {
      QBlock<TQ,TD> &b=QB[i]; d=b.S.len; mark.init(d);
      for (j=0; j<d; j++) {
         if (b.S[j]>=sref) mark[j]=1;
         else s2t+=b.S[j]*b.S[j];
      }

      mark.find(b.Ik); nk+=b.Ik.len; nt+=(mark.len-b.Ik.len);
   }

   if (nk+nt!=SV.len) wblog(FL, 
   "WRN Severe data inconsistency (%d+%d =? %d)", nk, nt, SV.len);
   if (!nk) wblog(FL, "WRN all data truncated !??");

   sref=sqrt(snorm*snorm-s2t);
   if (sref>0) {
      sfac=snorm/sref; if (fabs(sfac-1)>1E-2) {
         unsigned ns=strlen(str); sprintf(str+ns,
         "%sWRN readjusting norm by factor %.4g",ns ? "\n":"", sfac);
         flag++;
      }

      if (sfac!=1 && s2t!=0)
      for (i=0; i<nq; i++) QB[i].S*=sfac;
   }
   else wblog(FL,
   "WRN All states discarded !?? (%.4g, %.4g)",snorm,sqrt(s2t));

   if (S) {
      mxAddField2Scalar(FL, S, "stol",   numtoMx(stol));
      mxAddField2Scalar(FL, S, "svd",    SV.toMx());
      mxAddField2Scalar(FL, S, "svd2tr", numtoMx(s2t/(snorm*snorm)));
      mxAddField2Scalar(FL, S, "sfac",   numtoMx(sfac));
      mxAddField2Scalar(FL, S, "Nkeep",  numtoMx(nk));
      mxAddField2Scalar(FL, S, "Ntot",   numtoMx(SV.len));
      mxAddField2Scalar(FL, S, "flag",   numtoMx(flag));
      mxAddField2Scalar(FL, S, "info",   mxCreateString(str));
   }
   else if (disp) {
      double dbl=SV.norm2(); wblog(FL,
      "<i> sum(SVD) = %g = %g%+.2g", dbl, round(dbl), dbl-round(dbl));

      wblog(FL,"TST keep %d/%d states (%d)",nk, SV.len, Nkeep0);
      if (nt) wblog(FL,":x: Truncate %d/%d states.", nt, SV.len);
   }

   A1.clearQSpace();
   A2.clearQSpace();

   for (i=0; i<nq; i++) {
      QBlock<TQ,TD> &b=QB[i]; if (b.Ik.isEmpty()) continue;

      b.U.Select0(b.Ik,1);
      b.S.Select(b.Ik);
      b.V.Select0(b.Ik,1);

      wbDMatProd(b.V,b.S,VS);

      b.Q1.blockSum(QDIM,QS);
         if (!QS.recAllEqual()) wblog(FL,
         "ERR %s() all recs in QS must be equal",__FUNCTION__);
      QS.Resize(1,QS.dim2);

      if (b.U.rank()!=2) wblog(FL,"ERR invalid rank-%d",b.U.rank());
      DS.init(1,1); DS[0]=b.U.SIZE[1];

      A.initFromBlockMatrix(
         b.U, b.Q1, QS, QDIM,
         b.S1, DS
      );

      if (A.Append2AndDestroy(FL,A1,'u')) {
         sprintf(str+strlen(str), "\nQtot spaces must not overlap!\n\n");
         wberror(FL,str);
      }

      if (VS.rank()!=2) wblog(FL,"ERR invalid rank-%d",VS.rank());
      DS[0]=VS.SIZE[1];

      A.initFromBlockMatrix(
         VS, b.Q2, QS, QDIM,
         b.S2, DS
      );

      A.SkipZeroData(stol2,'b');

      if (A.Append2AndDestroy(FL,A2,'u')) {
         sprintf(str+strlen(str),"\nQtot spaces must not overlap!\n\n");
         wberror(FL,str);
      }
   }

   A2.PermuteLastTo(0);

   return flag;
}


template <class TQ, class TD>
void SVDData<TQ,TD>::getQSpaceUSV(
   QSpace<TQ,TD> &U,
   QSpace<TQ,double> &S,
   QSpace<TQ,TD> &VT,
   mxArray *I
) const {

   unsigned i,j,i0,d,n, nq=QB.len, QDIM;

   wbMatrix<TQ> QS;
   wbvector <double> SV;
   wbMatrix<unsigned> DS(1,1);
   QSpace<TQ,TD> A;
   QSpace<TQ,double> X;
   wbperm P;

   double snorm,eps,s2t=0;

   QDIM=QB[0].AS.QDIM; if (!QDIM) wblog(FL,"ERR QDIM=%d ???", QDIM);
   str[0]=0;

   U.clearQSpace();
   S.clearQSpace();
   VT.clearQSpace();

   for (n=i=0; i<nq; i++) n+=QB[i].S.len;
   SV.init(n);

   for (i0=i=0; i<nq; i++, i0+=d) {
      const wbvector<double> &sb=QB[i].S; d=sb.len;
      memcpy(SV.data+i0, sb.data, d*sizeof(double));
   }

   SV.Sort(P,-1);

   snorm=SV.norm(); eps=snorm*1E-10;

   for (i=0; i<nq; i++) {
      QBlock<TQ,TD> &b=QB[i];
      if (b.S.norm()<eps) { s2t+=b.S.norm2(); continue; }

      b.Q1.blockSum(QDIM,QS);
      if (!QS.recAllEqual()) wblog(FL, 
      "ERR All recs in QS should be equal !??");

      QS.Resize(1,QS.dim2);

      X.init(1,2,QDIM);
      X.DATA[0]->initp(1,b.S.len,b.S.data);
      for (j=0; j<2; j++) X.SetQ(0,j,QS.data);

      if (X.Append2AndDestroy(FL,S,'u')) {
         sprintf(str+strlen(str), "\nQtot spaces must not overlap!\n\n");
         wberror(FL,str);
      }

      DS[0]=b.U.dim2;
      A.initFromBlockMatrix(b.U, b.Q1, QS, QDIM, b.S1, DS);

      if (A.Append2AndDestroy(FL,U,'u')) {
         sprintf(str+strlen(str), "\nQtot spaces must not overlap!\n\n");
         wberror(FL,str);
      }

      DS[0]=b.VT.dim1;
      A.initFromBlockMatrix(b.VT, QS, b.Q2, QDIM, DS, b.S2);

      if (A.Append2AndDestroy(FL,VT,'u')) {
         sprintf(str+strlen(str),"\nQtot spaces must not overlap!\n\n");
         wberror(FL,str);
      }
   }

   if (I) {
      mxAddField2Scalar(FL, I, "svd",    SV.toMx());
      mxAddField2Scalar(FL, I, "snorm",  numtoMx(snorm));
      mxAddField2Scalar(FL, I, "eps",    numtoMx(eps));
      mxAddField2Scalar(FL, I, "s2t",    numtoMx(s2t));
      mxAddField2Scalar(FL, I, "info",   mxCreateString(str));
   }
}


char mxGetLRDir(
   const char *F, int L,
   const mxArray *a, unsigned k"k-th argument" (for log purposes)
){
   char ldir=0, e=0;

   if (mxIsChar(a)) {
      if (mxGetString(a,str)) wblog(F,L,
      "ERR failed to read lrdir string (arg #%d)",k);

      if (!strcmp(str,"LR") || !strcmp(str,"12") || !strcmp(str,">>"))
      ldir=+1; else

      if (!strcmp(str,"RL") || !strcmp(str,"21") || !strcmp(str,"<<"))
      ldir=-1; else


      wblog(F,L,
      "ERR argument #%d >%s< must be of type `LR' or `RL'",k,str);
   }
   else 
   if (mxIsNumber(a)) {
      if (mxGetNumber(a,ldir)) e++;
      else {
         if (ldir==12 || ldir==+1) ldir=+1; else
         if (ldir==21 || ldir==-1) ldir=-1; else e++;
      }

      if (e) wblog(F,L,
      "ERR argument #%d must be of type `LR' or `RL' (%d).",k,ldir);
   }
   else wblog(F,L,"ERR invalid argument #%d (lrdir)",k); 

   return ldir;
}


template <class TQ, class T1, class T2, class T3>
void twoSiteInit(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &PSI, unsigned &K,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
){
   const unsigned r1=Psi1.rank(), r2=Psi2.rank();
   wbindex I1, I2;
   wbperm P;

   if (ldir!=+1 && ldir!=-1) wblog(FL,"WRN ldir=%c<%d> ???",ldir,ldir);
   if ((ldir>0 && ic1>=r1) || (ldir<0 && ic2>=r2)) wblog(FL,
   "ERR contraction index out of bounds (%d/%d, %d/%d)",ic1,r1,ic2,r2);
   
   if (twoSite) {
      if (ldir>0) {
         Psi1.contract(ic1+1, Psi2, ic2+1, PSI);
         I1.Index(r1-1); I2.Index(r2-1)+=(r1-1);
      }
      else {
         Psi2.contract(ic2+1, Psi1, ic1+1, PSI);
         I1.Index(r2-1); I2.Index(r1-1)+=(r2-1);
      }
   }
   else {
      if (ldir>0)
           { PSI=Psi1; I1.Index(r1).Skip(ic1); I2.init(1,&ic1); }
      else { PSI=Psi2; I1.Index(r2).Skip(ic2); I2.init(1,&ic2); }
   }

   if (PSI.isEmpty()) wblog(FL,
   "WRN %s yields empty QSpace !??", __FUNCTION__);

   P.Cat(I1,I2);
   PSI.Permute(P); K=I1.len;
}


template <class TQ, class T1, class T2, class T3>
void twoSiteFinal(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &A1,
   QSpace<TQ,T3> &A2,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
){
   if (ldir!=+1 && ldir!=-1) wblog(FL,"WRN ldir=%d ???", ldir);
   QSpace<TQ,T3> AX;

   if (twoSite) {
      if (ldir>0) {
         A1.PermuteLastTo (ic1);
         A2.PermuteFirstTo(ic2);
      }
      else {
         A2.permuteFirstTo(ic1, AX);
         A1.permuteLastTo (ic2, A2); AX.save2(A1);
      }
   }
   else {
      if (ldir>0) {
         A1.PermuteLastTo(ic1);
         A2.contract(2, Psi2, ic2+1, AX);
         AX.permuteFirstTo(ic2, A2);
      }
      else {
         A2.contract(2, Psi1, ic1+1, AX);
         A1.permuteLastTo(ic2, A2);
         AX.permuteFirstTo(ic1,A1);
      }
   }
}


template <class TQ, class TD>
mxArray* mpsOrthoQS(
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &A1,
   QSpace<TQ,TD> &A2,
   unsigned K, unsigned Nkeep, double stol,
   const char *info
){
   SVDData<TQ,TD> SVD;
   wbvector<double> SV;

   mxArray *S=NULL;
   char dflag=0, sflag=0, Sflag=0;
   unsigned i;

   if (info && info[0]) {
      wbstring istr(info);

      dflag=istr.count('d');
      sflag=istr.count('s');
      Sflag=istr.count('S');

      if (dflag>1 || sflag>1 || Sflag>1 ||
          dflag+sflag+Sflag!=(int)istr.length())
      wblog(FL,"ERR Invalid flag set >%s<", info);
   }

   if (sflag || Sflag)
   S=mxCreateStructMatrix(1,1,0,NULL);

   SVD.blockSVD(PSI,K);
   i=SVD.dmrgTruncate(Nkeep, stol, SV, A1,A2, S, dflag);


   if (S) {
      mxAddField2Scalar(FL, S, "K", numtoMx(K));
      mxAddField2Scalar(FL, S, "D1", A1.getDim().toMx());
      mxAddField2Scalar(FL, S, "DD",PSI.getDim().toMx());
   }
   else if (i) wblog(FL,str);

   if (Sflag)
   mxAddField2Scalar(FL,S,"SVD", SVD.toMx());

   return S;
}


template <class TQ, class TD>
mxArray* mpsGetSVD(
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &U,
   QSpace<TQ,double> &S,
   QSpace<TQ,TD> &VT,
   unsigned K,
   const char *info
){
   SVDData<TQ,TD> SVD;

   mxArray *I=NULL;
   char sflag=0, Sflag=0;

   if (info && info[0]) {
      wbstring istr(info);

      sflag=istr.count('s');
      Sflag=istr.count('S');

      if (sflag>1 || Sflag>1 || sflag+Sflag!=(int)istr.length())
      wblog(FL,"ERR Invalid flag set >%s<", info);
   }

   if (sflag || Sflag)
   I=mxCreateStructMatrix(1,1,0,NULL);

   SVD.blockSVD(PSI,K);
   SVD.getQSpaceUSV(U,S,VT,I);

   if (I) {
      mxAddField2Scalar(FL,I,"K", numtoMx(K));
      mxAddField2Scalar(FL,I,"D1",  U.getDim().toMx());
      mxAddField2Scalar(FL,I,"DD",PSI.getDim().toMx());
   }

   if (Sflag)
   mxAddField2Scalar(FL,I,"SVD", SVD.toMx());

   return I;
}


template <class TQ, class TD>
mxArray* SVDData<TQ,TD>::toMx() const {

   const char *fn[]={"Q1","Q2","S1","S2","U","S","V","AS","I0","Ik"};

   mxArray *S=mxCreateStructMatrix(QB.len,1,10,fn);

   for (unsigned i=0; i<QB.len; i++) {
      mxSetFieldByNumber(S,i,0, QB[i].Q1.toMx() );
      mxSetFieldByNumber(S,i,1, QB[i].Q2.toMx() );
      mxSetFieldByNumber(S,i,2, QB[i].S1.toMx() );
      mxSetFieldByNumber(S,i,3, QB[i].S2.toMx() );
      mxSetFieldByNumber(S,i,4, QB[i].U .toMx() );
      mxSetFieldByNumber(S,i,5, QB[i].S .toMx() );
      mxSetFieldByNumber(S,i,6, QB[i].V .toMx() );
      mxSetFieldByNumber(S,i,7, QB[i].AS.toMx() );
      mxSetFieldByNumber(S,i,8, QB[i].I0.toMx() );
      mxSetFieldByNumber(S,i,9, QB[i].Ik.toMx() );
   }

   return S;
}


#endif

