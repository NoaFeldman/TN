#ifndef __WB_QSPACE_AUX_CC__
#define __WB_QSPACE_AUX_CC__

//====================================================================//
// auxilliary map for QSpace => CGC contraction
//====================================================================//

// NB! use "#pragma omp critical"
// when accessing matlab workspace (eg. through MXPut())
// after catching an exception; e.g. see ~/Source/MPI/tst_openmp.cc
// Wb,Aug06,14

template <class TQ>
template <class TA, class TB>
unsigned OrthoCGS<TQ>::contractCGS_ortho(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   unsigned isym
){
   unsigned ra=A.rank(F_L), rb=B.rank(F_L), rc=ra+rb-2*ica.len;

   if (A.qtype!=B.qtype || A.QDIM!=B.QDIM) wblog(F_L,
      "ERR %s() qtype mismatch (%s <> %s)",
      FCT,A.qStr().data,B.qStr().data);
   if (ica.len!=icb.len || ica.len>ra || icb.len>rb) wblog(F_L,
      "ERR %s() invalid contraction indices (%d/%d; %d %d)",
      FCT,ica.len,icb.len,ra,rb);
   if (Ia.len!=Ib.len || !Ia.len) wblog(F_L,"ERR %s()\n"
      "invalid input data (%d/%d)", FCT, Ia.len, Ib.len);
   if (isym>=A.qtype.len) wblog(F_L,"ERR %s() "
      "symmetry index out bounds (%d/%d)",FCT,isym+1,A.qtype.len);
   if (A.CGR.isEmpty() ^ B.CGR.isEmpty()) wblog(FL,
      "ERR %s() got mixed empty/non-empty CGR data (%s <> %s)",
      FCT, A.CGR.sizeStr().data, B.CGR.sizeStr().data);

   unsigned i, N=Ia.len; int gotC=(A.CGR.data ? 1 : 0), e=0;

   if (gotC) {
   if (!A.CGR(Ia[0],isym).cgr || !B.CGR(Ib[0],isym).cgr) { gotC=-1;
      for (i=0; i<N; ++i) {
      if (!A.CGR(Ia[i],isym).isAbelian() ||
          !B.CGR(Ib[i],isym).isAbelian()) { e=gotC; break; }
      }
   }
   else if (A.qtype[isym].isAbelian()) {
      for (i=0; i<N; ++i) {
      if (!A.CGR(Ia[i],isym).isAbelian() || 
          !B.CGR(Ib[i],isym).isAbelian()) { e=gotC; break; }
      }
   }
   else { gotC=2;
      for (i=0; i<N; ++i) {
      if (A.CGR(Ia[i],isym).isAbelian() || 
          B.CGR(Ib[i],isym).isAbelian()) { e=gotC; break; }
      }
   }
   if (e) wblog(FL,
      "ERR %s() got mixed CGR_ABELIAN data\n(%d: 0x%lX, 0x%lX; e=%d)",
      FCT,i+1, A.CGR(Ia[i],isym).cgr, B.CGR(Ib[i],isym).cgr,e
   ); }

   if (A.CGR.dim2!=A.qtype.len || B.CGR.dim2!=B.qtype.len) wblog(F_L,
      "ERR %s() invalid CGR data (nsym=%d,%d/%d,%d)",
      FCT,A.CGR.dim2,B.CGR.dim2,A.qtype.len,B.qtype.len
   );

   wbvector<unsigned> qdc;

   A.qtype.Qpos(qdc);

   const TQ *qA=A.QIDX.data+qdc[isym], *qB=B.QIDX.data+qdc[isym];
   unsigned m=0, n=A.qtype[isym].qlen(), dA=A.QIDX.dim2, dB=B.QIDX.dim2;
   wbvector< CRef<TQ> > cc(N);

   Q.init(F_L,
      A.CGR(Ia[0],isym), ica,
      B.CGR(Ib[0],isym), icb
   );

   rtype=CGR_DEFAULT;

   if (!Q.isEmpty() && Q.rank()!=rc) wblog(FL,
      "ERR %s() invalid QSet r=%d/%d !??",FCT,Q.rank(),rc);
   else Q.t=A.qtype[isym];

#ifndef NSAFEGUARDS
   QSet<TQ> qsA_(Q.t,A.itags,ra), qsB_(Q.t,B.itags,rb), qsA, qsB, qsC;


   for (i=0; i<N; ++i) {
      const CRef<TQ> &a=A.CGR(Ia[i],isym), &b=B.CGR(Ib[i],isym);
      
      qsA=qsA_; qsA.qs.initStride( qA+Ia[i]*dA, n,ra,A.QDIM);
      qsB=qsB_; qsB.qs.initStride( qB+Ib[i]*dB, n,rb,B.QDIM);
      if (a.cgr) { a.adapt(qsA); a.cgr->checkQ(F_L,qsA); }
      if (b.cgr) { b.adapt(qsB); b.cgr->checkQ(F_L,qsB); }

      if (i) {
         qsC.init(F_L,a,ica,b,icb);
         if (qsC!=Q) wblog(F_L,"ERR %s() "
            "got non-unique output symmetry sector\n%d: %s <> 1: %s",
            FCT,i+1,qsC.toStr().data, Q.toStr().data
         );
      }
      else if (gotC<0) { Q=qsC; }
   }
#endif

   if (gotC<=0) {
      Q.init(A.qtype[isym]); cgr=NULL;
      cgQ.init(1,1); cgQ[0]=1;
      cgR.init(1,Ia.len).set(RTD(1));

      cgr=NULL; rtype=CGR_ABELIAN;
      cgp.init(); conj=0;

      return 1;
   }


      for (i=0; i<N; ++i) {
         gCX.contract(FL,
            A.CGR(Ia[i],isym), ica,
            B.CGR(Ib[i],isym), icb, cc[i]
         );


         if (cc[i].rtype!=CGR_DEFAULT) {
            if (cc[i].rtype==CGR_CTR_ZERO) { if (cc[i].cgw.len) e=1; } else
            if (cc[i].rtype==CGR_CTR_SCALAR) { if (cc[i].cgw.len!=1) e=2; }
            else e=3;
         } else if (!cc[i].cgw.len) { e=4; }
         if (e) wblog(FL,"ERR %s() got %s",FCT,cc[i].toStr().data);

         if (i) { if (rtype!=cc[i].rtype) wblog(FL,
            "ERR %s() got inconsistent contraction type (%d/%d)\n%s",
            rtype.toStr().data, cc[i].rtype.toStr().data, cc[i].toStr().data);
         }
         else rtype=cc[i].rtype;

         if (m<cc[i].cgw.len) m=cc[i].cgw.len;

         if (i) {
            if (cc[i].cgp!=cc[0].cgp || cc[i].conj!=cc[0].conj)
            wblog(FL,"ERR %s() got different CRef (%d):\n"
               "cgp =%s <> %s; conj=%d/%d !?", FCT,i+1,
               cc[i].cgp.toStr().data, cc[0].cgp.toStr().data,
               cc[i].conj, cc[0].conj
            );
         }

      }

   cgr =(CDATA_TQ*)(cc[0].cgr);
   cgp =cc[0].cgp;
   conj=cc[0].conj;

   wbarray<RTD> cgw(m,N);
   if (m) { for (i=0; i<N; ++i) {
      MEM_CPY<RTD>(cgw.col(i), cc[i].cgw.len, cc[i].cgw.data);
   }}


   unsigned rval=cgw.QRdecomp(FL,cgQ,cgR);

   cgQ.SkipTiny(CG_SKIP_EPS2);
   cgR.SkipTiny(CG_SKIP_EPS2);

   if (rc<3) {
      RTD cfac=cc[0].normExt();
      if (cfac!=1) { cgQ*=cfac; cgR*=(RTD(1)/cfac); }
   }

#ifndef NSAFEGUARDS
   if (rval>m) wblog(FL,"ERR rval=%d/%d !?",rval,m);

#endif

   return rval;
};


template <class TQ>
template <class TD>
unsigned OrthoCGS<TQ>::getCGS_ortho(
   const char *F, int L,
   QSpace<TQ,TD> &A, C_UVEC &Ia, unsigned isym
){
   unsigned ra=A.rank(F_L);

   if (A.gotCGS(F_L)<=1) wblog(FL,
      "WRN %s() not expecting OM (%d)",FCT,A.gotCGS(F_L));

   if (!Ia.len) wblog(FL,"ERR %s() got empty Ia !??",FCT);
   if (A.CGR.dim2!=A.qtype.len) wblog(F_L,
      "ERR %s() invalid CGR data (%d/%d)",FCT,A.CGR.dim2,A.qtype.len);
   if (isym>=A.qtype.len) wblog(F_L,"ERR %s() "
      "symmetry index out bounds (%d/%d)",FCT,isym+1,A.qtype.len); 

   unsigned i, m=0, N=Ia.len;
   wbvector<unsigned> qdc;

   A.qtype.Qpos(qdc);

   const TQ *qA=A.QIDX.data+qdc[isym];
   unsigned n=A.qtype[isym].qlen(), dA=A.QIDX.dim2;
   CDATA_TQ X;

   const CRef<TQ> &R=A.CGR(Ia[0],isym);
   init(R);

   QSet<TQ> qsA(A.qtype[isym],A.itags,ra);
   for (i=0; i<N; ++i) {
      const CRef<TQ> &ar=A.CGR(Ia[i],isym);
      qsA.qs.initStride(qA+Ia[i]*dA, n, ra, A.QDIM);

      if (!ar.cgr) {
         if (ar.isAbelian()) { if (!m) m=1; }
         else wblog(FL,"ERR %s() got %s",FCT,ar.toStr().data);
      }
#ifndef NSAFEGUARDS
      else if (!ar.sameQSet(qsA)) {
         MXPut X(FL,"a"); X.add(qsA,"Q").add(ar,"R").add(i+1,"i");
            if (ar.cgr) { X.add(*ar.cgr,"R"); }; X.put();
         wblog(FL,"ERR %s(%d/%d,%d)) got QSet mismatch\n%s\n%s",
         FCT,i+1,N,isym+1, ar.toStr().data, qsA.toStr().data);
      }
#endif
      if (i) { if (qsA!=Q) wblog(F_L,"ERR %s() "
         "got non-unique symmetry sector\n%2d: %s\n 1: %s",
         FCT,i+1,qsA.toStr().data, Q.toStr().data
      ); }
      else Q=qsA;

      if (m<ar.cgw.len) m=ar.cgw.len;
   }

   if (!m || (R.cgr && m>(R.cgr->getOM())) || (!R.cgr && m>1)) wblog(FL,
      "ERR %s() OM out of bounds (%d/%d) !?",FCT,m,R.getOM());

   wbarray<RTD> cgw(m,N);
   for (i=0; i<N; ++i) {
      const CRef<TQ> &ar=A.CGR(Ia[i],isym);
      if (ar.cgr) {
         if (!ar.cgw.len) wblog(FL,"ERR %s() %s",FCT,ar.toStr().data);
         MEM_CPY<RTD>(cgw.col(i), ar.cgw.len, ar.cgw.data);
      }
      else {
         if (!ar.isAbelian()) wblog(FL,"ERR %s()",FCT,ar.toStr().data);
         cgw(0,i)=1;
      }
   }

   unsigned rval=cgw.QRdecomp(FL,cgQ,cgR);

   double e2 =
      cgQ.SkipTiny(CG_SKIP_EPS2)
    + cgR.SkipTiny(CG_SKIP_EPS2);

   if (CG_VERBOSE>1 && e2>1E-4*CG_SKIP_EPS2) {
      double qmin=cgQ.aMin(1);
      double rmin=cgR.aMin(1); wblog(FL,
     "TST %s() skipping %.3g (QR @ %.3g, %.3g)",FCT,SQRT(e2),qmin,rmin);
      MXPut(FL,"Iqr").add(cgQ,"Q").add(cgR,"R")
       .add(e2,"e2").add(qmin,"qmin").add(rmin,"rmin"); 
   }


   return rval;
};


template <class TQ, class TA, class TB, class TC>
int contractDATA_ortho(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   const wbvector< OrthoCGS<TQ> > &MC,
   QSpace<TQ,TC> &C, unsigned ic
){
   wbvector<unsigned> S(MC.len);
   unsigned i, ra=A.rank(F_L), rb=B.rank(F_L), rc=ra+rb-2*ica.len;

   if (MC.len!=A.qtype.len || A.qtype.len!=B.qtype.len) wblog(FL,
      "ERR %s() size mistmatch (%d/%d/%d)",
      FCT, MC.len, A.qtype.len, B.qtype.len);
   if (Ia.len!=Ib.len || !Ia.len) wblog(F_L,"ERR %s()\n"
      "invalid input data (%d/%d)", FCT, Ia.len, Ib.len);
   if (C.CGR.dim1!=C.DATA.len || C.CGR.dim2!=C.qtype.len) wblog(F_L,
      "ERR %s() invalid CGR data (C: %dx%d <> %dx%d)",
      FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len);
   if (rc!=C.rank(F_L)) wblog(FL,
      "ERR %s() inconsistent rc=%d/%d !??",FCT,rc,C.rank(F_L));

   for (i=0; i<MC.len; ++i) {
      if (!MC[i].cgQ.isMatrix()) wblog(FL,
         "ERR %s() got S=%s !??",FCT,MC[i].cgQ.sizeStr().data);
      S[i]=MC[i].cgQ.SIZE[1];
   }

   wbIndex I(S);
   wbarray<TC> Ck;
   unsigned j,k,l=-1, nx=0; INDEX_T *idx=I.data, m=S.prod();
   RTD x, one=1, eps=1E-99;

   if (!m) wblog(FL,"ERR %s() got total OM=%d !??",FCT,m);
   if (ic+m>C.DATA.len) wblog(FL,
      "ERR %s() C.DATA out of bounds (%d/%d)",FCT,ic,m,C.DATA.len);

   while (++I) { ++l;
      for (j=0; j<I.len; ++j) {
         CRef<TQ> &cj = C.CGR(ic+l,j);
         cj.cgr   = MC[j].cgr;
         cj.cgp   = MC[j].cgp;
         cj.conj  = MC[j].conj;
         cj.rtype = MC[j].rtype;

         if (cj.cgr || cj.rtype==CGR_CTR_SCALAR) {
            MC[j].cgQ.getCol(idx[j], cj.cgw);
         }
         else if (cj.rtype!=CGR_ABELIAN || cj.cgp.len || cj.cgw.len)
            wblog(FL,"ERR %s() invalid abelian CRef\n%s",
            FCT,cj.toStr().data
         );

         cj.SortDegQ(F_L);

         if (cj.cgr) cj.cgr->checkValidStat(FL);

      }

      C.DATA[ic+l]->init();
   }

   if ((++l)!=m) wblog(FL,"ERR %s() %d/%d !??",FCT,l,m);

   try {
      for (k=0; k<Ia.len; ++k) { I.reset(); l=-1;

         if (k) Ck.init();
         A.DATA.el(Ia[k])->contract(FL,ica, *B.DATA.el(Ib[k]),icb, Ck);

         while (++I) { ++l;
            for (x=one, j=0; j<MC.len; ++j) {
               x*=MC[j].cgR(idx[j],k); if (ABS(x)<eps) break;
            }
            if (j<MC.len) { ++nx;
               continue;
            }
            C.DATA.el(ic+l)->Plus(Ck,TC(x),'i');
         }
      }
   }
   catch (...) {
   #pragma omp critical
    { unsigned ia=Ia[k], ib=Ib[k];
      const wbarray<TA> &a=(*A.DATA[ia]);
      const wbarray<TB> &b=(*B.DATA[ib]);
      const wbarray<TC> &c=(*C.DATA[ic]);

      MXPut S(FL); S.add(A,"A").add(B,"B")
         .add(a,"a").add(ia+1,"ia").add(ica+1,"ica")
         .add(b,"b").add(ib+1,"ib").add(icb+1,"icb").add(Ck,"c");

      wblog(F,L,"ERR %s::%s\n"
         "A(%2d) %6s [ %s ] @ [%s]\n"
         "B(%2d) %6s [ %s ] @ [%s]\n"
         "C(%2d) %6s [ %s ]", shortFL(FL),FCT,
      ia+1, a.sizeStr().data, A.qrec2Str(ia).data, (ica+1).toStr().data,
      ib+1, b.sizeStr().data, B.qrec2Str(ib).data, (icb+1).toStr().data,
      ic+1, c.sizeStr().data, C.qrec2Str(ic).data);
   }}

   return nx;
};


template <class TQ, class TD>
unsigned combineDATA_ortho(
   const char *F, int L,
   const QSpace<TQ,TD> &A, C_UVEC &Ia,
   const wbvector< OrthoCGS<TQ> > &MC, QSpace<TQ,TD> &C, unsigned ic
){
   wbvector<unsigned> S(MC.len);
   unsigned i,j,k, ra=A.rank(F_L);

   if (!Ia.len) wblog(F_L,"ERR %s() got empty ia",FCT);
   if (MC.len!=A.qtype.len || A.qtype!=C.qtype) wblog(FL,
      "ERR %s() size/qtype mistmatch\n(%d/%d; %s <> %s)",
      FCT, MC.len, A.qtype.len, A.qStr().data, C.qStr().data);
   if (C.CGR.dim1!=C.DATA.len || C.CGR.dim2!=C.qtype.len) wblog(F_L,
      "ERR %s() invalid CGR data (C: %dx%d <> %dx%d)",
      FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len);
   if (ra!=C.rank(F_L)) wblog(FL,
      "ERR %s() inconsistent r=%d/%d !??",FCT,ra,C.rank(F_L));

   for (i=0; i<MC.len; ++i) {
      const wbvector<unsigned> &Sq=MC[i].cgQ.SIZE, &Sr=MC[i].cgR.SIZE;

      if (Sq.len!=2 || Sr.len!=2 || Sq[0]<Sq[1] || Sq[1]!=Sr[0])
         wblog(FL,"ERR %s() got QR = %s * %s !?",FCT,
         MC[i].cgQ.sizeStr().data, MC[i].cgR.sizeStr().data
      );
      S[i]=Sq[1];
   }

   wbIndex I(S);
   INDEX_T l=-1, m=S.prod(), *idx=I.data;
   wbvector<unsigned> mark(m);
   RTD x, one=1, eps=1E-99;

   if (!m) wblog(FL,"ERR %s() got total OM=%d !??",FCT,m);
   if (ic+m>C.DATA.len) wblog(FL,
      "ERR %s() C.DATA out of bounds (%d/%d)",FCT,ic,m,C.DATA.len);

   while (++I) { ++l;
      for (j=0; j<I.len; ++j) {
         CRef<TQ> &cj = C.CGR(ic+l,j);

         cj.cgr   = MC[j].cgr;
         cj.cgp   = MC[j].cgp;
         cj.conj  = MC[j].conj;
         cj.rtype = MC[j].rtype;

         if (cj.cgr) {
            MC[j].cgQ.getCol(idx[j], cj.cgw);
         }
         else if (cj.rtype!=CGR_ABELIAN || cj.cgp.len || cj.cgw.len)
            wblog(FL,"ERR %s() invalid abelian CRef\n%s",
            FCT,cj.toStr().data
         );
      }

      C.DATA[ic+l]->init();
   }

   if ((++l)!=m) wblog(FL,"ERR %s() %d/%d !??",FCT,l,m);

   try {
      for (k=0; k<Ia.len; ++k) { I.reset(); l=-1;
         const wbarray<TD> &Ak=*A.DATA.el(Ia[k]);
         while (++I) { ++l;
            for (x=one, j=0; j<MC.len; ++j) {
               x*=MC[j].cgR(idx[j],k); if (ABS(x)<eps) break;
            }
            if (j==MC.len) { ++mark[l];
               C.DATA.el(ic+l)->Plus(Ak,TD(x),'i');
            }
         }
      }
   }
   catch (...) {
   #pragma omp critical
    { unsigned ia=Ia[k];
      MXPut S(FL); S.add(A,"A").add(ia+1,"ia")
         .add(*A.DATA[ia],"a").add(*C.DATA[ic+l],"c");

      wblog(F,L,"ERR %s::%s\n"
         "A(%d): %6s  [ %s ]\n"
         "C(%d): %6s  [ %s ]", shortFL(FL), FCT,
      ia+1, A.DATA[ia]->sizeStr().data, A.qrec2Str(ia).data,
      ic+1, C.DATA[ic]->sizeStr().data, C.qrec2Str(ic).data);
   }}

   return mark.numZeros();
};


template <class TQ, class TA, class TB, class TC>
void contractDATA_plain(const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic
){

   wbvector<unsigned> S(A.qtype.len);
   unsigned i,j, ra=A.rank(F_L), rb=B.rank(F_L), rc=ra+rb-2*ica.len;

   if (A.qtype!=B.qtype || !A.isAbelian()) wblog(FL,
      "ERR %s() got non-abelian %s <> %s",
      FCT, A.qStr().data, B.qStr().data);
   if (Ia.len!=Ib.len || !Ia.len) wblog(F_L,"ERR %s()\n"
      "invalid input data (%d/%d)", FCT, Ia.len, Ib.len);
   if (C.QIDX.dim1!=C.DATA.len) wblog(F_L,
      "ERR %s() QSpace not yet setup !?? (C: %dx%d <> %dx%d)",
      FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len);
   if (rc!=C.rank(F_L)) wblog(FL,
      "ERR %s() inconsistent rc=%d/%d !??",FCT,rc,C.rank(F_L));

   if (ic>=C.DATA.len) wblog(FL,
      "ERR %s() C.DATA out of bounds (%d/%d)",FCT,ic,C.DATA.len);
   if (!C.DATA[ic]) wblog(FL,
      "ERR %s() got null DATA !??  (ic=%d/%d)",FCT,ic,C.DATA.len);

   if (!C.CGR.isEmpty()) {
      wblog(FL,"WRN %s() "
         "got non-empty CGR for qtype=%s",FCT,C.qtype.toStr().data);
      if (C.CGR.dim1!=C.DATA.len || C.CGR.dim2!=C.qtype.len) wblog(F_L,
         "ERR %s() invalid CGR data (C: %dx%d <> %dx%d)",
         FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len
      );
      for (j=0; j<C.CGR.dim2; ++j) {
         if (!C.CGR(ic,j).isAbelian()) wblog(FL,
        "ERR %s() got non-abelian cgref data !??",FCT);
      }
   }

   C.DATA[ic]->init();

   try {
      for (i=0; i<Ia.len; ++i) {
         A.DATA.el(Ia[i])->contract(
            FL,ica, *B.DATA.el(Ib[i]),icb, *C.DATA[ic]
         );
      }
   }
   catch (...) {
   #pragma omp critical
    { unsigned ia=Ia[i], ib=Ib[i];
      const wbarray<TA> &a=(*A.DATA[ia]);
      const wbarray<TB> &b=(*B.DATA[ib]);
      const wbarray<TC> &c=(*C.DATA[ic]);

      MXPut(FL,"icd")
         .add(A,"A").add(ia+1,"ia").add(a,"a").add(ica+1,"ica")
         .add(B,"B").add(ib+1,"ib").add(b,"b").add(icb+1,"icb")
         .add(C,"C");

      wblog(F,L,"ERR %s::%s\n"
         "A(%d): %6s  [ %s ] @ [%s]\n"
         "B(%d): %6s  [ %s ] @ [%s]\n"
         "C(%d): %6s  [ %s ]", shortFL(FL),FCT,
      ia+1, a.sizeStr().data, A.qrec2Str(ia).data, (ica+1).toStr().data,
      ib+1, b.sizeStr().data, B.qrec2Str(ib).data, (icb+1).toStr().data,
      ic+1, c.sizeStr().data, C.qrec2Str(ic).data);
   }}
};



template<class TQ, class TD>
void getQall(const wbvector< QSpace<TQ,TD> > &F, wbMatrix<TQ> &QA){

   unsigned i,i0,k=0,QDIM,dim2;

   wbvector < wbMatrix<TQ>* > pQ(F.len);

   for (i=0; i<F.len; ++i) if (!F[i].isEmpty()) break;
   if (i>=F.len) { QA.init(); return; }

   pQ[k++]=&(F[i].QIDX); QDIM=F[i].QDIM; dim2=F[i].QIDX.dim2;

   if (QDIM==0 || dim2==0 || dim2 % QDIM) wblog(FL,
   "ERR invalid operator (empty; %d/%d)",dim2,QDIM);

   for (i0=i, i=i0+1; i<F.len; ++i) { if (F[i].isEmpty()) continue;

      if (F[i].qtype!=F[i0].qtype) wblog(FL,
         "ERR qtype inconsistency %d,%d/%d: %s <> %s",i0+1,i+1,F.len,
         F[i].qStr().data, F[i0].qStr().data
      );

      if (F[i].QDIM!=QDIM || F[i].QIDX.dim2!=dim2)
      wblog(FL,"ERR dimensional inconsistency %d/%d: %d/%d <> %d/%d",
      i0+1, i+1, F.len, F[i].QIDX.dim2, F[i].QDIM, dim2, QDIM);

      pQ[k++]=&(F[i].QIDX);
   }

   pQ.len=k;

   QA.CAT(1,pQ);
   QA.Reshape(QA.dim1*QDIM, QA.dim2/QDIM);
   QA.makeUnique();
};


int mxIsQSpace(
   const char *F, int L,
   const mxArray *a,
   unsigned &rank,
   char cflag,
   unsigned k,
   const unsigned *rmin,
   const unsigned *rmax,
   const char *istr
){
   mxArray *Q, *D, *q;
   unsigned i,l,n,k2, M=0, N=0,
     r, r1=(rmin ? *rmin : rank), r2=(rmax ? *rmax : r1);

   int fidQ,fidD, isq=1;

   if (!a || !mxIsStruct(a) || !mxIsVector(a)) {
      if (F || istr) wblog(F_L,
         "ERR %sinvalid vector QSpace structure (%s)",
          istr ? wbstring(istr,"\n").data : "",
          mxTypeSize2Str(a).data);
      return 0;
   }

   n=mxGetNumberOfElements(a);
   if (int(k)>=int(n)) {
      if (F || istr) wblog(F_L,
         "ERR %sstructure index out of bounds (%s <> %d)",
          istr ? wbstring(istr,"\n").data : "",
          mxTypeSize2Str(a).data, int(k)+1);
      return 0;
   }

   if (int(k)==-1 && n>1) {
      if (F || istr) wblog(F_L,
         "ERR %sexpecting single QSpace (got %d)",
          istr ? wbstring(istr,"\n").data : "", n);
      return 0;
   }


   fidQ=mxGetFieldNumber(a,"Q");
   fidD=mxGetFieldNumber(a,"data");

   if (fidQ<0 || fidD<0) {
      if (F || istr) wblog(F_L,
         "ERR %smxArray is not of type {Q,data,...}",
          istr ? wbstring(istr,"\n").data : "");
      return 0;
   }

   if (int(k)<0) { k=0; k2=n; } else k2=k+1;
   for (; k<k2; ++k) {

      Q=mxGetFieldByNumber(a,k,fidQ); 
      D=mxGetFieldByNumber(a,k,fidD);

      if (!Q || !D) {
         if (!Q && !D) continue;
         if ((Q && !mxIsEmpty(Q)) || (D && !mxIsEmpty(D))) {
            if (F || istr) wblog(F_L,
               "ERR %sQSpace data only set partially (k=%d)",
                istr ? wbstring(istr,"\n").data : "", k+1);
            return 0;
         }
      }

      if ((!Q || mxIsEmpty(Q)) && (!D || mxIsEmpty(D))) continue;

      if (Q) { r=mxGetNumberOfElements(Q); } else { r=0; }
      if (r) {
         if (!mxIsVector(Q) || !mxIsCell(Q)) {
            if (F || istr) wblog(F_L,
               "ERR %sQ-field not a vector cell (%d: %s)",
                istr ? wbstring(istr,"\n").data : "",
                k+1, mxTypeSize2Str(Q).data);
            return 0;
         }
      }

      if (int(rank)<0) {
         rank=r; if (int(r1)<0) r1=rank; if (int(r2)<0) r2=r1;
         if (r1>r2) wblog(F_L,
         "ERR %s() invalid rank-range [%d %d]",FCT,r1,r2); 
      }
      else if (r<r1 || r>r2) {
         if (F || istr) wblog(F_L,
            "ERR %sinvalid rank (k=%d: r=%d not in [%d %d])",
             istr ? wbstring(istr,"\n").data : "",
             k+1,r,r1,r2);
         return 0;
      }

      for (i=0; i<r; ++i) {
         q=mxGetCell(Q,i);
         if (mxMatSafeGuard(0,L,q)) {
            if (F || istr) wblog(F_L,
               "ERR %s() invalid input data (%s%s%d)\n%s",
               FCT, istr ? istr : "", istr ? ": ": "",
               q ? int(mxGetNumberOfElements(q)) : -1, str);
            return 0;
         }
         if (i==0) { M=mxGetM(q); N=mxGetN(q); }
         if ((i && (M!=mxGetM(q) || N!=mxGetN(q)))
            || mxGetNumberOfDimensions(q)>2
          ){ if (F || istr) wblog(F_L,
               "ERR %sQ{%d} dimension mismatch (k=%d)",
                istr ? wbstring(istr,"\n").data : "", i+1, k+1);
             return 0;
         }
      }

      if (!D) {
         if (mxGetFieldNumber(a,"data")<0) {
            if (F || istr) wblog(F_L,
               "ERR %smissing field `data' (k=%d)",
                istr ? wbstring(istr,"\n").data : "", k+1);
            return 0;
         }
      }
      else {
         l=mxGetNumberOfElements(D);

         if (l!=M && (M || l!=1)) {
            if (F || istr) wblog(F_L,
               "ERR %sQSpace inconsistency data: %d/%s (%d)",
                istr ? wbstring(istr,"\n").data : "",
                M,mxSize2Str(D).data,k+1);
            return 0;
         }
         if (!l) continue;

         if (!mxIsVector(D) || !mxIsCell(D)) {
            if (F || istr) wblog(F_L,
               "ERR A(%d).data{} field not a cell vector (%s)",
                istr ? wbstring(istr,"\n").data : "",
                k+1,mxTypeSize2Str(D).data);
            return 0;
         }

         for (i=0; i<l; i++)
         if (mxMatSafeGuard(0,L, mxGetCell(D,i), r, cflag)) {
            if (F || istr) wblog(F_L,
               "ERR %sinvalid field A(%d).data{%d}\n%s",
                istr ? istr : "", k+1,i+1, str);
            return 0;
         }
      }

      if (Q && mxIsEmpty(Q) && D && mxIsEmpty(D)) isq=-1;
   }

   return isq;
};


bool mxIsEmptyQSpace(const mxArray *a, unsigned k) {

   mxArray *Q=NULL, *D=NULL;
   unsigned i1,i2;
   int fidQ, fidD;


   if (!a || mxIsEmpty(a)) return 1;

   if (k>=(unsigned)mxGetNumberOfElements(a)) wblog(FL,
      "ERR index out of bounds (%d/%d)",k,mxGetNumberOfElements(a));
   
   fidQ=mxGetFieldNumber(a,"Q");
   fidD=mxGetFieldNumber(a,"data");

   if (fidQ<0 || fidD<0) wblog(FL,
      "ERR mxArray is not of type {Q,data,...}", FL);

   Q=mxGetFieldByNumber(a,k,fidQ);  i1=(Q && !mxIsEmpty(Q));
   D=mxGetFieldByNumber(a,k,fidD);  i2=(D && !mxIsEmpty(D));

   if (i2 && !i1) {
      if (mxGetNumberOfElements(D)==1) {
         const mxArray *a=mxGetCell(D,0);
         if (mxGetNumberOfElements(a)==1) { return 0; }
      }
   }
   if (i1 && !i2) {
      unsigned i=0, n=mxGetNumberOfElements(Q);
      for (; i<n; ++i) {
         const mxArray *a=mxGetCell(Q,i);
         if (mxGetNumberOfElements(a)) { break; }
      }
      if (i==n) { return 1; }
   }

   if (i1 || i2) {
      if (i1 ^ i2) wblog(FL,
         "WRN Q or data set in QSpace, but not both !??");
      return 0;
   }

   return 1;
}


bool mxIsScalarQSpace(const mxArray *a, unsigned k) {

   int fidQ, fidD, fidI;
   mxArray *Q=NULL, *D=NULL, *I=NULL;

   if (!a || mxIsEmpty(a)) return 0;

   if (k>=(unsigned)mxGetNumberOfElements(a)) wblog(FL,
      "ERR index out of bounds (%d/%d)",k,mxGetNumberOfElements(a));
   
   fidQ=mxGetFieldNumber(a,"Q");
   fidD=mxGetFieldNumber(a,"data");
   fidI=mxGetFieldNumber(a,"info");

   if (fidQ<0 || fidD<0 || fidI<0) wblog(FL,
      "ERR mxArray is not of type {Q,data,info}", FL);

   Q=mxGetFieldByNumber(a,k,fidQ);
   if (Q && mxGetNumberOfElements(Q)) return 0;

   D=mxGetFieldByNumber(a,k,fidD);
   if (!mxIsCell(D) || mxGetNumberOfElements(D)!=1) return 0;
   else {
      const mxArray *a=mxGetCell(D,0);
      if (mxGetNumberOfElements(a)!=1) return 0;
   }

   I=mxGetFieldByNumber(a,k,fidI);
   if (I && (mxGetM(I)>1 || mxGetNumberOfDimensions(I)>2)) return 0;

   return 1;
};


bool mxIsQSpaceArr(
    const char *F, int L,
    const mxArray *S, unsigned rank, int arrdim,
    const unsigned *rmin,
    const unsigned *rmax,
    char cflag
){
    unsigned e=0, m=mxGetM(S), n=mxGetN(S);

    if (!m || !n) {
       if (F) wblog(F,L,"WRN %s received empty array (%dx%d)",FCT,m,n);
       return 1;
    }

    if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2) {
       if (F) wblog(F,L,
          "ERR %s() invalid structure array (%s)",
           FCT,mxTypeSize2Str(S).data);
       return 0;
    }

    switch (arrdim) {
       case 0: if (m!=1 || n!=1) e++; break;
       case 1: if (m!=1 && n!=1) e++; break;
    }
    if (e) {
        if (F) wblog(F,L,"ERR invalid QSpace %s (%dx%d)\n%s",
           arrdim==0 ? "scalar":"vector",m,n,str);
        return 0;
    }

    if (!mxIsQSpace(F,L,S,rank,cflag,-2,rmin,rmax)) {
        return 0;
    }

    return 1;
};


bool mxIsQSpaceVecVec(const mxArray *C, unsigned rank, char cflag) {

    unsigned i,m=0,n=0;
    const mxArray *a;

    if (C) { m=mxGetM(C); n=mxGetN(C); }

    if (!C || !m || !n) { 
       wblog(FL,"WRN got empty object (%s; %dx%d)\n--> %s",
       mxGetClassName(C),m,n, shortFL(FL));
       return 1;
    }

    if (!mxIsCell(C) || (m>1 && n>1)) {
       wblog(FL,"ERR got invalid object (%s)\n",mxTypeSize2Str(C).data);
       return 0;
    }

    n*=m;

    for (i=0; i<n; i++) { a=mxGetCell(C,i);
       mxIsQSpace(FL,a,rank,cflag,-2);
       return 0;
    }

    return 1;
}


bool mxsIsQSpaceVec(const mxArray *S, int fid, unsigned rank, char cflag) {

    unsigned m=mxGetM(S), n=mxGetN(S);
    mxArray *a;

    if ((m!=1 && n!=1) || !mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2) {
       sprintf(str, "%s:%d need vector structure (%dx%d)", FL,m,n);
       return 0;
    }

    n*=m;

    if (!n || fid<0) { sprintf(str,
      "%s:%d Invalid or empty data set (%d;%d)", FL, n, fid);
       return 0;
    }

    for (unsigned k=0; k<n; k++) {
       a=mxGetFieldByNumber(S,k,fid); if (!a || mxIsEmpty(a)) continue;
       mxIsQSpace(FL,a,rank,cflag,0);
    }

    return 1;
}


bool mxsIsQSpaceVEC(
    const mxArray *S, int fid,
    unsigned rank, char cflag
){
    unsigned i,k, d=0, m=mxGetM(S), n=mxGetN(S);
    mxArray *a;

    if ((m!=1 && n!=1) || !mxIsStruct(S) ||
        mxGetNumberOfDimensions(S)!=2) { sprintf(str,
      "%s:%d need vector structure of QSpace vectors \n"
      "containing {Q,data,...} structures (%d)", FL, m); return 0;
    }

    n*=m;

    if (!n || fid<0) { sprintf(str,
      "%s:%d Invalid or empty data set (%d;%d)", FL, n, fid);
       return 0;
    }

    for (i=0; i<n; i++) {
       a=mxGetFieldByNumber(S,i,fid); if (!a || mxIsEmpty(a)) continue;

       m=mxGetM(a); k=mxGetN(a);

       if (d==0) d=m*k;

       if ((m!=1 && k!=1) || n!=d || mxGetNumberOfDimensions(a)!=2 ||
          !mxIsCell(a)) { sprintf(str,
          "%s:%d field must be %d-dim QSpace vector structure (%dx%d)",
           FL, d, m, k); return 0;
       }

       for (k=0; k<d; k++) {
           mxIsQSpace(FL,a,rank,cflag,k);
           return 0;
       }
    }

    return 1;
}


int mxIsQSpaceVec(
   const char *F, int L,
   const char *fname0, const char *vname,
   unsigned rank,
   const unsigned *rmin,
   const unsigned *rmax,
   char cflag,
   unsigned N
){
    unsigned i, imax=1000, r1=-1,r2=-1, Nflag=0;
    mxArray *a;

    if (int(N)>0) {
       if (N>imax) wblog(FL,"ERR %s() expected chain length N=%d !??",FCT,N);
       Nflag=2;
    }

    wbstring fname(strlen(fname0)+12);

    for (i=0; i<imax; ++i) {
       if (Nflag && i>=Nflag) { i=N-Nflag; Nflag=imax; }
       snprintf(fname.data,fname.len,"%s_%02d.mat", fname0, i);

       Wb::matFile fid;
       if (!fid.open(FL,fname.data,"r")) break;

       printf("\r !X! "); fflush(0);

       a=matGetVariableInfo(fid.mfp,vname);

       printf("\r     %s|%s|%s(%s%s) %d ...   \r",
          shortFL(F_L), shortFL(FL), FCT,
          shortFL(fname0 ? fname0:"<data>"), (vname ? vname:""), i
       ); fflush(0);

       if (!a || mxIsEmpty(a)) continue;

       if (!mxIsQSpace(F,L,a,rank,cflag)) return -1;

       if (int(r1)<0 || int(r2)<0) {
          if (int(r1)<0) { r1=(rmin ? *rmin : rank); }
          if (int(r2)<0) { r2=(rmax ? *rmax : r1); }
          if (r1>r2) wblog(FL,
             "ERR %s() invalid rank-range [%d %d]",FCT,r1,r2); 
       }

       if (rank<r1 || rank>r2) {
          if (F) wblog(F,L,
             "ERR invalid QSpaceVec (%d: rank=%d [%d..%d])\n%s '%s'",
              i,rank,r1,r2,shortFL(FL),str);
          return -1;
       }
       if (r1!=r2) rank=-1;

       mxDestroyArray(a);
    }
    printf("\r%85s\r","");

    if (Nflag && i!=N) wblog(FL,
       "ERR  inconsistent NRG data (chain length %d/%d) !??", i,N);
    if (i>=imax) wblog(FL,"ERR more than %d files !??", imax);

    return i;
};


int mxIsQSpaceVEC(
   const char *F, int L,
   const char *fname0, const char *vname,
   unsigned rank,
   const unsigned *rmin,
   const unsigned *rmax,
   char cflag
){
    unsigned i,k,m,n,d=0,imax=1000, r1=-1,r2=-1;
    mxArray *a;

    wbstring fname(strlen(fname0)+12);

    for (i=0; i<imax; i++) {
       snprintf(fname.data,fname.len,"%s_%02d.mat", fname0, i);

       Wb::matFile fid;
       if (!fid.open(F_L,fname.data,"r")) break;

       a=matGetVariableInfo(fid.mfp,vname);

       if (!a || mxIsEmpty(a)) continue;

       m=mxGetM(a); n=mxGetN(a);

       if (d==0) d=m*n;

       if ((m!=1 && n!=1) || n!=d || !mxIsCell(a) ||
           mxGetNumberOfDimensions(a)!=2) { 
           if (F) wblog(F,L,"ERR sequence of %d-dim vectors of QSpace\n"
              "structures required (%dx%d; %d)",d,m,n,i);
           return -1;
       }

       for (k=0; k<d; k++) {
          if (!mxIsQSpace(F,L,a,rank,cflag,k)) return -1;
          if (i==0 && k==0) {
             r1=(rmin ? *rmin : rank);
             r2=(rmax ? *rmax : r1);
             if (r1>r2) wblog(F,L,
                "ERR %s() invalid rank-range [%d %d]",FCT,r1,r2); 
          }

          if (rank<r1 || rank>r2) {
             if (F) wblog(F,L,
                "ERR invalid QSpaceVec (%d: rank=%d [%d..%d])\n%s",
                 k+1,rank,r1,r2,str);
             return -1;
          }
          if (r1!=r2) rank=-1;
       }

       mxDestroyArray(a);
    }

    if (i>=imax) wblog(F,L,
    "WRN More than %d files - skip !??", imax);

    return i;
}


template<class TQ, class TD>
void mxInitQSpaceVec(const char *F, int L,
   const mxArray* S,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref,
   const unsigned *rmin=NULL,
   const unsigned *rmax=NULL
){
   unsigned k,m,n,l=0;

   if (!S) { FF.init(); return; }

   m=mxGetM(S); n=mxGetN(S);
   if (m==0 && n==0) { FF.init(); return; }

   if (!mxIsQSpaceVec(F_L,S,-1,rmin,rmax)) return;

   if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2)
       wblog(F,L,"ERR %s needs QSpace vector on input\n(%d, %dx%d)",
       FCT, mxIsStruct(S), m, n);

   n*=m;

   if (FF.len!=n) FF.initDef(n);

   for (k=0; k<n; k++) {
      FF[k].init(F,L,S,ref,k);
      if (l) FF[k].checkQ(F,L,FF[l-1]); else
      if (!FF[k].isEmpty()) { l=k+1; }
   }
}

template<class TQ, class TD>
void mxInitQSpaceVecR23(const char *F, int L,
   const mxArray* S,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref=0
){
   unsigned rmin=2, rmax=3;
   mxInitQSpaceVec(F,L,S,FF,ref,&rmin,&rmax);
};


template<class TQ, class TD>
void mxcInitQSpaceVec(const char *F, int L,
   const mxArray* C,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref
){
   unsigned i,m,n;

   if (!C) { FF.init(); return; }

   m=mxGetM(C); n=mxGetN(C);

   if (!mxIsCell(C) || mxGetNumberOfDimensions(C)!=2 || (m!=1 && n!=1))
   wblog(F,L,"ERR %s needs QSpace vector on input",FCT);

   n*=m;

   if (FF.len!=n) FF.initDef(n);

   for (i=0; i<n; i++)
   FF[i].init(F,L,mxGetCell(C,i),ref);
}


template<class TQ, class TD>
void mxInitQSpaceVecVec(const char *F, int L,
   const mxArray* C,
   wbvector< wbvector< QSpace<TQ,TD> > > &FF,
   const char ref
){
   unsigned i,m=0,n=0;
   const mxArray *a;

   if (C) { m=mxGetM(C); n=mxGetN(C); }

   if (!C || m==0 || n==0) {
   FF.init(); return; }

   if (!mxIsCell(C) || (m!=1 && n!=1)) {
      wblog(F,L,"ERR invalid QSpace cell vector: got %s(%d,%d)\n"
      "--> %s:%d", mxGetClassName(C),m,n, basename(__FILE__),__LINE__);
   }

   n*=m; if (FF.len!=n) FF.initDef(n);

   for (i=0; i<n; i++) { a=mxGetCell(C,i);
       if (mxIsQSpaceVec(FL,a))
          mxInitQSpaceVec(F,L,a,FF[i],ref);
       else {
          wbstring istr(str); sprintf(str,
         "Invalid QSpace cell vector (%d: %s)",i+1,mxGetClassName(C));
          wblog(F,L,"%s\n%s",str,istr.data);
       }
   }
}


template<class TQ, class TD>
void mxInitQSpaceMat(const char *F, int L,
   const mxArray* S,
   wbMatrix< QSpace<TQ,TD> > &FF,
   const char ref
){
   unsigned i,j,m,n;

   if (!S) { FF.init(); return; }

   if (mxIsCell(S)) { mxcInitQSpaceVec(F,L,S,FF,ref); return; }

   m=mxGetM(S); n=mxGetN(S);
   if (m==0 && n==0) { FF.init(); return; }

   if (!mxIsQSpaceMat(FL,S)) { wbstring istr(str);
      sprintf(str,"Invalid QSpace (%s)", mxGetClassName(S));
      wblog(F,L,"%s\n%s",str,istr.data);
   }

   if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2)
       wblog(F,L,"ERR %s needs QSpace array on input\n(%d, %dx%d)",
       FCT, mxIsStruct(S), m, n);

   
   if (FF.dim1!=m || FF.dim2!=n) FF.initDef(m,n);

   for (i=0; i<m; i++)
   for (j=0; j<n; j++) FF(i,j).init(F,L,S,ref,i+j*m);
}


template<class TQ, class TD>
void mxcInitQSpaceMat(const char *F, int L,
   const mxArray* C,
   wbMatrix< QSpace<TQ,TD> > &FF,
   const char ref
){
   unsigned i,j,m,n;
   const mxArray *a;

   if (!C) { FF.init(); return; }

   m=mxGetM(C); n=mxGetN(C);

   if (!mxIsCell(C) || mxGetNumberOfDimensions(C)!=2 || (m!=1 && n!=1))
   wblog(F,L,"ERR %s needs QSpace cell vector on input.",FCT);

   m*=n;

   for (i=0; i<m; i++) {
      a=mxGetCell(C,i); mxIsQSpaceVec(F,L,a);

      if (i>0) {
         if (n!=mxGetNumberOfElements(a)) { wbstring istr(str);
            sprintf(str,"cell QSpace dimensions mismatch (%d/%d)",
            n,mxGetNumberOfElements(a));
            wblog(F,L,"%s\n%s",str,istr.data);
            wblog(FL,"ERR");
         }
      }
      else {
         n=mxGetNumberOfElements(a);
         if (FF.dim1!=m || FF.dim2!=n) FF.initDef(m,n);
      }

      for (j=0; j<n; j++)
      FF(i,j).init(F,L,a,ref,j);
   }
}


template<class TQ, class TD>
mxArray* QSpaceVecVec2Mx(
   const wbvector< wbvector< QSpace<TQ,TD> > > &F
){
   mxArray *S;

   S=mxCreateCellMatrix(1,F.len);

   for (unsigned i=0; i<F.len; i++)
   mxSetCell(S,i,F[i].toMx());

   return S;
}


template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &F,
   wbindex &I
){
   if (!I.len) { QSpace<TQ,TD> X; return X.toMx(); }

   wbvector<const QSpace<TQ,TD>*> qq(I.len);
   for (unsigned i=0; i<I.len; i++) qq[i]=&(F[I[i]]);

   return qq.toMxP();
}



#endif

