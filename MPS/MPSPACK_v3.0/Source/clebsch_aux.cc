#ifndef __WB_CLEBSCH_AUX_CC__
#define __WB_CLEBSCH_AUX_CC__

// ------------------------------------------------------------------ //
// See also header file clebsch.hh.
// Wb,Sep20,09 ; Wb,Dec17,14
// ------------------------------------------------------------------ //

// get sign of first non-zero element (-1 for negative, +1 otherwise)

template <class TD> inline
int CG::signFirstVal(const char *F, int L,
   const TD *d, SPIDX_T n, double eps1, double eps2
){
   SPIDX_T i; double a,x;

   for (i=0; i<n; ++i) { x=double(d[i]); a=fabs(x);
      if (a>eps1) return (x<0 ? -1 : +1);
      if (a>eps2 && F) wblog(F,L,
        "WRN %s() ignores small values (%.3g/%.3g)",FCT,x,eps2
      );
   }
   return +1;
};


template <class TD> inline
int CG::rangeSignConvention(const char *F, int L,
   TD *d, SPIDX_T n, double eps1, double eps2
){
   if (signFirstVal(F,L,d,n,eps1,eps2)<0) {
      for (SPIDX_T i=0; i<n; ++i) d[i]=-d[i];
      return -1;
   }
   return +1;
};


namespace CG {


template<class T>
double FixRational(const char *F, int L,
   T *d, SPIDX_T n,
   unsigned niter __attribute__ ((unused)),
   T eps1, T eps2 __attribute__ ((unused))
){

   double r2=0;
   if (!WbUtil<T>::isInt()) {
      double x,r; T q;


      for (SPIDX_T i=0; i<n; ++i) { q=round(d[i]);
         if (fabs(d[i]-q)<eps1) { x=d[i]-q; r2+=x*x;
            d[i]=q;
         }
      }

      if ((r=std::sqrt(r2))>CG_EPS2) {
         wblog(F_L,"WRN skipped %.2g (%d)", r,n);
      }
   }

   return r2;
};


template<>
double FixRational(const char *F, int L,
   double *d, SPIDX_T n, unsigned niter,
   double eps1, double eps2
){
   double rz,ra,r2;

   r2=Wb::FixRational(F_L,d,n,
      'r',
      niter,
      1024,
      eps1>0 ? eps1 : CG_SKIP_EPS1,
      eps2>0 ? eps2 : CG_SKIP_EPS2,
      &rz,
      &ra
   );

   if ((ra=std::sqrt(ra))>CG_EPS2) {
      wblog(F_L,"WRN skipped %.2g [%.2g, %.2g; %3.g; %d]",
      ra, std::sqrt(rz/n),std::sqrt(r2/n),CG_EPS2,n);
   }

   return r2;
};

};


template <class TD>
double cdata<TD>::SkipTiny(const char *F, int L){

   double r2=CG::FixRational(F_L, this->D.data, this->D.len, 12);
   this->Compress(F_L,CG_SKIP_EPS1);
   return r2;
};


template <class TQ, class TD>
double CData<TQ,TD>::SkipTiny(const char *F, int L){

   double e2=cgd.SkipTiny(F_L), e=std::sqrt(e2);
   if (e>CG_EPS1) wblog(F_L,"WRN %s() skipped %.3g",FCT,e);

   return e2;
};


template <class TQ, class TD>
double genRG_base<TQ,TD>::SkipTiny(const char *F, int L) {

    double e2=
        CG::FixRational(F_L,Z.data,Z.numel(),12);

    return e2;
};


template <class TQ>
SPIDX_T findMaxWeight(
   const QType &q, const wbMatrix<double> &Z, qset<TQ> *J, wbperm *P_
){
   unsigned i, r=(q.type ? q.qlen(): Z.dim2); SPIDX_T k;
   wbMatrix<double> z2(Z);
   wbperm P;

   if (Z.isEmpty()) wblog(FL,"ERR %s() "
      "got empty z-labels (%dx%d; %s)",FCT,Z.dim1,Z.dim2,q.toStr().data); 
   if (r && Z.dim2!=r) wblog(FL,
      "ERR %s()\ninconsistent z-labels for %s (%dx%d; %d)",
       FCT,q.toStr().data,Z.dim1,Z.dim2,r);

   z2.FlipCols();
   z2.sortRecs_float(P,-1);
   k=P[0];

   if (J) {
      qset<double> qm(Z.dim2,Z.rec(k));


      if (q.type==QTYPE_SUN) {

         if (r<1 || r>9) wblog(FL,
            "ERR %s() got sym='%s' !??",FCT,q.toStr().data);

         for (i=r-1; i>0; --i) {
            qm.data[i]=num2int(FL, (qm.data[i] - qm.data[i-1]) / double(i+1));
         }; qm.data[i]=num2int(FL, (qm.data[i] ));
      }
      else if (q.type==QTYPE_SpN) {





         if (r<2 || r>9) wblog(FL,
            "ERR %s() got sym='%s' !??",FCT,q.toStr().data);


         for (i=r-1; i>0; --i) {
            qm.data[i]=num2int(FL, (qm.data[i] - qm.data[i-1]) / double(i+1));
         }; qm.data[i]=num2int(FL, (qm.data[i] ));
      }

      J->initT(FL,qm);
   }

   for (i=1; i<z2.dim1; ++i) { if (z2.recDiff2(0,i)>1E-8) break; }
   if (i>1) {
      MXPut(FL,"ans").add(i+1,"i").add(z2,"z2").add(Z,"Z").add(P,"P");
      wblog(FL,"ERR %s()\nmaximum weight state not unique (%d)",FCT,i);
   }

   if (P_) P.save2(*P_);


  return k;
};


template <class TQ, class TD>
void get_SU2mat(
   TQ s2,
   wbvector<double> &sz,
   wbsparray<TD> &Sp, wbsparray<TD> &Sz,
   wbsparray<TD> &S2, wbsparray<TD> &E
){
   QType q("SU2");

   { double q=s2;
     if (q!=double(int(q)) || q<0 || q>2E3) wblog(FL,
       "ERR invalid spin S=%g (@ %.3g) !??",q/2,q-round(q)
     );
   }

   SPIDX_T i, D=SPIDX_T(s2+1);
   TD z, s=s2/2.0;
   wbsparray<TD> X1,X2,Sz2; 

   E.initIdentity(D);
   Sz.initDiag(D); Sz2.initDiag(D); Sp.initz(D,D,D-1); sz.init(D);

   for (i=0; i<D; i++)
   for (z=+s, i=0; i<D; i++, z-=1) {
       sz[i]=Sz[i]=z; Sz2[i]=z*z; if (i) {
       Sp.setRec(i-1, i-1, i, sqrt(s*(s+1)-z*(z+1))); }
   }

   wbMatProd(Sp,Sp,X1,'N','C');
   wbMatProd(Sp,Sp,X2,'C','N'); X1.Plus(FL,X2)*=0.5; X1.Plus(FL,Sz2);
   
   X1.save2(S2);
};


cgdStatus& cgdStatus::init_time(char flag) {

   if (t==CGD_DEFAULT || (!flag && t==CGD_ABELIAN)) {
      ctime=mtime=0; ID=0;
      return *this;
   }


   struct timespec now_;
   int e=clock_gettime(CLOCK_REALTIME,&now_);
   double now=now_.tv_sec + now_.tv_nsec/1E9;

   if (e<0) wblog(FL,"ERR %s() gettime returned error (e=%d) !?",FCT,e);

   if (!flag || flag=='a' || flag=='c') {
       mtime=ctime=now; setID();
   } else
   if (flag=='m') { mtime=now; 
      if (!ctime) ctime=now; else
      if (ctime>mtime) wblog(FL,"WRN %s() got ctime > mtime !?",FCT);
      if (!ID) setID();
   }
   else wblog(FL,"ERR %s() invalid flag=%c<%d>",FCT,flag,flag);

   return *this;
};


wbstring cgdStatus::to_vstr(double t, char vflag) const {

   if (t) {
      wbstring s(40);
      unsigned l=0, n=s.len; time_t tsec=t;

      struct tm *q=localtime(&tsec);
      l=snprintf(s.data,n,"%02d/%02d/%d %02d:%02d:%05.2f",
         q->tm_mon+1, q->tm_mday, q->tm_year+1900,
         q->tm_hour, q->tm_min, q->tm_sec+t-tsec
      );

      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
      return s;
   }
   return "(0)";
};


wbstring cgdStatus::toStr(char vflag) const {
   wbstring sout;

   if (!vflag) { sout=tstr(); return sout; }

   unsigned l,n;
   if (mtime<ctime) wblog(FL,"ERR %s() "
      "got invalid mtime=%.2f (ctime=%.2f)",FCT,mtime,ctime);

   if (vflag==1 || vflag=='v') {
      n=64; sout.init(n);
      l=snprintf(sout.data,n,"%.2f",ctime);
   }
   else {
      n=128; sout.init(n);
      l=snprintf(sout.data,n,"%s",to_vstr(ctime,'v').data);
   }

   if (l<sout.len) l+=snprintf(sout.data+l,n-l,
      " (%6.4g; %8x; %s)", mtime-ctime,ID,tstr());
   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return sout;
};


CGR_TYPE cgdStatus::init(const char *F, int L, const mxArray* a) {

   CGR_TYPE rt=CGR_DEFAULT;

   if (!a || mxIsEmpty(a)) {
      init(CGD_DEFAULT); return rt;
   }

   if (mxIsDouble(a)) {
      wbvector<double> cid(F_L,a); int l=3;
      if (cid.len<4 || cid.len>5) wblog(F_L,
         "ERR %s() invalid CRef::cgr data (%d)",FCT,cid.len);

      ctime  = cid[0];
      mtime  = cid[1];
      ID     = cid[2];

      if (cid.len==5) {
         rt=CGR_TYPE(cid[l++]);
         if (rt>CGR_NUM_TYPES) wblog(F_L,
            "ERR %s() rtype out of bounds (%d/%d)",FCT,rt,CGR_NUM_TYPES
         );
      }

      t=CGD_TYPE(cid[l]);
      if (t>=CGD_NUM_TYPES) wblog(F_L,
         "ERR %s() invalid CData type %g",FCT,double(cid[4])
      );
   }
   else if (mxIsStruct(a)) {
      unsigned n=mxGetNumberOfElements(a);
      int ic=mxGetFieldNumber(a,"ctime"), id=mxGetFieldNumber(a,"ID"),
          im=mxGetFieldNumber(a,"mtime"), iflag=mxGetFieldNumber(a,"flag");

      if (ic<0 || im<0 || id<0 || iflag<0 || n!=1) wblog(F_L,
         "ERR %s() unexpected cgdStatus structure (%d %d %d %d; %d)",
         FCT,ic,im,id,iflag,n
      );

      mxGetNumber(mxGetFieldByNumber(a,0,ic), ctime);
      mxGetNumber(mxGetFieldByNumber(a,0,im), mtime);
      mxGetNumber(mxGetFieldByNumber(a,0,id),  ID  );
      mxGetNumber(mxGetFieldByNumber(a,0,iflag), t );

      if (t>=CGD_NUM_TYPES) wblog(F_L,
         "ERR %s() invalid CData type %g",FCT,t
      ); 
   }
   else wblog(F_L,
     "ERR %s() unexpected cgdStatus (%s)",FCT,mxGetClassName(a)
   );

   return rt;
};


mxArray* cgdStatus::toMx() const {
   double cid[4]={ ctime, mtime, double(ID), double(t) };
   return wbvector<double>(4,cid).toMx();
};

mxArray* cgdStatus::toMx(CGR_TYPE rt) const {
   double cid[5]={ctime, mtime, double(ID), double(rt), double(t) };
   return wbvector<double>(5,cid).toMx();
};


mxArray* cgdStatus::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={"ctime","mtime","ID","flag"};
   return mxCreateStructMatrix(m,n,4,fields);
};


mxArray* cgdStatus::add2MxStruct(mxArray *S, unsigned k) const {

   mxSetFieldByNumber(S,k,0, numtoMx(ctime));
   mxSetFieldByNumber(S,k,1, numtoMx(mtime));
   mxSetFieldByNumber(S,k,2, numtoMx(ID   ));
   mxSetFieldByNumber(S,k,3,
     t==CGD_DEFAULT ? wbstring("").toMx() : toStr().toMx());

   return S;
};


template <class TQ>
QMap<TQ>& QMap<TQ>::init(const char *F, int L,
   const QVec &qv, const wbMatrix<TQ> &q1, const wbMatrix<TQ> &q2, 
   const wbMatrix< wbMatrix<TQ> > &QQ,
   const wbMatrix< wbMatrix< const wbvector< CRef<TQ> >* > > *SM
){
   unsigned i,j,k, l=0, n=0; 
   unsigned d1=q1.dim1, d2=q2.dim1, d=q1.dim2, mn=d1*d2;
   wbindex I; wbperm P; TQ *qq;

   if (!q1.dim1 || !q2.dim1) wblog(F,L,
      "ERR got empty QIDX (%d/%d)",q1.dim1, q2.dim1);
   if (q2.dim2!=d || QQ.dim1!=d1 || QQ.dim2!=d2) wblog(F,L,
      "ERR severe QDim mismatch (%d/%d, %d/%d, %d/%d)",
       q1.dim2, q2.dim2, QQ.dim1, d1, QQ.dim2, d2);
   if (SM && (SM->dim1!=QQ.dim1 || SM->dim2!=QQ.dim2)) wblog(FL,
      "ERR severe dimension mismatch in SS (%d/%d, %d/%d)",
       SM->dim1, QQ.dim1, SM->dim2, QQ.dim2
   );

   for (i=0; i<mn; ++i) { const wbMatrix<TQ> &q=QQ.data[i];
      n+=q.dim1;
      if (q.dim2!=d) wblog(F,L,
         "ERR inconsistency in qset record length (%d/%d)",q.dim2,d);
      if (SM) {
         const wbMatrix< const wbvector< CRef<TQ> >* > &S=SM->data[i];
         if (S.dim1!=q.dim1 || S.dim2!=qv.len) wblog(FL,
            "ERR severe dimension mismatch in SM (%d/%d, %d/%d)",
             S.dim1, q.dim1, S.dim2, qv.len
         ); 
      }
   }

   qvec=qv; Q1=q1; Q2=q2;
   Q.init(n,d);
   I1.init(n); I2.init(n); II.init(n,2); D.init();

   if (SM)
        cg3.init(n,qv.len);
   else cg3.init();

   for (l=i=0; i<d1; ++i)
   for (  j=0; j<d2; ++j) {
      const wbMatrix<TQ> &q=QQ(i,j); qq=q.data;
      for (k=0; k<q.dim1; ++k, ++l, qq+=d) {
         I1[l]=i; I2[l]=j;
         Q.recSetP(l,qq); if (SM) {
         cg3.recSetP(l, (*SM)(i,j).rec(k) ); }
      }
   }
 
   Q.groupRecs(P,D); I1.Select(P); I2.Select(P);
   if (SM) cg3.recPermute(P);

   for (l=i=0; i<D.len; ++i) { d=D[i];
   for (j=0; j<d; ++j,++l) {
      II(l,0)=i;
      II(l,1)=j;
   }}

   return *this;
};


template <class TQ>
void QMap<TQ>::Skip(const wbindex &Ix){

   if (!Ix.len) return;

   unsigned i;
   wbindex I,J;
   wbvector<char> m(D.len), mark(II.dim1);

   if (D.len!=Q.dim1) wblog(FL,
      "ERR %s() size inconsistencty (%d/%d)",FCT,D.len,Q.dim1); 
   if (I1.len!=II.dim1 || I2.len!=II.dim1) wblog(FL,
      "ERR %s() size inconsistencty (%d,%d/%d)",FCT,I1.len,I2.len,II.dim1); 
   if (!Ix.isUnique(D.len)) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,Ix.max()+1,D.len); 

   Ix.invert(D.len,I);
   Q.Set2Recs(I); D.Select(I);

   m.set((const wbvector<unsigned>&)Ix,1);
   mark.set(1);
   for (i=0; i<II.dim1; i++) { if (m[II(i,0)]) mark[i]=0; }
   mark.find(J);

   if (!cg3.isEmpty()) {
      if (cg3.dim1!=II.dim1) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,cg3.dim1,II.dim1); 
      cg3.Set2Recs(J);
   }

   I1.Select(J); I2.Select(J); II.Set2Recs(J);
};


template <class TQ>
template <class TD>
int QMap<TQ>::getCGZlist(const char *F, int L,
   wbMatrix<TQ> &qc,
   wbvector<TD> &dc,
   wbvector<INDEX_T>* IC,
   const wbvector<INDEX_T>* Ix,
   const qset<TQ> *q3,
   wbIndex *m3,
   const wbvector<double>* cc,
   double eps1,
   double eps2
) const {

   unsigned i,j,k,l,s,it,nt,ni,i0,iz,mQ,mq, N=0, n1=I1.len, n2=qvec.len;

   wbarray< const wbMatrix<double>* > Z(n1,n2,3);
   wbMatrix< wbMatrix<SPIDX_T> > IJ(n1,n2);
   wbMatrix< wbvector<TD> > DD(n1,n2);
   wbvector<unsigned> d0,dz;
   wbvector<char> mark(n1);
   wbindex I;
   TD a;

   CDATA_TQ X;
   wbperm p3(3);

   if (I2.len!=n1 || II.dim1!=n1 || n1!=cg3.dim1 || n2!=cg3.dim2) wblog(F,L,
      "ERR %s() severe size mismatch (%d,%d,%d/%d; %d/%d)",
       FCT, I1.len, I2.len, II.dim1, n1, n2, qvec.len
   );

   mark.set(1);
   if (Ix) {
      for (i=0; i<Ix->len; ++i) { k=Ix->data[i];
         if (k>=mark.len) wblog(FL,
            "ERR %s() index out of bounds (%d/%d)",FCT,k+1,mark.len); 
         mark[k]=0;
      }
   }

   if (q3) {
      if (q3->len!=Q.dim2 || II.dim1!=n1) wblog(FL,
         "ERR %s() size mismatch (%d/%d; %d/%d)",
          FCT, q3->len,Q.dim2, II.dim1,n1);
      for (k=0; k<Q.dim1; ++k) { if (!Q.recCompareP(k,q3->data)) break; }
      if (k==Q.dim1) {
          MXPut(FL).add(*this,"M").add(*q3,"q3");
          wblog(FL,"ERR %s()\nno match of Q=[%s] within *this QMap (%d)",
          FCT,q3->toStr().data,Q.dim1); 
      }
      for (i=0; i<II.dim1; ++i) { if (II(i,0)!=k) mark[i]=0; }

      if (m3) {
         if (mark.sum()!=1) wblog(FL,
            "ERR %s() single record expected (%d)",FCT,mark.sum());
         i=mark.find1(); if (int(i)<0) wblog(FL,"ERR %s() !??",FCT);

         wbvector<unsigned> S(n2);
         for (j=0; j<S.len; ++j) {
            S[j]=cg3(i,j)->data[0].cgr->getOM(FL);
         }

         if (m3->isEmpty()) { m3->init(S); }
         else {
            if (m3->SIZE!=S) wblog(FL,
               "ERR %s() multiplicity setting changed !?? [%s] <> [%s]",
               FCT,m3->SIZE.toStr().data, S.toStr().data);
         }
         if (!(++(*m3))) return 0;
      }
   }
   else if (m3)
   wblog(FL,"ERR %s() m3 may only be used together with q3",FCT); 

   if (cc) {
      if (!m3) wblog(FL,
         "ERR %s() cc may only be used together with q3 and m3",FCT);
      if (m3->SIZE.numelFindGT(1)!=1) wblog(FL,
         "ERR %s() cc only applies with outer multiplicity\n"
         "in one symmetry only",FCT);
   }

   if (!mark.any()) {
      wblog(FL,"WRN %s() got empty list into CGR",FCT);
      qc.init(); dc.init(); if (IC) IC->init(N);
      return 0;
   }

   qvec.Qlen(d0,dz);

   for (i=0; i<n1; ++i) { if (!mark[i]) continue; ni=1;
   for (j=0; j<n2; ++j) { const wbvector< CRef<TQ> > &c3=(*cg3(i,j));

       if (!c3.len) wblog(FL,"ERR %s() got empty cg3's !?",FCT);
       if (!c3[0].cgr) wblog(FL,"ERR %s() got empty cgr !?",FCT);

       const CDATA_TQ &S = (*c3[0].cgr);

       c3[0].getP(p3,'i');
       if (p3.len!=3) wblog(FL,
          "ERR %s() got invalid cgp.len=%d",FCT,p3.len);

       if (S.qs.len!=3*d0[j]) wblog(FL,"ERR CData "
          "length inconsistency in Q (%d != 3*%d)",S.qs.len,d0[j]); 

       if (dz[j]) { for (k=0; k<3; ++k) {
          Z(i,j,k) = &gRG.getR(FL,S,p3[k]).Z;

          if (Z(i,j,k)->dim2!=dz[j]) wblog(FL,
             "ERR CData length inconsistency (%s) Z(%d,%d,%d): "
             "%s @ dz=%d", S.qStr().data, i+1,j+1,k+1,
             Z(i,j,k)->sizeStr().data, dz[j]
          );
       }}

       if (!m3 || m3->SIZE[j]==1) {
          if (S.cgd.D.len) {
             ni*=S.cgd.nnz();
             S.cgd.IDX.getCols(p3,IJ(i,j));
             DD(i,j)=S.cgd.D;
          }
          else if (S.isAbelian()) { TD x=1;
             IJ(i,j).init(1,3); DD(i,j).init(1,&x);
          }
          else wblog(FL,
             "ERR %s() got invalid CData\n%s",FCT,S.toStr().data);

          if (c3.len!=1 || c3[0].cgw.len!=1) wblog(FL,"ERR %s() "
             "got invalid CRef (%d; %d)",FCT,c3.len,c3[0].cgw.len);
          DD(i,j)*=c3[0].cgw[0];
       }
       else if (cc) {
          cgsparray x; QSet<TQ> Q0,Q2;
          if (cc->len != m3->SIZE[j] || c3.len!=cc->len) wblog(FL,
             "ERR %s() length mismatch of coefficients having OM\n"
             "(%d/%d/%d)",FCT,cc->len,m3->SIZE[j],c3.len); 
          for (l=0; l<c3.len; ++l) {
             X.init(c3[l]); 
             X.cgd*=(cc->data[l]);
             if (l) {
                x+=X.cgd;
                if (Q2.init(c3[l])!=Q0) wblog(FL,
                   "ERR %s() got QSet mismatch in OM !?",FCT
                );
             }
             else {
                x=X.cgd; Q0.init(c3[l]);
             }
          }
          x.SkipTiny(1E-12);

          ni*=x.nnz();
          x.IDX.getCols(0,2,IJ(i,j));
          DD(i,j)=x.D;
       }
       else {
          X.init(c3.el((*m3)[j]));
          const cgsparray &x=X.cgd;
          ni*=x.nnz();
          x.IDX.getCols(0,2,IJ(i,j));
          DD(i,j)=x.D;

       }

       if (eps1<=0) {
       k=S.cgd.D.numZeros(1E-12); if (k) {
          wblog(FL,"WRN cg3(%dx%d): (%d,%d) got zero value%s (%d/%d) !??",
          n1, n2, i+1,j+1, k!=1 ? "s":"",k, S.cgd.D.len); 
       }}

   }; N+=ni; if (ni==0) wblog(FL,"ERR got ni=%d !?",ni); }

   mQ=d0.sum()+dz.sum(); l=n2-1;
   qc.init(N,3*mQ); dc.init(N); if (IC) IC->init(N);

   for (k=i=0; i<n1; ++i) { if (!mark[i]) continue;
      unsigned sz[n2], I[n2]; nt=1;
      for (j=0; j<n2; ++j) { I[j]=0; sz[j]=DD(i,j).len; nt*=sz[j]; }

      for (it=0; it<nt; ++it, ++k) {
         TQ *q=qc.rec(k);
         TD &c=dc[k]; c=1; if (IC) IC->data[k]=i;

         for (i0=iz=0, j=0; j<n2; ++j) {
            if (!cg3(i,j) || !cg3(i,j)->len || !cg3(i,j)->data[0].cgr)
               wblog(FL,"ERR %s() unexpected cg3 data",FCT);

            const QSet<TQ> &Q=(cg3(i,j)->data[0]);

            wbMatrix<SPIDX_T> &I1=IJ(i,j);
            mq=Q.t.qlen();

            c*=DD(i,j)[I[j]];

            if (d0[j]) { s=d0[j]*sizeof(TQ);
               memcpy(q     , Q.qs.data     , s);
               memcpy(q+  mQ, Q.qs.data+  mq, s);
               memcpy(q+2*mQ, Q.qs.data+2*mq, s);
               i0+=d0[j]; q+=d0[j];
            }
            if (dz[j]) {
               cpyRange(Z(i,j,0)->ref(I1(I[j],0)), q,      dz[j]);
               cpyRange(Z(i,j,1)->ref(I1(I[j],1)), q+mQ,   dz[j]);
               cpyRange(Z(i,j,2)->ref(I1(I[j],2)), q+2*mQ, dz[j]);
               iz+=dz[j]; q+=dz[j];
            }
         }

         a=fabs(c); if (a<eps1) { --k;
            if (a>eps2) wblog(FL,"WRN %s() got small CGC data "
              "(%.3g; %g %g)",FCT,c,eps1,eps2);
         }

         j=0; I[0]++;
         while(I[j]>=sz[j] && j<l) { I[j]=0; ++I[++j]; }
      }
   }

   if (k!=N) {
      if (!k || k>N) wblog(FL,
         "ERR %s() got no CG data (%d/%d) !??",FCT,k,N);
      qc.dim1=k; dc.len=k; if (IC) IC->len=k;
   }

   return 1;
};


template <class TQ>
template <class TD>
void QMap<TQ>::getIdentityQ(const char *F, int L,
   const wbvector<INDEX_T> &Sa, const wbvector<INDEX_T> &Sb,
   QSpace<TQ,TD> &A, char vflag
) const {

   unsigned i,j,l,k,i0,i1,im,mi,N, db,Db,
        n=I1.len, QDIM=qvec.Qlen(), nq=qvec.len;
   char isa=qvec.isAbelian();

   wbMatrix<INDEX_T> mm(n,nq);
   wbvector<unsigned> d0,dz;
   wbvector< wbvector<INDEX_T> > DD(D.len), MM(D.len);
   double sd=0;
   RTD cfac=1;

   qvec.Qlen(d0,dz);

   if (I2.len!=n || II.dim1!=n || cg3.dim1!=n || cg3.dim2!=nq)
       wblog(F,L,"ERR %s() severe size mismatch (%d,%d,%d/%d; %d/%d)",
       FCT, I1.len, I2.len, II.dim1, cg3.dim1, cg3.dim2, qvec.len);
   if (Q1.dim1!=Sa.len || Q1.dim2!=QDIM || Q.dim1!=D.len ||
       Q2.dim1!=Sb.len || Q2.dim2!=QDIM || Q.dim2!=QDIM) wblog(FL,
      "ERR QMap::%s() block size mismatch\n%dx%d %dx%d; "
      "%dx%d %dx%d; %dx%d %dx%d", Q1.dim1,Q1.dim2, Sa.len,QDIM,
       Q2.dim1,Q2.dim2, Sb.len,QDIM, Q.dim1, Q.dim2, D.len, QDIM
   );

   for (i=0; i<D.len; i++) { DD[i].init(D[i]); MM[i].init(D[i]); }
   for (i=0; i<n; i++) {
      INDEX_T &di=DD.el(II(i,0)).el(II(i,1)), si=Sa.el(I1[i])*Sb.el(I2[i]);
      if (di) wblog(FL,"ERR %s() block size not unique (%d/%d)",FCT,di,si);
      di=si;
   }

   for (N=i=0; i<n; i++, N+=mi) {
      for (mi=1, j=0; j<nq; j++) {
         const wbvector< CRef<TQ> > &c3=(*cg3(i,j));
         if (!c3.len || !c3[0].cgr) wblog(FL,
            "ERR %s() got invalid cg3 data",FCT);
         if (c3[0].getOM(FL)!=c3.len) wblog(FL,
            "ERR %s() got OM mismatch (%d/%d)",FCT,c3.len,c3[0].getOM(FL));
         mi*=(mm(i,j)=c3.len);
      
         if (c3[0].cgr->qs.len!=3*d0[j]) wblog(FL,
            "ERR CData length inconsistency in Q (%d != 3*%d)",
            c3[0].cgr->qs.len, d0[j]
         );
      }

      MM[II(i,0)][II(i,1)]=mi;

   }

   A.init(N,3,QDIM);
   if (!isa) { A.qtype=qvec; A.setupCGR(); }

   for (l=i=0; i<n; i++) {
      wbIndex I(nq,mm.rec(i)); im=0; i0=II(i,0); i1=II(i,1);

      Db=prodRange(DD[i0].data, MM[i0].data, DD[i0].len);
      k=prodRange(DD[i0].data, MM[i0].data, i1, INDEX_T(0));
      db=DD[i0][i1];

      while (++I) {
         A.QIDX.recSetB(l,0,QDIM,Q1.rec(I1[i]));
         A.QIDX.recSetB(l,1,QDIM,Q2.rec(I2[i]));
         A.QIDX.recSetB(l,2,QDIM,Q .rec(i0));

         if (!isa) {
            for (cfac=1, j=0; j<nq; ++j) {
               if (qvec[j].isAbelian()) {
                  A.CGR(l,j).initAbelian();
               }
               else {
                  A.CGR(l,j)=cg3(i,j)->data[I[j]];
                  cfac*=A.CGR(l,j).NormSignR();
               }
            }
         }
         
         A.DATA[l]->initIdentityB3(
            Sa[I1[i]], Sb[I2[i]], Db, k+im*db, TD(cfac));

         sd+=A.DATA[l]->numel();
         if (sd>2E8 && vflag=='V') wblog(FL,
             "WRN %s() occupied memory (%d/%d): %.3g GB",
              FCT,i+1,n,sd/(1<<27)
         );

         ++l; ++im;
      }
   }

   A.otype=QS_AMATRIX;
   A.itags.init3();

   { wbvector<INDEX_T> s, S;
     A.getDim(s,&S); if (S[2]!=S[0]*S[1]) {
        MXPut(FL).add(*this,"M").add(A,"A").add(s,"sd").add(S,"S");
        wblog(FL,"ERR %s() overall dimension mismatch\n[ %s ]: [ %s ]",
        FCT,s.toStr().data,S.toStr().data);
     }
   }

   if (sd>5E8 && vflag) wblog(FL,
   "WRN %s() using %.3g GB memory (%d blocks)",FCT,sd/(1<<27),n);
};


template <class TQ>
mxArray* QMap<TQ>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={
      "qvec","Q1","Q2","Q", "D","I1","I2","II","CG"
   };
   return mxCreateStructMatrix(m,n,9,fields);
}

template <class TQ>
void QMap<TQ>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      unsigned s=0; int k=0;
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (k=mxGetFieldNumber(S,"Q2"))<0) wblog(FL,
      "ERR %s must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,k);
   }

   mxSetFieldByNumber(S,i, 0, qvec.toStr().toMx());
   mxSetFieldByNumber(S,i, 1, Q1    .toMx());
   mxSetFieldByNumber(S,i, 2, Q2    .toMx());
   mxSetFieldByNumber(S,i, 3, Q     .toMx());
   mxSetFieldByNumber(S,i, 4, D     .toMx());
   mxSetFieldByNumber(S,i, 5,(I1+1) .toMx());
   mxSetFieldByNumber(S,i, 6,(I2+1) .toMx());
   mxSetFieldByNumber(S,i, 7,(II+1) .toMx());
   mxSetFieldByNumber(S,i, 8, cg3   .toMxP('c'));
};


template <class TQ>
QMap<TQ>& QMap<TQ>::init(const char *F, int L, const mxArray *S) {

   mxArray *a; if (!S || !mxIsStruct(S)) { init(); return *this; }

   a=mxGetField(S,0,"qvec"); if (a) qvec.init(F,L,a);
   a=mxGetField(S,0,"Q1");   if (a) Q1.init(F,L,a);
   a=mxGetField(S,0,"Q2");   if (a) Q2.init(F,L,a);
   a=mxGetField(S,0,"Q" );   if (a) Q .init(F,L,a);
   a=mxGetField(S,0,"D" );   if (a) D .init(F,L,a);
   a=mxGetField(S,0,"I1");   if (a) I1.init(F,L,a);
   a=mxGetField(S,0,"I2");   if (a) I2.init(F,L,a);
   a=mxGetField(S,0,"II");   if (a) II.init(F,L,a);

   return *this;
}


template <class T>
mxArray* pairPatternCG<T>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"i","j","k","fac"};
   return mxCreateStructMatrix(m,n,4,fields);
};

template <class T>
void pairPatternCG<T>::add2MxStruct(mxArray *S, unsigned l) const {

   mxSetFieldByNumber(S,l, 0, numtoMx(i+1));
   mxSetFieldByNumber(S,l, 1, numtoMx(j+1));
   mxSetFieldByNumber(S,l, 2, (k+1).toMx());
   mxSetFieldByNumber(S,l, 3, fac.toMx());
};



#endif

