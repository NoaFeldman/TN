#ifndef __WB_REC_HPSORT_CC__
#define __WB_REC_HPSORT_CC__

/* ------------------------------------------------------------------ *
 * Heap sort algorithm for vector array
 * Reference: Numerical Recipies C, p337f
 * AWb (C) Jan 2006
 * ------------------------------------------------------------------ *
 * mergesort (ranks #2, vs. #3 for hpsort)
    + on average faster than hpsort
    + stable (does not reverse order within degenerate subspaces)
    - requires worst case an additional space of size N
    + in-place merge-sort: implemented in C++ STL (!) [stable_sort]
    + easily parallelizes!
    * while quicksort is faster in typical cases, it is
      (i) not stable, and (ii) has a worst-case performance of N^2
 * hpsort, for comparsion
    * is NOT stable (ie. may reverse order within degenerate subspaces)
    * slower on average, while in-place
    * still performs MEM_CPY even on already sorted array (!?)
    * reverse sort is actually faster than sorting a sorted array (!?)
 * Wb,Dec29,11
 * ------------------------------------------------------------------ */

// by including P, this ensures to keep initial order of degenerate
// subspaces, in that what is actually is sorted is [ra,(1:N)']

template<class T>
void Wb::hpsort(
    wbvector<T> &ra,
    wbperm &P,
    char dir=+1
){
    size_t i,j,l,ir, N=ra.len; PERM_T ia,*p;
    T rra;

    char q, Q=(dir>0 ? +1 : -1);

    P.Index(N); p=P.data; if (N<2) return;

    l = (N>>1);
    ir= N-1;


    for (;;) {
        if(l>0) {
            rra=ra[--l];
            ia=p[l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[0];

            ia=p[ir]; p[ir]=p[0];

            if ((--ir)==0) {
                ra[0]=rra;
                p[0]=ia;
                break;
            }
        }

        i=l;

        j=l+l+1;
        while (j<=ir) {
            if (j<ir) {
               q=NUMCMP(ra[j],ra[j+1]); if (q==0) q=NUMCMP(p[j],p[j+1]);
               if (q!=Q) j++;
            }

            q=NUMCMP(rra,ra[j]); if (q==0) q=NUMCMP(ia,p[j]);
            if (q!=Q) {
                ra[i]=ra[j];
                p[i]=p[j];

                i = j;
                j = ((j+1)<<1) - 1;
            }
            else break;
        }
        ra[i]=rra;
        p[i]=ia;
    }

    return;
};



template<class T>
int Wb::hpsort(
    T *A,
    size_t lda,
    size_t N,
    wbperm &P,
    char dir=+1,
    char lex=+1"row-major" (lex>0), else "col-major"
){
    size_t i,j,l, ia,ir, *p;

    char q, Q=(dir>0 ? +1 : -1);

    P.Index(N); p=P.data; if (N<2) return 0;

    l=(N>>1); ir=N-1;
    T x[lda];


    for (;;) {
        if (l>0) {
           ia=p[--l]; MEM_CPY<T>(x,lda,A+l*lda);
        }
        else {
            ia=p[ir];   MEM_CPY<T>(x,lda,A+ir*lda);
            p[ir]=p[0]; MEM_CPY<T>(A+ir*lda,lda,A);

            if ((--ir)==0) {
                p[0]=ia; MEM_CPY<T>(A,lda,x);
                break;
            }
        }

        i=l; j=2*l+1;

        while (j<=ir) { T *aj=A+j*lda;
            if (j<ir) {
                q=Wb::recCompare(aj,aj+lda,lda,lex);
                if (q==0) q=NUMCMP(p[j],p[j+1]);
                if (q!=Q) { ++j; aj+=lda; }
            }

            q=Wb::recCompare(x,aj,lda,lex); if (q==0) q=NUMCMP(ia,p[j]);
            if (q!=Q) {
                MEM_CPY<T>(A+i*lda,lda,A+j*lda);
                p[i]=p[j];

                i=j;
                j=((j+1)<<1)-1;
            }
            else break;
        }
        p[i]=ia; MEM_CPY<T>(A+i*lda,lda,x);
    }

    return 0;
};



template<class T>
int Wb::hpsort(
    wbMatrix<T> &ra, wbperm &P, char dir, char lex
){
    return hpsort(ra.data, ra.dim2, ra.dim1, P, dir, lex);
};


template <class T>
size_t Wb::matchSortedIdx(
   const T *A, size_t lda, size_t na,
   const T *B, size_t ldb, size_t nb, wbindex &Ia, wbindex &Ib,
   size_t m,
   char lex,"row-major" (lex>0), else "col-major"
   INDEX_T *ma_, INDEX_T *mb_,
   T eps
){
   size_t i,j,k,ia,ib,ig,iga,igb,ng,m1,m2, ma=0, mb=0, mtot=0;
   size_t *p1, *p2; int nwrn=0;
   char c;

   if (int(m)<0) { m=lda;
      if (lda!=ldb) wblog(FL,
     "WRN %s() got different record length (%d/%d)",FCT,lda,ldb);
   }
   if (!m) wblog(FL,"ERR finite m required "
      "(m=%d; %dx%d <> %dx%d)",m,lda,na,lda,nb,ldb);
   if (m>lda || m>ldb) wblog(FL,
      "ERR m out of bounds (m=%d; %d/%d)",m,lda,ldb); 
   if (!A || !B) wblog(FL,
      "ERR %s() got null-input (0x%lX, 0x%lX)",FCT,A,B); 


    wbvector<size_t> gra(2*na), grb(2*nb);

    if (lex<=0) {
       A+=(lda-m);
       B+=(ldb-m);
    }

    for (ng=iga=igb=ia=ib=0; ia<na && ib<nb;) {
        const T *a=A+ia*lda, *b=B+ib*ldb;

        c=Wb::recCompare(a,b,m,lex,eps,&nwrn);
        if (c<0) { ia++; continue; } else
        if (c>0) { ib++; continue; }

        p1=gra.data+(iga++); gra[iga++]=ia++;
        p2=grb.data+(igb++); grb[igb++]=ib++;

        m1=m2=1; ++ng;
        while (ia<na) {
            c=Wb::recCompare(a,a+lda,m,lex,eps,&nwrn);
            if (c==0) { gra[iga++]=(ia++); ++m1; a+=lda; }
            else {
               if (c>0) wblog(FL, "ERR %s() "
                  "expects ascending input order (%d; %d)",FCT,c,ia);
               break;
            }
        }
        while (ib<nb) {
            c=Wb::recCompare(b,b+ldb,m,lex,eps,&nwrn);
            if (c==0) { grb[igb++]=(ib++); ++m2; b+=ldb; }
            else {
               if (c>0) wblog(FL, "ERR %s() expects "
                  "ascending input order (%d; %d)",FCT,c,ib);
               break;
            }
        }
        mtot+=m1*m2; (*p1)=m1; (*p2)=m2;
    }

    Ia.init(mtot); Ib.init(mtot);

    for (k=ig=iga=igb=0; ig<ng; ++ig) {
        m1=gra[iga]; p1=gra.ref(iga+1); iga+=(m1+1);
        m2=grb[igb]; p2=grb.ref(igb+1); igb+=(m2+1);

        if (m1>1) ma+=m1;
        if (m2>1) mb+=m2;

        for (i=0; i<m1; ++i)
        for (j=0; j<m2; ++j,++k) {
            Ia[k]=p1[i];
            Ib[k]=p2[j];
        }
    }

    if (ma || mb) i=(ma ? ma : 1)*(mb ? mb : 1); else i=0;
    if (ma_) (*ma_)=ma;
    if (mb_) (*mb_)=mb;


    return i;
};


#endif

