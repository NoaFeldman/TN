#ifndef __WB_CMAT_HH__
#define __WB_CMAT_HH__

/* ---------------------------------------------------------------- *
 * complex matrix build on wbMatrix<double>
 *
 * Real part (R) is always present
 * Imag part (I) is only present, when it is !=0 -- then it always
 * has the same dimensions as R.
 *
 * AWb (C) Mar 2006
 * ---------------------------------------------------------------- */

class wbCMat {

  public:

    wbCMat() : isdiag(0), nrefs(0) { memset(user,0,2); };

    wbCMat(const wbCMat &Z) { init(Z); };

    wbCMat(const mxArray *a) : isdiag(0), nrefs(0) {
       memset(user,0,2);
       init(a);
    };

    wbCMat(unsigned d1, unsigned d2) : isdiag(0), nrefs(0) {
       memset(user,0,2);
       R.init(d1,d2);
    };

    wbCMat(unsigned d1, unsigned d2, wbcomplex *z) : isdiag(0), nrefs(0) {
       unsigned i, s=d1*d2, iflag=0;

       memset(user,0,2);

       R.init(d1,d2);
       for (i=0; i<s; i++) R.data[i]=z[i].r;

       for (i=0; i<s; i++) if (z[i].i!=0.) { iflag=1; break; }

       if (!iflag) return;

       I.init(d1,d2);
       for (i=0; i<s; i++) I.data[i]=z[i].i;
    };

   ~wbCMat() {}; 

    void init(unsigned d1=0, unsigned d2=0) {
        isdiag=0; R.init(d1,d2); I.init();
    };
    void init(const mxArray *a);
    void init0(const mxArray *a);

    char isDiag() const;

    wbcomplex operator() (unsigned i, unsigned j) {
        if (i>=R.dim1 || j>=R.dim2)
        wblog(FL,"ERR Index out of bounds.");

        wbcomplex z(R(i,j),0);
        if (I.data) z.i=I(i,j);
        return z;
    }

    void mat2mxc(mxArray*, unsigned idx) const;

    void setZero() {
       memset(R.data, 0, R.dim1*R.dim2*sizeof(double));
       I.init();
    };

    void set(unsigned, unsigned, const double &);
    void set(unsigned, unsigned, const wbcomplex &);

    unsigned checkSize() const {
       if ((I.data!=NULL && (I.dim1!=R.dim1 || I.dim2!=R.dim2)) ||
           (I.data==NULL && (I.dim1!=0 || I.dim2!=0))
       ){
           wblog(FL,"WRN Severe data inconsistency! (%dx%d, %dx%d)",
           R.dim1, R.dim2, I.dim1, I.dim2);
           return 0;
       }
       return 1;
    }

    unsigned totSize() {
       unsigned n=R.dim1*R.dim2;
       if (I.data) n*=2;
       return n;
    };

    char sameSize(const wbCMat &Z) const {
       if (!checkSize() || !Z.checkSize()) 
       wblog(FL,"ERR Severe data inconsistency!");

       if (R.dim1!=Z.R.dim1 || R.dim2!=Z.R.dim2) return 0;

       return 1;
    };

    void initSize(const wbCMat &Z) {
       R.init(Z.R.dim1,Z.R.dim2);
       I.init(Z.I.dim1,Z.I.dim2);
    };

    wbCMat& init(const wbCMat &Z) {
       R=Z.R; I=Z.I; memset(user,0,2); isdiag=Z.isdiag; nrefs=0;
       return *this;
    };

    wbCMat& operator=(const wbCMat &Z) { return init(Z); };

    void operator+=(const wbCMat &Z);
    void operator*=(const wbcomplex &z);

    void conj() {
       if (!I.data) return;
       unsigned i,s=I.dim1*I.dim2;
       for (i=0; i<s; i++) I.data[i]=-I.data[i];
    };

    void transp() { R.Transpose(); if (I.data) I.Transpose(); };

    void toWbComplex(wbMatrix<wbcomplex> &Z) const;

    double maxRelDiff(const wbCMat &Z) const;

    void Contract(
       char idx, const wbCMat &M, wbCMat &Aout,
       char mflag='N') const;

    wbCMat& plus (const wbCMat &B, const wbcomplex &z=1.);
    wbCMat& minus(const wbCMat &B, const wbcomplex &z=1.);

    wbcomplex scalarProd(const wbCMat &Z) const;

    void put (const char *vname="ans", const char *ws="base") const;
    void info(const char *vname="ans") const;

    wbMatrix<double> R, I;
    char user[2];

    char isdiag;

    int nrefs;

  protected:
  private:

};


#endif

