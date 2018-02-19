#ifndef __WB_IDAT_HCC__
#define __WB_IDAT_HCC__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* class to manage iteration data of BiCG and the likes */

class IDAT {

  public:

    IDAT () { memset(this, 0, sizeof(IDAT)); };
    IDAT (unsigned maxit) { init(maxit); };

    void init(unsigned maxit=0) {
       resvec.init(); memset(this, 0, sizeof(IDAT));
       resvec.init(maxit);
    }

   ~IDAT () {}; 

    void put(const char *vname="ans", const char *ws="base") const;

    char flag;
    unsigned niter;
    double relres;
    wbvector<double> resvec;

  protected:
  private:

};


void IDAT::put(const char *vname, const char *ws) const {

    const char *fnames[] = { "flag", "niter", "relres", "resvec" };
    mxArray *S=mxCreateStructMatrix(1,1, 4,fnames);

    mxSetFieldToNumber(S,0,0, (double)flag);
    mxSetFieldToNumber(S,0,1, (double)niter);
    mxSetFieldToNumber(S,0,2, relres);

    resvec.mat2mxs(S,3,1);

#ifdef MATLAB_MEX_FILE
    unsigned i=mexPutVariable(ws, vname, S);
    if (i==1) wblog(FL,
       "ERR Could not put variable %s into workspace %s.",vname,ws);
#else
    wblog(FL,"ERR %s(%s,%s,%lX) not available outsite MatLab",
    FCT,vname?vname:"",ws?ws:"",&S);
#endif

    mxDestroyArray(S);
}


#endif

