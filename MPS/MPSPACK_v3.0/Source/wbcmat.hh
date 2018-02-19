#ifndef __WB_CMATRIX_HH__
#define __WB_CMATRIX_HH__

// ----------------------------------------------------------------- //
// complex matrix class
// Wb,May27,06
// ----------------------------------------------------------------- //

class wbcmat : public wbMatrix<wbcomplex> {

  public:


    void put (const char *vname="ans", const char *ws="base") const;
    void getReal(wbMatrix<double> &R) const;
    void getImag(wbMatrix<double> &I) const;
    void set(const wbMatrix<double> &R, const wbMatrix<double> &I);
    
  protected:
  private:

};


void wbcmat::put(const char *vname, const char* ws) const {

    unsigned i,j,k;
    mxArray *a=mxCreateDoubleMatrix(dim1, dim2, mxCOMPLEX);
    double *dd;

    dd=mxGetPr(a); k=0;
    for (j=0; j<dim2; j++)
    for (i=0; i<dim1; i++) dd[k++]=data[i*dim2+j].r;

    dd=mxGetPi(a); k=0;
    for (j=0; j<dim2; j++)
    for (i=0; i<dim1; i++) dd[k++]=data[i*dim2+j].i;

#ifdef MATLAB_MEX_FILE
    i=mexPutVariable(ws, vname, a);
    if (i==1){
        mexPrintf("Could not put variable %s into workspace %s.",vname,ws);
        wberror(FLINE,"");
    }
#else
    wblog(FL,"ERR %s(%s,%lX,%s) not available outsite MatLab",
    FCT,vname?vname:"",a,ws?ws:"");
#endif

    mxDestroyArray(a);
    return;
}


void wbcmat::getReal(wbMatrix<double> &R) const {
   R.init(dim1,dim2);
   for (unsigned s=dim1*dim2, i=0; i<s; i++)
   R.data[i]=data[i].r;
}

void wbcmat::getImag(wbMatrix<double> &I) const {
   I.init(dim1,dim2);
   for (unsigned s=dim1*dim2, i=0; i<s; i++)
   I.data[i]=data[i].i;
}

void wbcmat::set(const wbMatrix<double> &R, const wbMatrix<double> &I) {
   if (!R.hasSameSize(I)) wblog(FL,
   "ERR Dimension mismatch (%dx%d; %dx%d)!",R.dim1,R.dim2,I.dim1,I.dim2);

   init(R.dim1,R.dim2);

   for (unsigned s=dim1*dim2, i=0; i<s; i++)
   data[i].set(R.data[i], I.data[i]);
}


#endif

