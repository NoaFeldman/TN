#ifndef __WB_COMPLEX_HCC__
#define __WB_COMPLEX_HCC__

class wbcomplex;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
// complex number class
// Wb,Dec20,01


class wbcomplex {

  public:

    wbcomplex () : r(0.), i(0.) {};
    wbcomplex (const wbcomplex &z) : r(z.r), i(z.i) {};

    wbcomplex (double r0, double i0=0.) : r(r0), i(i0) {};

   ~wbcomplex () {};

    wbcomplex& operator= (const wbcomplex &z) {
       r=z.r; i=z.i; return *this;
    };

    wbcomplex& operator= (const double &dbl) {
       r=dbl; i=0.;  return *this;
    };

    operator double() const {
       if (i!=0) wblog(FL,
          "ERR %s() converting complex %g%+gi to double",FCT,r,i);
       return r;
    };


    void init(const double r0, const double i0) { r=r0; i=i0; };
    void set (const double r0, const double i0) { r=r0; i=i0; };

    void setRand(double fac=2., double shift=-1.);

    bool isinf() const { return (std::isinf(r) || std::isinf(i)); };
    bool isnan() const { return (std::isnan(r) || std::isnan(i)); };

    char operator== (const wbcomplex &) const;
    char operator!= (const wbcomplex &) const;

    void operator+= (const wbcomplex &);
    void operator-= (const wbcomplex &);
    void operator*= (const wbcomplex &);
    void operator/= (const wbcomplex &);

    char operator== (const double &) const;
    char operator!= (const double &) const;
    void operator*= (const double &);

    wbcomplex operator+  (const wbcomplex &) const;
    wbcomplex operator-  (const wbcomplex &) const;
    wbcomplex operator*  (const wbcomplex &) const;
    wbcomplex operator/  (const wbcomplex &) const;

    wbcomplex operator+  (const double &) const;
    wbcomplex operator-  (const double &) const;
    wbcomplex operator*  (const double &) const;
    wbcomplex operator/  (const double &) const;

    bool operator>  (const wbcomplex &) const;
    bool operator>= (const wbcomplex &) const;
    bool operator<  (const wbcomplex &) const;
    bool operator<= (const wbcomplex &) const;

    wbcomplex times   (const wbcomplex &) const;
    wbcomplex ConjProd(const wbcomplex &z)const;

    void Conj() { i=-i; };
    wbcomplex conj() const { return wbcomplex(r,-i); };


    double abs2 () const { return double(r*r+i*i); };
    double abs  () const {
       double ar=fabs(r), ai=fabs(i), x;

       if (ar>=ai) {
          if (ar==0.) return 0;
          x=i/r; return ar*std::sqrt(1.+x*x);
       }
       else {
          x=r/i; return ai*std::sqrt(1.+x*x);
       }
    };

    wbcomplex sqrt() const;
    double arg(char flag='R') const;

    wbcomplex round() const {
       return wbcomplex( ::round(r), ::round(i) );
    };

    wbstring toStr(const char *fmt="%.4g") const;
    void put (const char *vname="ans", const char *ws="base") const;

    double r,i;

  protected:

  private:
};


inline wbcomplex wbcomplex::operator+ (const wbcomplex &z) const
{  return wbcomplex(r+z.r, i+z.i); }

inline wbcomplex wbcomplex::operator- (const wbcomplex &z) const
{  return wbcomplex(r-z.r, i-z.i); }

inline wbcomplex wbcomplex::operator* (const wbcomplex &z) const
{  return wbcomplex( r*z.r-i*z.i, i*z.r+r*z.i); }

inline wbcomplex wbcomplex::operator+ (const double &dbl) const
{  return wbcomplex(r+dbl,i); }

inline wbcomplex wbcomplex::operator- (const double &dbl) const
{  return wbcomplex(r-dbl,i); }

inline wbcomplex wbcomplex::operator* (const double &dbl) const
{  return wbcomplex(r*dbl, i*dbl); }

inline wbcomplex wbcomplex::ConjProd(const wbcomplex &z) const
{  return wbcomplex(r*z.r+i*z.i, -i*z.r+r*z.i); }


inline wbcomplex wbcomplex::operator/ (const wbcomplex &b) const {
   double x, den;

   if (fabs(b.r)>=fabs(b.i)) {
      if (b.r==0) { wblog(FL,"WRN DIV/0");
      return wbcomplex(Inf,Inf); }

      x=b.i/b.r;
      den = b.r + x*b.i;

      return wbcomplex((r+x*i)/den, (i-x*r)/den);
   }
   else {
      x=b.r/b.i;
      den = b.i + x*b.r;

      return wbcomplex((x*r+i)/den, (x*i-r)/den);
   }
}

inline wbcomplex wbcomplex::operator/  (const double &dbl) const {
    if (dbl==0.) {
        wblog(FL, "ERR DIV/0! (%g+i%g)/%g", r, i, dbl);
        return wbcomplex(Inf,Inf);
    }
    return wbcomplex(r/dbl, i/dbl);
}


inline void wbcomplex::operator+= (const wbcomplex &z)
{  r+=z.r; i+=z.i; }

inline void wbcomplex::operator-= (const wbcomplex &z)
{  r-=z.r; i-=z.i; }

inline void wbcomplex::operator*= (const wbcomplex &z) {
    double rr;
    rr= r*z.r - i*z.i;
    i = r*z.i + i*z.r;
    r = rr;
}

inline void wbcomplex::operator*= (const double &dbl)
{  r*=dbl; i*=dbl; }

inline wbcomplex wbcomplex::times (const wbcomplex &z) const {
    double a,b; a=r*z.r; b=i*z.i;
    wbcomplex zz(a-b, (r+i)*(z.r+z.i)-a-b);
    return zz;
}

inline void wbcomplex::operator/= (const wbcomplex &z) {
    if (z==0.) {
        wblog(FL,"ERR Division by zero! (%g+i%g)/(%g+i%g)",r,i,z.r,z.i);
        return;
    }

    double dbl = 1. / z.abs2();
    double rr;

    rr= (r * z.r + i * z.i) * dbl;
    i = (i * z.r - r * z.i) * dbl;
    r = rr;
}


inline char wbcomplex::operator== (const wbcomplex &z) const {
    return ( float(r)==float(z.r) && float(i)==float(z.i) );
}

inline char wbcomplex::operator== (const double &dbl) const {
    return (i==0. && float(r)==float(dbl));
}

inline char wbcomplex::operator!= (const wbcomplex &z) const {
    return (!((*this)==z));
}

inline char wbcomplex::operator!= (const double &dbl) const {
    return (i!=0. || float(r)!=float(dbl));
}


inline bool wbcomplex::operator> (const wbcomplex &z) const {
   if (r>z.r || (r==z.r && i>z.i))  return 1; else return 0;
}
inline bool wbcomplex::operator>=(const wbcomplex &z) const {
   if (r>z.r || (r==z.r && i>=z.i)) return 1; else return 0;
}

inline bool wbcomplex::operator< (const wbcomplex &z) const {
   if (r<z.r || (r==z.r && i<z.i))  return 1; else return 0;
}
inline bool wbcomplex::operator<=(const wbcomplex &z) const {
   if (r<z.r || (r==z.r && i<=z.i)) return 1; else return 0;
}


inline void wbcomplex::setRand(double fac, double shift) {
   static char first_call=1;

   if (first_call) {
      srand((unsigned int)time((time_t*)NULL));
      first_call=0;
   }

   fac/=(double)RAND_MAX;

   if (shift==0.)
        { r=fac*rand();       i=fac*rand();       }
   else { r=fac*rand()+shift; i=fac*rand()+shift; }
}






inline wbcomplex wbcomplex::sqrt () const {
    if (r==0 && i==0)
    return wbcomplex(0.,0.);

    double ar=fabs(r), ai=fabs(i), w, x;

    if (ar>=ai) {
       x=ai/ar;
       w=std::sqrt(ar)*std::sqrt(0.5*(1.+std::sqrt(1+x*x)));
    }
    else {
       x=ar/ai;
       w=std::sqrt(ai)*std::sqrt(0.5*(x +std::sqrt(1+x*x)));
    }

    if (r>=0.) return wbcomplex( w, i/(2.*w) ); else
    if (i>=0.) return wbcomplex( i/(2.*w), w );
    else       return wbcomplex(-i/(2.*w),-w );
}

inline double wbcomplex::arg (const char flag) const {

    double phi=0.;

    if (r==0) {
        if (i==0)
             return  0.;
        else return (0.5 * (i>0 ? +M_PI :+M_PI));
    }

    phi=atan(i/r); if (r<0) phi+=M_PI;


    if (flag=='D') phi *= 180./M_PI; else
    if (flag!='R') wblog(FL,
    "ERR Invalid flag [D(eg), R(ad)] (%c<%d>)", flag, flag);

    return phi;
}


void wbcomplex::put(const char *vname, const char *ws) const {

    double *d;
    mxArray *a=mxCreateDoubleMatrix(1, 1, i==0. ? mxREAL : mxCOMPLEX);

    d=mxGetPr(a); d[0]=r;

    if (i!=0.) { d=mxGetPi(a); d[0]=i; }

#ifdef MATLAB_MEX_FILE
    int e=mexPutVariable(ws, vname, a);
    if (e==1) wblog(FL,
       "ERR Could not put variable %s into workspace %s.",vname,ws);
#else
    wblog(FL,"ERR %s(%s,%s,%lX) not available outsite MatLab",
    FCT,vname?vname:"",ws?ws:"",&a);
#endif

    mxDestroyArray(a);
}







inline wbcomplex operator+(const wbcomplex &z) { return z; }
inline wbcomplex operator-(const wbcomplex &z) {
   return wbcomplex(-z.r, -z.i);
}

inline wbcomplex operator*(double d, const wbcomplex &z) {
   return wbcomplex(z.r*d, z.i*d);
}

inline wbcomplex operator+(double d, const wbcomplex &z) {
   return wbcomplex(d+z.r, +z.i);
}

inline wbcomplex operator-(double d, const wbcomplex &z) {
   return wbcomplex(d-z.r, -z.i);
}

inline wbcomplex operator/(double d, const wbcomplex &z) {
   return wbcomplex(d,0)/z;
}

inline bool isnan(const wbcomplex &z) { return z.isnan(); };


template <class T> inline T real(const T &d) { return d; };
template <class T> inline T imag(const T &d) { return 0; };
template <class T> inline T conj(const T &d) { return d; };

template <class T> inline T CONJ(const T &d) { return d; };
template <class T> inline T ABS (const T &a) { return (a<0 ? -a : a); };
template <class T> inline T ABS2(const T &d) { return d*d; };
template <class T> inline T NORM (const T &a){ return (a<0 ? -a : a); };
template <class T> inline T NORM2(const T &d){ return d*d; };

inline double real(const wbcomplex &z) { return z.r; };
inline double imag(const wbcomplex &z) { return z.i; };
inline wbcomplex conj(const wbcomplex &z) {
   return wbcomplex(z.r,-z.i);
};


inline double ROUND(const double &d) { return ::round(d); }
inline wbcomplex ROUND(const wbcomplex &d) { return d.round(); }

#ifdef __WB_MPFR_HH__
template <unsigned P> inline
Wb::mpfr__<P> ROUND(const Wb::mpfr__<P> &q) { return q.round(); };
#endif


inline wbcomplex SQRT(const wbcomplex &d) { return d.sqrt(); }
inline double SQRT(const double &d) { return std::sqrt(d); }
inline long double SQRT(const long double &d) { return std::sqrt(d); }

#ifdef __WB_MPFR_HH__
template <unsigned P> inline
Wb::mpfr__<P> SQRT(const Wb::mpfr__<P> &q) { return q.sqrt(); };
#endif


#ifdef _COMPLEX_H
  #define ISREAL(x)    (typeid(x)!=typeid(wbcomplex) && typeid(x)!=typeid(complex))
  #define ISCOMPLEX(x) (typeid(x)==typeid(wbcomplex) || typeid(x)==typeid(complex))
#else                                     
  #define ISREAL(x)    (typeid(x)!=typeid(wbcomplex))
  #define ISCOMPLEX(x) (typeid(x)==typeid(wbcomplex))
#endif


template <> inline 
double ABS(const double &a) { return std::fabs(a); };
template <> inline
double NORM(const double &a) { return std::fabs(a); };

template <> inline
wbcomplex CONJ(const wbcomplex &a) { return a.conj(); };

inline double ABS  (const wbcomplex &a) { return a.abs();  };
inline double ABS2 (const wbcomplex &a) { return a.abs2(); };
inline double NORM (const wbcomplex &a) { return a.abs();  };
inline double NORM2(const wbcomplex &a) { return a.abs2(); };
 
inline bool ISINF(const double &a) { return std::isinf(a); }
inline bool ISINF(const wbcomplex &a) {
   return (std::isinf(a.r) || std::isinf(a.i)); }

inline bool ISNAN(const double &a) { return std::isnan(a); }
inline bool ISNAN(const wbcomplex &a) {
   return (std::isnan(a.r) || std::isnan(a.i)); }

inline double REAL(const double &a) { return a; }
inline double REAL(const wbcomplex &a) { return a.r; }

#ifdef __WB_MPFR_HH__

template <unsigned P> inline
Wb::mpfr__<P> ABS(const Wb::mpfr__<P> &a) { return a.abs(); };
template <unsigned P> inline
Wb::mpfr__<P> NORM(const Wb::mpfr__<P> &a) { return a.abs(); };

template <unsigned P> inline
Wb::mpfr__<P> ABS2(const Wb::mpfr__<P> &a) { return a.abs2(); };
template <unsigned P> inline
Wb::mpfr__<P> NORM2(const Wb::mpfr__<P> &a) { return a.abs2(); };

#endif


template <class T> inline double ABSDIFF_F(
   const T &d1, const T &d2
){ return double(fabs(float(d1)-float(d2))); };

template <> inline double ABSDIFF_F(
   const wbcomplex &d1, const wbcomplex &d2
){ 
   float dr=float(d1.r)-float(d2.r), di=float(d1.i)-float(d2.i);
   return std::sqrt(double(dr*dr+di*di));
};


inline wbcomplex exp(const wbcomplex &z) {
   double r=exp(z.r);
   return wbcomplex(r*cos(z.i), r*sin(z.i));
};

inline wbcomplex expi(const double &x) {
   return wbcomplex(cos(x), sin(x));
};

inline wbstring toStr(const wbcomplex &z, const char* fmt="%.4g");
inline wbstring toStr(const double &d, const char* fmt="%.4g");


#endif

