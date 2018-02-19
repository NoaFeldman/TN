#ifndef __WB_MPFR_CC__
#define __WB_MPFR_CC__

/**********************************************************************/
   namespace Wb {
/**********************************************************************/

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

template <unsigned P>
mpfr__<P>& mpfr__<P>::init_s(const char *F, int L,
   const char *s, unsigned base
){
   if (!s || !s[0]) wblog(FL,"ERR %s() got empty string",FCT);
   int e=mpfr_set_str(f,s,base,GMP_RNDD);
   if (e) wblog(F_L,
      "ERR '%s'\n=> mpfr_set_str() returned %d @ base=%d !??",s,e,base);
   return *this;
};


template <unsigned P>
mpfr__<P>& mpfr__<P>::init_d(
   const char *F, int L, double val, char dcheck
){
   if (isnan(val) || isinf(val))
      wblog(FL,"WRN %s() got val=%g",FCT,val);

   init(F_L,val,'l');

   if (val && dcheck) {
      double x=getval_d(), e=std::fabs((x-val)/val);
      if (e>1e-14) {
         fprintf(stderr,"\n  %30.20g\n  %30.20g\n\n",val,x);
         wblog(FL,"ERR %s() got conversion error @ %.3g (1E-14)",FCT,e);
      }
      else wblog(FL,"TST %s() init_d() @ %.3g",FCT,e); 
   }
   return *this;
};

template <unsigned P>
mpfr__<P>& mpfr__<P>::init(const char *F, int L, double val,
   char lflag
){
   int e=0;


   if (::round(val)==val && std::fabs(val)<1E16) {
      mpfr_set_si(f,::round(val),DEF_RND);
      return *this;
   }

   if (std::fabs(val)>1E-4) {
      long p,q; unsigned niter; double r,v0=val;
      e=Rational(val,p,q,&r,&niter,NULL,6,1000000,1E-10,1E-12);
      if (!e) {
         if (std::fabs((r*q)/p)<1E-14) {
            mpfr_t b; mpfr_init2(b,P);
            mpfr_set_si(f,p,DEF_RND);
            mpfr_set_si(b,q,DEF_RND);
            mpfr_div(f,f,b,DEF_RND); return *this;
         }
         else { val=v0; }
      }
   }

   double q=::exp10(6-std::floor(std::log10(std::fabs(val)))), x=val*q;
   if (std::fabs(::round(x)-x)<1E-7) {
      char s[20], l=sprintf(s,"%.12g",val);
         if (l>14) wblog(FL,"WRN %s() s=%s",FCT,s);
      init_s(FL,s); return *this;
   }

   if (lflag) {
      e=mpfr_set_d(f,val,DEF_RND);
      if (e) wblog(F_L,"WRN %s() got e=%d for val=%.12g",FCT,e,val);
   }
   else { wblog(F_L,"ERR %s() "
     "only integers or simple fractions accepted (%.8g)", FCT,val);
   }

   return *this;
};


template <unsigned P>
wbstring mpfr__<P>::toStr(int base) const {
   wbstring s(prec(base)+7);
   to_str(s.data,s.len,base);
   return s;
};


template <unsigned P>
size_t mpfr__<P>::to_str(char *s, size_t n, int base) const {

   size_t l=prec(base)+7;
   mp_exp_t e;
   if (l>n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,n,l);

   mpfr_get_str(s,&e,base,0,f,GMP_RNDU);

   l=strlen(s); if (l<n) {
   l+=snprintf(s+l,n-l," @%d",int(e)); }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);

   return l;
};









template <unsigned P>
wbstring& mpfr__<P>::toStr(wbstring &sout, int base, size_t N) const {

   size_t i=0, l=3, n=prec(base)+9; 
   char xflag=0; if (int(N)<0) { xflag=1; N=0; }

   mp_exp_t e;

   if (sout.len!=n+1) { sout.init(n+1); }
   char *q, *s=sout.data; s[0]=s[1]=' ';

   q=mpfr_get_str(s+l,&e,base,N,f,GMP_RNDU);
   if (q!=s+l) wblog(FL,
      "ERR %s() mpfr_get_str returned %lX",FCT,(unsigned long)q);


   if (s[l]!='@') {
      if (e==0) {"0."
         if (s[l]=='-')
              { s[1]=s[l]; s[2]='0'; s[3]='.'; }
         else {            s[1]='0'; s[2]='.'; }
         l=1;
      }
      else if (e==-1) { e=0;"0.0"
         if (s[l]=='-')
              { s[0]=s[l]; s[1]='0'; s[2]='.'; s[3]='0'; }
         else {            s[0]='0'; s[1]='.'; s[2]='0'; }
         l=0;
      }
      else if (e==+2) { i=(--l); e=0;"."
         s[i]=s[i+1]; if (s[i]=='-') { ++i;
         s[i]=s[i+1]; }; ++i;
         if (!s[i] || !s[i+1]) wblog(FL,
            "ERR %s() got ultrashort string !??  (%s)",FCT,s+l);
         s[i]=s[i+1]; s[i+1]='.';
      }
      else { i=(--l);
         s[i]=s[i+1]; if (s[i]=='-') { ++i;
         s[i]=s[i+1]; }; s[i+1]='.';
         --e;
      }
   }
   else {
      e=0;
   }

   if (l) {
      for (i=0; i<n; ++i) {
          s[i]=s[i+l]; if (!s[i]) { break; }
      }
   }

   if (e) { l=strlen(s);
      l+=snprintf(s+l,n-l,"%c%d",base<=10 ? 'e':'@', int(e));
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);
   }

   if (xflag) {
      mpfr__<P> x; 
      try {
         x.init_s(FL,s,base);
      }
      catch (...) {
         wblog(FL,"ERR %s() `%s'\n%s",FCT,s,toStr().data); 
      }

      int i=mpfr_cmp(x.f,f);
      if (i) { wblog(FL,
         "ERR %s() got difference in base-%d representation",FCT,base);
         mpfr_out_str(stderr,base,0,  f,GMP_RNDD); fprintf(stderr,"\n");
         mpfr_out_str(stderr,base,0,x.f,GMP_RNDD); fprintf(stderr,"\n");
      }
   }

   return sout;
};


};

#endif

