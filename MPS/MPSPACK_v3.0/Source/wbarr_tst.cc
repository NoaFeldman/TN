
// reference counting for entier wbarray object including SIZE
// -> cannot create ShapeAliases (!) without duplicating the whole thing
// -> use reference counting within data[] only.
// Wb,Sep16,09

/* *** SCRATCH / DEPRECATED *** */

class wbarray_header {
    unsigned nrefs, ref, rank, NH, numel;
    double fac, add; char flag;

    wbarray_header() 
     : nrefs(0), rank(0), NH(0), numel(0), fac(1.), add(0.), flag=0 {};

    void init(unsigned r, unsigned nh, unsigned ne) {
       nrefs=1;
       rank=r; NH=nh; numel=ne;
       fac=1.; add=0.; flag='N';
    };
};

template<class T>
class wbarray {

  public:

    wbarray(unsigned d1=0)
     : data(NULL) { unsigned s[]={d1}; INIT(s,1); };

    wbarray(const UVEC &d)
     : data(NULL) { INIT(d.data, d.len); };

    wbarray(const wbarray<T> A) {
       data=A.data; nrefs()++;
    }

   ~wbarray() {
       DELETE_DATA();
    };

    void instantiate() {
       wbarray_header &h=header();
       if (h.nrefs<=1) return;

       h.nrefs--;

       unsigned nt=h.NH+h.numel;
       T *d0=data-d.NH;

       data = new T[nt];
       memcpy(data,d0,nt*sizeof(T));
       data+=h.NH;
    }

    void init(unsigned d1=0) {
       unsigned s[]={ (d1==0) ? 0 : 1, d1 };
       RENEW(wbvector<unsigned>(2,s));
    };

    template<class T0>
    void init(const wbarray<T0> &a) { 
       RENEW(a.SIZE);
       for (unsigned n=SIZE.prod(), i=0; i<n; i++)
       data[i]=T(a.data[i]);
    };

    T& operator[] (unsigned i) const { return data[i]; }

    T& operator() (const wbvector<unsigned> &I) const {
       if (I.len!=rank()) wblog(FL,
       "ERR %s - invalid index set (%s; %d)",__FUNCTION__,
       I.toStr().data, rank());

       return operator()(I.data);
    }

    T& operator() (const unsigned *I) const {
       const unsigned *SIZE=size(), r=rank(), l=r-1;
       unsigned i=0,j=0;

       for (j=I[l], i=l-1; i<r; i--) j+=(I[i]+j*SIZE[i]);
       return data[j];
    }

  protected:

    wbarray_header& header() { return *(((wbarray_header*)data)-1); }

    unsigned* size () { return (unsigned*)(data - header().NH); }
    unsigned& nrefs() { return header().nrefs; }
    unsigned& rank () { return header().rank;  }

    void DELETE_DATA() {
       if ((--header().nrefs)>0) return;

       data -= header().NH;

#ifndef __WB_MEM_CHECK__
       delete [] data;
#else
       Wb::DELETE(FL,data);
#endif
       data=NULL;
    }

    void INIT(const wbvector<unsigned> &S) { INIT(S.data,S.len); };
    void INIT(const unsigned *S, unsigned r) {
        DELETE_DATA(); if (r==0) return;

        unsigned i,n,NH,s, sT=sizeof(T), *Sh,
           nh=sizeof(wbarray_header),
           ns=r*sizeof(unsigned);

        for (n=S[0]; i=0; i<r; i++) { n*=S[i]; }
        if (n==0) return;

        NH=(ns+nh)/sT+1; s=NH+n;
#ifndef __WB_MEM_CHECK__
        data = new T[s]; if (!data) MemErrMsg(FL,s);
#else
        Wb::NEW(FL,data,s);
#endif
        Sh=(unsigned*)data; memcpy(Sh,S,ns);

        data+=NH;

        header().init(r,NH,n);
    }




    T* data;

    char isref;

  private:

}


#endif

