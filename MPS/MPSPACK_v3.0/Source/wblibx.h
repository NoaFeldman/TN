#ifndef __WB_WBLIBX_HH__
#define __WB_WBLIBX_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

   void usage(const char *file="", int line=0, const char* estr="");
   int  isHelpIndicator(const mxArray *a);

   void wberror(const char *file, int line, const char* istr);

   void dbstop(const char* file, int line);


   wbvector<unsigned>
       Str2Idx(const char* s, unsigned offset=0, bool exitflag=1);

   template<class T>
   int Str2Idx(const char* s, wbvector<T> &idx, unsigned offset=0);

   char* memsize2Str(double x, char *s, unsigned l=-1);


template <class T>
class is_pointer_ { public: char q() { return 0; }; };

template <class T>
class is_pointer_<T*> { public: char q() { return 1; }; };

template <class T>
class is_pointer_<T**> { public: char q() { return 2; }; };

template <class T>
class is_pointer_<T***> { public: char q() { return 3; }; };

template <class T>
inline char ispointer(const T& x) { return is_pointer_<T>().q(); };


#endif

