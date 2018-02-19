#ifndef __WB_MEM_TRACK_HH__
#define __WB_MEM_TRACK_HH__

// Wb,Jul27,12

/* ------------------------------------------------------------------ */
#if ( defined __WBDEBUG__ || defined __WB_MEM_CHECK__ )
/* ------------------------------------------------------------------ */

#include "memtrack.cc"

#define WB_NEW(p,n) Wb::NEW(__FILE__,__LINE__,p,n);
#define WB_NEW_1(p,x) p = new x; Wb::gML.add_ptr(__FILE__,__LINE__,p);

#define WB_DELETE(p) \
   if (p) { Wb::gML.rm_ptr(__FILE__,__LINE__,p); delete [] p; p=NULL; }

#define WB_DELETE_1(p) \
   if (p) { Wb::gML.rm_ptr(__FILE__,__LINE__,p); delete p; p=NULL; }

#else

#define WB_NEW(p,n) Wb::NEW(p,n)
#define WB_NEW_1(p,x) { p = new x; \
  if (!p) wblog(__FILE__,__LINE__,"ERR failed to allocate instance"); \
}

#define WB_DELETE(p)   if (p) { delete [] p; p=NULL; }
#define WB_DELETE_1(p) if (p) { delete p; p=NULL; }

namespace Wb {

  template<class T>
  inline void NEW(T* &p, const size_t &n) {
     if (n) {
      { try { p = new T[n]; } catch (...) { p=NULL; } }

        if (!p) { printf("\n");
            wblog(FL,"ERR out of memory (%d * %d = %.6g GB) !!?",
            n, (int)sizeof(T), n*sizeof(T)/double(1<<30));
        }
     }
     else { p=NULL; }
  };



}

#endif

#endif

