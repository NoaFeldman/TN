#ifndef __WB_REC_HPSORT_HH__
#define __WB_REC_HPSORT_HH__

namespace Wb {

/* ------------------------------------------------------------------ */
// generalized matrix/record sort
// sorted by rows/records of length m, with actual leading length lda
/* ------------------------------------------------------------------ */

template<class T>
int hpsort(
    T *A,
    size_t lda,
    size_t N,
    wbperm &P,
    char dir=+1,
    char lex=+1"row-major" (lex>0), else "col-major"
);

template<class T>
int hpsort(
    wbMatrix<T> &ra,
    wbperm &P,
    char dir=+1,
    char lex=+1"row-major" (lex>0), else "col-major"
);

template<class T>
void hpsort(
    wbvector<T> &ra,
    wbperm &P,
    char dir=+1
);

template<class T>
void hpsort(wbvector<T> &ra) { wbperm P; hpsort(ra,P); }


template <class T>
size_t matchSortedIdx(
    const T *da, size_t lda, size_t na,
    const T *db, size_t ldb, size_t nb, wbindex &Ia, wbindex &Ib,
    size_t m,
    char lex,"row-major" (lex>0), else "col-major"
    INDEX_T *ma=NULL,
    INDEX_T *mb=NULL,
    T eps=0
);

};

#endif

