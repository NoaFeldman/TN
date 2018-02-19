// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

char USAGE[]="";

#include "wblib.h"

typedef struct __dummy_struct {
   int i[8];
   double d[8];
} MyStruct;


void tst_memcmp(double d1, double d2) {
    unsigned s=sizeof(double);

    wblog(FL, "memcmp(%4g, %4g) = %2d (%2d)", d1, d2,
      memcmp(&d1,&d2,s),
      memcmp(&d2,&d1,s)
    );
}


int main (int argc, char **argv) {






























































    return 0;
}

