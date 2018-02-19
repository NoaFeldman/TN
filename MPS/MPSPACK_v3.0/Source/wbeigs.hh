#ifndef __WB_EIGS_HH__
#define __WB_EIGS_HH__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

#define AUPDFUN  dsaupd
#define EUPDFUN  dseupd
#define AUPDSTR "DSAUPD"
#define EUPDSTR "DSEUPD"

extern "C" {


#define dsaupd dsaupd_
void dsaupd (

   const int    &,
   const char   &,
   const int    &,
   const char[2] ,
   const int    &,
   const double &,
         double[],
   const int    &,
         double[],
   const int    &,
   const int[11] ,
   const int[11] ,
         double[],
         double[],
   const int    &,

   const int    &
);



#define dseupd dseupd_
void dseupd (

   const int    &,
   const char   &,
   const int[]   ,
   const double[],
   const double[],
   const int    &,
   const double &,


   const char   &,
   const int    &,
   const char[2] ,
   const int    &,
   const double &,
         double[],
   const int    &,
         double[],
   const int    &,
   const int[11] ,
   const int[11] ,
         double[],
         double[],
   const int    &,
   const int    &
);


}"C"

char* info2str(int info) {
   switch (info) {
      case    -1: strcpy(str,"N<=0");              break;
      case    -2: strcpy(str,"NEV<=0");            break;
      case    -3: strcpy(str,"NCV<=NEV || NCV>N"); break;
      case    -4: strcpy(str,"MAXIT<0");           break;
      case    -5: strcpy(str,"invalid WHICH");     break;
      case    -6: strcpy(str,"invalid BMAT");      break;
      case    -7: strcpy(str,
                  "length of private work array WORKL not sufficient"); break;
      case    -8: strcpy(str,"error return from trid. eigenvalue "
                  "calculation (from lapack::dsteqr)");  break;
      case    -9: strcpy(str,"starting vector is zero"); break;
      case   -10: strcpy(str,"invalid MODE=IPARAM[6]");  break;
      case   -11: strcpy(str,"IPARAM[6]==1 && BMAT=='G' are incompatable"); break;
      case   -12: strcpy(str,"invalid IPARAM[0] (must be equal to 0 or 1)"); break;
      case   -13: strcpy(str,"NEV && WHICH=='BE' are incompatable."); break;
      case -9999: strcpy(str,"Could not build an Arnoldi factorization."); break;
      case   -15: strcpy(str,"HOWMNY must be one of 'A' or 'S' if RVEC==true"); break;
      case   -16: strcpy(str,"HOWMNY = 'S' not yet implemented"); break;
      case   -17: strcpy(str,"data inconsistency between DSAUPD and DSEUPD"); break;
      default: sprintf(str,"unknown info=%d",info);
   }

   return str;
}


#endif

