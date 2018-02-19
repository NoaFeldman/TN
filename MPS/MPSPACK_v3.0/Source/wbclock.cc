#ifndef __WB_CLOCK_CC__
#define __WB_CLOCK_CC__

// -------------------------------------------------------------------- //
// -------------------------------------------------------------------- //

mxArray* WbClock::toMx(const char *F, int L) {
   return MXPut(F,L).add(name,"istr")
     .add(gettime("int"),"int")
     .add(gettime("cpu"),"cpu")
     .add(gettime("asm"),"asm")
   .toMx();
};

#endif

