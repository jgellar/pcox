#include "survS.h"
#include "R_ext/Rdynload.h"

void R_init_pcox(DllInfo *info){
     R_RegisterCCallable("pcox",  "coxcount1",  (DL_FUNC) &coxcount1);
};
