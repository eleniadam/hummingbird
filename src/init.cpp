#include "R_ext/Rdynload.h"
#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

/* Function Registration */

static const R_CallMethodDef callMethods[] = {
	CALLDEF(hummingbirdPostAdjustmentInternal, 5),
	CALLDEF(hummingbirdEMinternal, 6),
	NULL
};

void R_init_myRoutines(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}

