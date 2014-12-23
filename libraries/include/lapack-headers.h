#if defined(WIN32) || defined(linux)
  #include "mkl_cblas.h"
  #include "mkl_types.h"
  #include "mkl_lapack.h"
  #include "mkl_blas.h"
#elif defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#endif

