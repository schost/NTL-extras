#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "mat_ZZ_pX_extra.h"
#include "magma_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(int opt){
  ZZ p = to_ZZ("234898888888888888888888888888888888888888881111111111");
  ZZ_p::init(p);

  Mat<ZZ_pX> A, B, C, C2;
  random_mat_ZZ_pX(A, 10, 10, 200);
  random_mat_ZZ_pX(B, 10, 200, 200);

  mul_CRT_CTFT(C, A, B);
  mul_direct(C2, A, B);

  for (long u = 0; u < C.NumRows(); u++)
    for (long v = 0; v < C.NumCols(); v++)
      cout << C[u][v]-C2[u][v];
  cout << endl;
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
