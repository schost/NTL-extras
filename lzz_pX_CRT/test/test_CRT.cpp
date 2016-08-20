#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a CRT + multimod                                      */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long j = 1; j < 20; j++){
    Vec<zz_pX> val, q, mmod;
    q.SetLength(j);
    val.SetLength(j);
    for (long i = 0; i < j; i++){
      q[i] = random_zz_pX(i*3+5);
      val[i] = random_zz_pX(i*3+5) % q[i];
    }
    zz_pX_CRT crt = zz_pX_CRT(q);
    zz_pX f;
    crt.combine(f, val);
    crt.multimod(mmod, f);
    for (long i = 0; i < j; i++){
      cout << val[i]-(f % q[i]) << ",";
      cout << mmod[i] - (f % q[i]) << "  ";
    }
    cout << endl;
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
