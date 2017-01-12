#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does an interpolation                                      */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long j = 1; j < 10; j++){
    Vec<zz_p> q, val;
    q.SetLength(j);
    val.SetLength(j);
    for (long i = 0; i < j; i++){
      q[i] = random_zz_p();
      val[i] = random_zz_p();
    }
    zz_pX_Multipoint_General evQ(q);
    zz_pX_Multipoint * ev = &evQ;
    
    zz_pX f;
    ev->interpolate(f, val);
    for (long i = 0; i < j; i++)
      cout << val[i]-eval(f, q[i]) << " ";
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
