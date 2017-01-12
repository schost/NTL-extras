#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long j = 1; j < 10; j++){
    Vec<zz_p> q;
    q.SetLength(j);
    for (long i = 0; i < j; i++){
      q[i] = random_zz_p();
    }
    zz_pX_Multipoint * ev;
    if (1){
      zz_pX_Multipoint_General evQ(q);
      ev = &evQ;
      cout << ev->length() << endl;
    }
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
