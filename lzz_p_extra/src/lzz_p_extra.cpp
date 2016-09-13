#include <NTL/lzz_p.h>

#include "lzz_p_extra.h"

NTL_CLIENT

/*-----------------------------------------------------------*/
/* finds q that has order > n                                */
/* q > q0                                                    */
/*-----------------------------------------------------------*/
void find_root(zz_p& q, long n, const zz_p& q0){

  long p = zz_p::modulus();
  if ((p-1) < n)
    Error("order too large with respect to p\n");

  q = q0;
  while(1){
    q += 1;
    zz_p tmp = q*q;
    zz_p a = to_zz_p(1);
    long too_small = 0;

    for (long i = 0; i < n; i++){
      a *= tmp;
      if (a == to_zz_p(1)){
	too_small = 1;
	break;
      }
    }
    if (too_small == 0)
      return;
  }
} 


/*-----------------------------------------------------------*/
/* finds q that has order > n                                */
/*-----------------------------------------------------------*/
void find_root(zz_p& q, long n){
  find_root(q, n, to_zz_p(1));
}
