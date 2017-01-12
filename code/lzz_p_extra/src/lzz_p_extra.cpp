#include <NTL/lzz_p.h>

#include "lzz_p_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a = 0                                                */
/*------------------------------------------------------------*/
long order(const zz_p& a){
  if (a == 0)
    return -1;
  long o = 1;
  zz_p ap = a;
  while (ap != 1){
    ap *= a;
    o++;
  }
  return o;
}

/*------------------------------------------------------------*/
/* multiplicative order of w                                  */
/* assumes it is a power of 2                                 */
/*------------------------------------------------------------*/
long order_dyadic(const zz_p& w){
  if (w == zz_p(1))
    return 1;
  else
    return 2 * order_dyadic(w * w);
}


/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord){

  long p = zz_p::modulus();
  if ((p-1) < ord)
    LogicError("order too large with respect to p\n");

  while (1){
    a = random_zz_p();
    long ok = 1;
    zz_p ap = a;
    for (long i = 1; i < ord; i++){
      if (ap == 1)
	ok = 0;
      ap *= a;
    }
    if (ok == 1)
      return;
  }
}
