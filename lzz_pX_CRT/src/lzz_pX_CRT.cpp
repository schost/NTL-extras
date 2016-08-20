#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CRT over zz_p                                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* multiple remainders                                        */
/*------------------------------------------------------------*/
void zz_pX_CRT::multimod(Vec<zz_pX>& val, const zz_pX& f) const {
  long n = this->n;
  val.SetLength(n);
  long d = this->tree.length();
  val.SetLength(n);
  val[0] = f % tree[d-1][0];

  for (long i = d-2; i >= 0; i--){
    long len = tree[i].length();
    if (len & 1)
      val[len-1] = val[(len-1)/2];
    for (long j = len/2-1; j >= 0; j--){
      val[2*j+1] = val[j] % tree[i][2*j+1];
      val[2*j] = val[j] % tree[i][2*j];
    }
  }
}

/*------------------------------------------------------------*/
/* Chinese remainder                                          */
/*------------------------------------------------------------*/
void zz_pX_CRT::combine(zz_pX& f, const Vec<zz_pX>& val) const {
  long n = this->n;
  long d = this->tree.length();
 
  Vec<zz_pX> combinations;
  combinations.SetLength(n);

  for (long i = 0; i < n; i++)
    combinations[i] = MulMod((val[i] % tree[0][i]), cofactors[i], tree[0][i]);
  for (long i = 0; i < d-1; i++){
    long len = tree[i].length();
    for (long j = 0; j < len/2; j++)
      combinations[j] = combinations[2*j]*tree[i][2*j+1] + combinations[2*j+1]*tree[i][2*j];
    if (len & 1)
      combinations[(len-1)/2] = combinations[len-1];
  }

  f = combinations[0];
}
