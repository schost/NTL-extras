#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* does a forward FFT                                         */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::evaluate(Vec<zz_p>& val, const zz_pX& f) const {
  fftRep frep(INIT_SIZE, k);
  TofftRep(frep, f, k);
  long *frept = &frep.tbl[0][0];

  val.SetLength(n);
  for (long i = 0; i < n; i++)
    val[i] = frept[i];
}

/*------------------------------------------------------------*/
/* vectorial multipoint evaluation                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const{
  val.SetLength(f.length());
  for (long i = 0; i < f.length(); i++)
    evaluate(val[i], f[i]);  
}

/*------------------------------------------------------------*/
/* does an inverse FFT                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::interpolate(zz_pX& f, const Vec<zz_p>& val) const{
  fftRep frep(INIT_SIZE, k);
  long *frept = &frep.tbl[0][0];

  for (long i = 0; i < n; i++)
    frept[i] = rep(val[i]);

  FromfftRep(f, frep, 0, n-1);
}


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
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord){
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
