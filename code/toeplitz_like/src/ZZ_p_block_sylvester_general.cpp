#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* is this useful?                                    */
/*----------------------------------------------------*/
void ZZ_p_block_sylvester_general::init(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec){
  this->type = type;
  this->prec = prec;
  f = fs;
  initialized = true;
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_block_sylvester_general::mul_right(const Vec<ZZ_p> &rhs){
  if (!initialized) 
    throw "must call init first";
  Vec<ZZ_pX> rhs_poly;
  create_lhs_list(rhs_poly,rhs);
  ZZ_pX result;
  for (long i = 0; i < f.length(); i++)
    result += trunc(f[i] * rhs_poly[i], prec);
  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(result, i);
  return v;
}

/*----------------------------------------------------*/
/* input: Vec of polynomials fs                       */
/*        type                                        */
/*        output precision                            */
/*----------------------------------------------------*/
ZZ_p_block_sylvester_general::ZZ_p_block_sylvester_general(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec):
  ZZ_p_block_sylvester(type, prec) {
  init(fs, type, prec);
}
