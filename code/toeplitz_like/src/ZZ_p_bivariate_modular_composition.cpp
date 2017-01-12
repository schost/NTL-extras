#include <NTL/tools.h>
#include <time.h>
#include <cmath>

#include "ZZ_p_block_sylvester.h"
#include "ZZ_pX_extra.h"
#include "mat_ZZ_pX_extra.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* is this useful?                                    */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::init (const ZZ_pX &f, const Vec<long> &type_new, long prec_new){
  initialized = true;
  type = type_new;
  f_field = f;
  
  prec = prec_new;
  sqrtP = ceil(sqrt(type.length()));
  ZZ_pX running = ZZ_pX{1};
  Vec<ZZ_pX> fs;
  
  fs.SetLength(sqrtP);
  ZZ_pX_poly_multiplier multF(f, prec);
  for (long i = 0; i < sqrtP; i++){
    fs[i] = running;
    multF.mul(running, running);
    trunc(running, running, prec);
  }
  
  F_field = running;
  create_rhs_matrix(B, fs);
}

/*----------------------------------------------------*/
/* given the result of create_lhs_list, returns       */
/* a matrix representation with dimension             */
/* ceil(len(type)/sqrtP) x sqrtP                      */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::create_lhs_matrix(Mat<ZZ_pX> &A, const Vec<ZZ_pX> &v){
  A.SetDims(ceil((type.length()*1.0) / sqrtP), sqrtP);
  for (long i = 0; i < A.NumRows(); i++)
    for (long j = 0; j < A.NumCols(); j++)
      if(A.NumCols() * i + j < v.length())
	A.put(i, j, v[A.NumCols() * i + j]);
}

/*----------------------------------------------------*/
/* creates a matrix of f^i for i = 0..sqrtP-1         */
/* partitions each f^i into len(type) slices          */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::create_rhs_matrix(Mat<ZZ_pX> &B,const Vec<ZZ_pX> &v){

  long cols = ceil( (1.0*prec) / (1.0*(max_of_type+1)));
  
  B.SetDims(v.length(), cols);
  for (long i = 0; i < B.NumRows(); i++){
    long acc = 0;
    for (long j = 0; j < B.NumCols(); j++){
      ZZ_pX p;
      for (long s = 0; s < max_of_type + 1; s++)
	SetCoeff(p, s, coeff(v[i], acc + s));
      acc += max_of_type + 1;
      B.put(i, j, p);
    }
  }
 }

/*----------------------------------------------------*/
/* converts the sliced up matrix into a vector        */
/*----------------------------------------------------*/
void ZZ_p_bivariate_modular_composition::deslice(Vec<ZZ_pX> &D, const Mat<ZZ_pX> &C){

  D.SetLength(C.NumRows());

  for (long i = 0; i < C.NumRows(); i++){
    ZZ_pX v;

    v.rep.SetLength(prec + max_of_type + 1);
    ZZ_p* cv = v.rep.elts();
    const ZZ_pX * Ci = C[i].elts();

    const ZZ_p * cC = Ci[0].rep.elts();
    long old_len = Ci[0].rep.length();
    for (long j = 0; j < old_len; j++)
      cv[j] = cC[j];
    cv += max_of_type + 1;

    for (long a = 1; a < C.NumCols(); a++){
      const ZZ_p * cC = Ci[a].rep.elts();
      long new_len = Ci[a].rep.length();
      for (long j = 0; j < min(old_len-(max_of_type+1), new_len); j++)
	cv[j] += cC[j];
      for (long j = max(0, old_len-(max_of_type+1)); j < new_len; j++)
	cv[j] = cC[j];
      old_len = new_len;
      cv += max_of_type + 1;
    }
    v.normalize();
    D[i] = trunc(v, prec);
  }

}

/*----------------------------------------------------*/
/* multiplies rhs using Horner's rule                 */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right_Horners(const Vec<ZZ_p> &rhs){
  if (!initialized) 
    throw "must init first";
  Vec<ZZ_pX> rhs_poly;
  create_lhs_list(rhs_poly, rhs);
  ZZ_pX result = rhs_poly[rhs_poly.length()-1];
  for (long i = rhs_poly.length()-2; i >= 0; i--)
    result = trunc((result * f_field), prec) + rhs_poly[i];

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(result, i);
  return v;
}

/*----------------------------------------------------*/
/* mult using the baby steps / giant steps algorithm  */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right_comp(const Vec<ZZ_p> &rhs){

  if (!initialized) 
    throw "must init first";
  Mat<ZZ_pX> A0, A;
  Vec<ZZ_pX> rhs_poly;

  create_lhs_list(rhs_poly, rhs);
  create_lhs_matrix(A0, rhs_poly);
  mul_CRT_CTFT(A, A0, B);
  Vec<ZZ_pX> B1; // TODO: change this name!

  deslice(B1, A);

  ZZ_pX p;
  p = B1[B1.length() - 1];
  ZZ_pX_poly_multiplier multF(F_field, prec);
  for (long i = B1.length()-2; i >= 0; i--){
    multF.mul(p, p);
    p = trunc(p, prec) + B1[i];
  }

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(p, i);

  return v;
}

/*----------------------------------------------------*/
/* header for multiplying                             */
/* TODO: thresholds                                   */
/*----------------------------------------------------*/
Vec<ZZ_p> ZZ_p_bivariate_modular_composition::mul_right(const Vec<ZZ_p> &rhs){
  return mul_right_comp(rhs);
}

/*----------------------------------------------------*/
/* default constructor; does nothing                  */
/*----------------------------------------------------*/
ZZ_p_bivariate_modular_composition::ZZ_p_bivariate_modular_composition(){
}

/*----------------------------------------------------*/
/* input: Vec of polynomials fs                       */
/*        type                                        */
/*        output precision                            */
/*----------------------------------------------------*/
ZZ_p_bivariate_modular_composition::ZZ_p_bivariate_modular_composition(const ZZ_pX& f, 
								       const Vec<long> &type, long prec): 
  ZZ_p_block_sylvester(type, prec) {
  init(f, type, prec);
}
