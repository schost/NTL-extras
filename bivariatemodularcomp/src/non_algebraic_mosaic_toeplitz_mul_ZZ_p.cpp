#include "non_algebraic_mosaic_toeplitz_mul_ZZ_p.h"
using namespace NTL;

void non_algebraic_mosaic_toeplitz_mul_ZZ_p::init(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec){
  this->type = type;
  this->prec = prec;
  fs_field = fs;
  initialized = true;
}

non_algebraic_mosaic_toeplitz_mul_ZZ_p::non_algebraic_mosaic_toeplitz_mul_ZZ_p
(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec):
mosaic_toeplitz_mul_ZZ_p(type,prec){
  init(fs,type,prec);
}

Vec<ZZ_p> non_algebraic_mosaic_toeplitz_mul_ZZ_p::mult_right(const Vec<ZZ_p> &rhs){
	if (!initialized) throw "must init first";
	Vec<ZZ_pX> rhs_poly;
	create_lhs_list(rhs_poly,rhs);
	ZZ_pX result;
	for (long i = 0; i < fs_field.length(); i++)
		result += trunc(fs_field[i] * rhs_poly[i],prec);
	Vec<ZZ_p> v;
	v.SetLength(prec);
	for (long i = 0; i < prec; i++)
		v[i] = coeff(result, i);
	return v;
}
