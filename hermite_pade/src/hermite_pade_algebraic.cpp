#include <cstdlib>
#include <ctime>
#include <cmath>

#include "hermite_pade_algebraic.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"

Vec<long> hermite_pade_algebraic::update_type(){
  if (rank + 1 == sizeY) return type;
	switch_context(0);
	
	cout << "sizeY: " << sizeY << " rank: " << rank << endl;
	
  Vec<ZZ_p> ex1; // mult with A to get a column
  ex1.SetLength(sizeY, ZZ_p(0));
  Vec<ZZ_p> add1;
  add1.SetLength(sizeY - rank);
  for (long i = 0; i < add1.length(); i++){
  	add1[i] = rand() % 100 + 1;
  	ex1[rank+i] = add1[i];
  }
  
  Vec<ZZ_p> ex2; // mult with A to get a column
  Vec<ZZ_p> add2;
  add2.SetLength(sizeY-rank);
  ex2.SetLength(sizeY, ZZ_p(0));
  for (long i = 0; i < add2.length(); i++){
  	add2[i] = rand() % 100 + 1;
  	ex2[rank+i] = add2[i];
  }
	
	ex1 = mulA_right(ex1);
	ex2 = mulA_right(ex2);
	Vec<ZZ> sol1,sol2;
	DAC(sol1,conv<Vec<ZZ>>(ex1),0);
	DAC(sol2,conv<Vec<ZZ>>(ex2),0);
	ex1.kill();
	ex2.kill();
	conv(ex1,sol1);
	conv(ex2,sol2);
	ex1.SetLength(sizeY, ZZ_p(0));
  ex2.SetLength(sizeY, ZZ_p(0));
  
	for (long i = 0; i < add1.length(); i++){
  	ex1[rank+i] = -add1[i];
  }
  for (long i = 0; i < add2.length(); i++){
  	ex2[rank+i] = -add2[i];
  }
	ex1 = flip_on_type(mul_Y_right(ex1));
	ex2 = flip_on_type(mul_Y_right(ex2));

	Vec<zz_p> v1 = conv<Vec<zz_p>>(conv<Vec<ZZ>>(ex1));
	Vec<zz_p> v2 = conv<Vec<zz_p>>(conv<Vec<ZZ>>(ex2));
	
	zz_pXY b1(split_on_type(v1));
	zz_pXY b2(split_on_type(v2));
	zz_pXY gcd;
	GCD(gcd,b1,b2);
	cout << "GCD: " << gcd << endl;
	// extracting the type
	Vec<long> t;
	long total = 0;
	for (long i = 0; i < gcd.coeffX.length(); i++){
	  long degree = deg(gcd.coeffX[i]);
	  t.append(degree);
	  total += degree+1;
	}
  return t;
}

Vec<zz_pX> hermite_pade_algebraic::split_on_type(const Vec<zz_p> &v){
	long acc = 0;
	Vec<zz_pX> result;
	for (long i = 0; i < type.length(); i++){
		Vec<zz_p> f;
		for (long j = 0; j < type[i]+1; j++){
		  f.append(v[acc+j]);
		}
		acc += type[i] + 1;
		result.append(conv<zz_pX>(f));
	}
	return result;
}


void hermite_pade_algebraic::init(const ZZX &f, const Vec<long> &type, long prec_inp, long fft_index){
  set_up_field(fft_index);
  
  long prec = deg(f) + 1;
  f_full_prec = f;
  this->type = type;
  zz_pX f_field;
  conv(f_field, f);
  this->prec = prec;
  level = 0;
  
  // setting up the mosaic Hankel matrix
  Vec<hankel> vec_H;
  zz_pX running {1};
  for (long i = 0; i < type.length(); i++){
    Vec<zz_p> v;
    Vec<zz_p> inp_vec;
    conv(v,running);
    v.SetLength(prec, zz_p{0});
    for (int j = 0; j < v.length(); j++)
      inp_vec.append(v[v.length() - 1 - j]);
    for (int j = 0; j < type[i]; j++)
        inp_vec.append(zz_p{0});
    vec_H.append(hankel(inp_vec, prec, type[i] + 1));
    running = running * f_field;
  }
  Vec<Vec<hankel>> hankel_matrices;
  hankel_matrices.append(vec_H);
  mosaic_hankel MH(hankel_matrices);
  
  cauchy_like_geometric_special CL; // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int; // preconditioners
  
  // find the rank of the original matrix
  Vec<zz_p> e_zz_p, f_zz_p; // the diagonal matrices
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA, CL); // inverting M mod p

  sizeX = X_int.length();
  sizeY = Y_int.length();
	
  // initializing the bivariate modular comp
  ZZ_pX f_ZZ_pX;
  conv(f_ZZ_pX, f);
  BivariateModularComp M(f_ZZ_pX, type, sizeX); // could pass in the precomputed stuff
  // initializing the pointer variables and vectors
  vec_M.kill();
  vec_M.SetLength(1);
  vec_M[0] = M;

  // converting the preconditioners that do not change
  this->e = conv<Vec<ZZ>>(e_zz_p);
  this->f = conv<Vec<ZZ>>(f_zz_p);
  zz_p c_zz, d_zz; 
  X_int.point(c_zz, 0);
  Y_int.point(d_zz, 0); 
  c = conv<ZZ>(c_zz);
  d = conv<ZZ>(d_zz);

  // initializing the X_int and Y_int stuff
  zz_p w_zz_p, w2;
  X_int.point(w_zz_p, 1);
  Y_int.point(w2, 1);
  w_zz_p = w_zz_p / c_zz;
  this->w = w_zz_p.LoopHole();
  ZZ_p w_p = conv<ZZ_p>(this->w);
  // find the order of w
  order = order_dyadic(w_zz_p);

  ZZ_pX_Multipoint_FFT X_int_ZZ_p(w_p, conv<ZZ_p>(this->c), sizeX);
  ZZ_pX_Multipoint_FFT Y_int_ZZ_p(w_p, conv<ZZ_p>(this->d), sizeY);

	vec_X_int.kill();
	vec_Y_int.kill();
	vec_X_int.SetLength(1);
	vec_Y_int.SetLength(1);
  vec_X_int[0] = X_int_ZZ_p;
  vec_Y_int[0] = Y_int_ZZ_p;
}

// todo: f (the input poly) != f (the ZZ attribute)
hermite_pade_algebraic::hermite_pade_algebraic(const ZZX &f, const Vec<long>& type, long prec_inp, long fft_index){
  init(f,type,prec_inp,fft_index);
}

void hermite_pade_algebraic::set_up_bmc(){
	ZZ_pX f_p;
  conv(f_p, f_full_prec);
  BivariateModularComp m_new(f_p, type, rank);
  vec_M.append(m_new);
}

























