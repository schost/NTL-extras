#include "hermite_pade_non_algebraic.h"
#include "hermite_pade_algebraic.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"
using namespace NTL;

Vec<hankel> hermite_pade_non_algebraic::increase_rank(){
	return {};
}


hermite_pade_non_algebraic::hermite_pade_non_algebraic
(const Vec<ZZX> &fs, const Vec<long> &type, long fft_init): hermite_pade(fft_init){
	this->type = type;
	vec_fs = fs;
	
	// find the highest degree
	long max_deg = 0;
	for (long i = 0; i < fs.length(); i++)
	  if (max_deg < deg(fs[i])) max_deg= deg(fs[i]);
	
	// setting up the mosaic Hankel Matrix
	Vec<hankel> vec_H;
	for (long i = 0; i < fs.length(); i++){
		Vec<ZZ> v_ZZ = conv<Vec<ZZ>>(fs[i]);
	  Vec<zz_p> v = conv<Vec<zz_p>>(v_ZZ);
	  v.SetLength(max_deg+1);
	  Vec<zz_p> inp;
	  for (long j = 0; j < v.length(); j++)
	    inp.append(v[v.length() - 1 - j]);
	  for (long j = 0; j< type[i]; j++)
	    inp.append(zz_p{0});
	  vec_H.append(hankel(inp, max_deg+1, type[i]+1)); 
	}
	Vec<Vec<hankel>> hankel_matrices;
	hankel_matrices.append(vec_H);
	mosaic_hankel MH(hankel_matrices);
	Mat<zz_p> mat;
	to_dense(mat,MH);
	cout << mat << endl;
	
	cauchy_like_geometric_special CL; // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int; // preconditioners
  Vec<zz_p> e_zz_p, f_zz_p; // the diagonal matrices
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA, CL); // inverting M mod p
  sizeX = X_int.length();
  sizeY = Y_int.length();
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
  
  vec_M.kill();
  set_up_bmc();
  Vec<ZZ_p> ex;
  ex.SetLength(sizeY,ZZ_p{0});
  ex[sizeY-1] = 1;
  cout << "rank: " << rank << endl;
  to_dense(mat,CL);
  cout << mat << endl;
  cout << ex << endl;
  cout << mulA_right(ex) << endl;
}


void hermite_pade_non_algebraic::set_up_bmc(){
	Vec<ZZ_pX> fs_p;
	conv(fs_p, vec_fs);
	BivariateModularComp m_new(fs_p, type, sizeX);
	vec_M.append(m_new);
}

Vec<ZZ_p> hermite_pade_non_algebraic::mul_M_right(const Vec<ZZ_p> &b){
	return hermite_pade::mul_M_right(b);
}
