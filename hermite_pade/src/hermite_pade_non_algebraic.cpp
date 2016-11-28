#include "hermite_pade_non_algebraic.h"
#include "hermite_pade_algebraic.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"
#include "non_algebraic_mosaic_toeplitz_mul_ZZ_p.h"
#include <stdlib.h>
#include <time.h>
using namespace NTL;

Vec<hankel> hermite_pade_non_algebraic::increase_rank(long add){
  cout << "adding " << add << endl;
  rows_added = add;
  Vec<hankel> result;
  for (long i = 0; i < type.length(); i++){
    cout << "type[i]: " << type[i] << endl;
    Vec<long> adding;
    Vec<long> values;
    long to_add = max(1, add - type[i]);
    adding.SetLength(add+type[i],0);
    cout << "to add: " << to_add << endl;
    for (long j = 0; j < to_add; j++){
      values.append(rand()%100+1);
      cout << "added: " << values[j] << endl;
      adding[j+add-to_add] = values[j];
    }
    Mat<zz_p> mat;
    Vec<zz_p> add_p = conv<Vec<zz_p>>(adding);
    hankel h{add_p, add, type[i]+1};
    to_dense(mat,h);
    cout << mat << endl;
    vec_added.append(values);
    result.append(h);
  }
	return result;
}


hermite_pade_non_algebraic::hermite_pade_non_algebraic
(const Vec<ZZX> &fs, const Vec<long> &type, long fft_init): hermite_pade(fft_init){
	srand(time(NULL));
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
  original_sizeX = sizeX;
  
  cout << "rank: " << rank << endl;
  
  // check if we need to add more rows
  if (sizeY -rank != 1){
    cout << "ADD MORE ROWS" << endl;
    hankel_matrices.append(increase_rank(sizeY-rank-1));
    MH = mosaic_hankel(hankel_matrices);
    to_dense(mat,MH);
		cout << mat << endl;
  }
  
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
	vec_M.append(new non_algebraic_mosaic_toeplitz_mul_ZZ_p(fs_p, type, sizeX));
}

Vec<ZZ_p> hermite_pade_non_algebraic::mul_M_right(const Vec<ZZ_p> &b){
	Vec<ZZ_p> upper = hermite_pade::mul_M_right(b);
	Vec<Vec<ZZ_p>> b_split = split_on_type(flip_on_type(b));
	ZZ_pX lower;
	if (rows_added != 0){
	for (long i = 0; i < type.length(); i++){
		  Vec<ZZ_p> values;
		  for (long j = vec_added[i].length()-1; j >= 0; j--)
		    values.append(ZZ_p(vec_added[i][j]));
		  lower += conv<ZZ_pX>(values) * conv<ZZ_pX>(b_split[i]);
		}
		for (long i = 0; i < deg(lower)+1; i++){
	 	 upper[i+original_sizeX] = lower[i];
		}
	}
	return upper;
}


















