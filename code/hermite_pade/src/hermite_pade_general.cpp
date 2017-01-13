#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"
#include "ZZ_hermite_pade.h"
#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
/* hermite pade for general polynomials                           */
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

/*---------------------------------------------------------*/
/* tries to increase the rank of M by adding random blocks */
/*---------------------------------------------------------*/
Vec<hankel> hermite_pade_general::increase_rank(long add){
  // cout << "adding " << add << endl;
  rows_added = add;
  Vec<hankel> result;
  for (long i = 0; i < type.length(); i++){
    Vec<long> adding;
    Vec<long> values;
    long to_add = max(1, add - type[i]);
    adding.SetLength(add+type[i],0);
    for (long j = 0; j < to_add; j++){
      values.append(rand() % 100);
      adding[j+add-to_add] = values[j];
    }
    Vec<zz_p> add_p = conv<Vec<zz_p>>(adding);
    hankel h{add_p, add, type[i]+1};
    vec_added.append(values);
    result.append(h);
  }
  return result;
}

/*----------------------------------------------------------------*/
/* creates a new block Sylvester matrix                           */
/*----------------------------------------------------------------*/
SmartPtr<ZZ_p_block_sylvester> hermite_pade_general::create_bmc(){
  Vec<ZZ_pX> fs_p;
  conv(fs_p, vec_fs);
  return MakeSmart<ZZ_p_block_sylvester_general>(ZZ_p_block_sylvester_general(fs_p, type, original_sizeX));
}

/*----------------------------------------------------------------*/
/* does the product by M and the extra rows                       */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade_general::mul_M_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> upper = hermite_pade::mul_M_right(b);
  if (rows_added != 0){
    Vec<Vec<ZZ_p>> b_split = split_on_type(flip_on_type(b));
    ZZ_pX lower;
    for (long i = 0; i < type.length(); i++){
      Vec<ZZ_p> values;
      for (long j = vec_added[i].length()-1; j >= 0; j--)
				values.append(ZZ_p(vec_added[i][j]));
      lower += conv<ZZ_pX>(values) * conv<ZZ_pX>(b_split[i]);
    }
    upper.SetLength(original_sizeX);
    //cout << "upper: " << upper << endl;
    //cout << "lower: " << lower << endl;
    for (long i = 0; i < deg(lower)+1; i++){
    	//cout << "i: " << i << endl;
      upper.append(lower[i]);
      //cout << "upper: " << upper[i+original_sizeX] << endl;
    }
  }
  return upper;
}

/*----------------------------------------------------------------*/
/* fs: the power series                                           */
/* type: the type of the approximant                              */
/* sigma: requested precision                                     */
/*----------------------------------------------------------------*/
hermite_pade_general::hermite_pade_general(const Vec<ZZX> &fs, const Vec<long> &type, long sigma, long fft_init): 
  hermite_pade(fft_init){

  this->type = type;
  vec_fs = fs;
  
  // setting up the mosaic Hankel Matrix
  Vec<hankel> vec_H;
  for (long i = 0; i < fs.length(); i++){
    Vec<ZZ> v_ZZ = conv<Vec<ZZ>>(fs[i]);
    Vec<zz_p> v = conv<Vec<zz_p>>(v_ZZ);
    v.SetLength(sigma);
    Vec<zz_p> inp;
    for (long j = 0; j < v.length(); j++)
      inp.append(v[v.length() - 1 - j]);
    for (long j = 0; j< type[i]; j++)
      inp.append(zz_p{0});
    vec_H.append(hankel(inp, sigma, type[i]+1)); 
  }
  Vec<Vec<hankel>> hankel_matrices;
  hankel_matrices.append(vec_H);
  mosaic_hankel MH(hankel_matrices);
  
  lzz_p_cauchy_like_geometric CL; // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int; // preconditioners
  Vec<zz_p> e_zz_p, f_zz_p; // the diagonal matrices
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA, CL); // inverting M mod p
  cout << "original rank: " << rank << endl;
  sizeX = X_int.length();
  sizeY = Y_int.length();
  original_sizeX = sizeX;
  
  //Mat<zz_p> mat;
  
  if (sizeY -rank != 1){
    hankel_matrices.append(increase_rank(sizeY-rank-1));
    MH = mosaic_hankel(hankel_matrices);
    //to_dense(mat, MH);
    //cout << "new mat: " << mat << endl;
  }
  else
    rows_added = 0;

  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA, CL); // inverting M mod p
  cout << "new rank: " << rank << endl;
  sizeX = X_int.length();
  sizeY = Y_int.length();
  //CL.to_dense(mat);
  //cout << "CL: " << mat << endl;

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

  vec_M.SetLength(0);
  vec_X_int.SetLength(1);
  vec_Y_int.SetLength(1);
  vec_X_int[0] = X_int_ZZ_p;
  vec_Y_int[0] = Y_int_ZZ_p;
  set_up_bmc();
}




