#include "hermite_pade.h"

hermite_pade::hermite_pade(const ZZX &f, const Vec<long>& type){
  long prec = deg(f) + 1;
  f_full_prec = f;
  this->type = type;
  zz_pX f_field;
  conv(f_field, f);
    
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
    vec_H.append(hankel(inp_vec,prec,type[i]+1));
    running = running * f_field;
  }
  Vec<Vec<hankel>> hankel_matrices;
  hankel_matrices.append(vec_H);
  mosaic_hankel MH(hankel_matrices);

  // setting up the Cauchy matrix
  cauchy_like_geometric_special CL; // the cauchy matrix
  zz_pX_Multipoint_Geometric X_int, Y_int; // preconditioners
  Vec<zz_p> e_zz_p, f_zz_p; // the diagonal matrices
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert_fast(invM,CL); // inverting M mod p

  // initializing the bivariate modular comp
  ZZ_pX f_ZZ_pX;
  conv(f_ZZ_pX, f);
  BivariateModularComp M(f_ZZ_pX, type, rank); // could pass in the precomputed stuff

  // initialize the pointer variables and vectors
  vec_M.append(M);
  Vec<ZZ> e_ZZ, f_ZZ;
  conv(e_ZZ, e_zz_p);
  conv(f_ZZ, f_zz_p);
  conv(this->e, e_ZZ);
  conv(this->f, f_ZZ);

}










