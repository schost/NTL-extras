#include "hermite_pade.h"

hermite_pade::hermite_pade(const ZZX &f, const Vec<long>& type, long prec_inp){
  long prec = deg(f) + 1;
  f_full_prec = f;
  this->type = type;
  zz_pX f_field;
  conv(f_field, f);
  this->prec = prec_inp;
    
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
  sizeX = X_int.length();
  sizeY = Y_int.length();

  // initializing the bivariate modular comp
  ZZ_pX f_ZZ_pX;
  conv(f_ZZ_pX, f);
  BivariateModularComp M(f_ZZ_pX, type, prec); // could pass in the precomputed stuff

  // initializing the pointer variables and vectors
  vec_M.append(M);
  this->M = &M;

  // coverting the preconditioners that donot change
  Vec<ZZ> e_ZZ, f_ZZ;
  conv(e_ZZ, e_zz_p);
  conv(f_ZZ, f_zz_p);
  conv(this->e, e_ZZ);
  conv(this->f, f_ZZ);
  zz_p c, d;
  X_int.point(c,0);
  Y_int.point(d,0);
  conv(this->c,c);
  conv(this->d,d);

  // initializing the X_int and Y_int stuff
  zz_p w_zz_p;
  X_int.point(w_zz_p,1);
  w_zz_p = w_zz_p / c;
  ZZ w;
  ZZ_p w_p;
  conv(w,w_zz_p);
  conv(w_p,w);
  conv(this->w,w);
  // find the order of w
  order = find_order(w_zz_p);
  ZZ_pX_Multipoint_FFT X_int_ZZ_p(w_p,conv<ZZ_p>(this->c), sizeX);
  ZZ_pX_Multipoint_FFT Y_int_ZZ_p(w_p,conv<ZZ_p>(this->d), sizeY);
  this->vec_X_int.append(X_int_ZZ_p);
  this->vec_Y_int.append(Y_int_ZZ_p);
  this->X_int = &X_int_ZZ_p;
  this->Y_int = &Y_int_ZZ_p;
}


void hermite_pade::switch_context(long n){
  if (n + 1 < vec_M.length()){ // it has already been computed
    n++; // since when n = 0, its actually p^1
    ZZ_p::init(p_powers[n]);
    M = &vec_M[n];
    X_int = &vec_X_int[n];
    Y_int = &vec_Y_int[n]; 
  } else{
    // calculating the new power of p
    ZZ p_new(p);
    long pow2 = power_long(2,n++);
    power(p_new, pow2); // 2^n isn't going to be very large
    ZZ_p::init(p_new);
    // creating the new bivariate modular comp
    ZZ_pX f_p;
    conv(f_p, f_full_prec);
    BivariateModularComp m_new(f_p, type, prec);
    // computing w mod p^2^n
    ZZ new_w;
    lift_root_of_unity(new_w, this->w, order, p, pow2);
    ZZ_p w_p, c_p, d_p;
    ZZ_pX_Multipoint_FFT X_new (w_p, c_p, sizeX);
    ZZ_pX_Multipoint_FFT Y_new (w_p, d_p, sizeY);
    // update
    vec_M.append(m_new);
    vec_X_int.append(X_new);
    vec_Y_int.append(Y_new);
    M = &vec_M[n];
    X_int = &vec_X_int[n];
    Y_int = &vec_Y_int[n];
  }
}


long hermite_pade::find_order(zz_p w){
  if (w == zz_p(1))
    return 1;
  else
    return 2 * find_order(power(w,2));
}

