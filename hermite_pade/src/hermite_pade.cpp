#include <cstdlib>
#include <ctime>
#include <cmath>

#include "hermite_pade.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_p_extra.h"

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::flip_on_type (const Vec<ZZ_p> &v){
  Vec<ZZ_p> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<Vec<ZZ>> hermite_pade::flip_on_type (const Vec<Vec<ZZ>> &v){
  Vec<Vec<ZZ>> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* returns a block of bi-sub-diagonal hankel matrices             */
/* size is determined by the type                                 */
/*----------------------------------------------------------------*/
Vec<hankel> hermite_pade::create_random_hankel_on_type(){
  Vec<hankel> mh;
  for (long i = 0; i < type.length(); i++){
    long cols = type[i] + 1;
    Vec<zz_p> running;
    running.SetLength(added + cols - 1, zz_p(0));
    long random1 = std::rand() % 1000;
    long random2 = std::rand() % 1000;
    running[added-1] = random1;
    if (added-1 != 0)
      running[added-2] = random2;
    diagonals1.append(random1);
    diagonals2.append(random2);
    mh.append(hankel(running, added, cols));
  }

  return mh;
}

/*----------------------------------------------------------------*/
/* computes Bv, where B is made of anti-diagonal matrices         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_bottom_mosaic_diagonal(const Vec<ZZ_p> &v){
  if (added == 0){ 
    Vec<ZZ_p> empty;		
    return empty;
  }
  long acc = 0;
  Vec<ZZ_p> x;
  x.SetLength(added, ZZ_p(0));
  for (long i = 0; i < type.length(); i++){
    for (long j = 0; j < std::min(added, type[i]+1); j++)
      x[j] += v[acc+j] * diagonals1[i];
    for (long j = 0; j < std::min(type[i]+1, added-1); j++)
      x[j+1] += v[acc+j] * diagonals2[i];
    acc += type[i] + 1;
  }
  return x;
}

// todo: f (the input poly) != f (the ZZ attribute)
hermite_pade::hermite_pade(const ZZX &f, const Vec<long>& type, long prec_inp, long fft_index){
  std::srand(std::time(0));
  zz_p::FFTInit(fft_index);
  ctx = zz_pContext(INIT_FFT, fft_index);

  p = zz_p::modulus();
  ZZ_p::init(ZZ(p));
  p_powers.append(ZZ(p));
  cout << "p: " << p << endl;
  long prec = deg(f) + 1;
  f_full_prec = f;
  this->type = type;
  zz_pX f_field;
  conv(f_field, f);
  this->prec = prec;
  level = 0;

  long type_sum = 0;
  for (long i = 0; i < type.length(); i++)
    type_sum += type[i] + 1;
  
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
  cout << "original rank: " << rank << endl;

  added = 0;
  if (rank < type_sum-1){
    added = this->prec - rank;
    hankel_matrices.append(create_random_hankel_on_type());
    MH = mosaic_hankel(hankel_matrices);
  }

  // setting up the Cauchy matrix
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA, CL); // inverting M mod p
  cout << "new rank: " << rank << endl;
  sizeX = X_int.length();
  sizeY = Y_int.length();
  cout << "number of cols: " << sizeY << endl;

  // initializing the bivariate modular comp
  ZZ_pX f_ZZ_pX;
  conv(f_ZZ_pX, f);
  BivariateModularComp M(f_ZZ_pX, type, rank); // could pass in the precomputed stuff
  // initializing the pointer variables and vectors
  vec_M.append(M);

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

  vec_X_int.append(X_int_ZZ_p);
  vec_Y_int.append(Y_int_ZZ_p);
}

void hermite_pade::switch_context(long n){
  if (n < vec_M.length()){ // it has already been computed
    ZZ_p::init(p_powers[n]);
    level = n;
  } 
  else{
    // calculating the new power of p
    ZZ p_new(p);
    long pow2 = power_long(2,n);
    power(p_new, p_new, pow2); // 2^n isn't going to be very large
    p_powers.append(p_new);
    ZZ_p::init(p_new);
    // creating the new bivariate modular comp
    ZZ_pX f_p;
    conv(f_p, f_full_prec);
    BivariateModularComp m_new(f_p, type, rank);
    // computing w mod p^2^n
    ZZ new_w;
    lift_root_of_unity(new_w, this->w, order, p, pow2);
    ZZ_p w_p, c_p, d_p;
    conv(w_p, new_w);
    conv(c_p, c);
    conv(d_p, d);
    ZZ_pX_Multipoint_FFT X_new (w_p, c_p, sizeX);
    ZZ_pX_Multipoint_FFT Y_new (w_p, d_p, sizeY);
    // update
    vec_M.append(m_new);
    vec_X_int.append(X_new);
    vec_Y_int.append(Y_new);
    level = n;
  }
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix Y_int^t D_f                         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::find_original_sol(const Vec<ZZ_p> &b){
  Vec<ZZ_p> x;
  Vec<ZZ_p> f = conv<Vec<ZZ_p>>(this->f); // TODO: e and f should be called e_ZZ and f_ZZ?
  mul_diagonal(x, f, b);
  vec_Y_int[level].mul_left(x, x);
  return x;
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f       */
/* (CL is cauchy-geometric-like)                                  */
/* b need not have size CL.NumCols()                              */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mulA_right(const Vec<ZZ_p>& b){
  Vec<ZZ_p> b_loc = b;
  b_loc.SetLength(sizeY, ZZ_p(0)); // pad it
  Vec<ZZ_p> x = find_original_sol(b_loc);

  // TODO: should there be a function for the product by M?
  Vec<ZZ_p> flipped = flip_on_type(x);
  x = vec_M[level].mult(flipped);
  Vec<ZZ_p> res = mul_bottom_mosaic_diagonal(flipped);
  for (long i = 0; i < res.length(); i++)
    x.append(res[i]);

  vec_X_int[level].mul_right(x, x);
  Vec<ZZ_p> e = conv<Vec<ZZ_p>>(this->e);
  mul_diagonal(x, e, x);
  return x;
}

/*----------------------------------------------------------------*/
/* if Mx = b mod p^(2^{n-1}), updates x so that Mx = b mod p^(2^n)*/
/*----------------------------------------------------------------*/
void hermite_pade::update_solution(Vec<ZZ>& x, const Vec<ZZ_p> &b, long n){
  Vec<ZZ_p> r = mulA_right(conv<Vec<ZZ_p>>(x)) - b; // mod p^(2^n)
  Vec<ZZ> r_ZZ = conv<Vec<ZZ>>(r);
  for (long i = 0; i < r.length(); i++)
    r_ZZ[i] = r_ZZ[i] / p_powers[n-1];
  Vec<ZZ> x_1;
  DAC(x_1, r_ZZ, n-1);
  x = x - p_powers[n-1] * x_1;
}

/*----------------------------------------------------------------*/
/* solves for Mx = b mod p^(2^n)                                  */
/*----------------------------------------------------------------*/
void hermite_pade::DAC(Vec<ZZ> &x, const Vec<ZZ>& b_in, long n){

  if (n == 0){ // we are mod p
    zz_pContext push;
    ctx.restore();
    Vec<zz_p> x_zz_p;
    Vec<zz_p> b_zz_p = conv<Vec<zz_p>>(b_in);
    b_zz_p.SetLength(rank); // TODO: why do we need this?
    mul_right(x_zz_p, invA, b_zz_p);
    x = conv<Vec<ZZ>>(x_zz_p);
    return;
  }

  long old_n = level;
  switch_context(n); 

  Vec<ZZ_p> b_n = conv<Vec<ZZ_p>>(b_in);
  Vec<ZZ> b = conv<Vec<ZZ>>(b_n); // this reduces b_in mod p^(2^n).
  DAC(x, b, n-1); 

  update_solution(x, b_n, n);
  switch_context(old_n); 
}

bool hermite_pade::can_reconstruct(const Vec<ZZ_p> &v, long n){
  if (n > 0) 
    return false;
  for (long i = 0; i < v.length(); i++){
    ZZ a,b;
    try{
      long result = ReconstructRational(a,b,conv<ZZ>(v[i]),p_powers[n], p_powers[n-1], p_powers[n-1]);
      if (result == 0) return false;
    }
    catch(...){
      return false;
    }
  }
  return true;
}

void hermite_pade::reconstruct(Vec<Vec<ZZ>> &sol, const Vec<ZZ_p> &v, long n){
  for (long i = 0; i < v.length(); i++){
    ZZ a, b;
    ReconstructRational(a, b, conv<ZZ>(v[i]), p_powers[n], p_powers[n-1], p_powers[n-1]);
    Vec<ZZ> temp;
    temp.append(a);
    temp.append(b);
    sol.append(temp);
  }
}

void hermite_pade::find_rand_sol(Vec<Vec<ZZ>> &sol){
  zz_pContext zz_p_push;
  zz_pContext ZZ_p_push;
  ctx.restore();
  long n = 0; // start at p^2^n
  switch_context(n);

  Vec<ZZ_p> extractor; // mult with A to get a column
  extractor.SetLength(sizeY, ZZ_p(0));
  extractor[rank] = 1; // just for now, take the last column
  Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A
  Vec<ZZ> x_ZZ, b_ZZ;
  b_ZZ = conv<Vec<ZZ>> (b);
  DAC(x_ZZ, b_ZZ, n); // solution mod p

  Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
  x.SetLength(sizeY, ZZ_p(0));
  x[rank] = -1;
  Vec<ZZ_p> soln = find_original_sol(x);

  // loop until we get enough prec
  while(!can_reconstruct(soln, n) && n < 4){
    switch_context(++n);
    
    Vec<ZZ_p> extractor; // mult with A to get a column
    extractor.SetLength(sizeY, ZZ_p(0));
    extractor[rank] = 1; // just for now, take the last column
    Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A

    update_solution(x_ZZ, b, n);

    Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
    x.SetLength(sizeY,ZZ_p(0));
    x[rank] = -1;
    soln.kill();
    soln = find_original_sol(x);

    ZZ_p first;
    for (long i = 0; i < soln.length(); i++)
      if (soln[i]._ZZ_p__rep % p_powers[0] != ZZ(0)){
        first = 1/soln[i];
        break;
      }
    soln *= first;
  }
  
  reconstruct(sol, soln, n);
  sol = flip_on_type(sol);
  for (long times = 0; times < 10; times++){
    long prime = RandomPrime_long(32);
    ZZ_p::init(ZZ(prime));
    Vec<ZZ_p> soln2;
    for (long i = 0; i < sol.length(); i++){
      ZZ_p a = conv<ZZ_p>(sol[i][0]);
      ZZ_p b = conv<ZZ_p>(sol[i][1]);
      soln2.append(a/b);
    }
    ZZ_pX f_p;
    conv(f_p, f_full_prec);
    BivariateModularComp m_new(f_p, type, rank);
    cout << "Double check with p = " << ZZ_p::modulus() << ": " << m_new.mult(soln2) << endl;
  }
}





