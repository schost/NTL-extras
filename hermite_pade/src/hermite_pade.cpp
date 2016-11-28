#include "hermite_pade.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"
using namespace NTL;

hermite_pade::~hermite_pade(){
  for (long i = 0; i < vec_M.length(); i++) delete vec_M[i];
}

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
    set_up_bmc();
    
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
    
    vec_X_int.append(X_new);
    vec_Y_int.append(Y_new);
    level = n;
  }
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix Y_int^t D_f                         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_Y_right(const Vec<ZZ_p> &b){
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
  Vec<ZZ_p> x = mul_Y_right(b_loc);

  x = mul_M_right(x);

  return mul_X_right(x);
}


Vec<ZZ_p> hermite_pade::mul_M_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> flipped = flip_on_type(b);
  return vec_M[level]->mult_right(flipped);
}

Vec<ZZ_p> hermite_pade::mul_X_right(Vec<ZZ_p> b){
	vec_X_int[level].mul_right(b, b);
  Vec<ZZ_p> e = conv<Vec<ZZ_p>>(this->e);
  Vec<ZZ_p> x;
  mul_diagonal(x, e, b);
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

void hermite_pade::set_up_field(long fft_index){
	zz_p::FFTInit(fft_index);
  ctx = zz_pContext(INIT_FFT, fft_index);
	srand(time(NULL));
  p = zz_p::modulus();
  ZZ_p::init(ZZ(p));
  p_powers.append(ZZ(p));
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
  if (n == 0) 
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

Vec<Vec<ZZ_p>> hermite_pade::split_on_type(const Vec<ZZ_p> &v){
	long acc = 0;
	Vec<Vec<ZZ_p>> result;
	for (long i = 0; i < type.length(); i++){
		Vec<ZZ_p> f;
		for (long j = 0; j < type[i]+1; j++){
		  f.append(v[acc+j]);
		}
		acc += type[i] + 1;
		result.append(f);
	}
	return result;
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
  Vec<ZZ_p> soln = mul_Y_right(x);

  // loop until we get enough prec
  while(!can_reconstruct(soln, n)){
    switch_context(++n);
    cout << "n: " << n << endl;
    
    Vec<ZZ_p> extractor; // mult with A to get a column
    extractor.SetLength(sizeY, ZZ_p(0));
    extractor[rank] = 1; // just for now, take the last column
    Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A

    update_solution(x_ZZ, b, n);

    Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
    x.SetLength(sizeY,ZZ_p(0));
    x[rank] = -1;
    soln.kill();
    soln = mul_Y_right(x);

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
  
  /***********************************************
  ZZ_p::init(ZZ(13));
  Vec<ZZ_p> check;
  for (long i = 0; i < sol.length(); i++){
    check.append(conv<ZZ_p>(sol[i][0]) / conv<ZZ_p>(sol[i][1]));	
  }
  
  cout << "check: " << mul_M_right(check) << endl;
  
  
  ***********************************************/
  
}
























