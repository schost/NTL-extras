#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* points in geometric progression                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* computes the polynomial                                   */
/*     f = (x-1)(x-a)(x-a^2)...(x-a^(n-1))                   */
/*-----------------------------------------------------------*/
void fan_in(zz_pX& f, const zz_p& a, long n){
  if (n <= 0){
    cerr << "too few points : n <= 0 in fanIn.\n";
    exit(-1);
  }
  long j, i = 1;

  while (i <= n/2)
    i<<=1;
  i>>=1;
  
  SetX(f);
  SetCoeff(f, 0, to_zz_p(-1));

  zz_pX tmp;
  zz_p b, c, d;
  b = a;
  c = a;

  // the algorithm is a basic divide-and-conquer
  while (i >= 1){
    tmp.rep.SetLength(deg(f)+1);
    d = c;
    for (long k = 0; k < tmp.rep.length(); k++){
      mul(tmp.rep[k], f.rep[k], d);
      div(d, d, b);
    }
    mul(f, f, tmp); // the cost is here: 1 polynomial multiplication
    mul(b, b, b);
    mul(c, c, c);
    mul(c, c, c);

    // accomodate the odd cases
    j = n&i;
    if (j != 0){
      long de = f.rep.length();
      f.rep.SetLength(de+1);
      f.rep[de] = f.rep[de-1];
      for (long k = de-1; k > 0; k--){
	mul(f.rep[k], f.rep[k], b);
	sub(f.rep[k], f.rep[k-1], f.rep[k]);
      }
      mul(f.rep[0], f.rep[0], -b);
      mul(c, c, b);
      mul(c, c, b);
      mul(c, c, a);
      mul(b, b, a);
    }
    i >>= 1;
  }
}

/*-----------------------------------------------------------*/
/* Let P=(x-1)(x-a)(x-a^2)...(x-a^(n-1))                     */
/* Evaluates P' on 1,a,...,a^(n-1)                           */
/* Not so clever: we should certainly not use divisions.     */
/*-----------------------------------------------------------*/
void derivative_fan_in(vec_zz_p& der, const zz_p& a, long n){
  vec_zz_p P, Q, powA;
  zz_p tmp;

  P.SetLength(n);
  Q.SetLength(n);
  powA.SetLength(n);
  der.SetLength(n);

  // powers of A
  powA[0] = to_zz_p(1);
  for (long i = 1; i < n; i++)
    mul(powA[i], powA[i-1], a);

  // small recurrence formulas give the result
  P[0] = to_zz_p(1);
  for (long i = 1; i < n; i++){
    mul(P[i], P[i-1], powA[i-1]);
    sub(tmp, powA[i], to_zz_p(1));
    mul(P[i], P[i], tmp);
  }
    
  Q[n-1] = to_zz_p(1);
  for (long i = n-2; i >= 0; i--){
    sub(tmp, powA[i], powA[n-1]);
    mul(Q[i], Q[i+1], tmp);
    div(Q[i], Q[i], powA[n-2-i]);
  }

  for (long i = 0; i < n; i++)
    mul(der[i], P[i], Q[i]);
}


/*-----------------------------------------------------------*/
/* initializes all quantities                                */
/*-----------------------------------------------------------*/
zz_pX_Multipoint_Geometric::zz_pX_Multipoint_Geometric(const zz_p& q, long eval_n, long interp_n, long i0){

  zz_p rho = q*q;
  this->q = rho;
  this->i0 = i0;
  
  // ----------------------------------
  // precomputation for evaluation
  // ----------------------------------
  this->eval_n = eval_n;
  this->eval_inverse_powers_square_q.SetLength(eval_n);
  this->eval_inverse_powers_square_q_shifted.SetLength(eval_n);

  // sets S=sum_{0 <= i < 2n-1} q^{i^2} x^i
  // sets inverse_powers_square_q
  this->eval_S.rep.SetLength(2*eval_n-1);
  zz_p power_q = to_zz_p(1);
  eval_S.rep[0] = to_zz_p(1);
  this->eval_inverse_powers_square_q[0] = to_zz_p(1);
  for (long i = 1; i < 2*eval_n-1; i++){
    zz_p tmp;
    tmp = eval_S.rep[i-1]*power_q;
    power_q =  q*power_q;
    eval_S.rep[i] = tmp*power_q;
    if (i < eval_n)
      this->eval_inverse_powers_square_q[i] = eval_S.rep[i];
  }
  this->eval_S.normalize();
  inv(this->eval_inverse_powers_square_q);

  // sets inverse_powers_square_q_shifted
  zz_p shift = power(rho, i0);
  zz_p tmp_shift = to_zz_p(1);
  for (long i = 0; i < eval_n; i++){
    this->eval_inverse_powers_square_q_shifted[i] = this->eval_inverse_powers_square_q[i]*tmp_shift;
    tmp_shift = tmp_shift*shift;
  }

  // finally, does the FFT
  long k = NextPowerOfTwo(2*eval_n-1);
  TofftRep(this->eval_S_FFT, this->eval_S, k);

  // ----------------------------------
  // precomputation for interpolation
  // ----------------------------------
  this->n = interp_n;
  this->inverse_powers_square_q.SetLength(interp_n);
  this->inverse_powers_square_q_shifted.SetLength(interp_n);
  this->inverse_derivative.SetLength(interp_n);

  // sets M as the master polynomial of rho^i = q^(2*i)
  // sets revM as the reverse of it
  fan_in(this->M, rho, interp_n);
  reverse(this->revM, this->M);

  // sets S=sum_{0 <= i < 2n-1} q^{i^2} x^i
  // sets inverse_powers_square_q
  this->S.rep.SetLength(2*interp_n-1);
  power_q = to_zz_p(1);
  S.rep[0] = to_zz_p(1);
  this->inverse_powers_square_q[0] = to_zz_p(1);
  for (long i = 1; i < 2*interp_n-1; i++){
    zz_p tmp;
    tmp = S.rep[i-1]*power_q;
    power_q =  q*power_q;
    S.rep[i] = tmp*power_q;
    if (i < interp_n)
      this->inverse_powers_square_q[i] = S.rep[i];
  }
  this->S.normalize();
  inv(this->inverse_powers_square_q);

  // sets inverse_derivative
  derivative_fan_in(this->inverse_derivative, rho, interp_n);
  inv(this->inverse_derivative);
  shift = power(rho, i0);
  zz_p inv_shift = 1/shift;
  zz_p tmp_inv_shift = to_zz_p(1);
  for (long i = 0; i < interp_n; i++){
    zz_p tmp;
    tmp = this->inverse_derivative[i]*this->inverse_powers_square_q[i];
    this->inverse_derivative[i] = tmp*tmp_inv_shift;
    tmp_inv_shift = tmp_inv_shift*inv_shift;
  }

  // does the FFT_s
  k = NextPowerOfTwo(2*interp_n-1);
  TofftRep(this->S_FFT, S, k);
  TofftRep(this->revM_FFT, revM, k);

  // prepares the matrix if size not too large and the ratio justifies it
  if (eval_n < 200 && n > 2*eval_n){
    Mat<zz_p> mat1, mat2;
    mat1.SetDims(eval_n, eval_n);
    mat2.SetDims(eval_n, eval_n);
    zz_p tmp1, tmp2;
    tmp1 = power(rho, i0);
    tmp2 = power(rho, i0+eval_n);
    
    for (long i = 0; i < eval_n; i++){
      zz_p ump1 = to_zz_p(1), ump2 = to_zz_p(1);
      for (long j = 0; j < eval_n; j++){
	mat1[i][j] = ump1;
	ump1 *= tmp1;
	mat2[i][j] = ump2;
	ump2 *= tmp2;
      }
      tmp1 *= rho;
      tmp2 *= rho;
    }
    inv(mat1, mat1);
    mul(shift_matrix, mat2, mat1);
  }
}

/*-----------------------------------------------------------*/
/* initializes all quantities                                */
/*-----------------------------------------------------------*/
zz_pX_Multipoint_Geometric::zz_pX_Multipoint_Geometric(const zz_p& q, long n, long i0) : zz_pX_Multipoint_Geometric(q, n, n, i0) { 
}

/*-----------------------------------------------------------*/
/* values[i] = f(q^(2(i+i0))), i = 0..eval_n-1               */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate_one(Vec<zz_p>& val, const zz_pX& f) const{
  val.SetLength(eval_n);
  
  zz_pX T;
  T.rep.SetLength(eval_n);
  
  for (long i = 0; i < eval_n; i++)
    T.rep[i] = coeff(f, eval_n-1-i) * eval_inverse_powers_square_q_shifted[eval_n-1-i];
  T.normalize();
  
  if (eval_n > NTL_zz_pX_MUL_CROSSOVER){
   long k = NextPowerOfTwo(2*eval_n-1);
    fftRep T_FFT(INIT_SIZE, k);
    TofftRep(T_FFT, T, k);
    mul(T_FFT, T_FFT, eval_S_FFT);
    FromfftRep(T, T_FFT, eval_n-1, 2*eval_n-2);
  }
  else{
    T = middle_product(eval_S, T, eval_n);
  }
  
  for (long i = 0; i < eval_n; i++)
    val[i] = coeff(T, i) * eval_inverse_powers_square_q[i];
}

/*-----------------------------------------------------------*/
/* values[i] = f(q^(2(i+i0))), i = 0..n-1                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate(Vec<zz_p>& val, const zz_pX& f_init) const {
  zz_pX f = f_init;

  if (n == 1){
    val.SetLength(1);
    val[0] = coeff(f, 0);
    return;
  }

  zz_p qpow = power(q, eval_n);
  val.SetLength(0);
  while(val.length() < n){
    Vec<zz_p> tmp;
    evaluate_one(tmp, f);
    val.append(tmp);
    zz_p shift = to_zz_p(1);
    for (long i = 0; i <= deg(f); i++){
      f.rep[i] *= shift;
      shift *= qpow;
    }    
  }

  val.SetLength(n);
}

/*-----------------------------------------------------------*/
/* values[i] = f(q^(2(i+i0))), i = 0..n-1                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f_init) const{
  
  cerr << "geom eval: ";
  double t = GetTime();

  if (shift_matrix.NumRows() == 0){
    val.SetLength(f_init.length());
    for (long i = 0; i < f_init.length(); i++)
      evaluate(val[i], f_init[i]);  
    cerr << GetTime()-t << " using shift " << endl;
    return;
  }

  long d = f_init.length();
  Mat<zz_p> coeff_f;
  coeff_f.SetDims(eval_n, d);
  val.SetLength(d);

  Vec<zz_p> tmp;
  tmp.SetLength(eval_n);
  for (long j = 0; j < d; j++){
    val[j].SetLength(n);
    evaluate_one(tmp, f_init[j]);
    for (long i = 0; i < eval_n; i++){
      val[j][i] = tmp[i];
      coeff_f[i][j] = tmp[i];;
    }
  }
    
  long idx = eval_n;
  while (idx < n){
    mul(coeff_f, shift_matrix, coeff_f);
    for (long i = 0; i < eval_n; i++)
      if (i+idx < n)
	for (long j = 0; j < d; j++)
	  val[j][i+idx] = coeff_f[i][j];
    idx += eval_n;
  }

  cerr << GetTime()-t << " using matrix " << endl;
}

/*------------------------------------------------------------*/
/* finds P of degree < n such that                            */
/* values[i] = P(q^(2(i+i0))), i = 0..n-1                     */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::interpolate(zz_pX& f, const vec_zz_p& val) const {

  if (n == 1){
    f.rep.SetLength(1);
    f.rep[0] = val[0];
    f.normalize();
    return;
  }

  zz_pX PV;
  PV.rep.SetLength(n);
  for (long i = 0; i < n; i++)
    PV.rep[i] = val[i];
  PV.normalize();

  if (n > NTL_zz_pX_MUL_CROSSOVER){
     long k = NextPowerOfTwo(2*n-1);
     fftRep PV_FFT(INIT_SIZE, k);
     TofftRep(PV_FFT, PV, k);
     mul(PV_FFT, PV_FFT, revM_FFT);
     FromfftRep(PV, PV_FFT, 0, n-1);
  }
  else{
    MulTrunc(PV, PV, revM, n);
  }

  zz_pX R;
  R.rep.SetLength(n);
  for (long i = 0; i < n; i++)
    R.rep[i] = coeff(PV, i) * inverse_powers_square_q[n-1-i];
  R.normalize();

  if (n > NTL_zz_pX_MUL_CROSSOVER){
     long k = NextPowerOfTwo(2*n-1);
     fftRep R_FFT(INIT_SIZE, k);
     TofftRep(R_FFT, R, k);
     mul(R_FFT, R_FFT, S_FFT);
     FromfftRep(R, R_FFT, n-1, 2*n-2);
  }
  else{
    R = middle_product(S, R, n);
  }

  f.rep.SetLength(n);
  for (long i = 0; i < n; i++)
    f.rep[i] = coeff(R, i) * inverse_derivative[i];
  f.normalize();
}

/*------------------------------------------------------------*/
/* a naive conversion to a dense matrix                       */
/* maybe promote to all multipoint matrices?                  */
/*------------------------------------------------------------*/
void to_dense(Mat<zz_p>& M, const zz_pX_Multipoint_Geometric& X){
   long n = X.n;
   zz_pX f;
   set(f);
   Vec<zz_p> v;
   v.SetLength(n);
   M.SetDims(n, n);
   
   for (long i = 0; i < n; i++){
     X.evaluate(v, f);
     for (long j = 0; j < n; j++)
       M[j][i] = v[j];
     f <<= 1;
  }
}
