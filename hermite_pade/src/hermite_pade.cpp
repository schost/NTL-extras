#include "hermite_pade.h"
#include <cstdlib>
#include <ctime>
#include <cmath>

Vec<ZZ_p> construct_diagonal (ZZ c, long n){
  ZZ_p c_p;
  conv(c_p,c);
  ZZ_p running(1);
  Vec<ZZ_p> v;
  for (long i = 0; i < n; i++){
    v.append(running);
    running = c_p * running;
  }
  return v;
}

Vec<ZZ_p> invert_diagonal (Vec<ZZ_p> c){
  Vec<ZZ_p> v;
  for (long i = 0; i < c.length(); i++)
    v.append((ZZ_p(1) / c[i]));
  return v;
}

Vec<ZZ_p> hermite_pade::flip_on_type (const Vec<ZZ_p> &v){
  Vec<ZZ_p> r;
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--){
      r.append(v[j+acc]);
    }
    acc += type[i] + 1;
  }
  return r;
}

Vec<Vec<ZZ>> hermite_pade::flip_on_type (const Vec<Vec<ZZ>> &v){
  Vec<Vec<ZZ>> r;
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--){
      r.append(v[j+acc]);
    }
    acc += type[i] + 1;
  }
  return r;
}

hankel hermite_pade::create_random_hankel (long rows, long cols){
  Vec<zz_p> running;
  cout << "rows: " << rows << " cols: " << cols << endl;
  cout << "R: " << rows+cols-1 << endl;
  running.SetLength(rows+cols-1, zz_p(0));
  long random1 = std::rand()%1000;
  long random2 = std::rand()%1000;
  running[rows-1] = random1;
  if (rows-1 != 0)
    running[rows-2] = random2;
  diagonals1.append(random1);
  diagonals2.append(random2);
  Mat<zz_p> m;
  to_dense(m, hankel(running,rows,cols));
  cout << running << endl;
  cout << "Hankel: " << endl;
  cout << m << endl;
  return hankel(running, rows, cols);
}

Vec<hankel> hermite_pade::create_random_hankel_on_type(long rows, const Vec<long> &type){
  Vec<hankel> mh;
  cout << "needed rows: " << rows << endl;
  for (long i = 0; i < type.length(); i++)
    mh.append(create_random_hankel(rows,type[i]+1));
  return mh;
}

// multiplication of left diagonal matrices
Vec<ZZ_p> hermite_pade::mul_bottom_mosaic_diagonal(const Vec<ZZ_p> &v, long rows){
	if (added == 0){
		Vec<ZZ_p> empty;		
		return empty;
	}
  Vec<ZZ_p> diag_ZZ_p1, diag_ZZ_p2;
  Vec<Vec<ZZ_p>> partial_results;
  conv(diag_ZZ_p1, diagonals1);
  conv(diag_ZZ_p2, diagonals2);
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
  	Vec<ZZ_p> x1, x2;
  	x1.SetLength(rows, ZZ_p(0));
  	x2.SetLength(rows, ZZ_p(0));
  	for (long j = 0; j < std::min(rows,type[i]+1); j++)
  	  x1[j] = v[acc+j] * diag_ZZ_p1[i];
  	for (long j = 0; j < std::min(type[i],rows-1); j++)
  	  x2[j+1] = v[acc+j] * diag_ZZ_p2[i];
  	acc += type[i] + 1;
  	Vec<ZZ_p> x;
  	for (long j = 0; j < rows; j++)
  	  x.append(x1[j] + x2[j]);
  	partial_results.append(x);
  }
  Vec<ZZ_p> result;
  for (long i = 0; i < rows; i++){
  	ZZ_p sum;
  	for(long j = 0; j < type.length(); j++){
  		sum += partial_results[j][i];
  	}
  	result.append(sum);
  }
  return result;
}

hermite_pade::hermite_pade(const ZZX &f, const Vec<long>& type, long prec_inp, long fft_index){
  std::srand(std::time(0));
  zz_p::FFTInit(fft_index);
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
  
  long type_sum = 0;
  for (long i = 0; i < type.length(); i++){
    cout << "type: " << type[i] << endl;
    type_sum += 1+type[i];
  }
  
  Mat<zz_p> mat;
    
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
    to_dense(mat,hankel(inp_vec,prec,type[i]+1));
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
  rank = invert(invA,CL); // inverting M mod p
  cout << "original rank: " << rank << endl;

  cout << "type_sum: " << type_sum << endl;
  added = 0;
  if (rank < type_sum-1){
  	long diff = this->prec - rank;
  	added = diff;
  	hankel_matrices.append(create_random_hankel_on_type(diff, type));
  	MH = mosaic_hankel(hankel_matrices);
  	cout << "diags: " << diagonals1 << endl;
  	cout << "diags: " << diagonals2 << endl;
  }

  // setting up the Cauchy matrix
  to_cauchy_grp(CL, X_int, Y_int, e_zz_p, f_zz_p, MH); // converting from Hankel to Cauchy
  rank = invert(invA,CL); // inverting M mod p
  cout << "new rank: " << rank << endl;
  sizeX = X_int.length();
  sizeY = Y_int.length();
  cout << "sizeY: " << sizeY << endl;

  // initializing the bivariate modular comp
  ZZ_pX f_ZZ_pX;
  conv(f_ZZ_pX, f);
  BivariateModularComp M(f_ZZ_pX, type, rank); // could pass in the precomputed stuff

  // initializing the pointer variables and vectors
  vec_M.append(M);

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
  zz_p w_zz_p, w2;
  X_int.point(w_zz_p,1);
  Y_int.point(w2,1);
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
  this->M = &vec_M[0];
  this->X_int = &vec_X_int[0];
  this->Y_int = &vec_Y_int[0];
  this->CL = CL;
  
  //Mat<zz_p> mat;
  to_dense(mat, MH);
  cout << mat << endl;
  //cout << "RANK: " << rank << endl;
}

void hermite_pade::switch_context(long n){
  if (n < vec_M.length()){ // it has already been computed
    ZZ_p::init(p_powers[n]);
    M = &vec_M[n];
    X_int = &vec_X_int[n];
    Y_int = &vec_Y_int[n]; 
  } else{
    // calculating the new power of p
    ZZ p_new(p);
    long pow2 = power_long(2,n);
    power(p_new,p_new, pow2); // 2^n isn't going to be very large
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
    conv(w_p,new_w);
    conv(c_p,c);
    conv(d_p,d);
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

void hermite_pade::dostuff(){
  Vec<ZZ_p> extractor; // mult with A to get a column
  extractor.SetLength(sizeY, ZZ_p(0));
  extractor[rank] = 1; // just for now, take the last column
  for (long i = 0; i < 6; i++){
  	switch_context(i);
  	Vec<ZZ_p> b = conv<Vec<ZZ_p>>(mulA_right(extractor));
  	Vec<ZZ_p> x;
  	DAC(x,b,i);
  	x.append(ZZ_p(-1));
  	x = conv<Vec<ZZ_p>>(find_original_sol(x));
  	cout << x << endl;
  	
  	if(i == 5){
  		Vec<Vec<ZZ>> sol;
  		reconstruct(sol,x,5);
  		cout << sol << endl;
  		sol = flip_on_type(sol);
  		cout << sol << endl;
  	}
  }
  
}

Vec<ZZ> hermite_pade::mulA_right(Vec<ZZ_p> b){
  // Y_int = D_d * Y_int * D_d^(-1)
  // D_e X_int M Y_int^t D_f
  b.SetLength(sizeY, ZZ_p(0)); // padding it
  Vec<ZZ_p> x,e(this->e),f(this->f);
  ZZ_pX temp;
  Vec<ZZ_p> D_d = construct_diagonal(this->d, sizeY);
  Vec<ZZ_p> invD_d = invert_diagonal(D_d);
  mul_diagonal_right(x,f,b);
  mul_diagonal_right(x,invD_d,x);
  conv(temp, x);
  Y_int->evaluate(x,temp);
  mul_diagonal_right(x,D_d,x);
  Vec<ZZ_p> flipped = flip_on_type(x);
  x = M->mult(flipped);
  Vec<ZZ_p> res = mul_bottom_mosaic_diagonal(flipped,added);
  for (long i = 0; i < res.length(); i++)
    x.append(res[i]);
  conv(temp, x);
  X_int->evaluate(x,temp);
  mul_diagonal_right(x,e,x);
  return conv<Vec<ZZ>>(x);
}

Vec<ZZ> hermite_pade::find_original_sol(const Vec<ZZ_p> &b){
  ZZ_pX temp;
  Vec<ZZ_p> x,f(this->f);
  Vec<ZZ_p> D_d = construct_diagonal(this->d, sizeY);
  Vec<ZZ_p> invD_d = invert_diagonal(D_d);
  mul_diagonal_right(x,f,b);
  mul_diagonal_right(x,invD_d,x);
  conv(temp, x);
  Y_int->evaluate(x,temp);
  mul_diagonal_right(x,D_d,x);
  return conv<Vec<ZZ>>(x);
}

void hermite_pade::mul_diagonal_right(Vec<ZZ_p> &x, const Vec<ZZ_p> &b, const Vec<ZZ_p> &a){
  if (b.length() != a.length()) throw "size mismatch";
  x.SetLength(b.length());
  for (long i = 0; i < b.length(); i++)
    x[i] = b[i] * a[i];
}

void hermite_pade::DAC(Vec<ZZ_p> &x, const Vec<ZZ_p>& b_in, long n){
  Vec<ZZ> b_ZZ;
  conv(b_ZZ,b_in);
  switch_context(n);
  Vec<ZZ_p> b;
  conv(b,b_ZZ);
  if (n == 0){ // since p^2^0 == p 
    Vec<zz_p> x_zz_p;
    Vec<zz_p> b_zz_p;
    Vec<ZZ> b_ZZ;
    conv(b_ZZ,b);
    conv(b_zz_p, b_ZZ);
    b_zz_p.SetLength(rank);
    mul_right(x_zz_p, invA, b_zz_p);
    Vec<ZZ> x_ZZ;
    conv(x_ZZ, x_zz_p);
    conv(x, x_ZZ);
    return;
  }
  Vec<ZZ_p> x_0,x_1;
  DAC(x_0, b, n-1); // first recursive call
  switch_context(n);
  Vec<ZZ_p> r; // the error
  r = conv<Vec<ZZ_p>>(mulA_right(x_0));
  r = r - b;
  Vec<ZZ> r_ZZ;
  conv(r_ZZ, r);
  for (long i = 0; i < r.length(); i++){
    r_ZZ[i] = r_ZZ[i] / p_powers[n-1];
    }	
  conv(r,r_ZZ);
  DAC(x_1,r,n-1);
  switch_context(n);
  ZZ_p p_pow;
  conv(p_pow, p_powers[n-1]);
  x = x_0 - p_pow * x_1;
  return;
}

bool hermite_pade::can_reconstruct(const Vec<ZZ_p> &v, long n){
  if (n <= 2) return false;
  for (long i = 0; i < v.length(); i++){
    ZZ a,b;
    try{
    	long result = ReconstructRational(a,b,conv<ZZ>(v[i]),p_powers[n], p_powers[n-1], p_powers[n-1]);
        if (result == 0) return false;
    }catch(...){return false;}
  }
  /*
  Vec<ZZ> v2;
  conv(v2, v);
  long prime = RandomPrime_long(16);
  ZZ_p::init(ZZ(prime));
  Vec<ZZ_p> v3;
  conv(v3,v2);
  ZZ_pX f_p;
  conv(f_p, f_full_prec);
  BivariateModularComp M(f_p, type, rank);
  Vec<ZZ_p> check = M.mult(flip_on_type(v3));
  cout << "diff prime" << check << endl;
  for (long i = 0; i < check.length(); i++)
    if (check[i] != ZZ_p(0)){
      switch_context(n);
      return false;
    }
  switch_context(n);
  */
  return true;
}

void hermite_pade::reconstruct(Vec<Vec<ZZ>> &sol, const Vec<ZZ_p> &v, long n){
  for (long i = 0; i < v.length(); i++){
    ZZ a,b;
    ReconstructRational(a,b,conv<ZZ>(v[i]),p_powers[n], p_powers[n-1], p_powers[n-1]);
    Vec<ZZ> temp;
    temp.append(a);
    temp.append(b);
    sol.append(temp);
  }
}

void hermite_pade::find_rand_sol(Vec<Vec<ZZ>> &sol){
  Vec<ZZ_p> b; // rhs of the equation Ax = b
  Vec<ZZ_p> extractor; // mult with A to get a column
  extractor.SetLength(sizeY, ZZ_p(0));
  extractor[rank] = 1; // just for now, take the last column
  b = conv<Vec<ZZ_p>>(mulA_right(extractor)); // b is the last column of A
  long n = 0; // start at p^2^n
  Vec<ZZ_p> x,x_1,soln;
  DAC(x,b,n); // solution mod p
  Mat<zz_p> mat;
  to_dense(mat,CL);
  
  // padding x
  long x_length = x.length();
  x.SetLength(sizeY,ZZ_p(0));
  x[rank] = -1;
  soln = conv<Vec<ZZ_p>>(find_original_sol(x));
  x.SetLength(x_length);
  cout << "soln: " << soln << endl;
  cout << "init check: " << M->mult(flip_on_type(soln)) << endl;
  
  Vec<ZZ> soln_ZZ,x_ZZ;
  conv(soln_ZZ,soln);
  conv(x_ZZ, x);
  
  // loop until we get enough prec
  while(!can_reconstruct(conv<Vec<ZZ_p>>(soln_ZZ),n)){
    switch_context(++n);
  	Vec<ZZ_p> b; // rhs of the equation Ax = b
  	Vec<ZZ_p> extractor; // mult with A to get a column
  	extractor.SetLength(sizeY, ZZ_p(0));
  	extractor[rank] = 1; // just for now, take the last column
  	b = conv<Vec<ZZ_p>>(mulA_right(extractor)); // b is the last column of A
  	Vec<ZZ_p> x,x_1,soln;
  	conv(x,x_ZZ);
  	x.SetLength(x_length);
    Vec<ZZ_p> r; // the error
    r = conv<Vec<ZZ_p>>(mulA_right(x));
    r = r - b;
    Vec<ZZ> r_ZZ;
    conv(r_ZZ,r);
  	for (long i = 0; i < r.length(); i++){
  		cout << "check R: " << r_ZZ[i] % p_powers[0] << endl;
   		r_ZZ[i] = r_ZZ[i] / p_powers[n-1];
   	}	
    conv(r,r_ZZ);
    DAC(x_1,r,n-1);
    switch_context(n);
    ZZ_p p_pow;
    conv(p_pow, p_powers[n-1]);
    x = x - p_pow * x_1;
    
    // padding x
    x.SetLength(sizeY,ZZ_p(0));
  	x[rank] = -1;
  	//cout << "MUL A: " << mulA_right(x) << endl;
  	soln = conv<Vec<ZZ_p>>(find_original_sol(x));
  	ZZ_p first;
    for (long i = 0; i < soln.length(); i++)
      if (soln[i]._ZZ_p__rep % p_powers[0] != ZZ(0)){
        first = soln[i];
        cout << "i: " << i << endl;
        break;
      }
   // cout << "first: " << first   << endl;
    for (long i = 0; i < soln.length(); i++){
      soln[i] = soln[i] / first;
    }
    cout << "check: " << M->mult(flip_on_type(soln)) << endl;
  	x.SetLength(x_length);
    conv(soln_ZZ, soln);
    conv(x_ZZ,x);
  }
  reconstruct(sol,conv<Vec<ZZ_p>>(soln_ZZ),n);
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
  cout << "soln2: " << soln2 << endl;
  }
}















