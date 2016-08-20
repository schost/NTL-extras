#include <sstream> 
#include <algorithm> 
#include <NTL/tools.h>
#include <NTL/vector.h>

#include "ZZ_CRT.h"
#include "ratrecon.h"
#include "ZZXY.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"
#include "magma_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over ZZ           */
/* a ZZXY is simply a vector of ZZX                           */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* total degree                                               */
/*------------------------------------------------------------*/
long ZZXY::tdeg() const {
  long res = -1;
  for (long i = 0; i < coeffX.length(); i++)
    res = max(res, deg(coeffX[i]) + i);
  return res;
}

/*------------------------------------------------------------*/
/* degree in Y                                                */
/*------------------------------------------------------------*/
long ZZXY::degY() const {
  return coeffX.length()-1;
}

/*------------------------------------------------------------*/
/* degree in X                                                */
/*------------------------------------------------------------*/
long ZZXY::degX() const {
  long d = -1;
  for (long i = 0; i < coeffX.length(); i++)
    d = max(d, deg(coeffX[i]));
  return d;
}

/*------------------------------------------------------------*/
/* resizes the array of coefficients                          */
/* to remove the trailing entries that are zero, if any       */
/*------------------------------------------------------------*/
void ZZXY::normalize(){
  long idx = coeffX.length()-1;
  while (idx >= 0 && coeffX[idx] == 0)
    idx--;
  coeffX.SetLength(idx+1);
}

/*------------------------------------------------------------*/
/* builds from a file                                         */
/* format: [[a0,a1,..,aN],...,[x0,x1,..,xM]]                  */
/* where [a0,a1,..,aN] = coeff(f,Y^0) \in Z[X], etc.          */
/*------------------------------------------------------------*/
void ZZXY::read_from_file(const string& filename){
  ifstream input;
  input.open(filename);
  string line;
  getline(input, line);
  replace(line.begin(), line.end(), ',', ' ');
  stringstream line_stream;
  line_stream.str(line);

  line_stream >> coeffX;
  input.close();
}


/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
istream& operator>>(istream& s, ZZXY& x){
  NTL_INPUT_CHECK_RET(s, s >> x.coeffX);
  x.normalize();
  return s;
}

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const ZZXY& a){
  return s << a.coeffX;
}

/*------------------------------------------------------------*/
/* initializes M<y,x>=QQ[y,x] with lex order y > x            */
/*------------------------------------------------------------*/
void magma_init_bi_QQ(){
  cout << "M<YYY,XXX> := PolynomialRing(Rationals(),2);\n";
}


/*------------------------------------------------------------*/
/* prints a poly in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZXY & a){
  cout << "(M!(0)";
  long d = a.degY();
  for (long i = 0; i <= d; i++){
    cout << "+(";
    magma_output(a.coeffX[i], "XXX");
    cout << ")*YYY^" << i;
   }
  cout << ")";
}

/*------------------------------------------------------------*/
/* assigns a poly in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZXY & a, const string & name){
  cout << name << " := ";
  magma_output(a);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(ZZXY & dyF, const ZZXY & F){
  long dy = F.degY();
  Vec<ZZX> coeffs_dyF;

  if (dy <= 0){
    coeffs_dyF.SetLength(0);
  }
  else{
    coeffs_dyF.SetLength(dy);
    for (long i = 0; i < dy; i++)
      coeffs_dyF[i] = (i+1)*F.coeffX[i+1];
  }
  dyF = ZZXY(coeffs_dyF);
}



/*------------------------------------------------------------*/
/* prints a pair in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZ_bivariate_regular_chain& T){
  cout << "[";
  magma_output(T.T1, "XXX");
  cout << ", ";
  magma_output(T.T2);
  cout << "]";
}

/*------------------------------------------------------------*/
/* assigns a pair in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_bivariate_regular_chain& T, const string& name){
  cout << name << " := ";
  magma_output(T);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* prints a vector of pairs in M = QQ[y > x]                  */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_bivariate_regular_chain>& T){
  cout << "[";
  for (long i = 0; i < T.length(); i++){
    magma_output(T[i]);
    if (i < T.length()-1)
      cout << ", ";
  }
  cout << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector of pairs in M = QQ[y > x]                 */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_bivariate_regular_chain>& T, const string& name){
  cout << name << " := ";
  magma_output(T);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* utilities for solve                                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* conversion Vec<zz_p> to Vec<long>                          */
/*------------------------------------------------------------*/
static inline
void conv(Vec<long>& along, const Vec<zz_p>& a){
  along.SetLength(a.length());
  for (long i = 0; i < a.length(); i++)
    along[i] = a[i]._zz_p__rep;
}

/*------------------------------------------------------------*/
/* sets the lengths of the entries of array ar                */
/* (we need this several below)                               */
/*------------------------------------------------------------*/
template <class T>
static inline void set_all_lengths(Vec<Vec<Vec<T>>>& ar, long sz, const Vec<long>& deg_X, const Vec<long>& deg_Y){
  ar.SetLength(sz);
  for (long j = 0; j < sz; j++){
    long dX = deg_X[j];
    long dY = deg_Y[j];
    ar[j].SetLength(dY+2);
    ar[j][0].SetLength(dX+1);
    for (long ell = 0; ell <= dY; ell++)
      ar[j][ell+1].SetLength(dX);
  }
}

/*------------------------------------------------------------*/
/* takes all coefficients of sols, appends them to images_T2  */
/*------------------------------------------------------------*/
static inline 
void add_to(Vec<Vec<Vec<Vec<long>>>> & images_T2, const Vec<zz_pEX_augmented>& sols, const Vec<long>& deg_X, const Vec<long>& deg_Y){
  Vec<long> tmp;
  long sz = sols.length();
  for (long j = 0; j < sz; j++){
    long dX = deg_X[j];
    long dY = deg_Y[j];
    conv(tmp, sols[j].T1X.rep);
    for (long k = 0; k <= dX; k++)
      images_T2[j][0][k].append(tmp[k]);
    for (long ell = 0; ell <= dY; ell++){
      conv(tmp, coeff(sols[j].T2, ell)._zz_pE__rep.rep);
      for (long k = 0; k < tmp.length(); k++)
    	images_T2[j][ell+1][k].append(tmp[k]);
      for (long k = tmp.length(); k < dX; k++)
    	images_T2[j][ell+1][k].append(0);
    }
  }
}

/*------------------------------------------------------------*/
/* takes all coefficients of sols, puts them in images_T2     */
/*------------------------------------------------------------*/
static inline 
void set(Vec<Vec<Vec<long>>> & one_T2, const Vec<zz_pEX_augmented>& sols, const Vec<long>& deg_X, const Vec<long>& deg_Y){
  long sz = sols.length();
  for (long j = 0; j < sz; j++){
    long dX = deg_X[j];
    long dY = deg_Y[j];
    conv(one_T2[j][0], sols[j].T1X.rep);
    for (long ell = 0; ell <= dY; ell++){
      conv(one_T2[j][ell+1], coeff(sols[j].T2, ell)._zz_pE__rep.rep);
      long a = one_T2[j][ell+1].length();
      one_T2[j][ell+1].SetLength(dX);
      for (long k = a; k < dX; k++)
    	one_T2[j][ell+1][k] = 0;
    }
  }
}


/*------------------------------------------------------------*/
/* computes a family of regular chains                        */
/* whose solutions is V(FZ, GZ)                               */
/* assumes V(FZ, GZ) finite, not empty                        */
/* otherwise, sols.length is set to zero                      */
/* the decomposition may not be the equiprojectable one       */
/*------------------------------------------------------------*/
void solve(Vec<ZZ_bivariate_regular_chain> & sols, const ZZXY & FZ, const ZZXY & GZ){

  long p = 288230376151711813;

  /*------------------------------------------------------------*/
  /* A first run with 5 primes                                  */
  /* we will keep the ones with most frequent degree pattern    */
  /*------------------------------------------------------------*/
  Vec<long> primes_5;   // our primes
  Vec<Vec<zz_pEX_augmented>> images_sols_5; // the solutions modulo these primes 
  Vec<Vec<long>> list_deg_X; // the X-degrees modulo these primes
  Vec<Vec<long>> list_deg_Y; // the X-degrees modulo these primes
  
  for (long i = 0; i < 5; i++){
    p = NextPrime(p+1);
    
    zz_p::init(p);
    zz_pXY F, G;
    conv(F, FZ);
    conv(G, GZ);

    Vec<zz_pEX_augmented> sols;
    solve(sols, F, G);
    
    long sz = sols.length();
    for (long j = 0; j < sz; j++){ // makes it a rational parametrization
      zz_pEPush push(sols[j].T1); 
      zz_pX dT1 = diff(sols[j].T1X);
      for (long k = 0; k <= deg(sols[j].T2); k++){
	zz_pE tmp;
	conv(tmp, MulMod(dT1, rep(coeff(sols[j].T2, k)), sols[j].T1X));
	SetCoeff(sols[j].T2, k, tmp);
      }
    }

    images_sols_5.append(sols);
    primes_5.append(p);
    Vec<long> deg_X, deg_Y;
    deg_X.SetLength(sz);
    deg_Y.SetLength(sz);
    for (long j = 0; j < sz; j++){
      deg_X[j] = deg(sols[j].T1X);
      deg_Y[j] = deg(sols[j].T2);
    }
    list_deg_X.append(deg_X);
    list_deg_Y.append(deg_Y);
  }
  
  /*------------------------------------------------------------*/
  /* find candidates by majority                                */
  /* nbs[i] = number of primes p s.t. degrees[p] = degrees[i].  */
  /*------------------------------------------------------------*/
  Vec<long> nbs;
  nbs.SetLength(5);
  for (long i = 0; i < 5; i++){
    nbs[i] = 0;
    for (long j = 0; j < 5; j++)
      if ((list_deg_X[i] == list_deg_X[j]) && (list_deg_Y[i] == list_deg_Y[j]))
	  nbs[i]++;
  }
  // finds the max of the nbs
  int mx = 0;
  for (long i = 1; i < 5; i++)
    if (nbs[i] > nbs[mx])
      mx = i;
  // deg_X, deg_Y are the common degree sequences, sz is their length
  Vec<long> deg_X = list_deg_X[mx];
  Vec<long> deg_Y = list_deg_Y[mx];
  long sz = deg_X.length();

  // no solution, or infinitely many solutions
  if (sz == 0){
    sols.SetLength(0);
    return;
  }

  /*------------------------------------------------------------*/
  /* images_T2 stores all the remainders we will do CRT to      */
  /* 1st dimension: selects the trig. set. we work with         */
  /* 2nd dimension: selcts polynomial; 0=T1, 1..dY+1 = T2       */
  /* 3rd dimension: selects coefficient in X                    */
  /* for any such i,j,k, images_T2[i][j][k]=vector of remainders*/
  /*------------------------------------------------------------*/
  Vec<Vec<Vec<Vec<long>>>> images_T2;
  set_all_lengths(images_T2, sz, deg_X, deg_Y);

  // copy the coefficients we know in images_T2
  // todo: convert them to rational representation
  // 
  // primes = array of primes known so far; modulus = &*primes
  Vec<long> primes;
  ZZ modulus = to_ZZ(1);
  for (long i = 0; i < 5; i++)
    if ((deg_X == list_deg_X[i]) && (deg_Y == list_deg_Y[i])){
      add_to(images_T2, images_sols_5[i], deg_X, deg_Y);
      primes.append(primes_5[i]);
      modulus *= primes_5[i];
    }

  /*------------------------------------------------------------*/
  /* a loop to find a witness solution, using a new prime       */
  /* the witness solution is a Vec<Vec<long>>                   */
  /*------------------------------------------------------------*/
  long p_witness;
  Vec<Vec<Vec<long>>> T_witness;
  T_witness.SetLength(sz);
  set_all_lengths(T_witness, sz, deg_X, deg_Y);

  bool found_witness = false;
  while (! found_witness){
    p = NextPrime(p+1);

    zz_p::init(p);
    zz_pXY F, G;
    conv(F, FZ);
    conv(G, GZ);
    Vec<zz_pEX_augmented> sols_witness;
    solve(sols_witness, F, G);

    long sz_local = sols_witness.length();
    if (sz_local != sz)
      continue;
    for (long j = 0; j < sz; j++){
      if (deg(sols_witness[j].T1X) != deg_X[j])
	continue;
      if (deg(sols_witness[j].T2) != deg_Y[j])
	continue;
    }

    for (long j = 0; j < sz; j++){ // makes it a rational parametrization
      zz_pEPush push(sols_witness[j].T1); 
      zz_pX dT1 = diff(sols_witness[j].T1X);
      for (long k = 0; k <= deg(sols_witness[j].T2); k++){
	zz_pE tmp;
	conv(tmp, MulMod(dT1, rep(coeff(sols_witness[j].T2, k)), sols_witness[j].T1X));
	SetCoeff(sols_witness[j].T2, k, tmp);
      }
    }
   
    found_witness = true;
    p_witness = p;
    set(T_witness, sols_witness, deg_X, deg_Y);
  }

  // we will reconstruct T2 by CRT
  Vec<Vec<Vec<ZZ>>> T2;
  T2.SetLength(sz);
  set_all_lengths(T2, sz, deg_X, deg_Y);

  // we will reconstruct T2_rec and T2_den by rational reconstruction
  Vec<Vec<Vec<ZZ>>> T2_rec;
  T2_rec.SetLength(sz);
  Vec<ZZ> T2_den;
  T2_den.SetLength(sz);

  /*------------------------------------------------------------*/
  /* the last loop                                              */
  /* does CRT and rational reconstruction until finished        */
  /*------------------------------------------------------------*/
  bool finished = false;
  long nextnum = ceil( 1.5*(double)primes.length() );
  while (!finished){

    double t;

    t = GetTime();
    p = NextPrime(p+1);
    zz_p::init(p);
    zz_pXY F, G;
    conv(F, FZ);
    conv(G, GZ);
    cerr << "time1: init " << GetTime()-t << endl;

    t = GetTime();
    Vec<zz_pEX_augmented> sols;
    solve(sols, F, G);
    cerr << "time2: solve " << GetTime()-t << endl;

    t = GetTime();
    long sz_local = sols.length();
    if (sz_local != sz)
      continue;
    for (long j = 0; j < sz; j++){
      if (deg(sols[j].T1X) != deg_X[j])
	continue;
      if (deg(sols[j].T2) != deg_Y[j])
	continue;
    }
    for (long j = 0; j < sz; j++){ // makes it a rational parametrization
      zz_pEPush push(sols[j].T1); 
      zz_pX dT1 = diff(sols[j].T1X);
      for (long k = 0; k <= deg(sols[j].T2); k++){
	zz_pE tmp;
	conv(tmp, MulMod(dT1, rep(coeff(sols[j].T2, k)), sols[j].T1X));
	SetCoeff(sols[j].T2, k, tmp);
      }
    }
    primes.append(p);
    modulus *= p;
    add_to(images_T2, sols, deg_X, deg_Y);
    cerr << "time3: convert " << GetTime()-t << endl;

    if (primes.length() < nextnum)
      continue;

    nextnum = ceil( 1.5*(double)primes.length() );

    t = GetTime();
    ZZ_CRT crt(primes);
    for (long j = 0; j < sz; j++){
      long dX = deg_X[j];
      long dY = deg_Y[j];
      for (long k = 0; k <= dX; k++)
    	crt.crt(T2[j][0][k], images_T2[j][0][k]);
      for (long ell = 0; ell <= dY; ell++){
      	for (long k = 0; k < dX; k++)
      	  crt.crt(T2[j][ell+1][k], images_T2[j][ell+1][k]);
      }
    }
    cerr << "time4: crt " << GetTime()-t << endl;

    t = GetTime();
    long s1 = 1;
    for (long j = 0; j < sz; j++)
      s1 = s1 * ReconstructRational(T2_rec[j], T2_den[j], T2[j], modulus, T_witness[j], p_witness);
    if (s1 != 0)
      finished = true;
    cerr << "time5: ratrecon " << GetTime()-t << endl;

  }
  cerr << "-----number of primes: " << primes.length() << endl;
  sols.SetLength(sz);

  /*------------------------------------------------------------*/
  /* assigns the output                                         */
  /*------------------------------------------------------------*/
  for (long j = 0; j < sz; j++){
    long dX = deg_X[j];
    long dY = deg_Y[j];
    sols[j].T1 = 0;
    for (long k = 0; k <= dX; k++)
      SetCoeff(sols[j].T1, k, T2_rec[j][0][k]);
    sols[j].T2.coeffX.SetLength(dY+1);
    for (long ell = 0; ell <= dY; ell++){
      for (long k = 0; k < dX; k++)
	SetCoeff(sols[j].T2.coeffX[ell], k, T2_rec[j][ell+1][k]);
    }
  }
}
