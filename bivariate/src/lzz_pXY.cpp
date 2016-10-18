#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/lzz_pEX.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"
#include "lzz_pXY.h"
#include "ZZXY.h"
#include "lzz_pX_extra.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over zz_p         */
/* a zz_pXY is simply a vector of zz_pX                       */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* total degree                                               */
/*------------------------------------------------------------*/
long zz_pXY::tdeg() const {
  long res = -1;
  for (long i = 0; i < coeffX.length(); i++)
    res = max(res, deg(coeffX[i]) + i);
  return res;
}

/*------------------------------------------------------------*/
/* degree in Y                                                */
/*------------------------------------------------------------*/
long zz_pXY::degY() const {
  return coeffX.length()-1;
}

/*------------------------------------------------------------*/
/* degree in X                                                */
/*------------------------------------------------------------*/
long zz_pXY::degX() const {
  long d = -1;
  for (long i = 0; i < coeffX.length(); i++)
    d = max(d, deg(coeffX[i]));
  return d;
}

/*------------------------------------------------------------*/
/* resizes the array of coefficients                          */
/* to remove the trailing entries that are zero, if any       */
/*------------------------------------------------------------*/
void zz_pXY::normalize(){
  long idx = coeffX.length()-1;
  while (idx >= 0 && coeffX[idx] == 0)
    idx--;
  coeffX.SetLength(idx+1);
}

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const zz_pXY& a){
  return s << a.coeffX;
}

/*------------------------------------------------------------*/
/* initializes V<y>=GF(p)[x][y]                               */
/*------------------------------------------------------------*/
void magma_init_Y(){
  cout << "V<y> := PolynomialRing(U);\n";
}

/*------------------------------------------------------------*/
/* initializes M<y,x>=GF(p)[y,x] with lex order y > x         */
/*------------------------------------------------------------*/
void magma_init_bi(){
  cout << "M<Y,X> := PolynomialRing(k,2);\n";
}

/*------------------------------------------------------------*/
/* prints a poly in V = GF(p)[x][y]                           */
/*------------------------------------------------------------*/
void magma_output(const zz_pXY & a){
  cout << "V!([";
  long d = a.degY();
  if (d == -1){
    cout << "])";
    return;
  }
  for (long i = 0; i < d; i++){
    magma_output(a.coeffX[i]);
    cout << ", ";
  }
  magma_output(a.coeffX[d]);
  cout << "])";
}

/*------------------------------------------------------------*/
/* assigns a poly in V = GF(p)[x][y]                          */
/*------------------------------------------------------------*/
void magma_assign(const zz_pXY & a, const string & n){
  cout << n << " := ";
  magma_output(a);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* prints a poly in M = GF(p)[y > x]                          */
/*------------------------------------------------------------*/
void magma_output_bi(const zz_pXY & a){
  cout << "(M!(0)";
  long d = a.degY();
  for (long i = 0; i <= d; i++){
    cout << "+(";
    magma_output(a.coeffX[i], "X");
    cout << ")*Y^" << i;
   }
  cout << ")";
}

/*------------------------------------------------------------*/
/* assigns a poly in M = GF(p)[y > x]                         */
/*------------------------------------------------------------*/
void magma_assign_bi(const zz_pXY & a, const string & name){
  cout << name << " := ";
  magma_output_bi(a);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* random element f with deg(f,x) < dx and deg(f,y) < dy      */
/*------------------------------------------------------------*/
void random_zz_pXY(zz_pXY & f, long dx, long dy){
  f.coeffX.SetLength(dy+1);
  for (long i = 0; i <= dy; i++)
    f.coeffX[i] = random_zz_pX(dx);
}

/*------------------------------------------------------------*/
/* evaluation with respect to X at a point                    */
/*------------------------------------------------------------*/
void evaluate(zz_pX & value, const zz_pXY & f, const zz_p & point){
  value = to_zz_pX(0);
  for (long i = 0; i <= f.degY(); i++)
    SetCoeff(value, i, eval(f.coeffX[i], point));
  value.normalize();
}

/*------------------------------------------------------------*/
/* multipoint evaluation with respect to X                    */
/* result is a vector of polynomials                          */
/*------------------------------------------------------------*/
void evaluate(Vec<zz_pX> & values, const zz_pXY & f, const zz_pX_Multipoint & ev){
  Vec<Vec<zz_p>> val;
  ev.evaluate(val, f.coeffX);

  long d = f.degY();
  long n = ev.length();
  values.SetLength(n);
  for (long j = 0; j < n; j++){
    values[j].rep.SetLength(d+1);
    zz_p *cf = values[j].rep.elts();
    for (long i = 0; i <= d; i++)
      cf[i] = val[i][j];
    values[j].normalize();
  }
}


/*------------------------------------------------------------------------*/
/* returns the resultant and the last non-trivial sub-resultant           */
/*   rres     is the resultant                                            */
/*   subres1  is the last subresultant (of degree <= 1)                   */
/*   subres2  is filled if the subres1 is degenerated (of degree 0). Then */
/*            subres2 is the last non-trivial subresultant.               */
/* All of this might fail, in the sense that subres1 and subres2 could    */
/* be 0 instead of having meaningful values. In particular, if rres = 0,  */
/* no further information is obtained. However, if a subres is <> 0, then */
/* it is exact.                                                           */
/* Warning: subres1 and subres2 are not allowed to be aliases of the      */
/* input parameters.                                                      */
/*------------------------------------------------------------------------*/
void resultant_subresultant(zz_p& rres, zz_pX& subres1, zz_pX& subres2, const zz_pX& a, const zz_pX& b){
  zz_p res, saved_res;
   
  clear(subres1);
  clear(subres2);
 
  if (IsZero(a) || IsZero(b))
    clear(res);
  else if (deg(a) == 0 && deg(b) == 0) 
    set(res);
  else {
    long d0, d1, d2;
    long sg = 1;
    zz_p lc;
    set(res);

    long n = max(deg(a),deg(b)) + 1;
    zz_pX u(INIT_SIZE, n), v(INIT_SIZE, n), saved_u(INIT_SIZE, n);

    u = a;
    v = b;

    for (;;) {
      d0 = deg(u);
      d1 = deg(v);
      lc = LeadCoeff(v);

      saved_u = u;
      PlainRem(u, u, v);
      swap(u, v);
 
      d2 = deg(v);
      if (d2 >= 0) {
	power(lc, lc, d0-d2);
	saved_res = res;
	mul(res, res, lc);
	if (d0 & d1 & 1) 
	  NTL::negate(res, res);
	if (!((d0-d1) & 1)) 
	  sg = -sg;
      }
      else { // d2 < 0 so v = 0
	if (d1 == 0) { // previous one = non-zero constant
	  subres1 = sg*saved_res*saved_u;
	  if (deg(subres1) > 1) { // degenerate case
	    subres2 = subres1;
	    if (deg(subres1) > 2)
	      clear(subres1);
	    else
	      subres1 = saved_res*lc*power(LeadCoeff(saved_u),d0);
	  }
	  power(lc, lc, d0);
	  mul(res, res, lc);
	}
	else { // previous one = positive degree
	  subres1 = sg*res*u;
	  if (deg(subres1) > 1) { // *very* degenerate case, not reliable
	    subres2 = subres1;
	    if (deg(subres1) > 2)
	      clear(subres1);
	    else
	      subres1 = res*lc*power(LeadCoeff(u),d0);
	  }
	  clear(res);
	}
	break;
      }
    }

    rres = res;
  }
}



/*------------------------------------------------------------*/
/* computes the resultant of f and g                          */
/* together with the subresultant of degree 1, if possible    */
/*------------------------------------------------------------*/
void resultant(zz_pX & res, zz_pX& subres0, zz_pX& subres1, const zz_pXY& f, const zz_pXY& g){
  long dyf = f.degY();
  long dyg = g.degY();

  long df = f.tdeg();
  long dg = g.tdeg();

  long xtra = 10;
  long dr = df*dg + xtra;

  Vec<zz_p> values;
  values.SetLength(dr);

  zz_pX_Multipoint * ev;
  zz_pX_Multipoint_General evQ;
  zz_pX_Multipoint_Geometric evG;


  zz_p a0 = random_zz_p();
  zz_p a = a0*a0;
  zz_p pow_a = to_zz_p(1);
  long OK = 1;

  if (eval(f.coeffX[dyf], pow_a) == 0 || eval(g.coeffX[dyg], pow_a) == 0)
    OK = 0;
  for (long i = 0; i < dr; i++){
    pow_a *= a;
    if (pow_a == to_zz_p(1))
      OK = 0;
    if (eval(f.coeffX[dyf], pow_a) == 0 || eval(g.coeffX[dyg], pow_a) == 0)
      OK = 0;
  }

  if (OK == 1){
    evG = zz_pX_Multipoint_Geometric(a0, df+1, dr, 1);
    ev = &evG;
  }
  else{
    Vec<zz_p> points;
    points.SetLength(dr);
    for (long i = 0; i < dr; i++){
      zz_p pt;
      do{
  	pt = random_zz_p();
      } while (eval(f.coeffX[dyf], pt) == 0 || eval(g.coeffX[dyg], pt) == 0);
      points[i] = pt;
    }
    evQ = zz_pX_Multipoint_General(points);
    ev = &evQ;
  }

  cerr << ((OK == 1) ? "geom " : "general ") << "in length " << ev->length() << endl;

  double tint0, tint1, tres, t;

  t = GetTime();
  Vec<zz_pX> valuesF, valuesG;
  evaluate(valuesF, f, *ev);
  evaluate(valuesG, g, *ev);
  tint0 = GetTime()-t;

  Vec<zz_p> vsubres0, vsubres1;
  vsubres0.SetLength(dr);
  vsubres1.SetLength(dr);

  t = GetTime();
  long OK_res = 1;
  for (long i = 0; i < dr; i++){
    zz_pX tmp1, tmp2;
    resultant_subresultant(values[i], tmp1, tmp2, valuesF[i], valuesG[i]);
    if (values[i] == 0 || deg(tmp1) != 1)
      OK_res = 0;
    else{
      vsubres0[i] = coeff(tmp1, 0);
      vsubres1[i] = coeff(tmp1, 1);
    }
  }
  tres = GetTime()-t;

  t = GetTime();
  ev->interpolate(res, values);
  if (OK_res == 1){
    ev->interpolate(subres0, vsubres0);
    ev->interpolate(subres1, vsubres1);
  }
  else{
    subres0 = 0;
    subres1 = 0;
  }

  tint1 = GetTime()-t;

  cerr << "evaluation: " << tint0 << ", resultant: " << tres << ", interpolation: " << tint1 << endl;
}

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(zz_pXY & dyF, const zz_pXY & F){
  long dy = F.degY();
  Vec<zz_pX> coeffs_dyF;

  if (dy <= 0){
    coeffs_dyF.SetLength(0);
  }
  else{
    coeffs_dyF.SetLength(dy);
    for (long i = 0; i < dy; i++)
      coeffs_dyF[i] = (i+1)*F.coeffX[i+1];
  }
  dyF = zz_pXY(coeffs_dyF);
}

/*------------------------------------------------------------*/
/* conversion zz_pXY to zz_pEX                                */
/*------------------------------------------------------------*/
void conv(zz_pEX & fEX, const zz_pXY & f){
  Vec<zz_pE> coeffs;
  conv(coeffs, f.coeffX);
  conv(fEX, coeffs);
}

/*------------------------------------------------------------*/
/* conversion zz_pEX to zz_pXY                                */
/*------------------------------------------------------------*/
void conv(zz_pXY & f, const zz_pEX & fEX){
  Vec<zz_pE> coeffsE;
  Vec<zz_pX> coeffs;
  conv(coeffsE, fEX);
  conv(coeffs, coeffsE);
  f = zz_pXY(coeffs);
}

/*------------------------------------------------------------*/
/* conversion zzXY to zz_pXY                                  */
/*------------------------------------------------------------*/
void conv(zz_pXY & f, const ZZXY & fZ){
  Vec<zz_pX> coeffs;
  conv(coeffs, fZ.coeffX);
  f = zz_pXY(coeffs);
}

/*------------------------------------------------------------*/
/* reduces P(x,y) modulo (R(x), y-S(x))                       */
/*------------------------------------------------------------*/
void reduce_naive(zz_pX & rem, const zz_pX& R, const zz_pX& S, const zz_pXY & P){
  zz_pXModulus R_mod(R);
  zz_pXMultiplier S_mul(S, R_mod);

  rem = 0;
  long n = P.degY();
  for (long i = n; i >= 0; i--)
    rem = MulMod(rem, S_mul, R_mod) + P.coeffX[i];
}

/*------------------------------------------------------------*/
/* reduces P(x,y) modulo (R(x), y-S(x))                       */
/*------------------------------------------------------------*/
void reduce(zz_pX & rem, const zz_pX& R, const zz_pX& S, const zz_pXY & P){
  zz_pXModulus R_mod(R);
  zz_pXMultiplier S_mul(S, R_mod);

  long d = deg(R);
  long m = P.degX();
  long n = P.degY();

  long n1 = SqrRoot(n);
  long n2 = ceil( ((double)(n+1)) / ((double)n1) );
  long c = ceil( ((double)(d+1)) / ((double)(m+1)) );

  Mat<zz_pX> coeff_P, powers_S;

  double t0 = 0, t1 = 0, t;

  t = GetTime();
  powers_S.SetDims(n1, c);
  zz_pX pS;
  set(pS);
  for (long i = 0; i < n1; i++){
    for (long j = 0; j < c; j++){
      powers_S[i][j].rep.SetLength(m+1);
      long shift = j*(m+1);
      for (long k = 0; k <= m; k++)
        powers_S[i][j].rep[k] = coeff(pS, shift+k);
      powers_S[i][j].normalize();
    }
    MulMod(pS, pS, S_mul, R_mod);
  }
  
  coeff_P.SetDims(n2, n1);
  for (long i = 0; i < n2; i++){
    for (long j = 0; j < n1; j++){
      long idx = i*n1+j;
      if (idx <= n)
	coeff_P[i][j] = P.coeffX[idx];
    }
  }
  t0 += GetTime()-t;

  t = GetTime();
  Mat<zz_pX> res;
  multiply_waksman(res, coeff_P, powers_S);
  t1 += GetTime()-t;

  t = GetTime();
  clear(rem);
  zz_pXMultiplier pS_mul(pS, R_mod);
  for (long i = n2-1; i >= 0; i--){
    zz_pX tmp;
    tmp.rep.SetLength(d+2*(m+1));
    for (long j = 0; j < c; j++){
      long ell = res[i][j].rep.length();
      for (long k = 0; k < ell; k++)
        tmp.rep[j*(m+1)+k] += res[i][j].rep[k];
    }
    tmp.normalize();
    tmp %= R;
    MulMod(rem, rem, pS_mul, R_mod);
    rem += tmp;
  }
  t0 += GetTime()-t;
}

/*------------------------------------------------------------*/
/* naive multiplication                                       */
/*------------------------------------------------------------*/
void mul_naive(zz_pXY& c, const zz_pXY& a, const zz_pXY& b){
  zz_pXY c_tmp;
  long da = a.degY();
  long db = b.degY();

  c_tmp.coeffX.SetLength(da + db + 1, to_zz_pX(0));
  for (long i = 0; i <= da; i++)
    for (long j = 0; j <= db; j++)
      c_tmp.coeffX[i+j] += a.coeffX[i]*b.coeffX[j];
  c_tmp.normalize();
  c = c_tmp;
}


/*------------------------------------------------------------*/
/* contents, as a polynomial in x                             */
/*------------------------------------------------------------*/
void contents(zz_pX& c, const zz_pXY& a){
  long da = a.degY();
  if (da < 0){
    c = 0;
    return;
  }
  
  c = a.coeffX[0];
  for (long i = 1; i <= da; i++)
    c = GCD(c, a.coeffX[i]);
}


/*------------------------------------------------------------*/
/* computes b = a(x+c, y)                                     */
/*------------------------------------------------------------*/
void shift_X(zz_pXY& b, const zz_pXY& a, const zz_p& c){
  long dX = a.degX(), dY = a.degY();
  
  b.coeffX.SetLength(dY+1);
  zz_pX_shift s(c, dX);
  for (long i = 0; i <= dY; i++)
    s.shift(b.coeffX[i], a.coeffX[i]);

}

/*------------------------------------------------------------*/
/* computes c = to_zz_XY(c0*lc), makes it primitive           */
/*------------------------------------------------------------*/
static void reconstruct(zz_pXY & c, const zz_pEX& c0, const zz_pX & lc){
  Vec<zz_pE> coeffs;
  long d = deg(c0);
  coeffs.SetLength(d+1);
  zz_pE lc_E = to_zz_pE(lc);
  for (long i = 0; i <= d; i++)
    coeffs[i] = coeff(c0, i) * lc_E;
  zz_pEX pol_coeffs;
  conv(pol_coeffs, coeffs);
  conv(c, pol_coeffs);

  zz_pX cts;
  contents(cts, c);
  for (long i = 0; i <= d; i++)
    c.coeffX[i] /= cts;
}

/*------------------------------------------------------------*/
/* tests whether c(d, Y) divides a(d, Y) for a random d       */
/*------------------------------------------------------------*/
static bool random_trial_division(const zz_pXY& c, const zz_pXY& a){
  zz_p d = random_zz_p();
  zz_pX cd, ad;
  evaluate(cd, c, d);
  evaluate(ad, a, d);
  return ((ad % cd) == 0);
}

/*------------------------------------------------------------*/
/* GCD                                                        */
/* assumes that the GCD has multiplicity 1 in a               */
/*------------------------------------------------------------*/
void GCD(zz_pXY& c, const zz_pXY& a, const zz_pXY& b){
  long da = a.degY(), db = b.degY();
  if (da < 0){
    c = b;
    return;
  }
  if (db < 0){
    c = a;
    return;
  }

  zz_p e = random_zz_p();
  zz_pXY a_shifted, b_shifted, c_shifted;
  shift_X(a_shifted, a, e);
  shift_X(b_shifted, b, e);
  zz_pX lc_a_shifted = a_shifted.coeffX[da];
  
  zz_pX a0, b0, diff_a0, inv_diff_a0, c0;
  evaluate(a0, a_shifted, to_zz_p(0));
  evaluate(b0, b_shifted, to_zz_p(0));
  diff_a0 = diff(a0);
  c0 = GCD(a0, b0);
  inv_diff_a0 = InvMod(diff_a0, c0);

  zz_pEContext push;
  long prec = 1;
  zz_pX powX = to_zz_pX(1);
  powX = LeftShift(powX, 1);
  zz_pE::init(powX);

  Vec<zz_p> coeff_c0 = conv<Vec<zz_p>> (c0);
  Vec<zz_pX> coeff_c0_X = conv<Vec<zz_pX>> (coeff_c0);
  zz_pEX c0_E;
  Vec<zz_pE> coeff_c0_E = conv<Vec<zz_pE>> (coeff_c0_X);
  conv(c0_E, coeff_c0_E);
  
  Vec<zz_p> coeff_inv = conv<Vec<zz_p>> (inv_diff_a0);
  Vec<zz_pX> coeff_inv_X = conv<Vec<zz_pX>> (coeff_inv);
  zz_pEX inv_E;
  Vec<zz_pE> coeff_inv_E = conv<Vec<zz_pE>> (coeff_inv_X);
  conv(inv_E, coeff_inv_E);

  reconstruct(c_shifted, c0_E, lc_a_shifted);

  while (! random_trial_division(c_shifted, a_shifted)){
    powX = LeftShift(powX, prec);
    prec = 2 * prec;

    Vec<zz_pX> c0_X = conv<Vec<zz_pX>> (c0_E.rep);
    Vec<zz_pX> inv_X = conv<Vec<zz_pX>> (inv_E.rep);

    zz_pE::init(powX);

    Vec<zz_pE> coeff_c0_E = conv<Vec<zz_pE>> (c0_X);
    conv(c0_E, coeff_c0_E);
    Vec<zz_pE> coeff_inv_E = conv<Vec<zz_pE>> (inv_X);
    conv(inv_E, coeff_inv_E);

    zz_pEX a_E, diff_E;
    conv(a_E, a_shifted);
    diff_E = diff(a_E);
    a_E = a_E % c0_E;
    diff_E = diff_E % c0_E;
    inv_E = 2*inv_E - MulMod(MulMod(diff_E, inv_E, c0_E), inv_E, c0_E);
    c0_E = c0_E + MulMod(MulMod(diff(c0_E) % c0_E, inv_E, c0_E), a_E, c0_E);

    reconstruct(c_shifted, c0_E, lc_a_shifted);

  }

  shift_X(c, c_shifted, -e);
  
  zz_pX contents_a, contents_b;
  contents(contents_a, a);
  contents(contents_b, b);
  zz_pX g = GCD(contents_a, contents_b);

  for (long i = 0; i <= c.degY(); i++)
    c.coeffX[i] *= g;
}

/*------------------------------------------------------------*/
/* GCD                                                        */
/* assumes that the GCD has multiplicity 1 in a               */
/*------------------------------------------------------------*/
void GCD(Vec<zz_pX> & c, const Vec<zz_pX>& a, const Vec<zz_pX>& b){
  zz_pXY aY(a);
  zz_pXY bY(b);
  zz_pXY cY;
  GCD(cY, aY, bY);
  c = cY.coeffX;
}
