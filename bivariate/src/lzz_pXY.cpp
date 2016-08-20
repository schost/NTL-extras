#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/lzz_pEX.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"
#include "lzz_pXY.h"
#include "ZZXY.h"

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

 
