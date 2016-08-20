#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pXFactoring.h>

#include "magma_output.h"
#include "lzz_pX_CRT.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* squafree decomposition using Yun's algorithm               */
/* todo: check the characteristic                             */
/*------------------------------------------------------------*/
void sqfree_decomposition(Vec<zz_pX>& factors, const zz_pX& f){
  zz_pX df = diff(f);
  zz_pX u = GCD(f, df);
  zz_pX v = f / u;
  zz_pX w = df / u;

  factors.SetLength(0);
  do{
    zz_pX b = w - diff(v);
    zz_pX h = GCD(v, b);
    v = v / h;
    w = b / h;
    if (deg(h) > 0)
      factors.append(h);
  } while (deg(v) > 0);

}

/*------------------------------------------------------------*/
/* helper function for splitting modulus into two factors:    */
/* one "mod_unit" where elt is a unit (we compute its inverse)*/
/* one "mod_zero" where elt is zero                           */
/* assumes modulus to be squarefree                           */
/*------------------------------------------------------------*/
void split(zz_pX & mod_zero, zz_pX & mod_unit, zz_pX & inverse, const zz_pX & elt, const zz_pX & modulus){
  mod_zero = GCD(elt, modulus);
  if (deg(mod_zero) < deg(modulus)){
    mod_unit = modulus / mod_zero;
    MakeMonic(mod_unit);
    inverse = InvMod(elt % mod_unit, mod_unit);
  }
  else{
    mod_unit = 1;
    inverse = 0;
  }
}

/*------------------------------------------------------------*/
/* Let the "degree" of a vector of zz_p be the index of the   */
/* last non-zero entry.                                       */
/* This function splits "modulus" into factors f[i]           */
/* such that for any two roots x,x' of f[i],                  */
/*  deg(elts(x))=deg(elts(x'))=:d[i]                           */
/* and d[i] != d[i'] for i != i'                              */
/* assumes modulus to be squarefree                           */
/*------------------------------------------------------------*/
void split_for_regularize(Vec<zz_pX> & factor_modulus, const Vec<zz_pX> & elts, const zz_pX & modulus){
  factor_modulus.SetLength(deg(modulus));
  zz_pX tmp_modulus = modulus;
  long n = elts.length()-1;
  long i = 0;

  while(n >= 0 && deg(tmp_modulus) > 0){
    zz_pX gcd = GCD(elts[n], tmp_modulus);
    if (deg(gcd) < deg(tmp_modulus)){
      factor_modulus[i] = tmp_modulus / gcd;
      tmp_modulus = gcd;
      i++;
    }
    n--;
  }
  
  if (deg(tmp_modulus) > 0){
    factor_modulus[i] = tmp_modulus;
    i++;
  }

  factor_modulus.SetLength(i);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Stuff for zz_pEX_augmented                                 */
/* a class similar to zz_pEX, with a copy of a zz_pEContext   */
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* outputs as a sequence of M<y,x>                            */
/*------------------------------------------------------------*/
void magma_output(const zz_pEX_augmented & F){
 zz_pEPush push(F.T1); 
  cout << "[";
  magma_output(F.T1X, "X");
  cout << ", ";
  zz_pXY t2;
  conv(t2, F.T2);
  magma_output_bi(t2);
  cout << "]";
}

/*------------------------------------------------------------*/
/* assigns to "name" as a sequence of M<y,x>                  */
/*------------------------------------------------------------*/
void magma_assign(const zz_pEX_augmented & F, const string & name){
  cout << name << " := ";
  magma_output(F);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* random element with deg(T1,x) = dX and deg(T2,y) = dY      */
/*------------------------------------------------------------*/
void random_monic_zz_pEX_augmented(zz_pEX_augmented & res, long dX, long dY){
  zz_pX f = random_zz_pX(dX);
  SetCoeff(f, dX, 1);
  
  Vec<zz_pX> g;

  g.SetLength(dY+1);
  if (dY >= 0){
    for (long j = 0; j < dY; j++)
      g[j] = random_zz_pX(dX);
    do {
      g[dY] = random_zz_pX(dX);
    } while (deg(GCD(g[dY], f)) != 0);
  }
  res = zz_pEX_augmented(f, g);
}

/*------------------------------------------------------------*/
/* replaces T1 by new_mod and reduces T2 by it                */
/*------------------------------------------------------------*/
void mod(zz_pEX_augmented & F_reduced, const zz_pEX_augmented & F, const zz_pX & new_mod){
  if (F.T1X % new_mod != 0){
    cout << "wrong parameters in mod\n";
    exit(-1);
  }

  Vec<zz_pX> F_red_coeffs;
  long d = deg(F.T2);
  F_red_coeffs.SetLength(d+1);
  for (long i = 0; i <= d; i++) 
    F_red_coeffs[i] = rep(coeff(F.T2, i)) % new_mod;

  F_reduced = zz_pEX_augmented(new_mod, F_red_coeffs);
}

/*------------------------------------------------------------*/
/* applies the above with all entries of new_mod              */
/*------------------------------------------------------------*/
void multimod(Vec<zz_pEX_augmented> & reduced, const zz_pEX_augmented & F, const zz_pX_CRT & new_mod){
  if (F.T1X % new_mod.master() != 0){
    cout << "wrong parameters in mod\n";
    exit(-1);
  }

  Vec<Vec<zz_pX>> F_red_coeffs;
  long n = new_mod.length();
  reduced.SetLength(n);
  F_red_coeffs.SetLength(n);

  long d = deg(F.T2);
  for (long k = 0; k < n; k++)
    F_red_coeffs[k].SetLength(d+1);
  for (long i = 0; i <= d; i++) {
    Vec<zz_pX> tmp;
    new_mod.multimod(tmp, rep(coeff(F.T2, i)));
    for (long k = 0; k < n; k++)
      F_red_coeffs[k][i] = tmp[k];
  }
  for (long i = 0; i < n; i++)
    reduced[i] = zz_pEX_augmented(new_mod.moduli(i), F_red_coeffs[i]);
}


/*------------------------------------------------------------*/
/* inverse of the previous one                                */
/*------------------------------------------------------------*/
void combine(zz_pEX_augmented & result, const Vec<zz_pEX_augmented> & rems){

  Vec<zz_pX> mods;
  long n = rems.length();
  mods.SetLength(n);
  for (long i = 0; i < n; i++)
    mods[i] = rems[i].T1X;
  zz_pX_CRT crt(mods);

  Vec<zz_pX> coeffs;
  long dmax = -1;
  for (long k = 0; k < n; k++)
    dmax = max(dmax, deg(rems[k].T2));
  coeffs.SetLength(dmax+1);
  if (dmax != -1){
    Vec<zz_pX> list;
    list.SetLength(n);
    for (long i = 0; i <= dmax; i++){
      for (long k = 0; k < n; k ++)
	list[k] = rep(coeff(rems[k].T2, i));
      crt.combine(coeffs[i], list);
    }
  }
  
  result = zz_pEX_augmented(crt.master(), coeffs);
}

/*------------------------------------------------------------*/
/* splits F into polys where T2 is monic modulo T1            */
/*------------------------------------------------------------*/
void make_monic_and_split(Vec<zz_pEX_augmented> & monic_F, const zz_pEX_augmented & F){
  Vec<zz_pX> split_modulus;
  Vec<zz_pX> coeffs;
  long d = deg(F.T2);
  coeffs.SetLength(d+1);
  for (long i = 0; i <= d; i++)
    coeffs[i] = rep(coeff(F.T2, i));
  split_for_regularize(split_modulus, coeffs, F.T1X);
  zz_pX_CRT crt(split_modulus);
  multimod(monic_F, F, split_modulus);
}

/*------------------------------------------------------------*/
/* computes the squarefree part of f                          */
/* splits if needed, to ensure that the result is monic       */
/*------------------------------------------------------------*/
void squarefree(Vec<zz_pEX_augmented> & sqfree, const zz_pEX_augmented & f){

  pair_zz_pEX_augmented pairs;
  {
    zz_pEPush push(f.T1);
    pairs = pair_zz_pEX_augmented(f.T2, diff(f.T2));
  }

  Vec<zz_pEX_augmented> gcds;
  regular_gcd(gcds, pairs);

  Vec<zz_pX> mods;
  mods.SetLength(gcds.length());
  for (long i = 0; i < gcds.length(); i++)
    mods[i] = gcds[i].T1X;
  zz_pX_CRT crt(mods);
  multimod(sqfree, f, mods);

  for (long i = 0; i < gcds.length(); i++){
    zz_pEPush push(gcds[i].T1);    
    sqfree[i].T2 = sqfree[i].T2 / gcds[i].T2;
    MakeMonic(sqfree[i].T2);
  }
}

/*------------------------------------------------------------*/
/* computes the squarefree part of all entries of f           */
/* splits if needed, to ensure that the result is monic       */
/*------------------------------------------------------------*/
void squarefree(Vec<zz_pEX_augmented> & sqfree, const Vec<zz_pEX_augmented> & f){
  sqfree.SetLength(0);

  for (long i = 0; i < f.length(); i++){
    Vec<zz_pEX_augmented> tmp;
    squarefree(tmp, f[i]);
    sqfree.append(tmp);
  }
}

/*------------------------------------------------------------*/
/* computes a decomposition of the zeros of f, g              */
/* assumes their resultant is non-zero                        */
/*------------------------------------------------------------*/
void solve(Vec<zz_pEX_augmented> & sols, const zz_pXY & f, const zz_pXY & g){
  zz_pX res, sres0, sres1;
  Vec<zz_pEX_augmented> sols_mult;

  double t;
  t = GetTime();
  resultant(res, sres0, sres1, f, g);
  if (res == 0){
    sols.SetLength(0);
    return;
  }
  cerr << "resultant: " << GetTime()-t << ", degree=" << deg(res) << endl;

  t = GetTime();
  Vec<zz_pX> uu;
  sqfree_decomposition(uu, res);
  cerr << "squarefree: " << GetTime()-t << ", size=" << uu.length() << endl;
  
  t = GetTime();
  sols_mult.SetLength(0);
  for (long i = 0; i < uu.length(); i++){
    zz_pX u0, u1, inv;
    split(u0, u1, inv, sres1, uu[i]);
    if (deg(u1) > 0){
      Vec<zz_pX> tmp_inv2;
      tmp_inv2.SetLength(2);
      tmp_inv2[1] = 1;
      tmp_inv2[0] = (inv*sres0) % u1;
      zz_pEX_augmented tmp_inv(u1, tmp_inv2);
      sols_mult.append(tmp_inv);
    }
    if (deg(u0) > 0){
      pair_zz_pEX_augmented fg_EX = pair_zz_pEX_augmented(u0, f, g);
      Vec<zz_pEX_augmented> tmp;
      regular_gcd(tmp, fg_EX);
      sols_mult.append(tmp);
    }
  }
  squarefree(sols, sols_mult);
  cerr << "GCD: " << GetTime()-t << endl;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Stuff for pair_zz_pEX_augmented                            */
/* similar to zz_pEX_augmented, but with two polynomials in Y */
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* outputs as a sequence of M<y,x>                            */
/*------------------------------------------------------------*/
void magma_output(const pair_zz_pEX_augmented & F){
  zz_pEPush push(F.T1); 
  zz_pXY t2;

  cout << "[";
  magma_output(F.T1X, "X");
  cout << ", ";
  conv(t2, F.T20);
  magma_output_bi(t2);
  cout << ", ";
  conv(t2, F.T21);
  magma_output_bi(t2);
  cout << "]";
}

/*------------------------------------------------------------*/
/* assigns to "name" as a sequence of M<y,x>                  */
/*------------------------------------------------------------*/
void magma_assign(const pair_zz_pEX_augmented & F, const string & name){
  cout << name << " := ";
  magma_output(F);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* random, with deg(T1,x)=dX, deg(T2,y)=dY0, deg(T21,y)=dY1   */
/*------------------------------------------------------------*/
void random_pair_monic_zz_pEX_augmented(pair_zz_pEX_augmented & res, long dX, long dY0, long dY1){
  zz_pX f = random_zz_pX(dX);
  SetCoeff(f, dX, 1);
  
  Vec<zz_pX> g0, g1;

  g0.SetLength(dY0+1);
  if (dY0 >= 0){
    for (long j = 0; j < dY0; j++)
      g0[j] = random_zz_pX(dX);
    do {
      g0[dY0] = random_zz_pX(dX);
    } while (deg(GCD(g0[dY0], f)) != 0);
  }

  g1.SetLength(dY1+1);
  if (dY1 >= 0){
    for (long j = 0; j < dY1; j++)
      g1[j] = random_zz_pX(dX);
    do {
      g1[dY1] = random_zz_pX(dX);
    } while (deg(GCD(g1[dY1], f)) != 0);
  }

  res = pair_zz_pEX_augmented(f, g0, g1);
}

/*------------------------------------------------------------*/
/* replaces T1 by new_mod and reduces T20 and T21 by it       */
/*------------------------------------------------------------*/
void mod(pair_zz_pEX_augmented & F_reduced, const pair_zz_pEX_augmented & F, const zz_pX & new_mod){
  if (F.T1X % new_mod != 0){
    cout << "wrong parameters in mod\n";
    exit(-1);
  }

  Vec<zz_pX> F_red_coeffs0, F_red_coeffs1;

  long d0 = deg(F.T20);
  F_red_coeffs0.SetLength(d0+1);
  for (long i = 0; i <= d0; i++) 
    F_red_coeffs0[i] = rep(coeff(F.T20, i)) % new_mod;

  long d1 = deg(F.T21);
  F_red_coeffs1.SetLength(d1+1);
  for (long i = 0; i <= d1; i++) 
    F_red_coeffs1[i] = rep(coeff(F.T21, i)) % new_mod;

  F_reduced = pair_zz_pEX_augmented(new_mod, F_red_coeffs0, F_red_coeffs1);
}

/*------------------------------------------------------------*/
/* applies the above with all entries of new_mod              */
/*------------------------------------------------------------*/
void multimod(Vec<pair_zz_pEX_augmented> & reduced, const pair_zz_pEX_augmented & F, const zz_pX_CRT & new_mod){
  if (F.T1X % new_mod.master() != 0){
    cout << "wrong parameters in mod\n";
    exit(-1);
  }

  Vec<Vec<zz_pX>> F_red_coeffs0, F_red_coeffs1;

  long n = new_mod.length();
  reduced.SetLength(n);
  F_red_coeffs0.SetLength(n);
  F_red_coeffs1.SetLength(n);

  long d0 = deg(F.T20);
  for (long k = 0; k < n; k++)
    F_red_coeffs0[k].SetLength(d0+1);
  for (long i = 0; i <= d0; i++) {
    Vec<zz_pX> tmp;
    new_mod.multimod(tmp, rep(coeff(F.T20, i)));
    for (long k = 0; k < n; k++)
      F_red_coeffs0[k][i] = tmp[k];
  }

  long d1 = deg(F.T21);
  for (long k = 0; k < n; k++)
    F_red_coeffs1[k].SetLength(d1+1);
  for (long i = 0; i <= d1; i++) {
    Vec<zz_pX> tmp;
    new_mod.multimod(tmp, rep(coeff(F.T21, i)));
    for (long k = 0; k < n; k++)
      F_red_coeffs1[k][i] = tmp[k];
  }

  for (long i = 0; i < n; i++)
    reduced[i] = pair_zz_pEX_augmented(new_mod.moduli(i), F_red_coeffs0[i], F_red_coeffs1[i]);
}

/*------------------------------------------------------------*/
/* inverse of the previous one                                */
/*------------------------------------------------------------*/
void combine(pair_zz_pEX_augmented & result, const Vec<pair_zz_pEX_augmented> & rems){

  Vec<zz_pX> mods;
  long n = rems.length();
  mods.SetLength(n);
  for (long i = 0; i < n; i++)
    mods[i] = rems[i].T1X;
  zz_pX_CRT crt(mods);

  Vec<zz_pX> coeffs0;
  long dmax0 = -1;
  for (long k = 0; k < n; k++)
    dmax0 = max(dmax0, deg(rems[k].T20));
  coeffs0.SetLength(dmax0+1);
  if (dmax0 != -1){
    Vec<zz_pX> list;
    list.SetLength(n);
    for (long i = 0; i <= dmax0; i++){
      for (long k = 0; k < n; k ++)
	list[k] = rep(coeff(rems[k].T20, i));
      crt.combine(coeffs0[i], list);
    }
  }

  Vec<zz_pX> coeffs1;
  long dmax1 = -1;
  for (long k = 0; k < n; k++)
    dmax1 = max(dmax1, deg(rems[k].T21));
  coeffs1.SetLength(dmax1+1);
  if (dmax1 != -1){
    Vec<zz_pX> list;
    list.SetLength(n);
    for (long i = 0; i <= dmax1; i++){
      for (long k = 0; k < n; k ++)
	list[k] = rep(coeff(rems[k].T21, i));
      crt.combine(coeffs1[i], list);
    }
  }
  
  result = pair_zz_pEX_augmented(crt.master(), coeffs0, coeffs1);
}

/*------------------------------------------------------------*/
/* splits F into polys where T2_idx is monic modulo T1        */
/*------------------------------------------------------------*/
void make_monic_and_split(Vec<pair_zz_pEX_augmented> & monic_F, const pair_zz_pEX_augmented & F, long idx){
  Vec<zz_pX> split_modulus;
  Vec<zz_pX> coeffs;

  if (idx == 0){
    long d = deg(F.T20);
    coeffs.SetLength(d+1);
    for (long i = 0; i <= d; i++)
      coeffs[i] = rep(coeff(F.T20, i));
  }
  else if (idx == 1){
    long d = deg(F.T21);
    coeffs.SetLength(d+1);
    for (long i = 0; i <= d; i++)
      coeffs[i] = rep(coeff(F.T21, i));
  }
  else{
    cout << "wrong parameter to make_monic_and_split\n";
    exit(-1);
  }

  split_for_regularize(split_modulus, coeffs, F.T1X);
  zz_pX_CRT crt(split_modulus);
  multimod(monic_F, F, split_modulus);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class of stacks of pairs (copied from NTL)        */
/* used to hold the remaining tasks in regular GCD            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
struct _zz_pEX_local_stack {
  long top;
  Vec<pair_zz_pEX_augmented> data;

  _zz_pEX_local_stack() { 
    top = -1; 
  }

  void pop(pair_zz_pEX_augmented & a) { 
    a = data[top];
    top--;
  }

  long size() const {
    return top+1;
  }

  long empty() const { 
    return (top == -1); 
  }

  void push(const pair_zz_pEX_augmented & a){
    if (top+1 >= data.length()) 
      data.SetLength(max(32, long(1.414*data.length())));
    top++;
    data[top] = a;
  }

  void push(const zz_pEContext & ctx, const zz_pEX & a, const zz_pEX & b){
    if (top+1 >= data.length()) {
      data.SetLength(max(32, long(1.414*data.length())));
      top++;
    }
    zz_pEPush push(ctx);
    data[top] = pair_zz_pEX_augmented(a, b);
  }
};

/*------------------------------------------------------------*/
/* computes a regular gcd, splitting if needed                */
/*------------------------------------------------------------*/
void regular_gcd(Vec<zz_pEX_augmented> & h, const pair_zz_pEX_augmented & a_b){
  _zz_pEX_local_stack stack;

  Vec<pair_zz_pEX_augmented> a_monic_b;
  make_monic_and_split(a_monic_b, a_b, 0); // splits along the 1st entry

  for (long i = 0; i < a_monic_b.length(); i++){
    Vec<pair_zz_pEX_augmented> a_monic_b_monic;
    make_monic_and_split(a_monic_b_monic, a_monic_b[i], 1);
    for (long j = 0; j < a_monic_b_monic.length(); j++)
      stack.push(a_monic_b_monic[j]);
  }

  long nb = 0;
  h.SetLength(0);
  
  while (! stack.empty()){
    pair_zz_pEX_augmented current;
    stack.pop(current);
    zz_pEPush push(current.T1);

    long d0 = deg(current.T20);
    long d1 = deg(current.T21);
    
    if (d0 == -1){
      if (nb == h.length())
  	h.SetLength(max(32, long(1.414*h.length())));
      h[nb] = zz_pEX_augmented(current.T21);
      nb++;
    }
    else if (d1 == -1){
      if (nb == h.length())
  	h.SetLength(max(32, long(1.414*h.length())));
      h[nb] = zz_pEX_augmented(current.T20);
      nb++;
    }
    else if (d0 < d1){
      Vec<pair_zz_pEX_augmented> next;
      make_monic_and_split(next, pair_zz_pEX_augmented(current.T21 % current.T20, current.T20), 0);
      for (long i = 0; i < next.length(); i++)
  	stack.push(next[i]);
    }
    else{
      Vec<pair_zz_pEX_augmented> next;
      make_monic_and_split(next, pair_zz_pEX_augmented(current.T20 % current.T21, current.T21), 0);
      for (long i = 0; i < next.length(); i++)
  	stack.push(next[i]);
    }
  }

  h.SetLength(nb);

  for (long i = 0; i < nb; i++){
    zz_pEPush push(h[i].T1);
    MakeMonic(h[i].T2);
  }
}



