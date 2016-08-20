#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* builds the subproduct tree naively (no Huffman)            */
/*------------------------------------------------------------*/
void build_subproduct_tree(Vec<Vec<zz_pX>> & tree, const Vec<zz_pX> & q){

  long len_tree = 1;
  long len = q.length();
  while (len != 1){
    len_tree++;
    len = (len+1) / 2;
  }

  tree.SetLength(len_tree);
  tree[0] = q;
  
  long idx = 1;
  len = q.length();
  while (len != 1){
    long len_new = (len+1) / 2;
    tree[idx].SetLength(len_new);
    for (long j = 0; j < len / 2; j++)
      tree[idx][j] = tree[idx-1][2*j] * tree[idx-1][2*j+1];
    if (len & 1)
      tree[idx][len_new-1] = tree[idx-1][len-1];
    len = len_new;
    idx++;
  }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* uses subproduct tree.                                      */
/* TODO: thresholds                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* multipoint evaluation                                      */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::evaluate(Vec<zz_p>& val, const zz_pX& f) const {

  // cerr << "general\n";

  long n = this->n;
  val.SetLength(n);
  long d = this->tree.length();
  Vec<zz_pX> remainders;
  remainders.SetLength(n);
  remainders[0] = f % tree[d-1][0];

  for (long i = d-2; i >= 0; i--){
    long len = tree[i].length();
    if (len & 1)
      remainders[len-1] = remainders[(len-1)/2];
    for (long j = len/2-1; j >= 0; j--){
      remainders[2*j+1] = remainders[j] % tree[i][2*j+1];
      remainders[2*j] = remainders[j] % tree[i][2*j];
    }
  }
  for (long i = 0; i < n; i++)
    val[i] = coeff(remainders[i], 0);

}

/*------------------------------------------------------------*/
/* vectorial multipoint evaluation                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const{
  val.SetLength(f.length());
  for (long i = 0; i < f.length(); i++)
    evaluate(val[i], f[i]);  
}

/*------------------------------------------------------------*/
/* multipoint interpolation                                   */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::interpolate(zz_pX& f, const Vec<zz_p>& val) const {
  long n = this->n;
  long d = this->tree.length();
 
 Vec<zz_pX> combinations;
  combinations.SetLength(n);

  for (long i = 0; i < n; i++)
    SetCoeff(combinations[i], 0, val[i]*this->cofactors[i]);
  for (long i = 0; i < d-1; i++){
    long len = tree[i].length();
    for (long j = 0; j < len/2; j++)
      combinations[j] = combinations[2*j]*tree[i][2*j+1] + combinations[2*j+1]*tree[i][2*j];
    if (len & 1)
      combinations[(len-1)/2] = combinations[len-1];
  }

  f = combinations[0];
}

