#include "vec_ZZ_p_extra.h"

NTL_CLIENT

/*---------------------------------------------*/
/* random vector of length d                   */
/*---------------------------------------------*/
void random(Vec<ZZ_p>& A, long d){
  A.SetLength(d);
  for (long i = 0; i < A.length(); i++)
    A[i] = random_ZZ_p();
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv_naive(Vec<ZZ_p>& invA, const Vec<ZZ_p>& A){
  long n = A.length();
  invA.SetLength(n);

  for (long i = 0; i < n; i++)
    invA[i] = 1/A[i];
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<ZZ_p>& invA, const Vec<ZZ_p>& A){
  long n = A.length();
  Vec<ZZ_p> tmp;
  tmp.SetLength(n);
  invA.SetLength(n);

  if (n == 0)
    return;

  tmp[0] = A[0];
  for (long i = 1; i < n; i++)
    tmp[i] = tmp[i-1]*A[i];
  ZZ_p aux = 1/tmp[n-1];
  for (long i = n-1; i >= 1; i--){
    invA[i] = aux*tmp[i-1];
    aux *= A[i];
  }
  invA[0] = aux;
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<ZZ_p>& A){
  Vec<ZZ_p> tmp = A;
  inv(A, tmp);
}

/*---------------------------------------------*/
/* pairwise product of the entries of a and b  */
/*---------------------------------------------*/
void mul_diagonal(Vec<ZZ_p> &x, const Vec<ZZ_p> &b, const Vec<ZZ_p> &a){
  if (b.length() != a.length()) 
    throw "size mismatch";
  x.SetLength(b.length());
  for (long i = 0; i < b.length(); i++)
    x[i] = b[i] * a[i];
}

