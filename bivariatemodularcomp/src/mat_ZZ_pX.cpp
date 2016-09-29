#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;

void add(Mat<ZZ_pX>& X, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B){
  long n = A.NumRows();
  long m = A.NumCols();

  if (B.NumRows() != n || B.NumCols() != m)
    throw "matrix add: dimension mismatch";

  X.SetDims(n,m);

  long i,j;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= m; j++)
      add(X(i,j),A(i,j),B(i,j));
}

void mul_aux(Mat<ZZ_pX>& X, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B){
  long n = A.NumRows();
  long l = A.NumCols();
  long m = B.NumCols();

  if (l != B.NumRows())
    throw "matrix mul: dimension mismatch";

  X.SetDims(n,m);

  long i, j, k;
  ZZ_pX acc,tmp;

  for (i = 1; i <= n; i++){
    for (j = 1; j <= m; j++){
      clear(acc);
      for (k = 1; k <= l; k++){
        mul(tmp,A(i,k),B(k,j));
	add(acc,acc,tmp);
      }
      X(i,j) = acc;
    }
  }
}

void mul(Mat<ZZ_pX>& X, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B){
  if (&X == &A || &X == &B){
    Mat<ZZ_pX> tmp;
    mul_aux(tmp,A,B);
    X = tmp;
  }
  else mul_aux(X,A,B);
}