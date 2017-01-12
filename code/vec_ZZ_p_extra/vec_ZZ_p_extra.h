#ifndef __VEC_ZZ_P_EXTRA_H
#define __VEC_ZZ_P_EXTRA_H

#include <NTL/vector.h>
#include <NTL/ZZ_p.h>

NTL_CLIENT

/*---------------------------------------------*/
/* random vector of length d                   */
/*---------------------------------------------*/
void random(Vec<ZZ_p>& A, long d);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv_naive(Vec<ZZ_p>& invA, const Vec<ZZ_p>& A);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<ZZ_p>& invA, const Vec<ZZ_p>& A);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<ZZ_p>& A);

/*---------------------------------------------*/
/* pairwise product of the entries of a and b  */
/*---------------------------------------------*/
void mul_diagonal(Vec<ZZ_p> &x, const Vec<ZZ_p> &b, const Vec<ZZ_p> &a);



#endif
