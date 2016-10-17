#ifndef __LZZ_P_EXTRA__H
#define __LZZ_P_EXTRA__H

NTL_CLIENT

/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a = 0                                                */
/*------------------------------------------------------------*/
long order(const zz_p& a);

/*------------------------------------------------------------*/
/* multiplicative order of w                                  */
/* assumes it is a power of 2                                 */
/*------------------------------------------------------------*/
long order_dyadic(const zz_p& w);

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord);

#endif
