#ifndef __LZZ_P_EXTRA__H
#define __LZZ_P_EXTRA__H

NTL_CLIENT

/*-----------------------------------------------------------*/
/* finds q that has order > n                                */
/* q > q0                                                    */
/*-----------------------------------------------------------*/
void find_root(zz_p& q, long n, const zz_p &q0);

/*-----------------------------------------------------------*/
/* finds q that has order > n                                */
/*-----------------------------------------------------------*/
void find_root(zz_p& q, long n);

#endif
