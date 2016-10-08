#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){

  for (long i = 1; i < 1000; i+=1){

    Vec<long> src_vec, tmp_vec, dest_vec;
    src_vec.SetLength(2*i);
    tmp_vec.SetLength(2*i);
    dest_vec.SetLength(2*i);

    long *src = src_vec.elts();
    long *tmp = tmp_vec.elts();
    long *dest = dest_vec.elts();

    for (long j = 0; j < i; j++){
      src[j] = random_zz_p().LoopHole();
      dest[j] = src[j];
    }
    for (long j = i; j < 2*i; j++)
      src[j] = 0;

    zz_pX_Multipoint_CTFT c(i);
    c.reduce(src);
    for (long j = 0; j < i; j++)   // warning, src is in [0,2p)
      src[j] = src[j] % zz_p::modulus();
    
    c.CRT(src, tmp);
    
    for (long j = 0; j < i; j++)
      assert (dest[j] == src[j]);
    
    cout << i << endl;
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
  zz_p::FFTInit(0);

  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
