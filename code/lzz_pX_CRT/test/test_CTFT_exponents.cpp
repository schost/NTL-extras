#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes the exponents of a few integers                   */
/*------------------------------------------------------------*/
void check(int opt){

  for (long i = 0; i < 20; i+=1){
    cout << i << ":   ";
    Vec<long> degrees;
    CTFT_exponents(degrees, i);
    long j = 0;
    while (degrees[j] != -1){
      cout << (1 << degrees[j]) << " ";
      j++;
    }
    cout << endl;
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
