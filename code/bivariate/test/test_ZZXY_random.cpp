#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "ZZXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* reads from a file and outputs the result                   */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  ZZXY Fpoly;
  random(Fpoly, 10, 10, 10);
  cout << Fpoly << endl;

  ZZX G;
  random(G, 10, 10);
  cout << G << endl;
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
