#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  for (long i = 1; i < 100; i+=1){
    Mat<zz_p> val;
    val.SetDims(i, i);
    zz_pX f, x;
    SetCoeff(f, 0, 1);
    SetCoeff(x, 1, 1);
    zz_pX_Multipoint_CTFT c(i);

    for (long j = 0; j < i; j++){
      Vec<zz_p> val_i;      
      c.evaluate(val_i, f);
      for (long k = 0; k < i; k++){
	val[k][j] = val_i[k];
      }
      f *= x;
    }

    Mat<zz_p> val_t;
    val_t.SetDims(i, i);
    for (long j = 0; j < i; j++){
      Vec<zz_p> in, out;
      in.SetLength(i);
      in[j] = 1;
      c.evaluate_t(out, in);
      for (long k = 0; k < i; k++){
	val_t[k][j] = out[k];
      }
    }
    
    cout << i << " " << (val_t == transpose(val)) << " ";

    Mat<zz_p> val_inv_t;
    val_inv_t.SetDims(i, i);
    for (long j = 0; j < i; j++){
      Vec<zz_p> in, out;
      in.SetLength(i);
      in[j] = 1;
      c.interpolate_t(out, in);
      for (long k = 0; k < i; k++){
	val_inv_t[k][j] = out[k];
      }
    }

    cout << (val_inv_t == inv(val_t)) << endl;
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
