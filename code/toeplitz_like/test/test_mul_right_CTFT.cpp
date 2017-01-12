#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_p_toeplitz.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(0);
  
  for (long i1 = 700; i1 < 701; i1 += 1)
    for (long i2 = 700; i2 < 701; i2 += 1){
      
      Vec<zz_p> dat;
      random(dat, i1+i2-1);
      lzz_p_toeplitz h(dat, i1, i2);
      Vec<zz_p> input, output, output2;
      random(input, i2);
      
      if (opt == 1){
	h.mul_right(output, input);
	h.mul_right_CTFT(output2, input);
	assert (output == output2);
	cout << i1 << " " << i2 << endl;
      }
      else{
	cout << i1 << " " << i2 << " ";
	double t;
	
	t = GetTime();
	for (long j = 0; j < 40000; j++)
	  h.mul_right(output, input);
	t = GetTime() - t;
	cout << t << " ";
	
	t = GetTime();
	for (long j = 0; j < 40000; j++)
	  h.mul_right_CTFT(output, input);
	t = GetTime() - t;
	cout << t << " ";
	
	cout << endl;
      }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
