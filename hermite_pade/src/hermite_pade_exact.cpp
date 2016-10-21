#include "hermite_pade_exact.h"

hermite_pade_exact::hermite_pade_exact(const ZZX &f, const Vec<long> &type, long prec_inp, long fft_index): HP(f,type,prec_inp,fft_index){
  new_type = HP.update_type();
  long total = 0;
  for (long i = 0; i < new_type.length(); i++)
    total += new_type[i]+1;
  cout <<"t from exact: " << new_type << " " << total<< endl;
  cout << trunc(f,total-1) << endl;
  HP = hermite_pade(trunc(f,total-1),new_type,prec_inp,fft_index);
}

Vec<long> hermite_pade_exact::find_sol(Vec<Vec<ZZ>> &sol){
  HP.find_rand_sol(sol);
  return new_type;
}
