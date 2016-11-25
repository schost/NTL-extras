#include "hermite_pade_exact.h"

hermite_pade_exact::hermite_pade_exact(const ZZX &f, const Vec<long> &type, long prec_inp, long fft_index){
  hermite_pade_algebraic hp(f,type,prec_inp,fft_index);
  new_type = hp.update_type();
  cout << "new type: " << new_type << endl;
  for (long i  = 0; i < new_type.length(); i++)
    new_type[i] = max(0, new_type[i]);
  long total = 0;
  for (long i = 0; i < new_type.length(); i++)
    total += new_type[i]+1;
  cout <<"t from exact: " << new_type << " " << total<< endl;
  cout << trunc(f,total-1) << endl;
  HP = new hermite_pade_algebraic(trunc(f,total-1),new_type,prec_inp,fft_index);
}

hermite_pade_exact::hermite_pade_exact(const ZZX &f, const ZZ &denom, const Vec<long> &type, long prec_inp, long fft_index){
  hermite_pade_algebraic hp(f,denom,type,prec_inp,fft_index);
  new_type = hp.update_type();
  cout << "new type: " << new_type << endl;
  for (long i  = 0; i < new_type.length(); i++)
    new_type[i] = max(0, new_type[i]);
  long total = 0;
  for (long i = 0; i < new_type.length(); i++)
    total += new_type[i]+1;
  cout <<"t from exact: " << new_type << " " << total<< endl;
  cout << trunc(f,total-1) << endl;
  HP = new hermite_pade_algebraic(trunc(f,total-1),denom,new_type,prec_inp,fft_index);
}


hermite_pade_exact::hermite_pade_exact(const Vec<ZZX> &f, const Vec<long> &type, long prec_inp, long fft_index){

}

Vec<long> hermite_pade_exact::find_sol(Vec<Vec<ZZ>> &sol){
  HP->find_rand_sol(sol);
  return new_type;
}
