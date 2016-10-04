#include "hermite_pade.h"

int main(int argc, char* argv[]){
  ZZX f;
  int N = 8;
  for(long x = 1; x <= N; x++){
    SetCoeff(f,x-1,x);
  }
  cout << "f: " << f << endl;
  Vec<long> type;
  type.append(3);
  type.append(4);
  cout << "type: " << type << endl;
  hermite_pade hp(f,type,N,10);
  Vec<Vec<ZZ>> v;
  hp.find_rand_sol(v);
  cout << v << endl;
}
