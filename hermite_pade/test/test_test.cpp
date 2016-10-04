#include "hermite_pade.h"

int main(int argc, char* argv[]){
  ZZX f;
  ZZ x;
  int N = 0;
  while (cin >> x) {
    SetCoeff(f,N++,x);
  }
  //cout << f << endl;
  Vec<long> type;
  type.append(3);
  type.append(4);
  hermite_pade hp(f,type,N,stoi(argv[1]));
  Vec<Vec<ZZ>> v;
  hp.find_rand_sol(v);
  cout << v << endl;
}
