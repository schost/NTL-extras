#include <iostream>
#include <sstream>
#include "hermite_pade.h"
using namespace std;

int main(){
  ZZX f;
  Vec<long> type;
  string s;
  getline(cin,s);
  istringstream iss(s);
  long x;
  long i = 0;
  while(iss>>x){
    SetCoeff(f,i++,x);
  }
  cout << "f: " << f << endl;
  getline(cin,s);
  iss = istringstream(s);
  while (iss >> x){
    type.append(x);
  }
  cout << "type: " << type << endl;
  hermite_pade hp(f,type,deg(f)+1,10);
  Vec<Vec<ZZ>> sol;
  hp.find_rand_sol(sol);
  cout << sol << endl;
}
