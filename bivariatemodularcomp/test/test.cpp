#include "bivariatemodularcomp.h"
#include "mat_ZZX.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <iomanip>
#include <random>
#include <functional>
using namespace NTL;
using namespace std;

ostream& operator<< (ostream& out, const vector<long>& v){
  for (auto &i: v) out << i << " ";
  return out;
}

// prints the resulting vectors to result.txt
int main (){
  default_random_engine gen;
  uniform_int_distribution<long> dist;
  auto rand = bind(dist,gen);
  ofstream fout("result.txt");
  string s;
  ZZX f;
  long f_deg = 2500;
  long type_size = 50;
  long type_entry = 49;
  long prec = 2500;
  ZZ p_field;
  Vec<ZZ>rhs_full;
  Vec<ZZ_p>rhs;
  Vec<ZZ>rhs_full2;
  Vec<ZZ_p>rhs2;
  Vec<long> type;
  istringstream iss;
  //reading the field (p)
  cout << "Enter the field" << endl;
  cin >> p_field;
  ZZ_p::init(p_field);
  //reading the polynomial
  cout << "Enter the degree" << endl;
  cin >> f_deg;
  cout << "Enter the precision" << endl;
  cin >> prec;
  for (long i = 0; i <= f_deg; i++)
    SetCoeff(f,i,rand());
  //reading the type
  cout << "Enter the type length" << endl;
  cin >> type_size;
  cout << "Enter the type entry" << endl;
  cin >> type_entry;
  for (long i = 0; i < type_size; i++)
    type.append(type_entry);
  // reading the rhs
  long to = type_size * (type_entry+1);
  for (long i = 0; i < to; i++){
    rhs_full.append(ZZ(rand()));
    rhs_full2.append(ZZ(rand()));
  }
  conv(rhs,rhs_full);
  ZZ_pX f_p;
  conv(f_p,f);
  cout << "f: " << f_p << endl;
  cout << "rhs: " << rhs << endl;
//  conv(rhs2,rhs_full2);
  auto p = BivariateModularComp(f,type,prec);
  clock_t t;
  ofstream ofs("result.txt");
  t = clock();
  auto m1 = p.mult(rhs);
 //auto m3 = p.mult(rhs2);
  cout << "Matrix mult took: " << clock()-t << endl;
  fout << "Matrix: " << m1 << endl;
  t = clock();
  auto m2 = p.mult_Horners(rhs);
 //auto m4 = p.mult_Horners(rhs2);
  cout << "Horner's rule took: " << clock()-t <<  endl;
  fout << "Horner's: " << m2 << endl;
  cout << "equal? " << boolalpha << (m1 == m2) << endl;
}
