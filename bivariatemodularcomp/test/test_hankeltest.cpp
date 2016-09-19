#include "mosaic_hankel.h"
using namespace NTL;
using namespace std;

int main(){
  zz_p::init(13);
  Vec<zz_p> v;
  for (int i = 0; i < 5; i++)
    v.append(zz_p{i});
  hankel h{v,4,5};
  Mat<zz_p> M;
  to_dense(M,h);
  cout << M << endl;
}
