#include "hermite_pade.h"

class hermite_pade_exact{
	hermite_pade HP;
	Vec<long> new_type;
	
	public:
	hermite_pade_exact(const ZZX &f, const Vec<long> &type, long prec_inp, long fft_index);
	
	Vec<long> find_sol(Vec<Vec<ZZ>> &sol);
};
