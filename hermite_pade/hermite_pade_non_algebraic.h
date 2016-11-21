#include "hermite_pade.h"

class hermite_pade_non_algebraic : public hermite_pade{
	NTL::Vec<NTL::ZZX> vec_fs;
	
	void set_up_bmc() override;
	
	void init();
	
	void mul_M_right() override;
	
	public:
	hermite_pad_non_algebraic(const NTL::Vec<NTL::ZZX> &fs,
	                          const Vec<long>& type,
	                          long prec,
	                          long fft_init = 0);
};
