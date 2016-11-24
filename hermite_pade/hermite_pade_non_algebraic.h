#include "hermite_pade.h"

class hermite_pade_non_algebraic : public hermite_pade{
	NTL::Vec<NTL::ZZX> vec_fs;
	NTL::Vec<NTL::ZZX> vec_added;
	
	/*---------------------------------------------------------*/
	/* tries to increase the rank of M by adding random blocks */
	/*---------------------------------------------------------*/
	Vec<hankel> increase_rank();
	
	void set_up_bmc() override;
	
	Vec<ZZ_p> mul_M_right(const NTL::Vec<NTL::ZZ_p> &) override;
	
	public:
	hermite_pade_non_algebraic(const NTL::Vec<NTL::ZZX> &fs,
	                          const Vec<long>& type,
	                          long fft_init = 0);
};
