#ifndef __MAT_ZZpx_H__
#define __MAT_ZZpX_H__

#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>

void add(NTL::Mat<NTL::ZZ_pX>&, const NTL::Mat<NTL::ZZ_pX>&, const NTL::Mat<NTL::ZZ_pX>&);
void mul(NTL::Mat<NTL::ZZ_pX>&, const NTL::Mat<NTL::ZZ_pX>&, const NTL::Mat<NTL::ZZ_pX>&);
#endif
