#ifndef __MAT_ZZX_H__
#define __MAT_ZZX_H__

#include <NTL/matrix.h>
#include <NTL/ZZX.h>

void add(NTL::Mat<NTL::ZZX>&, const NTL::Mat<NTL::ZZX>&, const NTL::Mat<NTL::ZZX>&);
void mul(NTL::Mat<NTL::ZZX>&, const NTL::Mat<NTL::ZZX>&, const NTL::Mat<NTL::ZZX>&);
#endif
