#ifndef __hermite_pade_h__
#define __hermite_pade_h__

#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include "mosaic_hankel.h"

// abstract base class for hermite pade classes
class hermite_pade{
	protected:
	NTL::Vec<long> type; // the type of the approximant
  long rank; // the rank of M
  long w; // w mod p
  long order; // order of w
  long p, prec, sizeX, sizeY;
  long added; // number of rows added
  
	/*----------------------------------------------------------------*/
  /* applies a block reversal to v                                  */
  /* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
  /* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                     
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> flip_on_type (const Vec<ZZ_p> &v);

  /*----------------------------------------------------------------*/
  /* applies a block reversal to v                                  */
  /* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
  /* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                     
  /*----------------------------------------------------------------*/
  Vec<Vec<ZZ>> flip_on_type(const Vec<Vec<ZZ>> &v);

  public:
  virtual void find_rand_sol (NTL::Vec<NTL::Vec<NTL::ZZ>> &sol) = 0;
  virtual ~hermite_pade();
};

#endif
