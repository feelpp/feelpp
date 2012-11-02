// testrandom_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "base/random/cl_random_impl.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I testrandom_I (random_state& randomstate)
{
  var uint32 ran = random32(randomstate);
  var bool negative = (ran & 1);
  var bool algo = ((ran>>1) & 1);
  ran = ran >> 2;
  ran = ran & ((1<<8)-1);
  var uintC len =
    (ran == 0 ? 0 :
     ran <= 80 ? 1 :
     ran <= 128 ? 2 :
     ran <= 158 ? 3 :
     ran <= 172 ? 4 :
     ran <= 200 ? (ran-153)/4 : // 5..11
     ran-189 // 12..66
    );
  CL_ALLOCA_STACK;
  var uintD* MSDptr;
  num_stack_alloc_1(len,MSDptr=,);
  if (algo)
    { testrandom_UDS(randomstate,MSDptr,len); }
  else
    { random_UDS(randomstate,MSDptr,len); }
  var cl_I x = UDS_to_I(MSDptr,len);
  return (negative ? -x : x);
}

}  // namespace cln
