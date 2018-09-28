// Digit sequence level random number generator.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/random/cl_random_impl.h"


// Implementation.

#include "cln/random.h"
#include "base/digitseq/cl_DS.h"
#include "base/cl_low.h"

namespace cln {

void testrandom_UDS (random_state& randomstate, uintD* MSDptr, uintC len)
{
  // Idea from Torbjörn Granlund, see his "random2.c" file in gmp 2.0.
  var uintD* ptr = MSDptr mspop len;
  DS_clear_loop(MSDptr,len,ptr);
  var uintC bit_pos = 0;
  var uint32 ran = 0;
  var uintC ran_bits = 0;
  while (bit_pos < intDsize*len)
    { if (ran_bits < log2_intDsize+1)
        { ran = random32(randomstate); ran_bits = 32; }
      var uintL n_bits = (ran >> 1) % intDsize + 1; // number of bits
      if (ran & 1)
        { // put in a bit string of n_bits bits at position bit_pos.
          if (bit_pos + n_bits > intDsize*len)
            { n_bits = intDsize*len - bit_pos; }
          if (bit_pos / intDsize == (bit_pos + n_bits - 1) / intDsize)
            { // need to modify one digit
              lspref(ptr,bit_pos/intDsize) |= (((uintD)1 << n_bits) - 1) << (bit_pos%intDsize);
            }
            else
            { // need to modify two adjacent digits
              lspref(ptr,bit_pos/intDsize) |= ((uintD)(-1) << (bit_pos%intDsize));
              lspref(ptr,bit_pos/intDsize+1) |= (((uintD)1 << ((bit_pos+n_bits)%intDsize)) - 1);
            }
        }
      bit_pos = bit_pos + n_bits;
      ran = ran >> (log2_intDsize+1); ran_bits -= log2_intDsize+1;
    }
}

}  // namespace cln
