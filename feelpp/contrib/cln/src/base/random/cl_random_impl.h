// cl_random implementation

#ifndef _CL_RANDOM_IMPL_H
#define _CL_RANDOM_IMPL_H

#include "cln/random.h"

namespace cln {

// random_UDS(randomstate,MSDptr,len) füllt die UDS MSDptr/len/..
// mit len Zufallsdigits.
// > randomstate: ein Random-State, wird verändert
// > MSDptr/len/..: wo die Zufallsdigits abgelegt werden sollen
// > len: gewünschte Anzahl von Zufallsdigits
extern void random_UDS (random_state& randomstate, uintD* MSDptr, uintC len);

// testrandom_UDS(randomstate,MSDptr,len) füllt die UDS MSDptr/len/..
// mit len Zufallsdigits, mit relativ langen Blöcken von Nullen und Einsen.
// > randomstate: ein Random-State, wird verändert
// > MSDptr/len/..: wo die Zufallsdigits abgelegt werden sollen
// > len: gewünschte Anzahl von Zufallsdigits
extern void testrandom_UDS (random_state& randomstate, uintD* MSDptr, uintC len);

}  // namespace cln

#endif /* _CL_RANDOM_IMPL_H */
