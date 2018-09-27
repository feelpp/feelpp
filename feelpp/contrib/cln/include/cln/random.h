// Public random number operations.

#ifndef _CL_RANDOM_H
#define _CL_RANDOM_H

#include "cln/types.h"
#include "cln/modules.h"

namespace cln {

class random_state {
public:
	struct { uint32 hi; uint32 lo; } seed;
// Constructor:
	random_state ();
};

// random32(randomstate) liefert eine neue Zufallszahl.
// > randomstate: ein Random-State, wird verändert
// < ergebnis: eine 32-Bit-Zufallszahl
extern uint32 random32 (random_state& randomstate);

#if defined(HAVE_FAST_LONGLONG)
// random64(randomstate) liefert eine neue Zufallszahl.
// > randomstate: ein Random-State, wird verändert
// < ergebnis: eine 64-Bit-Zufallszahl
inline uint64 random64 (random_state& randomstate)
{
	return ((uint64)random32(randomstate) << 32)
	       | (uint64)random32(randomstate);
}
#endif

// Ein globaler Zufallszahlengenerator.
extern random_state default_random_state;
class cl_random_def_init_helper
{
	static int count;
public:
	cl_random_def_init_helper();
	~cl_random_def_init_helper();
};
static cl_random_def_init_helper cl_random_def_init_helper_instance;

// Das ist der Default-Generator.
inline uint32 random32 (void)
	{ return random32(default_random_state); }
#if defined(HAVE_FAST_LONGLONG)
inline uint64 random64 (void)
	{ return random64(default_random_state); }
#endif

}  // namespace cln

#endif /* _CL_RANDOM_H */
