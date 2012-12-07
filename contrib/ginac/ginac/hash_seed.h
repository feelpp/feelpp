#ifndef GINAC_HASH_SEED_H
#define GINAC_HASH_SEED_H
/**
 * We need a hash function which gives different values for objects of
 * different types. Hence we need some unique integer for each type.
 * Fortunately, standard C++ RTTI class `type_info' stores a pointer to
 * mangled type name. Normally this pointer is the same for all objects
 * of the same type (although it changes from run to run), so it can be
 * used for computing hashes. However, on some platforms (such as woe32)
 * the pointer returned by type_info::name() might be different even for
 * objects of the same type! Hence we need to resort to comparing string
 * representation of the (mangled) type names. This is quite expensive,
 * so we compare crc32 hashes of those strings. We might got more hash
 * collisions (and slower evaluation as a result), but being a bit slower
 * is much better than being wrong.
 */
#include <typeinfo>
#include <cstring>
#include "crc32.h"
#include "utils.h"
#ifdef _WIN32
#define GINAC_HASH_USE_MANGLED_NAME 1
#endif
namespace GiNaC
{
#ifndef GINAC_HASH_USE_MANGLED_NAME
static inline unsigned make_hash_seed(const std::type_info& tinfo)
{
	// this pointer is the same for all objects of the same type.
	// Hence we can use that pointer 
	const void* mangled_name_ptr = (const void*)tinfo.name();
	unsigned v = golden_ratio_hash((p_int)mangled_name_ptr);
	return v;
}
#else
static unsigned make_hash_seed(const std::type_info& tinfo)
{
	const char* mangled_name = tinfo.name();
	return crc32(mangled_name, std::strlen(mangled_name), 0);
}
#endif
} // namespace GiNaC
#endif /* GINAC_HASH_SEED_H */


