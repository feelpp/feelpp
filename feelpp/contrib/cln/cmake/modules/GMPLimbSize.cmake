
include(CheckCXXSourceCompiles)

set(_save_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
set(_save_required_includes ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES "${GMP_LIBRARIES} ${CMAKE_REQUIRED_LIBRARIES}")
set(CMAKE_REQUIRED_INCLUDES "${GMP_INCLUDE_DIR} ${CMAKE_REQUIRED_INCLUDES}")

CHECK_CXX_SOURCE_COMPILES("
	#include <gmp.h>
	template<bool COND> struct Static_Assert;
	template<> struct Static_Assert<true> { };
	#if defined(__GMP_BITS_PER_MP_LIB)
	Static_Assert<8*sizeof(mp_limb_t) == __GMP_BITS_PER_MP_LIB> check;
	int main() { return 0; }
	#endif
	"
	cl_gmp_has_nails
)
if (${cl_gmp_has_nails})
    message(SEND_ERROR "nails in MP limbs are not supported.")
endif()

CHECK_CXX_SOURCE_COMPILES("
	#include <gmp.h>
	template<bool COND> struct Static_Assert;
	template<> struct Static_Assert<true> { };
	Static_Assert<sizeof(mp_limb_t) == sizeof(long)> check;
	int main() { return 0; }
	"
	mp_limb_t_is_long
)

CHECK_CXX_SOURCE_COMPILES("
	#include <gmp.h>
	template<bool COND> struct Static_Assert;
	template<> struct Static_Assert<true> { };
	Static_Assert<sizeof(mp_limb_t) == sizeof(long long)> check;
	int main() { return 0; }
	"
	mp_limb_t_is_long_long
)
set(CMAKE_REQUIRED_INCLUDES ${_save_required_includes})
set(CMAKE_REQUIRED_LIBRARIES ${_save_required_libraries})

if (mp_limb_t_is_long)
    message(STATUS "sizeof(mp_limb_t) == sizeof(long)")
    set(GMP_DEMANDS_UINTD_LONG 1 CACHE INTERNAL "sizeof(mp_limb_t) == sizeof(long)")
else (mp_limb_t_is_long)
    if (mp_limb_t_is_long_long)
        message(STATUS "sizeof(mp_limb_t) == sizeof(long long)")
        set(GMP_DEMANDS_UINTD_LONG_LONG 1 CACHE INTERNAL "sizeof(mp_limb_t) == sizeof(long long)")
    else()
        message(SEND_ERROR "Don't know which C integer type has sizeof(mp_limb_t)")
    endif()
endif (mp_limb_t_is_long)
