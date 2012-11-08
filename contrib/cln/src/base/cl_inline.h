/*
 * Selectively inline a function in *some* translation units.
 * See cl_maybe_inline.h file for the explanation.
 */
#undef CL_INLINE
#undef CL_INLINE_DECL
#define CL_INLINE static inline
#define CL_INLINE_DECL(fcn) CL_INLINE_HINT fcn ## _inline
