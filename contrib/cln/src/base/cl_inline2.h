/*
 * Selectively inline a function in *some* translation units.
 * See cl_maybe_inline.h file for the explanation.
 */
#undef CL_INLINE2
#undef CL_INLINE2_DECL
#define CL_INLINE2 static inline
#define CL_INLINE2_DECL(fcn) CL_INLINE_HINT fcn ## _inline
