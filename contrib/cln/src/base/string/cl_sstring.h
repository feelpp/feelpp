// Simple strings.

#ifndef _CL_SSTRING_H
#define _CL_SSTRING_H

namespace cln {

// Liefert einen String.
// Mit malloc_hook() alloziert, mit free_hook() freizugeben.
extern char * cl_sstring (const char * ptr, uintC len);

}  // namespace cln

#endif /* _CL_SSTRING_H */
