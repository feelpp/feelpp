// cln/string.hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/string.h"


// Implementation.

namespace cln {

uintptr_t hashcode (const cl_string& str)
{
    var uintptr_t code = 0x61284AF3;
    // We walk through all characters. It may take some time for very
    // long strings, but it's better than completely ignoring some characters.
    var intptr_t len = str.size();
    var const char * ptr = str.asciz();
    for (; len > 0; len--) {
        var unsigned char c = *ptr++;
        code = (code << 5) | (code >> 27); // rotate left 5 bits
        code += (intptr_t)c << 16;
        code ^= (intptr_t)c;
        code &= 0xFFFFFFFF;
    }
    return code;
}

}  // namespace cln
