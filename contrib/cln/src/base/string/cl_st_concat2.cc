// cl_string concatenation.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/string.h"


// Implementation.
#include "base/string/cl_st_make0.h"

namespace cln {

const cl_string operator+ (const char* str1, const cl_string& str2)
{
    unsigned long len1 = ::strlen(str1);
    unsigned long len2 = strlen(str2);
    var cl_heap_string* str = cl_make_heap_string(len1+len2);
    var char * ptr = &str->data[0];
    {
        var const char * ptr1 = asciz(str1);
        for (var unsigned long count = len1; count > 0; count--)
            *ptr++ = *ptr1++;
    }
    {
        var const char * ptr2 = asciz(str2);
        for (var unsigned long count = len2; count > 0; count--)
            *ptr++ = *ptr2++;
    }
    *ptr++ = '\0';
    return str;
}

}  // namespace cln
