// CLN internal macros extra
// This file must be included after other system include files.

#ifndef _CL_XMACROS_H
#define _CL_XMACROS_H

// Swap the contents of two variables:  swap(int, x1, x2);
  #define swap(swap_type,swap_var1,swap_var2)  \
    { var swap_type swap_temp;                                             \
      swap_temp = swap_var1; swap_var1 = swap_var2; swap_var2 = swap_temp; \
    }

#endif /* _CL_XMACROS_H */
