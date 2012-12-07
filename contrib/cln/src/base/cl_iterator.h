// Abstract iterators.

#ifndef _CL_ITERATOR_H
#define _CL_ITERATOR_H

#include "cln/types.h"


// An iterator's typical use is a loop, but you have an abstraction over
// the loop's initialization, step and end-test.
// Example:
//    foo_iterator foo_loop = ...;
//    while (!foo_loop.endp()) {
//        foo element = foo_loop.next();
//        ...
//    }
// It is allowed to call endp() as many times as you want, and to terminate
// the loop any time you want.

template <class T>
class cl_abstract_iterator {
public:
	virtual bool endp () = 0;
	virtual T& next () = 0;
};


#endif /* _CL_ITERATOR_H */
