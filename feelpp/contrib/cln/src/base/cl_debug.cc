// Debugging support for dynamic typing.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/object.h"


// Implementation.

#include "cln/io.h"

namespace cln {

// The default printer function.
void cl_dprint_unknown (cl_heap* pointer)
{
	fprint(cl_debugout, "<unknown @0x");
	fprinthexadecimal(cl_debugout, (uintptr_t) pointer);
	fprint(cl_debugout, " refcount=");
	fprintdecimal(cl_debugout, pointer->refcount);
	fprint(cl_debugout, " type=");
	fprinthexadecimal(cl_debugout, (uintptr_t) pointer->type);
	fprint(cl_debugout, ">");
}

static void cl_dprint_unknown_immediate (cl_heap* pointer)
{
	fprint(cl_debugout, "<unknown @0x");
	fprinthexadecimal(cl_debugout, (uintptr_t) pointer);
	fprint(cl_debugout, ">");
}

// Print an object. This function is callable from the debugger.
extern "C" void* cl_print (cl_uint word);
void* cl_print (cl_uint word)
{
	var cl_heap* pointer = (cl_heap*)word;
	if (cl_pointer_p(word)) {
		var const cl_class* type = pointer->type;
		if (type->dprint)
			type->dprint(pointer);
		else
			cl_dprint_unknown(pointer);
	} else {
		var const cl_class* type = cl_immediate_classes[cl_tag(word)];
		if (type && type->dprint)
			type->dprint(pointer);
		else
			cl_dprint_unknown_immediate(pointer);
	}
	cl_debugout << std::endl; // newline and flush output
	return pointer;
}

void cl_gcobject::debug_print () const
{
	cl_print(word);
}

void cl_gcpointer::debug_print () const
{
	cl_print(word);
}

void cl_rcobject::debug_print () const
{
	cl_print(word);
}

void cl_rcpointer::debug_print () const
{
	cl_print(word);
}

}  // namespace cln
