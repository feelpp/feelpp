// cl_free_heap_object().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/object.h"


// Implementation.

#include "cln/malloc.h"

namespace cln {

void cl_free_heap_object (cl_heap* pointer)
{
	// This is invoked when pointer->refcount gets decremented to 0.
	var const cl_class* type = pointer->type;
	if (type->destruct)
		type->destruct(pointer);
	free_hook(pointer);
}


// The best place to put the free software license is cl_free.o.

// NB about #ident: To see the strings in the .comment section of an ELF
// executable, use "strings < executable", not "strings executable".
// Better put the license into the data section: it is more portable (some
// C++ compilers may not understand #ident), is not lost in object formats
// like a.out, and is taken into account by "intelligent" strings commands.

static const char * copyright_notice[] = {
  "                                                                    \n"
  "Copyright (c)      Bruno Haible 1988-2008                           \n"
  "Copyright (c)   Richard Kreckel 2000-2012                           \n"
  "Copyright (c) Alexei Sheplyakov 2008-2010                           \n"
  "                                                                    \n"
  "This program is free software; you can redistribute it and/or modify\n"
  "it under the terms of the GNU General Public License as published by\n"
  "the Free Software Foundation; either version 2, or (at your option) \n"
  "any later version.                                                  \n"
  "                                                                    \n"
  "This program is distributed in the hope that it will be useful, but \n"
  "WITHOUT ANY WARRANTY; without even the implied warranty of          \n"
  "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   \n"
  "General Public License for more details.                            \n"
  "                                                                    \n"
  "You should have received a copy of the GNU General Public License   \n"
  "along with this program; if not, write to the Free Software         \n"
  "Foundation, 51 Franklin Street, Fifth Floor, Boston, MA             \n"
  "02110-1301, USA.\n"
  "                                                                    ",
  (const char *) &copyright_notice
};

}  // namespace cln
