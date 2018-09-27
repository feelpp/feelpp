// cl_fget().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/string.h"


// Implementation.

#include "cln/io.h"
#include "base/string/cl_spushstring.h"

namespace cln {

const cl_string cl_fget (std::istream& stream, char delim)
{
	var cl_spushstring buffer;
	// Handling of eofp is tricky: EOF is reached when (!stream.good()) || (stream.eof()).
	while (stream.good()) {
		var int c = stream.get();
		if (stream.eof())
			break;
		if (c==delim) {
			stream.unget();
			break;
		}
		buffer.push(c);
	}
	return buffer.contents();
}

}  // namespace cln
