// (C) Copyright 2005 Matthias Troyer 

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Matthias Troyer

#include <boost/archive/impl/archive_pointer_oserializer.ipp>
#include <boost/mpi/detail/mpi_datatype_cache.hpp>

namespace boost { namespace mpi { namespace detail {

mpi_datatype_map mpi_datatype_cache;
  
} } }
