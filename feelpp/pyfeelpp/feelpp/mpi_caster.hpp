#ifndef FEELPP_PYBIND11_BOOST_MPI_COMM
#define FEELPP_PYBIND11_BOOST_MPI_COMM

#include <boost/mpi/communicator.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/pybind11.h>

// Import mpi4py on demand
#define VERIFY_MPI4PY(func)                                                    \
  if (!func)                                                                   \
  {                                                                            \
    int rc = import_mpi4py();                                                  \
    if (rc != 0)                                                               \
    {                                                                          \
      std::cout << "ERROR: could not import mpi4py!" << std::endl;             \
      throw std::runtime_error("Error when importing mpi4py");                 \
    }                                                                          \
  }

namespace pybind11
{
namespace detail
{
template <>
class type_caster<boost::mpi::communicator>
{
public:
  // Define this->value of type MPICommWrapper
  PYBIND11_TYPE_CASTER(boost::mpi::communicator, _("MPICommunicator"));

  // Python to C++
  bool load(handle src, bool)
  {
    // Simplified version of isinstance(src, mpi4py.MPI.Comm) - avoids
    // segfault when pybind11 tries to convert some other random type to
    // MPICommWrapper
    if (not hasattr(src, "Allgather"))
      return false;
    VERIFY_MPI4PY(PyMPIComm_Get);
    value = boost::mpi::communicator( *PyMPIComm_Get( src.ptr() ), boost::mpi::comm_attach );
    return true;
  }

  // C++ to Python
  static handle cast(boost::mpi::communicator src,
                     pybind11::return_value_policy policy, handle parent)
  {
    VERIFY_MPI4PY(PyMPIComm_New);
    return pybind11::handle(PyMPIComm_New(src));
  }

  operator boost::mpi::communicator() { return this->value; }
};
} // namespace detail
} // namespace pybind11
#endif