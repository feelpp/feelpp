#include <feel/feel.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
#include <boost/parameter/python.hpp>
#include <boost/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelfilters/loadmesh.hpp>
#include<feel/feelfilters/exporter.hpp>
#include<feel/feelfilters/detail/mesh.hpp>

#define BOOST_PYTHON_MAX_ARITY 20

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

struct loadMesh2_fwd
{
    template<class A0>
        boost::shared_ptr<Mesh<Simplex<2>>> operator() ( 
                boost::type<boost::shared_ptr<Mesh<Simplex<2>>>>, A0 const& a0)
        {
           return loadMesh2(a0);
        }
};

BOOST_PYTHON_MODULE(libPyFeelpp)
{
   
   if (import_mpi4py()<0) return ;

 class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>());

class_<Feel::Simplex<2>>("Simplex",init<>())
        .def("dim",&Feel::Simplex<2>::dimension);


    class_<Feel::Mesh<Feel::Simplex<2>>,boost::shared_ptr<Feel::Mesh<Feel::Simplex<2>>>,boost::noncopyable>("Mesh",init<>())
        .def("new",&Feel::Mesh<Simplex<2>>::New)
        .staticmethod("new")
        .def("clear",&Feel::Mesh<Simplex<2>>::clear);


    //class_<Exporter<Mesh<Simplex<2>>>,boost::shared_ptr<Exporter<Mesh<Simplex<2>>>>>("Exporter",init<>());
    
   def("loadMesh",loadMesh3);
  // def("exporter",exporter2);

/*
      def(
            "loadmesh2",
            py::function<
            loadMesh2_fwd,
            mpl::vector<
            boost::shared_ptr<Mesh<Simplex<2>>>,
            tag::mesh (Mesh<Simplex<2>>*)
            >
            >()
       );
*/

}
