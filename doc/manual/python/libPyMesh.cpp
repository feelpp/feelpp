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

#include<feel/feelcore/environment.hpp>
#include<feel/feelfilters/loadmesh.hpp>
#include<feel/feelfilters/exporter.hpp>
#include<feel/feelfilters/detail/mesh.hpp>


using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

// redefine exporter and loadMesh methods, created from a BOOST_PARAMETER_FUNCTION, into simple definition
 
    template<typename MeshType,int N>
void expo_w ( boost::shared_ptr<MeshType> m)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}

boost::shared_ptr<Mesh<Simplex<2>>> loadMesh_w (Mesh<Simplex<2>>* mesh)
{
    return loadMesh(_mesh=mesh);
}

// creation of the libPyMesh library, that we will use in the Python script 

BOOST_PYTHON_MODULE(libPyMesh)
{

    if (import_mpi4py()<0) return ;

// definition of the Environment object and methods and classes link to it 
    class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
        .def("worldComm",&Feel::detail::Environment::worldComm,return_value_policy<copy_non_const_reference>())
        .staticmethod("worldComm");

    class_<WorldComm>("WorldComm",init<>());


// definition of the geometrical object (Simplex and Hypercube) and of the Mesh class 
    class_<Feel::Simplex<2>>("Simplex",init<>())
        .def("dim",&Feel::Simplex<2>::dimension);

    class_<Feel::Hypercube<2>>("Hypercube",init<>())
        .def("dim",&Feel::Hypercube<2>::dimension);



    class_<Feel::Mesh<Feel::Simplex<2>>,boost::shared_ptr<Feel::Mesh<Feel::Simplex<2>>>,boost::noncopyable>("Mesh",init<>())
        .def("new",&Feel::Mesh<Simplex<2>>::New)
        .staticmethod("new")
        .def("clear",&Feel::Mesh<Simplex<2>>::clear);


//definition of the loadMesh and exporter methods with functions define before 
    def("loadMesh",loadMesh_w);

    class_<ExporterEnsightGold<Mesh<Simplex<2>>,1>>("Exporter",init<WorldComm>())
        .def("setMesh",&Exporter<Mesh<Simplex<2>>>::setMesh) 
        .def("addRegions",&Exporter<Mesh<Simplex<2>>>::addRegions)
        .def("save",&ExporterEnsightGold<Mesh<Simplex<2>>,1>::save);


    def("export",expo_w<Mesh<Simplex<2>>,1>);
}
