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


//#define BOOST_PYTHON_MAX_ARITY 20

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

/*
struct loadMesh2_fwd
{
    template<class A0>
        boost::shared_ptr<Mesh<Simplex<2>>> operator() ( 
                boost::type<boost::shared_ptr<Mesh<Simplex<2>>>>, A0 const& a0)
        {
            return loadMesh2(a0);
        }
};



   template<typename MeshType,int N>
   boost::shared_ptr<Exporter<MeshType,N> > (Exporter<MeshType,N>::*New1) (po::variables_map const&,std::string,WorldComm const&) = &Exporter<MeshType,N>::New;
 


   boost::shared_ptr<Exporter<Mesh<Simplex<2>>,1>> (Exporter<Mesh<Simplex<2>>,1>::*New1) (po::variables_map const&,std::string,WorldComm const&) = &Exporter<Mesh<Simplex<2>>,1>::New;

   boost::shared_ptr<Exporter<Mesh<Simplex<2>>,1>> (Exporter<Mesh<Simplex<2>>,1>::*New2) (std::string const&,std::string,WorldComm const&) = &Exporter<Mesh<Simplex<2>>,1>::New;
 


boost::shared_ptr<Exporter<Mesh<Simplex<2>>,1>> New1 (po::variables_map const& x,std::string y,WorldComm const& z) 
{
    return Exporter<Mesh<Simplex<2>>,1>::New(x,y,z);
}

boost::shared_ptr<Exporter<Mesh<Simplex<2>>,1>> New2 () 
{
    return Exporter<Mesh<Simplex<2>>,1>::New();
}

    template<typename MeshType, int N>
void setMesh1 (ExporterEnsightGold<Mesh<Simplex<2>>,1> e,boost::shared_ptr<MeshType> x)
{
    e.setMesh(x);

}



void expo ( boost::shared_ptr<Mesh<Simplex<2>>> m)
{
    auto x=Exporter<Mesh<Simplex<2>>,1>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}
*/

    template<typename MeshType,int N>
void expo2 ( boost::shared_ptr<MeshType> m)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}

/*
    template<typename MeshType,int N>
boost::shared_ptr<Exporter<MeshType,N>> expo3 ( boost::shared_ptr<MeshType> m,WorldComm w)
{
    auto x=ExporterEnsightGold<MeshType,N>::New();
    x->setMesh(m);
    return x;
}
*/

boost::shared_ptr<Mesh<Simplex<2>>> loadMesh3 (Mesh<Simplex<2>>* mesh)
{
    return loadMesh(_mesh=mesh);
}



BOOST_PYTHON_MODULE(libPyFeelpp)
{

    if (import_mpi4py()<0) return ;

    class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
       .def("worldComm",&Feel::detail::Environment::worldComm,return_value_policy<copy_non_const_reference>())
       .staticmethod("worldComm");
  
    /*
       .def("vm",&Feel::detail::Environment::about)
       .staticmethod("vm");
     */


    class_<Feel::Simplex<2>>("Simplex",init<>())
        .def("dim",&Feel::Simplex<2>::dimension);

    class_<Feel::Hypercube<2>>("Hypercube",init<>())
        .def("dim",&Feel::Hypercube<2>::dimension);



    class_<Feel::Mesh<Feel::Simplex<2>>,boost::shared_ptr<Feel::Mesh<Feel::Simplex<2>>>,boost::noncopyable>("Mesh",init<>())
        .def("new",&Feel::Mesh<Simplex<2>>::New)
        .staticmethod("new")
        .def("clear",&Feel::Mesh<Simplex<2>>::clear);



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


    class_<WorldComm>("WorldComm",init<>());

    class_<ExporterEnsightGold<Mesh<Simplex<2>>,1>>("Exporter",init<WorldComm>())
        //.def("new",New1<Mesh<Simplex<2,1,2>,double,0>,1>)
        //.def("new",&Exporter<Mesh<Simplex<2>>>::New1)
        .def("setMesh",&Exporter<Mesh<Simplex<2>>>::setMesh) 
        .def("addRegions",&Exporter<Mesh<Simplex<2>>>::addRegions)
        .def("save",&ExporterEnsightGold<Mesh<Simplex<2>>,1>::save);


   // def("new",New2);
   //def("setMesh1",setMesh1<Mesh<Simplex<2>>,1>); 
   // def("export0",expo);
    def("export",expo2<Mesh<Simplex<2>>,1>);
   //def("export2",expo3<Mesh<Simplex<2>>,1>);

    /////////////////
}
