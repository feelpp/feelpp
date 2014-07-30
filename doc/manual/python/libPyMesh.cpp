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


#define SIMPLEX(_,n,type) type<n+1,1>("Simplex");
#define HYPERCUBE(_,n,type) type<n+1,1>("Hypercube");

#define PCH1(_,n,type) type<1,n+1>();
#define PCH2(_,n,type) type<2,n+1>();
#define PCH3(_,n,type) type<3,n+1>();

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;


    template<typename MeshType, int N>
void expo_w ( boost::shared_ptr<MeshType> m)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}

    template<typename MeshType>
boost::shared_ptr<MeshType> loadMesh_w (MeshType* mesh)
{
    return loadMesh(_mesh=mesh);
}

template<int Order,
         int Tag = 0,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename MeshType>

typename meta::Pch<MeshType,Order,Tag,Pts>::ptrtype
Pch_w( boost::shared_ptr<MeshType> mesh)
{
    return Pch<Order,Tag,Pts,MeshType>(mesh);
}




    template <int n,int N>
void def_wrapper (std::string s)
{    
    std::ostringstream f;
    std::ostringstream g;

    if(s.compare("Simplex")==0)
    {
        f<<"Simplex"<<n;
        g<<"MeshS"<<n;
        class_<Feel::Simplex<n>>(f.str().c_str(),init<>());
        class_<Feel::Mesh<Simplex<n>>,boost::shared_ptr<Feel::Mesh<Simplex<n>>>,boost::noncopyable>(g.str().c_str(),init<>())
            .def("new",&Feel::Mesh<Simplex<n>>::New)
            .staticmethod("new");

        def("loadMesh",loadMesh_w<Mesh<Simplex<n>>>);
        def("export",expo_w<Mesh<Simplex<n>>,N>);

    }

    else if(s.compare("Hypercube")==0)
    {
        f<<"Hypercube"<<n;
        g<<"MeshH"<<n;   
        class_<Feel::Hypercube<n>>(f.str().c_str(),init<>());
        class_<Feel::Mesh<Hypercube<n>>,boost::shared_ptr<Feel::Mesh<Hypercube<n>>>,boost::noncopyable>(g.str().c_str(),init<>())
            .def("new",&Feel::Mesh<Hypercube<n>>::New)
            .staticmethod("new");

        def("loadMesh",loadMesh_w<Mesh<Hypercube<n>>>);
        def("export",expo_w<Mesh<Hypercube<n>>,N>);
    }

}

    template <int n,int k>
void def_wrapper_Pch ()
{
    std::ostringstream h;
    std::ostringstream i;
    std::ostringstream j;

    j<<"newPch"<<k;

    h<<"PchS"<<n<<k;
    i<<"FunctSpaceS"<<n<<k;
    
    class_<Feel::meta::Pch<Mesh<Simplex<n>>,k>,boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<n>>,k>>>(h.str().c_str(),no_init);

    class_<Feel::FunctionSpace<Mesh<Simplex<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>,boost::shared_ptr<Feel::FunctionSpace<Mesh<Simplex<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>>,boost::python::bases<Feel::FunctionSpaceBase>>(i.str().c_str(),no_init);

    def(j.str().c_str(),Pch_w<k,0,PointSetEquiSpaced,Mesh<Simplex<n>>>);

    h.str("");
    i.str(""); 

    h<<"PchH"<<n<<k;
    i<<"FunctSpaceH"<<n<<k;

    class_<Feel::meta::Pch<Mesh<Hypercube<n>>,k>,boost::shared_ptr<Feel::meta::Pch<Mesh<Hypercube<n>>,k>>>(h.str().c_str(),no_init);

    class_<Feel::FunctionSpace<Mesh<Hypercube<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>,boost::shared_ptr<Feel::FunctionSpace<Mesh<Hypercube<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>>,boost::python::bases<Feel::FunctionSpaceBase>>(i.str().c_str(),no_init);

    def(j.str().c_str(),Pch_w<k,0,PointSetEquiSpaced,Mesh<Hypercube<n>>>);

    h.str("");
    i.str(""); 
}




BOOST_PYTHON_MODULE(libPyMesh)
{

    if (import_mpi4py()<0) return ;

    class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
        .def("worldComm",&Feel::detail::Environment::worldComm,return_value_policy<copy_non_const_reference>())
        .staticmethod("worldComm");

    class_<WorldComm>("WorldComm",init<>());

    class_<Feel::FunctionSpaceBase>("FunctionSpaceBase",no_init);

    /*
    def_wrapper<1,1>("Simplex");
    def_wrapper<2,1>("Simplex");
    def_wrapper<3,1>("Simplex");
    def_wrapper<1,1>("Hypercube");
    def_wrapper<2,1>("Hypercube");
    def_wrapper<3,1>("Hypercube");
    
    
    def_wrapper_Pch<1,1>();
    def_wrapper_Pch<1,2>();
    def_wrapper_Pch<1,3>();

    def_wrapper_Pch<2,1>();
    def_wrapper_Pch<2,2>();
    def_wrapper_Pch<2,3>();

    def_wrapper_Pch<3,1>();
    def_wrapper_Pch<3,2>();
    def_wrapper_Pch<3,3>();
    */

    
    BOOST_PP_REPEAT(3,SIMPLEX,def_wrapper)
    BOOST_PP_REPEAT(3,HYPERCUBE,def_wrapper)


    BOOST_PP_REPEAT(3,PCH1,def_wrapper_Pch)
    BOOST_PP_REPEAT(3,PCH2,def_wrapper_Pch)
    BOOST_PP_REPEAT(3,PCH3,def_wrapper_Pch)
}

