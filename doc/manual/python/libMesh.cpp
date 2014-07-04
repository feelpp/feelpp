#include <feel/feel.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/python.hpp>
#include <boost/python.hpp>
#include <boost/mpl/vector.hpp>

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

/*
template<uint16_type Dim,uint16_type Order,uint16_type RDim>
Feel::Simplex<Dim,Order,RDim>* ConstruSimplex ()
{
    return new Feel::Simplex<Dim,Order,RDim>();
}
*/


struct loadMesh_fwd
{
    template<class A0,class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8,class A9,class A10,class A11,class A12,class A13,class A14,class A15>
void operator() ( 
        boost::type<Feel::detail::mesh<Mesh<Simplex<2>>>::ptrtype>, &self, A0 const& a0, A1 const& a1,A2 const& a2,A3 const& a3,A4 const& a4,A5 const& a5,A6 const& a6,A7 const& a7,A8 const& a8,A9 const& a9,A10 const& a10,A11 const& a11,A12 const& a12,A13 const& a13,A14 const& a14,A15 const& a15
        )
    {
        self.loadMesh(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15);
    }
};


//template<uint16_type Dim,uint16_type Order,uint16_type RDim>
BOOST_PYTHON_MODULE(libMesh)
{
   
   if (import_mpi4py()<0) return ;
       
   class_<Feel::Simplex<2>>("Simplex",init<>())
    .def("dim",&Feel::Simplex<2>::dimension);
   
   
   class_<Feel::Mesh<Feel::Simplex<2>>,boost::noncopyable>("Mesh",init<>());


    def(
        "loadMesh",py::function<
                loadMesh_fwd
                ,mpl::vector<
                    Feel::detail::mesh<Mesh<Simplex<2>>>::ptrtype
                    , tag::mesh(*)
                    , tag::filename*(*(boost::is_convertible<mpl::_,std::string>))
                    , tag::desc*(*)
                    , tag::h*(*(boost::is_arithmetic<mpl::_>))
                    , tag::straighten*(bool)
                    , tag::refine*(*(boost::is_integral<mpl::_>))
                    , tag::update*(*(boost::is_integral<mpl::_>))
                    , tag::physical_are_elementary_regions*(bool)
                    , tag::worldcomm*(WorldComm)
                    , tag::force_rebuild*(*(boost::is_integral<mpl::_>))
                    , tag::respect_partition*(bool)
                    , tag::rebuild_partitions*(bool)
                    , tag::partitions*(*(boost::is_integral<mpl::_>))
                    , tag::partitioner*(*(boost::is_integral<mpl::_>))
                    , tag::partition_file*(*(boost::is_integral<mpl::_>))
                    , tag::depends*(*(boost::is_convertible<mpl::_,std::string>))
                    >
           >()
     );

                  
  // def("ConstruSimplex",ConstruSimplex<1,1,1>,return_value_policy<manage_new_object>()); 
               
   
  
  
       

}

