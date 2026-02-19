/**
 * @file plot.cpp
 * @author christophe prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief 
 * @version 0.1
 * @date 2023-03-18
 * 
 * @copyright Copyright (c) 2023 Université de Strasbourg
 * 
 */

#include <feel/feelpython/pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <fmt/core.h>
#include <mpi4py/mpi4py.h>

#include <feel/feeldiscr/mesh.hpp>

namespace py = pybind11;

namespace Feel
{
/**
 * @brief exporter to plotly
 * 
 * it is required to renumber the points and elements of the mesh from 0 with continguous numbering
 * we have to hold a mapping between the Feel++ numbering and the plotly numbering
 * 
 * @tparam Dim topology dimension
 */
template<int Dim>
class ExporterPlotly
{
public:

    ExporterPlotly( std::shared_ptr<Mesh<Simplex<Dim, 1>>> const& mesh )
        :
        mesh_( mesh )
    {
        int i = 0;
        for( auto const& [id,p] : mesh_->points() )
            mapping_p_[id] = i++;
        LOG(INFO) << "number of points: " << i;

        p_.resize(i,Dim);

        i = 0;
        for ( auto const& e_ptr : elements(mesh_) )
        {
            auto const& e = boost::unwrap_ref( e_ptr );
            mapping_e_[e.id()] = i++;
        }
        LOG( INFO ) << "number of elements: " << i;
        c_.resize(i,Dim+1);

        for ( auto const& e_ptr : elements( mesh_ ) )
        {
            auto const& e = boost::unwrap_ref( e_ptr );
            auto nodesG = e.G();
            em_matrix_col_type<double> G( const_cast<double*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );

            for ( int i = 0; i < e.numVertices; i++ )
            {
                c_( mapping_e_[e.id()], i ) = mapping_p_[e.point( i ).id()];
                p_.row( mapping_p_[e.point( i ).id()] ) = G.col( i );
            }
        }
    }

    auto& getPoints() const { return p_; }
    auto& getConnectivity() const { return c_; }

#if 0
    void add( std::string const& name, std::shared_ptr<element_type<FunctionSpace<Mesh<Simplex<Dim, 1>>, bases<Lagrange<1, Scalar,Continuous>>>>> const& u )
    {
        
    }
#endif

private:
    //! mesh
    std::shared_ptr<Mesh<Simplex<Dim, 1>>> mesh_;
    //! prefix
    std::string prefix_;
    //! map of point ids
    std::unordered_map<int,int> mapping_p_,mapping_e_;
    //! map of point coordinates
    eigen_matrix_xx_col_type<double> p_;
    //! connectivity
    eigen_matrix_xx_col_type<int> c_;
    //!
    std::unordered_map<std::string, eigen_matrix_xx_col_type<double>> data_;
};
}
template <int Dim>
void _inst( py::module& m )
{
    using namespace Feel;
    using mesh_t = Mesh<Simplex<Dim, 1>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    constexpr int RealDim = Dim;
    using plotly_t = ExporterPlotly<Dim>;
    py::class_<plotly_t, std::shared_ptr<plotly_t>>( m, fmt::format( "ExporterPlotly{}", Dim ).c_str(), py::dynamic_attr() )
        .def( py::init<mesh_ptr_t const&>(), py::arg( "mesh" ), "Construct a new plotly mesh" )
        .def( "getPoints", &plotly_t::getPoints, "get Mesh coordinates", py::return_value_policy::reference )
        .def( "getConnectivity", &plotly_t::getConnectivity, "get Mesh connectivity", py::return_value_policy::reference )
        ;
}
PYBIND11_MODULE( _plot, m )
{
    using namespace Feel;

    if ( import_mpi4py() < 0 ) return;
    _inst<2>( m );
    _inst<3>( m );
}