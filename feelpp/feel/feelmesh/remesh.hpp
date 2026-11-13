/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2020-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#pragma once

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/metric.hpp>
#include <feel/feelmesh/remesher.hpp>
#include <fmt/core.h>
namespace Feel
{
/**
 * @brief list of arguments type for remesh named parameters functions
 *
 * @tparam T
 */
template <typename T>
using args_remesh_type = NA::arguments<
    typename na::mesh::template required_as_t<std::shared_ptr<T>>,
    typename na::metric::template required_as_t<std::string>,
    typename na::required_elts::template required_as_t<std::vector<std::string>>,
    typename na::required_facets::template required_as_t<std::vector<std::string>>,
    typename na::parent::template required_as_t<std::shared_ptr<T>>,
    typename na::vm::template required_as_t<po::variables_map const&>,
    typename na::verbose::template required_as_t<int>>;


/**
 * @brief given a metric and required elements/facets, remesh a mesh
 *
 * @tparam MeshT type of the mesh to be remeshed
 * @param r mesh to be remeshed
 * @param metric_expr metric expression
 * @param req_elts required elements set
 * @param req_facets  required facets set
 * @param parent parent mesh (possibly nullptr)
 * @return std::shared_ptr<MeshT>
 */
template <typename MeshT>
std::tuple<std::shared_ptr<MeshT>, int>
remeshImpl( std::shared_ptr<MeshT> const& r,
        std::string const& metric_expr,
        std::vector<std::string> const& req_elts,
        std::vector<std::string> const& req_facets,
        std::shared_ptr<MeshT> const& parent,
        nl::json const& params )
{
    static int cpt = 0;
    using mesh_t = MeshT;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;

    auto setMetric = [&metric_expr, &req_facets]( auto const& P1h, auto& R ){ 
        if ( metric_expr.find("gradedls") != std::string::npos )
        {
            bool created = createGradedLevelsetMetric( P1h, metric_expr, 
                                                        [&R](auto const& met ){ R->setMetric(met); }, req_facets );
            if ( !created )
            {
                LOG( INFO ) << fmt::format( "[remesh] invalid metric expression={}\n", metric_expr ) << std::endl;
                throw std::invalid_argument( fmt::format( "[remesh] invalid metric expression={}\n", metric_expr ) );
            }
        }
        else 
        {
            LOG( INFO ) << fmt::format( "[remesh] metric={}\n", metric_expr ) << std::endl;
            R->setMetric( expr( P1h, expr( metric_expr ) ) );
        }
    };
    if ( r->worldComm().localSize() == 1 )
    {
        auto R = std::make_shared<Remesh<mesh_t>>( r, req_elts, req_facets, parent, std::string{}, params );
        auto P1h = Pch<1>( r );
        setMetric( P1h, R );
        auto m = R->execute();
        return std::tuple{ m, ++cpt };
    }
    else // parallel
    {
        std::string imeshParaName = fmt::format( "themesh{}_i.json", cpt );
        std::string omeshParaName = fmt::format( "themesh{}_o.json", cpt );
        std::string imeshParaPath = ( fs::path( Environment::appRepository() ) / imeshParaName ).string();
        std::string omeshParaPath = ( fs::path( Environment::appRepository() ) / imeshParaName ).string();
        int nPartition = r->worldComm().localSize();
        r->saveHDF5( imeshParaPath );
        r->worldComm().barrier();
        if ( r->worldComm().isMasterRank() )
        {
            auto meshSeqInput = [&]() {
                if constexpr ( mesh_t::is_order_dynamic )
                    return new mesh_t(RuntimeOrder{r->order()}, "", Environment::worldCommSeqPtr());
                else
                    return new mesh_t(Environment::worldCommSeqPtr());
            }();
            auto meshSeq = loadMesh( _mesh = meshSeqInput, _savehdf5 = 0,
                                     _filename = imeshParaPath,
                                     //_update=update_,  TODO test FACE_MINIMAL
                                     _straighten = false );
            auto P1h = Pch<1>( meshSeq );
            auto R = std::make_shared<Remesh<mesh_t>>( meshSeq, req_elts, req_facets, parent );
            setMetric( P1h, R );
            auto m = R->execute();
            using io_t = PartitionIO<mesh_t>;
            io_t io( omeshParaPath );
            io.write( partitionMesh( m, nPartition /*, partitionByRange, partconfig*/ ) );
        }
        r->worldComm().barrier();
        auto meshAdaptInput = [&]() {
            if constexpr ( mesh_t::is_order_dynamic )
                return new mesh_t(RuntimeOrder{r->order()}, "", r->worldCommPtr());
            else
                return new mesh_t(r->worldCommPtr());
        }();
        auto meshAdaptPara = loadMesh( _mesh = meshAdaptInput, _savehdf5 = 0,
                                       _filename = omeshParaPath /*, _update=update_*/ );
        return std::tuple{ meshAdaptPara, ++cpt };
    }
}

/**
 * @brief remesh with named arguments
 * @ingroup Mesh
 *
 * @tparam MeshType
 * @param _mesh mesh to be adapted
 * @param _required_elts list of markers of required elements
 * @param _required_facets list of markers of required facets
 * @param _parent parent mesh if relation is built
 *
 * @code {.cpp}
 * // adapt mesh with constant mesh size
 * auto adapted_mesh = remesh( _mesh=mesh, _metric="0.1" );
 * // adapt mesh with constant mesh size with required elements marked omega
 * auto adapted_mesh = remesh( _mesh=mesh, _metric="0.1", _required_elts=std::vector{{"omega"}} );
 * @endcode
 *
 * @return std::shared_ptr<MeshType>
 */
template <typename... Ts>
auto remesh( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& mesh = args.get( _mesh );
    auto&& metric = args.get( _metric );
    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype( mesh )>>>;

    auto args0 = std::move( args )
                     .add_default_arguments( NA::make_default_argument( _vm, Environment::vm() ),
                                             NA::make_default_argument( _parent, std::shared_ptr<mesh_type>{} ),
                                             NA::make_default_argument( _required_elts, std::vector<std::string>{} ),
                                             NA::make_default_argument( _required_facets, std::vector<std::string>{} ),
                                             NA::make_default_argument( _params, nl::json{} ) );

    po::variables_map const& vm = args0.get( _vm );
    auto required_elts = args0.get( _required_elts );
    auto required_facets = args0.get( _required_facets );
    auto&& parent = args0.get( _parent );
    auto&& params = args0.get( _params );
    return remeshImpl<mesh_type>( mesh, metric, required_elts, required_facets, parent, params );
}
} // namespace Feel
