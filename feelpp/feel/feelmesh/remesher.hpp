/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Feb 2020

 Copyright (C) 2020 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#pragma once
//#include <boost/container_hash/hash.hpp>
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelmesh/concatenate.hpp>
#include <feel/feelmesh/submeshdata.hpp>
#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
#include <mmg/libmmg.h>
#include <parmmg/libparmmg.h>
#endif
#include <variant>

namespace Feel
{
/**
 * @brief  a pairing function is a process to uniquely encode two natural numbers into a single natural number.
 *
 * @param a
 * @param b
 * @return constexpr int
 */
constexpr int cantor_pairing( int a, int b )
{
    return ( a + b ) * ( a + b + 1 ) / 2 + b;
}
/**
 * @brief compute cantor pairing recursively for number of natural numbers > 2
 *
 * @tparam T integrer type
 * @param a integer to be paired with
 * @param b integer to be paired with
 * @param ts remaining pack that will be treated recursively
 * @return constexpr int
 */
template <typename... T>
constexpr int cantor_pairing( int a, int b, T&&... ts )
{
    return cantor_pairing( cantor_pairing( a, b ), ts... );
}

enum class RemeshMode
{
    MMG = 1,
    PARMMG = 2
};

/**
 * @brief build command line options for remesh
 *
 * @param prefix
 * @return Feel::po::options_description
 */
Feel::po::options_description
remesh_options( std::string const& prefix );

#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )

/**
 * @brief set MMG options from command line
 * @ingroup Mesh
 * @tparam topoDim topological dimension
 * @tparam realDim real dimension
 *
 * @param prefix prefix to discriminate between different remeshers
 * @param sol mmg solution
 * @param mesh mmg mesh
 */
template <int topoDim, int realDim>
void setMMGOptions( nl::json const& j, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh, MMG5_pSol sol );

/**
 * @brief Class that handles remeshing in sequential using mmg and parallel using parmmg
 *
 * @tparam MeshType mesh type to be remeshed
 */
template <typename MeshType>
class Remesh
{
  public:
    using mmg_mesh_t = std::variant<MMG5_pMesh, PMMG_pParMesh>;
    using mesh_t = MeshType;
    using mesh_ptrtype = std::shared_ptr<MeshType>;
    using mesh_ptr_t = mesh_ptrtype;
    using scalar_metric_t = typename Pch_type<mesh_t, 1>::element_type;
    using marker_type = typename mesh_t::element_type::marker_type;
    typedef SubMeshData<index_type> smd_type;
    typedef std::shared_ptr<smd_type> smd_ptrtype;

    Remesh()
        : Remesh( nullptr )
    {
    }
    Remesh( std::shared_ptr<MeshType> const& mesh )
        : Remesh( mesh, std::vector<std::string>{}, std::vector<std::string>{} ) {}
    Remesh( std::shared_ptr<MeshType> const& mesh, nl::json const& params )
        : Remesh( mesh, std::vector<std::string>{}, std::vector<std::string>{}, nullptr, std::string{}, params ) {}
    Remesh( std::shared_ptr<MeshType> const& mesh,
            boost::any const& required_element_markers )
        : Remesh( mesh, required_element_markers, std::vector<std::string>{} ) {}
    Remesh( std::shared_ptr<MeshType> const& mesh,
            boost::any const& required_element_markers,
            boost::any const& required_facet_markers,
            std::shared_ptr<MeshType> const& parent = nullptr,
            std::string const& prefix = {},
            nl::json const& params = {} );
    Remesh( Remesh const& r ) = default;
    Remesh( Remesh&& r ) = default;

    ~Remesh()
    {
        if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Free_all( MMG5_ARG_start,
                                MMG5_ARG_ppMesh, &std::get<MMG5_pMesh>( M_mmg_mesh ), MMG5_ARG_ppMet, &M_mmg_sol,
                                MMG5_ARG_end );
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
                MMGS_Free_all( MMG5_ARG_start,
                               MMG5_ARG_ppMesh, &std::get<MMG5_pMesh>( M_mmg_mesh ), MMG5_ARG_ppMet, &M_mmg_sol,
                               MMG5_ARG_end );
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
                MMG2D_Free_all( MMG5_ARG_start,
                                MMG5_ARG_ppMesh, &std::get<MMG5_pMesh>( M_mmg_mesh ), MMG5_ARG_ppMet, &M_mmg_sol,
                                MMG5_ARG_end );
        }
        else if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
        {
            PMMG_Free_all( PMMG_ARG_start,
                           PMMG_ARG_ppParMesh, &std::get<PMMG_pParMesh>( M_mmg_mesh ),
                           PMMG_ARG_end );
        }
    }

    /**
     * set scalar metric
     */
    void setMetric( scalar_metric_t const& );
    void setMetricLs( scalar_metric_t const& );

    /**
     * execute remesh task
     */
    mesh_ptrtype execute( bool run = true )
    {
        if ( run )
        {
            if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
            {
                auto mesh = std::get<MMG5_pMesh>( M_mmg_mesh );
                auto fname = soption( _name = "remesh.save", _prefix = M_prefix );
                auto save = [&fname, &mesh, this]( auto fn1, auto fn2 )
                {
                    if ( !fname.empty() )
                    {
                        auto mesh_name = fmt::format( "{}.mesh", fname );
                        int ier = fn1( mesh, mesh_name.c_str() );
                        if ( ier == MMG5_STRONGFAILURE )
                            throw std::logic_error( "Unable to save the mesh to mesh format." );
                        auto name_sol = fmt::format( "{}.sol", fname );
                        ier = fn2( mesh, this->M_mmg_sol, name_sol.c_str() );
                        if ( ier == MMG5_STRONGFAILURE )
                            throw std::logic_error( "Unable to save the mesh to sol format." );
                    }
                };

                if constexpr ( dimension_v<MeshType> == 3 )
                {
                    save( MMG3D_saveMesh, MMG3D_saveSol );
                    MMG3D_mmg3dlib( mesh, M_mmg_sol );
                }
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
                    MMGS_mmgslib( mesh, M_mmg_sol );
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
                {
                    save( MMG2D_saveMesh, MMG2D_saveSol );
                    MMG2D_mmg2dlib( mesh, M_mmg_sol );
                }
            }
            else if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
            {
                auto mesh = std::get<PMMG_pParMesh>( M_mmg_mesh );
                if ( M_mesh->worldCommPtr()->localSize() > 1 )
                {
                    int ier = PMMG_parmmglib_distributed( mesh );
                    if ( ier == PMMG_STRONGFAILURE )
                        throw std::logic_error( "Unable to adapt mesh in distributed mode." );
                }
                else
                {
                    int ier = PMMG_parmmglib_centralized( mesh );
                    if ( ier == PMMG_STRONGFAILURE )
                        throw std::logic_error( "Unable to adapt mesh in centralized mode." );
                }
            }
        }
        auto r = this->mmg2Mesh();
        return r;
    }

    /**
     * transform a Feel++ \p mesh into an Mmg mesh
     */
    mmg_mesh_t mesh2Mmg( std::shared_ptr<MeshType> const& m_in );
    mmg_mesh_t mesh2Mmg() { return mesh2Mmg( M_mesh ); }

    /**
     * convert a Mmg mesh into a Feel++ \p mesh
     */
    mesh_ptrtype mmg2Mesh( mmg_mesh_t const& m_in );
    mesh_ptrtype mmg2Mesh() { return mmg2Mesh( M_mmg_mesh ); }

    /**
     * @brief Set the Keep Relation For Required Entities
     *
     * @param keep 0 do not keep relation, 1 keep
     */
    void setKeepRelationForRequiredEntities( bool keep ) {
        keep_relation_required_ = keep;
    }

    /**
     * @brief Set the Preserve Floating Interface
     *
     * @param preserve 0 or 1
     */
    void setPreserveFloatingInterface( bool preserve )
    {
        if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
        {
            auto mesh = std::get<MMG5_pMesh>( M_mmg_mesh );
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                MMG3D_Set_iparameter( mesh, M_mmg_sol, MMG3D_IPARAM_opnbdy, preserve );
            }
            if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
            {
                MMG2D_Set_iparameter( mesh, M_mmg_sol, MMG2D_IPARAM_opnbdy, preserve );
            }
        }
        else
        {}
    }

    /**
     * Insert object
    */
    mesh_ptrtype mmgls2Mesh( mmg_mesh_t const& m_in);
    std::shared_ptr<MeshType> mesh_from_ls( scalar_metric_t const&, std::vector<int>, std::vector<int>);
    std::shared_ptr<MeshType> mesh_from_ls_met( scalar_metric_t const&, scalar_metric_t const&, std::vector<int>, std::vector<int>);

  private:
    void setParameters();
    void setCommunicatorAPI();
    using neighbor_ent_t = std::map<int, std::pair<int, int>>;
    using local_face_index_t = std::vector<std::vector<int>>;
    using comm_t = std::tuple<neighbor_ent_t, local_face_index_t, std::map<int, int>>;
    comm_t getCommunicatorAPI();

  private:
    std::shared_ptr<MeshType> M_mesh, M_parent_mesh;
    std::string M_prefix;
    nl::json M_params;
    RemeshMode M_mode = RemeshMode::PARMMG;
    mmg_mesh_t M_mmg_mesh;
    MMG5_pSol M_mmg_sol;
    MMG5_pSol M_mmg_met;

    std::unordered_map<int, int> pt_id;
    std::unordered_map<int, std::pair<int, int>> face_id;

    bool keep_relation_required_;
    boost::any M_required_element_markers;
    boost::any M_required_facet_markers;
    std::mutex mutex_;
    // sub mesh data
    smd_ptrtype M_smd;
    std::map<int, std::pair<int, int>> M_id2tag, M_id2tag_face;

    
    std::map<int,std::tuple<ElementsType,int,marker_type>> M_mapMmgFragmentIdToMeshFragementDesc;
};

template <typename MeshType>
Remesh<MeshType> remesher( std::shared_ptr<MeshType> const& m,
                           boost::any required_element_markers = std::vector<int>{},
                           boost::any required_facet_markers = std::vector<int>{},
                           std::shared_ptr<MeshType> const& parent = {},
                           std::string const& prefix = {},
                           nl::json const& params = {} )
{
    return Remesh<MeshType>{ m, required_element_markers, required_facet_markers, parent, prefix, params };
}

template <typename MeshType>
Remesh<MeshType>::Remesh( std::shared_ptr<MeshType> const& mesh,
                          boost::any const& required_element_markers,
                          boost::any const& required_facet_markers,
                          std::shared_ptr<MeshType> const& parent,
                          std::string const& prefix,
                          nl::json const& params )
    : M_mesh( mesh ),
      M_parent_mesh( parent ),
      M_prefix( prefix ),
      M_params( params ),
      M_mode( RemeshMode::MMG ),
      M_mmg_mesh(),
      M_mmg_sol( nullptr ),
      M_mmg_met( nullptr ),
      keep_relation_required_( false ),
      M_required_element_markers( required_element_markers ),
      M_required_facet_markers( required_facet_markers ),
      M_smd( new smd_type( M_parent_mesh ? M_parent_mesh : M_mesh ) )
{
    if ( params.contains("remesh") )
    {
        if ( params["remesh"].contains("required" ) )
        {
            if ( M_params["remesh"]["required"].contains("elements") )
            {
                M_required_element_markers = M_params["remesh"]["required"]["elements"].template get<std::vector<std::string>>();
            }
            if ( M_params["remesh"]["required"].contains( "facets" ) )
            {
                M_required_facet_markers = M_params["remesh"]["required"]["facets"].template get<std::vector<std::string>>();
            }
            if ( M_params["remesh"]["required"].contains( "keep_relation" ) )
                keep_relation_required_ = M_params["remesh"]["required"]["keep_relation"].template get<bool>();
            if ( !M_params["remesh"]["required"]["facets"].template get<std::vector<std::string>>().empty() &&
                 keep_relation_required_ )
            {
                // this is necessary otherwise the hack to keep the relation and require facets will fail
                M_params["remesh"]["opnbdy"] = 1;
            }
        }
    }
    if ( M_parent_mesh && M_parent_mesh->isRelatedTo( M_mesh ) == false )
        throw std::logic_error( "invalid parent mesh" );
    else if ( M_parent_mesh && M_parent_mesh->isRelatedTo( M_mesh ) )
        VLOG(2) << fmt::format("-- remesh use parent mesh") << std::endl;

    if ( mesh->worldCommPtr()->localSize() == 1 && M_mode == RemeshMode::MMG )
    {
        MMG5_pMesh m = nullptr;
        if constexpr ( dimension_v<MeshType> == 3 )
            MMG3D_Init_mesh( MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &m, MMG5_ARG_ppLs,&M_mmg_met, MMG5_ARG_ppMet, &M_mmg_sol,
                             MMG5_ARG_end );
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
            MMGS_Init_mesh( MMG5_ARG_start,
                            MMG5_ARG_ppMesh, &m, MMG5_ARG_ppLs,&M_mmg_met, MMG5_ARG_ppMet, &M_mmg_sol,
                            MMG5_ARG_end );
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
            MMG2D_Init_mesh( MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &m, MMG5_ARG_ppLs,&M_mmg_met, MMG5_ARG_ppMet, &M_mmg_sol,
                             MMG5_ARG_end );
        M_mmg_mesh = m;
        setParameters();
    }
    else
    {

        PMMG_pParMesh m = nullptr;
        PMMG_Init_parMesh( PMMG_ARG_start,
                           PMMG_ARG_ppParMesh, &m,
                           PMMG_ARG_pMesh, PMMG_ARG_pMet,
                           PMMG_ARG_dim, mesh->realDimension(), PMMG_ARG_MPIComm, static_cast<MPI_Comm>( mesh->worldCommPtr()->comm() ),
                           PMMG_ARG_end );

        M_mmg_mesh = m;
        PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_verbose, 0 );
        PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_mmgVerbose, 1 );
        PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_debug, 5 );
        //PMMG_Set_dparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_DPARAM_hmin, 1e-2 );
    }

    this->mesh2Mmg();
    this->setParameters();
    this->setCommunicatorAPI();
}

template <typename MeshType>
void Remesh<MeshType>::setMetric( scalar_metric_t const& m )
{
    if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
    {
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_sol,
                                    MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
        {
            if ( MMGS_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_sol,
                                   MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
        {
            if ( MMG2D_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_sol,
                                    MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }

        for ( auto const& welt : M_mesh->elements() )
        {
            auto const& [key, elt] = boost::unwrap_ref( welt );
            for ( auto const& ldof : m.functionSpace()->dof()->localDof( elt.id() ) )
            {
                size_type index = ldof.second.index();
                uint16_type local_dof = ldof.first.localDof();
                auto s = m( index );
                int pos = pt_id[elt.point( local_dof ).id()];
                if constexpr ( dimension_v<MeshType> == 3 )
                {
                    if ( MMG3D_Set_scalarSol( M_mmg_sol, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
                {
                    if ( MMGS_Set_scalarSol( M_mmg_sol, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
                {
                    if ( MMG2D_Set_scalarSol( M_mmg_sol, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
            }
        }
    }
    if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
    {
        auto mesh = std::get<PMMG_pParMesh>( M_mmg_mesh );
        LOG( INFO ) << fmt::format( " - setting metric  size: {}...", m.nDof() );
        if ( PMMG_Set_metSize( mesh, MMG5_Vertex, m.nLocalDof(), MMG5_Scalar ) != 1 )
        {
            fmt::print( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                        "Unable to allocate the metric array.\n" );
            throw std::logic_error( "Unable to set metric size" );
        }
        for ( auto const& welt : M_mesh->elements() )
        {
            auto const& [key, elt] = boost::unwrap_ref( welt );
            for ( auto const& ldof : m.functionSpace()->dof()->localDof( elt.id() ) )
            {
                size_type index = ldof.second.index();
                uint16_type local_dof = ldof.first.localDof();
                auto s = m( index );
                int pos = pt_id[elt.point( local_dof ).id()];
                if constexpr ( dimension_v<MeshType> == 3 )
                {
                    if ( PMMG_Set_scalarMet( mesh, s, pos ) != 1 )
                    {
                        fmt::print( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                    "Unable to set metric {} at pos {}.\n", s, pos );
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
            }
        }
        double* sol = new double[m.nLocalDof()];
        for ( int k = 0; k < m.nLocalDof(); k++ )
        {
            /* Vertex by vertex */
            if ( PMMG_Get_scalarMet( mesh, &sol[k] ) != 1 )
            {
                LOG( ERROR ) << fmt::format( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                             "Unable to get metrics {} \n", k );
            }
            if ( sol[k] <= 0.0 )
            {
                LOG( ERROR ) << fmt::format( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                             "Invalid metrics {} at {} \n", sol[k], k );
            }
        }
        delete[] sol;
    }
}


template <typename MeshType>
void Remesh<MeshType>::setMetricLs( scalar_metric_t const& m )
{
    if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
    {
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_met,
                                    MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
        {
            if ( MMGS_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_met,
                                   MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
        {
            if ( MMG2D_Set_solSize( std::get<MMG5_pMesh>( M_mmg_mesh ), M_mmg_met,
                                    MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar ) != 1 )
            {
                throw std::logic_error( "Unable to allocate the metric array." );
            }
        }

        for ( auto const& welt : M_mesh->elements() )
        {
            auto const& [key, elt] = boost::unwrap_ref( welt );
            for ( auto const& ldof : m.functionSpace()->dof()->localDof( elt.id() ) )
            {
                size_type index = ldof.second.index();
                uint16_type local_dof = ldof.first.localDof();
                auto s = m( index );
                int pos = pt_id[elt.point( local_dof ).id()];
                if constexpr ( dimension_v<MeshType> == 3 )
                {
                    if ( MMG3D_Set_scalarSol( M_mmg_met, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
                {
                    if ( MMGS_Set_scalarSol( M_mmg_met, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
                else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
                {
                    if ( MMG2D_Set_scalarSol( M_mmg_met, s, pos ) != 1 )
                    {
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
            }
        }
    }
    if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
    {
        auto mesh = std::get<PMMG_pParMesh>( M_mmg_mesh );
        LOG( INFO ) << fmt::format( " - setting metric  size: {}...", m.nDof() );
        if ( PMMG_Set_metSize( mesh, MMG5_Vertex, m.nLocalDof(), MMG5_Scalar ) != 1 )
        {
            fmt::print( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                        "Unable to allocate the metric array.\n" );
            throw std::logic_error( "Unable to set metric size" );
        }
        for ( auto const& welt : M_mesh->elements() )
        {
            auto const& [key, elt] = boost::unwrap_ref( welt );
            for ( auto const& ldof : m.functionSpace()->dof()->localDof( elt.id() ) )
            {
                size_type index = ldof.second.index();
                uint16_type local_dof = ldof.first.localDof();
                auto s = m( index );
                int pos = pt_id[elt.point( local_dof ).id()];
                if constexpr ( dimension_v<MeshType> == 3 )
                {
                    if ( PMMG_Set_scalarMet( mesh, s, pos ) != 1 )
                    {
                        fmt::print( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                    "Unable to set metric {} at pos {}.\n", s, pos );
                        throw std::logic_error( "Unable to set metric" );
                    }
                }
            }
        }
        double* sol = new double[m.nLocalDof()];
        for ( int k = 0; k < m.nLocalDof(); k++ )
        {
            /* Vertex by vertex */
            if ( PMMG_Get_scalarMet( mesh, &sol[k] ) != 1 )
            {
                LOG( ERROR ) << fmt::format( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                             "Unable to get metrics {} \n", k );
            }
            if ( sol[k] <= 0.0 )
            {
                LOG( ERROR ) << fmt::format( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                                             "Invalid metrics {} at {} \n", sol[k], k );
            }
        }
        delete[] sol;
    }
}


template <typename MeshType>
typename Remesh<MeshType>::mmg_mesh_t
Remesh<MeshType>::mesh2Mmg( std::shared_ptr<MeshType> const& m_in )
{
    int nVertices = nelements( points( m_in ) );
    int nTetrahedra = ( dimension_v<MeshType> == 3 ) ? nelements( elements( m_in ) ) : 0;
    int nPrisms = 0;
    auto boundary_and_marked_faces = concatenate( boundaryfaces( m_in ), markedfaces( m_in, M_required_facet_markers ) );
    int nTriangles = ( dimension_v<MeshType> == 3 ) ? m_in->numFaces() : nelements( elements( m_in ) ); //( dimension_v<MeshType> == 3 ) ? ( ( m_in->worldCommPtr()->localSize() > 1 ) ? nelements( boundary_and_marked_faces ) + nelements( interprocessfaces( m_in ) ) : m_in->numFaces() ) : nelements( elements( m_in ) );
    int nQuadrilaterals = 0;
    int nEdges = ( dimension_v<MeshType> == 2 ) ? m_in->numFaces() : 0;


    // create mmg marker fragmentation and the inverse mapping of mesh fragmentation
    // we only propagte the marker type 1
    M_mapMmgFragmentIdToMeshFragementDesc.clear();
    std::map<ElementsType,std::map<marker_type,int>> mapMarkerToMmgFragmentId;
    int mmgFragmentId = 0;
    for ( auto const& [et,fragMapping] : m_in->meshFragmentationByMarkerByEntity() )
    {
        auto & mapUp = mapMarkerToMmgFragmentId[et];
        for ( auto const& [fragId,fragMarker] : fragMapping )
        {
            mapUp.emplace(fragMarker,mmgFragmentId );
            M_mapMmgFragmentIdToMeshFragementDesc.emplace( mmgFragmentId, std::make_tuple(et,fragId,fragMarker) );
            ++mmgFragmentId;
        }
    }

    if ( std::holds_alternative<MMG5_pMesh>( M_mmg_mesh ) )
    {
        auto mesh = std::get<MMG5_pMesh>( M_mmg_mesh );
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Set_meshSize( mesh, nVertices, nTetrahedra, nPrisms, nTriangles,
                                     nQuadrilaterals, nEdges ) != 1 )
            {
                throw std::logic_error( "Error in MMG3D_Set_meshSize" );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
        {
            if ( MMGS_Set_meshSize( mesh, nVertices, nTriangles, nEdges ) != 1 )
            {
                throw std::logic_error( "Error in MMGS_Set_meshSize" );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
        {
            if ( MMG2D_Set_meshSize( mesh, nVertices, nTriangles, nQuadrilaterals, nEdges ) != 1 )
            {
                throw std::logic_error( "Error in MMG2D_Set_meshSize" );
            }
        }

        pt_id.reserve( nVertices );
        int k = 1;
        auto const& mapMarkerToMmgFragmentId_points = mapMarkerToMmgFragmentId.at( ElementsType::MESH_POINTS );
        for ( auto const& wpt : points( m_in ) )
        {
            auto const& pt = boost::unwrap_ref( wpt );
            int ptMarker = pt.hasMarker() ? mapMarkerToMmgFragmentId_points.at( pt.marker() ) : -1;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Set_vertex( std::get<MMG5_pMesh>( M_mmg_mesh ), pt( 0 ), pt( 1 ), pt( 2 ), ptMarker, k ) != 1 )
                {
                    throw std::logic_error( "Error in MMG3D_Set_vertex" );
                }
                if ( pt.hasMarker() && MMG3D_Set_requiredVertex( std::get<MMG5_pMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in MMG3D_Set_requiredVertex" );
                }
                pt_id[pt.id()] = k++;
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
            {
                if ( MMGS_Set_vertex( std::get<MMG5_pMesh>( M_mmg_mesh ), pt( 0 ), pt( 1 ), pt( 2 ), ptMarker, k ) != 1 )
                {
                    throw std::logic_error( "Error in MMGS_Set_vertex" );
                }
                if ( pt.hasMarker() && MMGS_Set_requiredVertex( std::get<MMG5_pMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in MMGS_Set_requiredVertex" );
                }
                pt_id[pt.id()] = k++;
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
            {
                if ( MMG2D_Set_vertex( std::get<MMG5_pMesh>( M_mmg_mesh ), pt( 0 ), pt( 1 ), ptMarker, k ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_vertex" );
                }
                if ( pt.hasMarker() && MMG2D_Set_requiredVertex( std::get<MMG5_pMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_requiredVertex" );
                }
                pt_id[pt.id()] = k++;
            }
        }
        auto required_element_ids = m_in->markersId( M_required_element_markers );
        int req_elts = nelements( markedelements( m_in, M_required_element_markers ) );

        int next_free_req = 1;
        int next_free_nreq = req_elts + 1;
        k = 1;
        bool required = false;
        auto const& mapMarkerToMmgFragmentId_elements = mapMarkerToMmgFragmentId.at( ElementsType::MESH_ELEMENTS );
        for ( auto const& welt : m_in->elements() )
        {
            auto const& [key, elt] = boost::unwrap_ref( welt );
            int eltMarker = elt.hasMarker() ? mapMarkerToMmgFragmentId_elements.at( elt.marker() ) : -1;
            required = elt.hasMarker()? elt.marker().hasOneOf( required_element_ids ) : false;
            int& id_elt = ( required ) ? next_free_req : next_free_nreq;
            M_id2tag[id_elt].first = eltMarker;
            int lab_or_id = ( keep_relation_required_ && required ) ? id_elt : eltMarker;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Set_tetrahedron( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                            pt_id[elt.point( 0 ).id()],
                                            pt_id[elt.point( 1 ).id()],
                                            pt_id[elt.point( 2 ).id()],
                                            pt_id[elt.point( 3 ).id()],
                                            lab_or_id, id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG3D_Set_tetrahedron" );
                }
                if ( required &&
                     MMG3D_Set_requiredTetrahedron( std::get<MMG5_pMesh>( M_mmg_mesh ), id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG3D_Set_requiredTetrahedron" );
                }
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
            {
                if ( MMGS_Set_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                        pt_id[elt.point( 0 ).id()],
                                        pt_id[elt.point( 1 ).id()],
                                        pt_id[elt.point( 2 ).id()],
                                        lab_or_id, id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMGS_Set_triangle" );
                }
                if ( required &&
                     MMGS_Set_requiredTriangle( std::get<MMG5_pMesh>( M_mmg_mesh ), id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMGS_Set_requiredTriangle" );
                }
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
            {
                if ( MMG2D_Set_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         pt_id[elt.point( 0 ).id()],
                                         pt_id[elt.point( 1 ).id()],
                                         pt_id[elt.point( 2 ).id()],
                                         lab_or_id, id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_triangle" );
                }
                if ( required &&
                     MMG2D_Set_requiredTriangle( std::get<MMG5_pMesh>( M_mmg_mesh ), id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_requiredTriangle" );
                }
            }
            if ( keep_relation_required_ && required )
            {
                if ( M_parent_mesh )
                {
                    // id_in_parent_mesh
                    if ( M_mesh->isSubMeshFrom( M_parent_mesh ) )
                    {
                        // get parent id from submesh id
                        M_id2tag[id_elt].second = M_mesh->subMeshToMesh( M_parent_mesh, elt.id() );
                    }
                    else if ( M_parent_mesh->isSubMeshFrom( M_mesh ) )
                    {
                        // get submesh id from parent id
                        M_id2tag[id_elt].second = M_parent_mesh->meshToSubMesh( M_mesh, elt.id() );
                    }
                }
                else
                {
                    M_id2tag[id_elt].second = elt.id();
                }

                //std::cout << fmt::format( "[feelpp->mmg] added required element {}, key:{}, tag: {}", id_elt, key, elt.markerOr( 0 ).value() ) << std::endl;
            }
            k++;
            id_elt++;
        }

        auto required_facet_ids = m_in->markersId( M_required_facet_markers );
        int f_req_elts = nelements( markedfaces( m_in, M_required_facet_markers ) );

        int f_next_free_req = 1;
        int f_next_free_nreq = f_req_elts + 1;
        k = 1;
        required = false;
        auto const& mapMarkerToMmgFragmentId_faces = mapMarkerToMmgFragmentId.at( ElementsType::MESH_FACES );
        for ( auto const& wface : m_in->faces() )
        {
            auto const& [key, face] = boost::unwrap_ref( wface );
            int faceMarker = face.hasMarker() ? mapMarkerToMmgFragmentId_faces.at( face.marker() ) : -1;
            required = face.hasMarker()? face.marker().hasOneOf( required_facet_ids ) : false;
            int& id_elt = required ? f_next_free_req : f_next_free_nreq;
            M_id2tag_face[id_elt].first = faceMarker;
            int lab_or_id = ( keep_relation_required_ && required ) ? id_elt : faceMarker;
            if constexpr ( dimension_v<MeshType> == 3 )
            {

                if ( MMG3D_Set_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         pt_id[face.point( 0 ).id()],
                                         pt_id[face.point( 1 ).id()],
                                         pt_id[face.point( 2 ).id()],
                                         lab_or_id, id_elt ) != 1 )
                {
                    using namespace std::string_literals;
                    throw std::logic_error( "Error in MMG3D_Set_triangle "s + std::to_string( face.id() ) + " " + std::to_string( faceMarker ) );
                }

                if ( required &&
                     MMG3D_Set_requiredTriangle( std::get<MMG5_pMesh>( M_mmg_mesh ), id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_requiredTriangle" );
                }
                face_id[face.id()] = std::pair{ id_elt, id_elt };
            }
            else if constexpr ( dimension_v<MeshType> == 2 )
            {

                if ( MMG2D_Set_edge( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                     pt_id[face.point( 0 ).id()],
                                     pt_id[face.point( 1 ).id()],
                                     lab_or_id, id_elt ) != 1 )
                {
                    using namespace std::string_literals;
                    throw std::logic_error( "Error in MMG2D_Set_edge "s + std::to_string( face.id() ) + " " + std::to_string( faceMarker ) );
                }
                if ( required && MMG2D_Set_requiredEdge( std::get<MMG5_pMesh>( M_mmg_mesh ), id_elt ) != 1 )
                {
                    throw std::logic_error( "Error in MMG2D_Set_requiredEdge" );
                }
            }
            if ( keep_relation_required_ && required )
            {
                if ( M_parent_mesh )
                {
                    // id_in_parent_mesh
                    M_id2tag_face[id_elt].second = M_mesh->subMeshToMesh( M_parent_mesh, face.id() );
                }
                else
                {
                    M_id2tag_face[id_elt].second = face.id();
                }

                //std::cout << fmt::format( "[feelpp->mmg] added required element {}, key:{}, tag: {}", id_elt, key, elt.markerOr( 0 ).value() ) << std::endl;
            }
#if 0
            if ( required )
            {
                 auto key = std::tuple { pt_id[face.point( 0 ).id()], pt_id[face.point( 1 ).id()], pt_id[face.point( 2 ).id()] };
                if ( M_parent_mesh )
                {
                    int id_in_parent_mesh = M_mesh->subMeshToMesh( face.id() );
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( id_elt, id_in_parent_mesh ) );
                }
                else
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( id_elt, face.id() ) );
                }
                std::cout << fmt::format( "[feelpp->mmg] added required facet {}, key:{}, tag: {}", id_elt, key, faceMarker ) << std::endl;

            }
#endif
            // update facet index
            k++;
            id_elt++;
        }
    }
    else if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
    {
        LOG( INFO ) << fmt::format( "nvertices: {}, nTriangles:{} b:{} i:{} t:{}, ntetra: {}", nVertices, nTriangles, nelements( boundaryfaces( M_mesh ) ),
                                    nelements( interprocessfaces( M_mesh ) ), nelements( boundaryfaces( M_mesh ) ) + nelements( interprocessfaces( M_mesh ) ), nTetrahedra );
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( PMMG_Set_meshSize( std::get<PMMG_pParMesh>( M_mmg_mesh ), nVertices, nTetrahedra, nPrisms, nTriangles,
                                    nQuadrilaterals, nEdges ) != 1 )
            {
                throw std::logic_error( "Error in PMMG_Set_meshSize" );
            }
        }
        pt_id.reserve( nVertices );
        int k = 1;
        for ( auto const& wpt : points( m_in ) )
        {
            auto const& pt = boost::unwrap_ref( wpt );
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                int ptMarker = pt.hasMarker() ? pt.marker().value() : 0;
                if ( PMMG_Set_vertex( std::get<PMMG_pParMesh>( M_mmg_mesh ), pt( 0 ), pt( 1 ), pt( 2 ), ptMarker, k ) != 1 )
                {
                    throw std::logic_error( "Error in PMMG_Set_vertex" );
                }
                if ( pt.hasMarker() && PMMG_Set_requiredVertex( std::get<PMMG_pParMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in PMMG_Set_requiredVertex" );
                }
                pt_id[pt.id()] = k;
                k++;
            }
        }
        assert( k == nVertices + 1 );
        auto required_element_ids = m_in->markersId( M_required_element_markers );
        k = 1;
        for ( auto const& welt : elements( m_in ) )
        {
            auto const& elt = boost::unwrap_ref( welt );
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                int eltMarker = elt.hasMarker() ? elt.marker().value() : 0;
                if ( PMMG_Set_tetrahedron( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                           pt_id[elt.point( 0 ).id()],
                                           pt_id[elt.point( 1 ).id()],
                                           pt_id[elt.point( 2 ).id()],
                                           pt_id[elt.point( 3 ).id()],
                                           eltMarker, k ) != 1 )
                {
                    throw std::logic_error( "Error in PMMG_Set_tetrahedron" );
                }
                if ( required_element_ids.count( eltMarker ) &&
                     PMMG_Set_requiredTetrahedron( std::get<PMMG_pParMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in PMMG_Set_requiredTetrahedron" );
                }
                k++;
            }
        }
        assert( k == nTetrahedra + 1 );
        auto required_facet_ids = m_in->markersId( M_required_facet_markers );

        k = 1;
        for ( auto const& wface : boundary_and_marked_faces )
        {
            auto const& face = boost::unwrap_ref( wface );
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                int faceMarker = face.hasMarker() ? face.marker().value() : 0;
                if ( PMMG_Set_triangle( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                        pt_id[face.point( 0 ).id()],
                                        pt_id[face.point( 1 ).id()],
                                        pt_id[face.point( 2 ).id()],
                                        faceMarker, k ) != 1 )
                {
                    using namespace std::string_literals;
                    throw std::logic_error( "Error in PMMG_Set_triangle "s + std::to_string( face.id() ) + " " + std::to_string( faceMarker ) );
                }
                if ( required_facet_ids.count( faceMarker ) &&
                     PMMG_Set_requiredTriangle( std::get<PMMG_pParMesh>( M_mmg_mesh ), k ) != 1 )
                {
                    throw std::logic_error( "Error in PMMG_Set_requiredTriangle" );
                }

                auto s = cantor_pairing( face.point( 0 ).id(), face.point( 1 ).id(), face.point( 2 ).id() );
                DVLOG( 3 ) << fmt::format( "face: {} pt: {} {} {} renumber: {} pt: {} {} {}",
                                           face.id(), face.point( 0 ).id(), face.point( 1 ).id(), face.point( 2 ).id(),
                                           k, pt_id[face.point( 0 ).id()], pt_id[face.point( 1 ).id()], pt_id[face.point( 2 ).id()] );
                face_id[face.id()] = std::make_pair( k, s );
                if ( M_parent_mesh )
                {
                    int id_in_parent_mesh = M_parent_mesh->subMeshToMesh( face.id() );
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, id_in_parent_mesh ) );
                }
                else
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, face.id() ) );
                }

                k++;
            }
        }
        std::vector<int> pointIdInElt( 4 ), pointIdInFace( 3 );
        for ( auto const& wface : interprocessfaces( M_mesh ) )
        {
            auto const& face = boost::unwrap_ref( wface );
            auto const& elt0 = face.element0();
            int idInElt0 = face.idInElement0();
            auto const& elt1 = face.element1();
            int idInElt1 = face.idInElement1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& elt = ( elt0isGhost ) ? elt1 : elt0;
            int faceIdInElt = ( elt0isGhost ) ? idInElt1 : idInElt0;
            for ( uint16_type f = 0; f < elt.numVertices; ++f )
                pointIdInElt[f] = elt.point( f ).id();
            for ( uint16_type f = 0; f < 3; ++f )
            {
                pointIdInFace[f] = pointIdInElt[elt.fToP( faceIdInElt, f )];
            }
            DVLOG( 3 ) << fmt::format( "interprocess face: {} {} {}, vids:{}, {} {} {}", face.id(), faceIdInElt, k, pointIdInFace,
                                       pt_id[pointIdInFace[0]], pt_id[pointIdInFace[1]], pt_id[pointIdInFace[2]] );
            int faceMarker = face.hasMarker() ? face.marker().value() : 0;
            if ( PMMG_Set_triangle( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                    pt_id[pointIdInFace[0]],
                                    pt_id[pointIdInFace[1]],
                                    pt_id[pointIdInFace[2]],
                                    faceMarker, k ) != 1 )
            {
                using namespace std::string_literals;
                throw std::logic_error( "Error in PMMG_Set_triangle "s + std::to_string( face.id() ) + " " + std::to_string( faceMarker ) );
            }
            face_id[face.id()] = std::make_pair( k, 0 );
            k++;
        }
        assert( k == nTriangles + 1 );
    }
    return M_mmg_mesh;
}

template <typename MeshType>
std::shared_ptr<MeshType>
Remesh<MeshType>::mmg2Mesh( mmg_mesh_t const& mesh )
{
    int ier;

    int nVertices = 0;
    int nTetrahedra = 0;
    int nTriangles = 0;
    int nEdges = 0;

    std::shared_ptr<MeshType> out = std::make_shared<MeshType>(M_mesh->worldCommPtr());

    if ( std::holds_alternative<MMG5_pMesh>( mesh ) )
    {
        if ( MMG3D_Get_meshSize( std::get<MMG5_pMesh>( mesh ), &nVertices, &nTetrahedra, NULL, &nTriangles, NULL,
                                 &nEdges ) != 1 )
        {
            ier = MMG5_STRONGFAILURE;
        }
        int corner, required, tag = 0;
        node_type n( mesh_t::nRealDim );
        for ( int k = 1; k <= nVertices; k++ )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                       &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format("Unable to get mesh vertex {} ", k) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
            {
                if ( MMGS_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                      &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format( "Unable to get mesh vertex {} ", k ) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            if constexpr ( dimension_v<MeshType> == 2 )
            {

                if ( MMG2D_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ),
                                       &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format( "Unable to get mesh vertex {} ", k ) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            using point_type = typename mesh_t ::point_type;
            point_type pt( k, n );
            pt.setProcessIdInPartition( 0 );
            pt.setProcessId( 0 );
            if ( tag >= 0 )
            {
                auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( tag );
                if ( et == ElementsType::MESH_POINTS ) // only marker come from point marker
                    pt.setMarker( marker );
            }
            out->addPoint( pt );
        }

        if constexpr ( dimension_v<MeshType> == 3 )
        {
            int req_elts = nelements( markedelements( M_mesh, M_required_element_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required elements: {}", req_elts ) << std::endl;
            int next_free_req = 1;
            int next_free_nreq = req_elts + 1;
            for ( int k = 1; k <= nTetrahedra; k++ )
            {
                int iv[4], lab;
                if ( MMG3D_Get_tetrahedron( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                            &( iv[0] ), &( iv[1] ),
                                            &( iv[2] ), &( iv[3] ),
                                            &( lab ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh tetra " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                int& id_elt = k;//required ? next_free_req : next_free_nreq;
                using element_type = typename mesh_t::element_type;
                int mmgFragmentId = ( keep_relation_required_ &&  required ) ? M_id2tag[lab].first : lab;
                element_type newElem;
                newElem.setId( id_elt );
                if ( mmgFragmentId >= 0 )
                {
                    auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( mmgFragmentId );
                    if ( et == ElementsType::MESH_ELEMENTS ) // only marker come from element marker
                        newElem.setMarker( marker );
                }
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 4; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( newElem, false );
                if ( keep_relation_required_ && required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( id_elt, M_id2tag[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
                //++id_elt;
            }
        }
        int f_req_elts = 0;
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            f_req_elts = nelements( markedfaces( M_mesh, M_required_facet_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required facets: {}", f_req_elts ) << std::endl;
        }
        else  if constexpr ( dimension_v<MeshType> == 2 )
        {
            f_req_elts = nelements( markedelements( M_mesh, M_required_element_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required elements: {}", f_req_elts ) << std::endl;
        }
        int f_next_free_req = 1;
        int f_next_free_nreq = f_req_elts + 1;
        for ( int k = 1; k <= nTriangles; k++ )
        {
            int iv[3], lab;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Get_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                         &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }

                using face_type = typename mesh_t::face_type;
                face_type newElem;
                int& id_elt = k;//required ? f_next_free_req : f_next_free_nreq;
                int mmgFragmentId = ( keep_relation_required_ && required ) ? M_id2tag_face[lab].first : lab;
                newElem.setId( id_elt );
                if ( mmgFragmentId >= 0 )
                {
                    auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( mmgFragmentId );
                    if ( et == ElementsType::MESH_FACES ) // only marker come from face marker
                        newElem.setMarker( marker );
                }
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addFace( newElem );
                //id_elt++;
            }       
            if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 ) 
            {
                if ( MMG2D_Get_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                         &( lab ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using element_type = typename mesh_t::element_type;
                element_type newElem;
                int mmgFragmentId = ( keep_relation_required_ && required ) ? M_id2tag[lab].first : lab;
                newElem.setId( k );       
                if ( mmgFragmentId >= 0 )
                {
                    auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( mmgFragmentId );
                    if ( et == ElementsType::MESH_ELEMENTS ) // only marker come from element marker
                        newElem.setMarker( marker );    
                }
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( newElem, false );
                if ( keep_relation_required_ && required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, M_id2tag[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
            }
        }
        if constexpr ( dimension_v<MeshType> == 2 )
        {
            int f_req_elts = nelements( markedfaces( M_mesh, M_required_facet_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required facets: {}", f_req_elts ) << std::endl;
        }
        for ( int k = 1; k <= nEdges; k++ )
        {
            int iv[2], lab, isridge = 0;
            if constexpr ( dimension_v<MeshType> == 2 )
            {
                if ( MMG2D_Get_edge( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                     &( iv[0] ), &( iv[1] ),
                                     &( lab ), &( isridge ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using face_type = typename mesh_t::face_type;
                int mmgFragmentId = ( keep_relation_required_ && required ) ? M_id2tag_face[lab].first : lab;
                face_type newElem;
                newElem.setId( k );
                if ( mmgFragmentId >= 0 )
                {
                    auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( mmgFragmentId );
                    if ( et == ElementsType::MESH_FACES ) // only marker come from face marker
                        newElem.setMarker( marker ); 
                }
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 2; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                auto [it, ins] = out->addFace( newElem );
#if 0
                if ( required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, M_id2tag_face[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
#endif
            }
        }
        out->setMarkerNames( M_mesh->markerNames() );
        out->updateForUse();
    }
    if ( std::holds_alternative<PMMG_pParMesh>( mesh ) )
    {
        if ( PMMG_Get_meshSize( std::get<PMMG_pParMesh>( mesh ), &nVertices, &nTetrahedra, NULL, &nTriangles, NULL,
                                &nEdges ) != 1 )
        {
            ier = MMG5_STRONGFAILURE;
        }
        int corner, required, tag = 0;
        node_type n( mesh_t::nRealDim );
        for ( int k = 1; k <= nVertices; k++ )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( PMMG_Get_vertex( std::get<PMMG_pParMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                      &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh vertex " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            using point_type = typename mesh_t::point_type;
            point_type pt( k, n );
            pt.setProcessIdInPartition( out->worldComm().localRank() );
            pt.setProcessId( out->worldComm().localRank() );
            pt.setMarker( tag );
            out->addPoint( std::move( pt ) );
        }

        if constexpr ( dimension_v<MeshType> == 3 )
        {
            for ( int k = 1; k <= nTetrahedra; k++ )
            {
                int iv[4], lab;
                if ( PMMG_Get_tetrahedron( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                           &( iv[0] ), &( iv[1] ),
                                           &( iv[2] ), &( iv[3] ),
                                           &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh tetra " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }

                using element_type = typename mesh_t::element_type;
                element_type newElem;
                newElem.setId( k );
                newElem.setMarker( lab );
                newElem.setProcessIdInPartition( out->worldComm().localRank() );
                newElem.setProcessId( out->worldComm().localRank() );
                for ( int i = 0; i < 4; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( std::move( newElem ) );
            }
        }
        auto [neigh_interface, local_faces, face_to_interface] = getCommunicatorAPI();
        for ( int k = 1; k <= nTriangles; k++ )
        {
            int iv[3], lab;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( PMMG_Get_triangle( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                        &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                        &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using face_type = typename mesh_t::face_type;
                face_type newElem;
                newElem.setId( k );

                if ( auto fit = face_to_interface.find( k ); fit != face_to_interface.end() )
                {
                    LOG( INFO ) << fmt::format( "interface {} face {} with pid {}", fit->second, k, neigh_interface[fit->second].first );
                    lab = 1234567;
                }
                newElem.setMarker( lab );
                newElem.setProcessIdInPartition( out->worldComm().localRank() );
                newElem.setProcessId( out->worldComm().localRank() );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addFace( std::move( newElem ) );
            }
        }
        out->setMarkerNames( M_mesh->markerNames() );
        out->updateForUse();
#if 0
        int current_pid = M_mesh->worldComm().localRank();
        int interf = 0;
        std::vector<mesh_ptr_t> ghost_out( neigh_interface.size() );
        std::vector<mesh_ptr_t> ghost_in( neigh_interface.size() );
        std::vector<mpi::request> reqs( 2*neigh_interface.size() );
        for ( auto [pid, sz] : neigh_interface )
        {
            ext_elements_t<mesh_t> range_ghost;
            range_element_ptr_t<mesh_t> r( new range_element_t<mesh_t>() );
#if 0
            // now extract elements shared elements
            for ( auto const& welt : elements( out ) )
            {
                auto const& elt = boost::unwrap_ref( welt );
                if ( elt.hasFaceWithMarker( 1234567 ) )
                    r->push_back( boost::cref( elt ) );
            }
#else
            auto rfaces = markedfaces(out, 1234567 );
            LOG(INFO) << fmt::format("interface with {} nfaces interface {} nfaces marked {} ",
                                     pid, sz, nelements(rfaces));
            for( auto const& wf : rfaces )
            {
                auto const& f = boost::unwrap_ref( wf );
                CHECK( f.isConnectedTo0() ) << f;
                r->push_back( boost::cref( f.element0() ) );
            }
#endif
            range_ghost = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                             r->begin(),
                                             r->end(),
                                             r );
            ghost_out[interf] = createSubmesh( out, range_ghost );
            // send mesh to pid and get back the ghosts from the other side
            ghost_in[interf] = std::make_shared<mesh_t>();
            int  tag_send = (current_pid < pid)?0:1;
            int  tag_recv = (current_pid < pid)?1:0;

            reqs[2 * interf] = out->worldCommPtr()->localComm().isend( pid, tag_send, *ghost_out[interf] );
            reqs[2 * interf + 1] = out->worldCommPtr()->localComm().irecv( pid, tag_recv, *ghost_in[interf] );
            interf++;
        }
        LOG( INFO ) << fmt::format( "done sending ghosts - n requests: {}", reqs.size() );
        mpi::wait_all( reqs.begin(), reqs.end() );
#endif
    }
    // check if the bimap is empty. If so, there is no submesh link between meshes
    if ( !M_smd->bm.empty() && M_parent_mesh )
    {
        out->setSubMeshData( M_smd );
    }
    else
    {
    }
    return out;
}

template <typename MeshType>
void Remesh<MeshType>::setParameters()
{
    setMMGOptions<dimension_v<MeshType>, real_dimension_v<MeshType>>( M_params, M_mmg_mesh, M_mmg_sol );
    if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
    {
#if 0
        /* Set number of iterations */
        if ( !PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_niter, 1 ) )
        {
            throw std::logic_error( "Unable to set the number of iteration of the adaptation process." );
        }
#endif
        /* Don't remesh the surface */
        if ( !PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_nosurf, 1 ) )
        {
            throw std::logic_error( "Unable to disable surface adaptation." );
        }

        /* Compute output nodes and triangles global numbering */
        if ( !PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_globalNum, 1 ) )
        {
            throw std::logic_error( "Unable to compute output globa ids for nodes and triangles." );
        }
    }
}

template <typename MeshType>
void Remesh<MeshType>::setCommunicatorAPI()
{
    if ( std::holds_alternative<PMMG_pParMesh>( M_mmg_mesh ) )
    {
        if ( !PMMG_Set_iparameter( std::get<PMMG_pParMesh>( M_mmg_mesh ), PMMG_IPARAM_APImode, PMMG_APIDISTRIB_faces ) )
        {
            throw std::logic_error( "Unable to set the faces API mode." );
        }
        auto [ifaces_beg, ifaces_end, ifaces] = M_mesh->interProcessFaces();
        std::map<int, int> neighbor_ent;
        std::map<std::pair<int, int>, int> ids_other;
        std::vector<int> loc;
        for ( auto it = ifaces_beg; it != ifaces_end; ++it )
        {
            auto const& faceip = boost::unwrap_ref( *it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = ( elt0isGhost ) ? elt0 : elt1;

            if ( !neighbor_ent.count( eltOffProc.processId() ) )
                neighbor_ent[eltOffProc.processId()] = 0;
            ++neighbor_ent[eltOffProc.processId()];
            ids_other[std::pair{ faceip.id(), eltOffProc.processId() }] = faceip.idInOthersPartitions( eltOffProc.processId() );
        }
        /* Set the number of interfaces */
        if ( !PMMG_Set_numberOfFaceCommunicators( std::get<PMMG_pParMesh>( M_mmg_mesh ), neighbor_ent.size() ) )
        {
            throw std::logic_error( "Unable to set the number of local parallel faces." );
        }
        LOG( INFO ) << "number of interfaces: " << neighbor_ent.size() << "\n";
        int lcomm = 0;
        int current_pid = M_mesh->worldComm().localRank();
        for ( auto [pid, sz] : neighbor_ent )
        {
            LOG( INFO ) << "setting ParMmg communicator with " << pid << " (" << sz << " shared faces)" << std::endl;
            /* Set nb. of entities on interface and rank of the outward proc */
            if ( !PMMG_Set_ithFaceCommunicatorSize( std::get<PMMG_pParMesh>( M_mmg_mesh ), lcomm, pid, sz ) )
            {
                throw std::logic_error( "Unable to set internode data." );
            }
            /* Set local and global index for each entity on the interface */
            std::vector<int> local( sz );
            std::vector<std::pair<int, int>> local_sorted( sz );
            std::vector<int> local_wo_ren( sz );
            std::vector<int> global( sz, 0 );
            std::vector<int> other( sz );
            int fid = 0;
            for ( auto const& wf : interprocessfaces( M_mesh, pid ) )
            {
                auto const& f = boost::unwrap_ref( wf );
                // if current  pid is greater than pid then sort the local faces with respect to the id
                // in the partition pid so that the faces are set in the same order
                if ( current_pid < pid )
                    local_sorted[fid] = { f.id(), ids_other[std::pair{ f.id(), pid }] };
                else
                    local_sorted[fid] = { ids_other[std::pair{ f.id(), pid }], f.id() };
#if 0
                local[fid] = face_id[f.id()].first;

                global[fid] = 0;//face_id[f.id()].second;

#endif
                local_wo_ren[fid] = f.id();
                other[fid] = ids_other[std::pair{ f.id(), pid }];
                ++fid;
            }
            std::sort( local_sorted.begin(), local_sorted.end() );
            int i = 0;
            for ( auto [id1, id2] : local_sorted )
            {
                int id = ( current_pid < pid ) ? id1 : id2;
                local[i] = face_id[id].first;
                ++i;
            }
            // we have sorted the faces so that we don't have to provide a global index
            if ( !PMMG_Set_ithFaceCommunicator_faces( std::get<PMMG_pParMesh>( M_mmg_mesh ), lcomm,
                                                      local.data(),
                                                      global.data(),
                                                      0 ) )
            {
                throw std::logic_error( "Unable to set interprocess face data." );
            }
            ++lcomm;
        }
    }
}
template <typename MeshType>
typename Remesh<MeshType>::comm_t
Remesh<MeshType>::getCommunicatorAPI()
{
    LOG( INFO ) << "getCommunicatorAPI - set interface data with neighbors...";
    int n_neigh;
    if ( !PMMG_Get_numberOfFaceCommunicators( std::get<PMMG_pParMesh>( M_mmg_mesh ), &n_neigh ) )
    {
        throw std::logic_error( "Unable to get the number of local parallel faces." );
    }
    //using neighbor_ent_t = std::map<int, int>;
    //using local_face_index_t = std::vector<std::vector<int>>;
    //using comm_t = std::tuple<neighbor_ent_t, local_face_index_t>;
    neighbor_ent_t neighbor_ent;

    int** local_faces = new int*[n_neigh];
    for ( int s = 0; s < n_neigh; ++s )
    {
        int pid, sz;
        /* Set nb. of entities on interface and rank of the outward proc */
        if ( !PMMG_Get_ithFaceCommunicatorSize( std::get<PMMG_pParMesh>( M_mmg_mesh ), s, &pid, &sz ) )
        {
            throw std::logic_error( "Unable to get interface data." );
        }
        LOG( INFO ) << fmt::format( "setting ParMmg communicator {} with {} and {} shared faces)", s, pid, sz );
        neighbor_ent[s] = std::pair{ pid, sz };
        local_faces[s] = new int[sz];
    }

    if ( !PMMG_Get_FaceCommunicator_faces( std::get<PMMG_pParMesh>( M_mmg_mesh ), local_faces ) )
    {
        throw std::logic_error( "Unable to get interface faces data." );
    }

    std::map<int, int> face_to_interface;
    local_face_index_t local_face_index( n_neigh );
    for ( int s = 0; s < n_neigh; ++s )
    {
        local_face_index[s].resize( neighbor_ent[s].second );
        std::copy( local_faces[s], local_faces[s] + neighbor_ent[s].second, local_face_index[s].begin() );
        for ( int f : local_face_index[s] )
        {
            face_to_interface[f] = s;
        }
        LOG( INFO ) << fmt::format( "interface {} local face with pid {} : {}", s, neighbor_ent[s].first, local_face_index );
        delete[] local_faces[s];
    }
    delete[] local_faces;
    return std::tuple{ neighbor_ent, local_face_index, face_to_interface };
}

template <typename MeshType>
std::shared_ptr<MeshType>
Remesh<MeshType>::mmgls2Mesh( mmg_mesh_t const& mesh )
{
    int ier;
    int nVertices = 0;
    int nTetrahedra = 0;
    int nTriangles = 0;
    int nEdges = 0;

    std::shared_ptr<MeshType> out = std::make_shared<MeshType>(M_mesh->worldCommPtr());

    if ( std::holds_alternative<MMG5_pMesh>( mesh ) )
    {
        if ( MMG3D_Get_meshSize( std::get<MMG5_pMesh>( mesh ), &nVertices, &nTetrahedra, NULL, &nTriangles, NULL,
                                 &nEdges ) != 1 )
        {
            ier = MMG5_STRONGFAILURE;
        }
        int corner, required, tag = 0;
        node_type n( mesh_t::nRealDim );
        for ( int k = 1; k <= nVertices; k++ )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                       &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format("Unable to get mesh vertex {} ", k) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
            {
                if ( MMGS_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                      &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format( "Unable to get mesh vertex {} ", k ) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            if constexpr ( dimension_v<MeshType> == 2 )
            {

                if ( MMG2D_Get_vertex( std::get<MMG5_pMesh>( mesh ), &( n[0] ), &( n[1] ),
                                       &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << fmt::format( "Unable to get mesh vertex {} ", k ) << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            using point_type = typename mesh_t ::point_type;
            point_type pt( k, n );
            pt.setProcessIdInPartition( 0 );
            pt.setProcessId( 0 );
            out->addPoint( pt );
        }

        if constexpr ( dimension_v<MeshType> == 3 )
        {
            int req_elts = nelements( markedelements( M_mesh, M_required_element_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required elements: {}", req_elts ) << std::endl;
            int next_free_req = 1;
            int next_free_nreq = req_elts + 1;
            for ( int k = 1; k <= nTetrahedra; k++ )
            {
                int iv[4], lab;
                if ( MMG3D_Get_tetrahedron( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                            &( iv[0] ), &( iv[1] ),
                                            &( iv[2] ), &( iv[3] ),
                                            &( lab ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh tetra " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                int& id_elt = k;//required ? next_free_req : next_free_nreq;
                using element_type = typename mesh_t::element_type;
                element_type newElem;
                newElem.setId( id_elt );
                

                newElem.setMarker( lab );
                
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 4; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( newElem, false );
                if ( keep_relation_required_ && required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( id_elt, M_id2tag[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
                //++id_elt;
            }
        }
        int f_req_elts = 0;
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            f_req_elts = nelements( markedfaces( M_mesh, M_required_facet_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required facets: {}", f_req_elts ) << std::endl;
        }
        else  if constexpr ( dimension_v<MeshType> == 2 )
        {
            f_req_elts = nelements( markedelements( M_mesh, M_required_element_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required elements: {}", f_req_elts ) << std::endl;
        }
        int f_next_free_req = 1;
        int f_next_free_nreq = f_req_elts + 1;
        for ( int k = 1; k <= nTriangles; k++ )
        {
            int iv[3], lab;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Get_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                         &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }

                using face_type = typename mesh_t::face_type;
                face_type newElem;
                int& id_elt = k;//required ? f_next_free_req : f_next_free_nreq;
                newElem.setId( id_elt );
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addFace( newElem );
                //id_elt++;
            }       
            if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 ) 
            {
                if ( MMG2D_Get_triangle( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                         &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                         &( lab ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using element_type = typename mesh_t::element_type;
                element_type newElem;
                newElem.setId( k );
                newElem.setMarker( lab );                
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( newElem, false );
                if ( keep_relation_required_ && required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, M_id2tag[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
            }
        }
        if constexpr ( dimension_v<MeshType> == 2 )
        {
            int f_req_elts = nelements( markedfaces( M_mesh, M_required_facet_markers ) );
            LOG(INFO) << fmt::format( "[mmg->feelpp] number of required facets: {}", f_req_elts ) << std::endl;
        }
        for ( int k = 1; k <= nEdges; k++ )
        {
            int iv[2], lab, isridge = 0;
            if constexpr ( dimension_v<MeshType> == 2 )
            {
                if ( MMG2D_Get_edge( std::get<MMG5_pMesh>( M_mmg_mesh ),
                                     &( iv[0] ), &( iv[1] ),
                                     &( lab ), &( isridge ), &( required ) ) != 1 )
                {
                    LOG( ERROR ) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using face_type = typename mesh_t::face_type;
                face_type newElem;
                newElem.setId( k );
                LOG(INFO) << "lab : " << lab  ;
                newElem.setProcessIdInPartition( 0 );
                newElem.setProcessId( 0 );
                for ( int i = 0; i < 2; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                auto [it, ins] = out->addFace( newElem );
#if 0
                if ( required )
                {
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( k, M_id2tag_face[lab].second ) );
                    //std::cout << fmt::format( "[mmg->feelpp] required elt  {} tag,id: {}", id_elt, M_id2tag[lab] ) << std::endl;
                }
#endif
            }
        }
        out->updateForUse();
    }




    // Parallel
    if ( std::holds_alternative<PMMG_pParMesh>( mesh ) )
    {
        if ( PMMG_Get_meshSize( std::get<PMMG_pParMesh>( mesh ), &nVertices, &nTetrahedra, NULL, &nTriangles, NULL,
                                &nEdges ) != 1 )
        {
            ier = MMG5_STRONGFAILURE;
        }
        int corner, required, tag = 0;
        node_type n( mesh_t::nRealDim );
        for ( int k = 1; k <= nVertices; k++ )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( PMMG_Get_vertex( std::get<PMMG_pParMesh>( mesh ), &( n[0] ), &( n[1] ), &( n[2] ),
                                      &( tag ), &( corner ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh vertex " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
            }
            using point_type = typename mesh_t::point_type;
            point_type pt( k, n );
            pt.setProcessIdInPartition( out->worldComm().localRank() );
            pt.setProcessId( out->worldComm().localRank() );
            out->addPoint( std::move( pt ) );
        }

        if constexpr ( dimension_v<MeshType> == 3 )
        {
            for ( int k = 1; k <= nTetrahedra; k++ )
            {
                int iv[4], lab;
                if ( PMMG_Get_tetrahedron( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                           &( iv[0] ), &( iv[1] ),
                                           &( iv[2] ), &( iv[3] ),
                                           &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh tetra " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }

                using element_type = typename mesh_t::element_type;
                element_type newElem;
                newElem.setId( k );
                newElem.setProcessIdInPartition( out->worldComm().localRank() );
                newElem.setProcessId( out->worldComm().localRank() );
                for ( int i = 0; i < 4; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addElement( std::move( newElem ) );
            }
        }
        auto [neigh_interface, local_faces, face_to_interface] = getCommunicatorAPI();
        for ( int k = 1; k <= nTriangles; k++ )
        {
            int iv[3], lab;
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( PMMG_Get_triangle( std::get<PMMG_pParMesh>( M_mmg_mesh ),
                                        &( iv[0] ), &( iv[1] ), &( iv[2] ),
                                        &( lab ), &( required ) ) != 1 )
                {
                    LOG(ERROR) << "Unable to get mesh triangle " << k << std::endl;
                    ier = MMG5_STRONGFAILURE;
                }
                using face_type = typename mesh_t::face_type;
                face_type newElem;
                newElem.setId( k );

                if ( auto fit = face_to_interface.find( k ); fit != face_to_interface.end() )
                {
                    LOG( INFO ) << fmt::format( "interface {} face {} with pid {}", fit->second, k, neigh_interface[fit->second].first );
                    lab = 1234567;
                    newElem.setMarker( lab );
                }
                
                newElem.setProcessIdInPartition( out->worldComm().localRank() );
                newElem.setProcessId( out->worldComm().localRank() );
                for ( int i = 0; i < 3; i++ )
                    newElem.setPoint( i, out->point( iv[i] ) );
                out->addFace( std::move( newElem ) );
            }
        }
        out->setMarkerNames( M_mesh->markerNames() );
        out->updateForUse();
#if 0
        int current_pid = M_mesh->worldComm().localRank();
        int interf = 0;
        std::vector<mesh_ptr_t> ghost_out( neigh_interface.size() );
        std::vector<mesh_ptr_t> ghost_in( neigh_interface.size() );
        std::vector<mpi::request> reqs( 2*neigh_interface.size() );
        for ( auto [pid, sz] : neigh_interface )
        {
            ext_elements_t<mesh_t> range_ghost;
            range_element_ptr_t<mesh_t> r( new range_element_t<mesh_t>() );
#if 0
            // now extract elements shared elements
            for ( auto const& welt : elements( out ) )
            {
                auto const& elt = boost::unwrap_ref( welt );
                if ( elt.hasFaceWithMarker( 1234567 ) )
                    r->push_back( boost::cref( elt ) );
            }
#else
            auto rfaces = markedfaces(out, 1234567 );
            LOG(INFO) << fmt::format("interface with {} nfaces interface {} nfaces marked {} ",
                                     pid, sz, nelements(rfaces));
            for( auto const& wf : rfaces )
            {
                auto const& f = boost::unwrap_ref( wf );
                CHECK( f.isConnectedTo0() ) << f;
                r->push_back( boost::cref( f.element0() ) );
            }
#endif
            range_ghost = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                             r->begin(),
                                             r->end(),
                                             r );
            ghost_out[interf] = createSubmesh( out, range_ghost );
            // send mesh to pid and get back the ghosts from the other side
            ghost_in[interf] = std::make_shared<mesh_t>();
            int  tag_send = (current_pid < pid)?0:1;
            int  tag_recv = (current_pid < pid)?1:0;

            reqs[2 * interf] = out->worldCommPtr()->localComm().isend( pid, tag_send, *ghost_out[interf] );
            reqs[2 * interf + 1] = out->worldCommPtr()->localComm().irecv( pid, tag_recv, *ghost_in[interf] );
            interf++;
        }
        LOG( INFO ) << fmt::format( "done sending ghosts - n requests: {}", reqs.size() );
        mpi::wait_all( reqs.begin(), reqs.end() );
#endif
    }
    // check if the bimap is empty. If so, there is no submesh link between meshes
    if ( !M_smd->bm.empty() && M_parent_mesh )
    {
        out->setSubMeshData( M_smd );
    }
    else
    {
    }
    return out;
}

template <typename MeshType>
std::shared_ptr<MeshType> Remesh<MeshType>::mesh_from_ls(scalar_metric_t const& m, std::vector<int> refMarkers, std::vector<int> fluidMarker )
{ 
    // Set levelset M_mmg_sol
    this->setMetric(m);
    
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "initmesh.o.mesh");
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "initmesh.o.mesh");
    
    // Add multiphysics references
    int npar = refMarkers.size() + fluidMarker.size();
    
    if (npar > 0)
    {
        if constexpr ( dimension_v<MeshType> == 2 )
            MMG2D_Set_iparameter(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,MMG2D_IPARAM_numberOfMat,npar);
        if constexpr ( dimension_v<MeshType> == 3 )
            MMG3D_Set_iparameter(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,MMG3D_IPARAM_numberOfMat,npar);
    }

    for (int i = 0;i<npar;i++)
    {
        auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( i );
            
        if (std::find(refMarkers.begin(), refMarkers.end(), marker.value()) != refMarkers.end())
        {
            if constexpr ( dimension_v<MeshType> == 2 )
                MMG2D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_NoSplit,i,i);
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_NoSplit,i,i);
        }
        else if (std::find(fluidMarker.begin(), fluidMarker.end(), marker.value()) != fluidMarker.end())
        {
            if constexpr ( dimension_v<MeshType> == 2 )
                MMG2D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_Split,3000,2000);
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_Split,3000,2000);
        }

    }

    // Execute remeshing from levelset
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_mmg2dls(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,NULL);
           
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_mmg3dls(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,NULL);
    
    // Save final mesh in mmg format
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "finalmesh.o.mesh");
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "finalmesh.o.mesh");

    // Get and export final mesh in mesh format
    auto finalMesh = mmgls2Mesh(this->M_mmg_mesh);

    return finalMesh;
}

template <typename MeshType>
std::shared_ptr<MeshType> Remesh<MeshType>::mesh_from_ls_met(scalar_metric_t const& m, scalar_metric_t const& m_ls, std::vector<int> refMarkers, std::vector<int> fluidMarker )
{ 
    // Set levelset M_mmg_sol
    this->setMetric(m);
    this->setMetricLs(m_ls);
    
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "initmesh.o.mesh");
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "initmesh.o.mesh");
    
    // Add multiphysics references
    int npar = refMarkers.size() + fluidMarker.size();
    
    if (npar > 0)
    {
        if constexpr ( dimension_v<MeshType> == 2 )
            MMG2D_Set_iparameter(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,MMG2D_IPARAM_numberOfMat,npar);
        if constexpr ( dimension_v<MeshType> == 3 )
            MMG3D_Set_iparameter(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,MMG3D_IPARAM_numberOfMat,npar);
    }

    for (int i = 0;i<npar;i++)
    {
        auto const& [et,fragId,marker] = M_mapMmgFragmentIdToMeshFragementDesc.at( i );
            
        if (std::find(refMarkers.begin(), refMarkers.end(), marker.value()) != refMarkers.end())
        {
            if constexpr ( dimension_v<MeshType> == 2 )
                MMG2D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_NoSplit,i,i);
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_NoSplit,i,i);
        }
        else if (std::find(fluidMarker.begin(), fluidMarker.end(), marker.value()) != fluidMarker.end())
        {
            if constexpr ( dimension_v<MeshType> == 2 )
                MMG2D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_Split,3000,2000);
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Set_multiMat(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,i,MMG5_MMAT_Split,3000,2000);
        }

    }

    // Execute remeshing from levelset
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_mmg2dls(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,this->M_mmg_met);
           
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_mmg3dls(std::get<MMG5_pMesh>(this->M_mmg_mesh),this->M_mmg_sol,this->M_mmg_met);
    
    // Save final mesh in mmg format
    if constexpr ( dimension_v<MeshType> == 2 )
        MMG2D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "finalmesh.o.mesh");
    if constexpr ( dimension_v<MeshType> == 3 )
        MMG3D_saveMesh(std::get<MMG5_pMesh>(this->M_mmg_mesh), "finalmesh.o.mesh");

    // Get and export final mesh in mesh format
    auto finalMesh = mmgls2Mesh(this->M_mmg_mesh);

    return finalMesh;
}


#endif // FEELPP_HAS_MMG && FEELPP_HAS_PARMMG

} // namespace Feel
