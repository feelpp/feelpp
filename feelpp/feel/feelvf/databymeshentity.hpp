/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 22 Oct 2020

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

namespace Feel
{
namespace vf
{

template <typename DataByMeshEntityType>
class EvaluateDataByMeshEntity
{
public :
    static const size_type context = 0;
    static const bool is_terminal = false;
    template <typename Func>
    struct HasTestFunction { static const bool result = false; };
    template <typename Func>
    struct HasTrialFunction { static const bool result = false; };
    template <typename Func>
    static const bool has_test_basis = false;
    template <typename Func>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    using this_type = EvaluateDataByMeshEntity<DataByMeshEntityType>;
    using data_by_mesh_entity_type = DataByMeshEntityType;
    using value_type = typename data_by_mesh_entity_type::value_type;

    EvaluateDataByMeshEntity( DataByMeshEntityType const& dataOnMeshEntities )
        :
        M_dataOnMeshEntities( dataOnMeshEntities )
        {}

    DataByMeshEntityType const& dataOnMeshEntities() const { return M_dataOnMeshEntities; }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using value_type = this_type::value_type;
        //using key_type = key_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_value( 0 )
        {}

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_value( 0 )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_value( 0 )
        {}
        template <typename IM>
        void init( IM const& im ) {}

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            this->update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            this->update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_value = 0;
            bool hasValue = false;
            gmc_ptrtype gmc = vf::detail::ExtractGm<Geo_t>::get( geom );
            if ( !gmc->isOnSubEntity() )
            {
                auto const& elt = gmc->element();
                bool isSameMesh = M_expr.dataOnMeshEntities().mesh()->isSameMesh( elt.mesh() );
                bool gmcEltIsSubMesh = elt.mesh()->isSubMeshFrom( M_expr.dataOnMeshEntities().mesh() );
                if ( M_expr.dataOnMeshEntities().entityType() == ElementsType::MESH_EDGES )
                {
                    if ( isSameMesh )
                    {
                        for ( uint16_type p=0;p<elt.nEdges();++p )
                        {
                            index_type edgeId = elt.edge(p).id();
                            if ( auto valOpt = M_expr.dataOnMeshEntities().valueAtEntityIdIfExists( edgeId ) )
                            {
                                M_value = *valOpt;
                                hasValue = true;
                                break;
                            }
                        }
                    }
                    else if ( gmcEltIsSubMesh )
                    {
                        if constexpr ( gmc_type::element_type::nDim == 1 )
                        {
                            index_type edgeId = elt.mesh()->subMeshToMesh( M_expr.dataOnMeshEntities().mesh(),elt.id() );
                            if ( auto valOpt = M_expr.dataOnMeshEntities().valueAtEntityIdIfExists( edgeId ) )
                            {
                                M_value = *valOpt;
                                hasValue = true;
                            }
                        }
                        else
                            CHECK( false ) << "TODO";
                    }
                    else
                    {
                        CHECK( false ) << "Something wrong";
                    }
                }
                else
                    CHECK( false ) << "TODO";
                CHECK( hasValue ) << "value not found in data";
            }
            else if ( gmc->isOnFace() )
            {
                CHECK( false ) << "TODO data on faces";
            }
            else if constexpr ( gmc_type::nDim == 3 && gmc_type::subEntityCoDim == 2 )
            {
                auto const& elt = gmc->element();
                CHECK( false ) << "TODO data on edges";
            }
            else
            {
                CHECK( false ) << "TODO data on points";
            }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            this->update( geom );
        }

        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_value;
        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_value;
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_value;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_value;
        }
    private :
        this_type const& M_expr;
        value_type M_value;
    };

private :
    DataByMeshEntityType const& M_dataOnMeshEntities;
};



template <typename MeshType>
auto dataByMeshEntityExpr( std::shared_ptr<CollectionOfDataByMeshEntity<MeshType>> const& dcome, std::string const& field )
{
    using _expr_type = EvaluateDataByMeshEntity< DataByMeshEntity<MeshType> >;
    return  Expr<_expr_type>( _expr_type( dcome->get( field ) ) );
}

} // vf
} // Feel
