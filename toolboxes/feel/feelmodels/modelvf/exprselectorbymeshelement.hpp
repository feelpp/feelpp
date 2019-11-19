/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2019-11-16

  Copyright (C) 2019 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_MODELS_VF_ExprSelectorByMeshElement_H
#define FEELPP_MODELS_VF_ExprSelectorByMeshElement_H 1

namespace Feel
{
namespace vf
{

template <typename IndexType>
class ExprSelectorByMeshElementMapping
{
public :
    using tag_type = uint16_type;
    ExprSelectorByMeshElementMapping() = default;

    template <typename MeshType>
    void updateForUse( std::map<std::string, elements_reference_wrapper_t<MeshType> > const& data )
        {
            uint16_type cpt=0;
            for ( auto const& [name,rangeElt] : data )
            {
                for ( auto const& eltWrap : rangeElt )
                    M_eltIdToTag[ unwrap_ref( eltWrap ).id() ] = cpt;
                M_nameToTag[name] = cpt++;
            }
        }

    tag_type idToTag( IndexType id ) const
        {
            auto itFindId = M_eltIdToTag.find( id );
            if ( itFindId != M_eltIdToTag.end() )
                return itFindId->second;
            return invalid_v<tag_type>;
        }

    tag_type nameToTag( std::string const& name ) const
        {
            auto itFindName = M_nameToTag.find( name );
            if ( itFindName != M_nameToTag.end() )
                return itFindName->second;
            return invalid_v<tag_type>;
        }
private :
    std::unordered_map<index_type,tag_type> M_eltIdToTag;
    std::map<std::string,tag_type> M_nameToTag;
};

template <typename IndexType, typename ExprType>
class ExprSelectorByMeshElement
{
public :
    using this_type = ExprSelectorByMeshElement<IndexType,ExprType>;
    using mapping_type = ExprSelectorByMeshElementMapping<IndexType>;
    using mapping_ptrtype = std::shared_ptr<mapping_type>;
    using mapping_tag_type = typename mapping_type::tag_type;
    using expression_type = ExprType;

    static const size_type context = expression_type::context;
    static const bool is_terminal = false;
    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = expression_type::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = expression_type::template HasTrialFunction<Func>::result;
    };

    template<typename Func>
    static const bool has_test_basis = expression_type::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = expression_type::template has_trial_basis<Func>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;
    typedef typename expression_type::value_type value_type;
    typedef typename expression_type::evaluate_type evaluate_type;

    ExprSelectorByMeshElement( mapping_ptrtype mapping, std::vector<std::pair<std::string,ExprType>> const& exprs )
        :
        M_mapping( mapping ),
        M_exprs( exprs )
        {}

    ExprSelectorByMeshElement( ExprSelectorByMeshElement const& ) = default;

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            uint16_type res = 0;
            for ( auto const& [name,expr] : M_exprs )
                res = std::max( res, expr.polynomialOrder() );
            return res;
        }

    //! expression is polynomial?
    bool isPolynomial() const
        {
            bool res = true;
            for ( auto const& [name,expr] : M_exprs )
                res = res && expr.isPolynomial();
            return res;
        }

    mapping_type const& mapping() const { return *M_mapping; }
    std::vector<std::pair<std::string,ExprType>> const& expressions() const { return M_exprs; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape expr_shape;
        typedef expr_shape shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_mapping( expr.mapping() ),
            M_currentTensorExpr( nullptr )
            {
                for ( auto const& [name,theexpr] : expr.expressions() )
                    M_tensorExprs.emplace( M_mapping.nameToTag( name ), tensor_expr_type( theexpr, geom, fev, feu ) );
            }
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_mapping( expr.mapping() ),
            M_currentTensorExpr( nullptr )
            {
                for ( auto const& [name,theexpr] : expr.expressions() )
                    M_tensorExprs.emplace( M_mapping.nameToTag( name ), tensor_expr_type( theexpr, geom, fev ) );
            }
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_mapping( expr.mapping() ),
            M_currentTensorExpr( nullptr )
            {
                for ( auto const& [name,theexpr] : expr.expressions() )
                    M_tensorExprs.emplace( M_mapping.nameToTag( name ), tensor_expr_type( theexpr, geom ) );
            }

        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            this->selectSubTensor( geom );
            if ( M_currentTensorExpr )
                M_currentTensorExpr->update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            this->selectSubTensor( geom );
            if ( M_currentTensorExpr )
                M_currentTensorExpr->update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            this->selectSubTensor( geom );
            if ( M_currentTensorExpr )
                M_currentTensorExpr->update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            this->selectSubTensor( geom );
            if ( M_currentTensorExpr )
                M_currentTensorExpr->update( geom, face );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( M_currentTensorExpr )
                return M_currentTensorExpr->evalijq( i, j, c1, c2, q );
            else
                return value_type(0);
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( M_currentTensorExpr )
                return M_currentTensorExpr->evaliq( i, c1, c2, q );
            else
                return value_type(0);
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( M_currentTensorExpr )
                return M_currentTensorExpr->evalq( c1, c2, q );
            else
                return value_type(0);
        }

    private :
        void selectSubTensor( Geo_t const& geom )
            {
                IndexType eid = vf::detail::ExtractGm<Geo_t>::get( geom )->id();
                mapping_tag_type tag = M_mapping.idToTag( eid );
                if ( tag != invalid_v<mapping_tag_type> )
                {
                    auto itFindTensorExpr = M_tensorExprs.find( tag );
                    if ( itFindTensorExpr !=  M_tensorExprs.end() )
                    {
                        M_currentTensorExpr = &itFindTensorExpr->second;
                        return;
                    }
                }
                M_currentTensorExpr = nullptr;
            }
    private :
        mapping_type const& M_mapping;
        std::map<uint16_type,tensor_expr_type> M_tensorExprs;
        tensor_expr_type* M_currentTensorExpr;
    };

private :
    mapping_ptrtype M_mapping;
    std::vector<std::pair<std::string,ExprType>> M_exprs;
};


template <typename IndexType, typename ExprType>
inline
Expr<ExprSelectorByMeshElement<IndexType,ExprType>>
    expr( typename ExprSelectorByMeshElement<IndexType,ExprType>::mapping_ptrtype mapping, std::vector<std::pair<std::string,ExprType>> const& exprs )
{
    typedef ExprSelectorByMeshElement<IndexType,ExprType> esbme_t;
    return Expr< esbme_t >( esbme_t( mapping,exprs ) );
}

} // namespace vf
} // namespace Feel

#endif
