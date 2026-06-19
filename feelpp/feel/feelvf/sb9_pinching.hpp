/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Feel++ Consortium

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
/**
   \file sb9_pinching.hpp
   \brief SB9 pinch shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_PINCHING_HPP
#define FEELPP_VF_SB9_PINCHING_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects the SB9 membrane or bending coefficient family.
 *
 * The values are used as compile-time tags by \ref SB9BendingKernelCache and
 * \ref SB9VectorOperator.
 */
enum class SB9PinchingKind
{
    /// Mid-surface membrane coefficient matrix, usually denoted Bm0.
    Bpc,
    /// Linear-through-thickness bending coefficient matrix, usually denoted Bb0.
    Bpz,  // Bpzeta ?  complètement à part Bpw ?  // mais lui ne sera pas de taille nombre de noeuds comme les autres, mais ça sera juste un scalaire donc faut coder à part ?
    Bpw
};

/**
 * \brief Element-local cache for SB9 membrane and bending coefficients.
 *
 * The cache reuses \ref SB9KernelBase for frame and Jacobian data and adds the
 * bending derivative coefficients derived from the shell geometry `vgamma`
 * modes. The coefficients are later projected into Mandel symmetric storage by
 * \ref fillVectorCoefficients.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9PinchingKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by SB9 bending operators.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build membrane/bending coefficient data for one element.
     *
     * The membrane part is available directly from the geometry cache through
     * `bx` and `by`; this constructor precomputes the two bending derivative
     * coefficients for each node.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9PinchingKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        // Formulation doc A. Gravouil
        for ( uint16_type node = 0; node < node_count/2; ++node )
        {
            M_pinching[node] = value_type( -2 ) * this->M_data.bz( node );
        }
        for ( uint16_type node = node_count/2; node < node_count; ++node )
        {
            M_pinching[node] = value_type( 2 ) * this->M_data.bz( node );
        }


        double sbpz = 0.0;
        std::cout << "DEBUG: Bpz = ";
        for (auto const& v : M_pinching) 
        {
            sbpz += v;
            std::cout << v << " ";
        }
        std::cout << "DEBUG: sum Bpz = " << sbpz << std::endl;


        double sbpc = 0.0;
        std::cout << "DEBUG: Bpc = ";
        for (int i = 0; i < 8; i++)
        {
            auto hall_zi = this->M_data.bz( i );
            sbpc += hall_zi;
            std::cout << hall_zi << " ";

        } 
        std::cout << "DEBUG: sum Bpc = " << sbpc << std::endl;

        std::cout << "DEBUG: Bpw = " << - value_type( 4 ) / this->M_data.thickness << std::endl;
        std::cout << "DEBUG: Bpw bis = " << value_type( -4 ) / this->M_data.thickness << std::endl;



        // // Formulation code Matlab
        // for ( uint16_type node = 0; node < node_count; ++node )
        // {
        //     M_pinching[node] = this->M_data.bz( node ) 
        //                        + this->M_data.vgamma( node, 0 )*this->M_invJ0( 2, 1 ) 
        //                        + this->M_data.vgamma( node, 1 )*this->M_invJ0( 2, 0 );
        // }

        // // Formulation Thèse Dia
        // for ( uint16_type node = 0; node < node_count; ++node )
        // {
        //     M_pinching[node] = value_type(0.5) / this->M_data.thickness;
        // }

    }

    /**
     * \brief Fill one SB9 bending-family coefficient vector.
     *
     * \tparam Kind Compile-time selector, either \ref SB9BendingKind::Bm0 or
     *         \ref SB9BendingKind::Bb0.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in Feel++ symmetric storage order.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    // template <SB9PinchingKind Kind, typename VectorType>
    // void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    // {
    //     static_assert( Kind == SB9PinchingKind::Bpc || Kind == SB9PinchingKind::Bpz,
    //                    "unsupported SB9 pinching vector kind" );

    //     if constexpr ( Kind == SB9PinchingKind::Bpc )
    //         this->fillPinchingCoefficients( coeff, component, this->M_data.bz( node ) );
    //     else // if constexpr ( Kind == SB9PinchingKind::Bpz )
    //         this->fillPinchingCoefficients( coeff, component, M_pinching[node] );
    // }


    template <SB9PinchingKind Kind, typename VectorType>
    void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9PinchingKind::Bpc || Kind == SB9PinchingKind::Bpz || Kind == SB9PinchingKind::Bpw,
                       "unsupported SB9 pinching vector kind" );

        if constexpr ( Kind == SB9PinchingKind::Bpc )
            this->fillPinchingCoefficients( coeff, component, this->M_data.bz( node ) );
        else // if constexpr ( Kind == SB9PinchingKind::Bpz )
            this->fillPinchingCoefficients( coeff, component, M_pinching[node] );
        // si je considère Bpw comme étant un VectorCoefficients j'ai une erreur de segmentation (test qui tourne en boucle)
        // else
        //     coeff( 2 ) =  - value_type( 4 ) / this->M_data.thickness;
            
    }
    /*
    template <SB9PinchingKind Kind, typename VectorType>
    void fillVectorw9Coefficients( VectorType& coeffw9 ) const
    {
        static_assert( Kind == SB9PinchingKind::Bpw,
                       "unsupported SB9 pinching scalar kind" );

        // auto bz = - value_type( 4 ) / this->M_data.thickness;
        // this->fillPinchingW9Coefficients( coeff, bz );
        coeffw9(2) =  - value_type( 4 ) / this->M_data.thickness;   // pareil le "-" n'est pas pris en compte
        // std::cout << "DEBUG coeff(2) Bpw = " << coeff(2) << std::endl;
    }
    */
    template <SB9PinchingKind Kind, typename VectorType>
    void fillScalarCoefficients( VectorType& coeff ) const
    {
        static_assert( Kind == SB9PinchingKind::Bpw,
                       "unsupported SB9 pinching scalar kind" );

        // auto bz = - value_type( 4 ) / this->M_data.thickness;
        // this->fillPinchingW9Coefficients( coeff, bz );
        coeff(2) =  - value_type( 4 ) / this->M_data.thickness;
        // std::cout << "DEBUG coeff(2) Bpw = " << coeff(2) << std::endl;
    }


    // // ça fonctionne avec SB9ScalarCoefficient mais faudrait plutot le considérer comme un vectorCoefficient ? de taille stockage_size*1
    // template <SB9PinchingKind Kind>
    // void fillScalarCoefficients( value_type& coeff ) const
    // {
    //     static_assert( Kind == SB9PinchingKind::Bpw,
    //                    "unsupported SB9 pinching scalar kind" );

    //     coeff =  - value_type( 4 ) / this->M_data.thickness;
    // }

    // template <SB9PinchingKind Kind, typename VectorType>
    // void fillScalarCoefficients( VectorType& coeff ) const  // VectorType -> le même que les autres ?
    // {
    //     static_assert( Kind == SB9PinchingKind::Bpw,
    //                    "unsupported SB9 pinching scalar kind" );

    //     coeff(2) =  - value_type( 4 ) / this->M_data.thickness;
    // }

// // -----------------

//     template <SB9PinchingKind Kind, typename VectorType>
//     void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
//     {
//         static_assert( Kind == SB9PinchingKind::Bpc || Kind == SB9PinchingKind::Bpz || Kind == SB9PinchingKind::Bpw,
//                        "unsupported SB9 pinching vector kind" );

//         if constexpr ( Kind == SB9PinchingKind::Bpc )
//             this->fillPinchingCoefficients( coeff, component, this->M_data.bz( node ) );
//         else if constexpr ( Kind == SB9PinchingKind::Bpz )
//             this->fillPinchingCoefficients( coeff, component, M_pinching[node] );
//         else 
//             this->fillPinchingW9Coefficients( coeff, - value_type( 4 ) / this->M_data.thickness );
//     }
//     // mais faut un coeff2 ou un truc du genre faut pas que les deux coeff se regroupe - c'est pour ça faut un sb9ScalarCoefficient



private:
    /// Precomputed bending derivative coefficients per node and in-plane axis.
    std::array<value_type, node_count> M_pinching{};
};


// //============================================
// // POUR PRENDRE EN COMPTE LA GEO DE 1 NOEUD ET 1 COMPOSANTE
// template <typename GeometryDataType>
// class SB9Pinchingw9KernelCache : public SB9KernelBase<GeometryDataType> // public SB9w9KernelBase<GeometryDataType>
// {
// public:
//     /// Shared SB9 geometry cache base.
//     using base_type = SB9w9KernelBase<GeometryDataType>;
//     /// Scalar type used by the geometry and generated coefficients.
//     using value_type = typename base_type::value_type;
//     /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
//     static constexpr uint16_type node_count = 1;// base_type::node_count;
//     /// Number of displacement components handled by SB9 bending operators.
//     static constexpr uint16_type component_count = 1;// base_type::component_count;

//     /**
//      * \brief Build membrane/bending coefficient data for one element.
//      *
//      * The membrane part is available directly from the geometry cache through
//      * `bx` and `by`; this constructor precomputes the two bending derivative
//      * coefficients for each node.
//      *
//      * \param data Element-local shell geometry data.
//      */
//     explicit SB9Pinchingw9KernelCache( GeometryDataType const& data )
//         :
//         base_type( data )
//     {
//         M_pinchingw9 = - value_type( 4 ) / this->M_data.thickness;
//     }

//     template <SB9PinchingKind Kind, typename VectorType>
//     void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
//     {
//         static_assert( Kind == SB9PinchingKind::Bpw,
//                        "unsupported SB9 pinching vector kind" );

//         this->fillPinchingw9Coefficients( coeff, component, M_pinchingw9 );
//         // coeff( 2 ) = M_pinchingw9;      // - value_type( 4 ) / this->M_data.thickness;
//     }

//     // template <SB9PinchingKind Kind, typename VectorType>
//     // void fillScalarCoefficients( VectorType& coeff ) const  // VectorType -> le même que les autres ?
//     // {
//     //     static_assert( Kind == SB9PinchingKind::Bpw,
//     //                    "unsupported SB9 pinching scalar kind" );

//     //     coeff(2) =  - value_type( 4 ) / this->M_data.thickness;
//     // }



// private:
//     /// Precomputed bending derivative coefficients per node and in-plane axis.
//     value_type M_pinchingw9;
// };

// //============================================






} // namespace detail

/**
 * \brief Build the SB9 mid-surface membrane coefficient expression.
 *
 * The returned expression evaluates the `Bm0` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bm0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bpc( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9PinchingKernelCache,
                                                detail::SB9PinchingKind::Bpc>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 bending coefficient expression.
 *
 * The returned expression evaluates the `Bb0` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bb0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bpz( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9PinchingKernelCache,
                                                detail::SB9PinchingKind::Bpz>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bpw( ProxyType const& proxy )
{
    // using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    // using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,   // SB9StabilizationOperator ou créer complètement un SB9ScalarOperator
    //                                             proxy_type::role,
    //                                             detail::SB9Pinchingw9KernelCache,
    //                                             detail::SB9PinchingKind::Bpw>;
    // return Expr<expr_type>( expr_type( proxy.element() ) );

    // using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    // using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,   // SB9StabilizationOperator ou créer complètement un SB9ScalarOperator
    //                                             proxy_type::role,
    //                                             detail::SB9PinchingKernelCache,
    //                                             detail::SB9PinchingKind::Bpw>;
    // return Expr<expr_type>( expr_type( proxy.element() ) );

    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9ScalarOperator<typename proxy_type::element_type,   // SB9StabilizationOperator ou créer complètement un SB9ScalarOperator
                                                proxy_type::role,
                                                detail::SB9PinchingKernelCache,
                                                detail::SB9PinchingKind::Bpw>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 membrane-plus-bending Mandel strain expression.
 *
 * This helper combines the membrane and bending coefficient expressions as
 * `Bm0 + zeta * Bb0` in the in-plane symmetric-storage entries and returns a
 * six-component Mandel vector with out-of-plane components set to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \return Feel++ Mandel-vector expression for the SB9 membrane/bending strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9Pinching( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    auto bpc = sb9Bpc( proxy );
    auto bpz = sb9Bpz( proxy );

    return mandel_vec<3>( cst( 0.0 ),
                        cst( 0.0 ),
                        component<2,0>( bpc ) + zetaExpr * component<2,0>( bpz ),
                        cst( 0.0 ),
                        cst( 0.0 ),
                        cst( 0.0 ) );
}

template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9PinchingW9( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    auto bpw = sb9Bpw( proxy );

    return mandel_vec<3>( cst( 0.0 ),  
                          cst( 0.0 ),
                          zetaExpr * component<2,0>( bpw ),
                          cst( 0.0 ),
                          cst( 0.0 ),
                          cst( 0.0 ) );

    // return mandel_vec<3>( cst( 0.0 ),
    //                       cst( 0.0 ),
    //                       zetaExpr * bpw,
    //                       cst( 0.0 ),
    //                       cst( 0.0 ),
    //                       cst( 0.0 ) );

    // return mandel_vec<3>( cst( 0.0 ),  
    //                       cst( 0.0 ),
    //                       - cst(4.0) / shellThickness()* zetaExpr, // * proxy,
    //                       cst( 0.0 ),
    //                       cst( 0.0 ),
    //                       cst( 0.0 ) );
}
} // namespace vf
} // namespace Feel

#endif
