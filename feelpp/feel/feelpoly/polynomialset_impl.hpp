/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 11:18:41 2014

   Copyright (C) 2014-2016 Feel++ Consortium

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

namespace Feel {

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
resizeAndSet( rank_t<0> )
{
#if 0
    Eigen::Tensor<value_type,3> i_grad( 1, nRealDim, 1 );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_first_derivative_normal_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_dn( 1,1 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ||
         ( vm::has_dynamic_basis_function_v<context> && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
    {
        Eigen::Tensor<value_type,3> i_hessian( nRealDim, nRealDim, 1 );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
        if ( vm::has_laplacian_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasLAPLACIAN( this->dynamicContext() ) ) )
        {
            Eigen::Tensor<value_type,2> i_lap( 1, 1 );
            std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
        }
    }
#endif
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
resizeAndSet( rank_t<1> )
{
#if 0
    Eigen::Tensor<value_type,3> i_grad( nComponents1, nRealDim, 1 );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_symm_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasSYMM( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_symm( nComponents1, nRealDim );
        std::fill( M_symm_grad.data(), M_symm_grad.data()+M_symm_grad.num_elements(), i_symm.constant(0.) );
    }
    if ( vm::has_first_derivative_normal_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_dn( nComponents1,1 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_div_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasDIV( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_div( 1, 1 );
        std::fill( M_div.data(), M_div.data()+M_div.num_elements(), i_div.constant(0.) );
    }

    if ( vm::has_curl_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasCURL( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_curl( (nComponents1==3)?nComponents1:1, 1 );
        std::fill( M_curl.data(), M_curl.data()+M_curl.num_elements(), i_curl.constant(0.) );
    }

    if ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ||
         ( vm::has_dynamic_basis_function_v<context> && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) )
         )
    {
        Eigen::Tensor<value_type,3> i_hessian( nComponents1, nRealDim, nRealDim );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
        if ( vm::has_laplacian_v<context> || vm::has_second_derivative_v<context> ||
             ( vm::has_dynamic_basis_function_v<context> && ( hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
        {
            Eigen::Tensor<value_type,2> i_lap( nComponents1, 1 );
            std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
        }
    }
#endif
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
resizeAndSet( rank_t<2> )
{
#if 0
    Eigen::Tensor<value_type,3> i_grad( nComponents1, nComponents2, nRealDim );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_first_derivative_normal_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_dn( nComponents1, nComponents2 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_div_v<context> || ( vm::has_dynamic_basis_function_v<context> && hasDIV( this->dynamicContext() ) ) )
    {
        Eigen::Tensor<value_type,2> i_div( nComponents1, 1 );
        std::fill( M_div.data(), M_div.data()+M_div.num_elements(), i_div.constant(0.) );
    }
#endif
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
update( geometric_mapping_context_ptrtype const& __gmc,
        precompute_ptrtype const& __pc )
{
    M_pc = __pc;
    M_gmc = __gmc ;

    if ( ( M_npoints != M_gmc->nPoints() ) ||
         ( M_npoints != __pc->nPoints() ) ||
         ( M_npoints != M_phi.shape()[1] ) )
    {
        //std::cout << "M_npoints = "  << M_npoints << "\n";
        //std::cout << "pc->npoints = "  << __pc->nPoints() << "\n";
        //std::cout << "phi->npoints = "  << M_phi.shape()[1] << "\n";
        //std::cout << "gmc->npoints = "  << M_gmc->nPoints() << "\n";
        M_npoints = __pc->nPoints();

        //const int ntdof = nDof*nComponents1;
        const int ntdof = this->nDofs();
        M_phi.resize( boost::extents[ntdof][M_npoints] );
        //M_gradphi.resize( boost::extents[ntdof][M_npoints] );

        int npoints_firstderivative = do_optimization_p1? 1 : M_npoints;
        int npoints_secondderivative = do_optimization_p2? 1 : M_npoints;

        // normal component
        if constexpr ( rank >=1 )
        {
            if constexpr ( vm::has_normal_component_v<context> )
            {
                M_normal_component.resize( boost::extents[ntdof][M_npoints] );
            }
            else if constexpr( vm::has_dynamic_basis_function_v<context> )
            {
                if ( hasNORMAL_COMPONENT( this->dynamicContext() ) )
                    M_normal_component.resize( boost::extents[ntdof][M_npoints] );
            }
        }

        // trace
        if constexpr ( rank == 2 )
        {
            if constexpr ( vm::has_trace_v<context> )
            {
                M_trace.resize( boost::extents[ntdof][M_npoints] );
            }
            else if constexpr ( vm::has_dynamic_basis_function_v<context> )
            {
                if ( hasTRACE( this->dynamicContext() ) )
                    M_trace.resize( boost::extents[ntdof][M_npoints] );
            }
        }

        // gradient
        if constexpr ( vm::has_grad_v<context> || vm::has_first_derivative_v<context> ||
                       ( second_derivative_require_grad && ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ) ) )
        {
            M_grad.resize( boost::extents[ntdof][npoints_firstderivative] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasGRAD( this->dynamicContext() ) || hasFIRST_DERIVATIVE( this->dynamicContext() ) ||
                 ( second_derivative_require_grad && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
                M_grad.resize( boost::extents[ntdof][npoints_firstderivative] );
        }

        // symm
        if constexpr ( vm::has_symm_v<context>  )
        {
            M_symm_grad.resize( boost::extents[ntdof][npoints_firstderivative] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasSYMM( this->dynamicContext() ) )
                M_symm_grad.resize( boost::extents[ntdof][npoints_firstderivative] );
        }

        // div
        if constexpr ( vm::has_div_v<context> )
        {
            M_div.resize( boost::extents[ntdof][npoints_firstderivative] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasDIV( this->dynamicContext() ) )
                M_div.resize( boost::extents[ntdof][npoints_firstderivative] );
        }

        // curl
        if constexpr ( vm::has_curl_v<context> )
        {
            M_curl.resize( boost::extents[ntdof][npoints_firstderivative] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasCURL( this->dynamicContext() ) )
                M_curl.resize( boost::extents[ntdof][npoints_firstderivative] );
        }

        // first derivative normal
        if constexpr ( vm::has_first_derivative_normal_v<context> )
        {
            M_dn.resize( boost::extents[ntdof][M_npoints] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) )
                M_dn.resize( boost::extents[ntdof][M_npoints] );
        }

        // hessian
        if constexpr ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context>  )
        {
            M_hessian.resize( boost::extents[ntdof][M_npoints] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) )
                M_hessian.resize( boost::extents[ntdof][M_npoints] );
        }

        // laplacian
        if constexpr ( vm::has_laplacian_v<context> )
        {
            M_laplacian.resize( boost::extents[ntdof][M_npoints] );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasLAPLACIAN( this->dynamicContext() ) )
                M_laplacian.resize( boost::extents[ntdof][M_npoints] );
        }

        resizeAndSet( rank_t<rank>() );
    }

    M_phi = M_pc.get()->phi();
    M_gradphi = M_pc.get()->gradPtr();

    update( __gmc );
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateGrad( geometric_mapping_context_type* thegmc, rank_t<0> )
{
    tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> B( thegmc->B( 0 ).data(), gmc_type::NDim, gmc_type::PDim );
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = nDof;
    Eigen::array<int, 3> tensorGradShapeAfterContract{{1, 1, nRealDim}};
    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 1)}};
    for ( uint16_type i = 0; i < I; ++i )
    {
        auto const& g_phi_i = (*M_gradphi)[i];
        for ( uint16_type q = 0; q < Q; ++q )
        {
            if constexpr(!gmc_type::is_linear)
                            new (&B) tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type>(thegmc->B( q ).data(), gmc_type::NDim, gmc_type::PDim );
            // grad = (gradphi_1,...,gradphi_nRealDim) * B^T
            M_grad[i][q].reshape( tensorGradShapeAfterContract ) = ((*M_gradphi)[i][q].contract( B,dims ));
            //M_grad[i][q] = (g_phi_i[q].contract( B,dims ));
#if 0
            M_dx[i][q] = M_grad[i][q].col( 0 );
            if ( NDim == 2 )
                M_dy[i][q] = M_grad[i][q].col( 1 );

            if ( NDim == 3 )
                M_dz[i][q] = M_grad[i][q].col( 2 );
#endif
        }
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateFirstDerivativeNormal( geometric_mapping_context_type* thegmc, rank_t<0> )
{
    // const uint16_type I = M_ref_ele->nbDof()*nComponents;
    // const uint16_type Q = nPoints();
    const uint16_type Q = M_npoints;//do_optimization_p1?1:M_npoints;
    const uint16_type I = nDof;

    for ( uint16_type i = 0; i < I; ++i )
    {
        for ( uint16_type q = 0; q < Q; ++q )
        {
            M_dn[i][q].setZero();
            //M_dn[i][q]( 0,0 ) = 0;
            for ( uint16_type l = 0; l < gmc_type::NDim; ++l )
            {
                //M_dn[i][q]( 0,0 ) += M_grad[i][q]( 0,l,0 ) * thegmc->unitNormal( l, q );
                M_dn[i][q]( 0,0 ) += this->grad(i,q)( 0,l,0 ) * thegmc->unitNormal( l, q );
            }
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateHessian( geometric_mapping_context_type* thegmc, rank_t<0> )
{
    precompute_type* __pc = M_pc.get().get();
    //auto* __gmc_pc = thegmc->pc().get();
    const uint16_type Q = do_optimization_p2?1:M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = nDof; //M_ref_ele->nbDof();
    //hess_type L;
    //Eigen::array<int, 3> tensorHessShapeAfterContract{{nRealDim, 1, nRealDim}};
    Eigen::array<int, 3> tensorHessShapeAfterContract{{1, nRealDim, nRealDim}};
    Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(1, 1)}};
    Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(0, 1)}};
    Eigen::array<dimpair_t, 1> dimsh = {{dimpair_t(0, 0)}};
    Eigen::array<dimpair_t, 1> tracedims = {{dimpair_t(0, 1)}};
    for ( uint16_type q = 0; q < Q; ++q )
    {
        tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim ,value_type> B ( thegmc->B( q ).data(), gmc_type::NDim, gmc_type::PDim );
        for ( uint16_type i = 0; i < I; ++i )
        {
            //M_hessian[i][q] = B.contract( __pc->hessian(i,q).contract(B,dims1), dims2);

            if constexpr ( !geometric_mapping_context_type::is_linear )//Geo_t::nOrder > 1 || !convex_type::is_simplex )
            {
                // gmc hessian NxPxP
                // pc hessian PxP
                // grad N
                // B NxP
                // hessian: N x N
                //std::cout << "Grad=" << M_grad[i][q].chip(0,0).chip(0,1) << "\n";
                tensor2_fixed_size_t<PDim,PDim,value_type> H0 = thegmc->hessian(q).contract( M_grad[i][q].chip(0,0).chip(0,1), dimsh );
                //std::cout << "H0=" << H0 << std::endl;
                tensor2_fixed_size_t<PDim,PDim,value_type> H00 = __pc->hessian(i,q).chip(0,2);
                //std::cout << "H00=" << H00 << std::endl;
                tensor2_fixed_size_t<PDim,PDim,value_type> H1 = H00 - H0;
                //std::cout << "H1=" << H1 << std::endl;
                M_hessian[i][q].chip(0,2) = H1.contract( B, dims1 ).contract( B, dims2 );
                //std::cout << "H=" << M_hessian[i][q] << std::endl;
            }
            else
            {
                // hess = (hess_phi * B^T)*B
                auto H1 = __pc->hessian(i,q).contract(B,dims1);
                M_hessian[i][q].reshape( tensorHessShapeAfterContract ) = H1.contract( B, dims2 );
            }

        } // q
    } // i
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateLaplacian( geometric_mapping_context_type* thegmc, rank_t<0> )
{
    const uint16_type Q = do_optimization_p2?1:M_npoints;
    const uint16_type I = nDof;
    for ( uint16_type q = 0; q < Q; ++q )
    {
        for ( uint16_type i = 0; i < I; ++i )
        {
#if 0
            M_laplacian[i][q] = M_hessian[i][q].trace(tracedims);
#else
            M_laplacian[i][q].setZero();
            for( int c = 0; c < nRealDim; ++c )
                M_laplacian[i][q](0,0) += M_hessian[i][q]( c, c, 0 );
            //M_laplacian[i][q](0,0) = M_hessian[i][q].trace();
#endif
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0> )
{
    geometric_mapping_context_type* thegmc = __gmc.get();

    if constexpr ( vm::has_grad_v<context> || vm::has_first_derivative_v<context>  ||
                   ( second_derivative_require_grad && ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ) ) )
    {
        this->updateGrad( thegmc, rank_t<0>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasGRAD( this->dynamicContext() ) || hasFIRST_DERIVATIVE( this->dynamicContext() ) ||
             ( second_derivative_require_grad && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
            this->updateGrad( thegmc, rank_t<0>{} );
    }

    if constexpr ( vm::has_first_derivative_normal_v<context> )
    {
        this->updateFirstDerivativeNormal( thegmc, rank_t<0>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) )
            this->updateFirstDerivativeNormal( thegmc, rank_t<0>{} );
    }


    if constexpr ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> )
    {
        this->updateHessian( thegmc, rank_t<0>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) )
            this->updateHessian( thegmc, rank_t<0>{} );
    }

    if constexpr ( vm::has_laplacian_v<context> )
    {
        this->updateLaplacian( thegmc, rank_t<0>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasLAPLACIAN( this->dynamicContext() ) )
            this->updateLaplacian( thegmc, rank_t<0>{} );
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateNormalComponent( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = M_npoints;
    const uint16_type I = M_normal_component.shape()[0];
    tensor_map_fixed_size_matrix_t<gmc_type::NDim,1,value_type> N ( thegmc->unitNormal( 0 ).data(),gmc_type::NDim,1 );

    for ( uint16_type i = 0; i < I; ++i )
        for ( uint16_type q = 0; q < Q; ++q )
        {
            if constexpr(!gmc_type::is_linear)
                            new (&N) tensor_map_fixed_size_matrix_t<gmc_type::NDim, 1,value_type>(thegmc->unitNormal( q ).data(), gmc_type::NDim, 1 );
            auto &r = M_normal_component[i][q](0,0);
            r= 0;
            auto* p = M_phi[i][q].data();
            auto* n = N.data();
            for( int c = 0; c < gmc_type::NDim; ++c )
                r += *(p+c)*(*(n+c));
        }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateGrad( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = M_grad.shape()[0];

    typedef typename boost::multi_array<value_type,4>::index_range range;
    Eigen::array<int, 3> tensorGradShapeAfterContract{{nComponents1, 1, nRealDim}};
    Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(1, 1)}};
    Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(1, 0)}};

    for ( uint16_type q = 0; q < Q; ++q )
    {
        tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> K ( thegmc->K( q ).data(),gmc_type::NDim, gmc_type::PDim );
        tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> B ( thegmc->B( q ).data(),gmc_type::NDim, gmc_type::PDim );
        for ( uint16_type i = 0; i < I; ++i )
        {
            auto & gradiq = M_grad[i][q];

            if constexpr ( is_hdiv_conforming )
            {
                auto v = (*M_gradphi)[i][q].contract(B,dims1);
                gradiq = K.contract(v,dims2)/thegmc->J(q);
            }
            else if constexpr ( is_hcurl_conforming )
            {
                //auto v = (*M_gradphi)[i][q].contract(B,dims1);
                gradiq = B.contract((*M_gradphi)[i][q].contract(B,dims1),dims2);
            }
            else
            {
                //std::cout << "left : " << M_grad[i][q].dimensions() << std::endl;
                //std::cout << "right : " << (*M_gradphi)[i][q].contract(B,dims1).dimensions() << std::endl;
                //M_grad[i][q] = (*M_gradphi)[i][q].contract(B,dims1);
                gradiq.reshape( tensorGradShapeAfterContract ) = (*M_gradphi)[i][q].contract(B,dims1);
            }

            if constexpr ( vm::has_symm_v<context> )
            {
                em_fixed_size_matrix_t<nComponents1,NDim,value_type> sg( M_symm_grad[i][q].data() );
                em_fixed_size_cmatrix_t<nComponents1,NDim,value_type> g( gradiq.data() );
                sg = (g+g.transpose())/2;
            }

            // update divergence if needed
            if constexpr ( vm::has_div_v<context> )
            {
                M_div[i][q].setZero();
                if constexpr ( is_hdiv_conforming )
                {
                    for( int c = 0; c < nRealDim; ++c )
                        M_div[i][q]( 0,0 ) +=  (*M_gradphi)[i][q](c,c,0)/thegmc->J(q);
                }
                else if constexpr ( is_hcurl_conforming )
                    {
                        //M_div[i][0]( 0,0 ) =  ( Bt*((*M_gradphi)[i][0]*Bt.transpose()) ).trace();
                    }
                else
                {
                    for( int c = 0; c < nRealDim; ++c )
                        M_div[i][q]( 0,0 ) += gradiq(c,c,0);
                }
            }

            // update curl if needed
            if constexpr ( vm::has_curl_v<context> )
            {
                if constexpr ( NDim == 2 )
                {
                    M_curl[i][q]( 0 ) =  gradiq( 1,0,0 ) - gradiq( 0,1,0 );
                    //M_curl[i][q]( 1 ) =  M_curl[i][q]( 0 );
                    //M_curl[i][q]( 2 ) =  M_curl[i][q]( 0 );
                }
                else if constexpr ( NDim == 3 )
                {
                    M_curl[i][q]( 0 ) =  gradiq( 2,1,0 ) - gradiq( 1,2,0 );
                    M_curl[i][q]( 1 ) =  gradiq( 0,2,0 ) - gradiq( 2,0,0 );
                    M_curl[i][q]( 2 ) =  gradiq( 1,0,0 ) - gradiq( 0,1,0 );
                }
            }
        }
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateSymm( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = M_grad.shape()[0];
    for ( uint16_type q = 0; q < Q; ++q )
    {
        for ( uint16_type i = 0; i < I; ++i )
        {
            em_fixed_size_matrix_t<nComponents1,NDim,value_type> sg( M_symm_grad[i][q].data() );
            em_fixed_size_cmatrix_t<nComponents1,NDim,value_type> g( M_grad[i][q].data() );
            sg = (g+g.transpose())/2;
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateDiv( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = M_grad.shape()[0];
    for ( uint16_type q = 0; q < Q; ++q )
    {
        for ( uint16_type i = 0; i < I; ++i )
        {
            M_div[i][q].setZero();
            if constexpr ( is_hdiv_conforming )
            {
                for( int c = 0; c < nRealDim; ++c )
                    M_div[i][q]( 0,0 ) +=  (*M_gradphi)[i][q](c,c,0)/thegmc->J(q);
            }
            else if constexpr ( is_hcurl_conforming )
            {
                //M_div[i][0]( 0,0 ) =  ( Bt*((*M_gradphi)[i][0]*Bt.transpose()) ).trace();
            }
            else
            {
                for( int c = 0; c < nRealDim; ++c )
                    M_div[i][q]( 0,0 ) += M_grad[i][q](c,c,0);
            }
        }
    }
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateCurl( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = M_grad.shape()[0];
    for ( uint16_type q = 0; q < Q; ++q )
    {
        for ( uint16_type i = 0; i < I; ++i )
        {
            auto const& gradiq = M_grad[i][q];
            if constexpr ( NDim == 2 )
            {
                M_curl[i][q]( 0 ) =  gradiq( 1,0,0 ) - gradiq( 0,1,0 );
                //M_curl[i][q]( 1 ) =  M_curl[i][q]( 0 );
                //M_curl[i][q]( 2 ) =  M_curl[i][q]( 0 );
            }
            else if constexpr ( NDim == 3 )
            {
                M_curl[i][q]( 0 ) =  gradiq( 2,1,0 ) - gradiq( 1,2,0 );
                M_curl[i][q]( 1 ) =  gradiq( 0,2,0 ) - gradiq( 2,0,0 );
                M_curl[i][q]( 2 ) =  gradiq( 1,0,0 ) - gradiq( 0,1,0 );
            }
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateFirstDerivativeNormal( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = M_npoints;
    const uint16_type I = M_dn.shape()[0];

    for ( uint16_type i = 0; i < I; ++i )
        for ( uint16_type q = 0; q < Q; ++q )
        {
            M_dn[i][q].setZero();
            for ( uint16_type c1 = 0; c1 < gmc_type::NDim; ++c1 )
            {
                //M_dn[i][q]( c1,0 ) = 0;
                for ( uint16_type l = 0; l < gmc_type::NDim; ++l )
                {
                    //M_dn[i][q]( c1,0 ) += M_grad[i][q]( c1,l ) * thegmc->unitNormal( l, q );
                    M_dn[i][q]( c1,0 ) += this->grad(i,q)( c1,l,0 ) * thegmc->unitNormal( l, q );
                }
            }
        }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateHessian( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    precompute_type* __pc = M_pc.get().get();
    const uint16_type Q = do_optimization_p2?1:M_npoints;
    const uint16_type I = M_hessian.shape()[0];

    Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(2, 1)}};
    Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(1, 1)}};
    Eigen::array<dimpair_t, 1> dimsh = {{dimpair_t(0, 1)}};
    Eigen::array<ptrdiff_t, 3> myshuffle = { 2,0,1 };

    for ( uint16_type i = 0; i < I; ++i )
    {
        for ( uint16_type q = 0; q < Q; ++q )
        {
            tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> B ( thegmc->B( q ).data(),gmc_type::NDim, gmc_type::PDim );

            if constexpr ( !geometric_mapping_context_type::is_linear )//Geo_t::nOrder > 1 || !convex_type::is_simplex )
            {
                // gmc hessian NxPxP
                // pc hessian PxP
                // grad N
                // B NxP
                // hessian: N x N
                //tensor3_fixed_size_t<nComponents1,PDim,PDim,value_type> H0 = thegmc->hessian(q).contract( M_grad[i][q].chip(0,2), dimsh );
                tensor3_fixed_size_t<nComponents1,PDim,PDim,value_type> H0 = thegmc->hessian(q).contract( M_grad[i][q].chip(0,2), dimsh ).shuffle( myshuffle );
                //std::cout << "H0=" << H0 << std::endl;
                tensor3_fixed_size_t<nComponents1,PDim,PDim,value_type> H00 = __pc->hessian(i,q);
                //std::cout << "H00=" << H00 << std::endl;
                tensor3_fixed_size_t<nComponents1,PDim,PDim,value_type> H1 = H00 - H0;
                //std::cout << "H1=" << H1 << std::endl;
                M_hessian[i][q] = H1.contract( B, dims1 ).contract( B, dims2 );
                //std::cout << "H=" << M_hessian[i][q] << std::endl;
            }
            else
            {
                auto H1 = __pc->hessian(i,q).contract(B,dims1);
                M_hessian[i][q] = H1.contract( B, dims2 );
                //M_hessian[i][q] = (B.contract( H1, dims2)).shuffle( myperm );
            }

            //M_hessian[i][q] = B.contract( H1, dims2);
            if constexpr ( vm::has_laplacian_v<context> )
            {
                M_laplacian[i][q].setZero();
                for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < nRealDim; ++c2 )
                        M_laplacian[i][q](c1,0) += M_hessian[i][q](c1,c2,c2);
                } // c1 component of vector field
            }
        } // i
    } //q
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateLaplacian( geometric_mapping_context_type* thegmc, rank_t<1> )
{
    const uint16_type Q = do_optimization_p2?1:M_npoints;
    const uint16_type I = M_laplacian.shape()[0];

    for ( uint16_type i = 0; i < I; ++i )
    {
        for ( uint16_type q = 0; q < Q; ++q )
        {
            M_laplacian[i][q].setZero();
            for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
            {
                for ( uint16_type c2 = 0; c2 < nRealDim; ++c2 )
                    M_laplacian[i][q](c1,0) += M_hessian[i][q](c1,c2,c2);
            } // c1 component of vector field
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1> )
{
    geometric_mapping_context_type* thegmc = __gmc.get();
    const uint16_type I = M_phi.shape()[0];
    const int Qid = M_npoints;
    const int Q = do_optimization_p1?1:M_npoints;

    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
    if constexpr ( is_hdiv_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Qid; ++q )
            {
                tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> K ( thegmc->K( q ).data(), gmc_type::NDim, gmc_type::PDim );
                tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> Bt ( thegmc->B( q ).data(), gmc_type::NDim, gmc_type::PDim );

                // covariant piola transform
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
                M_phi[ii][q] = K.contract((*M_pc)->phi(ii,q), dims )/thegmc->J(q);
            }
        }
    }
    else if constexpr ( is_hcurl_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Qid; ++q )
            {
                tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> K( thegmc->K( q ).data(), gmc_type::NDim, gmc_type::PDim );
                tensor_map_fixed_size_matrix_t<gmc_type::NDim, gmc_type::PDim,value_type> Bt( thegmc->B( q ).data(), gmc_type::NDim, gmc_type::PDim );

                // piola transform
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
                M_phi[ii][q] = Bt.contract((*M_pc)->phi(ii,q), dims );
            }
        }
    }


    if constexpr ( gmc_type::is_on_face )
    {
        if constexpr ( vm::has_normal_component_v<context> )
        {
            this->updateNormalComponent( thegmc, rank_t<1>{} );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasNORMAL_COMPONENT( this->dynamicContext() ) )
                this->updateNormalComponent( thegmc, rank_t<1>{} );
        }
    }

    if constexpr ( vm::has_grad_v<context> || vm::has_first_derivative_v<context> ||
                   ( second_derivative_require_grad && ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ) ) )
    {
        this->updateGrad( thegmc, rank_t<1>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasGRAD( this->dynamicContext() ) || hasFIRST_DERIVATIVE( this->dynamicContext() ) ||
             ( second_derivative_require_grad && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
            this->updateGrad( thegmc, rank_t<1>{} );
    }

    if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( !vm::has_symm_v<context> && hasSYMM( this->dynamicContext() ) )
            this->updateSymm( thegmc, rank_t<1>{} );
        if ( !vm::has_div_v<context> && hasDIV( this->dynamicContext() ) )
            this->updateDiv( thegmc, rank_t<1>{} );
        if ( !vm::has_curl_v<context> && hasCURL( this->dynamicContext() ) )
            this->updateCurl( thegmc, rank_t<1>{} );
    }


    if constexpr ( vm::has_first_derivative_normal_v<context> )
    {
        this->updateFirstDerivativeNormal( thegmc, rank_t<1>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasFIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) )
            this->updateFirstDerivativeNormal( thegmc, rank_t<1>{} );
    }

    if constexpr ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> )
    {
        this->updateHessian( thegmc, rank_t<1>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) )
            this->updateHessian( thegmc, rank_t<1>{} );
    }

    if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( !vm::has_laplacian_v<context> && hasLAPLACIAN( this->dynamicContext() ) )
            this->updateLaplacian( thegmc, rank_t<1>{} );
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateNormalComponent( geometric_mapping_context_type* thegmc, rank_t<2> )
{
    const uint16_type Q = M_npoints;
    const uint16_type I = M_normal_component.shape()[0];
    for ( uint16_type i = 0; i < I; ++i )
        for ( uint16_type q = 0; q < Q; ++q )
        {
            M_normal_component[i][q].setZero();
            tensor_map_fixed_size_matrix_t<gmc_type::NDim,1,value_type> N ( thegmc->unitNormal( q ).data(), gmc_type::NDim,1 );
            for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                for( int c = 0; c < gmc_type::NDim; ++c )
                    M_normal_component[i][q]( c1,0 ) +=  M_phi[i][q]( c1, c ) * N( c );
        }
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateTrace( geometric_mapping_context_type* thegmc, rank_t<2> )
{
    const uint16_type Q = M_npoints;
    const uint16_type I = M_trace.shape()[0];
    //Eigen::array<dimpair_t, 1> dims = {{dimpair_t(0, 1)}};
    Eigen::Tensor<value_type, 0> res;
    for ( uint16_type i = 0; i < I; ++i )
        for ( uint16_type q = 0; q < Q; ++q )
        {
            res=M_phi[i][q].trace();
            M_trace[i][q]=res();
        }
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateGrad( geometric_mapping_context_type* thegmc, rank_t<2> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = M_grad.shape()[0];
    //typedef typename boost::multi_array<value_type,4>::index_range range;

    for ( uint16_type i = 0; i < I; ++i )
    {
        for ( uint16_type q = 0; q < Q; ++q )
        {
            tensor_map_fixed_size_matrix_t<gmc_type::NDim,gmc_type::PDim,value_type> B ( thegmc->B( q ).data(),gmc_type::NDim,gmc_type::PDim );
            Eigen::array<dimpair_t, 1> dims = {{dimpair_t(2, 1)}};
            auto & gradiq = M_grad[i][q];
            gradiq = (*M_gradphi)[i][q].contract( B,dims );
            // update divergence if needed
            if constexpr ( vm::has_div_v<context> )
            {
                auto & diviq = M_div[i][q];
                diviq.setZero();
                // div_i = sum_j \frac{\partial u_ij}{\partial x_j}
                for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                    for( uint16_type j = 0; j < nComponents2; ++j )
                        diviq( c1,0 ) +=  gradiq( c1, j, j );

                // div_j = sum_i \frac{\partial u_ij}{\partial x_i}
                // for( uint16_type j = 0; j < nComponents2; ++j )
                //     for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                //         M_div[i][q]( j,0 ) +=  M_grad[i][q]( c1, j, c1 );
            }
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
updateDiv( geometric_mapping_context_type* thegmc, rank_t<2> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;
    const uint16_type I = M_div.shape()[0];
    for ( uint16_type i = 0; i < I; ++i )
    {
        for ( uint16_type q = 0; q < Q; ++q )
        {
            auto const& gradiq = M_grad[i][q];
            auto & diviq = M_div[i][q];
            diviq.setZero();
            // div_i = sum_j \frac{\partial u_ij}{\partial x_j}
            for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                for( uint16_type j = 0; j < nComponents2; ++j )
                    diviq( c1,0 ) +=  gradiq( c1, j, j );
        }
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g, int SubEntityCoDim>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g,SubEntityCoDim>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<2> )
{
    geometric_mapping_context_type* thegmc = __gmc.get();

    if constexpr ( gmc_type::is_on_face )
    {
        if constexpr ( vm::has_normal_component_v<context> )
        {
            this->updateNormalComponent( thegmc, rank_t<2>{} );
        }
        else if constexpr ( vm::has_dynamic_basis_function_v<context> )
        {
            if ( hasNORMAL_COMPONENT( this->dynamicContext() ) )
                this->updateNormalComponent( thegmc, rank_t<2>{} );
        }
    }

    if constexpr ( vm::has_trace_v<context> )
    {
        this->updateTrace( thegmc, rank_t<2>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasTRACE( this->dynamicContext() ) )
            this->updateTrace( thegmc, rank_t<2>{} );
    }

    if constexpr ( vm::has_grad_v<context> || vm::has_first_derivative_v<context> ||
                   ( second_derivative_require_grad && ( vm::has_hessian_v<context> || vm::has_second_derivative_v<context> || vm::has_laplacian_v<context> ) ) )
    {
        this->updateGrad( thegmc, rank_t<2>{} );
    }
    else if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( hasGRAD( this->dynamicContext() ) || hasFIRST_DERIVATIVE( this->dynamicContext() ) ||
             ( second_derivative_require_grad && ( hasHESSIAN( this->dynamicContext() ) || hasSECOND_DERIVATIVE( this->dynamicContext() ) || hasLAPLACIAN( this->dynamicContext() ) ) ) )
            this->updateGrad( thegmc, rank_t<2>{} );
    }

    if constexpr ( vm::has_dynamic_basis_function_v<context> )
    {
        if ( !vm::has_div_v<context> && hasDIV( this->dynamicContext() ) )
            this->updateDiv( thegmc, rank_t<1>{} );
    }

#if 0
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value  || ( vm::has_dynamic_basis_function_v<context> && FIRST_DERIVATIVE_NORMAL( this->dynamicContext() ) ) )
        {
            const uint16_type I = nDof*nComponents1;
            const uint16_type Q = nPoints();

            for ( int i = 0; i < I; ++i )
                for ( uint16_type q = 0; q < Q; ++q )
                {
                    for ( uint16_type c1 = 0; c1 < NDim; ++c1 )
                    {
                        for ( uint16_type l = 0; l < NDim; ++l )
                        {
                            M_dn[i][q]( c1,0 ) += M_grad[i][q]( c1,l ) * thegmc->unitNormal( l, q );
                        }
                    }
                }
        }
#endif

}

} // Feel
