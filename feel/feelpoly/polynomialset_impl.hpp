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
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
resizeAndSet( rank_t<0> )
{
    Eigen::Tensor<value_type,3> i_grad( 1, nRealDim, 1 );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_first_derivative_normal<context>::value )
    {
        Eigen::Tensor<value_type,2> i_dn( 1,1 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        Eigen::Tensor<value_type,3> i_hessian( nRealDim, nRealDim, 1 );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
    }
    if ( vm::has_laplacian<context>::value || vm::has_second_derivative<context>::value  )
    {
        Eigen::Tensor<value_type,3> i_hessian( nRealDim, nRealDim, 1 );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
        Eigen::Tensor<value_type,2> i_lap( 1, 1 );
        std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
resizeAndSet( rank_t<1> )
{
    Eigen::Tensor<value_type,3> i_grad( nComponents1, nRealDim, 1 );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_first_derivative_normal<context>::value )
    {
        Eigen::Tensor<value_type,2> i_dn( nComponents1,1 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_div<context>::value )
    {
        Eigen::Tensor<value_type,2> i_div( 1, 1 );
        std::fill( M_div.data(), M_div.data()+M_div.num_elements(), i_div.constant(0.) );
    }

    if ( vm::has_curl<context>::value )
    {
        Eigen::Tensor<value_type,2> i_curl( nComponents1, 1 );
        std::fill( M_curl.data(), M_curl.data()+M_curl.num_elements(), i_curl.constant(0.) );
    }

    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        Eigen::Tensor<value_type,3> i_hessian( nComponents1, nRealDim, nRealDim );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
    }
    if ( vm::has_laplacian<context>::value || vm::has_second_derivative<context>::value  )
    {
        Eigen::Tensor<value_type,3> i_hessian( nComponents1, nRealDim, nRealDim );
        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), i_hessian.constant(0.) );
        Eigen::Tensor<value_type,2> i_lap( nComponents1, 1 );
        std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
resizeAndSet( rank_t<2> )
{
    Eigen::Tensor<value_type,3> i_grad( nComponents1, nComponents2, nRealDim );
    std::fill( M_grad.data(), M_grad.data()+M_grad.num_elements(), i_grad.constant(0.) );

    if ( vm::has_first_derivative_normal<context>::value )
    {
        Eigen::Tensor<value_type,2> i_dn( nComponents1, nComponents2 );
        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), i_dn.constant(0.) );
    }

    if ( vm::has_div<context>::value )
    {
        Eigen::Tensor<value_type,2> i_div( nComponents1, 1 );
        std::fill( M_div.data(), M_div.data()+M_div.num_elements(), i_div.constant(0.) );
    }

}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
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

                if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
                {
                    if ( do_optimization_p1 )
                        M_grad.resize( boost::extents[ntdof][1] );
                    else
                        M_grad.resize( boost::extents[ntdof][M_npoints] );

                        if ( vm::has_div<context>::value )
                        {
                            if ( do_optimization_p1 )
                                M_div.resize( boost::extents[ntdof][1] );
                            else
                                M_div.resize( boost::extents[ntdof][M_npoints] );
                        }

                        if ( vm::has_curl<context>::value )
                        {
                            if ( do_optimization_p1 )
                                M_curl.resize( boost::extents[ntdof][1] );
                            else
                                M_curl.resize( boost::extents[ntdof][M_npoints] );
                        }

                        if ( vm::has_first_derivative_normal<context>::value )
                        {
                                M_dn.resize( boost::extents[ntdof][M_npoints] );
                        }

                        if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
                        {
                                M_hessian.resize( boost::extents[ntdof][M_npoints] );
                        }
                        if ( vm::has_laplacian<context>::value || vm::has_second_derivative<context>::value  )
                        {
                            M_hessian.resize( boost::extents[ntdof][M_npoints] );
                            M_laplacian.resize( boost::extents[ntdof][M_npoints] );
                        }
                        resizeAndSet( rank_t<rank>() );
                }
        }

        M_phi = M_pc.get()->phi();
        M_gradphi = M_pc.get()->gradPtr();

        update( __gmc );
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0>, optimization_p1_t )
{
        //#pragma omp parallel
        //precompute_type* __pc = M_pc.get().get();
        geometric_mapping_context_type* thegmc = __gmc.get();

        if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
        {
                const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
                const uint16_type I = nDof; //M_ref_ele->nbDof();

                tensor_eigen_ublas_type B ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );

                for ( uint16_type i = 0; i < I; ++i )
                {
                    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 1)}};
                    M_grad[i][0] = (*M_gradphi)[i][0].contract( B,dims );
                }
#if 0
                // we need the normal derivative
                if ( vm::has_first_derivative_normal<context>::value )
                {
                        const uint16_type I = M_ref_ele->nbDof()*nComponents;
                        const uint16_type Q = nPoints();

                        for ( int i = 0; i < I; ++i )
                        {
                            for ( uint16_type l = 0; l < NDim; ++l )
                            {
                                value_type n = thegmc->unitNormal( l, 0 );
                                value_type gn = M_grad[i][0]( 0,l ) * n;
                                M_dn[i][0]( 0,0 ) += gn;
                            }
                        }
                }
#endif

        } // grad
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0>, no_optimization_p1_t )
{
    //#pragma omp parallel
    //precompute_type* __pc = M_pc.get().get();
    geometric_mapping_context_type* thegmc = __gmc.get();

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  ||
         vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        for ( uint16_type i = 0; i < I; ++i )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 1)}};
                M_grad[i][q] = (*M_gradphi)[i][q].contract( B,dims );
#if 0
                M_dx[i][q] = M_grad[i][q].col( 0 );

                if ( NDim == 2 )
                    M_dy[i][q] = M_grad[i][q].col( 1 );

                if ( NDim == 3 )
                    M_dz[i][q] = M_grad[i][q].col( 2 );
#endif
            }
        }
#if 0
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
        {
            const uint16_type I = M_ref_ele->nbDof()*nComponents;
            const uint16_type Q = nPoints();

            for ( int i = 0; i < I; ++i )
            {
                for ( uint16_type q = 0; q < Q; ++q )
                {
                    for ( uint16_type l = 0; l < NDim; ++l )
                    {
                        M_dn[i][q]( 0,0 ) += M_grad[i][q]( 0,l ) * thegmc->unitNormal( l, q );
                    }
                }
            }
        }
#endif
    } // grad

    //update( __gmc, rank_t<0>(), no_optimization_p1_t(), do_optimization_p2_t() );
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0>, no_optimization_p1_t, optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        // hessian only for P1 geometric mappings
        boost::multi_array<value_type,4> const& B3 = thegmc->B3();

        matrix_eigen_ublas_type B ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
        //#pragma omp for
        for ( uint16_type i = 0; i < I; ++i )
        {
            M_hessian[i][0] = B*__pc->hessian(i,0)*B.transpose();
            for( uint16_type q = 1; q < Q; ++q )
                M_hessian[i][q] = M_hessian[i][0];
            if ( vm::has_laplacian<context>::value  )
            {
                auto lap = M_hessian[i][0].trace();
                for( uint16_type q = 0; q < Q; ++q )
                    M_laplacian[i][q](0,0) = lap;

            }
        } // i
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0>, no_optimization_p1_t, no_optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {

        geometric_mapping_context_type* thegmc = __gmc.get();
        LOG(INFO) << "use no_optimization_p2_t element " << thegmc->id();
        precompute_type* __pc = M_pc.get().get();
        auto* __gmc_pc = thegmc->pc().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();
        hess_type L;
        for ( uint16_type q = 0; q < Q; ++q )
        {
            matrix_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
            for ( uint16_type i = 0; i < I; ++i )
            {
#if 0
                //L = __pc->hessian(i,q);
                //LOG(INFO) << "elt " << thegmc->id() << " L 1=" << L;
                for (int c1 = 0; c1 < gmc_type::NDim; ++c1 )
                {
                    //L -= __gmc_pc->hessian(i,c1,q)*M_grad[i][q](0,c1);
                    //LOG(INFO) << "elt " << thegmc->id() << " L c1 " << c1 << "="<< L;
                }
#endif
                M_hessian[i][q] = B*__pc->hessian(i,q)*B.transpose();
                if ( vm::has_laplacian<context>::value  )
                    M_laplacian[i][q](0,0) = M_hessian[i][q].trace();
            } // q
        } // i
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1>, no_optimization_p1_t )
{
        //precompute_type* __pc = M_pc.get().get();
        geometric_mapping_context_type* thegmc = __gmc.get();

        if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value )
        {
                const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
                const uint16_type I = nDof; //M_ref_ele->nbDof();

                typedef typename boost::multi_array<value_type,4>::index_range range;

                for ( uint16_type ii = 0; ii < I; ++ii )
                {
                       int ncomp= ( reference_element_type::is_product?nComponents1:1 );

                        for ( uint16_type c = 0; c < ncomp; ++c )
                        {
                                uint16_type i = I*c + ii;

                                for ( uint16_type q = 0; q < Q; ++q )
                                {
                                    tensor_eigen_ublas_type B ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );

                                    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(2, 1)}};
                                    M_grad[i][q] = (*M_gradphi)[i][q].contract( B,dims );

                                    // update divergence if needed
                                    if ( vm::has_div<context>::value )
                                    {
                                        M_div[i][q].setZero();
                                        for( int l = 0; l < nRealDim; ++l )
                                            M_div[i][q]( 0,0 ) +=  M_grad[i][q](l,l,0);
                                    }

                                    // update curl if needed
                                    if ( vm::has_curl<context>::value )
                                    {
                                        if ( NDim == 2 )
                                        {
                                            M_curl[i][q]( 0 ) =  M_grad[i][q]( 1,0,0 ) - M_grad[i][q]( 0,1,0 );
                                            M_curl[i][q]( 1 ) =  M_curl[i][q]( 0 );
                                            M_curl[i][q]( 2 ) =  M_curl[i][q]( 0 );
                                        }

                                        else if ( NDim == 3 )
                                        {
                                            M_curl[i][q]( 0 ) =  M_grad[i][q]( 2,1,0 ) - M_grad[i][q]( 1,2,0 );
                                            M_curl[i][q]( 1 ) =  M_grad[i][q]( 0,2,0 ) - M_grad[i][q]( 2,0,0 );
                                            M_curl[i][q]( 2 ) =  M_grad[i][q]( 1,0,0 ) - M_grad[i][q]( 0,1,0 );
                                        }
                                    }
                                }
                        }
                }
#if 0
                // we need the normal derivative
                if ( vm::has_first_derivative_normal<context>::value )
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
#if 0
        update( __gmc, rank_t<1>(), no_optimization_p1_t(), do_optimization_p2_t() );
#endif
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1>, no_optimization_p1_t, no_optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        auto* __gmc_pc = thegmc->pc().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        // hessian only for P1 geometric mappings
        boost::multi_array<value_type,4> const& B3 = thegmc->B3();
        hess_type L;
        //#pragma omp for
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            int ncomp= ( reference_element_type::is_product?nComponents1:1 );

            for ( uint16_type c = 0; c < ncomp; ++c )
            {
                uint16_type i = I*c + ii;
                for ( uint16_type q = 0; q < Q; ++q )
                {
                    matrix_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                    for ( uint16_type c1 = 0; c1 < NDim; ++c1 )
                    {
                        //L = __pc->hessian(i,c1,q);
                        //for (int c2 = 0; c2 < gmc_type::NDim; ++c2 )
                        //L -= __gmc_pc->hessian(i,c2,q)*M_grad[i][q](c1,c2);
                        L = B*__pc->hessian(i,c1,q)*B.transpose();
                        M_laplacian[i][q](c1,0) = L.trace();
                    } // c1 component of vector field

                } // q
            } // c
        } //ii
    }
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1>, no_optimization_p1_t, optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        // hessian only for P1 geometric mappings
        boost::multi_array<value_type,4> const& B3 = thegmc->B3();
        hess_type L;
        matrix_eigen_ublas_type B ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
        //#pragma omp for
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            int ncomp= ( reference_element_type::is_product?nComponents1:1 );

            for ( uint16_type c = 0; c < ncomp; ++c )
            {
                uint16_type i = I*c + ii;
                for ( uint16_type c1 = 0; c1 < NDim; ++c1 )
                {
                    L = B*__pc->hessian(i,c1,0)*B.transpose();
                    auto lap =  L.trace();
                    for( int q = 0; q < Q; ++q )
                        M_laplacian[i][q](c1,0) = lap;
                }// c1 component of vector field
            } // c
        } //ii
    }
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1>, optimization_p1_t )
{
    const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = nDof; //M_ref_ele->nbDof();

    //precompute_type* __pc = M_pc.get().get();
    geometric_mapping_context_type* thegmc = __gmc.get();

    tensor_eigen_ublas_type K ( thegmc->K( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
    tensor_eigen_ublas_type Bt ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
    if ( is_hdiv_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                // covariant piola transform
                M_phi[ii][q] = K.contract((*M_pc)->phi(ii,q), dims )/thegmc->J(q);
            }
        }
    }
    else if ( is_hcurl_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                // piola transform
                M_phi[ii][q] = Bt.contract((*M_pc)->phi(ii,q), dims );
            }
        }
    }

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
    {

        typedef typename boost::multi_array<value_type,4>::index_range range;

        //matrix_eigen_ublas_type Bt ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );

        matrix_eigen_grad_type grad_real = matrix_eigen_grad_type::Zero();
        //matrix_eigen_PN_type B=Bt.transpose();
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            int ncomp= ( reference_element_type::is_product?nComponents1:1 );

            for ( uint16_type c = 0; c < ncomp; ++c )
            {
                //std::cout << "component " << c << "\n";
                uint16_type i = I*c + ii;
                //uint16_type c1 = c;
                if ( is_hdiv_conforming )
                {
                    grad_real.noalias() = K*((*M_gradphi)[i][0]*Bt.transpose());
                    grad_real /= thegmc->J(0);
                }
                else if ( is_hcurl_conforming )
                {
                    grad_real.noalias() = Bt*((*M_gradphi)[i][0]*Bt.transpose());
                }
                else
                    grad_real.noalias() = (*M_gradphi)[i][0]*Bt.transpose();

                M_grad[i][0] = grad_real;

                // update divergence if needed
                if ( vm::has_div<context>::value )
                {
                    if ( is_hdiv_conforming )
                    {
                        M_div[i][0]( 0,0 ) =  (*M_gradphi)[i][0].trace()/thegmc->J(0);
                    }
                    else if ( is_hcurl_conforming )
                    {
                        M_div[i][0]( 0,0 ) =  ( Bt*((*M_gradphi)[i][0]*Bt.transpose()) ).trace();
                    }
                    else
                    {
                        M_div[i][0]( 0,0 ) =  M_grad[i][0].trace();
                    }
                }

                // update curl if needed
                if ( vm::has_curl<context>::value )
                {
                    if ( NDim == 2 )
                    {
#if 0
                        M_curl[i][0]( 0 ) = 0;
                        M_curl[i][0]( 1 ) = 0;
#else
                        M_curl[i][0]( 0 ) = M_grad[i][0]( 1,0 ) - M_grad[i][0]( 0,1 );
                        M_curl[i][0]( 1 ) = M_grad[i][0]( 1,0 ) - M_grad[i][0]( 0,1 );
#endif
                        M_curl[i][0]( 2 ) =  M_grad[i][0]( 1,0 ) - M_grad[i][0]( 0,1 );
                    }

                    else if ( NDim == 3 )
                    {
                        M_curl[i][0]( 0 ) =  M_grad[i][0]( 2,1 ) - M_grad[i][0]( 1,2 );
                        M_curl[i][0]( 1 ) =  M_grad[i][0]( 0,2 ) - M_grad[i][0]( 2,0 );
                        M_curl[i][0]( 2 ) =  M_grad[i][0]( 1,0 ) - M_grad[i][0]( 0,1 );
                    }
                }
            } // c
        } // ii

        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
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

    } // grad
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<2>, optimization_p1_t )
{
    const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = nDof; //M_ref_ele->nbDof();

    //precompute_type* __pc = M_pc.get().get();
    geometric_mapping_context_type* thegmc = __gmc.get();

    tensor_eigen_ublas_type B ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
    {

        typedef typename boost::multi_array<value_type,4>::index_range range;

        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                uint16_type i = I*c + ii;

                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(2, 1)}};
                M_grad[i][0] = (*M_gradphi)[i][0].contract( B,dims );
                // update divergence if needed
                if ( vm::has_div<context>::value )
                {
                    M_div[i][0].setZero();
                    for( uint16_type j = 0; j < nComponents2; ++j )
                        for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                            M_div[i][0]( j,0 ) +=  M_grad[i][0]( c1, j, c1 );
                }
            } // c
        } // ii
#if 0
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
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
    } // grad

}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<2>, no_optimization_p1_t )
{
    const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = nDof; //M_ref_ele->nbDof();

    //precompute_type* __pc = M_pc.get().get();
    geometric_mapping_context_type* thegmc = __gmc.get();



    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
    {

        typedef typename boost::multi_array<value_type,4>::index_range range;

        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                uint16_type i = I*c + ii;

                for ( uint16_type q = 0; q < Q; ++q )
                {
                    tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(2, 1)}};
                    M_grad[i][q] = (*M_gradphi)[i][q].contract( B,dims );
                    // update divergence if needed
                    if ( vm::has_div<context>::value )
                    {
                        M_div[i][q].setZero();
                        for( uint16_type j = 0; j < nComponents2; ++j )
                            for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                                M_div[i][q]( j,0 ) +=  M_grad[i][q]( c1, j, c1 );
                    }
                }
            } // c
        } // ii
#if 0
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
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
    } // grad

}

} // Feel
