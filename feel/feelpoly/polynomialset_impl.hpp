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
        if ( vm::has_laplacian<context>::value )
        {
            Eigen::Tensor<value_type,2> i_lap( 1, 1 );
            std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
        }
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
        if ( vm::has_laplacian<context>::value || vm::has_second_derivative<context>::value  )
        {
            Eigen::Tensor<value_type,2> i_lap( nComponents1, 1 );
            std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), i_lap.constant(0.) );
        }
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
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0> )
{
    geometric_mapping_context_type* thegmc = __gmc.get();

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  ||
         vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        const uint16_type Q = do_optimization_p1?1:M_npoints;
        const uint16_type I = nDof;
        Eigen::array<int, 3> tensorGradShapeAfterContract{{1, 1, nRealDim}};
        Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 1)}};
        for ( uint16_type i = 0; i < I; ++i )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                // grad = (gradphi_1,...,gradphi_nRealDim) * B^T
                M_grad[i][q].reshape( tensorGradShapeAfterContract ) = ((*M_gradphi)[i][q].contract( B,dims ));
#if 0
                M_dx[i][q] = M_grad[i][q].col( 0 );

                if ( NDim == 2 )
                    M_dy[i][q] = M_grad[i][q].col( 1 );

                if ( NDim == 3 )
                    M_dz[i][q] = M_grad[i][q].col( 2 );
#endif
            }
        }
#if 1
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
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
#endif
    } // grad
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        auto* __gmc_pc = thegmc->pc().get();
        const uint16_type Q = do_optimization_p2?1:M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();
        hess_type L;
        Eigen::array<int, 3> tensorHessShapeAfterContract{{nRealDim, 1, nRealDim}};
        Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(1, 1)}};
        Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(1, 0)}};
        for ( uint16_type q = 0; q < Q; ++q )
        {
            tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
            for ( uint16_type i = 0; i < I; ++i )
            {
                //M_hessian[i][q] = B.contract( __pc->hessian(i,q).contract(B,dims1), dims2);
                auto H1 = __pc->hessian(i,q).contract(B,dims1);
                M_hessian[i][q].reshape( tensorHessShapeAfterContract ) = B.contract( H1, dims2 );

                if ( vm::has_laplacian<context>::value  )
                {
                    M_laplacian[i][q].setZero();
                    for( int c = 0; c < nRealDim; ++c )
                        M_laplacian[i][q](0,0) += M_hessian[i][q]( c, c, 0 );
                        //M_laplacian[i][q](0,0) = M_hessian[i][q].trace();
                }
            } // q
        } // i
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1> )
{
    geometric_mapping_context_type* thegmc = __gmc.get();
    const uint16_type I = M_phi.shape()[0];
    const int Qid = M_npoints;
    const int Q = do_optimization_p1?1:M_npoints;

    Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
    if ( is_hdiv_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Qid; ++q )
            {
                tensor_eigen_ublas_type K ( thegmc->K( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                tensor_eigen_ublas_type Bt ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );

                // covariant piola transform
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
                M_phi[ii][q] = K.contract((*M_pc)->phi(ii,q), dims )/thegmc->J(q);
            }
        }
    }
    else if ( is_hcurl_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Qid; ++q )
            {
                tensor_eigen_ublas_type K ( thegmc->K( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                tensor_eigen_ublas_type Bt ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );

                // piola transform
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};
                M_phi[ii][q] = Bt.contract((*M_pc)->phi(ii,q), dims );
            }
        }
    }

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value )
    {
        const uint16_type Q = do_optimization_p1?1:M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = M_grad.shape()[0];

        typedef typename boost::multi_array<value_type,4>::index_range range;
        Eigen::array<int, 3> tensorGradShapeAfterContract{{nComponents1, 1, nRealDim}};
        Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(1, 1)}};
        Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(1, 0)}};

        for ( uint16_type q = 0; q < Q; ++q )
        {
            tensor_eigen_ublas_type K ( thegmc->K( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
            tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
            for ( uint16_type i = 0; i < I; ++i )
            {
                //M_grad[i][q] = (*M_gradphi)[i][q].contract( B,dims );

                if ( is_hdiv_conforming )
                {
                    auto v = (*M_gradphi)[i][q].contract(B,dims1);
                    M_grad[i][q] =  K.contract(v,dims2)/thegmc->J(q);
                }
                else if ( is_hcurl_conforming )
                {
                    //auto v = (*M_gradphi)[i][q].contract(B,dims1);
                    M_grad[i][q] =  B.contract((*M_gradphi)[i][q].contract(B,dims1),dims2);
                }
                else
                    M_grad[i][q].reshape( tensorGradShapeAfterContract ) = (*M_gradphi)[i][q].contract(B,dims1);

                // update divergence if needed
                if ( vm::has_div<context>::value )
                {
                    M_div[i][q].setZero();
                    if ( is_hdiv_conforming )
                    {
                        for( int c = 0; c < nRealDim; ++c )
                            M_div[i][q]( 0,0 ) +=  (*M_gradphi)[i][q](c,c,0)/thegmc->J(q);
                    }
                    else if ( is_hcurl_conforming )
                    {
                        //M_div[i][0]( 0,0 ) =  ( Bt*((*M_gradphi)[i][0]*Bt.transpose()) ).trace();
                    }
                    else
                    {
                        for( int c = 0; c < nRealDim; ++c )
                            M_div[i][q]( 0,0 ) +=  M_grad[i][q](c,c,0);
                    }
                }

                // update curl if needed
                if ( vm::has_curl<context>::value )
                {
                    if ( NDim == 2 )
                    {
                        M_curl[i][q]( 0 ) =  M_grad[i][q]( 1,0,0 ) - M_grad[i][q]( 0,1,0 );
                        M_curl[i][q]( 1 ) =  M_curl[i][q]( 0 );
                        //M_curl[i][q]( 2 ) =  M_curl[i][q]( 0 );
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
#if 1
        // we need the normal derivative
        if ( vm::has_first_derivative_normal<context>::value )
        {
            // const uint16_type I = nDof*nComponents1;
            // const uint16_type Q = nPoints();
            const uint16_type Q = M_npoints;//do_optimization_p1?1:M_npoints;
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
#endif
    }
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        auto* __gmc_pc = thegmc->pc().get();
        const uint16_type Q = do_optimization_p2?1:M_npoints;
        const uint16_type I = M_laplacian.shape()[0];

        hess_type L;
        for ( uint16_type i = 0; i < I; ++i )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                Eigen::array<dimpair_t, 1> dims1 = {{dimpair_t(2, 1)}};
                Eigen::array<dimpair_t, 1> dims2 = {{dimpair_t(1, 1)}};
                auto H1 = __pc->hessian(i,q).contract(B,dims1);
                Eigen::array<int, 3> myperm = {{1, 0, 2}};
                M_hessian[i][q] = (B.contract( H1, dims2)).shuffle( myperm );
                //M_hessian[i][q] = B.contract( H1, dims2);

                M_laplacian[i][q].setZero();
                for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < nRealDim; ++c2 )
                        M_laplacian[i][q](c1,0) += M_hessian[i][q](c1,c2,c2);
                } // c1 component of vector field
            } // i
        } //q
    }
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, rank_t<2> )
{
    const uint16_type Q = do_optimization_p1?1:M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = M_grad.shape()[0];
    geometric_mapping_context_type* thegmc = __gmc.get();

    if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
    {

        typedef typename boost::multi_array<value_type,4>::index_range range;

        for ( uint16_type i = 0; i < I; ++i )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                tensor_eigen_ublas_type B ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(2, 1)}};
                M_grad[i][q] = (*M_gradphi)[i][q].contract( B,dims );
                // update divergence if needed
                if ( vm::has_div<context>::value )
                {
                    M_div[i][q].setZero();
                    // div_i = sum_j \frac{\partial u_ij}{\partial x_j}
                    for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                        for( uint16_type j = 0; j < nComponents2; ++j )
                            M_div[i][q]( c1,0 ) +=  M_grad[i][q]( c1, j, j );

                    // div_j = sum_i \frac{\partial u_ij}{\partial x_i}
                    // for( uint16_type j = 0; j < nComponents2; ++j )
                    //     for( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                    //         M_div[i][q]( j,0 ) +=  M_grad[i][q]( c1, j, c1 );
                }
            }
        } // i
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
