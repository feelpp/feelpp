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
                                // in P0 and P1 hessian is 0
                                std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), hess_type::Zero() );
                        }
                        if ( vm::has_laplacian<context>::value || vm::has_second_derivative<context>::value  )
                        {
                            M_hessian.resize( boost::extents[ntdof][M_npoints] );
                            if ( do_optimization_p1 )
                            {
                                // in P0 and P1 hessian is 0
                                std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), hess_type::Zero() );
                            }

                            M_laplacian.resize( boost::extents[ntdof][M_npoints] );
                            if ( do_optimization_p1  )
                            {
                                // in P0 and P1 laplacian is 0
                                std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), laplacian_type::Zero() );
                            }
                        }
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, optimization_p1_t )
{
        //#pragma omp parallel
        //precompute_type* __pc = M_pc.get().get();
        geometric_mapping_context_type* thegmc = __gmc.get();

        if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
        {
                const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
                const uint16_type I = nDof; //M_ref_ele->nbDof();


                matrix_eigen_ublas_type Bt ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                matrix_eigen_grad_type grad_real = matrix_eigen_grad_type::Zero();

                for ( uint16_type i = 0; i < I; ++i )
                {
                    grad_real.noalias() = (*M_gradphi)[i][0]*Bt.transpose();
                    M_grad[i][0] = grad_real;
                }

                // we need the normal derivative
                if ( vm::has_first_derivative_normal<context>::value )
                {
                        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), dn_type::Zero() );
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

        } // grad
}

template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, no_optimization_p1_t )
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
                                matrix_eigen_ublas_type Bt ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                                M_grad[i][q] = (*M_gradphi)[i][q]*Bt.transpose();
#if 0
                                M_dx[i][q] = M_grad[i][q].col( 0 );

                                if ( NDim == 2 )
                                        M_dy[i][q] = M_grad[i][q].col( 1 );

                                if ( NDim == 3 )
                                        M_dz[i][q] = M_grad[i][q].col( 2 );
#endif
                        }
                }

                // we need the normal derivative
                if ( vm::has_first_derivative_normal<context>::value )
                {
                        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), dn_type::Zero() );
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

        } // grad

        update( __gmc, mpl::int_<0>(), no_optimization_p1_t(), do_optimization_p2_t() );
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, no_optimization_p1_t, optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        // hessian only for P1 geometric mappings
        boost::multi_array<value_type,4> const& B3 = thegmc->B3();

        std::fill( M_hessian.data(), M_hessian.data()+M_hessian.num_elements(), hess_type::Zero() );
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, no_optimization_p1_t, no_optimization_p2_t )
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, no_optimization_p1_t )
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
                                        //uint16_type c1 = c;
                                        matrix_eigen_ublas_type Bt ( thegmc->B( q ).data().begin(), gmc_type::NDim, gmc_type::PDim );
                                        //matrix_type const& Bq = thegmc->B( q );
                                        M_grad[i][q] = (*M_gradphi)[i][q]*Bt.transpose();
#if 0
                                        M_dx[i][q] = M_grad[i][q].col( 0 );

                                        if ( NDim == 2 )
                                                M_dy[i][q] = M_grad[i][q].col( 1 );

                                        if ( NDim == 3 )
                                                M_dz[i][q] = M_grad[i][q].col( 2 );
#endif
                                        // update divergence if needed
                                        if ( vm::has_div<context>::value )
                                        {
                                                M_div[i][q]( 0,0 ) =  M_grad[i][q].trace();
                                        }

                                        // update curl if needed
                                        if ( vm::has_curl<context>::value )
                                        {
                                                if ( NDim == 2 )
                                                {
                                                        M_curl[i][q]( 0 ) =  M_grad[i][q]( 1,0 ) - M_grad[i][q]( 0,1 );
                                                        M_curl[i][q]( 1 ) =  M_curl[i][q]( 0 );
                                                        M_curl[i][q]( 2 ) =  M_curl[i][q]( 0 );
                                                }

                                                else if ( NDim == 3 )
                                                {
                                                        M_curl[i][q]( 0 ) =  M_grad[i][q]( 2,1 ) - M_grad[i][q]( 1,2 );
                                                        M_curl[i][q]( 1 ) =  M_grad[i][q]( 0,2 ) - M_grad[i][q]( 2,0 );
                                                        M_curl[i][q]( 2 ) =  M_grad[i][q]( 1,0 ) - M_grad[i][q]( 0,1 );
                                                }
                                        }
                                }
                        }
                }

                // we need the normal derivative
                if ( vm::has_first_derivative_normal<context>::value )
                {
                        std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), dn_type::Zero() );
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
        }
        update( __gmc, mpl::int_<1>(), no_optimization_p1_t(), do_optimization_p2_t() );
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, no_optimization_p1_t, no_optimization_p2_t )
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
        std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), laplacian_type::Zero() );
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, no_optimization_p1_t, optimization_p2_t )
{
    if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value || vm::has_laplacian<context>::value  )
    {
        geometric_mapping_context_type* thegmc = __gmc.get();
        precompute_type* __pc = M_pc.get().get();
        const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
        const uint16_type I = nDof; //M_ref_ele->nbDof();

        // hessian only for P1 geometric mappings
        boost::multi_array<value_type,4> const& B3 = thegmc->B3();
        std::fill( M_laplacian.data(), M_laplacian.data()+M_laplacian.num_elements(), laplacian_type::Zero() );
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, optimization_p1_t )
{
    const uint16_type Q = M_npoints;//__gmc->nPoints();//M_grad.size2();
    const uint16_type I = nDof; //M_ref_ele->nbDof();

    //precompute_type* __pc = M_pc.get().get();
    geometric_mapping_context_type* thegmc = __gmc.get();

    matrix_eigen_ublas_type K ( thegmc->K( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
    matrix_eigen_ublas_type Bt ( thegmc->B( 0 ).data().begin(), gmc_type::NDim, gmc_type::PDim );
    if ( is_hdiv_conforming )
    {
        for ( uint16_type ii = 0; ii < I; ++ii )
        {
            for ( uint16_type q = 0; q < Q; ++q )
            {
                // covariant piola transform
                M_phi[ii][q].noalias() = K*(*M_pc)->phi(ii,q);
                M_phi[ii][q] /= thegmc->J(q);
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
                M_phi[ii][q].noalias() = Bt*(*M_pc)->phi(ii,q);
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
            std::fill( M_dn.data(), M_dn.data()+M_dn.num_elements(), dn_type::Zero() );
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
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<2>, mpl::bool_<true> )
{
}
template<typename Poly, template<uint16_type> class PolySetType>
template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g>
void
PolynomialSet<Poly,PolySetType>::Context<context_v, Basis_t,Geo_t,ElementType,context_g>::
update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<2>, mpl::bool_<false> )
{
}

} // Feel
