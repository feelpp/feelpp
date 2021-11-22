/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 2013-11-19

  Copyright (C) 2012 Universite Joseph Fourier (Grenoble I)

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
   \file solidmecresidual.hpp
   \author Vincent Chabannes
   \date 2013-11-19
 */
#ifndef __FEELPP_MODELS_VF_SOLIDMECSTVENANTKIRCHHOFF_H
#define __FEELPP_MODELS_VF_SOLIDMECSTVENANTKIRCHHOFF_H 1

namespace Feel
{
namespace FeelModels
{
/// \cond detail
/**
 * \class Det
 * \brief det of a matrix
 *
 * @author Vincent Chabannes
 * @see
 */
template<typename ElementType , typename ElementLameCoeffType, typename SpecificExprType>
class StressStVenantKirchhoff
{
public:

    /** @name Typedefs
     */
    //@{

    typedef StressStVenantKirchhoff<ElementType,ElementLameCoeffType,SpecificExprType> this_type;

    static const size_type context = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_lamecoeff = vm::JACOBIAN;

    typedef ElementType element_type;
    typedef ElementLameCoeffType element_lamecoeff_type;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename functionspace_type::reference_element_type* fe_ptrtype;
    typedef typename functionspace_type::reference_element_type fe_type;
    typedef typename functionspace_type::geoelement_type geoelement_type;
    typedef typename functionspace_type::gm_type gm_type;
    typedef typename functionspace_type::value_type value_type;
    typedef value_type evaluate_type;


    static const uint16_type rank = fe_type::rank;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type orderdisplacement = functionspace_type::basis_type::nOrder;
    //------------------------------------------------------------------------------//
    // lamecoeff functionspace
    typedef typename element_lamecoeff_type::functionspace_type functionspace_lamecoeff_type;
    typedef typename functionspace_lamecoeff_type::reference_element_type fe_lamecoeff_type;

    //------------------------------------------------------------------------------//

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<0> >,
                                            mpl::bool_<false>,
                                            typename mpl::if_<boost::is_same<Func,fe_type>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type >::type::value;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<2> >,
                                            typename mpl::if_<boost::is_same<Func,fe_type>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type,
                                            mpl::bool_<false> >::type::value;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    StressStVenantKirchhoff( element_type const & v, element_lamecoeff_type const& lambda, element_lamecoeff_type const& mu)
        :
        M_v( boost::cref(v) ),
        M_coeffLame1( boost::cref(lambda) ),
        M_coeffLame2( boost::cref(mu) )
    {}

    StressStVenantKirchhoff( StressStVenantKirchhoff const & op )
        :
        M_v( op.M_v ),
        M_coeffLame1( op.M_coeffLame1 ),
        M_coeffLame2( op.M_coeffLame2 )
    {}

    ~StressStVenantKirchhoff()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! polynomial order
    uint16_type polynomialOrder() const { return (SpecificExprType::value == 0)? (3*(orderdisplacement-1)) : (4*(orderdisplacement-1)); }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    element_type const& e() const { return M_v; }
    element_lamecoeff_type const& coeffLame1() const { return M_coeffLame1; }
    element_lamecoeff_type const& coeffLame2() const { return M_coeffLame2; }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename element_type::value_type value_type;
        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        // fe context
        typedef typename fe_type::PreCompute pc_type;
        typedef std::shared_ptr<pc_type> pc_ptrtype;
        typedef typename fe_type::template Context<context, fe_type, gm_type,geoelement_type,gmc_type::context,0, gmc_type::subEntityCoDim> ctx_type;
        typedef std::shared_ptr<ctx_type> ctx_ptrtype;
        // fe lamecoeff context
        typedef typename fe_lamecoeff_type::PreCompute pc_lamecoeff_type;
        typedef std::shared_ptr<pc_lamecoeff_type> pc_lamecoeff_ptrtype;
        typedef typename fe_lamecoeff_type::template Context<context_lamecoeff, fe_lamecoeff_type, gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_lamecoeff_type;
        typedef std::shared_ptr<ctx_lamecoeff_type> ctx_lamecoeff_ptrtype;


        //----------------------------------------------------------------------------------------------------//

        /*typedef Shape<gmc_type::nDim, Vectorial, false, false> shape_id;
        typedef Eigen::Matrix<value_type,shape_id::M,shape_id::N> loc_id_type;
        typedef boost::multi_array<loc_id_type,1> array_id_type;*/

        typedef Shape<gmc_type::nDim, Tensor2, false, false> shape_grad;
        typedef Eigen::Matrix<value_type,shape_grad::M,shape_grad::N> loc_grad_type;
        typedef boost::multi_array<loc_grad_type,1> array_grad_type;

        typedef Shape<gmc_type::nDim, Scalar, false, false> shape_scalar;
        typedef Eigen::Matrix<value_type,shape_scalar::M,shape_scalar::N> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type,1> array_scalar_type;

        typedef typename mpl::if_<boost::is_same<SpecificExprType,mpl::int_<0> >,
                                  shape_grad,
                                  shape_scalar >::type shape;


        typedef loc_grad_type loc_type;
        typedef array_grad_type array_type;

        //typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> locAssembly_LinearForm_type;

        //----------------------------------------------------------------------------------------------------//
        typedef Basis_i_t map_basis_fec_test_type;
        typedef Basis_j_t map_basis_fec_trial_type;
        //typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_test_type,vf::detail::gmc<0> >,
        //                        mpl::identity<vf::detail::gmc<0> >,
        //                         mpl::identity<vf::detail::gmc<1> > >::type::type key_fec_test_type;
        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_test_type,vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<1> > >::type::type basis_fec_test_key_type;

        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_trial_type,vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<1> > >::type::type basis_fec_trial_key_type;

        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type  basis_fec_test_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type* basis_fec_test_ptrtype;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type  basis_fec_trial_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type* basis_fec_trial_ptrtype;



        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_fecTest( fusion::at_key<basis_fec_test_key_type>( fev ).get() ),
            M_fecTrial( fusion::at_key<basis_fec_trial_key_type>( feu ).get() ),
            M_pc( new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.e().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_pcLameCoeff( new pc_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxLameCoeff( new ctx_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_lamecoeff_ptrtype const&)M_pcLameCoeff ) ),
            M_loc( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
        //M_locId( expr.e().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff1( expr.coeffLame1().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff2( expr.coeffLame2().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_zero( ret_type::Zero() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_fecTest( fusion::at_key<basis_fec_test_key_type>( fev ).get() ),
            M_pc( new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.e().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_pcLameCoeff( new pc_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxLameCoeff( new ctx_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_lamecoeff_ptrtype const&)M_pcLameCoeff ) ),
            M_loc( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
            //M_locId( expr.e().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff1( expr.coeffLame1().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff2( expr.coeffLame2().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_zero( ret_type::Zero() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pc( new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.e().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_pcLameCoeff( new pc_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxLameCoeff( new ctx_lamecoeff_type( expr.coeffLame1().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_lamecoeff_ptrtype const&)M_pcLameCoeff ) ),
            M_loc( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
            //M_locId( expr.e().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.e().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff1( expr.coeffLame1().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locLameCoeff2( expr.coeffLame2().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_zero( ret_type::Zero() )
        {}

        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom )
        {
            this->setGmc( geom );
            std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_type::Zero() );
            //std::fill( M_locId.data(), M_locId.data()+M_locId.num_elements(), loc_id_type::Zero() );
            std::fill( M_locGrad.data(), M_locGrad.data()+M_locGrad.num_elements(), loc_grad_type::Zero() );

            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            {
                M_pc->update( this->gmc()->pc()->nodes() );
                M_pcLameCoeff->update( this->gmc()->pc()->nodes() );
            }
            M_ctx->update( fusion::at_key<key_type>( geom ),  (pc_ptrtype const&) M_pc );
            M_ctxLameCoeff->update( fusion::at_key<key_type>( geom ),  (pc_lamecoeff_ptrtype const&) M_pcLameCoeff );
#if 1
            //M_expr.e().id( *M_ctx, M_locId );
            M_expr.e().grad( *M_ctx, M_locGrad );

            std::fill( M_locLameCoeff1.data(), M_locLameCoeff1.data()+M_locLameCoeff1.num_elements(), loc_scalar_type::Zero() );
            std::fill( M_locLameCoeff2.data(), M_locLameCoeff2.data()+M_locLameCoeff2.num_elements(), loc_scalar_type::Zero() );
            M_expr.coeffLame1().id( *M_ctxLameCoeff, M_locLameCoeff1 );
            M_expr.coeffLame2().id( *M_ctxLameCoeff, M_locLameCoeff2 );
#else

            const size_type elt_id = fusion::at_key<key_type>( geom )->id();
            for ( int l = 0; l < functionspace_type::basis_type::nDof; ++l )
                {
                    const int ncdof = functionspace_type::is_product?nComponents1:1;
                    for ( typename array_type::index c1 = 0; c1 < ncdof; ++c1 )
                        {
                            typename array_type::index ldof = functionspace_type::basis_type::nDof*c1+l;
                            const size_type gdof = M_expr.e().functionSpace()->dof()->localToGlobal( elt_id, l, c1 ).index();
                            const value_type v_ = M_expr.e().globalValue( gdof );

                            for ( uint16_type q = 0; q < M_geot->nPoints(); ++q )
                                {
                                    for ( typename array_type::index i = 0; i < nComponents1; ++i )
                                        {
                                            M_locId[q]( i,0 ) += v_*M_ctx->id( ldof, i, 0, q );

                                            for ( uint16_type j = 0; j < gmc_type::nDim/*nRealDim*/; ++j )
                                                {
                                                    M_locGrad[q]( i,j ) += v_*M_ctx->grad( ldof, i, j, q );
                                                }
                                        }
                                }
                        }
                }
#endif
            update(mpl::int_<gmc_type::nDim>(), mpl::int_<SpecificExprType::value>() );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalijq( i,j,c1,c2,q,mpl::int_<gmc_type::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ,mpl::int_<2> /*Dim*/ ) const
        {
            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q );

            const value_type du1dx = M_fecTest->grad( i, 0, 0, q ), du1dy = M_fecTest->grad( i, 0, 1, q );
            const value_type du2dx = M_fecTest->grad( i, 1, 0, q ), du2dy = M_fecTest->grad( i, 1, 1, q );

            // dFSv = dF*val(Sv)
            const value_type dFSv11 = du1tdx*M_loc[q](0,0) + du1tdy*M_loc[q](1,0);
            const value_type dFSv12 = du1tdx*M_loc[q](0,1) + du1tdy*M_loc[q](1,1);
            const value_type dFSv21 = du2tdx*M_loc[q](0,0) + du2tdy*M_loc[q](1,0);
            const value_type dFSv22 = du2tdx*M_loc[q](0,1) + du2tdy*M_loc[q](1,1);

            //auto dE = alpha_f*sym(gradt(u)) + alpha_f*alpha_f*0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            const value_type dE11 = du1tdx + /*0.5**/(du1vdx*du1tdx + du2vdx*du2tdx);
            const value_type dE12 = 0.5*( du1tdy + du2tdx + du1tdx*du1vdy + du2tdx*du2vdy + du1vdx*du1tdy + du2vdx*du2tdy );
            const value_type dE21 = 0.5*( du2tdx + du1tdy + du1tdy*du1vdx + du2tdy*du2vdx + du1vdy*du1tdx + du2vdy*du2tdx );
            const value_type dE22 = du2tdy + du1vdy*du1tdy + du2vdy*du2tdy;

            const value_type tracedE = dE11 + dE22;
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#else
            const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);
#endif

            // auto dS = idv(M_P0Coefflame1)*trace(dE)*Id + 2*idv(M_P0Coefflame2)*dE;
            const value_type dS11 = lame1*tracedE + 2*lame2*dE11;
            const value_type dS12 = 2*lame2*dE12;
            const value_type dS21 = 2*lame2*dE21;
            const value_type dS22 = lame1*tracedE + 2*lame2*dE22;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = (1.+du1vdx)*dS11 + du1vdy*dS21;
            const value_type FvdS12 = (1.+du1vdx)*dS12 + du1vdy*dS22;
            const value_type FvdS21 = du2vdx*dS11 + (1.+du2vdy)*dS21;
            const value_type FvdS22 = du2vdx*dS12 + (1.+du2vdy)*dS22;

            // dF*val(Sv) + val(Fv)*dS)*trans(grad(v))
            const value_type resTrial11 = dFSv11 + FvdS11;
            const value_type resTrial12 = dFSv12 + FvdS12;
            const value_type resTrial21 = dFSv21 + FvdS21;
            const value_type resTrial22 = dFSv22 + FvdS22;

            // trace( (dF*val(Sv) + val(Fv)*dS)*trans(grad(v)) )
            const value_type res11 = resTrial11*du1dx + resTrial12*du1dy;
            const value_type res22 = resTrial21*du2dx + resTrial22*du2dy;
            return res11+res22;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ,mpl::int_<3> /*Dim*/ ) const
        {
            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
            const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q ), du1tdz=M_fecTrial->grad( j, 0, 2, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q ), du2tdz=M_fecTrial->grad( j, 1, 2, q );
            const value_type du3tdx = M_fecTrial->grad( j, 2, 0, q ), du3tdy = M_fecTrial->grad( j, 2, 1, q ), du3tdz=M_fecTrial->grad( j, 2, 2, q );

            const value_type du1dx = M_fecTest->grad( i, 0, 0, q ), du1dy = M_fecTest->grad( i, 0, 1, q ), du1dz=M_fecTest->grad( i, 0, 2, q );
            const value_type du2dx = M_fecTest->grad( i, 1, 0, q ), du2dy = M_fecTest->grad( i, 1, 1, q ), du2dz=M_fecTest->grad( i, 1, 2, q );
            const value_type du3dx = M_fecTest->grad( i, 2, 0, q ), du3dy = M_fecTest->grad( i, 2, 1, q ), du3dz=M_fecTest->grad( i, 2, 2, q );

            // dFSv = dF*val(Sv)
            const value_type dFSv11 = du1tdx*M_loc[q](0,0) + du1tdy*M_loc[q](1,0) + du1tdz*M_loc[q](2,0);
            const value_type dFSv12 = du1tdx*M_loc[q](0,1) + du1tdy*M_loc[q](1,1) + du1tdz*M_loc[q](2,1);
            const value_type dFSv13 = du1tdx*M_loc[q](0,2) + du1tdy*M_loc[q](1,2) + du1tdz*M_loc[q](2,2);

            const value_type dFSv21 = du2tdx*M_loc[q](0,0) + du2tdy*M_loc[q](1,0) + du2tdz*M_loc[q](2,0);
            const value_type dFSv22 = du2tdx*M_loc[q](0,1) + du2tdy*M_loc[q](1,1) + du2tdz*M_loc[q](2,1);
            const value_type dFSv23 = du2tdx*M_loc[q](0,2) + du2tdy*M_loc[q](1,2) + du2tdz*M_loc[q](2,2);

            const value_type dFSv31 = du3tdx*M_loc[q](0,0) + du3tdy*M_loc[q](1,0) + du3tdz*M_loc[q](2,0);
            const value_type dFSv32 = du3tdx*M_loc[q](0,1) + du3tdy*M_loc[q](1,1) + du3tdz*M_loc[q](2,1);
            const value_type dFSv33 = du3tdx*M_loc[q](0,2) + du3tdy*M_loc[q](1,2) + du3tdz*M_loc[q](2,2);


            //auto dE = alpha_f*sym(gradt(u)) + alpha_f*alpha_f*0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            const value_type dE11 = du1tdx + /*0.5**/(du1vdx*du1tdx + du2vdx*du2tdx + du3vdx*du3tdx);
            const value_type dE12 = 0.5*( du1tdy + du2tdx + du1tdx*du1vdy + du2tdx*du2vdy + du3tdx*du3vdy + du1vdx*du1tdy + du2vdx*du2tdy + du3vdx*du3tdy );
            const value_type dE13 = 0.5*( du1tdz + du3tdx + du1tdx*du1vdz + du2tdx*du2vdz + du3tdx*du3vdz + du1vdx*du1tdz + du2vdx*du2tdz + du3vdx*du3tdz );
            const value_type dE21 = 0.5*( du2tdx + du1tdy + du1tdy*du1vdx + du2tdy*du2vdx + du3tdy*du3vdx + du1vdy*du1tdx + du2vdy*du2tdx + du3vdy*du3tdx );
            const value_type dE22 = du2tdy + du1vdy*du1tdy + du2vdy*du2tdy + du3vdy*du3tdy;
            const value_type dE23 = 0.5*( du2tdz + du3tdy + du1tdy*du1vdz + du2tdy*du2vdz + du3tdy*du3vdz + du1vdy*du1tdz + du2vdy*du2tdz + du3vdy*du3tdz );
            const value_type dE31 = 0.5*( du3tdx + du1tdz + du1tdz*du1vdx + du2tdz*du2vdx + du3tdz*du3vdx + du1vdz*du1tdx + du2vdz*du2tdx + du3vdz*du3tdx );
            const value_type dE32 = 0.5*( du3tdy + du2tdz + du1tdz*du1vdy + du2tdz*du2vdy + du3tdz*du3vdy + du1vdz*du1tdy + du2vdz*du2tdy + du3vdz*du3tdy );
            const value_type dE33 = du3tdz + du1vdz*du1tdz + du2vdz*du2tdz + du3vdz*du3tdz;

            const value_type tracedE = dE11 + dE22 + dE33;
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#else
            const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);
#endif
            // auto dS = idv(M_P0Coefflame1)*trace(dE)*Id + 2*idv(M_P0Coefflame2)*dE;
            const value_type dS11 = lame1*tracedE + 2*lame2*dE11;
            const value_type dS12 = 2*lame2*dE12;
            const value_type dS13 = 2*lame2*dE13;
            const value_type dS21 = 2*lame2*dE21;
            const value_type dS22 = lame1*tracedE + 2*lame2*dE22;
            const value_type dS23 = 2*lame2*dE23;
            const value_type dS31 = 2*lame2*dE31;
            const value_type dS32 = 2*lame2*dE32;
            const value_type dS33 = lame1*tracedE + 2*lame2*dE33;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = (1.+du1vdx)*dS11 + du1vdy*dS21 + du1vdz*dS31;
            const value_type FvdS12 = (1.+du1vdx)*dS12 + du1vdy*dS22 + du1vdz*dS32;
            const value_type FvdS13 = (1.+du1vdx)*dS13 + du1vdy*dS23 + du1vdz*dS33;
            const value_type FvdS21 = du2vdx*dS11 + (1.+du2vdy)*dS21 + du2vdz*dS31;
            const value_type FvdS22 = du2vdx*dS12 + (1.+du2vdy)*dS22 + du2vdz*dS32;
            const value_type FvdS23 = du2vdx*dS13 + (1.+du2vdy)*dS23 + du2vdz*dS33;
            const value_type FvdS31 = du3vdx*dS11 + du3vdy*dS21 + (1.+du3vdz)*dS31;
            const value_type FvdS32 = du3vdx*dS12 + du3vdy*dS22 + (1.+du3vdz)*dS32;
            const value_type FvdS33 = du3vdx*dS13 + du3vdy*dS23 + (1.+du3vdz)*dS33;

            // dF*val(Sv) + val(Fv)*dS)*trans(grad(v))
            const value_type resTrial11 = dFSv11 + FvdS11;
            const value_type resTrial12 = dFSv12 + FvdS12;
            const value_type resTrial13 = dFSv13 + FvdS13;
            const value_type resTrial21 = dFSv21 + FvdS21;
            const value_type resTrial22 = dFSv22 + FvdS22;
            const value_type resTrial23 = dFSv23 + FvdS23;
            const value_type resTrial31 = dFSv31 + FvdS31;
            const value_type resTrial32 = dFSv32 + FvdS32;
            const value_type resTrial33 = dFSv32 + FvdS33;

            // trace( (dF*val(Sv) + val(Fv)*dS)*trans(grad(v)) )
            const value_type res11 = resTrial11*du1dx + resTrial12*du1dy + resTrial13*du1dz;
            const value_type res22 = resTrial21*du2dx + resTrial22*du2dy + resTrial23*du2dz;
            const value_type res33 = resTrial31*du3dx + resTrial32*du3dy + resTrial33*du3dz;
            return res11+res22+res33;
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evaliq( i,c1,c2,q,mpl::int_<gmc_type::nDim>() );
        }
        value_type
        evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<2> /*Dim*/ ) const
        {
            const value_type res11 = M_loc[q](0,0)*M_fecTest->grad( i, 0, 0, q ) + M_loc[q](0,1)*M_fecTest->grad( i, 0, 1, q );
            const value_type res22 = M_loc[q](1,0)*M_fecTest->grad( i, 1, 0, q ) + M_loc[q](1,1)*M_fecTest->grad( i, 1, 1, q );
            return res11+res22;
        }
        value_type
        evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<3> /*Dim*/ ) const
        {
            const value_type res11 = M_loc[q](0,0)*M_fecTest->grad( i, 0, 0, q ) + M_loc[q](0,1)*M_fecTest->grad( i, 0, 1, q ) + M_loc[q](0,2)*M_fecTest->grad( i, 0, 2, q );
            const value_type res22 = M_loc[q](1,0)*M_fecTest->grad( i, 1, 0, q ) + M_loc[q](1,1)*M_fecTest->grad( i, 1, 1, q ) + M_loc[q](1,2)*M_fecTest->grad( i, 1, 2, q );
            const value_type res33 = M_loc[q](2,0)*M_fecTest->grad( i, 2, 0, q ) + M_loc[q](2,1)*M_fecTest->grad( i, 2, 1, q ) + M_loc[q](2,2)*M_fecTest->grad( i, 2, 2, q );
            return res11+res22+res33;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_loc[q]( c1,c2 );
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_loc[q];
        }

    private:

        void update( mpl::int_<2> /**/, mpl::int_<0> /*TypeOfExpr*/ )
        {
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#endif
            for ( uint16_type q = 0; q < M_geot->nPoints(); ++q )
                {
                    const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);

                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);

                    const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2));
                    const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2));
                    const value_type traceE = du1vdx+subtraceE1 + du2vdy+subtraceE2 ;

                    const value_type E11 = du1vdx + subtraceE1;
                    const value_type E12 = 0.5*( du1vdy + du2vdx + du1vdx*du1vdy + du2vdx*du2vdy );
                    const value_type E21 = 0.5*( du2vdx + du1vdy + du1vdy*du1vdx + du2vdy*du2vdx );
                    const value_type E22 = du2vdy + subtraceE2;

                    const value_type S11 = lame1*traceE + 2*lame2*E11;
                    const value_type S12 = 2*lame2*E12;
                    const value_type S21 = 2*lame2*E21;
                    const value_type S22 = lame1*traceE + 2*lame2*E22;

                    M_loc[q](0,0) /*const value_type FS11*/ = (1+du1vdx)*S11 + du1vdy*S21;
                    M_loc[q](0,1) /*const value_type FS12*/ = (1+du1vdx)*S12 + du1vdy*S22;
                    M_loc[q](1,0) /*const value_type FS21*/ = du2vdx*S11 + (1+du2vdy)*S21;
                    M_loc[q](1,1) /*const value_type FS22*/ = du2vdx*S12 + (1+du2vdy)*S22;

                }
        }

        void update( mpl::int_<2> /**/, mpl::int_<1> /*TypeOfExpr*/ )
        {
            update( mpl::int_<2>(), mpl::int_<0>() );
        }

        void update( mpl::int_<2> /**/, mpl::int_<2> /*TypeOfExpr*/ )
        {
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#endif

            for ( uint16_type q = 0; q < M_geot->nPoints(); ++q )
                {
                    const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);

                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);

                    const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2));
                    const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2));
                    const value_type traceE = du1vdx+subtraceE1 + du2vdy+subtraceE2 ;

                    const value_type E11 = du1vdx + subtraceE1;
                    const value_type E12 = 0.5*( du1vdy + du2vdx + du1vdx*du1vdy + du2vdx*du2vdy );
                    const value_type E21 = 0.5*( du2vdx + du1vdy + du1vdy*du1vdx + du2vdy*du2vdx );
                    const value_type E22 = du2vdy + subtraceE2;


                    M_loc[q](0,0) /*const value_type S11*/ = lame1*traceE + 2*lame2*E11;
                    M_loc[q](0,0) /*const value_type S12*/ = 2*lame2*E12;
                    M_loc[q](0,0) /*const value_type S21*/ = 2*lame2*E21;
                    M_loc[q](0,0) /*const value_type S22*/ = lame1*traceE + 2*lame2*E22;
                }

        }

        void update( mpl::int_<3> /*Dimension*/, mpl::int_<0> /*TypeOfExpr*/ )
        {
            //std::cout << "M_expr.e().coeffLame1() " << M_expr.coeffLame1()(0) << std::endl;
            //std::cout << "M_expr.e().coeffLame2() " << M_expr.coeffLame2()(0) << std::endl;
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#endif
            for ( uint16_type q = 0; q < M_geot->nPoints(); ++q )
                {
                    const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);

                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
                    const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);

                    const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2)+std::pow(du3vdx,2));
                    const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2)+std::pow(du3vdy,2));
                    const value_type subtraceE3 = 0.5*(std::pow(du1vdz,2)+std::pow(du2vdz,2)+std::pow(du3vdz,2));
                    const value_type traceE = du1vdx+subtraceE1 + du2vdy+subtraceE2 + du3vdz+subtraceE3 ;

                    const value_type E11 = du1vdx + subtraceE1;
                    const value_type E12 = 0.5*( du1vdy + du2vdx + du1vdx*du1vdy + du2vdx*du2vdy + du3vdx*du3vdy );
                    const value_type E13 = 0.5*( du1vdz + du3vdx + du1vdx*du1vdz + du2vdx*du2vdz + du3vdx*du3vdz );
                    const value_type E21 = 0.5*( du2vdx + du1vdy + du1vdy*du1vdx + du2vdy*du2vdx + du3vdy*du3vdx );
                    const value_type E22 = du2vdy + subtraceE2;
                    const value_type E23 = 0.5*( du2vdz + du3vdy + du1vdy*du1vdz + du2vdy*du2vdz + du3vdy*du3vdz );
                    const value_type E31 = 0.5*( du3vdx + du1vdz + du1vdz*du1vdx + du2vdz*du2vdx + du3vdz*du3vdx );
                    const value_type E32 = 0.5*( du3vdy + du2vdz + du1vdz*du1vdy + du2vdz*du2vdy + du3vdz*du3vdy );
                    const value_type E33 = du3vdz + subtraceE3;

                    const value_type S11 = lame1*traceE + 2*lame2*E11;
                    const value_type S12 = 2*lame2*E12;
                    const value_type S13 = 2*lame2*E13;
                    const value_type S21 = 2*lame2*E21;
                    const value_type S22 = lame1*traceE + 2*lame2*E22;
                    const value_type S23 = 2*lame2*E23;
                    const value_type S31 = 2*lame2*E31;
                    const value_type S32 = 2*lame2*E32;
                    const value_type S33 = lame1*traceE + 2*lame2*E33;

                    M_loc[q](0,0) /*const value_type FS11*/ = (1+du1vdx)*S11 + du1vdy*S21 + du1vdz*S31;
                    M_loc[q](0,1) /*const value_type FS12*/ = (1+du1vdx)*S12 + du1vdy*S22 + du1vdz*S32;
                    M_loc[q](0,2) /*const value_type FS13*/ = (1+du1vdx)*S13 + du1vdy*S23 + du1vdz*S33;
                    M_loc[q](1,0) /*const value_type FS21*/ = du2vdx*S11 + (1+du2vdy)*S21 + du2vdz*S31;
                    M_loc[q](1,1) /*const value_type FS22*/ = du2vdx*S12 + (1+du2vdy)*S22 + du2vdz*S32;
                    M_loc[q](1,2) /*const value_type FS23*/ = du2vdx*S13 + (1+du2vdy)*S23 + du2vdz*S33;
                    M_loc[q](2,0) /*const value_type FS31*/ = du3vdx*S11 + du3vdy*S21 + (1+du3vdz)*S31;
                    M_loc[q](2,1) /*const value_type FS32*/ = du3vdx*S12 + du3vdy*S22 + (1+du3vdz)*S32;
                    M_loc[q](2,2) /*const value_type FS33*/ = du3vdx*S13 + du3vdy*S23 + (1+du3vdz)*S33;
                }
        }


        void update( mpl::int_<3> /*Dimension*/, mpl::int_<1> /*TypeOfExpr*/ )
        {
            update( mpl::int_<3>(), mpl::int_<0>() );
        }

        void update( mpl::int_<3> /*Dimension*/, mpl::int_<2> /*TypeOfExpr*/ )
        {
#if 0
            const value_type lame1 = M_expr.coeffLame1()(0), lame2=M_expr.coeffLame2()(0);
#endif
            for ( uint16_type q = 0; q < M_geot->nPoints(); ++q )
                {
                    const value_type lame1 = M_locLameCoeff1[q](0,0), lame2 = M_locLameCoeff2[q](0,0);
                    //std::cout << " lame1 " << lame1 << "M_locLameCoeff1[q](0,0) "<< M_locLameCoeff1[q](0,0)  << std::endl;
                    //std::cout << " lame2 " << lame2 << "M_locLameCoeff2[q](0,0) "<< M_locLameCoeff2[q](0,0)  << std::endl;

                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
                    const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);

                    const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2)+std::pow(du3vdx,2));
                    const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2)+std::pow(du3vdy,2));
                    const value_type subtraceE3 = 0.5*(std::pow(du1vdz,2)+std::pow(du2vdz,2)+std::pow(du3vdz,2));
                    const value_type traceE = du1vdx+subtraceE1 + du2vdy+subtraceE2 + du3vdz+subtraceE3 ;

                    const value_type E11 = du1vdx + subtraceE1;
                    const value_type E12 = 0.5*( du1vdy + du2vdx + du1vdx*du1vdy + du2vdx*du2vdy + du3vdx*du3vdy );
                    const value_type E13 = 0.5*( du1vdz + du3vdx + du1vdx*du1vdz + du2vdx*du2vdz + du3vdx*du3vdz );
                    const value_type E21 = 0.5*( du2vdx + du1vdy + du1vdy*du1vdx + du2vdy*du2vdx + du3vdy*du3vdx );
                    const value_type E22 = du2vdy + subtraceE2;
                    const value_type E23 = 0.5*( du2vdz + du3vdy + du1vdy*du1vdz + du2vdy*du2vdz + du3vdy*du3vdz );
                    const value_type E31 = 0.5*( du3vdx + du1vdz + du1vdz*du1vdx + du2vdz*du2vdx + du3vdz*du3vdx );
                    const value_type E32 = 0.5*( du3vdy + du2vdz + du1vdz*du1vdy + du2vdz*du2vdy + du3vdz*du3vdy );
                    const value_type E33 = du3vdz + subtraceE3;

                    M_loc[q](0,0) /*const value_type S11*/ = lame1*traceE + 2*lame2*E11;
                    M_loc[q](0,0) /*const value_type S12*/ = 2*lame2*E12;
                    M_loc[q](0,0) /*const value_type S13*/ = 2*lame2*E13;
                    M_loc[q](0,0) /*const value_type S21*/ = 2*lame2*E21;
                    M_loc[q](0,0) /*const value_type S22*/ = lame1*traceE + 2*lame2*E22;
                    M_loc[q](0,0) /*const value_type S23*/ = 2*lame2*E23;
                    M_loc[q](0,0) /*const value_type S31*/ = 2*lame2*E31;
                    M_loc[q](0,0) /*const value_type S32*/ = 2*lame2*E32;
                    M_loc[q](0,0) /*const value_type S33*/ = lame1*traceE + 2*lame2*E33;
                }

        }



    private:
        this_type const& M_expr;
        gmc_ptrtype M_geot;
        basis_fec_test_ptrtype M_fecTest;
        basis_fec_trial_ptrtype M_fecTrial;
        //const uint16_type M_np;
        pc_ptrtype M_pc;
        ctx_ptrtype M_ctx;
        pc_lamecoeff_ptrtype M_pcLameCoeff;
        ctx_lamecoeff_ptrtype M_ctxLameCoeff;

        array_type M_loc;
        //array_id_type M_locId;
        array_grad_type M_locGrad;
        array_scalar_type M_locLameCoeff1,M_locLameCoeff2;

        //locAssembly_LinearForm_type M_locAssemblyLF;

        //ret_type M_zero;

        //tensor_expr_type M_tensor_expr;
        //vector_type M_det;
    };

private:
    boost::reference_wrapper<const element_type> M_v;
    boost::reference_wrapper<const element_lamecoeff_type> M_coeffLame1,M_coeffLame2;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */

template<class ElementType, class ElementLameCoeffType>
inline
Expr< StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<0> > >
stressStVenantKirchhoff( ElementType const& v, ElementLameCoeffType const&  lambda, ElementLameCoeffType const& mu )
{
    typedef StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<0> > stressStVenantKirchhoff_t;
    return Expr< stressStVenantKirchhoff_t >(  stressStVenantKirchhoff_t( v, lambda, mu ) );
}

template<class ElementType, class ElementLameCoeffType>
inline
Expr< StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<1> > >
stressStVenantKirchhoffResidual( ElementType const& v, ElementLameCoeffType const& lambda, ElementLameCoeffType const& mu )
{
    typedef StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<1>> stressStVenantKirchhoff_t;
    return Expr< stressStVenantKirchhoff_t >(  stressStVenantKirchhoff_t( v, lambda, mu ) );
}

template<class ElementType, class ElementLameCoeffType>
inline
Expr< StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<1> > >
stressStVenantKirchhoffResidual( std::shared_ptr<ElementType> const& v, ElementLameCoeffType const& lambda, ElementLameCoeffType const& mu )
{
    typedef StressStVenantKirchhoff<ElementType,ElementLameCoeffType, mpl::int_<1> > stressStVenantKirchhoff_t;
    return Expr< stressStVenantKirchhoff_t >(  stressStVenantKirchhoff_t( *v, lambda, mu ) );
}


template<class ElementType, class ElementLameCoeffType>
inline
Expr< StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<2> > >
stressStVenantKirchhoffJacobian( ElementType const& v, ElementLameCoeffType const& lambda, ElementLameCoeffType const& mu )
{
    typedef StressStVenantKirchhoff<ElementType,ElementLameCoeffType,mpl::int_<2> > stressStVenantKirchhoff_t;
    return Expr< stressStVenantKirchhoff_t >(  stressStVenantKirchhoff_t( v, lambda, mu ) );
}


} // namespace FeelModels
} // namespace Feel
#endif /* __SOLIDMECSTVENANTKIRCHHOFF_H */
