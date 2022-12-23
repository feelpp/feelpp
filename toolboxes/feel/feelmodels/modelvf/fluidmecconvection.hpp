/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2012-04-26

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
#ifndef __FEELPP_FSI_VF_FLUIDMECCONVECTION_H
#define __FEELPP_FSI_VF_FLUIDMECCONVECTION_H 1

namespace Feel {
namespace FeelModels {

/**
 * \class FluidMecConvectionImpl
 * \brief convection term in jac/res
 *
 * @author Vincent Chabannes
 * @see
 */
template<typename ElementVelocityType,typename SpecificExprType>
class FluidMecConvectionImpl
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FluidMecConvectionImpl<ElementVelocityType,SpecificExprType> this_type;

    static const size_type context = mpl::if_<mpl::or_<boost::is_same<SpecificExprType,mpl::int_<0> >,
                                                       boost::is_same<SpecificExprType,mpl::int_<1> > >,
                                              mpl::int_<vm::JACOBIAN|vm::KB>,
                                              mpl::int_<vm::JACOBIAN|vm::KB|vm::GRAD> >::type::value;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;

    typedef ElementVelocityType element_type;

    //------------------------------------------------------------------------------//
    // velocity functionspace
    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename functionspace_type::reference_element_type* fe_ptrtype;
    typedef typename functionspace_type::reference_element_type fe_type;
    //------------------------------------------------------------------------------//

    typedef typename functionspace_type::geoelement_type geoelement_type;
    typedef typename functionspace_type::value_type value_type;
    typedef value_type evaluate_type;

    static const uint16_type rank = fe_type::rank;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type ordervelocity = functionspace_type::basis_type::nOrder;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = mpl::if_<mpl::or_< boost::is_same<SpecificExprType,mpl::int_<0> >,
                                                      boost::is_same<SpecificExprType,mpl::int_<1> > >,
                                            mpl::bool_<false>,
                                            typename mpl::if_<boost::is_same<Func,fe_type>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type >::type::value;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<mpl::or_< boost::is_same<SpecificExprType,mpl::int_<2> >,boost::is_same<SpecificExprType,mpl::int_<3> > >,
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

    FluidMecConvectionImpl( element_type const & v )
        :
        M_v( boost::cref(v) )
    {}
    FluidMecConvectionImpl( FluidMecConvectionImpl const & op ) = default;
    ~FluidMecConvectionImpl()
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
    uint16_type polynomialOrder() const { return (SpecificExprType::value == 0 || SpecificExprType::value == 1)? (2*ordervelocity-1) : (3*ordervelocity-1); }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    element_type const& velocity() const { return M_v; }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename element_type::value_type value_type;


        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        // fe velocity context
        typedef typename fe_type::PreCompute pc_type;
        typedef std::shared_ptr<pc_type> pc_ptrtype;
        typedef typename fe_type::template Context<context_velocity, fe_type, gm_type,geoelement_type/*,gmc_type::context*/> ctx_type;
        typedef std::shared_ptr<ctx_type> ctx_ptrtype;


        // fe context for test and trial function
        typedef Basis_i_t map_basis_fec_test_type;
        typedef Basis_j_t map_basis_fec_trial_type;
        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_test_type,Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type basis_fec_test_key_type;

        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_trial_type,vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type basis_fec_trial_key_type;

        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type  basis_fec_test_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type* basis_fec_test_ptrtype;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type  basis_fec_trial_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type* basis_fec_trial_ptrtype;


        // output and useful container
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape_scalar;

        typedef Shape<gmc_type::nDim, Vectorial, false, false> shape_vectorial;
        //typedef Eigen::Matrix<value_type,shape_vectorial::M,shape_vectorial::N> loc_vectorial_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_vectorial::M,shape_vectorial::N>> loc_vectorial_type;
        typedef boost::multi_array<loc_vectorial_type,1> array_vectorial_type;

        typedef Shape<gmc_type::nDim, Tensor2, false, false> shape_tensor2;
        //typedef Eigen::Matrix<value_type,shape_tensor2::M,shape_tensor2::N> loc_tensor2_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_tensor2::M,shape_tensor2::N>> loc_tensor2_type;
        typedef boost::multi_array<loc_tensor2_type,1> array_tensor2_type;

        typedef loc_vectorial_type loc_id_type;
        typedef array_vectorial_type array_id_type;
        typedef loc_tensor2_type loc_grad_type;
        typedef array_tensor2_type array_grad_type;

        typedef Eigen::Matrix<value_type,shape_vectorial::M,shape_vectorial::N> loc_matrix_vectorial_type;
        typedef boost::multi_array<loc_matrix_vectorial_type,1> array_matrix_vectorial_type;

        using ret_type = Eigen::Map<loc_matrix_vectorial_type const>;

        typedef typename mpl::if_<mpl::or_< boost::is_same<SpecificExprType,mpl::int_<0> >,
                                            boost::is_same<SpecificExprType,mpl::int_<1> > >,
                                  shape_vectorial,
                                  shape_scalar >::type shape;

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
            M_pc( new pc_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_loc( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locId( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_zeroLocGrad(),
            M_zeroLocId()
            //M_zero( ret_type::Zero() )
            {
                M_zeroLocGrad.setZero();
                M_zeroLocId.setZero();
            }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_fecTest( fusion::at_key<basis_fec_test_key_type>( fev ).get() ),
            M_pc( new pc_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_loc( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locId( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_zeroLocGrad(),
            M_zeroLocId()
            //M_zero( ret_type::Zero() )
            {
                M_zeroLocGrad.setZero();
                M_zeroLocId.setZero();
            }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pc( new pc_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctx( new ctx_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ),
            M_loc( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locId( expr.velocity().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGrad( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_zeroLocGrad(),
            M_zeroLocId()
            //M_zero( ret_type::Zero() )
            {
                M_zeroLocGrad.setZero();
                M_zeroLocId.setZero();
            }

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& feu )
        {
            std::fill( M_locId.data(), M_locId.data()+M_locId.num_elements(), M_zeroLocId/*loc_id_type::Zero()*/ );
            std::fill( M_locGrad.data(), M_locGrad.data()+M_locGrad.num_elements(), M_zeroLocGrad/*loc_grad_type::Zero()*/ );

            const uint16_type nQuadPts = M_geot->nPoints();

            //M_ctx->update( fusion::at_key<key_type>( geom ),  (pc_ptrtype const&) M_pc );
#if 0
            M_expr.velocity().id( *M_ctx, M_locId );
            M_expr.velocity().grad( *M_ctx, M_locGrad );
#else
            const size_type elt_id = fusion::at_key<key_type>( geom )->id();
            for ( int l = 0; l < functionspace_type::basis_type::nDof; ++l )
                {
                    const int ncdof = functionspace_type::is_product?nComponents1:1;
                    for ( typename array_vectorial_type::index c1 = 0; c1 < ncdof; ++c1 )
                        {
                            typename array_vectorial_type::index ldof = functionspace_type::basis_type::nDof*c1+l;
                            const size_type gdof = M_expr.velocity().functionSpace()->dof()->localToGlobal( elt_id, l, c1 ).index();
                            const value_type v_ = M_expr.velocity().globalValue( gdof );

                            for ( uint16_type q = 0; q < nQuadPts; ++q )
                                {
                                    for ( typename array_vectorial_type::index i = 0; i < nComponents1; ++i )
                                        {
                                            M_locId[q]( i,0 ) += v_*M_fecTrial/*M_ctx*/->id( ldof, i, 0, q );

                                            for ( uint16_type j = 0; j < gmc_type::nDim/*nRealDim*/; ++j )
                                                {
                                                    M_locGrad[q]( i,j ) += v_*M_fecTrial/*M_ctx*/->grad( ldof, i, j, q );
                                                }
                                        }
                                }
                        }
                }
#endif
            //update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom )
        {
            std::fill( M_locId.data(), M_locId.data()+M_locId.num_elements(), M_zeroLocId/*loc_id_type::Zero()*/ );
            std::fill( M_locGrad.data(), M_locGrad.data()+M_locGrad.num_elements(), M_zeroLocGrad/*loc_grad_type::Zero()*/ );

            const uint16_type nQuadPts = M_geot->nPoints();

            M_ctx->update( fusion::at_key<key_type>( geom ),  (pc_ptrtype const&) M_pc );

#if 0
            M_expr.velocity().id( *M_ctx, M_locId );
            M_expr.velocity().grad( *M_ctx, M_locGrad );
#else
            const size_type elt_id = fusion::at_key<key_type>( geom )->id();
            for ( int l = 0; l < functionspace_type::basis_type::nDof; ++l )
                {
                    const int ncdof = functionspace_type::is_product?nComponents1:1;
                    for ( typename array_vectorial_type::index c1 = 0; c1 < ncdof; ++c1 )
                        {
                            typename array_vectorial_type::index ldof = functionspace_type::basis_type::nDof*c1+l;
                            const size_type gdof = M_expr.velocity().functionSpace()->dof()->localToGlobal( elt_id, l, c1 ).index();
                            const value_type v_ = M_expr.velocity().globalValue( gdof );

                            for ( uint16_type q = 0; q < nQuadPts; ++q )
                                {
                                    for ( typename array_vectorial_type::index i = 0; i < nComponents1; ++i )
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

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalijq( i, j, c1, c2, q, mpl::int_<gmc_type::nDim>(), mpl::int_<SpecificExprType::value>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> /*Dim*/, mpl::int_<2> /*TypeOfExpr*/ ) const
        {
            //auto const convecTerm = (trans(val(gradv(u)*idv(*M_P0Rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(*M_P0Rho)) ) )*id(v);
            const value_type u1t = M_fecTrial->id( j, 0, 0, q ), u2t = M_fecTrial->id( j, 1, 0, q );
            const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0);
            const value_type u1 = M_fecTest->id( i, 0, 0, q ), u2 = M_fecTest->id( i, 1, 0, q );

            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q );

            const value_type part1 = du1vdx*u1t + du1vdy*u2t + du1tdx*u1v + du1tdy*u2v;
            const value_type part2 = du2vdx*u1t + du2vdy*u2t + du2tdx*u1v + du2tdy*u2v;

            const value_type res = part1*u1 + part2*u2;

            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<3> /*Dim*/, mpl::int_<2> /*TypeOfExpr*/ ) const
        {
            const value_type u1t = M_fecTrial->id( j, 0, 0, q ), u2t = M_fecTrial->id( j, 1, 0, q ), u3t = M_fecTrial->id( j, 2, 0, q );
            const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0), u3v = M_locId[q](2,0);
            const value_type u1 = M_fecTest->id( i, 0, 0, q ), u2 = M_fecTest->id( i, 1, 0, q ), u3 = M_fecTest->id( i, 2, 0, q );

            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
            const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q ), du1tdz=M_fecTrial->grad( j, 0, 2, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q ), du2tdz=M_fecTrial->grad( j, 1, 2, q );
            const value_type du3tdx = M_fecTrial->grad( j, 2, 0, q ), du3tdy = M_fecTrial->grad( j, 2, 1, q ), du3tdz=M_fecTrial->grad( j, 2, 2, q );

            const value_type part1 = du1vdx*u1t + du1vdy*u2t + du1vdz*u3t + du1tdx*u1v + du1tdy*u2v + du1tdz*u3v;
            const value_type part2 = du2vdx*u1t + du2vdy*u2t + du2vdz*u3t + du2tdx*u1v + du2tdy*u2v + du2tdz*u3v;
            const value_type part3 = du3vdx*u1t + du3vdy*u2t + du3vdz*u3t + du3tdx*u1v + du3tdy*u2v + du3tdz*u3v;

            const value_type res = part1*u1 + part2*u2 + part3*u3;

            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> /*Dim*/, mpl::int_<3> /*TypeOfExpr*/ ) const
        {
            const value_type u1t = M_fecTrial->id( j, 0, 0, q ), u2t = M_fecTrial->id( j, 1, 0, q );
            const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0);
            const value_type u1 = M_fecTest->id( i, 0, 0, q ), u2 = M_fecTest->id( i, 1, 0, q );

            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q );

            const value_type part1 = du1vdx*u1t + du1vdy*u2t + du1tdx*u1v + du1tdy*u2v;
            const value_type part2 = du2vdx*u1t + du2vdy*u2t + du2tdx*u1v + du2tdy*u2v;

            const value_type divut = du1tdx + du2tdy;
            const value_type divuv = du1vdx + du2vdy;
            const value_type stab1 = 0.5*(divut*u1v + divuv*u1t);
            const value_type stab2 = 0.5*(divut*u2v + divuv*u2t);

            const value_type res = (part1+stab1)*u1 + (part2+stab2)*u2;

            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<3> /*Dim*/, mpl::int_<3> /*TypeOfExpr*/ ) const
        {
            const value_type u1t = M_fecTrial->id( j, 0, 0, q ), u2t = M_fecTrial->id( j, 1, 0, q ), u3t = M_fecTrial->id( j, 2, 0, q );
            const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0), u3v = M_locId[q](2,0);
            const value_type u1 = M_fecTest->id( i, 0, 0, q ), u2 = M_fecTest->id( i, 1, 0, q ), u3 = M_fecTest->id( i, 2, 0, q );

            const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
            const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
            const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);

            const value_type du1tdx = M_fecTrial->grad( j, 0, 0, q ), du1tdy = M_fecTrial->grad( j, 0, 1, q ), du1tdz=M_fecTrial->grad( j, 0, 2, q );
            const value_type du2tdx = M_fecTrial->grad( j, 1, 0, q ), du2tdy = M_fecTrial->grad( j, 1, 1, q ), du2tdz=M_fecTrial->grad( j, 1, 2, q );
            const value_type du3tdx = M_fecTrial->grad( j, 2, 0, q ), du3tdy = M_fecTrial->grad( j, 2, 1, q ), du3tdz=M_fecTrial->grad( j, 2, 2, q );

            const value_type divut = du1tdx + du2tdy + du3tdz;
            const value_type divuv = du1vdx + du2vdy + du3vdz;

            const value_type part1 = du1vdx*u1t + du1vdy*u2t + du1vdz*u3t + du1tdx*u1v + du1tdy*u2v + du1tdz*u3v;
            const value_type part2 = du2vdx*u1t + du2vdy*u2t + du2vdz*u3t + du2tdx*u1v + du2tdy*u2v + du2tdz*u3v;
            const value_type part3 = du3vdx*u1t + du3vdy*u2t + du3vdz*u3t + du3tdx*u1v + du3tdy*u2v + du3tdz*u3v;

            const value_type stab1 = 0.5*(divut*u1v + divuv*u1t);
            const value_type stab2 = 0.5*(divut*u2v + divuv*u2t);
            const value_type stab3 = 0.5*(divut*u3v + divuv*u3t);

            const value_type res = (part1+stab1)*u1 + (part2+stab2)*u2 + (part3+stab3)*u3;

            // trans(divt(u)*val(0.5*idv(*M_P0Rho)*idv(u))+val(0.5*idv(*M_P0Rho)*divv(u))*idt(u))*id(v),

            return res;
        }



        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return ret_type(M_loc[q].data());
        }

        value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_loc[q](c1,0);
        }
        
        ret_type
        evalq( uint16_type q ) const
        {
            return ret_type(M_loc[q].data());
        }

    private:

        /**
         *
         */
        void update( mpl::int_<2> /*Dim*/ , mpl::int_<0> /*TypeOfExpr*/ )
        {
            //std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_vectorial_type::Zero() );
            const uint16_type nQuadPts = M_geot->nPoints();
            for ( uint16_type q = 0; q < nQuadPts; ++q )
                {
                    const value_type locId0 = M_locId[q](0,0);
                    const value_type locId1 = M_locId[q](1,0);
                    M_loc[q](0,0) = M_locGrad[q](0,0)*locId0 + M_locGrad[q](0,1)*locId1;
                    M_loc[q](1,0) = M_locGrad[q](1,0)*locId0 + M_locGrad[q](1,1)*locId1;
                    //M_loc[q](0,0) = M_locGrad[q](0,0)*M_locId[q](0,0) + M_locGrad[q](0,1)*M_locId[q](1,0);
                    //M_loc[q](1,0) = M_locGrad[q](1,0)*M_locId[q](0,0) + M_locGrad[q](1,1)*M_locId[q](1,0);
                }
        }
        void update( mpl::int_<3> /*Dim*/, mpl::int_<0> /*TypeOfExpr*/)
        {
            //std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_vectorial_type::Zero() );
            const uint16_type nQuadPts = M_geot->nPoints();
            for ( uint16_type q = 0; q < nQuadPts; ++q )
                {
                    M_loc[q](0,0) = M_locGrad[q](0,0)*M_locId[q](0,0) + M_locGrad[q](0,1)*M_locId[q](1,0) + M_locGrad[q](0,2)*M_locId[q](2,0);
                    M_loc[q](1,0) = M_locGrad[q](1,0)*M_locId[q](0,0) + M_locGrad[q](1,1)*M_locId[q](1,0) + M_locGrad[q](1,2)*M_locId[q](2,0);
                    M_loc[q](2,0) = M_locGrad[q](2,0)*M_locId[q](0,0) + M_locGrad[q](2,1)*M_locId[q](1,0) + M_locGrad[q](2,2)*M_locId[q](2,0);
                }
        }
        /**
         *
         */
        void update( mpl::int_<2> /*Dim*/ , mpl::int_<1> /*TypeOfExpr*/ )
        {
            const uint16_type nQuadPts = M_geot->nPoints();
            for ( uint16_type q = 0; q < nQuadPts; ++q )
                {
                    const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0);
                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1);
                    const value_type divuv = du1vdx + du2vdy;

                    M_loc[q](0,0) = du1vdx*u1v + du1vdy*u2v + 0.5*divuv*u1v;
                    M_loc[q](1,0) = du2vdx*u1v + du2vdy*u2v + 0.5*divuv*u2v;
                }
        }
        void update( mpl::int_<3> /*Dim*/, mpl::int_<1> /*TypeOfExpr*/)
        {
            const uint16_type nQuadPts = M_geot->nPoints();
            for ( uint16_type q = 0; q < nQuadPts; ++q )
                {
                    const value_type u1v = M_locId[q](0,0), u2v = M_locId[q](1,0), u3v = M_locId[q](2,0);
                    const value_type du1vdx = M_locGrad[q](0,0), du1vdy = M_locGrad[q](0,1), du1vdz = M_locGrad[q](0,2);
                    const value_type du2vdx = M_locGrad[q](1,0), du2vdy = M_locGrad[q](1,1), du2vdz = M_locGrad[q](1,2);
                    const value_type du3vdx = M_locGrad[q](2,0), du3vdy = M_locGrad[q](2,1), du3vdz = M_locGrad[q](2,2);
                    const value_type divuv = du1vdx + du2vdy + du3vdz;

                    M_loc[q](0,0) = du1vdx*u1v + du1vdy*u2v + du1vdz*u3v + 0.5*divuv*u1v;
                    M_loc[q](1,0) = du2vdx*u1v + du2vdy*u2v + du2vdz*u3v + 0.5*divuv*u2v;
                    M_loc[q](2,0) = du3vdx*u1v + du3vdy*u2v + du3vdz*u3v + 0.5*divuv*u3v;
                }
        }

    private:
        this_type const& M_expr;

        gmc_ptrtype M_geot;

        basis_fec_test_ptrtype M_fecTest;
        basis_fec_trial_ptrtype M_fecTrial;

        pc_ptrtype M_pc;
        ctx_ptrtype M_ctx;

        array_matrix_vectorial_type M_loc;

        array_id_type M_locId;
        array_grad_type M_locGrad;

        //ret_type M_zero;
        loc_grad_type M_zeroLocGrad;
        loc_id_type M_zeroLocId;
    };

private:
    boost::reference_wrapper<const element_type> M_v;

};
/// \endcond

/**
 * \brief det of the expression tensor
 */
template<class ElementVelocityType>
inline
auto
fluidMecConvection( ElementVelocityType const& v )
{
    typedef FluidMecConvectionImpl<unwrap_ptr_t<ElementVelocityType>,mpl::int_<0> > convection_t;
    return Expr< convection_t >(  convection_t( unwrap_ptr(v) ) );
}

template<class ElementVelocityType>
inline
auto
fluidMecConvectionWithEnergyStab( ElementVelocityType const& v )
{
    typedef FluidMecConvectionImpl<unwrap_ptr_t<ElementVelocityType>,mpl::int_<1> > convection_t;
    return Expr< convection_t >(  convection_t( unwrap_ptr(v) ) );
}

template<class ElementVelocityType>
inline
auto
fluidMecConvectionJacobian( ElementVelocityType const& v )
{
    typedef FluidMecConvectionImpl<unwrap_ptr_t<ElementVelocityType>,mpl::int_<2> > convection_t;
    return Expr< convection_t >(  convection_t( unwrap_ptr(v) ) );
}

template<class ElementVelocityType>
inline
auto
fluidMecConvectionJacobianWithEnergyStab( ElementVelocityType const& v )
{
    typedef FluidMecConvectionImpl<unwrap_ptr_t<ElementVelocityType>,mpl::int_<3> > convection_t;
    return Expr< convection_t >(  convection_t( unwrap_ptr(v) ) );
}


} // namespace FeelModels
} // namespace Feel
#endif /* __FLUIDMECCONVECTION_H */
