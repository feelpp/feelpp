#ifndef _CAUCHY_GREEN_TENSOR_EXPR_HPP
#define _CAUCHY_GREEN_TENSOR_EXPR_HPP 1

namespace Feel {
namespace FeelModels {

template<typename ElementBackwardCharacteristicsType, typename ElementNormalType, int IMOrder>
class LeftCauchyGreenTensorExpr
{
public:
    typedef LeftCauchyGreenTensorExpr<ElementBackwardCharacteristicsType, ElementNormalType, IMOrder> this_type;

    typedef ElementBackwardCharacteristicsType element_backwardcharacteristics_type;
    typedef ElementNormalType element_normal_type;

    static const size_type context_backwardcharacteristics = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_normal = vm::JACOBIAN;
    static const size_type context = context_backwardcharacteristics;

    //------------------------------------------------------------------------------//
    // backward characteristics functionspace
    typedef typename element_backwardcharacteristics_type::functionspace_type functionspace_backwardcharacteristics_type;
    typedef typename functionspace_backwardcharacteristics_type::reference_element_type fe_backwardcharacteristics_type;
    //------------------------------------------------------------------------------//
    // normal functionspace
    typedef typename element_normal_type::functionspace_type functionspace_normal_type;
    typedef typename functionspace_normal_type::reference_element_type fe_normal_type;
    //------------------------------------------------------------------------------//
    // expression desc
    typedef typename functionspace_normal_type::geoelement_type geoelement_type;
    typedef typename functionspace_normal_type::value_type value_type;
    typedef value_type evaluate_type;
    static const uint16_type nDim = fe_normal_type::nDim;
    static const uint16_type nRealDim = fe_normal_type::nRealDim;
    static const uint16_type rank = fe_normal_type::rank;
    static const uint16_type nComponents1 = fe_normal_type::nComponents1;
    static const uint16_type nComponents2 = fe_normal_type::nComponents2;
    static const bool is_terminal = false;

    static const uint16_type orderbackwardcharacteristics = functionspace_backwardcharacteristics_type::basis_type::nOrder;
    static const uint16_type ordernormal = functionspace_normal_type::basis_type::nOrder;

    // imorder
    static const uint16_type imorderAuto = 2*(orderbackwardcharacteristics + ordernormal);
    static const uint16_type imorder = (IMOrder >= 0) ? IMOrder : imorderAuto;

    static const bool imIsPoly = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    //--------------------------------------------------------------------//
    // Constructors
    LeftCauchyGreenTensorExpr( element_backwardcharacteristics_type const& Y, element_normal_type const& N )
        : 
            M_backwardCharacteristics( boost::cref(Y) ),
            M_normal( boost::cref(N) )
    {}

    //! polynomial order
    uint16_type polynomialOrder() const { return imorder; }

    //! expression is polynomial?
    bool isPolynomial() const { return false; }

    //--------------------------------------------------------------------//
    // Accessors
    element_backwardcharacteristics_type const& backwardCharacteristics() const { return M_backwardCharacteristics; }
    element_normal_type const& normal() const { return M_normal; }

    //--------------------------------------------------------------------//
    // Expr tensor
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expr_type;
        typedef typename element_backwardcharacteristics_type::value_type value_type;

        struct is_zero { static const bool value = false; };

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        typedef typename expr_type::fe_normal_type fe_normal_type;
        typedef typename expr_type::fe_backwardcharacteristics_type fe_backwardcharacteristics_type;
        typedef typename expr_type::geoelement_type geoelement_type;

        // fe backwardcharacteristics context
        typedef typename fe_backwardcharacteristics_type::PreCompute pc_backwardcharacteristics_type;
        typedef std::shared_ptr<pc_backwardcharacteristics_type> pc_backwardcharacteristics_ptrtype;
        typedef typename fe_backwardcharacteristics_type::template Context<expr_type::context_backwardcharacteristics, fe_backwardcharacteristics_type, gm_type,geoelement_type> ctx_backwardcharacteristics_type;
        typedef std::shared_ptr<ctx_backwardcharacteristics_type> ctx_backwardcharacteristics_ptrtype;
        // fe normal context
        typedef typename fe_normal_type::PreCompute pc_normal_type;
        typedef std::shared_ptr<pc_normal_type> pc_normal_ptrtype;
        typedef typename fe_normal_type::template Context<expr_type::context_normal, fe_normal_type, gm_type,geoelement_type> ctx_normal_type;
        typedef std::shared_ptr<ctx_normal_type> ctx_normal_ptrtype;

        // Shapes types
        // res shape
        typedef Shape<expr_type::nRealDim, Tensor2, false, false> shape_type;
        typedef shape_type shape;
        typedef Eigen::Matrix<value_type,shape_type::M,shape_type::N> matrix_shape_type;
        typedef boost::multi_array<matrix_shape_type,1> array_shape_type;
        typedef Eigen::Matrix<value_type,shape_type::M,1> vector_shape_type;
        // scalar
        typedef Shape<shape_type::nDim, Scalar, false, false> shape_scalar;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_scalar::M,shape_scalar::N>> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type,1> array_scalar_type;
        // vectorial
        typedef Shape<shape_type::nDim, Vectorial, false, false> shape_vectorial;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_vectorial::M,shape_vectorial::N>> loc_vectorial_type;
        typedef boost::multi_array<loc_vectorial_type,1> array_vectorial_type;
        // vectorial transpose
        typedef Shape<shape_type::nDim, Vectorial, true, false> shape_vectorial_transpose;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_vectorial_transpose::M,shape_vectorial_transpose::N>> loc_vectorial_transpose_type;
        typedef boost::multi_array<loc_vectorial_transpose_type,1> array_vectorial_transpose_type;
        // tensor2
        typedef Shape<shape_type::nDim, Tensor2, false, false> shape_tensor2;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_tensor2::M,shape_tensor2::N>> loc_tensor2_type;
        typedef boost::multi_array<loc_tensor2_type,1> array_tensor2_type;

        tensor( expr_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
                M_expr( expr ),
                M_geot( fusion::at_key<key_type>( geom ) ),
                M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
                M_pcBackwardCharacteristics( new pc_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxBackwardCharacteristics( new ctx_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_backwardcharacteristics_ptrtype const&)M_pcBackwardCharacteristics ) ),
                M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locGradBackwardCharacteristics( expr.backwardCharacteristics().gradExtents(*fusion::at_key<key_type>( geom )) ),
                M_zeroLocScalar(),
                M_zeroLocVectorial(),
                M_zeroLocVectorialTranspose( ),
                M_zeroLocTensor2()
        {
            M_zeroLocScalar.setZero();
            M_zeroLocVectorial.setZero();
            M_zeroLocVectorialTranspose.setZero();
            M_zeroLocTensor2.setZero();
        }

        tensor( expr_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
                M_expr( expr ),
                M_geot( fusion::at_key<key_type>( geom ) ),
                M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
                M_pcBackwardCharacteristics( new pc_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxBackwardCharacteristics( new ctx_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_backwardcharacteristics_ptrtype const&)M_pcBackwardCharacteristics ) ),
                M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locGradBackwardCharacteristics( expr.backwardCharacteristics().gradExtents(*fusion::at_key<key_type>( geom )) ),
                M_zeroLocScalar(),
                M_zeroLocVectorial(),
                M_zeroLocVectorialTranspose(),
                M_zeroLocTensor2()
        {
            M_zeroLocScalar.setZero();
            M_zeroLocVectorial.setZero();
            M_zeroLocVectorialTranspose.setZero();
            M_zeroLocTensor2.setZero();
        }

        tensor( expr_type const& expr, Geo_t const& geom )
            :
                M_expr( expr ),
                M_geot( fusion::at_key<key_type>( geom ) ),
                M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
                M_pcBackwardCharacteristics( new pc_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
                M_ctxBackwardCharacteristics( new ctx_backwardcharacteristics_type( expr.backwardCharacteristics().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_backwardcharacteristics_ptrtype const&)M_pcBackwardCharacteristics ) ),
                M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
                M_locGradBackwardCharacteristics( expr.backwardCharacteristics().gradExtents(*fusion::at_key<key_type>( geom )) ),
                M_zeroLocScalar(),
                M_zeroLocVectorial(),
                M_zeroLocVectorialTranspose(),
                M_zeroLocTensor2()
        {
            M_zeroLocScalar.setZero();
            M_zeroLocVectorial.setZero();
            M_zeroLocVectorialTranspose.setZero();
            M_zeroLocTensor2.setZero();
        }

        template<typename IM>
        void init( IM const& im )
        {}

        expr_type const& expr() const { return M_expr; }

        void setGmc( Geo_t const& geom ) { M_geot = fusion::at_key<key_type>( geom ); }
        gmc_ptrtype const& gmc() const { return M_geot; }

        array_shape_type & locRes() { return M_locRes; }
        array_shape_type const& locRes() const { return M_locRes; }
        matrix_shape_type & locRes( uint16_type q ) { return M_locRes[q]; }
        matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }

        array_tensor2_type const& locGradBackwardCharacteristics() const { return M_locGradBackwardCharacteristics; }
        loc_tensor2_type const& locGradBackwardCharacteristics( uint16_type q ) const { return M_locGradBackwardCharacteristics[q]; }
        array_vectorial_type const& locNormal() const { return M_locNormal; }
        loc_vectorial_type const& locNormal( uint16_type q ) const { return M_locNormal[q]; }

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_ctxBackwardCharacteristics->update( fusion::at_key<key_type>( geom ),  (pc_backwardcharacteristics_ptrtype const&) M_pcBackwardCharacteristics );
            std::fill( M_locGradBackwardCharacteristics.data(), M_locGradBackwardCharacteristics.data()+M_locGradBackwardCharacteristics.num_elements(), this->M_zeroLocTensor2 );
            this->expr().backwardCharacteristics().grad( *M_ctxBackwardCharacteristics, M_locGradBackwardCharacteristics );

            M_ctxNormal->update( fusion::at_key<key_type>( geom ),  (pc_normal_ptrtype const&) M_pcNormal );
            std::fill( M_locNormal.data(), M_locNormal.data()+M_locNormal.num_elements(), this->M_zeroLocVectorial );
            this->expr().normal().id( *M_ctxNormal, M_locNormal );

            computeLeftCauchyGreenTensor( mpl::int_<shape::N>() );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            {
                M_pcBackwardCharacteristics->update( this->gmc()->pc()->nodes() );
                M_pcNormal->update( this->gmc()->pc()->nodes() );
            }
        this->update( geom );
        }

        matrix_shape_type const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return this->evalq( q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1,c2,q );
        }
        matrix_shape_type const&
        evaliq( uint16_type i, uint16_type q ) const
        {
            return this->evalq( q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_locRes[q](c1,c2);
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_locRes[q];
        }

    private:
        void computeLeftCauchyGreenTensor( mpl::int_<1> )
        {
        }
        void computeLeftCauchyGreenTensor( mpl::int_<2> )
        {
            /*
             * A = K - (KN)*(KN)^T / ( N^T*(KN) )
             * with K = (gradY)^-1 * (gradY)^-T = G^-1*G^-T = (G^T*G)^-1
             */
            for( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                vector_shape_type N = Eigen::Map<const vector_shape_type>( this->locNormal(q).data(), shape_type::M, 1 );
                matrix_shape_type gradY = Eigen::Map<const matrix_shape_type>( this->locGradBackwardCharacteristics(q).data(), shape_type::M, shape_type::N );

                matrix_shape_type K = ( gradY.transpose()*gradY ).inverse();

                this->locRes(q) = K - ( K * N * N.transpose() * K.transpose() ) / ( N.transpose() * K * N )(0);

            }
        }
        void computeLeftCauchyGreenTensor( mpl::int_<3> )
        {
            /*
             * A = K - (KN)*(KN)^T / ( N^T*(KN) )
             * with K = (gradY)^-1 * (gradY)^-T = G^-1*G^-T = (G^T*G)^-1
             */
            for( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                vector_shape_type N = Eigen::Map<const vector_shape_type>( this->locNormal(q).data(), shape_type::M, 1 );
                matrix_shape_type gradY = Eigen::Map<const matrix_shape_type>( this->locGradBackwardCharacteristics(q).data(), shape_type::M, shape_type::N );

                matrix_shape_type K = ( gradY.transpose()*gradY ).inverse();

                this->locRes(q) = K - ( K * N * N.transpose() * K.transpose() ) / ( N.transpose() * K * N )(0);
            }
        }

    private:
        expr_type const& M_expr;
        gmc_ptrtype M_geot;

        pc_normal_ptrtype M_pcNormal;
        ctx_normal_ptrtype M_ctxNormal;
        pc_backwardcharacteristics_ptrtype M_pcBackwardCharacteristics;
        ctx_backwardcharacteristics_ptrtype M_ctxBackwardCharacteristics;

        array_shape_type M_locRes;
        array_vectorial_type M_locNormal;
        array_tensor2_type M_locGradBackwardCharacteristics;

    protected:
        loc_scalar_type M_zeroLocScalar;
        loc_vectorial_type M_zeroLocVectorial;
        loc_vectorial_transpose_type M_zeroLocVectorialTranspose;
        loc_tensor2_type M_zeroLocTensor2;
    };


private:
    boost::reference_wrapper<const element_backwardcharacteristics_type> M_backwardCharacteristics;
    boost::reference_wrapper<const element_normal_type> M_normal;
};

template<int IMOrder = -1, typename ElementBackwardCharacteristicsType, typename ElementNormalType>
inline Expr< LeftCauchyGreenTensorExpr<ElementBackwardCharacteristicsType, ElementNormalType, IMOrder> >
leftCauchyGreenTensorExpr( ElementBackwardCharacteristicsType const& Y, ElementNormalType const& N )
{
    typedef LeftCauchyGreenTensorExpr< ElementBackwardCharacteristicsType, ElementNormalType, IMOrder > expr_leftcauchygreentensor_type;
    return Expr<expr_leftcauchygreentensor_type>( expr_leftcauchygreentensor_type( Y, N ) );
}



}
}
#endif
