#ifndef _LEVELSET_DELTA_EXPR_HPP
#define _LEVELSET_DELTA_EXPR_HPP 1

namespace Feel {

BOOST_PARAMETER_NAME( thickness )
BOOST_PARAMETER_NAME( use_adaptive_thickness )
BOOST_PARAMETER_NAME( impl )
BOOST_PARAMETER_NAME( imorder )

namespace FeelModels {

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorLevelsetDeltaBase
{
public:
    typedef ExprType expr_type;
    typedef typename expr_type::value_type value_type;

    typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
            mpl::identity<vf::detail::gmc<0> >,
            mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gmc_type::gm_type gm_type;

    typedef typename expr_type::fe_phi_type fe_phi_type;
    typedef typename expr_type::geoelement_type geoelement_type;

    // fe phi context
    typedef typename fe_phi_type::PreCompute pc_phi_type;
    typedef std::shared_ptr<pc_phi_type> pc_phi_ptrtype;
    typedef typename fe_phi_type::template Context<expr_type::context_phi, fe_phi_type, gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_phi_type;
    typedef std::shared_ptr<ctx_phi_type> ctx_phi_ptrtype;

    typedef typename matrix_node<value_type>::type matrix_node_type;

    static const bool isSameGeo = boost::is_same<typename gmc_type::element_type, geoelement_type>::value;

    // Shapes types
    // res shape
    typedef Shape<expr_type::nRealDim, Scalar, false, false> shape_type;
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


    tensorLevelsetDeltaBase( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcPhi( this->createPcPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_ctxPhi( this->createCtxPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_locRes( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPhi( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_isSameMesh( fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo(expr.phi().functionSpace()->mesh()) && isSameGeo ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocVectorialTranspose(),
            M_zeroLocTensor2()
    {
        if( !M_isSameMesh )
            expr.phi().functionSpace()->mesh()->tool_localization()->updateForUse();

        M_zeroLocScalar.setZero();
        M_zeroLocVectorial.setZero();
        M_zeroLocVectorialTranspose.setZero();
        M_zeroLocTensor2.setZero();
    }
    tensorLevelsetDeltaBase( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcPhi( this->createPcPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_ctxPhi( this->createCtxPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_locRes( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPhi( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_isSameMesh( fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo(expr.phi().functionSpace()->mesh()) && isSameGeo ),
            M_zeroLocScalar(),
            M_zeroLocVectorial( ),
            M_zeroLocVectorialTranspose(),
            M_zeroLocTensor2()
    {
        if( !M_isSameMesh )
            expr.phi().functionSpace()->mesh()->tool_localization()->updateForUse();

        M_zeroLocScalar.setZero();
        M_zeroLocVectorial.setZero();
        M_zeroLocVectorialTranspose.setZero();
        M_zeroLocTensor2.setZero();
    }
    tensorLevelsetDeltaBase( expr_type const& expr, Geo_t const& geom )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcPhi( this->createPcPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_ctxPhi( this->createCtxPhiIfSameGeo(expr, geom, mpl::bool_<isSameGeo>()) ),
            M_locRes( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPhi( expr.phi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_isSameMesh( fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo(expr.phi().functionSpace()->mesh()) && isSameGeo ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocVectorialTranspose(),
            M_zeroLocTensor2()
    {
        if( !M_isSameMesh )
            expr.phi().functionSpace()->mesh()->tool_localization()->updateForUse();

        M_zeroLocScalar.setZero();
        M_zeroLocVectorial.setZero();
        M_zeroLocVectorialTranspose.setZero();
        M_zeroLocTensor2.setZero();
    }

    virtual ~tensorLevelsetDeltaBase() {}

    expr_type const& expr() const { return M_expr; }

    void setGmc( Geo_t const& geom ) { M_geot = fusion::at_key<key_type>( geom ); }
    gmc_ptrtype const& gmc() const { return M_geot; }

    ctx_phi_ptrtype const& ctxPhi() const { return M_ctxPhi; }

    array_shape_type & locRes() { return M_locRes; }
    array_shape_type const& locRes() const { return M_locRes; }
    matrix_shape_type & locRes( uint16_type q ) { return M_locRes[q]; }
    matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }

    array_scalar_type const& locPhi() const { return M_locPhi; }
    loc_scalar_type const& locPhi( uint16_type q ) const { return M_locPhi[q]; }

    virtual void update( Geo_t const& geom )
    {
        this->updateCtxPhiIfSameGeo( geom, mpl::bool_<isSameGeo>() );
        std::fill( M_locPhi.data(), M_locPhi.data()+M_locPhi.num_elements(), this->M_zeroLocScalar );
        if( M_isSameMesh )
        {
            this->expr().phi().id( *M_ctxPhi, M_locPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<1>());
            auto setOfPts = ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) ?
                fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() :
                fusion::at_key<key_type>( geom )->element().vertices();
            this->expr().phi().idInterpolate( ptsreal, M_locPhi, this->expr().useInterpWithConfLoc(), setOfPts );
        }
    }

    virtual void update( Geo_t const& geom, uint16_type face )
    {
        this->updateCtxFacePhiIfSameGeo( geom, mpl::bool_<isSameGeo>() );
        std::fill( M_locPhi.data(), M_locPhi.data()+M_locPhi.num_elements(), this->M_zeroLocScalar );
        if( M_isSameMesh )
        {
            this->expr().phi().id( *M_ctxPhi, M_locPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<2>());
            this->expr().phi().idInterpolate( ptsreal, M_locPhi, this->expr().useInterpWithConfLoc(), fusion::at_key<key_type>(geom)->element().face( fusion::at_key<key_type>(geom)->faceId() ).vertices() );
        }
    }

    virtual matrix_shape_type const& evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return this->evalq( q );
    }
    virtual value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return this->evalq( c1,c2,q );
    }
    virtual value_type evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return this->evalq( c1,c2,q );
    }
    virtual matrix_shape_type const& evaliq( uint16_type i, uint16_type q ) const
    {
        return this->evalq( q );
    }
    virtual value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return M_locRes[q](0,0);
    }
    virtual matrix_shape_type const& evalq( uint16_type q ) const
    {
        return M_locRes[q];
    }

private:
    pc_phi_ptrtype createPcPhiIfSameGeo( expr_type const& expr, Geo_t const& geom, mpl::bool_<true> )
    {
        if( expr.phi().functionSpace() && expr.phi().functionSpace()->fe() )
            return pc_phi_ptrtype( new pc_phi_type( expr.phi().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) );
        else
            return pc_phi_ptrtype();
    }
    pc_phi_ptrtype createPcPhiIfSameGeo( expr_type const& expr, Geo_t const& geom, mpl::bool_<false> )
    {
        return pc_phi_ptrtype();
    }
    ctx_phi_ptrtype createCtxPhiIfSameGeo( expr_type const& expr, Geo_t const& geom, mpl::bool_<true> )
    {
        if( expr.phi().functionSpace() && expr.phi().functionSpace()->fe() )
            return ctx_phi_ptrtype( new ctx_phi_type( expr.phi().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_phi_ptrtype const&)M_pcPhi ) );
        else
            return ctx_phi_ptrtype();
    }
    ctx_phi_ptrtype createCtxPhiIfSameGeo( expr_type const& expr, Geo_t const& geom, mpl::bool_<false> )
    {
        return ctx_phi_ptrtype();
    }

    void updateCtxPhiIfSameGeo( Geo_t const& geom, mpl::bool_<true> )
    {
        if (fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcPhi->update(fusion::at_key<key_type>( geom )->pc()->nodes() );
        M_ctxPhi->update( fusion::at_key<key_type>( geom ),  (pc_phi_ptrtype const&) M_pcPhi );
    }
    void updateCtxPhiIfSameGeo( Geo_t const& geom, mpl::bool_<false> )
    {
    }
    void updateCtxFacePhiIfSameGeo( Geo_t const& geom, mpl::bool_<true> )
    {
        M_pcPhi->update(fusion::at_key<key_type>( geom )->pc()->nodes() );
        M_ctxPhi->update( fusion::at_key<key_type>( geom ), M_pcPhi );
    }
    void updateCtxFacePhiIfSameGeo( Geo_t const& geom, mpl::bool_<false> )
    {
    }

private:
    expr_type const& M_expr;

    gmc_ptrtype M_geot;

    pc_phi_ptrtype M_pcPhi;
    ctx_phi_ptrtype M_ctxPhi;

    array_shape_type M_locRes;
    array_scalar_type M_locPhi;

protected:
    bool M_isSameMesh;

    loc_scalar_type M_zeroLocScalar;
    loc_vectorial_type M_zeroLocVectorial;
    loc_vectorial_transpose_type M_zeroLocVectorialTranspose;
    loc_tensor2_type M_zeroLocTensor2;

    static Eigen::array<Eigen::IndexPair<int>, 1> S_tensorVectorialDotContraction;
    static Eigen::array<Eigen::IndexPair<int>, 1> S_tensorVectorialTransposeDotContraction;
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
Eigen::array<Eigen::IndexPair<int>, 1> 
tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>::S_tensorVectorialDotContraction = {{ Eigen::IndexPair<int>(0,0) }};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
Eigen::array<Eigen::IndexPair<int>, 1> 
tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>::S_tensorVectorialTransposeDotContraction = {{ Eigen::IndexPair<int>(1,1) }};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct NonAdaptiveThicknessPolicy
{
public:
    typedef tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> tensorbase_type;
    typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;

    typedef typename tensorbase_type::value_type value_type;
    typedef typename tensorbase_type::loc_scalar_type loc_scalar_type;

    NonAdaptiveThicknessPolicy( tensorbase_type const& tensor )
        :
            M_tensor( tensor ),
            M_locEpsilon()
    {
        M_locEpsilon(0,0) = tensor.expr().thickness();
    }

    void update( Geo_t const& geom ) {}
    void update( Geo_t const& geom, uint16_type face ) {}

    loc_scalar_type const& locEpsilon( uint16_type q ) const { return M_locEpsilon; }

private:
    tensorbase_type const& M_tensor;
    loc_scalar_type M_locEpsilon;
};

/** Generic case:
 * D = 1/(2*eps0) * ( 1 + cos(pi*phi/eps0) ) * |grad(phi)| if -eps0 <= phi <= eps0
 *   = 0 otherwise
 */
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType,
    template<typename, typename, typename, typename> class ThicknessPolicy>
struct tensorLevelsetDeltaGeneric
: 
    public tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>,
    public ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
    typedef ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType> thickness_policy_type;
public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::key_type key_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::shape_type shape_type;

    typedef typename super_type::matrix_node_type matrix_node_type;

    // array
    typedef typename super_type::loc_scalar_type loc_scalar_type;
    typedef typename super_type::array_scalar_type array_scalar_type;
    typedef typename super_type::loc_vectorial_type loc_vectorial_type;
    typedef typename super_type::array_vectorial_type array_vectorial_type;
    typedef typename super_type::loc_vectorial_transpose_type loc_vectorial_transpose_type;
    typedef typename super_type::array_vectorial_transpose_type array_vectorial_transpose_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::array_shape_type array_shape_type;
    typedef typename super_type::vector_shape_type vector_shape_type;

    using thickness_policy_type::locEpsilon;

    tensorLevelsetDeltaGeneric( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            super_type( expr, geom, fev, feu ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorLevelsetDeltaGeneric( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            super_type( expr, geom, fev ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorLevelsetDeltaGeneric( expr_type const& expr, Geo_t const& geom )
        :
            super_type( expr, geom ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}

    array_vectorial_transpose_type const& locGradPhi() const { return M_locGradPhi; }
    loc_vectorial_transpose_type const& locGradPhi( uint16_type q ) const { return M_locGradPhi[q]; }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        thickness_policy_type::update( geom );
        std::fill( M_locGradPhi.data(), M_locGradPhi.data()+M_locGradPhi.num_elements(), this->M_zeroLocVectorialTranspose );
        if( this->M_isSameMesh )
        {
            this->expr().phi().grad( *this->ctxPhi(), M_locGradPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<1>());
            auto setOfPts = ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) ?
                fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() :
                fusion::at_key<key_type>( geom )->element().vertices();
            this->expr().phi().gradInterpolate( ptsreal, M_locGradPhi, this->expr().useInterpWithConfLoc(), setOfPts );
        }

        this->updateImpl();
    }

    virtual void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );
        thickness_policy_type::update( geom, face );
        std::fill( M_locGradPhi.data(), M_locGradPhi.data()+M_locGradPhi.num_elements(), this->M_zeroLocVectorialTranspose );
        if( this->M_isSameMesh )
        {
            this->expr().phi().grad( *this->ctxPhi(), M_locGradPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<2>());
            this->expr().phi().gradInterpolate( ptsreal, M_locGradPhi, this->expr().useInterpWithConfLoc(), fusion::at_key<key_type>(geom)->element().face( fusion::at_key<key_type>(geom)->faceId() ).vertices() );
        }

        this->updateImpl();
    }

protected:
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            value_type const phi = this->locPhi(q)(0,0);
            value_type const eps = this->locEpsilon(q)(0,0);
            if( (phi < -eps) || (phi > eps) )
            {
                this->locRes(q)(0,0) = 0.;
            }
            else
            {
                auto const& gradPhi = this->locGradPhi(q);
                loc_scalar_type tensorGradPhi2 = gradPhi.contract( gradPhi, this->S_tensorVectorialTransposeDotContraction );
                value_type modGradPhi = math::sqrt( tensorGradPhi2(0) );
                this->locRes(q)(0,0) = 0.5 / eps * ( 1 + math::cos(M_PI*phi/eps) ) * modGradPhi;
            }
        }
    }

private:
    array_vectorial_transpose_type M_locGradPhi;
};

/** Localredist (delta(phi/modgradphi)) case:
 * D = 1/(2*eps0) * ( 1 + cos(pi*psi/eps0) ) if -eps0 <= psi <= eps0 with psi = phi / modgradphi
 *   = 0 otherwise
 */
template< typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType,
    template<typename, typename, typename, typename> class ThicknessPolicy>
struct tensorLevelsetDeltaLocalRedist
: 
    public tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>,
    public ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
    typedef ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType> thickness_policy_type;
public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::key_type key_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::shape_type shape_type;

    typedef typename super_type::matrix_node_type matrix_node_type;

    // array
    typedef typename super_type::loc_scalar_type loc_scalar_type;
    typedef typename super_type::array_scalar_type array_scalar_type;
    typedef typename super_type::loc_vectorial_type loc_vectorial_type;
    typedef typename super_type::array_vectorial_type array_vectorial_type;
    typedef typename super_type::loc_vectorial_transpose_type loc_vectorial_transpose_type;
    typedef typename super_type::array_vectorial_transpose_type array_vectorial_transpose_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::array_shape_type array_shape_type;
    typedef typename super_type::vector_shape_type vector_shape_type;

    using thickness_policy_type::locEpsilon;

    tensorLevelsetDeltaLocalRedist( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            super_type( expr, geom, fev, feu ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorLevelsetDeltaLocalRedist( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            super_type( expr, geom, fev ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorLevelsetDeltaLocalRedist( expr_type const& expr, Geo_t const& geom )
        :
            super_type( expr, geom ),
            thickness_policy_type( (super_type const&)(*this) ),
            M_locGradPhi( expr.phi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}

    array_vectorial_transpose_type const& locGradPhi() const { return M_locGradPhi; }
    loc_vectorial_transpose_type const& locGradPhi( uint16_type q ) const { return M_locGradPhi[q]; }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        thickness_policy_type::update( geom );
        std::fill( M_locGradPhi.data(), M_locGradPhi.data()+M_locGradPhi.num_elements(), this->M_zeroLocVectorialTranspose );
        if( this->M_isSameMesh )
        {
            this->expr().phi().grad( *this->ctxPhi(), M_locGradPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<1>());
            auto setOfPts = ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) ?
                fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() :
                fusion::at_key<key_type>( geom )->element().vertices();
            this->expr().phi().gradInterpolate( ptsreal, M_locGradPhi, this->expr().useInterpWithConfLoc(), setOfPts );
        }

        this->updateImpl();
    }

    virtual void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );
        thickness_policy_type::update( geom, face );
        std::fill( M_locGradPhi.data(), M_locGradPhi.data()+M_locGradPhi.num_elements(), this->M_zeroLocVectorialTranspose );
        if( this->M_isSameMesh )
        {
            this->expr().phi().grad( *this->ctxPhi(), M_locGradPhi );
        }
        else
        {
            matrix_node_type ptsreal = this->expr().phi().ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<2>());
            this->expr().phi().gradInterpolate( ptsreal, M_locGradPhi, this->expr().useInterpWithConfLoc(), fusion::at_key<key_type>(geom)->element().face( fusion::at_key<key_type>(geom)->faceId() ).vertices() );
        }

        this->updateImpl();
    }

protected:
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            value_type const eps = this->locEpsilon(q)(0,0);
            value_type const phi = this->locPhi(q)(0,0);
            auto const& gradPhi = this->locGradPhi(q);
            loc_scalar_type tensorGradPhi2 = gradPhi.contract( gradPhi, this->S_tensorVectorialTransposeDotContraction );
            value_type modGradPhi = math::sqrt( tensorGradPhi2(0) );
            value_type const psi = phi / modGradPhi;
            if( (psi < -eps) || (psi > eps) )
            {
                this->locRes(q)(0,0) = 0.;
            }
            else
            {
                this->locRes(q)(0,0) = 0.5 / eps * ( 1 + math::cos(M_PI*psi/eps) );
            }
        }
    }

private:
    array_vectorial_transpose_type M_locGradPhi;
};

/** Non-adaptive thickness distance (|grad(phi)| = 1) case:
 * D = 1/(2*eps0) * ( 1 + cos(pi*phi/eps0) ) if -eps0 <= phi <= eps0
 *   = 0 otherwise
 */
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType,
    template<typename, typename, typename, typename> class ThicknessPolicy>
struct tensorLevelsetDeltaDistance
: 
    public tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>,
    public ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
    typedef ThicknessPolicy<Geo_t, Basis_i_t, Basis_j_t, ExprType> thickness_policy_type;
public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::key_type key_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::shape_type shape_type;

    // array
    typedef typename super_type::loc_scalar_type loc_scalar_type;
    typedef typename super_type::array_scalar_type array_scalar_type;
    typedef typename super_type::loc_vectorial_type loc_vectorial_type;
    typedef typename super_type::array_vectorial_type array_vectorial_type;
    typedef typename super_type::loc_vectorial_transpose_type loc_vectorial_transpose_type;
    typedef typename super_type::array_vectorial_transpose_type array_vectorial_transpose_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::array_shape_type array_shape_type;
    typedef typename super_type::vector_shape_type vector_shape_type;

    using thickness_policy_type::locEpsilon;

    tensorLevelsetDeltaDistance( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            super_type( expr, geom, fev, feu ),
            thickness_policy_type( (super_type const&)(*this) )
    {}
    tensorLevelsetDeltaDistance( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            super_type( expr, geom, fev ),
            thickness_policy_type( (super_type const&)(*this) )
    {}
    tensorLevelsetDeltaDistance( expr_type const& expr, Geo_t const& geom )
        :
            super_type( expr, geom ),
            thickness_policy_type( (super_type const&)(*this) )
    {}

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        thickness_policy_type::update( geom );

        this->updateImpl();
    }

    void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );
        thickness_policy_type::update( geom, face );

        this->updateImpl();
    }

protected:
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            value_type const phi = this->locPhi(q)(0,0);
            value_type const eps = this->locEpsilon(q)(0,0);
            if( (phi < -eps) || (phi > eps) )
                this->locRes(q)(0,0) = 0.;
            else
                this->locRes(q)(0,0) = 0.5 / eps * ( 1 + math::cos(M_PI*phi/eps) );
        }
    }
};

template<typename ElementPhiType>
class LevelsetDeltaExpr
{
public:
    typedef LevelsetDeltaExpr<ElementPhiType> this_type;

    typedef ElementPhiType element_phi_type;

    static const size_type context_phi = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_phi;

    //------------------------------------------------------------------------------//
    // implementation
    enum class DeltaImpl {
        GENERIC, RENORMALISED, DISTANCE
    };

    //------------------------------------------------------------------------------//
    // phi functionspace
    typedef typename element_phi_type::functionspace_type functionspace_phi_type;
    typedef typename functionspace_phi_type::reference_element_type fe_phi_type;
    //------------------------------------------------------------------------------//
    // expression desc
    typedef typename functionspace_phi_type::geoelement_type geoelement_type;
    typedef typename functionspace_phi_type::value_type value_type;
    typedef value_type evaluate_type;
    static const uint16_type nDim = fe_phi_type::nDim;
    static const uint16_type nRealDim = fe_phi_type::nRealDim;
    static const uint16_type rank = fe_phi_type::rank;
    static const uint16_type nComponents1 = fe_phi_type::nComponents1;
    static const uint16_type nComponents2 = fe_phi_type::nComponents2;
    static const bool is_terminal = false;

    static const uint16_type orderphi = functionspace_phi_type::basis_type::nOrder;

    // imorder
    static const uint16_type imorderAuto = 6*orderphi - 4;

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
    LevelsetDeltaExpr( element_phi_type const& phi, double thickness ) : 
        M_phi( boost::cref(phi) ),
        M_thicknessDelta( thickness ),
        M_useAdaptiveThickness( false ),
        M_impl( DeltaImpl::GENERIC ),
        M_imOrder( -1 )
    {}

    ~LevelsetDeltaExpr() = default;

    //--------------------------------------------------------------------//
    //! polynomial order
    uint16_type polynomialOrder() const { return (M_imOrder>=0) ? M_imOrder : imorderAuto; }
    void setPolynomialOrder( int imorder ) { M_imOrder = imorder; }

    //! expression is polynomial?
    bool isPolynomial() const { return false; }

    //--------------------------------------------------------------------//
    bool useInterpWithConfLoc() const { return false; }

    //--------------------------------------------------------------------//
    // Accessors
    element_phi_type const& phi() const { return M_phi; }
    double thickness() const { return M_thicknessDelta; }

    //--------------------------------------------------------------------//
    // Options
    bool useAdaptiveThickness() const { return M_useAdaptiveThickness; }
    void setUseAdaptiveThickness( bool b ) { M_useAdaptiveThickness = b; }
    DeltaImpl const& impl() const { return M_impl; }
    void setImpl( DeltaImpl const& impl ) { M_impl = impl; }

    //--------------------------------------------------------------------//
    // Expr tensor
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expr_type;
        typedef typename element_phi_type::value_type value_type;

        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        // tensor base
        typedef tensorLevelsetDeltaBase<Geo_t, Basis_i_t, Basis_j_t, expr_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;
        // shape
        typedef typename tensorbase_type::shape shape;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        // is zero
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if( expr.useAdaptiveThickness() )
            {
                // TODO
            }
            else /* !expr.useAdaptiveThickness() */
            {
                switch( expr.impl() )
                {
                    case DeltaImpl::GENERIC:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom, fev, feu
                                    )
                                );
                    } break;
                    case DeltaImpl::RENORMALISED:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaLocalRedist< Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy >( 
                                    expr, geom, fev, feu
                                    )
                                );
                    } break;
                    case DeltaImpl::DISTANCE:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom, fev, feu
                                    )
                                );
                    } break;
                }
            }
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            if( expr.useAdaptiveThickness() )
            {
                // TODO
            }
            else /* !expr.useAdaptiveThickness() */
            {
                switch( expr.impl() )
                {
                    case DeltaImpl::GENERIC:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom, fev
                                    )
                                );
                    } break;
                    case DeltaImpl::RENORMALISED:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaLocalRedist< Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy >( 
                                    expr, geom, fev
                                    )
                                );
                    } break;
                    case DeltaImpl::DISTANCE:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom, fev
                                    )
                                );
                    } break;
                }
            }
        }

        tensor( this_type const& expr, Geo_t const& geom )
        {
            if( expr.useAdaptiveThickness() )
            {
                // TODO
            }
            else /* !expr.useAdaptiveThickness() */
            {
                switch( expr.impl() )
                {
                    case DeltaImpl::GENERIC:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom
                                    )
                                );
                    } break;
                    case DeltaImpl::RENORMALISED:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaLocalRedist< Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy >( 
                                    expr, geom
                                    )
                                );
                    } break;
                    case DeltaImpl::DISTANCE:
                    {
                        M_tensorbase.reset(
                                new tensorLevelsetDeltaDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type, NonAdaptiveThicknessPolicy>( 
                                    expr, geom
                                    )
                                );
                    } break;
                }
            }
        }

        template<typename IM>
        void init( IM const& im )
        {}

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
            M_tensorbase->update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensorbase->update( geom, face );
        }

        // eval
        matrix_shape_type const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,c1,c2,q );
        }
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i,c1,c2,q );
        }
        matrix_shape_type const&
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_tensorbase->evalq( q );
        }

    private:
        tensorbase_ptrtype M_tensorbase;
    };

private:
    boost::reference_wrapper<const element_phi_type> M_phi;
    double M_thicknessDelta;
    bool M_useAdaptiveThickness;
    DeltaImpl M_impl;
    int M_imOrder;
};

namespace detail {
template<typename Args>
struct compute_levelsetDelta_return
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<Args, Feel::tag::element>::type
    >::type
    >::type
    >::type element_type;
};
}

BOOST_PARAMETER_FUNCTION(
        ( Expr< LevelsetDeltaExpr<typename detail::compute_levelsetDelta_return<Args>::element_type> > ), // return type
        levelsetDelta, // function name
        tag, // tag types namespace
        ( required
          ( element, * )
          ( thickness, (double) )
        ) // required parameters
        ( optional
          ( use_adaptive_thickness, (bool), false )
          ( impl, (std::string), "generic" )
          ( imorder, (int), -1 )
        ) // optional parameters
    )
{
    typedef typename detail::compute_levelsetDelta_return<Args>::element_type ElementPhiType;
    typedef LevelsetDeltaExpr<ElementPhiType> lsdelta_t;
    lsdelta_t levelsetDeltaExpr( element, thickness );

    levelsetDeltaExpr.setUseAdaptiveThickness( use_adaptive_thickness );
    if ( impl == "generic" )
        levelsetDeltaExpr.setImpl( lsdelta_t::DeltaImpl::GENERIC );
    else if ( impl == "renorm" )
        levelsetDeltaExpr.setImpl( lsdelta_t::DeltaImpl::RENORMALISED );
    else if ( impl == "distance" )
        levelsetDeltaExpr.setImpl( lsdelta_t::DeltaImpl::DISTANCE );
    else
        CHECK( false ) << impl << " is not a valid DeltaImpl (available: generic, renorm, distance)";

    levelsetDeltaExpr.setPolynomialOrder( imorder );
    
    return Expr< lsdelta_t >( levelsetDeltaExpr );
}

}
}

#endif
