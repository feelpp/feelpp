#ifndef _HELFRICH_FORCE_EXPR_HPP
#define _HELFRICH_FORCE_EXPR_HPP 1

namespace Feel {
namespace FeelModels {

enum class HelfrichInnerDivImplementation
{
    GENERIC, DISTANCE
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
struct tensorHelfrichInnerDivTensorBase
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

    typedef typename expr_type::fe_normal_type fe_normal_type;
    typedef typename expr_type::fe_curvature_type fe_curvature_type;
    typedef typename expr_type::geoelement_type geoelement_type;

    // fe normal context
    typedef typename fe_normal_type::PreCompute pc_normal_type;
    typedef std::shared_ptr<pc_normal_type> pc_normal_ptrtype;
    typedef typename fe_normal_type::template Context<expr_type::context_normal, fe_normal_type, gm_type, geoelement_type, 0, gmc_type::subEntityCoDim> ctx_normal_type;
    typedef std::shared_ptr<ctx_normal_type> ctx_normal_ptrtype;
    // fe curvature context
    typedef typename fe_curvature_type::PreCompute pc_curvature_type;
    typedef std::shared_ptr<pc_curvature_type> pc_curvature_ptrtype;
    typedef typename fe_curvature_type::template Context<expr_type::context_curvature, fe_curvature_type, gm_type, geoelement_type, 0, gmc_type::subEntityCoDim> ctx_curvature_type;
    typedef std::shared_ptr<ctx_curvature_type> ctx_curvature_ptrtype;

    // Shapes types
    // res shape
    typedef Shape<expr_type::nRealDim, Vectorial, false, false> shape_type;
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


    tensorHelfrichInnerDivTensorBase( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
            M_pcCurvature( new pc_curvature_type( expr.curvature().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxCurvature( new ctx_curvature_type( expr.curvature().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_curvature_ptrtype const&)M_pcCurvature ) ),
            M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locCurvature( expr.curvature().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradCurvature( expr.curvature().gradExtents(*fusion::at_key<key_type>( geom )) ),
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
    tensorHelfrichInnerDivTensorBase( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
            M_pcCurvature( new pc_curvature_type( expr.curvature().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxCurvature( new ctx_curvature_type( expr.curvature().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_curvature_ptrtype const&)M_pcCurvature ) ),
            M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locCurvature( expr.curvature().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradCurvature( expr.curvature().gradExtents(*fusion::at_key<key_type>( geom )) ),
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
    tensorHelfrichInnerDivTensorBase( expr_type const& expr, Geo_t const& geom )
        :
            M_expr( expr ),
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_pcNormal( new pc_normal_type( expr.normal().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxNormal( new ctx_normal_type( expr.normal().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_normal_ptrtype const&)M_pcNormal ) ),
            M_pcCurvature( new pc_curvature_type( expr.curvature().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxCurvature( new ctx_curvature_type( expr.curvature().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_curvature_ptrtype const&)M_pcCurvature ) ),
            M_locRes( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locNormal( expr.normal().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locCurvature( expr.curvature().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradCurvature( expr.curvature().gradExtents(*fusion::at_key<key_type>( geom )) ),
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

    virtual ~tensorHelfrichInnerDivTensorBase() {}

    expr_type const& expr() const { return M_expr; }

    void setGmc( Geo_t const& geom ) { M_geot = fusion::at_key<key_type>( geom ); }
    gmc_ptrtype const& gmc() const { return M_geot; }

    array_shape_type & locRes() { return M_locRes; }
    array_shape_type const& locRes() const { return M_locRes; }
    matrix_shape_type & locRes( uint16_type q ) { return M_locRes[q]; }
    matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }

    array_vectorial_type const& locNormal() const { return M_locNormal; }
    loc_vectorial_type const& locNormal( uint16_type q ) const { return M_locNormal[q]; }
    array_scalar_type const& locCurvature() const { return M_locCurvature; }
    loc_scalar_type const& locCurvature( uint16_type q ) const { return M_locCurvature[q]; }
    array_vectorial_transpose_type const& locGradCurvature() const { return M_locGradCurvature; }
    loc_vectorial_transpose_type const& locGradCurvature( uint16_type q ) const { return M_locGradCurvature[q]; }

    virtual void update( Geo_t const& geom )
    {
        M_ctxNormal->update( fusion::at_key<key_type>( geom ),  (pc_normal_ptrtype const&) M_pcNormal );
        std::fill( M_locNormal.data(), M_locNormal.data()+M_locNormal.num_elements(), this->M_zeroLocVectorial );
        this->expr().normal().id( *M_ctxNormal, M_locNormal );

        M_ctxCurvature->update( fusion::at_key<key_type>( geom ),  (pc_curvature_ptrtype const&) M_pcCurvature );
        std::fill( M_locCurvature.data(), M_locCurvature.data()+M_locCurvature.num_elements(), this->M_zeroLocScalar );
        this->expr().curvature().id( *M_ctxCurvature, M_locCurvature );
        std::fill( M_locGradCurvature.data(), M_locGradCurvature.data()+M_locGradCurvature.num_elements(), this->M_zeroLocVectorialTranspose );
        this->expr().curvature().grad( *M_ctxCurvature, M_locGradCurvature );
    }

    virtual void update( Geo_t const& geom, uint16_type face )
    {
        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
        {
            M_pcNormal->update( this->gmc()->pc()->nodes() );
            M_pcCurvature->update( this->gmc()->pc()->nodes() );
        }
        this->update( geom );
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
        return M_locRes[q](c1,c2);
    }
    virtual matrix_shape_type const& evalq( uint16_type q ) const
    {
        return M_locRes[q];
    }

private:
    expr_type const& M_expr;

    gmc_ptrtype M_geot;

    pc_normal_ptrtype M_pcNormal;
    ctx_normal_ptrtype M_ctxNormal;
    pc_curvature_ptrtype M_pcCurvature;
    ctx_curvature_ptrtype M_ctxCurvature;

    array_shape_type M_locRes;
    array_vectorial_type M_locNormal;
    array_scalar_type M_locCurvature;
    array_vectorial_transpose_type M_locGradCurvature;

protected:
    loc_scalar_type M_zeroLocScalar;
    loc_vectorial_type M_zeroLocVectorial;
    loc_vectorial_transpose_type M_zeroLocVectorialTranspose;
    loc_tensor2_type M_zeroLocTensor2;
};

/**
 * Generic Helfrich force div argument:
 * f = - K^2 / 2 + 1/|gradPhi| * (1-NxN) grad( |gradPhi|*K )
 */
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
struct tensorHelfrichInnerDivTensorGeneric
: public tensorHelfrichInnerDivTensorBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorHelfrichInnerDivTensorBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::key_type key_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::shape_type shape_type;

    // fe modgradphi context
    typedef typename expr_type::fe_modgradphi_type fe_modgradphi_type;
    typedef typename fe_modgradphi_type::PreCompute pc_modgradphi_type;
    typedef std::shared_ptr<pc_modgradphi_type> pc_modgradphi_ptrtype;
    typedef typename fe_modgradphi_type::template Context<expr_type::context_modgradphi, fe_modgradphi_type, gm_type, geoelement_type, 0, gmc_type::subEntityCoDim> ctx_modgradphi_type;
    typedef std::shared_ptr<ctx_modgradphi_type> ctx_modgradphi_ptrtype;

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

    tensorHelfrichInnerDivTensorGeneric( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            super_type( expr, geom, fev, feu ),
            M_pcModGradPhi( new pc_modgradphi_type( expr.modgradphi().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxModGradPhi( new ctx_modgradphi_type( expr.modgradphi().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_modgradphi_ptrtype const&)M_pcModGradPhi ) ),
            M_locModGradPhi( expr.modgradphi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradModGradPhi( expr.modgradphi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorHelfrichInnerDivTensorGeneric( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            super_type( expr, geom, fev ),
            M_pcModGradPhi( new pc_modgradphi_type( expr.modgradphi().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxModGradPhi( new ctx_modgradphi_type( expr.modgradphi().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_modgradphi_ptrtype const&)M_pcModGradPhi ) ),
            M_locModGradPhi( expr.modgradphi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradModGradPhi( expr.modgradphi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}
    tensorHelfrichInnerDivTensorGeneric( expr_type const& expr, Geo_t const& geom )
        :
            super_type( expr, geom ),
            M_pcModGradPhi( new pc_modgradphi_type( expr.modgradphi().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxModGradPhi( new ctx_modgradphi_type( expr.modgradphi().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_modgradphi_ptrtype const&)M_pcModGradPhi ) ),
            M_locModGradPhi( expr.modgradphi().idExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradModGradPhi( expr.modgradphi().gradExtents(*fusion::at_key<key_type>( geom )) )
    {}

    array_scalar_type const& locModGradPhi() const { return M_locModGradPhi; }
    loc_scalar_type const& locModGradPhi( uint16_type q ) const { return M_locModGradPhi[q]; }
    array_vectorial_transpose_type const& locGradModGradPhi() const { return M_locGradModGradPhi; }
    loc_vectorial_transpose_type const& locGradModGradPhi( uint16_type q ) const { return M_locGradModGradPhi[q]; }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );

        M_ctxModGradPhi->update( fusion::at_key<key_type>( geom ),  (pc_modgradphi_ptrtype const&) M_pcModGradPhi );
        std::fill( M_locModGradPhi.data(), M_locModGradPhi.data()+M_locModGradPhi.num_elements(), this->M_zeroLocScalar );
        this->expr().modgradphi().id( *M_ctxModGradPhi, M_locModGradPhi );
        std::fill( M_locGradModGradPhi.data(), M_locGradModGradPhi.data()+M_locGradModGradPhi.num_elements(), this->M_zeroLocVectorialTranspose );
        this->expr().modgradphi().grad( *M_ctxModGradPhi, M_locGradModGradPhi );

        this->updateImpl();
    }

protected:
    void updateImpl()
    {
        Eigen::array<Eigen::IndexPair<int>, 1> dot_contraction = {{ Eigen::IndexPair<int>(0,1) }};
        Eigen::array<int, 2> transpose_shuffle = {{1, 0}};

        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            auto const& N = this->locNormal(q);
            value_type const K = this->locCurvature(q)(0,0);
            auto const& gradK = this->locGradCurvature(q);
            value_type const modGradPhi = this->locModGradPhi(q)(0,0);
            auto const& gradModGradPhi = this->locGradModGradPhi(q);

            value_type const K2_2 = 0.5*K*K;

            auto gradModGradPhiK_modGradPhi = gradModGradPhi*K/modGradPhi + gradK;
            loc_scalar_type tensorNdotGradModGradPhiK_modGradPhi = N.contract( gradModGradPhiK_modGradPhi, dot_contraction );
            value_type const NdotGradModGradPhiK_modGradPhi = tensorNdotGradModGradPhiK_modGradPhi(0);

            loc_vectorial_type res = -(NdotGradModGradPhiK_modGradPhi + K2_2)*N + gradModGradPhiK_modGradPhi.shuffle(transpose_shuffle);
            this->locRes(q) = Eigen::Map<const matrix_shape_type>( res.data(), shape_type::M, shape_type::N );

            //vector_shape_type N = Eigen::Map<const vector_shape_type>( this->locNormal(q).data(), shape_type::M, 1 );
            //value_type const K = this->locCurvature(q)(0,0);
            //value_type const K2_2 = 0.5*K*K;
            //loc_vectorial_type locGradK = this->locGradCurvature(q).shuffle( transpose_shuffle );
            //vector_shape_type gradK = Eigen::Map<const vector_shape_type>( locGradK.data(), shape_type::M, 1 );
            //value_type NdotGradK = ( N.transpose() * gradK )(0);
            //this->locRes(q) = - K2_2*N + (Eigen::Matrix<value_type, shape_type::M, shape_type::M>::Identity() - N*N.transpose()) * gradK;

        }
    }

private:
    pc_modgradphi_ptrtype M_pcModGradPhi;
    ctx_modgradphi_ptrtype M_ctxModGradPhi;
    array_scalar_type M_locModGradPhi;
    array_vectorial_transpose_type M_locGradModGradPhi;
};

/**
 * Distance-case (|gradPhi|=1) Helfrich force div argument:
 * f = - K^2 / 2 + 1/|gradPhi| * (1-NxN) grad( |gradPhi|*K )
 */
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
struct tensorHelfrichInnerDivTensorDistance
: public tensorHelfrichInnerDivTensorBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorHelfrichInnerDivTensorBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::key_type key_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::shape_type shape_type;

    // fe modgradphi context
    typedef typename expr_type::fe_modgradphi_type fe_modgradphi_type;
    typedef typename fe_modgradphi_type::PreCompute pc_modgradphi_type;
    typedef std::shared_ptr<pc_modgradphi_type> pc_modgradphi_ptrtype;
    typedef typename fe_modgradphi_type::template Context<expr_type::context_modgradphi, fe_modgradphi_type, gm_type, geoelement_type, 0, gmc_type::subEntityCoDim> ctx_modgradphi_type;
    typedef std::shared_ptr<ctx_modgradphi_type> ctx_modgradphi_ptrtype;

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

    tensorHelfrichInnerDivTensorDistance( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
            super_type( expr, geom, fev, feu )
    {}
    tensorHelfrichInnerDivTensorDistance( expr_type const& expr,
            Geo_t const& geom, Basis_i_t const& fev )
        :
            super_type( expr, geom, fev )
    {}
    tensorHelfrichInnerDivTensorDistance( expr_type const& expr, Geo_t const& geom )
        :
            super_type( expr, geom )
    {}

    void update( Geo_t const& geom )
    {
        super_type::update( geom );

        this->updateImpl();
    }

protected:
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            vector_shape_type N = Eigen::Map<const vector_shape_type>( this->locNormal(q).data(), shape_type::M, 1 );
            value_type const K = this->locCurvature(q)(0,0);
            value_type const K2_2 = 0.5*K*K;
            vector_shape_type gradK = Eigen::Map<const vector_shape_type>( this->locGradCurvature(q).data(), shape_type::M, 1 );
            this->locRes(q) = - K2_2*N + (Eigen::Matrix<value_type, shape_type::M, shape_type::M>::Identity() - N*N.transpose()) * gradK;

        }
    }
};

template<typename ElementNormalType, typename ElementCurvatureType, typename ElementModGradPhiType, int IMOrder>
class HelfrichInnerDivImpl
{
public:
    typedef HelfrichInnerDivImpl<ElementNormalType, ElementCurvatureType, ElementModGradPhiType, IMOrder> this_type;

    static const size_type context_normal = vm::JACOBIAN;
    static const size_type context_curvature = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_modgradphi = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_curvature;

    typedef ElementNormalType element_normal_type;
    typedef ElementCurvatureType element_curvature_type;
    typedef ElementModGradPhiType element_modgradphi_type;

    //------------------------------------------------------------------------------//
    // normal functionspace
    typedef typename element_normal_type::functionspace_type functionspace_normal_type;
    typedef typename functionspace_normal_type::reference_element_type fe_normal_type;
    //------------------------------------------------------------------------------//
    // curvature functionspace
    typedef typename element_curvature_type::functionspace_type functionspace_curvature_type;
    typedef typename functionspace_curvature_type::reference_element_type fe_curvature_type;
    //------------------------------------------------------------------------------//
    // modgradphi functionspace
    typedef typename element_modgradphi_type::functionspace_type functionspace_modgradphi_type;
    typedef typename functionspace_modgradphi_type::reference_element_type fe_modgradphi_type;
    //------------------------------------------------------------------------------//
    // expression desc
    typedef typename functionspace_normal_type::geoelement_type geoelement_type;
    typedef typename functionspace_normal_type::value_type value_type;
    static const uint16_type nDim = fe_normal_type::nDim;
    static const uint16_type nRealDim = fe_normal_type::nRealDim;
    static const uint16_type rank = fe_normal_type::rank;
    static const uint16_type nComponents1 = fe_normal_type::nComponents1;
    static const uint16_type nComponents2 = fe_normal_type::nComponents2;
    static const bool is_terminal = false;

    static const uint16_type ordernormal = functionspace_normal_type::basis_type::nOrder;
    static const uint16_type ordercurvature = functionspace_curvature_type::basis_type::nOrder;

    // imorder
    static const uint16_type imorder_K2N = 2*ordercurvature + ordernormal;
    static const uint16_type imorder_gradK = ordercurvature - 1;
    static const uint16_type imorder_NdotGradKN = 2*ordernormal + imorder_gradK;
    static const uint16_type imorderAuto = mpl::max< 
        mpl::int_<imorder_K2N>, mpl::int_<imorder_NdotGradKN>
        >::type::value;
    static const uint16_type imorder = (IMOrder >= 0) ? IMOrder : imorderAuto;

    static const bool imIsPoly = true;

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

    //------------------------------------------------------------------------------//
    HelfrichInnerDivImpl( element_normal_type const& N, element_curvature_type const& K,
            element_modgradphi_type const& modGradPhi, HelfrichInnerDivImplementation impl )
        :
            M_normal( boost::cref(N) ),
            M_curvature( boost::cref(K) ),
            M_modgradphi( boost::cref(modGradPhi) ),
            M_impl( impl )
    {}
    ~HelfrichInnerDivImpl() = default;

    //------------------------------------------------------------------------------//
    uint16_type polynomialOrder() const { return (IMOrder>=0)?IMOrder:imorderAuto; }
    bool isPolynomial() const { return true; }
    //------------------------------------------------------------------------------//
    element_normal_type const& normal() const { return M_normal; }
    element_curvature_type const& curvature() const { return M_curvature; }
    element_modgradphi_type const& modgradphi() const { return M_modgradphi; }

    HelfrichInnerDivImplementation const& impl() const { return M_impl; }

    //------------------------------------------------------------------------------//
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor 
    {
        typedef this_type expr_type;
        typedef typename element_normal_type::value_type value_type;

        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        // tensor base
        typedef tensorHelfrichInnerDivTensorBase<Geo_t, Basis_i_t, Basis_j_t, expr_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;
        // shape
        typedef typename tensorbase_type::shape shape;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        // is zero
        struct is_zero
        {
            static const bool value = false;
        };

        // constructor
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            switch( expr.impl() )
            {
                case HelfrichInnerDivImplementation::GENERIC:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom, fev, feu )
                            );
                    break;
                case HelfrichInnerDivImplementation::DISTANCE:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom, fev, feu )
                            );
                    break;
            }
        }
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            switch( expr.impl() )
            {
                case HelfrichInnerDivImplementation::GENERIC:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom, fev )
                            );
                    break;
                case HelfrichInnerDivImplementation::DISTANCE:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom, fev )
                            );
                    break;
            }
        }
        tensor( this_type const& expr, Geo_t const& geom )
        {
            switch( expr.impl() )
            {
                case HelfrichInnerDivImplementation::GENERIC:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorGeneric<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom )
                            );
                    break;
                case HelfrichInnerDivImplementation::DISTANCE:
                    M_tensorbase.reset( 
                            new tensorHelfrichInnerDivTensorDistance<Geo_t, Basis_i_t, Basis_j_t, expr_type>( expr, geom )
                            );
                    break;
            }
        }
        template<typename IM>
        void init( IM const& im )
        {}
        // update
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
    boost::reference_wrapper<const element_normal_type> M_normal;
    boost::reference_wrapper<const element_curvature_type> M_curvature;
    boost::reference_wrapper<const element_modgradphi_type> M_modgradphi;
    HelfrichInnerDivImplementation M_impl;
};

template<int IMOrder = -1, typename ElementNormalType, typename ElementCurvatureType, typename ElementModGradPhiType>
inline
Expr< HelfrichInnerDivImpl<ElementNormalType, ElementCurvatureType, ElementModGradPhiType, IMOrder> >
helfrichInnerDivExpr( ElementNormalType const& N, ElementCurvatureType const& K, ElementModGradPhiType const& modGradPhi, HelfrichInnerDivImplementation impl = HelfrichInnerDivImplementation::GENERIC )
{
    typedef HelfrichInnerDivImpl<ElementNormalType, ElementCurvatureType, ElementModGradPhiType, IMOrder> helfrich_inner_div_type;
    return Expr<helfrich_inner_div_type>( helfrich_inner_div_type(N, K, modGradPhi, impl) );
}

} // namespace FeelModels
} // namespace Feel

#endif // _HELFRICH_FORCE_EXPR_HPP
