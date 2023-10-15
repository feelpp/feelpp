#ifndef _GLOBAL_LEVELSET_EXPR_HPP
#define _GLOBAL_LEVELSET_EXPR_HPP 1

namespace Feel {
namespace FeelModels {

template<typename ElementLevelsetType>
class GlobalLevelsetImpl
{
public:
    typedef GlobalLevelsetImpl<ElementLevelsetType> this_type;

    //--------------------------------------------------------------------//
    typedef ElementLevelsetType element_levelset_type;
    typedef std::shared_ptr<element_levelset_type> element_levelset_ptrtype;
    typedef typename element_levelset_type::functionspace_type functionspace_levelset_type;
    typedef typename functionspace_levelset_type::reference_element_type fe_levelset_type;

    typedef typename functionspace_levelset_type::value_type value_type;
    typedef value_type evaluate_type;
    // dim
    static inline const uint16_type nDim = fe_levelset_type::nDim;
    static inline const uint16_type nRealDim = fe_levelset_type::nRealDim;

    //--------------------------------------------------------------------//
    static const size_type context = vm::JACOBIAN;
    static const size_type context_levelset = context;
    static const bool is_terminal = true;
    static inline const uint16_type orderlevelset = functionspace_levelset_type::basis_type::nOrder;
    static inline const uint16_type imorder = orderlevelset;
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

    //--------------------------------------------------------------------//
    explicit GlobalLevelsetImpl( std::vector<element_levelset_ptrtype> const& levelsets )
        : M_levelsets( levelsets )
    {}

    ~GlobalLevelsetImpl() = default;

    //--------------------------------------------------------------------//
    //! polynomial order
    uint16_type polynomialOrder() const { return imorder; }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    //--------------------------------------------------------------------//
    bool useInterpWithConfLoc() const { return false; }

    //--------------------------------------------------------------------//
    std::vector<element_levelset_ptrtype> const& levelsets() const { return M_levelsets; }
    element_levelset_ptrtype const& levelsets( size_type i ) const { return M_levelsets[i]; }

    //--------------------------------------------------------------------//
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expr_type;

        typedef typename expr_type::element_levelset_type::value_type value_type;
        typedef typename expr_type::functionspace_levelset_type::reference_element_type fe_levelset_type;
        typedef typename expr_type::functionspace_levelset_type::geoelement_type geoelement_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;
        // fe levelset context
        typedef typename fe_levelset_type::PreCompute pc_levelset_type;
        typedef std::shared_ptr<pc_levelset_type> pc_levelset_ptrtype;
        typedef typename fe_levelset_type::template Context<expr_type::context_levelset, fe_levelset_type, gm_type, geoelement_type, gmc_type::context> ctx_levelset_type;
        typedef std::shared_ptr<ctx_levelset_type> ctx_levelset_ptrtype;

        typedef typename matrix_node<value_type>::type matrix_node_type;

        static const bool isSameGeo = boost::is_same<typename gmc_type::element_type, geoelement_type>::value;

        // shape
        typedef Shape<expr_type::nRealDim, Scalar, false, false> shape_type;
        typedef Eigen::Matrix<value_type, shape_type::M, shape_type::N> matrix_shape_type;
        typedef boost::multi_array<matrix_shape_type, 1> array_shape_type;
        typedef shape_type shape;
        // scalar array
        typedef Shape<shape_type::nDim, Scalar, false, false> shape_scalar;
        typedef Eigen::TensorFixedSize<value_type, Eigen::Sizes<shape_scalar::M,shape_scalar::N>> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type, 1> array_scalar_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : 
                M_expr( expr ),
                M_pcLevelset( new pc_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom)->xRefs() ) ),
                M_ctxLevelset( new ctx_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom), M_pcLevelset ) ),
                M_hasRelationMesh( expr.levelsets().size() ),
                M_isSameMesh( expr.levelsets().size() ),
                M_zeroLocScalar(),
                M_locLevelsets( expr.levelsets().size(), array_scalar_type(expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom))) ),
                M_locRes( expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom)) )
        {
            for( size_type i = 0; i < expr.levelsets().size(); ++i )
            {
                M_hasRelationMesh[i] = fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo( expr.levelsets()[i]->functionSpace()->mesh() );
                M_isSameMesh[i] = M_hasRelationMesh[i] && isSameGeo;
                if( !M_isSameMesh[i] )
                    expr.levelsets()[i]->functionSpace()->mesh()->tool_localization()->updateForUse();
            }

            M_zeroLocScalar.setZero();
        }
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : 
                M_expr( expr ),
                M_pcLevelset( new pc_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom)->xRefs() ) ),
                M_ctxLevelset( new ctx_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom), M_pcLevelset ) ),
                M_hasRelationMesh( expr.levelsets().size() ),
                M_isSameMesh( expr.levelsets().size() ),
                M_zeroLocScalar(),
                M_locLevelsets( expr.levelsets().size(), array_scalar_type(expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom))) ),
                M_locRes( expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom)) )
        {
            for( size_type i = 0; i < expr.levelsets().size(); ++i )
            {
                M_hasRelationMesh[i] = fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo( expr.levelsets()[i]->functionSpace()->mesh() );
                M_isSameMesh[i] = M_hasRelationMesh[i] && isSameGeo;
                if( !M_isSameMesh[i] )
                    expr.levelsets()[i]->functionSpace()->mesh()->tool_localization()->updateForUse();
            }

            M_zeroLocScalar.setZero();
        }
        tensor( this_type const& expr, Geo_t const& geom )
            : 
                M_expr( expr ),
                M_pcLevelset( new pc_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom)->xRefs() ) ),
                M_ctxLevelset( new ctx_levelset_type( expr.levelsets(0)->functionSpace()->fe(), fusion::at_key<key_type>(geom), M_pcLevelset ) ),
                M_hasRelationMesh( expr.levelsets().size() ),
                M_isSameMesh( expr.levelsets().size() ),
                M_zeroLocScalar(),
                M_locLevelsets( expr.levelsets().size(), array_scalar_type(expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom))) ),
                M_locRes( expr.levelsets(0)->idExtents(*fusion::at_key<key_type>(geom)) )
        {
            for( size_type i = 0; i < expr.levelsets().size(); ++i )
            {
                M_hasRelationMesh[i] = fusion::at_key<key_type>(geom)->element().mesh()->isRelatedTo( expr.levelsets()[i]->functionSpace()->mesh() );
                M_isSameMesh[i] = M_hasRelationMesh[i] && isSameGeo;
                if( !M_isSameMesh[i] )
                    expr.levelsets()[i]->functionSpace()->mesh()->tool_localization()->updateForUse();
            }

            M_zeroLocScalar.setZero();
        }

        expr_type const& expr() const { return M_expr; }

        template<typename IM>
        void init( IM const& im )
        {}
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
            // Get levelsets values
            this->updateCtxIfSameGeo( geom, mpl::bool_<isSameGeo>() );
            for( size_type i = 0; i < M_locLevelsets.size(); ++i )
            {
                auto & locLevelset = M_locLevelsets[i];
                std::fill( locLevelset.data(), locLevelset.data()+locLevelset.num_elements(), this->M_zeroLocScalar );
                if( M_isSameMesh[i] )
                {
                    this->expr().levelsets(i)->id( *M_ctxLevelset, locLevelset );
                }
                else
                {
                    matrix_node_type ptsreal = this->expr().levelsets(i)->ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<1>());
                    auto setOfPts = ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) ?
                        fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() :
                        fusion::at_key<key_type>( geom )->element().vertices();
                    this->expr().levelsets(i)->idInterpolate( ptsreal, locLevelset, this->expr().useInterpWithConfLoc(), setOfPts );
                }
            }
            // Compute results: min( phi_i )
            this->updateGlobalLevelsetImpl( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            // Get levelsets values
            this->updateCtxFaceIfSameGeo( geom, mpl::bool_<isSameGeo>() );
            for( size_type i = 0; i < M_locLevelsets.size(); ++i )
            {
                auto & locLevelset = M_locLevelsets[i];
                std::fill( locLevelset.data(), locLevelset.data()+locLevelset.num_elements(), this->M_zeroLocScalar );
                if( M_isSameMesh[i] )
                {
                    this->expr().levelsets(i)->id( *M_ctxLevelset, locLevelset );
                }
                else
                {
                    matrix_node_type ptsreal = this->expr().levelsets(i)->ptsInContext(*fusion::at_key<key_type>(geom), mpl::int_<2>());
                    this->expr().levelsets(i)->idInterpolate( ptsreal, locLevelset, this->expr().useInterpWithConfLoc(), fusion::at_key<key_type>(geom)->element().face( fusion::at_key<key_type>(geom)->faceId() ).vertices() );
                }
            }
            // Compute results: min( phi_i )
            this->updateGlobalLevelsetImpl( geom );
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
            return M_locRes[q](0,0);
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_locRes[q];
        }

    private:
        void updateCtxIfSameGeo(Geo_t const& geom, mpl::bool_<true> )
        {   
            if (fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) /*face case*/
                M_pcLevelset->update(fusion::at_key<key_type>( geom )->pc()->nodes() );
            M_ctxLevelset->update( fusion::at_key<key_type>( geom ),  M_pcLevelset );
        }
        void updateCtxIfSameGeo(Geo_t const& geom, mpl::bool_<false> )
        {
        }

        void updateCtxFaceIfSameGeo(Geo_t const& geom, mpl::bool_<true> )
        {
            M_pcLevelset->update(fusion::at_key<key_type>( geom )->pc()->nodes() );
            M_ctxLevelset->update( fusion::at_key<key_type>( geom ), M_pcLevelset );
        }
        void updateCtxFaceIfSameGeo(Geo_t const& geom, mpl::bool_<false> )
        {
        }

        void updateGlobalLevelsetImpl( Geo_t const& geom )
        {
            // Compute results: min( phi_i )
            auto const& gmc = fusion::at_key<key_type>( geom );
            std::fill( M_locRes.data(), M_locRes.data()+M_locRes.num_elements(), matrix_shape_type::Zero() );
            for( uint16_type q = 0; q < gmc->nPoints(); ++q )
            {
                M_locRes[q](0,0) = M_locLevelsets[0][q](0,0);
                for( size_type i = 0; i < M_locLevelsets.size(); ++i )
                {
                    if( M_locLevelsets[i][q](0,0) < M_locRes[q](0,0) )
                        M_locRes[q](0,0) = M_locLevelsets[i][q](0,0);
                }
            }
        }
        
    private:
        expr_type const& M_expr;

        pc_levelset_ptrtype M_pcLevelset;
        ctx_levelset_ptrtype M_ctxLevelset;

        std::vector<bool> M_hasRelationMesh;
        std::vector<bool> M_isSameMesh;

        loc_scalar_type M_zeroLocScalar;
        std::vector<array_scalar_type> M_locLevelsets;

        array_shape_type M_locRes;
    };

private:
    std::vector<element_levelset_ptrtype> M_levelsets;
};

template<typename ElementLevelsetType>
inline Expr< GlobalLevelsetImpl<ElementLevelsetType> >
globalLevelsetExpr( std::vector< std::shared_ptr<ElementLevelsetType> > const& levelsets )
{
    typedef GlobalLevelsetImpl<ElementLevelsetType> expr_globallevelset_type;
    return Expr<expr_globallevelset_type>( expr_globallevelset_type(levelsets) );
}

} // namespace FeelModels
} // namespace Feel

#endif //_GLOBAL_LEVELSET_EXPR_HPP

