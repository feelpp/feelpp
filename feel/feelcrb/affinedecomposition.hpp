#ifndef __AFFINE_DECOMPOSITION_H
#define __AFFINE_DECOMPOSITION_H 1

#include <feel/feel.hpp>
#include <boost/shared_ptr.hpp>


namespace Feel
{



template <typename ModelType, typename TensorType, typename TestSpaceType, typename TrialSpaceType=TestSpaceType >
class BlockAD
{
public:
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename model_type::value_type value_type;
    typedef typename model_type::mesh_type mesh_type;
    //typedef typename model_type::mesh_ptrtype mesh_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename model_type::parameter_type parameter_type;

    static const bool use_block_structure = model_type::use_block_structure;
    static const bool is_matrix = std::is_same<TensorType,sparse_matrix_ptrtype>::value;
    static const int geo_dim = mesh_type::nDim;

    typedef typename model_type::space_type fullspace_type;
    typedef typename model_type::element_type element_type;

    typedef TestSpaceType testspace_type;
    typedef TrialSpaceType trialspace_type;

    typedef boost::shared_ptr<testspace_type> testspace_ptrtype;
    typedef boost::shared_ptr<trialspace_type> trialspace_ptrtype;

    typedef typename testspace_type::element_type testelement_type;
    typedef typename trialspace_type::element_type trialelement_type;

    typedef vf::detail::LinearForm<testspace_type,vector_type,vector_type> form1_type;
    typedef vf::detail::BilinearForm<testspace_type, trialspace_type, VectorUblas<value_type>> form2_type;

    typedef OperatorLinear< testspace_type , trialspace_type > operator_type;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;
    typedef FsFunctionalLinear< testspace_type > functional_type;
    typedef boost::shared_ptr<functional_type> functional_ptrtype;

    typedef OperatorLinearComposite< testspace_type , trialspace_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;
    typedef FsFunctionalLinearComposite< testspace_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef typename mpl::if_< mpl::bool_<is_matrix>,
                               operatorcomposite_ptrtype,
                               functionalcomposite_ptrtype >::type composite_ptrtype;


    typedef std::vector< std::vector< sparse_matrix_ptrtype >> affine_matrix_type;
    typedef std::vector< std::vector< affine_matrix_type >> block_affine_matrix_type;
    typedef std::vector< std::vector < vector_ptrtype >> affine_vector_type;
    typedef std::vector< affine_vector_type > affine_output_type;
    typedef std::vector< affine_output_type > block_affine_output_type;

    typedef std::vector< std::vector< value_type >> beta_vector_type;
    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > betaqm_type;

    typedef typename model_type::offline_merge_type offline_merge_type;
    typedef typename model_type::index_vector_type index_vector_type;

    typedef boost::tuple< affine_matrix_type, // Mqm
                          affine_matrix_type, // Aqm
                          affine_output_type > // Fqm
    affine_decomposition_type;


    typedef EIMFunctionBase< testspace_type, trialspace_type, parameterspace_type > fun_type;
    typedef boost::shared_ptr<fun_type> fun_ptrtype;


    BlockAD( model_ptrtype model, testspace_ptrtype testspace ) :
        M_model(model),
        M_Q( 0 ),
        M_use_operators_free( false ),
        M_use_ginac_expr( false ),
        Xtest( testspace ),
        Xtrial( testspace )
        {}

    BlockAD( model_ptrtype model, testspace_ptrtype testspace, trialspace_ptrtype trialspace ) :
        M_model(model),
        M_Q( 0 ),
        M_use_operators_free( false ),
        M_use_ginac_expr( false ),
        Xtest( testspace ),
        Xtrial( trialspace )
        {}


        template <typename BetaType>
    void addOpe( form2_type const& form2, BetaType const& beta, std::string filename )
        {
            add( form2.matrixPtr(), beta, filename );
        }

    template <typename ExprType, typename BetaType>
    void addOpe( ExprType const& expr, BetaType const& beta , std::string filename )
        {
            M_use_operators_free = true;
            auto ope = opLinearFree( _domainSpace=Xtrial,
                                     _imageSpace=Xtest,
                                     _expr=expr );

            if ( !M_composite )
                M_composite = opLinearComposite( _imageSpace=Xtest,
                                                 _domainSpace=Xtrial );

            newEntry( beta, filename );
            filename += std::to_string( M_Q );
            ope->setName( filename );
            M_composite->addElement( {M_Q-1,0}, ope );
        }

    template <typename BetaType>
    void addOpe( sparse_matrix_ptrtype mat, BetaType const& beta, std::string filename )
        {
            add( mat, beta, filename );
        }

    template <typename BetaType>
    void addFun( vector_ptrtype const& vec, BetaType const& beta, std::string filename )
        {
            add( vec, beta, filename );
        }

    template <typename BetaType>
    void addFun( form1_type const& form1, BetaType const& beta, std::string filename )
        {
            add( form1.vectorPtr(), beta, filename );
        }
    template <typename ExprType, typename BetaType>
    void addFun( ExprType const& expr, BetaType const& beta, std::string filename )
        {
            auto fun = functionalLinearFree( _space=Xtest,
                                             _expr=expr );

            M_use_operators_free = true;
            if ( !M_composite )
                M_composite = functionalLinearComposite( _space=Xtest );

            newEntry( beta, filename );
            filename += std::to_string( M_Q );
            fun->setName( filename );
            M_composite->addElement( {M_Q-1,0}, fun );
        }

    sparse_matrix_ptrtype ptrMat( int const& q, int const& m )
        {
            if ( q>=M_Q )
                resizeQ( q+1 );
            if ( m>=M_mMax[q] )
                resizeM( q, m+1 );
            if ( !M_Mqm[q][m] )
                M_Mqm[q][m] = backend()->newMatrix( _test=Xtest,
                                                    _trial=Xtrial );
            return M_Mqm[q][m];
        }
    vector_ptrtype ptrVec( int const& q, int const& m )
        {
            if ( q>=M_Q )
                resizeQ( q+1 );
            if ( m>=M_mMax[q] )
                resizeM( q, m+1 );
            if ( !M_Mqm[q][m] )
                M_Mqm[q][m] = backend()->newVector( Xtest );
            return M_Mqm[q][m];
        }

    void check()
        {
            CHECK( M_Q==M_mMax.size() )<<"M_Q="<<M_Q<<", M_mMax.size()="<<M_mMax.size()<<std::endl;
            CHECK( M_Q==M_beta.size() )<<"M_Q="<<M_Q<<", M_beta.size()="<<M_beta.size()<<std::endl;
            CHECK( M_Q==M_beta_mode.size() )<<"M_Q="<<M_Q
                                            <<", M_beta_mode.size()="<<M_beta_mode.size()<<std::endl;
            for ( int q=0; q<M_Q; q++ )
                CHECK( M_mMax[q]==M_beta[q].size() )<<"M_mMax[q]="<<M_mMax[q]
                                                    <<", M_beta[q].size()="<<M_beta[q].size()
                                                    <<", for q="<<q<<std::endl;
            if ( M_use_operators_free )
            {
                auto mMax = M_composite->countAllContributions();
                CHECK( M_Q==mMax.size() )<< "M_Q==mMax.size()\n";
                for ( int q=0; q<M_Q; q++ )
                    CHECK( mMax[q]==M_mMax[q] )<<"mMax[q]="<<mMax[q]
                                               <<", M_mMax[q]="<<M_mMax[q]
                                               <<", for q="<<q<<std::endl;
            }
            else
            {
                CHECK( M_Q==M_Mqm.size() )<<"M_Q="<<M_Q<<", M_Mqm.size()="<<M_Mqm.size()<<std::endl;
                for ( int q=0; q<M_Q; q++ )
                    CHECK( M_mMax[q]==M_Mqm[q].size() )<<"M_mMax[q]="<<M_mMax[q]
                                                       <<", M_Mqm[q].size()="<<M_Mqm[q].size()
                                                       <<", for q="<<q<<std::endl;
            }

        }
    void initializeMassMatrix()
        {
            initializeMassMatrix( mpl::bool_<use_block_structure>() );
        }
    void initializeMassMatrix( mpl::bool_<false> )
        {
            std::string filename = "GinacM-";

            if ( M_model->constructOperatorCompositeM() )
            {
                M_use_operators_free=true;
                M_composite = M_model->operatorCompositeM();
                M_mMax = M_composite->countAllContributions();
                M_Q = M_mMax.size();
            }
            else if ( boption("crb.stock-matrices") )
                assembleMassMatrix();
            else
                preAssembleMassMatrix();

            M_beta.resize( M_Q );
            M_beta_mode.resize( M_Q, BetaMode::SET );
            for( int q=0; q<M_Q; q++ )
            {
                M_beta[q].resize( M_mMax[q] );
                for( int m=0; m<M_mMax[q]; m++ )
                    M_beta[q][m]=1.;
            }
        }
    void initializeMassMatrix( mpl::bool_<true> )
        {
            std::string filename = "GinacM-";
            CHECK( !M_model->constructOperatorCompositeM() )<<"You can't construct the operator compositeM with block structure, set it to false and use default mass matrix\n";
            if ( boption("crb.stock-matrices") )
                assembleMassMatrix();
            else
                preAssembleMassMatrix();

            M_beta.resize( M_Q );
            M_beta_mode.resize( M_Q, BetaMode::SET );
            for( int q=0; q<M_Q; q++ )
            {
                M_beta[q].resize( M_mMax[q] );
                for( int m=0; m<M_mMax[q]; m++ )
                    M_beta[q][m]=1.;
            }
        }


    beta_vector_type computeBetaQm( parameter_type const& mu, double time )
        {
            if ( M_use_ginac_expr )
            {
                std::string symbol;
                for( int i=0; i<M_model->ParameterSpaceDimension; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;
            }

            for ( int q=0; q<M_Q; q++ )
            {
                switch ( M_beta_mode[q] )
                {
                case BetaMode::GINAC :
                    M_beta[q][0] = M_ginac[q].evaluate(M_symbols_map);
                    break;

                case BetaMode::SET :
                    break;

                case BetaMode::EIM :
                    for ( int m=0; m<M_mMax[q]; m++ )
                        M_beta[q][m] = M_eims[q]->beta( mu )(m);
                    break;

                default :
                    CHECK( false )<< "You did not fill betaMqm coefficient\n, you can use either ginac expressions, eim or betaMqm() function to do so"<<std::endl;
                    break;
                }
            }
            return M_beta;
        }
    beta_vector_type computeBetaQm( element_type const& T, parameter_type const& mu, double time )
        {
            if ( M_use_ginac_expr )
            {
                std::string symbol;
                for( int i=0; i<M_model->ParameterSpaceDimension; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;
            }

            for ( int q=0; q<M_Q; q++ )
            {
                switch ( M_beta_mode[q] )
                {
                case BetaMode::GINAC :
                    M_beta[q][0] = M_ginac[q].evaluate(M_symbols_map);
                    break;

                case BetaMode::SET :
                    break;

                case BetaMode::EIM :
                    for ( int m=0; m<M_mMax[q]; m++ )
                        M_beta[q][m] = M_eims[q]->beta( mu, T )(m);
                    break;

                default :
                    CHECK( false )<< "You did not fill betaMqm coefficient\n, you can use either ginac expressions, eim or betaMqm() function to do so"<<std::endl;
                    break;
                }
            }
            return M_beta;
        }

    beta_vector_type betaQm()
        {
            return M_beta;
        }

    double& setBeta( int const& q, int const& m )
        {
            if ( q>=M_Q )
                resizeQ( q+1 );

            if ( m>=M_mMax[q] )
                resizeM( q, m+1 );

            M_beta_mode[q] = BetaMode::SET;

            return M_beta[q][m];
        }

    double beta( int const& q, int const& m )
        {
            CHECK( q<M_Q ) << "Bad q index "<<q<<", max index is " << M_Q<<std::endl;
            CHECK( m<M_mMax[q] ) << "Bad m index "<<q<<", max index is " << M_mMax[q]<<std::endl;
            return M_beta[q][m];
        }

    std::vector< std::vector< TensorType >> compute()
        {
            return compute( mpl::bool_<is_matrix>() );
        }
    std::vector< std::vector< TensorType >> compute( mpl::bool_<true> )
        {
            if ( M_use_operators_free && M_Mqm.size()==0 )
            {
                M_Mqm.resize(M_Q);
                for ( int q=0; q<M_Q; q++ )
                {
                    int mMax = M_mMax[q];
                    M_Mqm[q].resize( mMax );
                    for( int m=0; m<mMax; m++ )
                    {
                        auto ope = M_composite->operatorlinear(q,m);
                        M_Mqm[q][0] = backend()->newMatrix( _trial=ope->domainSpace(),
                                                            _test=ope->dualImageSpace(),
                                                            _pattern=ope->pattern() );
                        ope->matPtr( M_Mqm[q][m] );
                    } //m
                } //q
            }  // if matrix not filled

            return M_Mqm;
        }
    std::vector< std::vector< TensorType >> compute( mpl::bool_<false> )
        {
            if ( M_use_operators_free && M_Mqm.size()==0 )
            {
                M_Mqm.resize(M_Q);
                for ( int q=0; q<M_Q; q++ )
                {
                    int mMax = M_mMax[q];
                    M_Mqm[q].resize( mMax );
                    for( int m=0; m<mMax; m++ )
                    {
                        auto fun = M_composite->functionallinear(q,m);
                        M_Mqm[q][0] = backend()->newVector( fun->space() );
                        fun->containerPtr( M_Mqm[q][m] );
                    } //m
                } //q
            }  // if matrix not filled

            return M_Mqm;
        }

    int Q()
        {
            return M_Q;
        }
    int mMax( int q )
        {
            int size = M_mMax.size();
            CHECK( q<size )<<"[mMax(q)] function called with bad q index : "<< q
                           << " and max possible index is "<< size << std::endl;

            return M_mMax[q];
        }

    sparse_matrix_ptrtype compute( int q, int m, bool transpose )
        {
            CHECK( q<M_Q ) << "compute called with bad q index : "<<q <<" and max index is "<< M_Q << std::endl;
            CHECK( m<M_mMax[q] ) << "compute called with bad m index : "<< m <<" and max index is "<< M_mMax[q];

            sparse_matrix_ptrtype matrix;
            if ( M_use_operators_free )
            {
                auto ope = M_composite->operatorlinear(q,m);
                matrix = backend()->newMatrix( _trial=ope->domainSpace(),
                                               _test=ope->dualImageSpace(),
                                               _pattern=ope->pattern() );
                ope->matPtr( matrix );
                if ( transpose )
                {
                    auto mt = backend()->newMatrix( _trial=ope->domainSpace(),
                                                           _test=ope->dualImageSpace(),
                                                           _pattern=ope->pattern() );
                    matrix->transpose( mt );
                    return mt;
                }
            }
            else
            {
                matrix  = backend()->newMatrix( _trial=Xtrial,
                                                _test=Xtest );
                if ( transpose )
                    M_Mqm[q][m]->transpose(matrix);
                else
                    matrix = M_Mqm[q][m];
            }
            return matrix;
        }
    vector_ptrtype compute( int q, int m )
        {
            CHECK( q<M_Q ) << "compute called with bad q index : "<<q <<" and max index is "<< M_Q << std::endl;
            CHECK( m<M_mMax[q] ) << "compute called with bad m index : "<< m <<" and max index is "<< M_mMax[q];

            vector_ptrtype vector;

            if ( M_use_operators_free )
            {
                auto fun = M_composite->functionallinear(q,m);
                vector = backend()->newVector( fun->space() );
                fun->containerPtr( vector );
            }
            else
                vector = M_Mqm[q][m];
            return vector;
        }

    TensorType merge( beta_vector_type const& beta )
        {
            return merge( beta, mpl::bool_<is_matrix>() );
        }
    TensorType merge( beta_vector_type const& beta, mpl::bool_<true> )
        {
            auto M = backend()->newMatrix( _test=Xtest, _trial=Xtrial );
            if ( M_use_operators_free )
            {
                M_composite->setScalars( beta );
                M_composite->sumAllMatrices( M );
            }
            else
            {
                M->zero();
                for ( int q=0; q<M_Q; q++ )
                    for( int m=0; m<M_mMax[q]; m++)
                        M->addMatrix( beta[q][m], M_Mqm[q][m] );
            }
            return M;
        }
    TensorType merge( beta_vector_type const& beta, mpl::bool_<false> )
        {
            auto F = backend()->newVector( Xtest );
            if( M_use_operators_free )
            {
                M_composite->setScalars( beta );
                M_composite->sumAllVectors( F );
            }
            else
            {
                F->zero();
                for( int q=0; q<M_Q; q++ )
                    for( int m=0; m<M_mMax[q]; m++ )
                        F->add( beta[q][m], M_Mqm[q][m] );
            }
            return F;
        }



    sparse_matrix_ptrtype computeLinear( int const& q, int const& m, bool transpose )
        {
            CHECK( q<M_Q ) << "compute called with bad q index : "<<q <<" and max index is "<< M_Q << std::endl;
            CHECK( m<M_mMax[q] ) << "compute called with bad m index : "<< m <<" and max index is "<< M_mMax[q];

            auto matrix  = backend()->newMatrix( _trial=Xtrial,
                                                 _test=Xtest );
            if ( transpose )
                M_Mqm[q][m]->transpose(matrix);
            else
                matrix = M_Mqm[q][m];

            return matrix;
        }

    std::vector< std::vector< sparse_matrix_ptrtype >> linearMqm()
        {
            if ( M_MqmLinear.size()==0 )
            {
                M_MqmLinear = M_model->computeLinearDecompositionA();
                if ( M_MqmLinear.size()==0 )
                {
                    for ( int q=0; q<M_Q; q++ )
                    {
                        if ( M_beta_mode[q]!=BetaMode::EIM )
                            M_MqmLinear.push_back( M_Mqm[q] );
                    }
                }

                M_QLinear = M_MqmLinear.size();
                M_mMaxLinear.resize( M_QLinear );
                M_betaLinear.resize( M_QLinear );
                for ( int q=0; q<M_QLinear; q++ )
                {
                    M_mMaxLinear[q] = M_MqmLinear[q].size();
                    M_betaLinear[q].resize( M_mMaxLinear[q], 0 );
                }
            }

            return M_MqmLinear;
        }
    beta_vector_type computeBetaLinear( parameter_type const& mu, double time=0 )
        {
            if ( M_betaLinear.size()==0 )
                linearMqm();

            if ( M_use_ginac_expr )
            {
                std::string symbol;
                for( int i=0; i<M_model->ParameterSpaceDimension; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;
            }

            for ( int q=0; q<M_Q; q++ )
            {
                switch ( M_beta_mode[q] )
                {
                case BetaMode::GINAC :
                    M_betaLinear[q][0] = M_ginac[q].evaluate(M_symbols_map);
                    break;

                case BetaMode::SET :
                    for ( int m=0; q<M_mMax[q]; m++)
                        M_betaLinear[q][m]=M_beta[q][m];
                    break;

                default :
                    CHECK( false )<< "You did not fill betaMqm coefficient\n, you can use either ginac expressions, eim or betaMqm() function to do so"<<std::endl;
                    break;
                }
            }

            return M_betaLinear;
        }

    int QLinear()
        {
            return M_QLinear;
        }
    int mMaxLinear( int const q )
        {
            return M_mMaxLinear[q];
        }

private:
    enum BetaMode { NONE=0, GINAC=1, EIM=2, SET=3 };

    void add( TensorType const& tensor, std::string symbol, std::string filename )
        {
            newEntry( symbol, filename );
            M_Mqm[M_Q-1][0] = tensor;
        }

    void add( sparse_matrix_ptrtype const& mat, fun_ptrtype eim, std::string filename )
        {
            int mMax = newEntry( eim );
            auto eimMat = backend()->newMatrix( _trial=Xtrial,
                                                _test=Xtrial );
            auto mesh = Xtrial->mesh();
            auto u = Xtrial->element();
            auto v = Xtrial->element();
            for ( int m=0; m<mMax; m++ )
            {
                M_Mqm[M_Q-1][m] = backend()->newMatrix( _trial=Xtrial,
                                                        _test=Xtest );
                form2( _trial=Xtrial, _test=Xtrial, _matrix=eimMat )
                    = integrate( elements(mesh), idt(u)*id(v)*idv( eim->q(m)) );
                mat->matMatMult( *eimMat, *M_Mqm[M_Q-1][m] );
            }
        }

    void add( vector_ptrtype const& vec, fun_ptrtype eim, std::string filename )
        {
            int mMax = newEntry( eim );
            auto eimMat = backend()->newMatrix( _trial=Xtest,
                                                _test=Xtest );
            auto mesh = Xtest->mesh();
            auto u = Xtest->element();
            auto v = Xtest->element();
            for ( int m=0; m<mMax; m++ )
            {
                M_Mqm[M_Q-1][m] = backend()->newVector( Xtest );
                form2( _trial=Xtest, _test=Xtest, _matrix=eimMat )
                    = integrate( elements(mesh), idt(u)*id(v)*idv( eim->q(m)) );
                M_Mqm[M_Q-1][m]->zero();
                M_Mqm[M_Q-1][m]->addVector( *vec, * eimMat );
            }
        }
    void resizeQ( int const& q )
        {
            M_Q = q;
            M_beta.resize( M_Q );
            M_beta_mode.resize( M_Q, BetaMode::NONE );
            M_mMax.resize( M_Q, 0 );
            M_ginac.resize( M_Q, expr<geo_dim>( "0", Symbols( {"x"} ), "ginacDefault" ) );
            if ( !M_use_operators_free )
                M_Mqm.resize( M_Q );
        }
    void resizeM( int const& q, int const& m )
        {
            M_mMax[q] = m;
            M_beta[q].resize( M_mMax[q], 0 );
            if ( !M_use_operators_free )
                M_Mqm[q].resize( M_mMax[q] );
        }

    void newEntry( std::string const& symbol, std::string filename )
        {
            if ( !M_use_ginac_expr )
            {
                M_use_ginac_expr=true;
                buildSymbolsVector();
                buildSymbolsMap();
            }

            filename = "Ginac" + filename + std::to_string( M_Q );

            resizeQ( M_Q+1 );
            resizeM( M_Q-1, 1 );

            M_ginac[M_Q-1]= expr<geo_dim>( symbol, Symbols( M_symbols_vec ), filename );
            M_beta_mode[M_Q-1]= BetaMode::GINAC;
        }
    int newEntry( fun_ptrtype eim )
        {
            M_model->setHasEim( true );
            resizeQ( M_Q+1 );
            bool error;
            int mMax = eim->mMax(error);
            M_mMax[M_Q-1]=mMax;
            if ( error )
                mMax++;

            M_beta[M_Q-1].resize( M_mMax[M_Q-1], 0 );

            if ( !M_use_operators_free )
                M_Mqm[M_Q-1].resize( mMax );
            else
                CHECK( false )<<"You cannot use eim with operators free, for now\n";

            M_eims.resize( M_Q );
            M_eims[M_Q-1]=eim;
            M_beta_mode[M_Q-1] = BetaMode::EIM;

            return mMax;
        }


    void assembleMassMatrix()
        {
            const bool is_composite = trialspace_type::is_composite || testspace_type::is_composite;
            assembleMassMatrix( mpl::bool_<is_composite>() );
        }
    void assembleMassMatrix( mpl::bool_<false> )
        {
            auto mesh = Xtrial->mesh();
            auto u = Xtrial->element();
            auto v = Xtest->element();

            M_Q = (int) 1 ;
            M_mMax.resize(M_Q);
            M_mMax[0] = 1;

            M_Mqm.resize( M_Q );
            M_Mqm[0].resize( 1 );
            M_Mqm[0][0] = backend()->newMatrix( _test=Xtest , _trial=Xtrial );

            form2( _test=Xtest, _trial=Xtrial, _matrix=M_Mqm[0][0] ) =
                integrate( _range=elements( mesh ), _expr= inner( idt(u), id(v) )  ) ;
            M_Mqm[0][0]->close();
        }
    void assembleMassMatrix( mpl::bool_<true> )
        {
            CHECK( !use_block_structure )<<"You cannot use composite spaces and block structure for now";
            index_vector_type index_vector;
            int n_spaces = fullspace_type::nSpaces;
            auto Xh = M_model->functionSpace();
            auto u = Xh->element();
            auto v = Xh->element();

            M_Q = n_spaces;
            M_Mqm.resize( M_Q );
            for ( int q=0; q<M_Q; q++ )
            {
                M_mMax.push_back(1);
                M_Mqm[q].resize(1);
                M_Mqm[q][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
            }
            AssembleMassMatrixInCompositeCase assemble_mass_matrix_in_composite_case ( u , v , this );
            fusion::for_each( index_vector, assemble_mass_matrix_in_composite_case );
        }
    struct AssembleMassMatrixInCompositeCase
    {
        AssembleMassMatrixInCompositeCase( element_type const u ,
                                           element_type const v ,
                                           std::vector< std::vector< sparse_matrix_ptrtype >> Mqm)
            :
            M_composite_u ( u ),
            M_composite_v ( v ),
            M_Mqm ( Mqm )
            {}

        template< typename T >
        void operator()( const T& t ) const
            {
                auto u = this->M_composite_u.template element< T::value >();
                auto v = this->M_composite_v.template element< T::value >();
                auto Xh = M_composite_u.functionSpace();
                auto mesh = Xh->mesh();

                int q = T::value;
                form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[q][0] ) +=
                    integrate( _range=elements( mesh ), _expr=trans( idt( u ) )*id( v ) );
                M_Mqm[q][0]->close();
            }

        element_type  M_composite_u;
        element_type  M_composite_v;
        std::vector< std::vector< sparse_matrix_ptrtype >> M_Mqm;
    };

    void preAssembleMassMatrix()
        {
            M_use_operators_free = true;
            const bool is_composite = trialspace_type::is_composite || testspace_type::is_composite;
            preAssembleMassMatrix( mpl::bool_<is_composite>() );

            M_mMax = M_composite->countAllContributions();
            M_Q = M_mMax.size();
        }
    void preAssembleMassMatrix( mpl::bool_<false> )
        {
            auto u = Xtrial->element();
            auto v = Xtest->element();

            auto mesh = Xtest->mesh();
            auto expr=integrate( _range=elements( mesh ) , _expr= trans(idt( u ))*id( v ) );

            auto M_composite = opLinearComposite( _domainSpace=Xtrial, _imageSpace=Xtest  );
            auto opfree = opLinearFree( _domainSpace=Xtrial, _imageSpace=Xtest, _expr=expr );

            opfree->setName("mass operator (automatically created)");
            M_composite->addElement( boost::make_tuple(0,0) , opfree );
        }
    void preAssembleMassMatrix( mpl::bool_<true> )
        {
            CHECK( !use_block_structure )<<"You cannot use composite spaces and block structure for now";
            index_vector_type index_vector;
            auto Xh = M_model->functionSpace();
            auto u = Xh->element();
            auto v = Xh->element();

            PreAssembleMassMatrixInCompositeCase preassemble_mass_matrix_in_composite_case ( u , v );
            fusion::for_each( index_vector, preassemble_mass_matrix_in_composite_case );

            M_composite = preassemble_mass_matrix_in_composite_case.opmass();
        }
    struct PreAssembleMassMatrixInCompositeCase
    {
        PreAssembleMassMatrixInCompositeCase( element_type const u ,
                                              element_type const v )
            :
            M_composite_u ( u ),
            M_composite_v ( v )
            {
                auto Xh = M_composite_u.functionSpace();
                op_mass = opLinearComposite( _imageSpace=Xh, _domainSpace=Xh );
            }

        template< typename T >
        void
        operator()( const T& t ) const
            {
                auto u = M_composite_u.template element< T::value >();
                auto v = M_composite_v.template element< T::value >();
                auto Xh = M_composite_u.functionSpace();
                auto mesh = Xh->mesh();

                auto expr = integrate( _range=elements( mesh ) , _expr=trans( idt( u ) )*id( v ) ) ;
                auto opfree = opLinearFree( _imageSpace=Xh, _domainSpace=Xh, _expr=expr );

                //each composant of the affine decomposition
                //is the subspace contribution
                int q=T::value;
                int m=0;
                auto tuple = boost::make_tuple( q , m );
                op_mass->addElement( tuple , opfree );
            }
        operatorcomposite_ptrtype opmass()
            {
                return op_mass;
            }

        element_type  M_composite_u;
        element_type  M_composite_v;
        operatorcomposite_ptrtype op_mass;
    };


    void buildSymbolsVector()
        {
            for( int i=0; i<M_model->ParameterSpaceDimension; i++)
            {
                std::string symbol = ( boost::format("mu%1%") %i ).str();
                M_symbols_vec.push_back( symbol );
            }
            M_symbols_vec.push_back( "t" );
        }
    void buildSymbolsMap()
        {
            for( int i=0; i<M_model->ParameterSpaceDimension; i++)
            {
                std::string symbol = ( boost::format("mu%1%") %i ).str();
                M_symbols_map.insert( std::pair< std::string,double> (symbol,0) );
            }
            M_symbols_map.insert( std::pair< std::string,double> ("t",0) );
        }


    /// beta coefficients
    beta_vector_type M_beta;
    std::vector< BetaMode > M_beta_mode;
    /// ginac expressions
    std::vector< Expr< GinacEx<geo_dim> >> M_ginac;

    model_ptrtype M_model;
    /// number of terms
    int M_Q;
    /// number of eim terms
    std::vector<int> M_mMax;
    /// affine decomposition
    std::vector< std::vector< TensorType >> M_Mqm;
    /// operators free
    composite_ptrtype M_composite;


    /// linear part
    std::vector< std::vector< TensorType >> M_MqmLinear;
    int M_QLinear;
    std::vector< int > M_mMaxLinear;
    beta_vector_type M_betaLinear;

    bool M_use_operators_free;

    /// Ginac tools
    std::map<std::string,double> M_symbols_map;
    std::vector< std::string > M_symbols_vec;

    bool M_use_ginac_expr;
    std::vector< fun_ptrtype > M_eims;

    testspace_ptrtype Xtest;
    trialspace_ptrtype Xtrial;

};


template< typename ModelType >
class AffineDecompositionVector
{
public:
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::vector_ptrtype vector_ptrtype;

    static const bool use_block_structure = model_type::use_block_structure;

    typedef typename mpl::range_c< int, 0, space_type::nSpaces > rangespace_type;

    template< int T >
    using sub_space = typename space_type::template sub_functionspace<T>::type::element_type;

    template< int Row >
    using oneblock_type = typename mpl::if_< mpl::bool_<use_block_structure>,
                                             BlockAD< model_type, vector_ptrtype, sub_space<Row> >,
                                             BlockAD< model_type, vector_ptrtype, space_type> >::type;
    template< int Row >
    using oneblock_ptrtype = boost::shared_ptr< oneblock_type<Row> >;

    template < typename  T >
    struct CreateVecBlock
    {
        typedef BlockAD< model_type, vector_ptrtype, sub_space<T::value> > block_type;
        typedef boost::shared_ptr<block_type> type;
    };
    typedef typename mpl::if_<mpl::bool_<use_block_structure>,
                              typename mpl::transform< rangespace_type,
                                                       CreateVecBlock<mpl::_1>,
                                                       mpl::back_inserter<fusion::vector<> > >::type,
                              fusion::vector<boost::shared_ptr<BlockAD<model_type,
                                                                       vector_ptrtype,
                                                                       space_type > > > >::type blockvector_ad_type;

    AffineDecompositionVector()
        {}

    AffineDecompositionVector( model_ptrtype model ) :
        M_model( model )
        {}

    AffineDecompositionVector( AffineDecompositionVector const& o ) :
        M_model( o.M_model ),
        M_blockAD( o.M_blockAD )
        {}

    template< int Row=0 >
    void createBlock()
        {
            return createBlock<Row>( mpl::bool_<use_block_structure>() );
        }
    template< int Row >
    void createBlock( mpl::bool_<true> )
        {
            if ( !fusion::at_c<Row>( M_blockAD ) )
                fusion::at_c<Row>( M_blockAD ) = oneblock_ptrtype<Row>(
                    new oneblock_type<Row>( M_model, M_model->functionSpace()->template functionSpace<Row>() ) );
        }

    template< int Row >
    void createBlock( mpl::bool_<false> )
        {
            CHECK( Row==0 ) << "Error : you want to create block "<< Row <<" in the affine decomposition, but you are not using block structure";
            if ( !fusion::at_c<0>( M_blockAD ) )
                fusion::at_c<0>( M_blockAD ) = oneblock_ptrtype<0>(
                    new oneblock_type<0>( M_model, M_model->functionSpace() ) );
        }

    template< int Row=0 >
    oneblock_ptrtype<Row> get() const
        {
            CHECK( fusion::at_c<Row>( M_blockAD ) )<< "Error trying to access block "<< Row <<" when it was not created yet.";

            return fusion::at_c<Row>( M_blockAD );
        }

    struct checkBlock
    {
        template <typename T>
        void operator() ( T& t ) const
            {
                if( t )
                    t->check();
            }
    };
    void check()
        {
            fusion::for_each( M_blockAD, checkBlock() );
        }

private :
    model_ptrtype M_model;
    blockvector_ad_type M_blockAD;

};

template <typename ModelType >
class AffineDecompositionMatrix :
        public boost::enable_shared_from_this< AffineDecompositionMatrix<ModelType> >
{
public :
    typedef AffineDecompositionMatrix<ModelType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    static const bool use_block_structure = model_type::use_block_structure;

    typedef typename mpl::range_c< int, 0, space_type::nSpaces > rangespace_type;

    template< int T >
    using sub_space = typename space_type::template sub_functionspace<T>::type::element_type;

    template< int Row, int Col >
    using oneblock_type = typename mpl::if_< mpl::bool_<use_block_structure>,
                                             BlockAD< model_type, sparse_matrix_ptrtype,
                                                      sub_space<Row>, sub_space<Col> >,
                                             BlockAD< model_type, sparse_matrix_ptrtype,
                                                      space_type, space_type > >::type;

    template< int Row, int Col >
    using oneblock_ptrtype = boost::shared_ptr< oneblock_type<Row,Col> >;

    template< typename Row>
    struct CreateMatBlock
    {
        template< typename Col >
        struct CreateRow
        {
            typedef BlockAD< model_type, sparse_matrix_ptrtype,
                             sub_space<Row::value>, sub_space<Col::value> > row_type;
            typedef boost::shared_ptr< row_type > type;
        };
        typedef typename mpl::transform< rangespace_type,
                                         CreateRow<mpl::_1>,
                                         mpl::back_inserter< fusion::vector<> > >::type type;
    };

    typedef typename mpl::if_< mpl::bool_<use_block_structure>,
                              typename mpl::transform<rangespace_type,
                                                      CreateMatBlock<mpl::_1>,
                                                      mpl::back_inserter< fusion::vector<> > >::type,
                               fusion::vector<fusion::vector<boost::shared_ptr< BlockAD<model_type,
                                                                                        sparse_matrix_ptrtype,
                                                                                        space_type, space_type >>>>
                               >::type blockmatrix_ad_type;


    AffineDecompositionMatrix()
        {}

    AffineDecompositionMatrix( model_ptrtype model ) :
        M_model( model )
        {}

    AffineDecompositionMatrix( AffineDecompositionMatrix const& o ) :
        M_model( o.M_model ),
        M_blockAD( o.M_blockAD )
        {}

    template< int Row=0, int Col=0 >
    void createBlock()
        {
            return createBlock<Row,Col>( mpl::bool_<use_block_structure>() );
        }
    template< int Row, int Col >
    void createBlock( mpl::bool_<true> )
        {
            if ( !fusion::at_c<Col>( fusion::at_c<Row>( M_blockAD ) ) )
                fusion::at_c<Col>( fusion::at_c<Row>( M_blockAD ) ) = oneblock_ptrtype<Row,Col>(
                    new oneblock_type<Row,Col>( M_model,
                                                M_model->functionSpace()->template functionSpace<Row>(),
                                                M_model->functionSpace()->template functionSpace<Col>() ) );
        }

    template< int Row, int Col >
    void createBlock( mpl::bool_<false> )
        {
            CHECK( Row==0 && Col==0 ) << "Error : you want to create block "<< Row<<","<<Col <<" in the affine decomposition, but you are not using block structure";

            if ( !fusion::at_c<0>( fusion::at_c<0>( M_blockAD ) ) )
                 fusion::at_c<0>( fusion::at_c<0>( M_blockAD ) ) = oneblock_ptrtype<0,0>(
                     new oneblock_type<0,0>( M_model, M_model->functionSpace() ) );
        }

    template< int Row=0, int Col=0 >
    oneblock_ptrtype< Row, Col > get() const
        {
            //     CHECK( fusion::at_c<Col>( fusion::at_c<Row>( M_blockAD ) )  )<< "Error trying to access block "<< Row<<","<<Col <<" when it was not created yet.";

            return fusion::at_c<Col>( fusion::at_c<Row>( M_blockAD ) ) ;
        }

    void check()
        {
            fusion::for_each( M_blockAD, checkRow() );
        }

    void initializeMassMatrix()
        {
            rangespace_type range;
            initializeMassRow init( this->shared_from_this() );
            fusion::for_each( range, init );
        }

private :
    struct initializeMassRow
    {
        initializeMassRow( self_ptrtype AD) :
            M_AD( AD )
            {}

        template<typename T>
        void operator() ( T& t) const
            {
                initializeMassBlock<T> block( M_AD );
                rangespace_type range;
                fusion::for_each( range, block );
            }

        self_ptrtype M_AD;
    };
    template <typename Row>
    struct initializeMassBlock
    {
        initializeMassBlock( self_ptrtype AD) :
            M_AD( AD )
            {}

        template<typename Col>
        void operator() ( Col& t) const
            {
                M_AD->template createBlock<Row::value,Col::value>();
                M_AD->template get<Row::value,Col::value>()->initializeMassMatrix();
            }
        self_ptrtype M_AD;
    };
    struct checkBlock
    {
        template< typename T>
        void operator() ( T& t) const
            {
                if ( t )
                    t->check();
            }
    };
    struct checkRow
    {
        template <typename T>
        void operator() ( T& t ) const
            {
                fusion::for_each( t, checkBlock() );
            }
    };

    model_ptrtype M_model;
    blockmatrix_ad_type M_blockAD;
};


} //namespace Feel

#endif // define __AFFINE_DECOMPOSITION_H
