#ifndef __AFFINE_DECOMPOSITION_H
#define __AFFINE_DECOMPOSITION_H 1

#include <feel/feel.hpp>
#include <boost/shared_ptr.hpp>

namespace Feel
{
template <typename EimDefinition, typename ModelType>
class AffineDecomposition
{
    template <typename TensorType, typename CompositeType>
    struct OneBlockAD;

public:
    typedef double value_type;
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename model_type::element_type element_type;
    typedef typename model_type::mesh_type mesh_type;
    typedef typename model_type::mesh_ptrtype mesh_ptrtype;

    /// Backend and Tensors
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename model_type::parameter_type parameter_type;

    typedef typename model_type::operator_type operator_type;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;
    typedef typename model_type::functional_type functional_type;
    typedef boost::shared_ptr<functional_type> functional_ptrtype;

    typedef typename model_type::operatorcomposite_type operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;
    typedef typename model_type::functionalcomposite_type functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef OneBlockAD< sparse_matrix_ptrtype,operatorcomposite_ptrtype > oneblockmatrix_type;
    typedef boost::shared_ptr< oneblockmatrix_type > oneblockmatrix_ptrtype;
    typedef OneBlockAD< vector_ptrtype,functionalcomposite_ptrtype > oneblockvector_type;
    typedef boost::shared_ptr< oneblockvector_type > oneblockvector_ptrtype;

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

    /**
     * Constructor
     *
     * @param pspace_dimension The dimension of the parameter space
     */
    AffineDecomposition( model_ptrtype model ) :
        M_model( model ),
        M_Nl(0)
        {
            init();
        }




    /**
     * Initialization of the affine decomposition
     */
    void init()
        {}

    template < typename OperatorType >
    void addMass( OperatorType const& ope, std::string const& symbol,
                  int const row=1, int const col=1 )
        {
            // we resize the block AD and create the concerned block if necessary
            if ( row>M_M.size() )
                M_M.resize( row );
            if ( col>M_M[row-1].size() )
                M_M[row-1].resize( col );
            if ( !M_M[row-1][col-1] )
                M_M[row-1][col-1] = oneblockmatrix_ptrtype ( new oneblockmatrix_type(M_model,row,col) );

            std::string filename = ( boost::format( "M%1%%2%-" ) %row %col ) .str();
            M_M[row-1][col-1]->add( ope, symbol, filename );

        }

    template < typename OperatorType >
    void addLhs( OperatorType const& ope, std::string const& symbol,
                 int const row=1, int const col=1 )
        {
            // we initialize the ginac tools if not already done
            if ( M_symbols_map.size() ==0)
                buildSymbolsMap();

            // we resize the block AD and create the concerned block if necessary
            if ( row>M_A.size() )
                M_A.resize( row );
            if ( col>M_A[row-1].size() )
                M_A[row-1].resize( col );
            if ( !M_A[row-1][col-1] )
                M_A[row-1][col-1] = oneblockmatrix_ptrtype ( new oneblockmatrix_type(M_model,row,col) );

            std::string filename = ( boost::format( "A%1%%2%-" ) %row %col ) .str();
            M_A[row-1][col-1]->add( ope, symbol, filename );
        }

    template < typename FunctionalType >
    void addOutput( FunctionalType const& fun, std::string const& symbol, int output,
                    int const row=1 )
        {
            // we initialize the ginac tools if not already done
            if ( M_symbols_map.size() ==0)
                buildSymbolsMap();

            M_Nl = std::max( M_Nl, output+1 );

            if ( row>M_F.size() )
                M_F.resize( row );
            if ( output>=M_F[row-1].size() )
                M_F[row-1].resize( output+1 );
            if ( !M_F[row-1][output] )
                M_F[row-1][output] = oneblockvector_ptrtype ( new oneblockvector_type(M_model,row) );
            std::string filename = ( boost::format( "F%1%-%2%-" ) %row %output ) .str();
            M_F[row-1][output]->add( fun, symbol, filename );
        }

    void initializeMassMatrix( int row=1, int col=1 )
        {
            if ( row>M_M.size() )
                M_M.resize( row );
            if ( col>M_M[row-1].size() )
                M_M[row-1].resize( col );
            if ( !M_M[row-1][col-1] )
                M_M[row-1][col-1] = oneblockmatrix_ptrtype ( new oneblockmatrix_type(M_model,row,col) );

            M_M[row-1][col-1]->initializeMassMatrix();
        }

    affine_decomposition_type compute( int row=1, int col=1 )
        {
            fatalBlockCheck( M_A, row, col, "compute");
            fatalBlockCheck( M_F, row, "compute");
            fatalBlockCheck( M_M, row, col, "compute");

            affine_output_type Fqm;
            for ( int output=0; output<M_Nl; output++)
                Fqm.push_back( M_F[row-1][output]->compute() );

            return boost::make_tuple( M_M[row-1][col-1]->compute(),
                                      M_A[row-1][col-1]->compute(),
                                      Fqm );
        }

    void count()
        {
            /* TODO ? */
        }


    affine_matrix_type Mqm( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_M, row, col, "Mqm");
            return M_M[row-1][col-1]->compute();
        }
    affine_matrix_type Aqm( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_A, row, col, "Aqm");
            return M_A[row-1][col-1]->compute();
        }
    affine_vector_type Fqm( int output, int row=1 ) const
        {
            fatalBlockCheck(M_F, row, output+1, "Fqm" );
            return M_F[row-1][output]->compute();
        }


    int Qm( int row=1, int col=1 ) const
        {
            fatalBlockCheck( M_M, row, col, "Qm");
            return M_M[row-1][col-1]->Q() ;
        }
    int Qa( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_A, row, col, "Qa");
            return M_A[row-1][col-1]->Q();
        }
    int Ql( int output, int row=1 ) const
        {
            fatalBlockCheck(M_F, row, output+1, "Ql" );
            return M_F[row-1][output]->Q();
        }
    int Nl() const
        {
            return M_Nl;
        }


    int mMaxM( int q, int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_M, row, col, "mMaxM" );
            return M_M[row-1][col-1]->mMax(q);
        }
    int mMaxA( int q, int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_A, row, col, "mMaxA" );
            return M_A[row-1][col-1]->mMax(q);
        }

    int mMaxF( int output, int q, int row=1 ) const
        {
            fatalBlockCheck(M_F, row, output+1, "mMaxF");
            return M_F[row-1][output]->mMax(q);
        }


    sparse_matrix_ptrtype M( uint16_type q, uint16_type m, bool transpose = false,
                                     int row=1, int col=1) const
        {
            fatalBlockCheck(M_M, row, col, "M");
            return M_M[row-1][col-1]->compute( q, m, transpose );
        }
    sparse_matrix_ptrtype A( uint16_type q, uint16_type m, bool transpose = false,
                                     int row=1, int col=1) const
        {
            fatalBlockCheck(M_A, row, col, "A");
            return M_A[row-1][col-1]->compute( q, m, transpose );
        }
    vector_ptrtype F( uint16_type output, uint16_type q, int m, int row=1 ) const
        {
            fatalBlockCheck(M_F, row, output+1, "F" );
            return M_F[row-1][output]->compute(q,m);
        }

    beta_vector_type betaMqm( parameter_type const& mu, double time=0, int row=1, int col=1 )
        {
            return betaQm( M_M, mu, time, row, col );
        }
    beta_vector_type betaAqm( parameter_type const& mu, double time=0, int row=1, int col=1 )
        {
            return betaQm( M_A, mu, time, row, col );
        }
    std::vector< beta_vector_type > betaFqm( parameter_type const& mu, double time=0, int row=1 )
        {
            std::vector< beta_vector_type > betaF;

            for ( int output=0; output<M_Nl; output++ )
                betaF.push_back( betaQm( M_F, mu, time, row, output+1 ) );

            return betaF;
        }
    template <typename AdType>
    beta_vector_type betaQm( AdType AD, parameter_type const& mu, double time, int row, int col )
        {
            fatalBlockCheck(AD, row, col, "betaQm");
            if ( AD[row-1][col-1]->useGinacExpr() )
            {
                std::string symbol;
                for( int i=0; i<M_model->ParameterSpaceDimension; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;

                for ( int i=0; i<AD[row-1][col-1]->Q(); i++ )
                    AD[row-1][col-1]->M_beta[i][0]
                        = AD[row-1][col-1]->M_ginac[i].evaluate(M_symbols_map);

                return AD[row-1][col-1]->M_beta;
            }
            else
            {
                CHECK( AD[row-1][col-1]->M_beta.size()>0 )<< "You did not fill betaMqm coefficient, you use either ginac expression or betaMqm() function to do so"<<std::endl;

                return AD[row-1][col-1]->M_beta;
            }
        }

    offline_merge_type offlineMerge( betaqm_type const& all_beta, bool only_time_dependent_terms,
                                     int row=1, int col=1 )
        {

            fatalBlockCheck( M_A, row, col, "offlineMerge" );
            fatalBlockCheck( M_M, row, col, "offlineMerge" );

            sparse_matrix_ptrtype A;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F( M_Nl );
            if ( !only_time_dependent_terms )
            {
                M = M_M[row-1][col-1]->merge( all_beta.template get<0>() );
                A = M_A[row-1][col-1]->merge( all_beta.template get<1>() );
            }

            for ( int output =0; output<M_Nl; output++ )
            {
                fatalBlockCheck( M_F, row, output+1, "offlineMerge" );
                F[output] = M_F[row-1][output]->merge( all_beta.template get<2>()[output] );
            }

            return boost::make_tuple( M, A, F );
        }


private:

    template <typename TensorType, typename CompositeType>
    struct OneBlockAD
    {
    public:
        static const bool is_matrix = std::is_same<TensorType,sparse_matrix_ptrtype>::value;

        OneBlockAD( model_ptrtype model, int row, int col=0 ) :
            M_model(model),
            M_Q( 0 ),
            M_use_operators_free( false ),
            M_row( row ),
            M_col( col ),
            M_use_ginac_expr( false )
            {}

        void add( TensorType const& mat, std::string const& symbol , std::string filename )
            {
                newEntry( symbol, filename );
                M_Mqm.resize( M_Q );
                M_Mqm[M_Q-1].push_back( mat );
            }
        void add( operator_ptrtype const& ope, std::string const& symbol , std::string filename )
            {
                M_use_operators_free = true;
                if ( !M_composite )
                    M_composite = opLinearComposite( _imageSpace=M_model->functionSpace(M_row-1),
                                                     _domainSpace=M_model->functionSpace(M_col-1) );

                newEntry( symbol, filename );
                filename += std::to_string( M_Q );
                ope->setName( filename );
                M_composite->addElement( {M_Q-1,0}, ope );
            }
        void add( functional_ptrtype const& fun, std::string const& symbol , std::string filename )
            {
                M_use_operators_free = true;
                if ( !M_composite )
                    M_composite = functionalLinearComposite( _space=M_model->functionSpace(M_row-1) );

                newEntry( symbol, filename );
                filename += std::to_string( M_Q );
                fun->setName( filename );
                M_composite->addElement( {M_Q-1,0}, fun );
            }


        void initializeMassMatrix()
            {
                std::string filename = ( boost::format("GinacM%1%%2%-") %M_row %M_col).str() ;

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
                for( int q=0; q<M_Q; q++ )
                {
                    M_beta[q].resize( M_mMax[q] );
                    for( int m=0; m<M_mMax[q]; m++ )
                        M_beta[q][m]=1.;
                }
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
                }
                else
                    matrix = M_Mqm[q][m];

                if ( transpose )
                    return matrix->transpose();
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
                auto M = M_model->newMatrix();
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
                        {
                            M->addMatrix( beta[q][m], M_Mqm[q][m] );
                        }
                }
                return M;
            }
        TensorType merge( beta_vector_type const& beta, mpl::bool_<false> )
            {
                auto F = M_model->newVector();
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


        bool useOperatorsFree() { return M_use_operators_free; }
        bool useGinacExpr() { return M_use_ginac_expr; }

        /// beta coefficients
        beta_vector_type M_beta;
        /// ginac expressions
        std::vector< Expr< GinacEx<2> >> M_ginac;

    private:
        void newEntry( std::string symbol, std::string filename )
            {
                M_use_ginac_expr = true;
                if ( M_symbols_vec.size()==0 )
                    buildSymbolsVector();
                filename = "Ginac" + filename + std::to_string( M_Q );
                M_Q++;
                M_mMax.push_back( 1 );
                M_ginac.push_back( expr( symbol,
                                         Symbols( M_symbols_vec ),
                                         filename ) );
                // create zero betaQm
                M_beta.resize( M_Q );
                M_beta[M_Q-1].push_back(0);
            }

        void assembleMassMatrix()
            {
                const bool is_composite = model_type::functionspace_type::is_composite;
                assembleMassMatrix( mpl::bool_<is_composite>() );
            }
        void assembleMassMatrix( mpl::bool_<false> )
            {
                auto Xh = M_model->functionSpace();
                auto mesh = Xh->mesh();
                auto u = Xh->element();
                auto v = Xh->element();

                M_Q = (int) 1 ;
                M_mMax.resize(M_Q);
                M_mMax[0] = 1;

                M_Mqm.resize( M_Q );
                M_Mqm[0].resize( 1 );
                M_Mqm[0][0] = backend()->newMatrix( _test=Xh , _trial=Xh );

                form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[0][0] ) =
                    integrate( _range=elements( mesh ), _expr=idt( u )*id( v )  );
                M_Mqm[0][0]->close();
            }
        void assembleMassMatrix( mpl::bool_<true> )
            {
                index_vector_type index_vector;
                int n_spaces = model_type::functionspace_type::nSpaces;
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
                const bool is_composite = model_type::functionspace_type::is_composite;
                preAssembleMassMatrix( mpl::bool_<is_composite>() );

                M_mMax = M_composite->countAllContributions();
                M_Q = M_mMax.size();
            }
        void preAssembleMassMatrix( mpl::bool_<false> )
            {
                auto Xh = M_model->functionSpace();
                auto u = Xh->element();
                auto v = Xh->element();

                auto mesh = Xh->mesh();
                auto expr=integrate( _range=elements( mesh ) , _expr=idt( u )*id( v ) );

                auto M_composite = opLinearComposite( _domainSpace=Xh, _imageSpace=Xh  );
                auto opfree = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr );

                opfree->setName("mass operator (automatically created)");
                M_composite->addElement( boost::make_tuple(0,0) , opfree );
            }
        void preAssembleMassMatrix( mpl::bool_<true> )
            {
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
                    mesh_ptrtype mesh = Xh->mesh();

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

        model_ptrtype M_model;
        /// number of terms
        int M_Q;
        /// number of eim terms
        std::vector<int> M_mMax;
        /// affine decomposition
        std::vector< std::vector< TensorType >> M_Mqm;
        /// operators free
        CompositeType M_composite;

        bool M_use_operators_free;
        int M_row;
        int M_col;
        std::vector< std::string > M_symbols_vec;

        bool M_use_ginac_expr;
    };

    /**
     * Initialize Ginac tools
     */
    void buildSymbolsMap()
        {
            for( int i=0; i<M_model->ParameterSpaceDimension; i++)
            {
                std::string symbol = ( boost::format("mu%1%") %i ).str();
                M_symbols_map.insert( std::pair< std::string,double> (symbol,0) );
            }
            M_symbols_map.insert( std::pair< std::string,double> ("t",0) );
        }
    template <typename ContainerType>
    void fatalBlockCheck( ContainerType& vec, int row, int col, std::string name ) const
        {
            int maxRow = vec.size();
            int maxCol = vec[maxRow-1].size();
            CHECK( row<=maxRow && col<= maxCol )
                << "[AD."<< name <<"()] Invalid block entries. You asked for "<< row <<","<<col
                << " and actual size is "<<maxRow<<","<<maxCol<<std::endl;
        }
    template <typename ContainerType>
    void fatalBlockCheck( ContainerType& vec, int row, std::string name ) const
        {
            int maxRow = vec.size();
            CHECK( row<=maxRow )<< "[AD."<< name <<"()] Invalid row number="<< row
                   << " and max row number="<< maxRow << std::endl;
        }

    model_ptrtype M_model;

    int M_Nl; ///< Number of Outputs

    std::vector< std::vector< oneblockmatrix_ptrtype >> M_M; // rowXcol
    std::vector< std::vector< oneblockmatrix_ptrtype >> M_A; // rowXcol
    std::vector< std::vector< oneblockvector_ptrtype >> M_F; // rowXoutput

    /// Ginac tools
    std::map<std::string,double> M_symbols_map;
}; // class AffineDecompostion


} // namespace Feel

#endif // define __AFFINE_DECOMPOSITION_H
