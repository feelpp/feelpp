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
        M_Nl(0),
        M_use_ginac_expr( false ),
        M_use_operators_free( false )
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
            if ( !M_use_ginac_expr )
            {
                buildSymbolsMap();
                M_use_ginac_expr=true;
            }

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


    affine_decomposition_type compute( int row=1, int col=1 )
        {
            blockCheck( M_M, row, col, "compute");
            fatalBlockCheck( M_A, row, col, "compute");
            fatalBlockCheck( M_F, row, "compute");

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
            fatalBlockCheck(M_F, row, output, "Fqm" );
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
            fatalBlockCheck(M_F, row, output, "Ql" );
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
            fatalBlockCheck(M_F, row, output, "mMaxF");
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
            fatalBlockCheck(M_F, row, output, "F" );
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
            if ( M_use_ginac_expr )
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
            M_col( col )
            {}

        void add( TensorType const& mat, std::string const& symbol , std::string filename )
            {
                newEntry( symbol, filename );
                M_Mqm.resize( M_Q );
                M_Mqm[M_Q-1].push_back( mat );
            }
        void add(  operator_ptrtype const& ope, std::string const& symbol , std::string filename )
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

        /// beta coefficients
        beta_vector_type M_beta;
        /// ginac expressions
        std::vector< Expr< GinacEx<2> >> M_ginac;

    private:
        void newEntry( std::string symbol, std::string filename )
            {
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
    void blockCheck( ContainerType& vec, int row, int col, std::string name )
        {
            if ( vec.size()<row)
            {
                LOG(INFO)<<  "[AD."<< name <<"()] Call of nonexistant row="<< row
                         << ". Creation of an empty one";
                vec.resize(row);
            }
            if ( vec[row-1].size()<col )
            {
                LOG(INFO)<< "[AD."<< name <<"()] Call of nonexistant col="<< col
                         << ". Creation of an empty one for row="<< row;
                vec[row-1].resize(col);
            }
            vec[row-1][col-1] = oneblockmatrix_ptrtype( new oneblockmatrix_type( M_model, row, col ) );
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
    bool M_use_ginac_expr;
    std::map<std::string,double> M_symbols_map;

    bool M_use_operators_free;
}; // class AffineDecompostion


} // namespace Feel

#endif // define __AFFINE_DECOMPOSITION_H
