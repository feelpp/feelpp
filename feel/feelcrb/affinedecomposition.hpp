#ifndef __AFFINE_DECOMPOSITION_H
#define __AFFINE_DECOMPOSITION_H 1

#include <feel/feel.hpp>
#include <boost/shared_ptr.hpp>

namespace Feel
{
template <typename EimDefinition, typename ModelType>
class AffineDecomposition
{
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
    AffineDecomposition( model_ptrtype model, int pspace_dimension=1 ) :
        M_pspace_dim( pspace_dimension ),
        M_model( model ),
        M_Nl(0),
        M_Qa(0),
        M_Qm(0),
        M_use_ginac_expr( false ),
        M_use_operators_free( false )
        {
            init();
        }

    /**
     * Initialization of the affine decomposition
     */
    void init()
        {
        }

    /**
     * Add terms in the Rhs affine decomposition using Ginac expression
     *
     * @param part the conerned part of the affine decomposition : can be A, M, F or L
     * @param tuple contain the matrix (a bilinear form) and the expression of the coefficient
     * @param row when using block structure : the row of the conserned block
     * @param col when using block structure : the column of the conserned block
     */
    void addLhs( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple,
                 int const row=1, int const col=1 )
        {
            // we initialize the ginac tools if not already done
            if ( M_use_ginac_expr==false )
            {
                buildGinacSymbols();
                M_use_ginac_expr=true;
            }
            if ( row>M_Aqm.size() )
            {
                M_Aqm.resize( row );
                M_ginacAq.resize( row );
                M_Qa.resize( row );
                M_mMaxA.resize( row );
                M_betaAqm.resize(row);
            }
            if ( col>M_Aqm[row-1].size() )
            {
                M_Aqm[row-1].resize( col );
                M_ginacAq[row-1].resize( col );
                M_Qa[row-1].resize( col );
                M_Qa[row-1][col-1]=0;
                M_mMaxA[row-1].resize(col);
                M_betaAqm[row-1].resize(col);
            }

            // add the matrix/vector to the concerned affine decomposition
            int size = M_Aqm[row-1][col-1].size();
            M_Aqm[row-1][col-1].resize( size +1 );
            M_Aqm[row-1][col-1][size].push_back( tuple.template get<0>() );

            // increment the counter of terms
            M_Qa[row-1][col-1]++;
            M_mMaxA[row-1][col-1].push_back( 1 );

            // add the ginac expression
            std::string filename = ( boost::format("GinacA%1%") %size ).str();
            M_ginacAq[row-1][col-1].push_back( expr( tuple.template get<1>(),
                                                       Symbols( M_symbols_vec ),
                                                       filename ) );
            // create zero betaAqm
            M_betaAqm[row-1][col-1].resize(size+1);
            M_betaAqm[row-1][col-1][size].push_back( 0 );
        }

    void addMass( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple,
                 int const row=1, int const col=1 )
        {
            // we initialize the ginac tools if not already done
            if ( M_use_ginac_expr==false )
            {
                buildGinacSymbols();
                M_use_ginac_expr=true;
            }
            if ( row>M_Mqm.size() )
            {
                M_Mqm.resize( row );
                M_ginacMq.resize( row );
                M_Qm.resize( row );
                M_mMaxM.resize( row );
                M_betaMqm.resize(row);
            }
            if ( col>M_Mqm[row-1].size() )
            {
                M_Mqm[row-1].resize( col );
                M_ginacMq[row-1].resize( col );
                M_Qm[row-1].resize( col );
                M_Qm[row-1][col-1]=0;
                M_mMaxM[row-1].resize(col);
                M_betaMqm[row-1].resize(col);
            }

            // add the matrix/vector to the concerned affine decomposition
            int size = M_Mqm[row-1][col-1].size();
            M_Mqm[row-1][col-1].resize( size +1 );
            M_Mqm[row-1][col-1][size].push_back( tuple.template get<0>() );

            // increment the counter of terms
            M_Qm[row-1][col-1]++;
            M_mMaxM[row-1][col-1].push_back( 1 );

            // add the ginac expression
            std::string filename = ( boost::format("GinacM%1%") %size ).str();
            M_ginacMq[row-1][col-1].push_back( expr( tuple.template get<1>(),
                                                       Symbols( M_symbols_vec ),
                                                       filename ) );
            // create zero betaMqm
            M_betaMqm[row-1][col-1].resize(size+1);
            M_betaMqm[row-1][col-1][size].push_back( 0 );
        }
    void addOutput( boost::tuple< vector_ptrtype , std::string > const & tuple, int output=0,
                    int const row=1 )
        {
            // we initialize the ginac tools if not already done
            if ( M_use_ginac_expr==false )
            {
                buildGinacSymbols();
                M_use_ginac_expr=true;
            }

            M_Nl = std::max( M_Nl, output+1);
            if ( row>M_Fqm.size() )
            {
                M_Fqm.resize( row );
                M_ginacFq.resize( row );
                M_Ql.resize( row );
                M_mMaxF.resize( row );
                M_betaFqm.resize(row);
            }
            if ( output>=M_Fqm[row-1].size() )
            {
                M_Fqm[row-1].resize( output+1 );
                M_ginacFq[row-1].resize( output+1 );
                M_Ql[row-1].resize( output+1 );
                M_Ql[row-1][output]=0;
                M_mMaxF[row-1].resize( output+1 );
                M_betaFqm[row-1].resize( output+1 );
            }

            // add the matrix/vector to the concerned affine decomposition
            int size = M_Fqm[row-1][output].size();
            M_Fqm[row-1][output].resize( size +1 );
            M_Fqm[row-1][output][size].push_back( tuple.template get<0>() );

            // increment the counter of terms
            M_Ql[row-1][output]++;
            M_mMaxF[row-1][output].push_back( 1 );

            // add the ginac expression
            std::string filename = ( boost::format("GinacL-%1%%2%") %output %size ).str();
            M_ginacFq[row-1][output].push_back( expr( tuple.template get<1>(),
                                                      Symbols( M_symbols_vec ),
                                                      filename ) );
            // create zero betaFqm
            M_betaFqm[row-1][output].resize(size+1);
            M_betaFqm[row-1][output][size].push_back( 0 );
        }
    /*void addOutput( boost::tuple< functional_ptrtype, std::string > const & tuple, int output=0,
                    int const row=1, int const col=1 )
        {

     }*/

    affine_decomposition_type compute( int row=1, int col=1 )
        {
            blockCheck( M_Mqm, row, col, "compute");
            fatalBlockCheck( M_Aqm, row, col, "compute");
            fatalBlockCheck( M_Fqm, row, "compute");

            return boost::make_tuple( M_Mqm[row-1][col-1],
                                      M_Aqm[row-1][col-1],
                                      M_Fqm[row-1] );
        }

    void count()
        {
            /* TODO */
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
        }

    template <typename ContainerType>
    void fatalBlockCheck( ContainerType& vec, int row, int col, std::string name ) const
        {
            CHECK( row<=vec.size() )<< "[AD."<< name <<"()] Invalid row number="<< row
                                    << " and max row number="<< vec.size()<< std::endl;
            CHECK( col<=vec[row-1].size() )<< "[AD."<< name <<"()] Invalid col number="<< col
                                           << " and max col number="<< vec[row-1].size()<< std::endl;
        }
    template <typename ContainerType>
    void fatalBlockCheck( ContainerType& vec, int row, std::string name ) const
        {
            CHECK( row<=vec.size() )<< "[AD."<< name <<"()] Invalid row number="<< row
                                    << " and max row number="<< vec.size()<< std::endl;
        }

    template <typename ContainerType>
    void fatalBlockCheck( ContainerType& vec, int row, std::string name, int output ) const
        {
            CHECK( row<=vec.size() )<< "[AD."<< name <<"()] Invalid row number="<< row
                                    << " and max row number="<< vec.size()<< std::endl;
            CHECK( output<vec[row-1].size() )<< "[AD."<< name
                                             << "()] Invalid output number="<< output
                                             << " and max output number="<< vec[row-1].size()-1<< std::endl;
        }


    affine_matrix_type Mqm( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_Mqm, row, col, "Mqm");
            return M_Mqm[row-1][col-1];
        }

    affine_matrix_type Aqm( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_Aqm, row, col, "Aqm");
            return M_Aqm[row-1][col-1];
        }

    affine_vector_type Fqm( int output, int row=1 ) const
        {
            fatalBlockCheck(M_Fqm, row, "Fqm", output );
            return M_Fqm[row-1][output];
        }

    int Qm( int row=1, int col=1 ) const
        {
            fatalBlockCheck( M_Qm, row, col, "Qm");
            return M_Qm[row-1][col-1];
        }

    int Qa( int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_Qa, row, col, "Qa");
            return M_Qa[row-1][col-1];
        }

    int Ql( int output, int row=1 ) const
        {
            fatalBlockCheck(M_Ql, row, "Ql", output );
            return M_Ql[row-1][output];
        }

    int Nl() const
        {
            return M_Nl;
        }


    int mMaxM( int q, int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_mMaxM, row, col, "mMaxM");
            double size = M_mMaxM[row-1][col-1].size();
            CHECK( q < size ) << "[AD.mMaxM] called with the bad q index "<<q<<" and max is "<<size<<"\n";
            return M_mMaxM[row-1][col-1][q];
        }

    int mMaxA( int q, int row=1, int col=1 ) const
        {
            fatalBlockCheck(M_mMaxA, row, col, "mMaxA");
            double size = M_mMaxA[row-1][col-1].size();
            CHECK( q < size ) << "[AD.mMaxA] called with the bad q index "<<q<<" and max is "<<size<<"\n";
            return M_mMaxA[row-1][col-1][q];
        }

    int mMaxF( int output, int q, int row=1 ) const
        {
            fatalBlockCheck(M_mMaxF, row, "mMaxF", output);
            CHECK( q<M_mMaxF[row-1][output].size())<< "[AD.mMaxF()] Invalid q index="<< q
                                                   << "and max is "<< M_mMaxF[row-1][output].size()<< std::endl;

            return M_mMaxF[row-1][output][q];
        }


    sparse_matrix_ptrtype M( uint16_type q, uint16_type m, bool transpose = false,
                                     int row=1, int col=1) const
        {
            fatalBlockCheck(M_Mqm, row, col, "M");
            if ( transpose )
                return M_Mqm[row-1][col-1][q][m]->transpose();
            else
                return M_Mqm[row-1][col-1][q][m];
        }
    sparse_matrix_ptrtype A( uint16_type q, uint16_type m, bool transpose = false,
                                     int row=1, int col=1) const
        {
            fatalBlockCheck(M_Aqm, row, col, "A");
            if ( transpose )
                return M_Aqm[row-1][col-1][q][m]->transpose();
            else
                return M_Aqm[row-1][col-1][q][m];
        }
    vector_ptrtype F( uint16_type output, uint16_type q, int m, int row=1 ) const
        {
            fatalBlockCheck(M_Fqm, row, "F", output);
            return M_Fqm[row-1][output][q][m];
        }



    beta_vector_type betaMqm( parameter_type const& mu, double time=0, int row=1, int col=1 )
        {
            fatalBlockCheck(M_betaMqm, row, col, "betaMqm");
            fatalBlockCheck(M_ginacMq, row, col, "betaMqm");
            if ( M_use_ginac_expr )
            {
                std::string symbol;
                for( int i=0; i<M_pspace_dim; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;

                CHECK( M_betaMqm[row-1][col-1].size()==M_ginacMq[row-1][col-1].size() )
                    << "[AD.betaMqm()] Error : not same size : "
                    << "M_betaMqm["<<row-1<<"]["<<col-1<<"]="<<M_betaMqm[row-1][col-1].size()
                    << "and M_ginacMq["<<row-1<<"]["<<col-1<<"]="<<M_ginacMq[row-1][col-1].size()
                    << std::endl;

                for ( int i=0; i<M_Qm[row-1][col-1]; i++ )
                    M_betaMqm[row-1][col-1][i][0]=M_ginacMq[row-1][col-1][i].evaluate(M_symbols_map);

                return M_betaMqm[row-1][col-1];
            }
            else
            {
                CHECK( M_betaMqm[row-1][col-1].size()>0 )<< "You did not fill betaMqm coefficient, you use either ginac expression or betaMqm() function to do so"<<std::endl;
                return M_betaMqm[row-1][col-1];
            }
        }

    beta_vector_type betaAqm( parameter_type const& mu, double time=0, int row=1, int col=1 )
        {
            fatalBlockCheck(M_betaAqm, row, col, "betaAqm");
            fatalBlockCheck(M_ginacAq, row, col, "betaAqm");
            if ( M_use_ginac_expr )
            {
                std::string symbol;
                for( int i=0; i<M_pspace_dim; i++ )
                {
                    symbol = symbol = ( boost::format("mu%1%") %i ).str();
                    M_symbols_map[symbol]=mu(i);
                }
                M_symbols_map["t"]=time;

                CHECK( M_betaAqm[row-1][col-1].size()==M_ginacAq[row-1][col-1].size() )
                    << "[AD.betaAqm()] Error : not same size : "
                    << "M_betaAqm["<<row-1<<"]["<<col-1<<"]="<<M_betaAqm[row-1][col-1].size()
                    << "and M_ginacAq["<<row-1<<"]["<<col-1<<"]="<<M_ginacAq[row-1][col-1].size()
                    << std::endl;

                for ( int i=0; i<M_Qa[row-1][col-1]; i++ )
                    M_betaAqm[row-1][col-1][i][0]=M_ginacAq[row-1][col-1][i].evaluate(M_symbols_map);

                return M_betaAqm[row-1][col-1];
            }
            else
            {
                CHECK( M_betaAqm[row-1][col-1].size()>0 )<< "You did not fill betaAqm coefficient, you can use either ginac expressions or betaAqm() function to do so"<<std::endl;
                return M_betaAqm[row-1][col-1];
            }
        }

    std::vector< beta_vector_type > betaFqm( parameter_type const& mu, double time=0, int row=1 )
        {
            if ( M_use_ginac_expr )
            {
                for ( int output=0; output<M_Nl; output++)
                {
                    fatalBlockCheck( M_betaFqm, row, "betaFqm", output );
                    fatalBlockCheck( M_ginacFq, row, "betaFqm", output );
                    std::string symbol;
                    for( int i=0; i<M_pspace_dim; i++ )
                    {
                        symbol = symbol = ( boost::format("mu%1%") %i ).str();
                        M_symbols_map[symbol]=mu(i);
                    }
                    M_symbols_map["t"]=time;

                    CHECK( M_betaFqm[row-1][output].size()==M_ginacFq[row-1][output].size() )
                        << "[AD.betaFqm()] Error : not same size : "
                        << "M_betaFqm["<<row-1<<"]="<<M_betaFqm[row-1][output].size()
                        << "and M_ginacFq["<<row-1<<"]="<<M_ginacFq[row-1][output].size()
                        << std::endl;

                    for ( int i=0; i<M_ginacFq[row-1][output].size(); i++ )
                        M_betaFqm[row-1][output][i][0]
                            =M_ginacFq[row-1][output][i].evaluate(M_symbols_map);
                }
            }
            else
            {
                CHECK( M_betaFqm[row-1].size()>0 )<< "You did not fill betaFqm coefficients, you can use either ginac expressions or betaFqm() function to do so"<<std::endl;
            }
            return M_betaFqm[row-1];
        }


private:
    /**
     * Initialize Ginac tools
     */
    void buildGinacSymbols()
        {
            for( int i=0; i<M_pspace_dim; i++)
            {
                std::string symbol = ( boost::format("mu%1%") %i ).str();
                M_symbols_vec.push_back( symbol );
                M_symbols_map.insert( std::pair< std::string,double> (symbol,0) );
            }
            M_symbols_vec.push_back( "t" );
            M_symbols_map.insert( std::pair< std::string,double> ("t",0) );
        }


    int M_pspace_dim;
    model_ptrtype M_model;

    int M_Nl; ///< Number of Outputs

    /// number of terms in the AD
    std::vector< std::vector<int>> M_Qa; // row x col
    std::vector< std::vector<int>> M_Qm; // row x col
    std::vector< std::vector<int>> M_Ql; // row x output

    /// number of eim term for each term of the AD
    std::vector< std::vector< std::vector< int >>> M_mMaxA; // row x col x q
    std::vector< std::vector< std::vector< int >>> M_mMaxM; // row x col x q
    std::vector< std::vector< std::vector< int >>> M_mMaxF; // row x output x q

    /// affine decompositions
    block_affine_matrix_type M_Aqm; // row x col x q x m
    block_affine_matrix_type M_Mqm; // row x col x q x m
    block_affine_output_type M_Fqm; // row x output x q x m

    /// beta coefficients
    std::vector< std::vector< beta_vector_type >> M_betaMqm; // row x col x q x m
    std::vector< std::vector< beta_vector_type >> M_betaAqm; // row x col x q x m
    std::vector< std::vector< beta_vector_type >> M_betaFqm; // row x output x q x m

    /// ginac expressions of the beta coefficient
    std::vector< std::vector< std::vector< Expr<GinacEx<2> >>>> M_ginacAq; // row x col x q
    std::vector< std::vector< std::vector< Expr<GinacEx<2> >>>> M_ginacMq; // row x col x q
    std::vector< std::vector< std::vector< Expr<GinacEx<2> >>>> M_ginacFq; // row x output x q

    /// Ginac tools
    bool M_use_ginac_expr;
    std::vector< std::string > M_symbols_vec;
    std::map<std::string,double> M_symbols_map;

    bool M_use_operators_free;
}; // class AffineDecompostion


} // namespace Feel

#endif // define __AFFINE_DECOMPOSITION_H
