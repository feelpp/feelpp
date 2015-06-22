/* -*- mode: c++ -*-

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2009-11-24

 Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 * @file   crbsaddlepoint.hpp
 * @date   Mon Nov 24 11:25:55 2014
 *
 * @brief
 *
 *
 */

#ifndef __CRBSADDLEPOINT_H
#define __CRBSADDLEPOINT_H 1

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>

#define QA(R,C) this->M_model->template Qa<R,C>()
#define QL(R,O) this->M_model->template Ql<R>(O)

namespace Feel
{

po::options_description crbSaddlePointOptions( std::string const& prefix = "" );

/**
 * \class CRBSaddlePoint
 * \brief Certfified Reduced BAsis for Saddle Point Problems
 *
 * \author JB Wahl
 */


template<typename TruthModelType>
class CRBSaddlePoint :
        public CRB<TruthModelType>,
        public boost::enable_shared_from_this<CRBSaddlePoint<TruthModelType>>
{
    typedef CRB<TruthModelType> super_crb;

    typedef CRBSaddlePoint<TruthModelType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
public:
    //@{ // Truth Model
    typedef TruthModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename model_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //@}

    //@{ /// Parameter Space
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    //@}

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;
    typedef double value_type;
    typedef typename convergence_type::value_type convergence;

    //@{ /// Function Space and Elements
    typedef typename model_type::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;
    //@}

    //@{ Backend and Matrix
    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;
    //@}

    //@{ /// Eigen Objects
    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic, 1> > map_dense_vector_type;
    typedef boost::tuple< std::vector<vectorN_type>,
                          std::vector<vectorN_type>,
                          std::vector<vectorN_type>,
                          std::vector<vectorN_type> > solutions_tuple;
    typedef boost::tuple< double,double,double,
                          std::vector< std::vector< double > >,
                          std::vector< std::vector< double > > > upper_bounds_tuple;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant
    //@}

    //@{ /// Exporter
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    //@}

    template<int T>
    using subspace_type = typename model_type::template subspace_type<T>;
    template<int T>
    using subspace_ptrtype = boost::shared_ptr<subspace_type<T>>;

    template <int T>
    using rbspace_ptrtype = typename model_type::template rbspace_ptrtype<T>;

    template<int T>
        using crbelementsdb_type = CRBElementsDB<subspace_type<T>,model_type>;
    template<int T>
    using crbelementsdb_ptrtype = boost::shared_ptr<crbelementsdb_type<T>>;

    template<int T>
    using subelement_type = typename subspace_type<T>::element_type;

    typedef fusion::vector< mpl::vector_c<int,0,0>, mpl::vector_c<int,0,1>,
                            mpl::vector_c<int,1,0>, mpl::vector_c<int,1,1> > matrixblockrange_type;
    typedef mpl::vector_c<int,0,1> vectorblockrange_type;

    typedef std::vector< std::vector< std::vector< std::vector< matrixN_type >>>> blockmatrixN_type;
    typedef std::vector< std::vector< std::vector< vectorN_type >>> blockvectorN_type;

    typedef typename model_type::matrix_ad_type matrix_ad_type;
    typedef typename model_type::vector_ad_type vector_ad_type;

    typedef boost::tuple<double, parameter_type, double, double> max_error_type;
    typedef boost::tuple<double, std::vector<double> > residual_error_type;

    /// Default constructor
    CRBSaddlePoint() :
        super_crb()
    {}


    /// constructor from command line options
    CRBSaddlePoint( std::string  name, model_ptrtype const & model ) :
        super_crb( name, model ),
        M_elements_database0(
            ( boost::format( "%1%" ) %ioption("crb.error-type") ).str(),
            name,
            ( boost::format( "%1%-%2%-%3%-elements0" )
              %name % ioption("crb.output-index") %ioption("crb.error-type") ).str(),
            model ),
        M_elements_database1(
            ( boost::format( "%1%" ) %ioption("crb.error-type") ).str(),
            name,
            ( boost::format( "%1%-%2%-%3%-elements1" )
              %name % ioption("crb.output-index") %ioption("crb.error-type") ).str(),
            model ),
        XN0( model->template rBFunctionSpace<0>() ),
        XN1( model->template rBFunctionSpace<1>() )
        {
            this->M_N=0;
            M_N0.push_back( 0 );
            M_N1.push_back( 0 );
        }


    //! copy constructor
    CRBSaddlePoint( CRBSaddlePoint const & o ) :
        super_crb( o )
    {}

    virtual void init()
    {
        this->setTruthModel();
        if ( this->loadDB() )
            LOG(INFO) << "Database " << this->lookForDB() << " available and loaded\n";

        // this will be in the offline step
        // (it's only when we enrich or create the database that we want to
        // have access to elements of the RB)
        M_elements_database0.setMN( M_N0[this->M_N] );
        M_elements_database1.setMN( M_N1[this->M_N] );

        bool load_elements_db= boption(_name="crb.load-elements-database");
        if( load_elements_db )
        {
            bool found = false;
            if( this->M_elements_database0.loadDB() )
            {
                if (Environment::isMasterRank() )
                    std::cout << "database for basis functions 0 "
                              << this->M_elements_database0.lookForDB() << " available and loaded\n";
                LOG(INFO) << "database for basis functions 0 "
                          << this->M_elements_database0.lookForDB() << " available and loaded\n";
                auto basis_functions = this->M_elements_database0.wn();
                XN0->setBasis( basis_functions );
                found = true;
            }
            if( this->M_elements_database1.loadDB() )
            {
                if (Environment::isMasterRank() )
                    std::cout << "database for basis functions 1 "
                              << this->M_elements_database1.lookForDB() << " available and loaded\n";
                LOG(INFO) << "database for basis functions 1 "
                          << this->M_elements_database1.lookForDB() << " available and loaded\n";
                auto basis_functions = this->M_elements_database1.wn();
                XN1->setBasis( basis_functions );
                found = true && found;
            }
            if ( !found )
            {
                if( Environment::isMasterRank() )
                    std::cout<<"Saddle Point : Warning ! No database for basis functions loaded. Start from the begining"<<std::endl;
                LOG( INFO ) <<"no database for basis functions loaded. Start from the begining";
            }
        }
    }

    WorldComm const& worldComm() const { return Environment::worldComm(); }


    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     *
     *\return compute online the lower bound
     *\and also condition number of matrix A
     */
    boost::tuple<std::vector<double>,matrix_info_tuple>
    lb( size_type N, parameter_type const& mu,
        std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ,
        std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold,
        bool print_rb_matrix=false, int K=0 ) const;

    virtual element_type expansionPrimal( vectorN_type const& U_coeff, int const N ) const;
    virtual element_type expansionDual( vectorN_type const& U_coeff, int const N ) const;
    element_type expansionSaddlePoint( vectorN_type const& U_coeff, int const N,
                                       std::vector<subelement_type<0>> WN0,
                                       std::vector<subelement_type<1>> WN1 ) const;

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    convergence_type offline();

    //@{ /// Database
    virtual void saveDB();
    virtual bool loadDB();
    //@}

private :
    void exportBasisFunctions();
    template<int Row>
        void offlineResidual( int N0, int N1, int Nadded0, int Nadded1 );
    template<int Row>
        void initResidualVectors();
    max_error_type maxResidual( int N );
    double onlineResidual( int Ncur, parameter_type const& mu, vectorN_type Un ) const;
    template <int Row>
        double onlineResidual( int Ncur, parameter_type const& mu, vectorN_type Un, bool test=false ) const;
    residual_error_type steadyPrimalResidual( int Ncur, parameter_type const& mu,  vectorN_type const& Un, double time=0 ) const
    {
        double rez = onlineResidual( Ncur, mu, Un );
        std::vector<double> coeffs(3,0);
        return boost::make_tuple( rez, coeffs );
    }

    struct ComputeADElements
    {
        ComputeADElements( self_ptrtype crb ) : M_crb( crb ) {}

        template <typename T>
        void operator() ( T& t ) const
            {
                using row = typename mpl::at_c<T,0>::type;
                using col = typename mpl::at_c<T,1>::type;
                int r = row::value;
                int c = col::value;
                fusion::vector<rbspace_ptrtype<0>,rbspace_ptrtype<1>> XN( M_crb->XN0, M_crb->XN1 );
                auto XNr = fusion::at_c<row::value>( XN );
                auto XNc = fusion::at_c<col::value>( XN );
                std::array<int,2> N = {{M_crb->M_N0[M_crb->M_N],
                                        M_crb->M_N1[M_crb->M_N]}};
                std::array<int,2> Nadded={{M_crb->M_Nadded0, M_crb->M_Nadded1}};

                for( size_type q=0; q<M_crb->M_model->template Qa<row::value,col::value>(false); q++ )
                {
                    M_crb->M_blockAqm_pr[r][c][q][0].conservativeResize( N[r], N[c] );

                    for( size_type i = N[r] - Nadded[r]; i<N[r]; i++)
                        for( size_type j=0; j<N[c]; j++ )
                            M_crb->M_blockAqm_pr[r][c][q][0]( i, j )
                                = M_crb->M_model->template Aqm<row::value,col::value>( q, 0,
                                                                        XNc->primalBasisElement(j),
                                                                        XNr->primalBasisElement(i) );

                    for( size_type j=N[c]-Nadded[c]; j<N[c]; j++ )
                        for( size_type i=0; i<N[r]; i++ )
                            M_crb->M_blockAqm_pr[r][c][q][0]( i, j )
                                = M_crb->M_model->template Aqm<row::value,col::value>( q, 0,
                                                                     XNc->primalBasisElement(j),
                                                                     XNr->primalBasisElement(i) );
                }

                if( c==0 )
                {
                    for( size_type q=0; q<M_crb->M_model->template Ql<row::value>(0); q++)
                    {
                        int m=0;
                        M_crb->M_blockFqm_pr[r][q][m].conservativeResize( N[r] );
                        for( size_type l=1; l<=Nadded[r]; l++ )
                        {
                            int index = N[r]-l;
                            M_crb->M_blockFqm_pr[r][q][m](index) = M_crb->M_model->template Fqm<row::value>( 0,q,m,XNr->primalBasisElement(index) );
                        }
                    }
                    for( size_type q=0;
                         q<M_crb->M_model->template Ql<row::value>( M_crb->M_output_index);
                         q++ )
                    {
                        int m=0;
                        M_crb->M_blockLqm_pr[r][q][m].conservativeResize( N[r] );
                        for( size_type l=1; l<=Nadded[r]; l++ )
                        {
                            int index = N[r]-l;
                            M_crb->M_blockLqm_pr[r][q][m](index) = M_crb->M_model->template Fqm<row::value>( M_crb->M_output_index,q,m,XNr->primalBasisElement(index) );
                        }
                    }

                }
            }
        self_ptrtype M_crb;
    };

    struct ComputeRBElements
    {
        typedef boost::shared_ptr<const self_type> crb_ptrtype;
        ComputeRBElements( crb_ptrtype crb, parameter_type const& mu, int N ) :
            M_crb( crb),
            M_mu( mu ),
            M_N( N ),
            M_A( M_crb->M_N0[M_N] + M_crb->M_N1[M_N],
                 M_crb->M_N0[M_N] + M_crb->M_N1[M_N] ),
            M_F( M_crb->M_N0[M_N] + M_crb->M_N1[M_N] ),
            M_L( M_crb->M_N0[M_N] + M_crb->M_N1[M_N] )
            {
                M_A.setZero( M_crb->M_N0[M_N] + M_crb->M_N1[M_N],
                             M_crb->M_N0[M_N] + M_crb->M_N1[M_N] );
                M_F.setZero( M_crb->M_N0[M_N] + M_crb->M_N1[M_N] );
                M_L.setZero( M_crb->M_N0[M_N] + M_crb->M_N1[M_N] );
            }

        template<typename T>
        void operator() (T& t) const
            {
                using row = typename mpl::at_c<T,0>::type;
                using col = typename mpl::at_c<T,1>::type;
                std::array<int,2> blockSize = {{ M_crb->M_N0[M_N], M_crb->M_N1[M_N] }};
                std::array<int,2> blockStart = {{ 0, M_crb->M_N0[M_N] }};
                int r = row::value;
                int c = col::value;
                int Qa = M_crb->M_model-> template Qa<row::value,col::value>(false);
                int N = M_crb->M_N0[M_N] + M_crb->M_N1[M_crb->M_N];

                if ( Qa!=0 )
                {
                    auto betaA = M_crb->M_model->template computeBetaAqm<row::value,col::value>( M_mu );
                    for( int q=0; q<Qa; q++ )
                    {
                        int mMax = M_crb->M_model->template mMaxA<row::value,col::value>(q,false);
                        for( int m=0; m<mMax; m++ )
                        {
                            M_A.block( blockStart[r], blockStart[c], blockSize[r], blockSize[c] )
                                += betaA[q][m]*M_crb->M_blockAqm_pr[r][c][q][m].block(0, 0, blockSize[r], blockSize[c] );
                        }
                    }
                }

                if ( c==0 )
                {
                    int Qf = M_crb->M_model->template Ql<row::value>(0);
                    if ( Qf!=0 )
                    {
                        auto betaF = M_crb->M_model->template computeBetaFqm<row::value>( 0, M_mu );
                        for( int q=0; q<Qf; q++ )
                        {
                            int mMax = M_crb->M_model->template mMaxF<row::value>( 0, q );
                            for( int m=0; m<mMax; m++ )
                            {
                                M_F.segment( blockStart[r], blockSize[r] )
                                    += betaF[q][m]*M_crb->M_blockFqm_pr[r][q][m].head( blockSize[r]);
                            }
                        }
                    }

                    int Ql = M_crb->M_model->template Ql<row::value>(M_crb->M_output_index);
                    if ( Ql!= 0 )
                    {
                        auto betaL = M_crb->M_model->template computeBetaFqm<row::value>( M_crb->M_output_index, M_mu );
                        for( int q=0; q<Ql; q++ )
                        {
                            int mMax = M_crb->M_model->template mMaxF<row::value>( M_crb->M_output_index, q );
                            for( int m=0; m<mMax; m++ )
                            {
                                M_L.segment( blockStart[r], blockSize[r] )
                                    += betaL[q][m]*M_crb->M_blockLqm_pr[r][q][m].head( blockSize[r]);
                            }
                        }
                    }
                }
            }

        boost::tuple<matrixN_type,vectorN_type,vectorN_type> getElements()
            {
                return boost::make_tuple( M_A, M_F, M_L );
            }


        crb_ptrtype M_crb;
        parameter_type M_mu;
        int M_N;
        mutable matrixN_type M_A;
        mutable vectorN_type M_F;
        mutable vectorN_type M_L;
    };

    void generateSuperSampling();
    bool buildSampling();

    crbelementsdb_type<0> M_elements_database0;
    crbelementsdb_type<1> M_elements_database1;

    int M_Nadded0;
    int M_Nadded1;
    std::vector<int> M_N0;
    std::vector<int> M_N1;

    rbspace_ptrtype<0> XN0;
    rbspace_ptrtype<1> XN1;

    blockmatrixN_type M_blockAqm_pr;
    blockvectorN_type M_blockFqm_pr;
    blockvectorN_type M_blockLqm_pr;

    std::vector< std::vector< std::vector< std::vector< std::vector< double >>>>> M_R_RhsRhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< vectorN_type >>>>> M_R_Lhs0Rhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< vectorN_type >>>>> M_R_Lhs1Rhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs0Lhs0;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs0Lhs1;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs1Lhs1;


    friend class boost::serialization::access;
    template<class Archive>
        void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
        void load( Archive & ar, const unsigned int version ) ;
    BOOST_SERIALIZATION_SPLIT_MEMBER()

}; // class CRBSaddlePoint


template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::convergence_type
CRBSaddlePoint<TruthModelType>::offline()
{
    bool rebuild_database = boption( "crb.rebuild-database" );
    bool orthonormalize_primal = boption(_name="crb.orthonormalize-primal") ;
    parameter_type mu ( this->M_Dmu );
    double delta_pr = 0;
    double delta_du = 0;
    size_type index =0;
    int output_index = this->M_output_index;

    element_ptrtype U( new element_type( this->M_model->functionSpace() ) );
    int proc_number = this->worldComm().globalRank();


    if ( Environment::isMasterRank() )
        std::cout << "Offline CRBSaddlePoint starts, this may take a while until Database is computed.." << std::endl;
    LOG(INFO) << "[CRBSaddlePoint::offline] Starting offline for output " << this->M_output_index << std::endl;
    LOG(INFO) << "[CRBSaddlePoint::offline] initialize underlying finite element model\n";

    if ( rebuild_database || this->M_N==0 )
    {
        // Generate sampling depending of the "crb.sampling-mode" option
        generateSuperSampling();

        this->M_N=0;
        M_N0.clear();
        M_N0.push_back( 0 );
        M_N1.clear();
        M_N1.push_back( 0 );
        this->M_maxerror = 1e10;

        LOG(INFO) << "[CRBSaddlePoint::offline] allocate reduced basis data structures\n";

        M_blockAqm_pr.resize(2);
        M_blockAqm_pr[0].resize(2);
        M_blockAqm_pr[1].resize(2);

        M_blockFqm_pr.resize(2);
        M_blockLqm_pr.resize(2);

        M_blockAqm_pr[0][0].resize( this->M_model->template Qa<0,0>(false) );
        for ( int q=0; q<M_blockAqm_pr[0][0].size(); q++ )
            M_blockAqm_pr[0][0][q].resize( this->M_model->template mMaxA<0,0>(q,false) );

        M_blockAqm_pr[0][1].resize( this->M_model->template Qa<0,1>(false) );
        for ( int q=0; q<M_blockAqm_pr[0][1].size(); q++ )
            M_blockAqm_pr[0][1][q].resize( this->M_model->template mMaxA<0,1>(q,false) );

        M_blockAqm_pr[1][0].resize( this->M_model->template Qa<1,0>(false) );
        for ( int q=0; q<M_blockAqm_pr[1][0].size(); q++ )
            M_blockAqm_pr[1][0][q].resize( this->M_model->template mMaxA<1,0>(q,false) );

        M_blockAqm_pr[1][1].resize( this->M_model->template Qa<1,1>(false) );
        for ( int q=0; q<M_blockAqm_pr[1][1].size(); q++ )
            M_blockAqm_pr[1][1][q].resize( this->M_model->template mMaxA<1,1>(q,false) );

        M_blockFqm_pr[0].resize( this->M_model->template Ql<0>(0) );
        for ( int q=0; q<M_blockFqm_pr[0].size(); q++ )
            M_blockFqm_pr[0][q].resize( this->M_model->template mMaxF<0>( 0, q) );

        M_blockFqm_pr[1].resize( this->M_model->template Ql<1>(0) );
        for ( int q=0; q<M_blockFqm_pr[1].size(); q++ )
            M_blockFqm_pr[1][q].resize( this->M_model->template mMaxF<1>( 0, q ) );

        M_blockLqm_pr[0].resize( this->M_model->template Ql<0>( output_index ) );
        for ( int q=0; q<M_blockLqm_pr[0].size(); q++ )
            M_blockLqm_pr[0][q].resize( this->M_model->template mMaxF<0>( output_index, q ) );

        M_blockLqm_pr[1].resize( this->M_model->template Ql<1>( output_index ) );
        for ( int q=0; q<M_blockLqm_pr[1].size(); q++ )
            M_blockLqm_pr[1][q].resize( this->M_model->template mMaxF<1>( output_index, q ) );
    } // if ( rebuild_database )

    bool use_predefined_WNmu = buildSampling();

    LOG( INFO )<<"[CRBSaddlePoint offline] M_error_type = "<<this->M_error_type;

    if ( use_predefined_WNmu )
    {
        this->M_maxerror = 1e10;
        this->M_iter_max = this->M_WNmu->size();
        mu = this->M_WNmu->at( std::min( this->M_N, this->M_iter_max-1 ) ); // first element
    }
    else if ( this->M_error_type != CRB_NO_RESIDUAL )
    {
        boost::tie( mu, index ) = this->M_Xi->min(true);
    }


    while( this->M_maxerror>this->M_tolerance && this->M_N<this->M_iter_max )
    {
        CRB_COUT<<"Construction of "<<this->M_N+1<<"/"<<this->M_iter_max<<" basis, mu = [";
        for ( int i=0; i<mu.size(); i++ )
            CRB_COUT<<mu( i )<<",";
        CRB_COUT<<"]"<<std::endl;
        M_Nadded0=0;
        M_Nadded1=0;

        LOG(INFO) <<"========================================"<<"\n";
        LOG(INFO) << "N=" << this->M_N << "/"  << this->M_iter_max << "( nb proc : "<<worldComm().globalSize()<<")";

        U->setName( ( boost::format( "fem-primal-N%1%-proc%2%" ) % (this->M_N)  % proc_number ).str() );
        mu.check();
        U->zero();

        LOG(INFO) << "[CRB::offline] solving primal" << "\n";

        tic();
        *U = this->M_model->solve( mu );
        toc( " -- probleme resolution");

        auto u = U->template elementPtr<0>();
        auto p = U->template elementPtr<1>();

        M_N0.push_back( M_N0[this->M_N] );
        M_N1.push_back( M_N1[this->M_N] );
        this->M_N++;

        XN1->addPrimalBasisElement( p );
        XN1->addDualBasisElement( p );
        M_Nadded1=1;
        M_N1[this->M_N] += M_Nadded1;

        tic();
        this->orthonormalize( M_N1[this->M_N], XN1->primalRB(), M_Nadded1, 1,
                              ioption("crb.saddlepoint.orthonormalize1") );
        toc(" -- orthonormalization of RBSpace #1", ioption("crb.saddlepoint.orthonormalize1"));

        XN0->addPrimalBasisElement( u );
        XN0->addDualBasisElement( u );
        M_Nadded0=1;

        if ( boption("crb.saddlepoint.add-false-supremizer") )
        {
            tic();
            auto us = this->M_model->falseSupremizer( mu, p );
            XN0->addPrimalBasisElement( us );
            XN0->addDualBasisElement( us );
            M_Nadded0++;
            toc(" -- false supremizer computation");
        }
        if ( boption("crb.saddlepoint.add-supremizer") )
        {
            tic();
            //auto us = this->M_model->supremizer( mu, p );
            auto us = this->M_model->supremizer( mu, XN1->primalRB().back() );
            XN0->addPrimalBasisElement( us );
            XN0->addDualBasisElement( us );
            M_Nadded0++;
            toc(" -- supremizer computation");
        }

        M_N0[this->M_N] += M_Nadded0;
        tic();
        this->orthonormalize( M_N0[this->M_N], XN0->primalRB(), M_Nadded0, 0,
                              ioption("crb.saddlepoint.orthonormalize0") );
        toc(" -- orthonormalization of RBSpace #0", ioption("crb.saddlepoint.orthonormalize0") );

        matrixblockrange_type matrixrange;
        ComputeADElements compute_elements( this->shared_from_this() );
        fusion::for_each( matrixrange, compute_elements );


        if ( use_predefined_WNmu )
        {
            //remmber that in this case M_iter_max = sampling size
            if( this->M_N < this->M_iter_max )
                mu = this->M_WNmu->at( this->M_N );
        }
        else
        {
            tic();
            offlineResidual<0>( M_N0[this->M_N], M_N1[this->M_N], M_Nadded0, M_Nadded1 );
            toc(" -- offline residual 0 computation");
            tic();
            offlineResidual<1>( M_N0[this->M_N], M_N1[this->M_N], M_Nadded0, M_Nadded1 );
            toc(" -- offline residual 1 computation");

            this->M_WNmu->push_back( mu, index );
            tic();
            this->M_WNmu_complement = this->M_WNmu->complement();
            toc(" -- WNmu complement creation",  FLAGS_v > 0);
            tic();
            boost::tie( this->M_maxerror, mu, delta_pr, delta_du ) = maxResidual( this->M_N );
            toc( " -- max residual computation ");
        }
        this->M_current_mu = mu;

        this->M_rbconv.insert( convergence(this->M_N, boost::make_tuple(this->M_maxerror,
                                                                        delta_pr,delta_du) ) );

        saveDB();
        M_elements_database0.setWn( boost::make_tuple(XN0->primalRB(), XN0->dualRB()) );
        M_elements_database1.setWn( boost::make_tuple(XN1->primalRB(), XN1->dualRB()) );
        M_elements_database0.saveDB();
        M_elements_database1.saveDB();

        CRB_COUT << "============================================================\n";
    } // while( this->M_maxerror>this->M_tolerance && this->M_N<this->M_iter_max )

    /*std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
     std::vector< vectorN_type > uNduold;*/

    vectorN_type uN(M_N0[this->M_N]+M_N1[this->M_N]);
    if (boption("crb.saddlepoint.test-residual"))
    {
        CRB_COUT << "\n TEST RESIDUAL \n";
        CRB_COUT << "WNmu.size()="<<this->M_WNmu->size()<<std::endl;
        for (int k=0; k<this->M_N; k++)
        {
            CRB_COUT << "====================================================\n";
            parameter_type mu_test = this->M_WNmu->at(k);
            uN.setZero();
            uN(k)=1;
            uN(M_N0[this->M_N]+k)=1;

            //lb( this->M_N, mu_test, uN, uNdu, uNold, uNduold );
            CRB_COUT<< "mu_test = [";
            for ( int i=0; i<mu.size(); i++ )
                CRB_COUT<<mu_test( i )<<",";
            CRB_COUT<<"]"<<std::endl;
            double rez_test0 = onlineResidual<0>( this->M_N, mu_test, uN, true );
            double rez_test1 = onlineResidual<1>( this->M_N, mu_test, uN, true );
            CRB_COUT << "rez_test0="<< std::sqrt(rez_test0)
                     << ", rez_test1="<< std::sqrt(rez_test1) <<std::endl;
        }
        CRB_COUT << "====================================================\n";
    }
    if ( boption("crb.visualize-basis") )
    {
        tic();
        exportBasisFunctions();
        toc(" -- export basis functions");
    }
    CRB_COUT << "Number of elements, in first base : "<< M_N0[this->M_N]<<", in second base : "
             <<M_N1[this->M_N]<<std::endl;

    return this->M_rbconv;
} // offline()


template<typename TruthModelType>
typename boost::tuple<std::vector<double>,
                      typename CRBSaddlePoint<TruthModelType>::matrix_info_tuple >
CRBSaddlePoint<TruthModelType>::lb( size_type N, parameter_type const& mu,
                                  std::vector< vectorN_type >& uN,
                                  std::vector< vectorN_type >& uNdu,
                                  std::vector<vectorN_type> & uNold,
                                  std::vector<vectorN_type> & uNduold,
                                  bool print_rb_matrix, int K ) const
{
    LOG(INFO) << "CRBSaddlePoint : online resolution for mu="<<mu;
    uN.resize(1);
    N = std::min( N, this->M_N );

    matrixN_type A ( M_N0[N] + M_N1[N], M_N0[N] + M_N1[N] );
    vectorN_type F ( M_N0[N] + M_N1[N] );
    vectorN_type L ( M_N0[N] + M_N1[N] );

    matrixblockrange_type matrixrange;
    ComputeRBElements compute_elements( this->shared_from_this(), mu, (int) N );
    fusion::for_each( matrixrange, compute_elements );

    boost::tie( A, F, L ) = compute_elements.getElements();
    if ( boption("crb.saddlepoint.transpose") )
    {
        int Qa01 = this->M_model-> template Qa<0,1>(false);
        int Qa10 = this->M_model-> template Qa<1,0>(false);
        if ( Qa01==0 && Qa10!=0 )
            A.block( 0, M_N0[N], M_N0[N], M_N1[N] )=A.block( M_N0[N], 0, M_N1[N], M_N0[N] ).transpose();
        else if( Qa01!=0 && Qa10==0 )
            A.block( M_N0[N], 0, M_N1[N], M_N0[N] )=A.block( 0, M_N0[N], M_N0[N], M_N1[N] ).transpose();
    }

    uN[0]=A.lu().solve(F);

    std::vector<double>output_vector(1);
    output_vector[0] = L.dot( uN[0] );

    auto matrix_info = boost::make_tuple( 0., 0. );

    return boost::make_tuple( output_vector, matrix_info);
} // lb()


template<typename TruthModelType>
template<int Row>
void
CRBSaddlePoint<TruthModelType>::offlineResidual( int N0, int N1, int Nadded0, int Nadded1 )
{
    bool transpose = boption("crb.saddlepoint.transpose");
    bool optimize = boption(_name="crb.optimize-offline-residual") ;

    if ( M_R_RhsRhs.size()==0 )
    {
        initResidualVectors<0>();
        initResidualVectors<1>();
    }

    int QRhs = this->M_model->template Ql<Row>( 0 );
    int QLhs0 = this->M_model->template Qa<Row,0>();
    int QLhs1 = this->M_model->template Qa<Row,1>();

    vector_ptrtype Z1 = this->M_model->template newL2Vector<Row>();
    vector_ptrtype Z2 = this->M_model->template newL2Vector<Row>();
    vector_ptrtype X0 = this->M_model->template newL2Vector<0>();
    vector_ptrtype X1 = this->M_model->template newL2Vector<1>();
    vector_ptrtype Y0 = this->M_model->template newL2Vector<0>();
    vector_ptrtype Y1 = this->M_model->template newL2Vector<1>();
    vector_ptrtype W = this->M_model->template newL2Vector<Row>();

    std::vector< std::vector< vector_ptrtype >> Fqm = this->M_model->template Fqm<Row>(0);
    std::vector< std::vector< sparse_matrix_ptrtype >> Lhs0 = this->M_model->template Aqm<Row,0>();
    std::vector< std::vector< sparse_matrix_ptrtype >> Lhs1 = this->M_model->template Aqm<Row,1>();

    if ( N0==Nadded0 )
    {
        LOG(INFO) << "[offlineResidual] Compute residual data\n";
        LOG(INFO) << "[offlineResidual] M_R_0x0\n";

        for ( int q1=0; q1<QRhs; q1++ )
        {
            for ( int m1=0; m1< this->M_model->template mMaxF<Row>(0,q1); m1++ )
            {
                this->M_model->l2solve( Z1, Fqm[q1][m1], Row );
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->template mMaxF<Row>(0,q2); m2++ )
                    {
                        this->M_model->l2solve( Z2, Fqm[q2][m2], Row );
                        M_R_RhsRhs[Row][q1][m1][q2][m2] = this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop
                } // q2 loop
            } // m1 loop
        } // q1 loop
    } //N==M_Nm


    // LHS0 LOOP ON I
    for ( int i=N0-Nadded0; i<N0; i++ )
    {
        *X0 = this->M_model->template rBFunctionSpace<0>()->primalBasisElement( i );
        for ( int q1=0; q1<QLhs0; q1++ )
        {
            for ( int m1=0; m1<this->M_model->template mMaxA<Row,0>(q1); m1++ )
            {
                Lhs0[q1][m1]->multVector( X0, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );

                // RHS LOOP
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->template mMaxF<Row>(0,q2); m2++ )
                    {
                        M_R_Lhs0Rhs[Row][q1][m1][q2][m2].conservativeResize(N0);
                        this->M_model->l2solve( Z2, Fqm[q2][m2], Row );
                        M_R_Lhs0Rhs[Row][q1][m1][q2][m2]( i ) = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop Rhs
                } // q2 loop Rhs

                // LHS0 LOOP ON J
                for ( int j=0; j<N0; j++ )
                {
                    *Y0 = this->M_model->template rBFunctionSpace<0>()->primalBasisElement(j);
                    if ( optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                            {
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].conservativeResize( N0, N0 );
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1].conservativeResize( N0, N0 );
                                Lhs0[q2][m2]->multVector( Y0, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0

                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1].conservativeResize( N0, N0 );
                        Lhs0[q1][m1]->multVector( Y0, W );
                        W->scale(-1.);
                        this->M_model->l2solve( Z2, W );
                        double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1]( i, j ) = prod;
                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1]( j, i ) = prod;
                    } //optimize
                    else
                    {
                        for ( int q2=0; q2<QLhs0; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                            {
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].conservativeResize( N0, N0 );
                                Lhs0[q2][m2]->multVector( Y0, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j )
                                    = this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    } // !optimize
                } // j loop Lhs0

                // LHS1 LOOP ON J
                for ( int j=0; j<N1; j++ )
                {
                    *Y1 = this->M_model->template rBFunctionSpace<1>()->primalBasisElement(j);
                    for ( int q2=0; q2<QLhs1; q2++ )
                    {
                        for( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                        {
                            M_R_Lhs0Lhs1[Row][q1][m1][q2][m2].conservativeResize( N0, N1 );
                            Lhs1[q2][m2]->multVector( Y1, W );
                            W->scale( -1. );
                            this->M_model->l2solve( Z2, W, Row );
                            M_R_Lhs0Lhs1[Row][q1][m1][q2][m2]( i, j )
                                = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                        } // m2 loop Lhs1
                    } // q2 loop Lhs1
                } // j loop Lhs1

            } // m1 loop Lhs0
        } // q1 loop Lhs0
    } // i loop Lhs0


    // LHS0 LOOP ON J
    for ( int j=N0-Nadded0; j<N0; j++ )
    {
        *X0 = this->M_model->template rBFunctionSpace<0>()->primalBasisElement( j );
        for ( int q1=0; q1<QLhs0; q1++ )
        {
            for ( int m1=0; m1<this->M_model->template mMaxA<Row,0>(q1); m1++ )
            {
                Lhs0[q1][m1]->multVector( X0, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );

                // LHS0 LOOP ON I
                for ( int i=0; i<N0; i++ )
                {
                    *Y0 = this->M_model->template rBFunctionSpace<0>()->primalBasisElement(i);
                    if ( optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                            {
                                Lhs0[q2][m2]->multVector( Y0, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    } //optimize
                    else
                    {
                        for ( int q2=0; q2<QLhs0; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                            {
                                Lhs0[q2][m2]->multVector( Y0, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j )
                                    = this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    } // !optimize
                } // i loop Lhs0

            } // m1 loop Lhs0
        } // q1 loop Lhs0
    } // j loop Lhs0


    // LHS1 LOOP ON I
    for ( int i=N1-Nadded1; i<N1; i++ )
    {
        *X1 = this->M_model->template rBFunctionSpace<1>()->primalBasisElement( i );
        for ( int q1=0; q1<QLhs1; q1++ )
        {
            for ( int m1=0; m1<this->M_model->template mMaxA<Row,1>(q1); m1++ )
            {
                Lhs1[q1][m1]->multVector( X1, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );

                // RHS LOOP
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->template mMaxF<Row>(0,q2); m2++ )
                    {
                        M_R_Lhs1Rhs[Row][q1][m1][q2][m2].conservativeResize(N1);
                        this->M_model->l2solve( Z2, Fqm[q2][m2], Row );
                        M_R_Lhs1Rhs[Row][q1][m1][q2][m2]( i ) = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop Rhs
                } // q2 loop Rhs

                // LHS1 LOOP ON J
                for ( int j=0; j<N1; j++ )
                {
                    *Y1 = this->M_model->template rBFunctionSpace<1>()->primalBasisElement(j);
                    if ( optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                            {
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].conservativeResize( N1, N1 );
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1].conservativeResize( N1, N1 );
                                Lhs1[q2][m2]->multVector( Y1, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0

                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1].conservativeResize( N1, N1 );
                        Lhs1[q1][m1]->multVector( Y1, W );
                        W->scale(-1.);
                        this->M_model->l2solve( Z2, W );
                        double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( i, j ) = prod;
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( j, i ) = prod;
                    } //optimize
                    else
                    {
                        for ( int q2=0; q2<QLhs1; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                            {
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].conservativeResize( N1, N1 );
                                Lhs1[q2][m2]->multVector( Y1, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j )
                                    = this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1
                    } // !optimize
                } // j loop Lhs1

            } // m1 loop Lhs1
        } // q1 loop Lhs1
    } // i loop Lhs1


    // LHS1 LOOP ON J
    for ( int j=N1-Nadded1; j<N1; j++ )
    {
        *X1 = this->M_model->template rBFunctionSpace<1>()->primalBasisElement( j );
        for ( int q1=0; q1<QLhs1; q1++ )
        {
            for ( int m1=0; m1<this->M_model->template mMaxA<Row,1>(q1); m1++ )
            {
                Lhs1[q1][m1]->multVector( X1, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );

                // LHS0 LOOP ON I
                for ( int i=0; i<N0; i++ )
                {
                    *Y0 = this->M_model->template rBFunctionSpace<0>()->primalBasisElement( i );
                    for ( int q2=0; q2<QLhs0; q2++ )
                    {
                        for ( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                        {
                            Lhs0[q2][m2]->multVector( Y0, W );
                            W->scale(-1.);
                            this->M_model->l2solve( Z2, W, Row );
                            M_R_Lhs0Lhs1[Row][q2][m2][q1][m1]( i, j )
                                = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                        }
                    } // q1 loop Lhs0
                } // i loop Lhs0

                // LHS1 LOOP ON I
                for ( int i=0; i<N1; i++ )
                {
                    *Y1 = this->M_model->template rBFunctionSpace<1>()->primalBasisElement(i);
                    if ( optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                            {
                                Lhs1[q2][m2]->multVector( Y1, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1

                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1].conservativeResize( N1, N1 );
                        Lhs1[q1][m1]->multVector( Y1, W );
                        W->scale(-1.);
                        this->M_model->l2solve( Z2, W );
                        double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( i, j ) = prod;
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( j, i ) = prod;
                    } //optimize
                    else
                    {
                        for ( int q2=0; q2<QLhs1; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                            {
                                Lhs1[q2][m2]->multVector( Y1, W );
                                W->scale( -1. );
                                this->M_model->l2solve( Z2, W, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j )
                                    = this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1
                    } // !optimize
                } // j loop Lhs1

            } // m1 loop Lhs1
        } // q1 loop Lhs1
    } // j loop Lhs1
} // offlineResidual

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::max_error_type
CRBSaddlePoint<TruthModelType>::maxResidual( int N )
{
    std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
    std::vector< vectorN_type > uNduold;

    double rez=0;
    parameter_type mu;

    // we evaluate residual for each parameter in the complement of W_MNmu
    for ( int k=0; k<this->M_WNmu_complement->size(); k++ )
    {
        parameter_type const& current_mu = this->M_WNmu_complement->at(k);
        lb( N, current_mu, uN, uNdu, uNold, uNduold );
        double current_rez = onlineResidual( N, current_mu, uN[0] );
        if ( current_rez>rez )
        {
            mu = current_mu;
            rez = current_rez;
        }
    } // loop on M_WNmu_complement

    // we find the proc which has the max residual
    int world_size = Environment::worldComm().globalSize();
    std::vector<double> max_world( world_size );
    mpi::all_gather( Environment::worldComm().globalComm(),
                     rez,
                     max_world );
    auto it_max = std::max_element( max_world.begin(), max_world.end() );
    int proc_having_good_mu = it_max - max_world.begin();

    // we broadcast the good parameter and the value of the max residual
    auto tuple = boost::make_tuple( mu, rez );
    boost::mpi::broadcast( Environment::worldComm(), tuple, proc_having_good_mu );
    mu = tuple.template get<0>();
    rez = tuple.template get<1>();


    CRB_COUT << std::setprecision(15) << "[CRBSaddlePoint] max residual="<< rez << std::endl;

    return boost::make_tuple( rez, mu, 0, 0 );
}

template<typename TruthModelType>
double
CRBSaddlePoint<TruthModelType>::onlineResidual( int Ncur, parameter_type const& mu,
                                                vectorN_type Un ) const
{
    double res0 = onlineResidual<0>( Ncur, mu, Un );
    double res1 = onlineResidual<1>( Ncur, mu, Un );

    return std::sqrt( res0 + res1 );
}


template<typename TruthModelType>
template<int Row>
double
CRBSaddlePoint<TruthModelType>::onlineResidual( int Ncur, parameter_type const& mu,
                                                vectorN_type Un, bool test ) const
{
    int N0 = M_N0[Ncur];
    int N1 = M_N1[Ncur];
    int N [2] = { M_N0[Ncur], M_N1[Ncur] };

    CHECK( Un.size() == N0 + N1 )
        << "invalide size of Un, vector can't be cut, Un.size="
        <<Un.size()<<", N0="<<N0<<", N1="<<N1<<std::endl;

    vectorN_type un = Un.head( N0 );
    vectorN_type pn = Un.tail( N1 );

    int QRhs = this->M_model->template Ql<Row>( 0 );
    int QLhs0 = this->M_model->template Qa<Row,0>();
    int QLhs1 = this->M_model->template Qa<Row,1>();

    beta_vector_type betaLhs0 = this->M_model->template computeBetaAqm<Row,0>(mu);
    beta_vector_type betaLhs1 = this->M_model->template computeBetaAqm<Row,1>(mu);;
    beta_vector_type betaRhs;
    if ( QRhs )
        betaRhs = this->M_model->template computeBetaFqm<Row>(0,mu);


    double RhsRhs = 0;
    double Lhs0Rhs = 0;
    double Lhs0Lhs0 = 0;
    double Lhs0Lhs1 =0;
    double Lhs1Rhs = 0;
    double Lhs1Lhs1 = 0;

    // LOOP ON RHS
    for ( int q1=0; q1<QRhs; q1++ )
    {
        for ( int m1=0; m1<this->M_model->template mMaxF<Row>( 0, q1 ); m1++ )
        {
            double beta_q1 = betaRhs[q1][m1];
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxF<Row>( 0, q2 ); m2++ )
                {
                    double beta_q2 = betaRhs[q2][m2];
                    RhsRhs += M_R_RhsRhs[Row][q1][m1][q2][m2]*beta_q1*beta_q2;
                }
            } // q2 loop rhs
        } // m1 loop rhs
    } // q1 loop rhs

    // LOOP ON LHS0
    for ( int q1=0; q1<QLhs0; q1++ )
    {
        for ( int m1=0; m1<this->M_model->template mMaxA<Row,0>(q1); m1++ )
        {
            double beta_q1 = betaLhs0[q1][m1];

            // SUBLOOP ON RHS
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxF<Row>( 0, q2 ); m2++ )
                {
                    double beta_q2 = betaRhs[q2][m2];
                    Lhs0Rhs += beta_q1*beta_q2*M_R_Lhs0Rhs[Row][q1][m1][q2][m2].head(N0).dot(un);
                } // m2 loop rhs
            } // q2 loop rhs

            // SUBLOOP ON LHS0
            for ( int q2=0; q2<QLhs0; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxA<Row,0>(q2); m2++ )
                {
                    double beta_q2 = betaLhs0[q2][m2];

                    auto m = M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].block(0,0,N0,N0)*un;
                    Lhs0Lhs0 += beta_q1*beta_q2*un.dot(m);
                } // m2 loop lhs0
            } // q2 loop lhs0

            // SUBLOOP ON LHS1
            for ( int q2=0; q2<QLhs1; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                {
                    double beta_q2 = betaLhs1[q2][m2];
                    auto m = M_R_Lhs0Lhs1[Row][q1][m1][q2][m2].block(0,0,N0,N1)*pn;
                    Lhs0Lhs1 += beta_q1*beta_q2*un.dot(m);
                } // m2 loop on lhs1
            } // q2 loop on lhs1

        } // m1 loop lhs0
    } // q1 loop lhs0

    // LOOP ON LHS1
    for ( int q1=0; q1<QLhs1; q1++ )
    {
        for ( int m1=0; m1<this->M_model->template mMaxA<Row,1>(q1); m1++ )
        {
            double beta_q1 = betaLhs1[q1][m1];

            // SUBLOOP ON RHS
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxF<Row>( 0, q2 ); m2++ )
                {
                    double beta_q2 = betaRhs[q2][m2];
                    Lhs1Rhs += beta_q1*beta_q2*M_R_Lhs1Rhs[Row][q1][m1][q2][m2].head(N1).dot(pn);
                } // m2 loop rhs
            } // q2 loop rhs

            // SUBLOOP ON LHS1
            for ( int q2=0; q2<QLhs1; q2++ )
            {
                for ( int m2=0; m2<this->M_model->template mMaxA<Row,1>(q2); m2++ )
                {
                    double beta_q2 = betaLhs1[q2][m2];
                    auto m = M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].block(0,0,N1,N1)*pn;
                    Lhs1Lhs1 += beta_q1*beta_q2*pn.dot(m);
                } // m2 loop on lhs1
            } // q2 loop on lhs1

        } // m1 loop lhs1
    } // q1 loop lhs1

    if ( test )
    {
        CRB_COUT << "test online Residual : \n RhsRhs="<<RhsRhs << ", Lhs0Rhs="<< Lhs0Rhs
                 <<", Lhs0Lhs0="<< Lhs0Lhs0 <<", Lhs0Lhs1="<<Lhs0Lhs1
                 <<", Lhs1Lhs1="<< Lhs1Lhs1 <<", Lhs1Rhs="<< Lhs1Rhs <<std::endl;

        vector_ptrtype F = this->M_model->template newL2Vector<Row>();
        vector_ptrtype F0 = this->M_model->template newL2Vector<Row>();
        vector_ptrtype F1 = this->M_model->template newL2Vector<Row>();
        vector_ptrtype Un = this->M_model->template newL2Vector<0>();
        vector_ptrtype Pn = this->M_model->template newL2Vector<1>();
        sparse_matrix_ptrtype A0 = this->M_model->template newL2Matrix<Row,0>();
        sparse_matrix_ptrtype A1 = this->M_model->template newL2Matrix<Row,1>();
        vector_ptrtype Rhs = this->M_model->template newL2Vector<Row>();
        vector_ptrtype Lhs0 = this->M_model->template newL2Vector<Row>();
        vector_ptrtype Lhs1 = this->M_model->template newL2Vector<Row>();

        auto truc0 = un(0)*this->M_model->template rBFunctionSpace<0>()->primalBasisElement(0);
        for ( int k=1; k<un.size(); k++ )
        {
            truc0 += un(k)*this->M_model->template rBFunctionSpace<0>()->primalBasisElement(k);
        }
        *Un = truc0;

        auto truc1 = pn(0)*this->M_model->template rBFunctionSpace<1>()->primalBasisElement(0);
        for ( int k=1; k<un.size(); k++ )
        {
            truc1 += pn(k)*this->M_model->template rBFunctionSpace<1>()->primalBasisElement(k);
        }
        *Pn = truc1;


        auto blockF = this->M_model->M_F[0].template get<Row>();
        if ( blockF )
        {
            auto beta = blockF->computeBetaQm( mu );
            auto Fqm = blockF->compute();
            for ( int q=0; q<blockF->Q(); q++ )
                for( int m=0; m<blockF->mMax(q); m++ )
                    F->add( beta[q][m], Fqm[q][m]);
        }
        this->M_model->l2solve( Rhs, F, Row );

        auto blockA0 = this->M_model->M_A.template get<Row,0>();
        if ( blockA0 )
        {
            auto beta = blockA0->computeBetaQm(mu);
            auto Aqm = blockA0->compute();
            for ( int q=0; q<blockA0->Q(); q++ )
                for( int m=0; m<blockA0->mMax(q); m++ )
                    A0->addMatrix( beta[q][m], Aqm[q][m]);
            A0->multVector( Un, F0 );
            F0->scale(-1);
        }
        else if ( boption("crb.saddlepoint.transpose") && this->M_model->M_A.template get<0,Row>() )
        {
            auto beta = this->M_model->M_A.template get<0,Row>()->computeBetaQm(mu);
            auto Aqm = this->M_model->M_A.template get<0,Row>()->compute();
            sparse_matrix_ptrtype truc= this->M_model->template newL2Matrix<0,Row>();
            for ( int q=0; q<this->M_model->M_A.template get<0,Row>()->Q(); q++ )
                for( int m=0; m<this->M_model->M_A.template get<0,Row>()->mMax(q); m++ )
                    truc->addMatrix( beta[q][m], Aqm[q][m]);
            truc->transpose( A0 );
            A0->multVector( Un, F0 );
            F0->scale(-1);
        }
        this->M_model->l2solve( Lhs0, F0, Row );

        auto blockA1 = this->M_model->M_A.template get<Row,1>();
        if ( blockA1 )
        {
            auto beta = blockA1->computeBetaQm(mu);
            auto Aqm = blockA1->compute();
            for ( int q=0; q<blockA1->Q(); q++ )
                for( int m=0; m<blockA1->mMax(q); m++ )
                    A1->addMatrix( beta[q][m], Aqm[q][m]);
            A1->multVector( Pn, F1 );
            F1->scale(-1);
        }
        else if ( boption("crb.saddlepoint.transpose") && this->M_model->M_A.template get<1,Row>() )
        {
            auto beta = this->M_model->M_A.template get<1,Row>()->computeBetaQm(mu);
            auto Aqm = this->M_model->M_A.template get<1,Row>()->compute();
            sparse_matrix_ptrtype truc= this->M_model->template newL2Matrix<1,Row>();
            for ( int q=0; q<this->M_model->M_A.template get<1,Row>()->Q(); q++ )
                for( int m=0; m<this->M_model->M_A.template get<1,Row>()->mMax(q); m++ )
                    truc->addMatrix( beta[q][m], Aqm[q][m]);
            truc->transpose( A1 );
            A1->multVector( Pn, F1 );
            F1->scale(-1);
        }
        this->M_model->l2solve( Lhs1, F1, Row );


        auto RhsRhs_test = this->M_model->scalarProduct(Rhs,Rhs, Row);
        auto Lhs0Rhs_test = 2.*this->M_model->scalarProduct(Lhs0, Rhs, Row);
        auto Lhs0Lhs0_test = this->M_model->scalarProduct(Lhs0, Lhs0, Row);
        auto Lhs0Lhs1_test = 2.*this->M_model->scalarProduct(Lhs0, Lhs1, Row);
        auto Lhs1Rhs_test = 2.*this->M_model->scalarProduct(Lhs1, Rhs, Row);
        auto Lhs1Lhs1_test = this->M_model->scalarProduct(Lhs1, Lhs1, Row);

        CRB_COUT << "test block residual : \n RhsRhs="<<RhsRhs_test << ", Lhs0Rhs="<< Lhs0Rhs_test
                 <<", Lhs0Lhs0="<< Lhs0Lhs0_test <<", Lhs0Lhs1="<<Lhs0Lhs1_test
                 <<", Lhs1Lhs1="<< Lhs1Lhs1_test <<", Lhs1Rhs="<< Lhs1Rhs_test <<std::endl;

    }

    double delta = math:: abs( RhsRhs + Lhs0Rhs + Lhs0Lhs0 + Lhs0Lhs1 + Lhs1Rhs + Lhs1Lhs1 );

    return delta;
} // onlineResidual


/**
 * Initialize the vectors to stock data of the residual evaluation
 */
template<typename TruthModelType>
template<int Row>
void
CRBSaddlePoint<TruthModelType>::initResidualVectors()
{
    int QRhs = this->M_model->template Ql<Row>( 0 );
    int QLhs0 = this->M_model->template Qa<Row,0>();
    int QLhs1 = this->M_model->template Qa<Row,1>();

    if ( M_R_RhsRhs.size()==0 )
    {
        M_R_RhsRhs.resize(2);
        M_R_Lhs0Rhs.resize(2);
        M_R_Lhs1Rhs.resize(2);
        M_R_Lhs0Lhs0.resize(2);
        M_R_Lhs0Lhs1.resize(2);
        M_R_Lhs1Lhs1.resize(2);
    }

    M_R_RhsRhs[Row].resize(QRhs);
    for ( int q1=0; q1<QRhs; q1++ )
    {
        int mMax1 = this->M_model->template mMaxF<Row>(0,q1);
        M_R_RhsRhs[Row][q1].resize(mMax1);
        for ( int m1=0; m1<mMax1; m1++ )
        {
            M_R_RhsRhs[Row][q1][m1].resize(QRhs);
            for ( int q2=0; q2<QRhs; q2++ )
                M_R_RhsRhs[Row][q1][m1][q2].resize( this->M_model->template mMaxF<Row>(0,q2) );
        }
    }

    M_R_Lhs0Lhs0[Row].resize(QLhs0);
    M_R_Lhs0Lhs1[Row].resize(QLhs0);
    M_R_Lhs0Rhs[Row].resize(QLhs0);
    for ( int q1=0; q1<QLhs0; q1++ )
    {
        int mMax1 = this->M_model->template mMaxA<Row,0>(q1);
        M_R_Lhs0Lhs0[Row][q1].resize(mMax1);
        M_R_Lhs0Lhs1[Row][q1].resize(mMax1);
        M_R_Lhs0Rhs[Row][q1].resize(mMax1);
        for ( int m1=0; m1<mMax1; m1++ )
        {
            M_R_Lhs0Lhs0[Row][q1][m1].resize(QLhs0);
            M_R_Lhs0Lhs1[Row][q1][m1].resize(QLhs1);
            M_R_Lhs0Rhs[Row][q1][m1].resize(QRhs);
            for ( int q2=0; q2<QLhs0; q2++ )
                M_R_Lhs0Lhs0[Row][q1][m1][q2].resize( this->M_model->template mMaxA<Row,0>(q2) );
            for ( int q2=0; q2<QLhs1; q2++ )
                M_R_Lhs0Lhs1[Row][q1][m1][q2].resize( this->M_model->template mMaxA<Row,1>(q2) );
            for ( int q2=0; q2<QRhs; q2++ )
                M_R_Lhs0Rhs[Row][q1][m1][q2].resize( this->M_model->template mMaxF<Row>(0,q2) );
        }
    }

    M_R_Lhs1Lhs1[Row].resize(QLhs1);
    M_R_Lhs1Rhs[Row].resize(QLhs1);
    for ( int q1=0; q1<QLhs1; q1++ )
    {
        int mMax1 = this->M_model->template mMaxA<Row,1>(q1);
        M_R_Lhs1Lhs1[Row][q1].resize(mMax1);
        M_R_Lhs1Rhs[Row][q1].resize(mMax1);

        for ( int m1=0; m1<mMax1; m1++ )
        {
            M_R_Lhs1Lhs1[Row][q1][m1].resize(QLhs1);
            M_R_Lhs1Rhs[Row][q1][m1].resize(QRhs);
            for ( int q2=0; q2<QLhs1; q2++ )
                M_R_Lhs1Lhs1[Row][q1][m1][q2].resize( this->M_model->template mMaxA<Row,1>(q2) );
            for ( int q2=0; q2<QRhs; q2++ )
                M_R_Lhs1Rhs[Row][q1][m1][q2].resize( this->M_model->template mMaxF<Row>(0,q2) );
        }
    }
}


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::exportBasisFunctions()
{
    auto u_vec = XN0->primalRB();
    auto p_vec = XN1->primalRB();

    auto e = exporter( _mesh=this->M_model->functionSpace()->mesh(), _name="basis-functions" );
    CHECK( u_vec.size()>0) << "[CRBSaddlePoint::exportBasisFunctions] Error : there are no element to export in first rb\n";
    CHECK( p_vec.size()>0) << "[CRBSaddlePoint::exportBasisFunctions] Error : there are no element to export in second rb\n";

    for( int index=0; index<this->M_N; index++)
    {
        std::string mu_str;
        auto mu=this->M_WNmu->at( index );

        for ( int i=0; i<mu.size(); i++)
            mu_str += ( boost::format( "_%1%" ) %mu[i] ).str() ;

        int u_count=0;
        for( int u_index=M_N0[index]; u_index<M_N0[index+1]; u_index++ )
        {
            u_count++;
            std::string basis_name = ( boost::format( "u_pr_%1%.%2%_param") %index %u_count ).str();
            std::string name = basis_name + mu_str;
            e->step(0)->add( name, u_vec[u_index] );
        }

        for( int p_index=M_N1[index]; p_index<M_N1[index+1]; p_index++ )
        {
            std::string basis_name = ( boost::format( "p_pr_%1%_param") %index ).str();
            std::string name = basis_name + mu_str;
            e->step(0)->add( name, p_vec[p_index] );
        }
    }
    e->save();
}


template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansionPrimal( vectorN_type const& U_coeff, int const N ) const
{
    auto WN0 = XN0->primalRB();
    auto WN1 = XN1->primalRB();
    return expansionSaddlePoint( U_coeff, N, WN0, WN1 );
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansionDual( vectorN_type const& U_coeff, int const N ) const
{
    auto WN0 = XN0->dualRB();
    auto WN1 = XN1->dualRB();
    return expansionSaddlePoint( U_coeff, N, WN0, WN1 );
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansionSaddlePoint( vectorN_type const& U_coeff, int const N,
                                                      std::vector<subelement_type<0>> WN0,
                                                      std::vector<subelement_type<1>> WN1 ) const
{
    int Nwn = N>0 ? N:this->M_N;

    CHECK( Nwn <= WN0.size() )<< "invalid expansion size\n";
    CHECK( M_N0[Nwn]+M_N1[Nwn] <= U_coeff.size() )<< "invalid expansion size, Nwn="
                                                  << M_N0[Nwn]+M_N1[Nwn]
                                                  << ", U_coeff.size="<<U_coeff.size()<<std::endl;
    CHECK( U_coeff.size() == M_N0[this->M_N] + M_N1[this->M_N] )
        << "invalide size of U_coeff, vector can't be cut\n";

    vectorN_type u_coeff = U_coeff.head( M_N0[this->M_N] );
    vectorN_type p_coeff = U_coeff.tail( M_N1[this->M_N] );

    element_type U = this->M_model->functionSpace()->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    u = Feel::expansion( WN0, u_coeff, M_N0[Nwn] ).container();
    p = Feel::expansion( WN1, p_coeff, M_N1[Nwn] ).container();

    return U;
}


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::generateSuperSampling()
{

    int proc_number = worldComm().globalRank();
    int total_proc = worldComm().globalSize();
    bool all_proc_same_sampling = boption( "crb.all-procs-have-same-sampling" );
    int sampling_size = ioption("crb.sampling-size");
    std::string sampling_mode = soption("crb.sampling-mode");
    std::string file_name = ( boost::format("M_Xi_%1%_"+sampling_mode+"-proc%2%on%3%")
                              %sampling_size
                              %proc_number
                              %total_proc ).str();

    if ( all_proc_same_sampling )
        file_name += "-all_proc_same_sampling";

    std::ifstream file ( file_name );

    if ( !file )
    {
        std::string supersamplingname =(boost::format("Dmu-%1%-generated-by-master-proc") %sampling_size ).str();

        if( sampling_mode == "log-random" )
            this->M_Xi->randomize( sampling_size , all_proc_same_sampling , supersamplingname );
        else if( sampling_mode == "log-equidistribute" )
            this->M_Xi->logEquidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else if( sampling_mode == "equidistribute" )
            this->M_Xi->equidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else
            throw std::logic_error( "[CRBSaddlePoint::offline] ERROR invalid option crb.sampling-mode, please select between log-random, log-equidistribute or equidistribute" );

        this->M_Xi->writeOnFile(file_name);
    }
    else
    {
        this->M_Xi->clear();
        this->M_Xi->readFromFile(file_name);
    }

    this->M_WNmu->setSuperSampling( this->M_Xi );
} //generateSuperSampling()


template<typename TruthModelType>
bool
CRBSaddlePoint<TruthModelType>::buildSampling()
{
    bool use_predefined_WNmu = boption("crb.use-predefined-WNmu");
    int N_log_equi = ioption("crb.use-logEquidistributed-WNmu");
    int N_equi = ioption("crb.use-equidistributed-WNmu");
    int N_random = ioption( "crb.use-random-WNmu" );

    std::string file_name = ( boost::format("SamplingWNmu") ).str();
    std::ifstream file ( file_name );
    this->M_WNmu->clear();

    if ( use_predefined_WNmu ) // In this case we want to read the sampling
    {
        if( ! file ) // The user forgot to give the sampling file
            throw std::logic_error( "[CRB::offline] ERROR the file SamplingWNmu doesn't exist so it's impossible to known which parameters you want to use to build the database" );
        else
        {
            int sampling_size = this->M_WNmu->readFromFile(file_name);
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Read WNmu ( sampling size : "
                         << sampling_size <<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Read WNmu ( sampling size : "
                       << sampling_size <<" )";
        }
    }
    else if ( this->M_error_type==CRB_NO_RESIDUAL )// We generate the sampling with choosen strategy
    {
        if ( N_log_equi>0 )
        {
            this->M_WNmu->logEquidistribute( N_log_equi , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Log-Equidistribute WNmu ( sampling size : "
                         <<N_log_equi<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Log-Equidistribute WNmu ( sampling size : "
                       <<N_log_equi<<" )";
        }
        else if ( N_equi>0 )
        {
            this->M_WNmu->equidistribute( N_equi , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Equidistribute WNmu ( sampling size : "
                         <<N_equi<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Equidistribute WNmu ( sampling size : "
                       <<N_equi<<" )";
        }
        else if ( N_random>0 )
        {
            this->M_WNmu->randomize( N_random , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Randomize WNmu ( sampling size : "
                         <<N_random<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Randomize WNmu ( sampling size : "
                       <<N_random<<" )";
        }
        else // In this case we don't know what sampling to use
            throw std::logic_error( "[CRB::offline] ERROR : You have to choose an appropriate strategy for the offline sampling : random, equi, logequi or predefined" );

        this->M_WNmu->writeOnFile(file_name);
        use_predefined_WNmu=true;
    } //build sampling

    return use_predefined_WNmu;
} //buildSampling()



template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    ar & boost::serialization::base_object<super_crb>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_N0 );
    ar & BOOST_SERIALIZATION_NVP( M_N1 );
}


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::load( Archive & ar, const unsigned int version )
{
    ar & boost::serialization::base_object<super_crb>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_N0 );
    ar & BOOST_SERIALIZATION_NVP( M_N1 );
}


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveDB()
{
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );

    if ( ofs )
    {
        boost::archive::text_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
} //saveDB()

template<typename TruthModelType>
bool
CRBSaddlePoint<TruthModelType>::loadDB()
{
    if( this->isDBLoaded() )
        return true;

    fs::path db = this->lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    fs::ifstream ifs( db );

    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        this->setIsLoaded( true );
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
} // loadDB()

} // namespace Feel


namespace boost
{
namespace serialization
{
template< typename T>
struct version< Feel::CRBSaddlePoint<T> >
{
    // at the moment the version of the CRBTrilinear DB is 0. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<0> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRBSaddlePoint<T> >::value;

} // namespace serialization
} // namespace boost



#endif // CRBSADDLEPOINT_H
