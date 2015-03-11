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

                for( size_type q=0; q<M_crb->M_model->template Qa<row::value,col::value>(); q++ )
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
                int Qa = M_crb->M_model-> template Qa<row::value,col::value>();
                int N = M_crb->M_N0[M_N] + M_crb->M_N1[M_crb->M_N];

                if ( Qa!=0 )
                {
                    auto betaA = M_crb->M_model->template computeBetaAqm<row::value,col::value>( M_mu );
                    for( int q=0; q<Qa; q++ )
                    {
                        int mMax = M_crb->M_model->template mMaxA<row::value,col::value>(q);
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
    boost::timer ti;
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
        ti.restart();
        // Generate sampling depending of the "crb.sampling-mode" option
        generateSuperSampling();
        LOG(INFO) << " -- super sampling init done in " << ti.elapsed() << "s";

        ti.restart();

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

        M_blockAqm_pr[0][0].resize( this->M_model->template Qa<0,0>() );
        for ( int q=0; q<M_blockAqm_pr[0][0].size(); q++ )
            M_blockAqm_pr[0][0][q].resize( this->M_model->template mMaxA<0,0>(q) );

        M_blockAqm_pr[0][1].resize( this->M_model->template Qa<0,1>() );
        for ( int q=0; q<M_blockAqm_pr[0][1].size(); q++ )
            M_blockAqm_pr[0][1][q].resize( this->M_model->template mMaxA<0,1>(q) );

        M_blockAqm_pr[1][0].resize( this->M_model->template Qa<1,0>() );
        for ( int q=0; q<M_blockAqm_pr[1][0].size(); q++ )
            M_blockAqm_pr[1][0][q].resize( this->M_model->template mMaxA<1,0>(q) );

        M_blockAqm_pr[1][1].resize( this->M_model->template Qa<1,1>() );
        for ( int q=0; q<M_blockAqm_pr[1][1].size(); q++ )
            M_blockAqm_pr[1][1][q].resize( this->M_model->template mMaxA<1,1>(q) );

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
        this->M_iter_max = this->M_WNmu->size();
        mu = this->M_WNmu->at( std::min( this->M_N, this->M_iter_max-1 ) ); // first element
    }

    if( this->M_error_type == CRB_NO_RESIDUAL || use_predefined_WNmu )
    {
        //in this case it makes no sens to check the estimated error
        this->M_maxerror = 1e10;
    }


    while( this->M_maxerror>this->M_tolerance && this->M_N<this->M_iter_max )
    {
        boost::timer timer;
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

        timer.restart();
        *U = this->M_model->solve( mu );
        if ( Environment::isMasterRank() )
            std::cout << "  -- primal problem solved in "<< timer.elapsed() << std::endl;
        timer.restart();

        if ( !use_predefined_WNmu)
            this->M_WNmu->push_back( mu, index );
        this->M_WNmu_complement = this->M_WNmu->complement();

        auto u = U->template elementPtr<0>();
        auto p = U->template elementPtr<1>();

        M_N0.push_back( M_N0[this->M_N] );
        M_N1.push_back( M_N1[this->M_N] );
        this->M_N++;

        XN1->addPrimalBasisElement( p );
        XN1->addDualBasisElement( p );
        M_Nadded1=1;
        M_N1[this->M_N] += M_Nadded1;

        this->orthonormalize( M_N1[this->M_N], XN1->primalRB(), M_Nadded1, 1,
                              ioption("crb.saddlepoint.orthonormalize1") );

        XN0->addPrimalBasisElement( u );
        XN0->addDualBasisElement( u );
        M_Nadded0=1;

        if ( boption("crb.saddlepoint.add-false-supremizer") )
        {
            timer.restart();
            auto us = this->M_model->falseSupremizer( mu, p );
            XN0->addPrimalBasisElement( us );
            XN0->addDualBasisElement( us );
            M_Nadded0++;
            CRB_COUT << "  -- false supremizer added in "<< timer.elapsed() << std::endl;
        }
        if ( boption("crb.saddlepoint.add-supremizer") )
        {
            timer.restart();
            //auto us = this->M_model->supremizer( mu, p );
            auto us = this->M_model->supremizer( mu, XN1->primalRB().back() );
            XN0->addPrimalBasisElement( us );
            XN0->addDualBasisElement( us );
            M_Nadded0++;
            CRB_COUT << "  -- supremizer added in "<< timer.elapsed() << std::endl;
        }

        M_N0[this->M_N] += M_Nadded0;
        this->orthonormalize( M_N0[this->M_N], XN0->primalRB(), M_Nadded0, 0,
                              ioption("crb.saddlepoint.orthonormalize0") );


        matrixblockrange_type matrixrange;
        ComputeADElements compute_elements( this->shared_from_this() );
        fusion::for_each( matrixrange, compute_elements );


        if ( ! use_predefined_WNmu )
        {
            bool already_exist;
            do
            {
                //initialization
                already_exist=false;
                //pick randomly an element
                mu = this->M_Dmu->element();
                //make sure that the new mu is not already is M_WNmu
                BOOST_FOREACH( auto _mu, *this->M_WNmu )
                {
                    if( mu == _mu )
                        already_exist=true;
                }
            }
            while( already_exist );
            this->M_current_mu = mu;
        }
        else
        {
            //remmber that in this case M_iter_max = sampling size
            if( this->M_N < this->M_iter_max )
            {
                mu = this->M_WNmu->at( this->M_N );
                this->M_current_mu = mu;
            }
        }

        this->M_rbconv.insert( convergence(this->M_N, boost::make_tuple(this->M_maxerror,
                                                                        delta_pr,delta_du) ) );

        saveDB();
        M_elements_database0.setWn( boost::make_tuple(XN0->primalRB(), XN0->dualRB()) );
        M_elements_database1.setWn( boost::make_tuple(XN1->primalRB(), XN1->dualRB()) );
        M_elements_database0.saveDB();
        M_elements_database1.saveDB();

        CRB_COUT << "============================================================\n";
    } // while( this->M_maxerror>this->M_tolerance && this->M_N<this->M_iter_max )

    if ( boption("crb.visualize-basis") )
    {
        exportBasisFunctions();
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
        int Qa01 = this->M_model-> template Qa<0,1>();
        int Qa10 = this->M_model-> template Qa<1,0>();
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
    else // We generate the sampling with choosen strategy
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
    } //build sampling

    return true;
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
