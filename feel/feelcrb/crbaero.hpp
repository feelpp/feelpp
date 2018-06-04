#ifndef CRBAERO_H
#define CRBAERO_H

#include <feel/feelcrb/crbblock.hpp>
#define POUT std::cout << "[" << Environment::worldComm().globalRank()<<"] "


namespace Feel
{
po::options_description
crbAeroOptions( std::string const& prefix="" );

template <typename TruthModelType>
class CRBAero :
        public CRBBlock<TruthModelType>
{
    typedef CRBBlock<TruthModelType> super;
    typedef CRBAero<TruthModelType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

public:
    typedef TruthModelType model_type;
    typedef typename super::blockmatrixN_type blockmatrixN_type;
    typedef typename super::truth_model_ptrtype truth_model_ptrtype;
    typedef typename super::element_type element_type;
    typedef typename super::parameter_type parameter_type;
    typedef typename super::element_ptrtype element_ptrtype;
    typedef typename super::matrix_info_tuple matrix_info_tuple;
    typedef typename super::map_dense_vector_type map_dense_vector_type;
    typedef typename super::map_dense_matrix_type map_dense_matrix_type;
    typedef typename super::space_type space_type;
    typedef typename super::beta_vector_type beta_vector_type;
    typedef typename super::matrixN_type matrixN_type;
    typedef typename super::max_error_type max_error_type;
    typedef typename super::error_estimation_type error_estimation_type;


    static self_ptrtype New( std::string const& name = "defaultname_crb",
                             crb::stage stage = crb::stage::online )
        {
            return New( name, boost::make_shared<model_type>(stage), stage );
        }

    static self_ptrtype New( std::string const& name,
                             truth_model_ptrtype const& model,
                             crb::stage stage = crb::stage::online,
                             std::string const& prefixElt = "")
        {
            auto crb = boost::shared_ptr<self_type>( new self_type(name, model, stage, prefixElt ));
            crb->init();
            return crb;
        }


protected:
    //! Default Constructor
    CRBAero( std::string const& name = "defaultname_crb",
              crb::stage stage = crb::stage::online,
              WorldComm const& worldComm = Environment::worldComm() ) :
        CRBAero( name, boost::make_shared<model_type>(stage), stage )
        {}

    //! constructor from command line options
    CRBAero( std::string const& name, truth_model_ptrtype const & model,
             crb::stage stage = crb::stage::online, std::string const& prefixExt = "" ) :
        super( name, model, stage, prefixExt )
        {
            CHECK( n_block==3 ) <<"CRBAero is supposed to work with a composite space u,p,T\n";
        }



public:
    //! Init function, specific to crb aero
    void init();
    void offlineSolve( element_type& u, element_type& udu, parameter_type& mu, element_ptrtype & dual_initial_field )override
        { u = this->M_model->offlineSolveAD(mu); }

    element_type solve( parameter_type const& mu )
    {
        return this->M_model->solve( mu );
    }

    matrix_info_tuple onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const override;
    double empiricalError( int N, parameter_type const& mu, std::vector<double> output_vec ) const;
    max_error_type maxErrorBounds( size_type N ) const override;
    error_estimation_type delta( size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN, std::vector< vectorN_type > const& uNdu, std::vector<vectorN_type> const& uNold, std::vector<vectorN_type> const& uNduold, int k=0 ) const override
        {
            std::vector< std::vector<double> > primal_residual_coeffs;
            std::vector< std::vector<double> > dual_residual_coeffs;
            std::vector<double> output_upper_bound;
            double delta_pr=0;
            double delta_du=0;
            primal_residual_coeffs.resize(1);

            output_upper_bound.resize(1);
            output_upper_bound[0]=-1;
            return boost::make_tuple( output_upper_bound ,primal_residual_coeffs,dual_residual_coeffs,delta_pr,delta_du );
        }

private:
    void updateJacobianOnline( const map_dense_vector_type& X, map_dense_matrix_type& J , parameter_type const& mu , int N ) const;
    void updateResidualOnline( const map_dense_vector_type& X, map_dense_vector_type& R , parameter_type const& mu , int N ) const;

    virtual void buildRbMatrixTrilinear( int number_of_added_elements, parameter_type& mu ) override;

    void updateSpecificTerms() override
        {
            int q_tri = this->M_model->QTri();
            if ( q_tri )
            {
                if ( this->M_blockTriqm_pr.size()==0 )
                {
                    this->M_blockTriqm_pr.resize( n_block );
                    for( int i = 0; i < n_block; i++ )
                        this->M_blockTriqm_pr[i].resize(n_block);
                }
                for (int r=0; r<n_block; r++ )
                    for ( int c=0; c<n_block; c++ )
                        this->M_blockTriqm_pr[r][c].resize(q_tri);
            }
        }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version );

    mutable matrixN_type M_Jbil;
    mutable vectorN_type M_Rli;

    using super::n_block;
    using super::subN;
    using super::blockAqm;
    using super::blockFqm;
    using super::blockLqm;
    using super::blockTriqm;
}; // class crbaero



template <typename TruthModelType>
void
CRBAero<TruthModelType>::init()
{
    LOG(INFO) << "CRBaero init begin\n";

    using Feel::cout;

    static_assert( space_type::nSpaces==3, "Invalid space type, CRBAero has to be used with composite space of 3 subspaces\n" );

    if ( !this->M_rebuild && this->loadDB() )
    {
        cout << "Database CRB Aero " << this->lookForDB() << " available and loaded with M_N0="
             << subN(0) << ", M_N1="<< subN(1) << ", M_N2="<< subN(2) <<", M_N="<<this->M_N <<std::endl;
        if( this->M_loadElementsDb )
        {
            if( this->M_elements_database.loadDB() )
            {
                auto size0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>()->primalRB().size();
                auto size1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>()->primalRB().size();
                auto size2 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<2>()->primalRB().size();
                cout<<"Database for basis functions " << this->M_elements_database.lookForDB() << " available and loaded with\n"
                    << size0 <<" primal basis functions in RBSpace0, "
                    << size1 << " primal basis functions in RBSpace1, "
                    << size2 << " primal basis functions in RBSpace2\n";
            }
            else
            {
                this->M_N=0;
            }
        }
    }
    else
    {
        this->M_scmM->setId( this->id() );
        this->M_scmA->setId( this->id() );
        this->M_elements_database.setId( this->id() );
    }
    cout << "Use DB id " << this->id() << std::endl;

    if ( this->M_N == 0 )
    {
        cout<< "Databases does not exist or incomplete -> Start from the begining\n";
        LOG( INFO ) <<"Databases does not exist or incomplete -> Start from the begining";
    }

    // fe vector is requiert in online : must not be TODO
    if ( this->M_use_newton && this->M_loadElementsDb && this->M_Rqm.empty() )
        boost::tie( boost::tuples::ignore, boost::tuples::ignore/*M_Jqm*/, this->M_Rqm ) = this->M_model->computeAffineDecomposition();

    LOG(INFO) << "CRBaero init end\n";
} // init


template <typename TruthModelType>
void
CRBAero<TruthModelType>::buildRbMatrixTrilinear( int number_of_added_elements, parameter_type& mu )
{
    tic();
    int N = this->WNmuSize();
    int N0 = this->subN(0,N);
    int N2 = this->subN(2,N);
    int n0 = N0 - this->subN(0,N-1);
    int n2 = N2 - this->subN(2,N-1);

    auto model = this->M_model;
    auto XN0 = model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN2 = model->rBFunctionSpace()->template rbFunctionSpace<2>();

    element_type Ut = model->functionSpace()->element();
    element_type Ur = model->functionSpace()->element();
    element_type Uc = model->functionSpace()->element();

    auto ut = Ut.template element<0>();
    auto ur = Ur.template element<0>();
    auto uc = Uc.template element<0>();

    auto Tt = Ut.template element<2>();
    auto Tr = Ur.template element<2>();
    auto Tc = Uc.template element<2>();

    // update the reduced operator for block 00
    int q_tri = model->QTri();
    for ( int q=0; q<q_tri; q++ )
    {
        this->blockTriqm(0,0)[q].resize( N0 );
        Ut.zero();
        Ur.zero();
        Uc.zero();

        for ( int k=0; k<N0; k++ )
        {
            ut = XN0->primalBasisElement(k);
            this->blockTriqm(0,0)[q][k]. conservativeResize( N0, N0 );
            auto trilinear_operator = model->assembleTrilinearOperator( Ut, q );
            for ( int i=0; i<N0; i++ )
            {
                ur = XN0->primalBasisElement(i);
                for ( int j=0; j<N0; j++ )
                {
                    if ( k>=N0-n0 || i>=N0-n0 || j>=N0-n0 )
                    {
                        uc = XN0->primalBasisElement(j);
                        this->blockTriqm(0,0)[q][k](i,j) = trilinear_operator->energy( Ur, Uc );
                    }
                }
            }
        }

        this->blockTriqm(2,0)[q].resize(N2);
        Ut.zero();
        Ur.zero();
        Uc.zero();

        for ( int k=0; k<N2; k++ )
        {
            ut = XN0->primalBasisElement(k);
            Tt = XN2->primalBasisElement(k);
            this->blockTriqm(2,0)[q][k]. conservativeResize( N2, N0 );
            auto trilinear_operator = model->assembleTrilinearOperator( Ut, q );
            for ( int i=0; i<N2; i++ )
            {
                Tr = XN2->primalBasisElement(i);
                for ( int j=0; j<N0; j++ )
                {
                    if ( k>=N2-n2 || i>=N2-n2 || j>=N0-n0 )
                    {
                        uc = XN0->primalBasisElement(j);
                        this->blockTriqm(2,0)[q][k](i,j) = trilinear_operator->energy( Uc, Ur );
                    }
                }
            }
        }
    }
    toc("CRBTrilinear : assemble reduced trilinear operator");
}


template <typename TruthModelType>
typename CRBAero<TruthModelType>::matrix_info_tuple
CRBAero<TruthModelType>::onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    tic();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    int sumN = N0+N1+N2;

    M_Jbil.resize( sumN, sumN );
    M_Rli.resize( sumN );

    matrixN_type J(sumN,sumN);
    vectorN_type R(sumN);
    vectorN_type L( sumN );
    uN[0].resize( sumN );
    uN[0].setZero();

    double *R_data = R.data();
    double *J_data = J.data();
    double *U_data = uN[0].data();

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( R_data, sumN );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_UN ( U_data, sumN );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( J_data, sumN , sumN );

    M_Jbil.setZero();
    M_Rli.setZero();

    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = model->computeBetaQm( map_UN, mu );
    for ( int q=0; q<model->sizeOfBilinearJ(); q++ )
    {
        for ( int m=0; m<model->mMaxA(q); m++ )
        {
            // first row
            M_Jbil.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
            M_Jbil.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
            M_Jbil.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);
            // second row
            M_Jbil.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
            M_Jbil.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
            M_Jbil.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);
            // third row
            M_Jbil.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
            M_Jbil.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
            M_Jbil.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
        }
    }

    for ( int q=0; q<model->sizeOfLinearR(); q++ )
    {
        for ( int m=0; m<model->mMaxF(0,q); m++ )
        {
            M_Rli.segment(    0,N0) += betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
            M_Rli.segment(   N0,N1) += betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
            M_Rli.segment(N0+N1,N2) += betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
        }
    }

    this->M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobianOnline, boost::ref( *this ), _1, _2  , mu , N );
    this->M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidualOnline, boost::ref( *this ), _1, _2  , mu , N );

    this->M_nlsolver->solve( map_J , map_UN , map_R, 1e-12, 100 );

    double conditioning = 0;
    double determinant =0;

    if ( computeOutput )
    {
        std::vector<beta_vector_type> betaFqm;
        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaFqm ) = model->computeBetaQm( uN[0], mu );
        L.setZero();
        int out_index = this->M_output_index;
        for ( size_type q=0; q<model->Ql(out_index); q++ )
        {
            for ( size_type m=0; m<model->mMaxF(out_index,q); m++ )
            {
                L.segment( 0,N0) += betaFqm[out_index][q][m]*blockLqm(0)[q][m].head(N0);
                L.segment(N0,N1) += betaFqm[out_index][q][m]*blockLqm(1)[q][m].head(N1);
                L.segment(N0+N1,N2) += betaFqm[out_index][q][m]*blockLqm(2)[q][m].head(N2);
            }
        }
        output_vector[0] = L.dot( uN[0] );
    }
    toc("CRBAERO onlineSolve", false);

    return boost::make_tuple(conditioning,determinant);
}


template <typename TruthModelType>
void
CRBAero<TruthModelType>::updateJacobianOnline( const map_dense_vector_type& X, map_dense_matrix_type& J , parameter_type const& mu , int N ) const
{
    tic();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    beta_vector_type betaJqm;

    J.setZero();
    boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = this->M_model->computeBetaQm( X, mu );
    J += M_Jbil;

    for ( int q=model->sizeOfBilinearJ(); q< model->Qa(); q++ )
    {
        for ( int m=0; m<model->mMaxA(q); m++ )
        {
            // first row
            J.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
            J.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
            J.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);
            // second row
            J.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
            J.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
            J.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);
            // third row
            J.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
            J.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
            J.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
        }
    }

    tic();
    auto betaTri = this->M_model->computeBetaTri( mu );
    for ( int q=0; q<model->QTri(); q++ )
    {
        for ( int k=0; k<N0; k++ )
        {
            for ( int i=0; i<N0; i++ )
            {
                J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].row(i).head(N)).dot(X.segment(0,N0));
                J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].col(i).head(N)).dot(X.segment(0,N0));
            }
        }
        for ( int k=0; k<N2; k++ )
        {
            for( int i=0; i<N2; i++ )
            {
                J(N0+N1+k,N0+N1+i) += betaTri[q]*(blockTriqm(2,0)[q][k].row(i)).dot(X.segment(0,N0));
            }

            for( int j=0; j<N0; j++ )
            {
                J(N0+N1+k,j) += betaTri[q]*(blockTriqm(2,0)[q][k].col(j)).dot(X.segment(N0+N1,N2));
            }
        }
    }
    toc("updateR trilinear", false);
    toc("CRBAERO updateJonline", false);
}

template <typename TruthModelType>
void
CRBAero<TruthModelType>::updateResidualOnline( const map_dense_vector_type& X, map_dense_vector_type& R , parameter_type const& mu , int N ) const
{
    tic();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    int sumN = N0+N1+N2;
    matrixN_type temp( sumN, sumN );

    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = model->computeBetaQm( X, mu );

    R.setZero();
    R += M_Rli;
    R += M_Jbil*X;

    for ( int q=model->sizeOfLinearR(); q<model->Ql(0); q++ )
    {
        for ( int m=0; m<model->mMaxF(0,q); m++ )
        {
            R.segment(    0,N0) += betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
            R.segment(   N0,N1) += betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
            R.segment(N0+N1,N2) += betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
        }
    }

    tic();
    temp.setZero();
    auto betaTri = this->M_model->computeBetaTri( mu );
    for ( int q=0; q<model->QTri(); q++ )
    {
        for ( int k=0; k<N0; k++ )
            for ( int i=0; i<N0; i++ )
                temp(k,i) = (blockTriqm(0,0)[q][k].row(i)).dot(X.segment(0,N0));
        for ( int k=0; k<N2; k++ )
            for ( int i=0; i<N2; i++ )
                temp(N0+N1+k,N0+N1+i) = (blockTriqm(2,0)[q][k].row(i)).dot(X.segment(0,N0));
        R += betaTri[q]*temp*X;
    }
    toc("updateR trilinear", false);
    toc("CRBAERO updateRonline", false);
}


template<typename TruthModelType>
typename CRBAero<TruthModelType>::max_error_type
CRBAero<TruthModelType>::maxErrorBounds( size_type N ) const
{
    std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
    std::vector< vectorN_type > uNduold;
    POUT <<"start maxerrobounds\n";
    double err=0;
    parameter_type mu;

    if( this->M_error_type == CRB_EMPIRICAL && this->M_WNmu->size() == 1 )
    {
        parameter_type mu( this->M_Dmu );
        size_type id;
        boost::tie( mu, id ) = this->M_Xi->max();
        return boost::make_tuple(1e5, mu, 0, 0);
    }

    // we evaluate the error (residual or empirical) for each parameter in the complement of W_MNmu
    for ( int k=0; k<this->M_WNmu_complement->size(); k++ )
    {
        parameter_type const& current_mu = this->M_WNmu_complement->at(k);
        auto o = this->lb( N, current_mu, uN, uNdu, uNold, uNduold );
        double current_err=0;
        if( this->M_error_type == CRB_EMPIRICAL )
            current_err = empiricalError( N, current_mu, o.template get<0>() );
        if ( current_err > err )
        {
            mu = current_mu;
            err = current_err;
        }
    } // loop on M_WNmu_complement
    Environment::worldComm().barrier();

    // we find the proc which has the max residual
    int world_size = Environment::worldComm().globalSize();
    std::vector<double> max_world( world_size );
    mpi::all_gather( Environment::worldComm().globalComm(),
                     err,
                     max_world );
    auto it_max = std::max_element( max_world.begin(), max_world.end() );
    int proc_having_good_mu = it_max - max_world.begin();

    // we broadcast the good parameter and the value of the max residual
    auto tuple = boost::make_tuple( mu, err );
    boost::mpi::broadcast( this->worldComm(), tuple, proc_having_good_mu );
    mu = tuple.template get<0>();
    err = tuple.template get<1>();

    Feel::cout << std::setprecision(15) << "[CRBSaddlePoint] max error="<< err << std::endl;
    POUT <<"finish maxerrobounds\n";
    return boost::make_tuple( err, mu, 0, 0 );
}

template<typename TruthModelType>
double
CRBAero<TruthModelType>::empiricalError( int N, parameter_type const& mu, std::vector<double> output_vec ) const
{
    double output = output_vec[0];
    std::vector<vectorN_type> uN2;
    std::vector<vectorN_type> uNdu2;//( nb_element );
    std::vector<vectorN_type> uNold2;
    std::vector<vectorN_type> uNduold2;
    auto o = this->lb(N-1, mu, uN2, uNdu2, uNold2, uNduold2);
    auto output_vec2 = o.template get<0>();
    auto output2 = output_vec2[0];

    return math::abs(output-output2);
}




template<typename TruthModelType>
template<class Archive>
void
CRBAero<TruthModelType>::serialize( Archive & ar, const unsigned int version )
{
    LOG(INFO) <<"[CRBAero::serialize] version : "<<version<<std::endl;

    ar & boost::serialization::base_object<super>( *this );

} //serialize( ... )



} // namespace Feel

#endif
