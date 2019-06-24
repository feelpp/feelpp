#ifndef CRBAERO_H
#define CRBAERO_H

#include <feel/feelcrb/crbblock.hpp>


#define VERBOSE false

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
    typedef std::shared_ptr<self_type> self_ptrtype;

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
            return New( name, std::make_shared<model_type>(stage), stage );
        }

    static self_ptrtype New( std::string const& name,
                             truth_model_ptrtype const& model,
                             crb::stage stage = crb::stage::online,
                             std::string const& prefixElt = "")
        {
            auto crb = std::shared_ptr<self_type>( new self_type(name, model, stage, prefixElt ));
            crb->init();
            return crb;
        }


protected:
    //! Default Constructor
    CRBAero( std::string const& name = "defaultname_crb",
              crb::stage stage = crb::stage::online,
              WorldComm const& worldComm = Environment::worldComm() ) :
        CRBAero( name, std::make_shared<model_type>(stage), stage )
        {}

    //! constructor from command line options
    CRBAero( std::string const& name, truth_model_ptrtype const & model,
             crb::stage stage = crb::stage::online, std::string const& prefixExt = "" ) :
        super( name, model, stage, prefixExt ),
        M_use_psit( boption(_prefix=this->M_prefix,_name="crb.aero.use-psit")),
        M_newton( boption(_prefix=this->M_prefix,_name="crb.aero.use-newton")),
        M_delta( doption(_prefix=this->M_prefix,_name="crb.aero.psit.delta0" ) ),
        M_rez(-1),
        M_store_rb_sol( boption(_prefix=this->M_prefix,_name="crb.aero.store-rb-sol"))
        {
            CHECK( n_block==3 ) <<"CRBAero is supposed to work with a composite space u,p,T\n";
        }



public:
    //! Init function, specific to crb aero
    void init();
    void offlineSolve( element_type& u, element_type& udu, parameter_type& mu, element_ptrtype & dual_initial_field )override
        {
            if ( boption(_prefix=this->M_prefix,_name="crb.solve-fem-monolithic") )
                u = this->M_model->solve(mu);
            else
                u = this->M_model->offlineSolveAD(mu);
        }

    element_type solve( parameter_type const& mu )
    {
        return this->M_model->solve( mu );
    }

    matrix_info_tuple onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const override;
    matrix_info_tuple onlineSolveNewton(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const;
    matrix_info_tuple onlineSolvePicard(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const;

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
    std::vector<bool> hasZeroMean() override
        {
            std::vector<bool> out( n_block,false );
            out[1] = this->model()->hasZeroMeanPressure();
            return out;
        }

    std::map<std::string,double> timerMap() const override
        {
            std::map<std::string,double> m;
            m["J"] = timerJ;
            m["Jtri"] = timerJtri;
            if ( timerJnl>1e-12 )
                m["Jnl"] = timerJnl;
            m["R"] = timerR;
            m["Rtri"] = timerRtri;
            if ( timerRnl>1e-12 )
                m["Rnl"] = timerRnl;
            m["Solve"] = timerSolve;
            if ( timerBeta>1e-12 )
                m["Beta"] = timerBeta;
            m["Niter"] = M_Niter;
            return m;
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

    double newR(const map_dense_vector_type& X, parameter_type const& mu , int N) const
        {
            int N0 = this->subN(0,N);
            int N1 = this->subN(1,N);
            int N2 = this->subN(2,N);
            int sumN = N0+N1+N2;
            vectorN_type R(sumN);
            double *R_data = R.data();
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( R_data, sumN );
            updateResidualOnline(X,map_R,mu,N);

            return R.norm();
        }
    void updatePsiT( const map_dense_vector_type& X, parameter_type const& mu, int N ) const
        {
            double new_rez = newR( X, mu, N );
            if ( M_rez==-1 )
                M_rez=new_rez;
            M_delta = M_delta*M_rez/new_rez;
            M_rez = new_rez;
        }


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version );

    bool M_use_psit, M_newton;
    mutable double M_delta, M_rez;

    mutable matrixN_type M_Jbil,M_M;
    mutable vectorN_type M_Rli;
    std::vector<matrixN_type> M_mass;

    mutable double timerJ, timerJtri, timerJnl, timerR, timerRtri, timerRnl, timerSolve, timerBeta, Niter;
    mutable int M_Niter;

    bool M_store_rb_sol;
    mutable std::map<std::string,std::vector<vectorN_type>> M_rb_sol;
    mutable std::map<std::string,std::vector<bool>> M_sol_done;

    using super::n_block;
    using super::subN;
    using super::blockAqm;
    using super::blockFqm;
    using super::blockLqm;
    using super::blockTriqm;
    using super::notEmptyAqm;
    using super::notEmptyFqm;
    using super::notEmptyLqm;
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
                        this->blockTriqm(0,0)[q][k](i,j) = trilinear_operator->energy( Uc, Ur );
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
            //ut = XN0->primalBasisElement(k);
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
    if ( M_use_psit )
    {
        tic();
        if ( M_mass.size()==0 )
            M_mass.resize(3);
        auto XN0 = this->model()->rBFunctionSpace()->template rbFunctionSpace<0>();
        auto XN2 = this->model()->rBFunctionSpace()->template rbFunctionSpace<2>();
        M_mass[0].conservativeResize(N0,N0);
        M_mass[2].conservativeResize(N2,N2);
        for ( int i=0; i<N0; i++ )
            for( int j=0; j<N0; j++ )
            {
                if ( i>=N0-n0 || j>=N0-n0 )
                {
                    if ( i==j )
                        M_mass[0](i,j) = this->model()->scalarProduct( XN0->primalBasisElement(i),
                                                                       XN0->primalBasisElement(j), 0 );
                    else
                        M_mass[0](i,j)=0;
                }
            }
        for ( int i=0; i<N2; i++ )
            for( int j=0; j<N2; j++ )
            {
                if ( i>=N2-n2 || j>=N2-n2 )
                {
                    if ( i==j )
                        M_mass[2](i,j) = this->model()->scalarProduct( XN2->primalBasisElement(i),
                                                                       XN2->primalBasisElement(j), 2 );
                    else
                        M_mass[2](i,j)=0;
                }
            }
        toc("CRBAero: update mass");
    }
    toc("CRBTrilinear : assemble reduced trilinear operator");
}

template <typename TruthModelType>
typename CRBAero<TruthModelType>::matrix_info_tuple
CRBAero<TruthModelType>::onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    timerJ=0; timerJtri=0; timerJnl=0; timerR=0; timerRtri=0; timerRnl=0; timerSolve=0; timerBeta=0;
    M_Niter=0;
    boost::mpi::timer t_solve;
    t_solve.restart();
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    int sumN = N0+N1+N2;
    typename CRBAero<TruthModelType>::matrix_info_tuple out =  boost::make_tuple(0.,0.);

    if ( M_store_rb_sol && !this->M_check_cvg )
    {
        if ( M_sol_done[mu.toString()].size() )
        {
            if ( M_sol_done[mu.toString()][N-1] )
            {
                uN[0] = M_rb_sol[mu.toString()][N-1];
                return out;
            }
        }
    }

    int number_of_neighbors = this->M_WNmu->size() - N+1;
    std::vector<int> index_vector;
    auto S = this->M_WNmu->searchNearestNeighbors( mu, number_of_neighbors, index_vector, false);
    int n_index=0;
    int index;
    do
    {
        index=index_vector[n_index];
        n_index++;
    }while( index>=N );

    auto mu_init = this->M_WNmu->at(index);
    auto u_init = this->projSol()[mu_init.toString()];
    if ( !(u_init.size()==3) )
    {
        Feel::cout << "u_init not stored for mu="<<mu_init.toString()<<std::endl;
        u_init=this->projSol().begin()->second;
    }

    auto u_prev = uN[0];
    int n0 = subN(0,N-1);
    int n1 = subN(1,N-1);
    int n2 = subN(2,N-1);
    int sumn = n0+n1+n2;

    uN[0].resize( sumN );
    uN[0].setZero();

    if ( this->M_check_cvg && u_prev.size()==sumn )
    {
        uN[0].segment(     0, n0 ) = u_prev.segment(     0, n0 );
        uN[0].segment(    N0, n1 ) = u_prev.segment(    n0, n1 );
        uN[0].segment( N0+N1, n2 ) = u_prev.segment( n0+n1, n2 );
    }
    else if ( boption(_prefix=this->M_prefix,_name="crb.aero.init-online"))
    {
        uN[0].segment(     0,N0 ) = u_init[0].head(N0);
        uN[0].segment(    N0,N1 ) = u_init[1].head(N1);
        uN[0].segment( N0+N1,N2 ) = u_init[2].head(N2);

    }

    if ( M_newton )
        onlineSolveNewton( N, mu, uN, uNdu, uNold, uNduold, output_vector, K, print_rb_matrix, computeOutput);
    else
        onlineSolvePicard( N, mu, uN, uNdu, uNold, uNduold, output_vector, K, print_rb_matrix, computeOutput);

    if ( !this->M_last_online_converged && ioption(_prefix=this->M_prefix,_name="crb.aero.online-continuation") )
    {
        Feel::cout << "CRBAero: use continuation for mu="<<mu.toString()<<std::endl;
        int n = ioption(_prefix=this->M_prefix,_name="crb.aero.online-continuation");
        uN[0].segment(     0,N0 ) = u_init[0].head(N0);
        uN[0].segment(    N0,N1 ) = u_init[1].head(N1);
        uN[0].segment( N0+N1,N2 ) = u_init[2].head(N2);

        for ( int i=0; i<=n; i++ )
        {
            parameter_type current_mu = (n-i)/(double)n*mu_init + i/(double)n*mu;
            if ( boption(_prefix=this->M_prefix,_name="crb.aero.log-continuation") )
                for ( int k=0; k<mu.size(); k++ )
                {
                    current_mu(k) = std::pow(mu_init(k),(n-i)/(double)n )*std::pow(mu(k),i/(double)n );
                }

            //Feel::cout << "CRBAero: online continuation current mu="<<current_mu.toString()<<std::endl;
            if ( M_newton )
                onlineSolveNewton( N, current_mu, uN, uNdu, uNold, uNduold, output_vector, K, print_rb_matrix, computeOutput);
            else
                onlineSolvePicard( N, current_mu, uN, uNdu, uNold, uNduold, output_vector, K, print_rb_matrix, computeOutput);
             if ( !this->M_last_online_converged )
                 Feel::cout << "CRBAero: continuation failed for i="<<i<<", mu="<<current_mu.toString()<<std::endl;
        }
    }

    // if ( !this->M_last_online_converged )
    //     Feel::cout << "CRBAero: Online Solver failed to converge, rez="<<M_rez <<std::endl;
    timerSolve = t_solve.elapsed();
    if ( M_store_rb_sol && !this->M_check_cvg )
    {
        if ( !M_sol_done[mu.toString()].size() )
        {
            M_sol_done[mu.toString()].resize(ioption(_prefix=this->M_prefix,_name="crb.dimension-max"));
            M_rb_sol[mu.toString()].resize(ioption(_prefix=this->M_prefix,_name="crb.dimension-max"));
        }
        M_sol_done[mu.toString()][N-1]=true;
        M_rb_sol[mu.toString()][N-1]=uN[0];
    }

    return out;
}

template <typename TruthModelType>
typename CRBAero<TruthModelType>::matrix_info_tuple
CRBAero<TruthModelType>::onlineSolvePicard(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    tic();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    if ( VERBOSE )
        Feel::cout << "CRBAero Online solve with N="<<N<<", N0="<<N0<<", N1="<<N1<<", N2="<<N2<<std::endl;
    int sumN = N0+N1+N2;
    double increment = this->M_fixedpointIncrementTol;

    M_Jbil.resize( sumN, sumN );
    M_Rli.resize( sumN );
    M_Jbil.setZero();
    M_Rli.setZero();

    matrixN_type J(sumN,sumN);
    vectorN_type R(sumN);
    vectorN_type L( sumN );
    //uN[0].resize( sumN );
    vectorN_type previous_uN( sumN );
    int fi=0;
    bool fixPointIsFinished = false;

    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = model->computeBetaQm( uN[0], mu );

    tic();
    for ( int q=0; q<model->sizeOfBilinearJ(); q++ )
    {
        for ( int m=0; m<model->mMaxA(q); m++ )
        {
            // first row
            if ( this->notEmptyAqm(0,0,q,m) )
                M_Jbil.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
            if ( this->notEmptyAqm(0,1,q,m) )
                M_Jbil.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
            if ( this->notEmptyAqm(0,2,q,m) )
                M_Jbil.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);

            // second row
            if ( this->notEmptyAqm(1,0,q,m) )
                M_Jbil.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
            if ( this->notEmptyAqm(1,1,q,m) )
                M_Jbil.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
            if ( this->notEmptyAqm(1,2,q,m) )
                M_Jbil.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);

            // third row
            if ( this->notEmptyAqm(2,0,q,m) )
                M_Jbil.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
            if ( this->notEmptyAqm(2,1,q,m) )
                M_Jbil.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
            if ( this->notEmptyAqm(2,2,q,m) )
                M_Jbil.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
        }
    }
    toc("Preassemble Bilinear terms", VERBOSE );
    tic();
    for ( int q=0; q<model->sizeOfLinearR(); q++ )
    {
        for ( int m=0; m<model->mMaxF(0,q); m++ )
        {
            if ( this->notEmptyFqm(0,q,m) )
                M_Rli.segment(    0,N0) += betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
            if ( this->notEmptyFqm(1,q,m) )
                M_Rli.segment(   N0,N1) += betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
            if ( this->notEmptyFqm(2,q,m) )
                M_Rli.segment(N0+N1,N2) += betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
        }
    }
    toc("Preassemble Linear terms", VERBOSE);

    do
    {
        previous_uN = uN[0];

        J.setZero();
        R.setZero();

        if ( model->sizeOfBilinearJ()!=this->model()->Qa() )
        {
            boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) =
                this->M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
            for ( int q=model->sizeOfBilinearJ(); q< model->Qa(); q++ )
            {
                for ( int m=0; m<model->mMaxA(q); m++ )
                {
                    // first row
                    if ( this->notEmptyAqm(0,0,q,m) )
                        J.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
                    if ( this->notEmptyAqm(0,1,q,m) )
                        J.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
                    if ( this->notEmptyAqm(0,2,q,m) )
                        J.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);

                    // second row
                    if ( this->notEmptyAqm(1,0,q,m) )
                        J.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
                    if ( this->notEmptyAqm(1,1,q,m) )
                        J.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
                    if ( this->notEmptyAqm(1,2,q,m) )
                        J.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);

                    // third row
                    if ( this->notEmptyAqm(2,0,q,m) )
                        J.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
                    if ( this->notEmptyAqm(2,1,q,m) )
                        J.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
                    if ( this->notEmptyAqm(2,2,q,m) )
                        J.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
                }
            }
            for ( int q=model->sizeOfLinearR(); q<model->Ql(0); q++ )
            {
                for ( int m=0; m<model->mMaxF(0,q); m++ )
                {
                    if ( this->notEmptyFqm(0,q,m) )
                        R.segment(    0,N0) += betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
                    if ( this->notEmptyFqm(1,q,m) )
                        R.segment(   N0,N1) += betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
                    if ( this->notEmptyFqm(2,q,m) )
                        R.segment(N0+N1,N2) += betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
                }
            }
        }

        auto betaTri = this->M_model->computeBetaTri( mu );
        for ( int q=0; q<model->QTri(); q++ )
        {
            for ( int k=0; k<N0; k++ )
                for ( int i=0; i<N0; i++ )
                    J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].row(i).head(N0)).dot(uN[0].segment(0,N0));
            for ( int k=0; k<N2; k++ )
                for( int i=0; i<N2; i++ )
                    J(N0+N1+k,N0+N1+i) += betaTri[q]*(blockTriqm(2,0)[q][k].row(i).head(N0)).dot(uN[0].segment(0,N0));
        }

        J += M_Jbil;
        R += M_Rli;



        uN[0] = J.fullPivLu().solve( R );

        increment = (uN[0]-previous_uN).norm();
        auto increment_abs = (uN[0]-previous_uN);
        fixPointIsFinished = increment < this->M_fixedpointIncrementTol || fi>=this->M_fixedpointMaxIterations;
        this->online_iterations_summary.first = fi;
        this->online_iterations_summary.second = increment;

        for ( int q=0; q<model->QTri(); q++ )
        {
            for ( int k=0; k<N0; k++ )
                for ( int i=0; i<N0; i++ )
                    J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].row(i).head(N0)).dot(uN[0].segment(0,N0));
            for ( int k=0; k<N2; k++ )
                for( int i=0; i<N2; i++ )
                    J(N0+N1+k,N0+N1+i) += betaTri[q]*(blockTriqm(2,0)[q][k].row(i).head(N0)).dot(uN[0].segment(0,N0));
        }
        M_rez = (J * uN[0] - R).norm() ;


        if( this->M_fixedpointVerbose )
            Feel::cout << "[CRBAero::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment << ", incrment_bas="<< M_rez<<std::endl;
        fi++;

    }while ( !fixPointIsFinished );

    this->M_last_online_converged = increment<this->M_fixedpointIncrementTol;

    // if ( !this->M_last_online_converged )
    //     Feel::cout << "CRBAero: Picard didn't converged, increment="<<increment<<", residual="<<residual_norm<<std::endl;

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
                if ( this->notEmptyLqm(0,q,m) )
                    L.segment( 0,N0) += betaFqm[out_index][q][m]*blockLqm(0)[q][m].head(N0);
                if ( this->notEmptyLqm(1,q,m) )
                    L.segment(N0,N1) += betaFqm[out_index][q][m]*blockLqm(1)[q][m].head(N1);
                if ( this->notEmptyLqm(2,q,m) )
                    L.segment(N0+N1,N2) += betaFqm[out_index][q][m]*blockLqm(2)[q][m].head(N2);
            }
        }
        output_vector[0] = L.dot( uN[0] );
    }

    toc("CRBAERO onlineSolve",VERBOSE);

    double conditioning = 0;
    double determinant = 0;

    return boost::make_tuple(conditioning,determinant);

}

template <typename TruthModelType>
typename CRBAero<TruthModelType>::matrix_info_tuple
CRBAero<TruthModelType>::onlineSolveNewton(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    tic();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    if ( VERBOSE )
        Feel::cout << "CRBAero Online solve with N="<<N<<", N0="<<N0<<", N1="<<N1<<", N2="<<N2<<std::endl;
    int sumN = N0+N1+N2;

    M_Jbil.resize( sumN, sumN );
    M_Rli.resize( sumN );

    matrixN_type J(sumN,sumN);
    vectorN_type R(sumN);
    vectorN_type L( sumN );
    //uN[0].resize( sumN );
    //uN[0].setZero();

    double *R_data = R.data();
    double *J_data = J.data();
    double *U_data = uN[0].data();

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( R_data, sumN );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_UN ( U_data, sumN );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( J_data, sumN , sumN );

    M_Jbil.setZero();
    M_Rli.setZero();


    boost::mpi::timer tBeta;
    tBeta.restart();
    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = model->computeBetaQm( map_UN, mu );
    timerBeta += tBeta.elapsed();

    tic();
    for ( int q=0; q<model->sizeOfBilinearJ(); q++ )
    {
        for ( int m=0; m<model->mMaxA(q); m++ )
        {
            // first row
            if ( this->notEmptyAqm(0,0,q,m) )
                M_Jbil.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
            if ( this->notEmptyAqm(0,1,q,m) )
                M_Jbil.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
            if ( this->notEmptyAqm(0,2,q,m) )
                M_Jbil.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);

            // second row
            if ( this->notEmptyAqm(1,0,q,m) )
                M_Jbil.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
            if ( this->notEmptyAqm(1,1,q,m) )
                M_Jbil.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
            if ( this->notEmptyAqm(1,2,q,m) )
                M_Jbil.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);

            // third row
            if ( this->notEmptyAqm(2,0,q,m) )
                M_Jbil.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
            if ( this->notEmptyAqm(2,1,q,m) )
                M_Jbil.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
            if ( this->notEmptyAqm(2,2,q,m) )
                M_Jbil.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
        }
    }
    toc("Preassemble Bilinear terms", VERBOSE );
    tic();
    for ( int q=0; q<model->sizeOfLinearR(); q++ )
    {
        for ( int m=0; m<model->mMaxF(0,q); m++ )
        {
            if ( this->notEmptyFqm(0,q,m) )
                M_Rli.segment(    0,N0) += -betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
            if ( this->notEmptyFqm(1,q,m) )
                M_Rli.segment(   N0,N1) += -betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
            if ( this->notEmptyFqm(2,q,m) )
                M_Rli.segment(N0+N1,N2) += -betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
        }
    }
    toc("Preassemble Linear terms", VERBOSE);

    if ( M_use_psit )
    {
        M_M.resize( sumN, sumN);
        M_M.setZero();
        M_M.block(0,0,N0,N0) = M_mass[0].block(0,0,N0,N0);
        M_M.block(N0+N1,N0+N1,N2,N2) = M_mass[2].block(0,0,N2,N2);
    }


    //Feel::cout << "Jbil=\n"<<M_Jbil<<"\n Rli=\n"<<M_Rli<<std::endl;
    this->M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobianOnline, boost::ref( *this ), _1, _2  , mu , N );
    this->M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidualOnline, boost::ref( *this ), _1, _2  , mu , N );

    this->M_nlsolver->setType( TRUST_REGION );
    this->M_nlsolver->setRelativeResidualTol( doption(_prefix=this->M_prefix,_name="crb.aero.snes.rtol"));
    this->M_nlsolver->setKspSolverType( PREONLY );
    auto out = this->M_nlsolver->solve( map_J , map_UN , map_R, 1e-12, 100 );

    this->M_last_online_converged = out.first>=0 && out.first<20;
    M_rez = newR( map_UN, mu, N );
    // if ( !this->M_last_online_converged )
    //     Feel::cout << "CRBAero Newton did not converged\n";

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
                if ( this->notEmptyLqm(0,q,m) )
                    L.segment( 0,N0) += betaFqm[out_index][q][m]*blockLqm(0)[q][m].head(N0);
                if ( this->notEmptyLqm(1,q,m) )
                    L.segment(N0,N1) += betaFqm[out_index][q][m]*blockLqm(1)[q][m].head(N1);
                if ( this->notEmptyLqm(2,q,m) )
                    L.segment(N0+N1,N2) += betaFqm[out_index][q][m]*blockLqm(2)[q][m].head(N2);
            }
        }
        output_vector[0] = L.dot( uN[0] );
    }
    toc("CRBAERO onlineSolve",VERBOSE);

    return boost::make_tuple(conditioning,determinant);
}


template <typename TruthModelType>
void
CRBAero<TruthModelType>::updateJacobianOnline( const map_dense_vector_type& X, map_dense_matrix_type& J , parameter_type const& mu , int N ) const
{
    tic();
    boost::mpi::timer tJ, tJtri, tJnl, tBeta;
    M_Niter++;
    tJ.restart();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    beta_vector_type betaJqm;

    J.setZero();
    tBeta.restart();
    boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = this->M_model->computeBetaQm( X, mu );
    timerBeta += tBeta.elapsed();
    J += M_Jbil;
    if ( M_use_psit )
    {
        updatePsiT( X, mu, N );
        J+=1./M_delta*M_M;
    }


    tic();
    for ( int q=model->sizeOfBilinearJ(); q< model->Qa(); q++ )
    {
        tJnl.restart();
        for ( int m=0; m<model->mMaxA(q); m++ )
        {
            // first row
            if ( this->notEmptyAqm(0,0,q,m) )
                J.block( 0,     0, N0, N0) += betaJqm[q][m]*blockAqm(0,0)[q][m].block(0,0,N0,N0);
            if ( this->notEmptyAqm(0,1,q,m) )
                J.block( 0,    N0, N0, N1) += betaJqm[q][m]*blockAqm(0,1)[q][m].block(0,0,N0,N1);
            if ( this->notEmptyAqm(0,2,q,m) )
                J.block( 0, N0+N1, N0, N2) += betaJqm[q][m]*blockAqm(0,2)[q][m].block(0,0,N0,N2);

            // second row
            if ( this->notEmptyAqm(1,0,q,m) )
                J.block( N0,     0, N1, N0) += betaJqm[q][m]*blockAqm(1,0)[q][m].block(0,0,N1,N0);
            if ( this->notEmptyAqm(1,1,q,m) )
                J.block( N0,    N0, N1, N1) += betaJqm[q][m]*blockAqm(1,1)[q][m].block(0,0,N1,N1);
            if ( this->notEmptyAqm(1,2,q,m) )
                J.block( N0, N0+N1, N1, N2) += betaJqm[q][m]*blockAqm(1,2)[q][m].block(0,0,N1,N2);

            // third row
            if ( this->notEmptyAqm(2,0,q,m) )
                J.block( N0+N1,     0, N2, N0) += betaJqm[q][m]*blockAqm(2,0)[q][m].block(0,0,N2,N0);
            if ( this->notEmptyAqm(2,1,q,m) )
                J.block( N0+N1,    N0, N2, N1) += betaJqm[q][m]*blockAqm(2,1)[q][m].block(0,0,N2,N1);
            if ( this->notEmptyAqm(2,2,q,m) )
                J.block( N0+N1, N0+N1, N2, N2) += betaJqm[q][m]*blockAqm(2,2)[q][m].block(0,0,N2,N2);
        }
        timerJnl += tJnl.elapsed();
    }
    toc("UpdateJ NL terms",VERBOSE);

    tic();
    tJtri.restart();
    auto betaTri = this->M_model->computeBetaTri( mu );
    for ( int q=0; q<model->QTri(); q++ )
    {
        if ( ioption(_prefix=this->M_prefix,_name="crb.aero.assemble-version")==1 )
        {
            for ( int k=0; k<N0; k++ )
            {
                J.block(k,0,1,N0) += betaTri[q]*((blockTriqm(0,0)[q][k].block(0,0,N0,N0))*(X.segment(0,N0))).transpose();
                J.block(k,0,1,N0) += betaTri[q]*((blockTriqm(0,0)[q][k].block(0,0,N0,N0).transpose())*(X.segment(0,N0))).transpose();
            }
            for ( int k=0; k<N2; k++ )
            {
                J.block(N0+N1+k,N0+N1,1,N2) += betaTri[q]*((blockTriqm(2,0)[q][k].block(0,0,N2,N0))*(X.segment(0,N0))).transpose();
                J.block(N0+N1+k,0,1,N0) += betaTri[q]*((blockTriqm(2,0)[q][k].block(0,0,N2,N0).transpose())*(X.segment(N0+N1,N2))).transpose();
            }

        }
        else
        {
            for ( int k=0; k<N0; k++ )
            {
                for ( int i=0; i<N0; i++ )
                {
                    J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].row(i).head(N0)).dot(X.segment(0,N0));
                    J(k,i) += betaTri[q]*(blockTriqm(0,0)[q][k].col(i).head(N0)).dot(X.segment(0,N0));
                }
            }
            for ( int k=0; k<N2; k++ )
            {
                for( int i=0; i<N2; i++ )
                {
                    J(N0+N1+k,N0+N1+i) += betaTri[q]*(blockTriqm(2,0)[q][k].row(i).head(N0)).dot(X.segment(0,N0));
                }

                for( int j=0; j<N0; j++ )
                {
                    J(N0+N1+k,j) += betaTri[q]*(blockTriqm(2,0)[q][k].col(j).head(N2)).dot(X.segment(N0+N1,N2));
                }
            }
        }
    }
    timerJtri += tJtri.elapsed();
    toc("updateJ trilinear", VERBOSE);
    timerJ += tJ.elapsed();
    toc("CRBAERO updateJonline", VERBOSE);
}

template <typename TruthModelType>
void
CRBAero<TruthModelType>::updateResidualOnline( const map_dense_vector_type& X, map_dense_vector_type& R , parameter_type const& mu , int N ) const
{
    tic();
    boost::mpi::timer tR, tRtri, tRnl, tBeta;
    tR.restart();
    auto model = this->M_model;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);
    int N2 = this->subN(2,N);
    int sumN = N0+N1+N2;
    matrixN_type temp( sumN, sumN );
    vectorN_type myR( sumN );

    tBeta.restart();
    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = model->computeBetaQm( X, mu );
    timerBeta += tBeta.elapsed();

    tic();
    R.setZero();
    myR.setZero();
    R += M_Rli;
    R += M_Jbil*X;
    toc("UpdateR some",VERBOSE);

    tic();
    for ( int q=model->sizeOfLinearR(); q<model->Ql(0); q++ )
    {
        tRnl.restart();
        for ( int m=0; m<model->mMaxF(0,q); m++ )
        {
            if ( this->notEmptyFqm(0,q,m) )
                R.segment(    0,N0) += betaRqm[0][q][m]*blockFqm(0)[q][m].head(N0);
            if ( this->notEmptyFqm(1,q,m) )
                R.segment(   N0,N1) += betaRqm[0][q][m]*blockFqm(1)[q][m].head(N1);
            if ( this->notEmptyFqm(2,q,m) )
                R.segment(N0+N1,N2) += betaRqm[0][q][m]*blockFqm(2)[q][m].head(N2);
        }
        timerRnl += tRnl.elapsed();
    }
    toc("UpdateR NL terms",VERBOSE);

    tic();
    tRtri.restart();
    temp.setZero();
    auto betaTri = this->M_model->computeBetaTri( mu );
    for ( int q=0; q<model->QTri(); q++ )
    {
        if ( ioption(_prefix=this->M_prefix,_name="crb.aero.assemble-version")==1 )
        {
            for ( int k=0; k<N0; k++ )
                temp.block(k,0,1,N0) = ((blockTriqm(0,0)[q][k].block(0,0,N0,N0))*(X.segment(0,N0))).transpose();
            for ( int k=0; k<N2; k++ )
                temp.block( N0+N1+k,N0+N1,1,N2) = ((blockTriqm(2,0)[q][k].block(0,0,N2,N0))*(X.segment(0,N0))).transpose();
        }
        else
        {
            for ( int k=0; k<N0; k++ )
                for ( int i=0; i<N0; i++ )
                    temp(k,i) = (blockTriqm(0,0)[q][k].row(i).head(N0)).dot(X.segment(0,N0));
            for ( int k=0; k<N2; k++ )
                for ( int i=0; i<N2; i++ )
                    temp(N0+N1+k,N0+N1+i) = (blockTriqm(2,0)[q][k].row(i).head(N0)).dot(X.segment(0,N0));
        }
        myR += betaTri[q]*temp*X;
    }
    R += myR;
    timerRtri += tRtri.elapsed();
    timerR += tR.elapsed();
    toc("updateR trilinear", VERBOSE);
    toc("CRBAERO updateRonline", VERBOSE);
}


template<typename TruthModelType>
typename CRBAero<TruthModelType>::max_error_type
CRBAero<TruthModelType>::maxErrorBounds( size_type N ) const
{
    std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
    std::vector< vectorN_type > uNduold;

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
