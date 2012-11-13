/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009-2012 UniversitÈ Joseph Fourier (Grenoble I)

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
   \file crb_trilinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Elisa Schenone
   \author Stéphane Veys
   \date 2012-11-10
 */
#ifndef __CRB_TRILINEAR_H
#define __CRB_TRILINEAR_H 1


#include <feel/feelcrb/crb.hpp>


namespace Feel
{
    
/**
 * \class CRB_TRILINEAR
 * \brief Certifed Reduced Basis for Trilinear forms class
 *
 * Implements the certified reduced basis method for treat trilinear forms
 * Herites from CRB class
 *
 * @author Elisa Schenone, Stéphane Veys
 * @see
 */
template<typename TruthModelType>
class CRB_TRILINEAR : public CRB
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRB_TRILINEAR()
        :
        CRB()
    {

    }

    //! constructor from command line options
    CRB( std::string  name, po::variables_map const& vm )
        :
        CRB( name, vm )
    {

    }

    //! constructor from command line options
    CRB_TRILINEAR_TR( std::string  name,
         po::variables_map const& vm,
         truth_model_ptrtype const & model )
        :
        CRB( name, vm, mode) )
    {

    }


    //! copy constructor
    CRB_TRILINEAR( CRB_TRILINEAR const & o )
        :
    CRB( o )
    {}

    //! destructor
    ~CRB_TRILINEAR()
    {}
    
    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     * \param K : index of time ( time = K*dt) at which we want to evaluate the output
     * Note : K as a default value for non time-dependent problems
     *
     *\return compute online the lower bound
     *\and also condition number of matrix A
     */

    boost::tuple<double,double> lb( size_type N, parameter_type const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu , std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, int K=0 ) const;

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    convergence_type offline();
    
    /**
     * Update the Jacobian Matrix for Newton Solver
     *
     */
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    
    /**
     * Update the Residual of the Newton Solver
     *
     */
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );

    //@}


private:
    std::vector < std::vector < std::vector<matrixN_type> > > M_Aqm_tril_pr;
    sparse_matrix_ptrtype M_bilinear_terms;
    vector_ptrtype M_linear_terms;    
};

template<typename TruthModelType>
typename CRB<TruthModelType>::convergence_type
CRB_TRILINEAR<TruthModelType>::offline()
{
    M_database_contains_variance_info = false;

    int proc_number = this->worldComm().globalRank();
    M_rbconv_contains_primal_and_dual_contributions = true;

    bool rebuild_database = this->vm()["crb.rebuild-database"].template as<bool>() ;
    orthonormalize_primal = this->vm()["crb.orthonormalize-primal"].template as<bool>() ;

    solve_dual_problem = false;
    orthonormalize_dual=false;

    M_Nm = this->vm()["crb.Nm"].template as<int>() ;
    bool seek_mu_in_complement = this->vm()["crb.seek-mu-in-complement"].template as<bool>() ;

    boost::timer ti;
    if( proc_number == 0 ) std::cout << "Offline CRB starts, this may take a while until Database is computed...\n";
    LOG(INFO) << "[CRB::offline] Starting offline for output " << M_output_index << "\n";
    LOG(INFO) << "[CRB::offline] initialize underlying finite element model\n";
    M_model->init();
    if( proc_number == 0 ) std::cout << " -- model init done in " << ti.elapsed() << "s\n";

    parameter_type mu( M_Dmu );

    double maxerror;
    double delta_pr;
    size_type index;
    int no_residual_index;
    //if M_N == 0 then there is not an already existing database
    if ( rebuild_database || M_N == 0)
    {

        ti.restart();
        //scm_ptrtype M_scm = scm_ptrtype( new scm_type( M_vm ) );
        //M_scm->setTruthModel( M_model );
        //    std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scm->offline();

        LOG(INFO) << "[CRB::offline] compute random sampling\n";


        int sampling_size = this->vm()["crb.sampling-size"].template as<int>();
        
        // random sampling
        M_Xi->randomize( sampling_size );
        M_Xi->writeOnFile(file_name);        
        M_WNmu->setSuperSampling( M_Xi );


        if( proc_number == 0 ) std::cout<<"[CRB offline] M_error_type = "<<M_error_type<<std::endl;

        std::cout << " -- sampling init done in " << ti.elapsed() << "s\n";
        ti.restart();

        // empty sets
        M_WNmu->clear();

        mu = M_Dmu->element();

        int size = mu.size();
        if( proc_number == 0 )
        {
            std::cout << "  -- start with mu = [ ";
            for ( int i=0; i<size-1; i++ ) std::cout<<mu( i )<<" ";
            std::cout<<mu( size-1 )<<" ]"<<std::endl;
        }
        //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

        // dimension of reduced basis space
        M_N = 0;

        maxerror = 1e10;
        delta_pr = 0;
        delta_du = 0;
        no_residual_index = 0;
        //boost::tie( maxerror, mu, index ) = maxErrorBounds( N );

        LOG(INFO) << "[CRB::offline] allocate reduced basis data structures\n";
        M_Aqm_pr.resize( M_model->Qa() );
        for(int q=0; q<M_model->Qa(); q++)
        {
            M_Aqm_pr[q].resize( 1 );
        }

        M_Aqm_tril_pr.resize( M_model->Qa_tril() );
        for(int q=0; q<M_model->Qa_tril(); q++)
        {
            M_Aqm_tril_pr[q].resize( 1 );
            M_Aqm_tril_pr[q][0].resize( M_model->mMaxA(q) );
        }

        M_Mqm_pr.resize( M_model->Qm() );
        for(int q=0; q<M_model->Qm(); q++)
        {
            M_Mqm_pr[q].resize( M_model->mMaxM(q) );
        }

        M_MFqm_pr.resize( M_model->Qmf() );
        for(int q=0; q<M_model->Qmf(); q++)
        {
            M_MFqm_pr[q].resize( 1 );
        }

        M_Fqm_pr.resize( M_model->Ql( 0 ) );

        for(int q=0; q<M_model->Ql( 0 ); q++)
        {
            M_Fqm_pr[q].resize( 1 );
        }
    }//end of if( rebuild_database )
#if 0
    else
    {
        mu = M_current_mu;
        no_residual_index=M_no_residual_index;
        if( proc_number == 0 )
        {
            std::cout<<"we are going to enrich the reduced basis"<<std::endl;
            std::cout<<"there are "<<M_N<<" elements in the database"<<std::endl;
        }
        LOG(INFO) <<"we are going to enrich the reduced basis"<<std::endl;
        LOG(INFO) <<"there are "<<M_N<<" elements in the database"<<std::endl;
    }//end of else associated to if ( rebuild_databse )
#endif
    
    sparse_matrix_ptrtype M,A; //,Adu,At;
    element_ptrtype InitialGuess;
    //vector_ptrtype MF;
    std::vector<vector_ptrtype> F,L;

    LOG(INFO) << "[CRB::offline] compute affine decomposition\n";
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm_tril;
    std::vector< std::vector<sparse_matrix_ptrtype> > Mqm;
    std::vector< std::vector<vector_ptrtype> > MFqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm;

    boost::tie( Mqm, Aqm, Fqm, MFqm ) = M_model->computeAffineDecomposition();

    element_ptrtype u( new element_type( M_model->functionSpace() ) );
    element_ptrtype uproj( new element_type( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );

    M_mode_number=1;

    LOG(INFO) << "[CRB::offline] starting offline adaptive loop\n";

    bool reuse_prec = this->vm()["crb.reuse-prec"].template as<bool>() ;

    sampling_ptrtype Sampling;
    int sampling_size=no_residual_index+1;

    LOG(INFO) << "[CRB::offline] strategy "<< M_error_type <<"\n";
    if( proc_number == 0 ) std::cout << "[CRB::offline] strategy "<< M_error_type <<"\n";

    while ( maxerror > M_tolerance && M_N < M_iter_max && no_residual_index<sampling_size )
    {

        boost::timer timer, timer2;
        LOG(INFO) <<"========================================"<<"\n";

        if( proc_number == 0 )
        {
            std::cout << "============================================================"<<std::endl;
            std::cout << "N=" << M_N << "/"  << M_iter_max <<" ( nb proc : "<<worldComm().globalSize()<<")"<<std::endl;
        }
        LOG(INFO) << "N=" << M_N << "/"  << M_iter_max << "\n";

        // for a given parameter \p mu assemble the left and right hand side

        u->setName( ( boost::format( "fem-primal-N%1%-proc%2%" ) % (M_N)  % proc_number ).str() );

        mu.check();
        u->zero();

//        boost::tie( boost::tuples::ignore, A, F, boost::tuples::ignore ) = M_model->update( mu , 1e30 );
        
        
        
        LOG(INFO) << "  -- updated model for parameter in " << timer2.elapsed() << "s\n";
        timer2.restart();
        
        
        LOG(INFO) << "[CRB::offline] solving primal" << "\n";
        M_model->solve( mu , *u );
        
        if(proc_number==0) std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n";
        timer2.restart();
        
        for ( size_type l = 0; l < M_model->Nl(); ++l )
        {
            F[l]->close();
            element_ptrtype eltF( new element_type( M_model->functionSpace() ) );
            *eltF = *F[l];
            LOG(INFO) << "u^T F[" << l << "]= " << inner_product( *u, *eltF ) << "\n";
        }


        LOG(INFO) << "[CRB::offline] energy = " << A->energy( *u, *u ) << "\n";

        M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        bool norm_zero = false;

        M_WN.push_back( *u );

	    number_of_added_elements=1;
        
        M_N+=number_of_added_elements;

        if ( orthonormalize_primal )
        {
            orthonormalize( M_N, M_WN, number_of_added_elements );
            orthonormalize( M_N, M_WN, number_of_added_elements );
            orthonormalize( M_N, M_WN, number_of_added_elements );
        }


        LOG(INFO) << "[CRB::offline] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";

        for  (size_type q = 0; q < M_model->Qa_tril(); ++q )
        {
            Aqm_tril[q][M_N] = M_model->compute_trilinear_form( *M_WN );
        }

        for  (size_type q = 0; q < M_model->Qa(); ++q )
        {
            M_Aqm_pr[q][0].conservativeResize( M_N, M_N );
            
            // only compute the last line and last column of reduced matrices
            for ( size_type i = M_N-number_of_added_elements; i < M_N; i++ )
            {
                for ( size_type j = 0; j < M_N; ++j )
                {
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }
            
            for ( size_type j=M_N-number_of_added_elements; j < M_N; j++ )
            {
                for ( size_type i = 0; i < M_N; ++i )
                {
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute Mq_pr, Mq_du, Mq_pr_du" << "\n";

        for ( size_type q = 0; q < M_model->Qm(); ++q )
        {
            M_Mqm_pr[q][0].conservativeResize( M_N, M_N );
            
            // only compute the last line and last column of reduced matrices
            for ( size_type i=M_N-number_of_added_elements ; i < M_N; i++ )
            {
                for ( size_type j = 0; j < M_N; ++j )
                {
                    M_Mqm_pr[q][0]( i, j ) = Mqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }
            for ( size_type j = M_N-number_of_added_elements; j < M_N ; j++ )
            {
                for ( size_type i = 0; i < M_N; ++i )
                {
                    M_Mqm_pr[q][0]( i, j ) = Mqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute MFqm" << "\n";

        for ( size_type q = 0; q < M_model->Qmf(); ++q )
        {
            M_MFqm_pr[q][0].conservativeResize( M_N );
            for ( size_type j = 0; j < M_N; ++j )
            {
                MFqm[q][0]->close();
                element_ptrtype eltMF( new element_type( M_model->functionSpace() ) );
                *eltMF = *MFqm[q][0];
                M_MFqm_pr[q][0]( j ) = inner_product( *eltMF , M_WN[j] );
            }
        }


        LOG(INFO) << "[CRB::offline] compute Fq_pr, Fq_du" << "\n";

        for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
        {
                M_Fqm_pr[q][0].conservativeResize( M_N );

                for ( size_type l = 1; l <= number_of_added_elements; ++l )
                {
                    int index = M_N-l;
                    M_Fqm_pr[q][0]( index ) = M_model->Fqm( 0, q, 0, M_WN[index] );
                }
        }//loop over q

        LOG(INFO) << "compute coefficients needed for the initialization of unknown in the online step\n";

        element_ptrtype primal_initial_field ( new element_type ( M_model->functionSpace() ) );
        element_ptrtype projection    ( new element_type ( M_model->functionSpace() ) );
        M_model->initializationField( primal_initial_field, mu ); //fill initial_field

        timer2.restart();

        M_compute_variance = this->vm()["crb.compute-variance"].template as<bool>();
        if ( M_database_contains_variance_info )
            throw std::logic_error( "[CRB::offline] ERROR : build variance is not actived" );

        maxerror=M_iter_max-M_N;
        
        bool already_exist;
        do
        {
            //initialization
            already_exist=false;
            //pick randomly an element
            mu = M_Dmu->element();
            //make sure that the new mu is not already is M_WNmu
            BOOST_FOREACH( auto _mu, *M_WNmu )
            {
                if( mu == _mu )
                    already_exist=true;
            }
        }
        while( already_exist );
        
        M_current_mu = mu;

        M_rbconv.insert( convergence( M_N, boost::make_tuple(maxerror,delta_pr,delta_du) ) );

        check( M_WNmu->size() );

        if ( this->vm()["crb.check.rb"].template as<int>() == 1 )std::cout << "  -- check reduced basis done in " << timer2.elapsed() << "s\n";

        timer2.restart();
        LOG(INFO) << "time: " << timer.elapsed() << "\n";
        if( proc_number == 0 ) std::cout << "============================================================\n";
        LOG(INFO) <<"========================================"<<"\n";

        //save DB after adding an element
        this->saveDB();
    }

    for  (size_type q = 0; q < M_model->Qa_tril(); ++q )
    {
        for (int k=0 ; k<M_N; k++) 
        {
            for ( int i = 0; i < M_N; ++i )
            {
                for ( int j = 0; j < M_N; ++j )
                {
                    M_Aqm_tril_pr[q][k]( i, j ) = Aqm_tril[q][k]->energy( M_WN[i], M_WN[j] );
                }
            }
        }
    }
    
    if( proc_number == 0 )
        std::cout<<"number of elements in the reduced basis : "<<M_N<<" ( nb proc : "<<worldComm().globalSize()<<")"<<std::endl;
    LOG(INFO) << " index choosen : ";
    BOOST_FOREACH( auto id, M_index )
    LOG(INFO)<<id<<" ";
    LOG(INFO)<<"\n";
    bool visualize_basis = this->vm()["crb.visualize-basis"].template as<bool>() ;

    if ( visualize_basis )
    {
        std::vector<wn_type> wn;
        std::vector<std::string> names;
        wn.push_back( M_WN );
        names.push_back( "primal" );
        exportBasisFunctions( boost::make_tuple( wn ,names ) );

        if ( orthonormalize_primal )
        {
            std::cout<<"[CRB::offline] Basis functions have been exported but warning elements have been orthonormalized"<<std::endl;
        }
    }

    if( proc_number == 0 ) std::cout << "Offline CRB is done\n";

    return M_rbconv;

}


template<typename TruthModelType>
boost::tuple<double,double>
CRB_TRILINEAR<TruthModelType>::lb( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu,  std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold,int K  ) const
{
    google::FlushLogFiles(google::GLOG_INFO);
        
    bool save_output_behavior = this->vm()["crb.save-output-behavior"].template as<bool>();
    
    //if K>0 then the time at which we want to evaluate output is defined by
    //time_for_output = K * time_step
    //else it's the default value and in this case we take final time
    double time_for_output = 1e30;
    
    double time_step = 1e30;
    int model_K = 1; //only one time step
    size_type Qm = 0;
    
    if ( N > M_N ) N = M_N;
    
    uN.resize( model_K );        
    
    int index=0;
    BOOST_FOREACH( auto elem, uN )
    {
        uN[index].resize( N );
        index++;
    }
    
    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    beta_vector_type betaMFqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;
    
    matrixN_type A ( ( int )N, ( int )N ) ;
    vectorN_type F ( ( int )N );
    
    //-- end of initialization step
    
    //vector containing outputs from time=time_step until time=time_for_output
    std::vector<double>output_time_vector;
    output_time_vector.resize( model_K );
    double output;
    int time_index=0;
    
    // init by 1, the model could provide better init
    uN[0].setOnes(N);
    
    boost::tie( betaMqm, betaAqm, betaFqm, betaMFqm ) = M_model->computeBetaQm( this->expansion( uN[0] , N ), mu ,0 );
    
    //LOG(INFO) << "betaMFqm = " << betaMFqm[0][0] <<"\n";//<< "," << betaMFqm[1][0] << "\n";
    //LOG(INFO) << "betaMqm = " << betaMqm[0][0] << "\n";
    //LOG(INFO) << "Qm = " << M_model->Qm() << "\n";
    //LOG(INFO) << "mMaxM = " << M_model->mMaxM(0) << "\n";
    
    google::FlushLogFiles(google::GLOG_INFO);
    
    /*
     ------> Here add Newton Method and continuous Gr,Pr
     */
    using namespace vf;
    Feel::ParameterSpace<2>::Element current_mu( mu );
    
    backend_ptrtype backend_primal_problem = backend_type::build( BACKEND_PETSC );
    
    vector_ptrtype R( backend_primal_problem->newVector( Xh ) );
    sparse_matrix_ptrtype J( backend_primal_problem->newMatrix( Xh,Xh ) );
    
    double gr = mu( 0 );
    double pr = mu( 1 );
    
    int Nmax=std::max( 1.0,std::max( std::ceil( std::log( gr ) ),std::ceil( std::log( pr )-std::log( 1.e-2 ) ) ) );
    
    for ( int i = 0; i < Nmax; ++i )
    {
        int denom = ( Nmax==1 )?1:Nmax-1;
        double current_Grashofs = math::exp( math::log( 1. )+i*( math::log( gr )-math::log( 1. ) )/denom );
        double current_Prandtl = math::exp( math::log( 1.e-2 )+i*( math::log( pr )-math::log( 1.e-2 ) )/denom );
        
        std::cout << "i/N = " << i+1 << "/" << Nmax <<std::endl;
        std::cout << " intermediary Grashof = " << current_Grashofs<<std::endl;
        std::cout<< " and Prandtl = " << current_Prandtl << "\n"<<std::endl;
        
        current_mu << current_Grashofs, current_Prandtl;
        
        M_model->computeBetaQm( current_mu );
        this->updateLinearTerms( current_mu );

        backend_primal_problem->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );
        backend_primal_problem->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );

        backend_primal_problem->nlSolve(_jacobian=J , _solution=T , _residual=R);
    }
    
    
    LOG(INFO) << "lb: solve with Newton done\n";
    
    vectorN_type previous_uN( M_N );
    
    google::FlushLogFiles(google::GLOG_INFO);
        
    
    if ( save_output_behavior )
    {
        time_index=0;
        std::ofstream file_output;
        std::string mu_str;
        
        for ( int i=0; i<mu.size(); i++ )
        {
            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
        }
        
        std::string name = "output_evolution" + mu_str;
        file_output.open( name.c_str(),std::ios::out );
        
        for ( double time=time_step; time<=time_for_output; time+=time_step )
        {
            file_output<<time<<"\t"<<output_time_vector[time_index]<<"\n";
            time_index++;
        }
        
        file_output.close();
    }
    
    int size=output_time_vector.size();
    return boost::make_tuple( output_time_vector[size-1], condition_number);
        
}
    
void CRB_TRILINEAR::updateLinearTerms( parameter_type const& mu, const int N )
{    
    LOG(INFO) << "compute reduce matrices\n";
    google::FlushLogFiles(google::GLOG_INFO);

    M_bilinear_terms.setZero( N,N );
    
    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        M_bilinear_terms += betaAqm[q][0]*M_Aqm_pr[q][0].block( 0,0,N,N );
    }
    M_bilinear_terms.close();
    
    M_linear_terms.setZero( N );
    
    for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
    {
        M_linear_terms += betaFqm[0][q][0]*M_Fqm_pr[q][0].head( N );
    }
    M_linear_terms.close();
}
    
void CRB_TRILINEAR ::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    J->zero();
    J->addMatrix(1.,M_bilinear_terms);
    matrixN_type M_NL;
    M_NL.setZero( N,N);

    for ( size_type q = 0; q < M_model->Qa_tril(); ++q )
    {
        
        for (int k = 0 ; k < N; ++k) 
        {
            for ( int i = 0; i < N; ++i )
            {
                M_NL( k, i ) += (M_Aqm_tril_pr[q][k].row(i)).dot(X);
                M_NL( k, i ) += ((M_Aqm_tril_pr[q][k].col(i)).transpose()).dot(X);
            }
        }
    }
    
    J->addMatrix(1.,M_NL);
    
}
    
void CRB_TRILINEAR ::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    R->zero();
    R = M_linear_terms;
    matrixN_type M_R;
    M_R.setZero( N, N);
    
    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        R += betaAqm[q][0]*M_Aqm_pr[q][0].block( 0,0,N,N ).dot(X);
    }
    for ( size_type q = 0; q < M_model->Qa_tril(); ++q )
    {
        
        for (int k = 0 ; k < N; ++k) 
        {
            for ( int i = 0; i < N; ++i )
            {
                M_R( k, i ) += X.dot(M_Aqm_tril_pr[q][k].row(i));
            }
        }
    }
    R += M_R.dot(R);
        
}


} // Feel

#endif /* __CRB_TRILINEAR_H */
