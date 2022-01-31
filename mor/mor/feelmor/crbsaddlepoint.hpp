/* -*- mode: c++ -*-

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2009-11-24

 Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)
 Copyright (C) 2011-present Feel++ Consortium
 
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
 */

#ifndef __CRBSADDLEPOINT_H
#define __CRBSADDLEPOINT_H 1

#include <feel/feelcore/environment.hpp>
#include <feel/feelcrb/crbblock.hpp>

namespace Feel
{
po::options_description crbSaddlePointOptions( std::string const& prefix="", int const& n_block=2 );

/**
 * \class CRBSaddlePoint
 * \brief Certfified Reduced BAsis for Saddle Point Problems
 *
 * \author JB Wahl
 */
template<typename TruthModelType>
class CRBSaddlePoint :
        public CRBBlock<TruthModelType>
{
    typedef CRBBlock<TruthModelType> super;


public:
    //@{ // Truth Model
    typedef double value_type;
    typedef TruthModelType model_type;
    typedef std::shared_ptr<model_type> truth_model_ptrtype;
    typedef typename model_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //@}

    //@{ /// Parameter Space
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    //@}

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;
    typedef typename convergence_type::value_type convergence;

    using self_type = CRBSaddlePoint;
    using self_ptrtype = std::shared_ptr<self_type>;

    //@{ /// Function Space and Elements
    typedef typename model_type::space_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;
    //@}

    //@{ Backend and Matrix
    typedef typename model_type::backend_type backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;
    //@}

    //@{ /// Eigen Objects
    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic, 1> > map_dense_vector_type;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant
    //@}

    //@{ /// Exporter
    typedef Exporter<mesh_type> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;
    //@}

    //@{ /// Database
    typedef CRBElementsDB<model_type> crb_elements_db_type;
    typedef std::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;
    //@}

    typedef std::vector< std::vector< std::vector< std::vector< matrixN_type >>>> blockmatrixN_type;
    typedef std::vector< std::vector< std::vector< vectorN_type >>> blockvectorN_type;
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
    CRBSaddlePoint( std::string const& name = "defaultname_crb",
                    crb::stage stage = crb::stage::online,
                    WorldComm const& worldComm = Environment::worldComm() ) :
        CRBSaddlePoint( name, std::make_shared<model_type>(stage), stage )
        {}

    //! constructor from command line options
    CRBSaddlePoint( std::string const& name, truth_model_ptrtype const & model,
                    crb::stage stage = crb::stage::online, std::string const& prefixExt = "" ) :
        super( name, model, stage, prefixExt )
        {
        }

public:
    void init()
        {
            using Feel::cout;
        if ( !this->M_rebuild && this->loadDB() )
        {
            cout << "Database CRB SP " << this->lookForDB() << " available and loaded with M_N0="
                 << subN(0) << ", M_N1="<< subN(1) <<", M_N="<<this->M_N <<std::endl;
            if( this->M_loadElementsDb )
            {
                if( this->M_elements_database.loadDB() )
                {
                    auto size0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>()->primalRB().size();
                    auto size1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>()->primalRB().size();
                    cout<<"Database for basis functions " << this->M_elements_database.lookForDB() << " available and loaded with\n"
                        << size0 <<" primal basis functions in RBSpace0 and "
                        << size1 << " primal basis functions in RBSpace0\n";
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
    }

    element_type solve( parameter_type const& mu )
    {
        return this->M_model->solve( mu );
    }

    void offlineResidual( int Ncur, int number_of_added_elements=1 ) override;

    max_error_type maxErrorBounds( size_type N ) const override;

    matrix_info_tuple fixedPointPrimal( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold, std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput=true ) const override;

    void fixedPointDual(  size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN, std::vector< vectorN_type > & uNdu,  std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K=0) const override;

    double correctionTerms(parameter_type const& mu, std::vector< vectorN_type > const & uN, std::vector< vectorN_type > const & uNdu,  std::vector<vectorN_type> const & /*uNold*/, int const k ) const override;


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


private :
    void initRezMatrix() override;

    template <int Row>
    void offlineResidualSP( int Ncur , int number_of_added_elements );

    double onlineResidual( int Ncur, parameter_type const& mu, vectorN_type Un ) const;
    template <int Row>
    double onlineResidualSP( int Ncur, parameter_type const& mu, vectorN_type Un, bool test=false ) const;
    void testResidual() override;
    double empiricalError( int N, parameter_type const& mu, std::vector<double> output_vec ) const;



    std::vector< std::vector< std::vector< std::vector< std::vector< double >>>>> M_R_RhsRhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< vectorN_type >>>>> M_R_Lhs0Rhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< vectorN_type >>>>> M_R_Lhs1Rhs;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs0Lhs0;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs0Lhs1;
    std::vector< std::vector< std::vector< std::vector< std::vector< matrixN_type >>>>> M_R_Lhs1Lhs1;

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    using super::subN;
    using super::M_blockAqm_pr;
    using super::M_blockAqm_du;
    using super::M_blockAqm_pr_du;
    using super::M_blockFqm_pr;
    using super::M_blockFqm_du;
    using super::M_blockLqm_pr;
    using super::M_blockLqm_du;
    using super::notEmptyAqm;
    using super::notEmptyFqm;
    using super::notEmptyLqm;
}; // class CRBSaddlePoint


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::initRezMatrix()
{
    int QLhs = this->M_model->Qa();
    int QRhs = this->M_model->Ql( 0 );

    if ( M_R_RhsRhs.size()==0 )
    {
        M_R_RhsRhs.resize(2);
        M_R_Lhs0Rhs.resize(2);
        M_R_Lhs1Rhs.resize(2);
        M_R_Lhs0Lhs0.resize(2);
        M_R_Lhs0Lhs1.resize(2);
        M_R_Lhs1Lhs1.resize(2);
    }

    for ( int r=0; r<2; r++ )
    {
        M_R_RhsRhs[r].resize(QRhs);
        for ( int q1=0; q1<QRhs; q1++ )
        {
            int mMax1 = this->M_model->mMaxF(0,q1);
            M_R_RhsRhs[r][q1].resize(mMax1);
            for ( int m1=0; m1<mMax1; m1++ )
            {
                M_R_RhsRhs[r][q1][m1].resize(QRhs);
                for ( int q2=0; q2<QRhs; q2++ )
                    M_R_RhsRhs[r][q1][m1][q2].resize( this->M_model->mMaxF(0,q2) );
            }
        }

        M_R_Lhs0Lhs0[r].resize(QLhs);
        M_R_Lhs0Lhs1[r].resize(QLhs);
        M_R_Lhs0Rhs[r].resize(QLhs);
        for ( int q1=0; q1<QLhs; q1++ )
        {
            int mMax1 = this->M_model->mMaxA(q1);
            M_R_Lhs0Lhs0[r][q1].resize(mMax1);
            M_R_Lhs0Lhs1[r][q1].resize(mMax1);
            M_R_Lhs0Rhs[r][q1].resize(mMax1);
            for ( int m1=0; m1<mMax1; m1++ )
            {
                M_R_Lhs0Lhs0[r][q1][m1].resize(QLhs);
                M_R_Lhs0Lhs1[r][q1][m1].resize(QLhs);
                M_R_Lhs0Rhs[r][q1][m1].resize(QRhs);
                for ( int q2=0; q2<QLhs; q2++ )
                    M_R_Lhs0Lhs0[r][q1][m1][q2].resize( this->M_model->mMaxA(q2) );
                for ( int q2=0; q2<QLhs; q2++ )
                    M_R_Lhs0Lhs1[r][q1][m1][q2].resize( this->M_model->mMaxA(q2) );
                for ( int q2=0; q2<QRhs; q2++ )
                    M_R_Lhs0Rhs[r][q1][m1][q2].resize( this->M_model->mMaxF(0,q2) );
            }
        }

        M_R_Lhs1Lhs1[r].resize(QLhs);
        M_R_Lhs1Rhs[r].resize(QLhs);
        for ( int q1=0; q1<QLhs; q1++ )
        {
            int mMax1 = this->M_model->mMaxA(q1);
            M_R_Lhs1Lhs1[r][q1].resize(mMax1);
            M_R_Lhs1Rhs[r][q1].resize(mMax1);

            for ( int m1=0; m1<mMax1; m1++ )
            {
                M_R_Lhs1Lhs1[r][q1][m1].resize(QLhs);
                M_R_Lhs1Rhs[r][q1][m1].resize(QRhs);
                for ( int q2=0; q2<QLhs; q2++ )
                    M_R_Lhs1Lhs1[r][q1][m1][q2].resize( this->M_model->mMaxA(q2) );
                for ( int q2=0; q2<QRhs; q2++ )
                    M_R_Lhs1Rhs[r][q1][m1][q2].resize( this->M_model->mMaxF(0,q2) );
            }
        }
    }
} // updateAffinedecompositionsize()





template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::matrix_info_tuple
CRBSaddlePoint<TruthModelType>::fixedPointPrimal(  size_type N, parameter_type const& mu,
                                                   std::vector< vectorN_type > & uN,
                                                   std::vector<vectorN_type> &,
                                                   std::vector< double > & output_vector,
                                                   int K, bool print_rb_matrix,
                                                   bool computeOutput ) const
{
    bool is_linear = this->M_model->isLinear();
    double output=0;
    double increment = this->M_fixedpointIncrementTol;
    int N0 = this->subN( 0,N );
    int N1 = this->subN( 1,N );

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;

    int Qa=this->M_model->Qa();
    int Ql=this->M_model->Ql(this->M_output_index);
    int Qf=this->M_model->Ql(0);
    std::vector<int> mMaxA(Qa);
    std::vector<int> mMaxL( Ql );
    std::vector<int> mMaxF( Qf );
    for ( size_type q = 0; q < Qa; ++q )
    {
        mMaxA[q]=this->M_model->mMaxA(q);
    }
    for ( size_type q = 0; q < Qf; ++q )
    {
        mMaxF[q]=this->M_model->mMaxF(0,q);
    }
    for ( size_type q = 0; q < Ql; ++q )
    {
        mMaxL[q]=this->M_model->mMaxF(this->M_output_index,q);
    }

    matrixN_type A( N0+N1, N0+N1 );
    vectorN_type F( N0+N1 );
    vectorN_type L( N0+N1 );

    if ( is_linear )
    {
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu );

        A.setZero( N0+N1,N0+N1);
        for ( size_type q=0; q<Qa; q++ )
        {
            for ( size_type m=0; m<mMaxA[q]; m++ )
            {
                if ( this->notEmptyAqm(0,0,q,m) )
                    A.block( 0,  0, N0, N0) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N0);
                if ( this->notEmptyAqm(0,1,q,m) )
                    A.block( 0, N0, N0, N1) += betaAqm[q][m]*M_blockAqm_pr[0][1][q][m].block(0,0,N0,N1);
                if ( this->notEmptyAqm(1,0,q,m) )
                    A.block(N0,  0, N1, N0) += betaAqm[q][m]*M_blockAqm_pr[1][0][q][m].block(0,0,N1,N0);
                if ( this->notEmptyAqm(1,1,q,m) )
                    A.block(N0, N0, N1, N1) += betaAqm[q][m]*M_blockAqm_pr[1][1][q][m].block(0,0,N1,N1);
            }
        }

        F.setZero(N0+N1);
        for ( size_type q=0; q<Qf; q++ )
        {
            for ( size_type m=0; m<mMaxF[q]; m++ )
            {
                if ( this->notEmptyFqm(0,q,m) )
                    F.head(N0) += betaFqm[0][q][m]*M_blockFqm_pr[0][q][m].head(N0);
                if ( this->notEmptyFqm(1,q,m) )
                    F.tail(N1) += betaFqm[0][q][m]*M_blockFqm_pr[1][q][m].head(N1);
            }
        }
        //Feel::cout << "A=\n"<<A<<"\n F=\n"<<F<<std::endl;

        uN[0] = A.fullPivLu().solve( F );
    }
    else //non-linear
    {
        vectorN_type previous_uN( N0+N1 );
        uN[0].setZero( N0+N1 ); // initial guess
        int fi=0;
        bool fixPointIsFinished = false;
        do
        {
            previous_uN = uN[0];
            if ( this->M_useRbSpaceContextEim && this->M_hasRbSpaceContextEim )
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                    this->M_model->computeBetaQm( uN[0], mu/*, N*/ );
            else
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                    this->M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
            A.setZero( N0+N1,N0+N1 );
            for ( size_type q=0; q<Qa; q++ )
            {
                for ( size_type m=0; m<mMaxA[q]; m++ )
                {
                    if ( this->notEmptyAqm(0,0,q,m) )
                        A.block( 0,  0, N0, N0) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N0);
                    if ( this->notEmptyAqm(0,1,q,m) )
                        A.block( 0, N0, N0, N1) += betaAqm[q][m]*M_blockAqm_pr[0][1][q][m].block(0,0,N0,N1);
                    if ( this->notEmptyAqm(1,0,q,m) )
                        A.block(N0,  0, N1, N0) += betaAqm[q][m]*M_blockAqm_pr[1][0][q][m].block(0,0,N1,N0);
                    if ( this->notEmptyAqm(1,1,q,m) )
                        A.block(N0, N0, N1, N1) += betaAqm[q][m]*M_blockAqm_pr[1][1][q][m].block(0,0,N1,N1);
                }
            }

            F.setZero(N0+N1);
            for ( size_type q=0; q<Qf; q++ )
            {
                for ( size_type m=0; m<mMaxF[q]; m++ )
                {
                    if ( this->notEmptyFqm(0,q,m) )
                        F.head(N0) += betaFqm[0][q][m]*M_blockFqm_pr[0][q][m].head(N0);
                    if ( this->notEmptyFqm(1,q,m) )
                        F.tail(N1) += betaFqm[0][q][m]*M_blockFqm_pr[1][q][m].head(N1);
                }
            }
            uN[0] = A.fullPivLu().solve( F );

            increment = (uN[0]-previous_uN).norm();
            auto increment_abs = (uN[0]-previous_uN).array().abs();
            fixPointIsFinished = increment < this->M_fixedpointIncrementTol || fi>=this->M_fixedpointMaxIterations;
            this->online_iterations_summary.first = fi;
            this->online_iterations_summary.second = increment;
            if( this->M_fixedpointVerbose  && this->worldComm().isMasterRank() )
            {
                DVLOG(2) << "iteration " << fi << " increment error: " << increment << "\n";
                VLOG(2)<<"[CRB::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment <<std::endl;
                double residual_norm = (A * uN[0] - F).norm() ;
                VLOG(2) << " residual_norm :  "<<residual_norm;
                Feel::cout << "[CRBSaddlePoint::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment << std::endl;
            }
            ++fi;
        }while ( !fixPointIsFinished );
    }


    if ( computeOutput )
    {
        L.setZero(N0+N1);
        for ( size_type q=0; q<Ql; q++ )
        {
            for ( size_type m=0; m<mMaxL[q]; m++ )
            {
                if ( this->notEmptyLqm(0,q,m) )
                    L.head(N0) += betaFqm[this->M_output_index][q][m]*M_blockLqm_pr[0][q][m].head(N0);
                if ( this->notEmptyLqm(1,q,m) )
                    L.tail(N1) += betaFqm[this->M_output_index][q][m]*M_blockLqm_pr[1][q][m].head(N1);
            }
        }
        output = L.dot( uN[0] );
        output_vector[0] = output;
    }


    double condition_number = 0;
    double determinant = 0;
    if( this->M_computeMatrixInfo )
    {
        condition_number = this->computeConditioning( A );
        determinant = A.determinant();
    }

    auto matrix_info = boost::make_tuple(condition_number,determinant);

    if( print_rb_matrix && !this->M_offline_step )
        this->printRBMatrix( A,mu );
    return matrix_info;
} //fixedPointPrimal()

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::fixedPointDual(  size_type N, parameter_type const& mu,
                                                 std::vector< vectorN_type > const& uN,
                                                 std::vector<vectorN_type> & uNdu,
                                                 std::vector<vectorN_type> & uNduold,
                                                 std::vector< double > & output_vector,
                                                 int K ) const
{
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);

    bool is_linear = this->M_model->isLinear();

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;

    matrixN_type Adu ( N0+N1, N0+N1 ) ;
    vectorN_type Fdu ( N0+N1 );

    int Qa=this->M_model->Qa();
    int Ql=this->M_model->Ql(this->M_output_index);

    std::vector<int> mMaxA(Qa);
    std::vector<int> mMaxF( Ql );
    for ( size_type q = 0; q < Qa; ++q )
    {
        mMaxA[q]=this->M_model->mMaxA(q);
    }
    for ( size_type q = 0; q < Ql; ++q )
    {
        mMaxF[q]=this->M_model->mMaxF(0,q);
    }

    if( is_linear )
    {
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu );
        Adu.setZero( N,N );
        for ( size_type q=0; q<Qa; q++ )
        {
            for ( size_type m=0; m<mMaxA[q]; m++ )
            {
                if ( this->notEmptyAqm(0,0,q,m) )
                    Adu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_du[0][0][q][m].block(0,0,N0,N0);
                if ( this->notEmptyAqm(0,1,q,m) )
                    Adu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_du[0][1][q][m].block(0,0,N0,N1);
                if ( this->notEmptyAqm(1,0,q,m) )
                    Adu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_du[1][0][q][m].block(0,0,N1,N0);
                if ( this->notEmptyAqm(1,1,q,m) )
                    Adu.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_du[1][1][q][m].block(0,0,N1,N1);
            }
        }
        Fdu.setZero( N0+N1 );
        for ( size_type q = 0; q < Ql ; ++q )
            {
                for(int m=0; m < mMaxF[q]; m++)
                {
                    if ( this->notEmptyFqm(0,q,m) )
                        Fdu.head(N0) -= betaFqm[this->M_output_index][q][m]*M_blockLqm_du[0][q][m].head( N0 );
                    if ( this->notEmptyFqm(1,q,m) )
                        Fdu.tail(N1) -= betaFqm[this->M_output_index][q][m]*M_blockLqm_du[1][q][m].head( N1 );
                }
            }
            uNdu[0] = Adu.fullPivLu().solve( Fdu );
    }
    else // non linear
    {
        double increment = this->M_fixedpointIncrementTol;
        vectorN_type next_uNdu( N0+N1 );
        uNdu[0].setZero( N0+N1 );
        int fi=0;
        do
        {
            // backup uNdu
            next_uNdu = uNdu[0];
            // update coefficients of affine decomposition
            // if ( this->M_useRbSpaceContextEim && this->M_hasRbSpaceContextEim )
            //     boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( uN[0], mu/*, N*/ );
            // else
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                    this->M_model->computeBetaQm( this->expansion( uN[0], N )/*dualRB*/, mu );
            // assemble rb matrix
            Adu.setZero( N0+N1, N0+N1 );
            for ( size_type q=0; q<Qa; q++ )
            {
                for ( size_type m=0; m<mMaxA[q]; m++ )
                {
                    if ( this->notEmptyAqm(0,0,q,m) )
                        Adu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_du[0][0][q][m].block(0,0,N0,N0);
                    if ( this->notEmptyAqm(0,1,q,m) )
                        Adu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_du[0][1][q][m].block(0,0,N0,N1);
                    if ( this->notEmptyAqm(1,0,q,m) )
                        Adu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_du[1][0][q][m].block(0,0,N1,N0);
                    if ( this->notEmptyAqm(1,1,q,m) )
                        Adu.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_du[1][1][q][m].block(0,0,N1,N1);
                }
            }

            Fdu.setZero( N0+N1 );
            for ( size_type q=0; q<Ql; q++ )
            {
                for ( size_type m=0; m<mMaxF[q]; m++ )
                {
                    if ( this->notEmptyFqm(0,q,m) )
                        Fdu.head(N0) += betaFqm[0][q][m]*M_blockFqm_du[0][q][m].head(N0);
                    if ( this->notEmptyFqm(1,q,m) )
                        Fdu.tail(N1) += betaFqm[0][q][m]*M_blockFqm_du[1][q][m].head(N1);
                }
            }
            uNdu[0] = Adu.fullPivLu().solve( Fdu );

            increment = (uNdu[0]-next_uNdu).norm();
            if( this->M_fixedpointVerbose  && this->worldComm().isMasterRank() )
            {
                VLOG(2)<<"[CRBSaddlePoint::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment;
                std::cout<<"[CRBSaddlePoint::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment << "\n";
            }
            fi++;
        } while ( increment > this->M_fixedpointIncrementTol && fi<this->M_fixedpointMaxIterations );
    } // if non linear
} // fixedPointDual

template< typename TruthModelType>
double
CRBSaddlePoint<TruthModelType>::correctionTerms(parameter_type const& mu, std::vector< vectorN_type > const & uN, std::vector< vectorN_type > const & uNdu,  std::vector<vectorN_type> const & /*uNold*/, int const k ) const
{
    int N=0;
    for ( int n=0; n<this->WNmuSize(); n++ )
        N = uN[0].size()==this->dimension(n) ? n:0;
    int N0 = this->subN(0,N);
    int N1 = this->subN(1,N);

    matrixN_type Aprdu ( N0+N1, N0+N1 ) ;
    vectorN_type Fdu ( N0+N1 );
    vectorN_type du ( N0+N1 );
    vectorN_type pr ( N0+N1 );

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;

    bool is_linear = this->M_model->isLinear();

    double correction=0;

    Aprdu.setZero( N0+N1 , N0+N1 );
    Fdu.setZero( N0+N1 );

    if ( is_linear )
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu/* ,time*/);
    else
    {
         if ( this->M_useRbSpaceContextEim && this->M_hasRbSpaceContextEim )
             boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( uN[0], mu/*, N*/ );
         else
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                this->M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
    }
    for(size_type q = 0;q < this->M_model->Ql(0); ++q)
    {
        for(int m=0; m < this->M_model->mMaxF(0,q); m++)
        {
            if ( this->notEmptyFqm(0,q,m) )
                Fdu.head(N0) += betaFqm[0][q][m]*M_blockFqm_du[0][q][m].head(N0);
            if ( this->notEmptyFqm(1,q,m) )
                Fdu.tail(N1) += betaFqm[0][q][m]*M_blockFqm_du[1][q][m].head(N1);
        }
    }
    for(size_type q = 0;q < this->M_model->Qa(); ++q)
    {
        for(int m=0; m < this->M_model->mMaxA(q); m++)
        {
            if ( this->notEmptyAqm(0,0,q,m) )
                Aprdu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_pr_du[0][0][q][m].block(0,0,N0,N0);
            if ( this->notEmptyAqm(0,1,q,m) )
                Aprdu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_pr_du[0][1][q][m].block(0,0,N0,N1);
            if ( this->notEmptyAqm(1,0,q,m) )
                Aprdu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_pr_du[1][0][q][m].block(0,0,N1,N0);
            if ( this->notEmptyAqm(1,1,q,m) )
                Aprdu.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_pr_du[1][1][q][m].block(0,0,N1,N1);
        }
    }

    du = uNdu[0];
    pr = uN[0];
    correction = -( Fdu.dot( du ) - du.dot( Aprdu*pr )  );

    return correction;
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::offlineResidual( int Ncur , int number_of_added_elements )
{
    tic();
    offlineResidualSP<0>( Ncur, number_of_added_elements );
    toc("OfflineResidual 0");
    tic();
    offlineResidualSP<1>( Ncur, number_of_added_elements );
    toc("OfflineResidual 1");
}

template<typename TruthModelType>
template <int Row>
void
CRBSaddlePoint<TruthModelType>::offlineResidualSP( int Ncur , int number_of_added_elements )
{
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();

    bool optimize = boption(_prefix=this->M_prefix,_name="crb.optimize-offline-residual") ;

    int N0 = this->subN(0,Ncur);
    int N1 = this->subN(1,Ncur);
    int n_added0 = N0 - this->subN(0,Ncur-1);
    int n_added1 = N1 - this->subN(1,Ncur-1);
    // in the case of SER we need to rebuild for the new EIM basis
    if( ioption(_prefix=this->M_prefix,_name="ser.eim-frequency") != 0 )
    {
        n_added0 = N0;
        n_added1 = N1;
    }

    int QLhs = this->M_model->Qa();
    int QRhs = this->M_model->Ql(0);

    auto Xh0 = this->M_model->functionSpace()->template functionSpace<0>();
    auto Xh1 = this->M_model->functionSpace()->template functionSpace<1>();

    vector_ptrtype Z1 = backend()->newVector( this->M_model->functionSpace()->template functionSpace<Row>() );
    vector_ptrtype Z2 = backend()->newVector( this->M_model->functionSpace()->template functionSpace<Row>() );
    vector_ptrtype X0 =backend()->newVector( this->M_model->functionSpace()->template functionSpace<0>() );
    vector_ptrtype X1 =backend()->newVector( this->M_model->functionSpace()->template functionSpace<1>() );
    vector_ptrtype W = backend()->newVector( this->M_model->functionSpace()->template functionSpace<Row>() );

    tic();
    this->M_model->initBlockMatrix();
    std::vector< std::vector< sparse_matrix_ptrtype >> Lhs0 = this->M_model->AqmBlock( Row, 0 );
    std::vector< std::vector< sparse_matrix_ptrtype >> Lhs1 = this->M_model->AqmBlock( Row, 1 );
    std::vector< std::vector< vector_ptrtype >> Fqm = this->M_model->FqmBlock( 0, Row );
    toc("initBlockmatrix");

    if ( N0==n_added0 )
    {
        if ( Row==0 && Environment::isMasterRank() )
        {
            boost::filesystem::path dir("Riesz");
            if ((boost::filesystem::exists(dir)))
                boost::filesystem::remove_all(dir);
            boost::filesystem::create_directory(dir);
        }

        tic();
        for ( int q=0; q<QRhs; q++ )
        {
            for ( int m=0; m< this->M_model->mMaxF(0,q); m++ )
            {
                this->M_model->l2solve( Z1, Fqm[q][m], Row );
                Z1->save( (boost::format("Riesz/%1%_%2%%3%") %Row %q %m).str()  );
            }
        }

        for ( int q1=0; q1<QRhs; q1++ )
        {
            for ( int m1=0; m1< this->M_model->mMaxF(0,q1); m1++ )
            {
                Z1->load( (boost::format("Riesz/%1%_%2%%3%") %Row %q1 %m1).str()  );
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->mMaxF(0,q2); m2++ )
                    {
                        Z2->load( (boost::format("Riesz/%1%_%2%%3%") %Row %q2 %m2).str()  );
                        M_R_RhsRhs[Row][q1][m1][q2][m2] = this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop
                } // q2 loop
            } // m1 loop
        } // q1 loop
        toc("rhsXrhs first update");
    } //N==M_Nm

    tic();
    for ( int i=N0-n_added0; i<N0; i++ )
    {
        *X0 = XN0->primalBasisElement( i );
        for ( int q=0; q<QLhs; q++ )
        {
            for ( int m=0; m<this->M_model->mMaxA(q); m++ )
            {
                Lhs0[q][m]->multVector( X0, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );
                Z1->save( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q %m %i).str()  );
            }
        }
    }
    toc("Precompute Lhs0");

    tic();
    for ( int i=N1-n_added1; i<N1; i++ )
    {
        *X1 = XN1->primalBasisElement( i );
        for ( int q=0; q<QLhs; q++ )
        {
            for ( int m=0; m<this->M_model->mMaxA(q); m++ )
            {
                Lhs1[q][m]->multVector( X1, W );
                W->scale(-1.);
                this->M_model->l2solve( Z1, W, Row );
                Z1->save( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q %m %i).str()  );
            }
        }
    }
    toc("Precompute Lhs1");

    // LHS0 LOOP ON I
    for ( int i=N0-n_added0; i<N0; i++ )
    {
        for ( int q1=0; q1<QLhs; q1++ )
        {
            for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
            {
                Z1->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q1 %m1 %i ).str()  );

                // RHS LOOP
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->mMaxF(0,q2); m2++ )
                    {
                        M_R_Lhs0Rhs[Row][q1][m1][q2][m2].conservativeResize(N0);
                        Z2->load( (boost::format("Riesz/%1%_%2%%3%") %Row %q2 %m2).str() );
                        M_R_Lhs0Rhs[Row][q1][m1][q2][m2]( i ) =
                            2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop Rhs
                } // q2 loop Rhs

                // LHS0 LOOP ON J
                for ( int j=0; j<N0; j++ )
                {
                    if ( Row==0 && optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].conservativeResize( N0, N0 );
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1].conservativeResize( N0, N0 );

                                Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q2 %m2 %j).str()  );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0

                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1].conservativeResize( N0, N0 );
                        Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q1 %m1 %j).str()  );

                        double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1]( i, j ) = prod;
                        M_R_Lhs0Lhs0[Row][q1][m1][q1][m1]( j, i ) = prod;
                    } // Row==0
                    else
                    {
                        for ( int q2=0; q2<QLhs; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].conservativeResize( N0, N0 );
                                Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q2 %m2 %j).str()  );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) =
                                    this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    }
                } // j loop Lhs0

                // LHS1 LOOP ON J
                for ( int j=0; j<N1; j++ )
                {
                    for ( int q2=0; q2<QLhs; q2++ )
                    {
                        for( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                        {
                            M_R_Lhs0Lhs1[Row][q1][m1][q2][m2].conservativeResize( N0, N1 );
                            Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q2 %m2 %j).str()  );
                            M_R_Lhs0Lhs1[Row][q1][m1][q2][m2]( i, j )
                                = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                        } // m2 loop Lhs1
                    } // q2 loop Lhs1
                } // j loop Lhs1

            } // m1 loop Lhs0
        } // q1 loop Lhs0
    } // i loop Lhs0

    // LHS0 LOOP ON J
    for ( int j=N0-n_added0; j<N0; j++ )
    {
        for ( int q1=0; q1<QLhs; q1++ )
        {
            for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
            {
                Z1->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q1 %m1 %j).str()  );

                // LHS0 LOOP ON I
                for ( int i=0; i<N0; i++ )
                {
                    if ( Row==0 && optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q2 %m2 %i).str()  );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs0Lhs0[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    } // Row==0
                    else
                    {
                        for ( int q2=0; q2<QLhs; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q2 %m2 %i).str()  );
                                M_R_Lhs0Lhs0[Row][q1][m1][q2][m2]( i, j ) =
                                    this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs0
                        } // q2 loop Lhs0
                    } // ! Row==0
                } // i loop Lhs0
            } // m1 loop Lhs0
        } // q1 loop Lhs0
    } // j loop Lhs0

    // LHS1 LOOP ON I
    for ( int i=N1-n_added1; i<N1; i++ )
    {
        for ( int q1=0; q1<QLhs; q1++ )
        {
            for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
            {
                Z1->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q1 %m1 %i).str()  );

                // RHS LOOP
                for ( int q2=0; q2<QRhs; q2++ )
                {
                    for ( int m2=0; m2<this->M_model->mMaxF(0,q2); m2++ )
                    {
                        M_R_Lhs1Rhs[Row][q1][m1][q2][m2].conservativeResize(N1);
                        Z2->load( (boost::format("Riesz/%1%_%2%%3%") %Row %q2 %m2).str() );
                        M_R_Lhs1Rhs[Row][q1][m1][q2][m2]( i ) =
                            2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                    } // m2 loop Rhs
                } // q2 loop Rhs

                // LHS1 LOOP ON J
                for ( int j=0; j<N1; j++ )
                {
                    if ( Row==1 && optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].conservativeResize( N1, N1 );
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1].conservativeResize( N1, N1 );
                                Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q2 %m2 %j).str()  );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1].conservativeResize( N1, N1 );
                        Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q1 %m1 %j).str()  );
                        double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( i, j ) = prod;
                        M_R_Lhs1Lhs1[Row][q1][m1][q1][m1]( j, i ) = prod;

                    } // if Row==1
                    else
                    {
                        for ( int q2=0; q2<QLhs; q2++ )
                        {
                            for( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].conservativeResize( N1, N1 );
                                Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q2 %m2 %j).str()  );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j )
                                    = this->M_model->scalarProduct( Z1, Z2, Row );
                            }
                        }
                    } // ! Row==1
                } // j loop Lhs1
            } // m1 loop Lhs1
        } // q1 loop Lhs1
    } // i loop Lhs1


    // LHS1 LOOP ON J
    for ( int j=N1-n_added1; j<N1; j++ )
    {
        for ( int q1=0; q1<QLhs; q1++ )
        {
            for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
            {
                Z1->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q1 %m1 %j).str()  );

                // LHS0 LOOP ON I
                for ( int i=0; i<N0; i++ )
                {
                    for ( int q2=0; q2<QLhs; q2++ )
                    {
                        for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                        {
                            Z2->load( (boost::format("Riesz/%1%0_%2%%3%_%4%") %Row %q2 %m2 %i).str()  );
                            M_R_Lhs0Lhs1[Row][q2][m2][q1][m1]( i, j )
                                = 2.0*this->M_model->scalarProduct( Z1, Z2, Row );
                        } // m2 loop Lhs0
                    } // q2 loop Lhs0
                } // i loop Lhs0

                // LHS1 LOOP ON I
                for ( int i=0; i<N1; i++ )
                {
                    if ( Row==1 && optimize )
                    {
                        for ( int q2=0; q2<q1; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q2 %m2 %i).str()  );
                                double prod = this->M_model->scalarProduct( Z1, Z2, Row );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j ) = prod;
                                M_R_Lhs1Lhs1[Row][q2][m2][q1][m1]( j, i ) = prod;
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1
                    }
                    else
                    {
                        for ( int q2=0; q2<QLhs; q2++ )
                        {
                            for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                            {
                                Z2->load( (boost::format("Riesz/%1%1_%2%%3%_%4%") %Row %q2 %m2 %i).str()  );
                                M_R_Lhs1Lhs1[Row][q1][m1][q2][m2]( i, j ) =
                                    this->M_model->scalarProduct( Z1, Z2, Row );
                            } // m2 loop Lhs1
                        } // q2 loop Lhs1
                    } // !Row==1
                } // j loop Lhs1
            } // m1 loop Lhs1
        } // q1 loop Lhs1
    } // j loop Lhs1

    this->M_model->clearBlockMatrix();

}

template<typename TruthModelType>
double
CRBSaddlePoint<TruthModelType>::onlineResidual( int Ncur, parameter_type const& mu,
                                                vectorN_type Un ) const
{
    tic();
    double res0 = onlineResidualSP<0>( Ncur, mu, Un );
    double res1 = onlineResidualSP<1>( Ncur, mu, Un );
    toc("online rez",false);
    return std::sqrt( res0 + res1 );
}

template<typename TruthModelType>
template<int Row>
double
CRBSaddlePoint<TruthModelType>::onlineResidualSP( int Ncur, parameter_type const& mu,
                                                vectorN_type Un, bool test ) const
{
    using Feel::cout;
    int N0 = this->subN(0,Ncur);
    int N1 = this->subN(1,Ncur);

    CHECK( Un.size() == N0 + N1 )
        << "invalide size of Un, vector can't be cut, Un.size="
        <<Un.size()<<", N0="<<N0<<", N1="<<N1<<std::endl;

    vectorN_type un = Un.head( N0 );
    vectorN_type pn = Un.tail( N1 );

    int QRhs = this->M_model->Ql( 0 );
    int QLhs0 = this->M_model->Qa();
    int QLhs1 = this->M_model->Qa();

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;
    if( this->M_model->isLinear() )
    {
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
            this->M_model->computeBetaQm( mu );
    }
    else
    {
        if ( this->M_useRbSpaceContextEim && this->M_hasRbSpaceContextEim )
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                this->M_model->computeBetaQm( Un, mu );
        else
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                this->M_model->computeBetaQm( this->expansion( Un, Ncur ), mu );
    }

    beta_vector_type betaLhs = betaAqm;
    beta_vector_type betaRhs= betaFqm[0];

    double RhsRhs = 0;
    double Lhs0Rhs = 0;
    double Lhs0Lhs0 = 0;
    double Lhs0Lhs1 =0;
    double Lhs1Rhs = 0;
    double Lhs1Lhs1 = 0;

    for ( int q1=0; q1<QRhs; q1++ )
    {
        for ( int m1=0; m1<this->M_model->mMaxF( 0, q1 ); m1++ )
        {
            double beta_q1 = betaRhs[q1][m1];
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxF( 0, q2 ); m2++ )
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
        for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
        {
            double beta_q1 = betaLhs[q1][m1];

            // SUBLOOP ON RHS
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxF( 0, q2 ); m2++ )
                {
                    double beta_q2 = betaRhs[q2][m2];
                    Lhs0Rhs += beta_q1*beta_q2*M_R_Lhs0Rhs[Row][q1][m1][q2][m2].head(N0).dot(un);
                } // m2 loop rhs
            } // q2 loop rhs

            // SUBLOOP ON LHS0
            for ( int q2=0; q2<QLhs0; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                {
                    double beta_q2 = betaLhs[q2][m2];

                    auto m = M_R_Lhs0Lhs0[Row][q1][m1][q2][m2].block(0,0,N0,N0)*un;
                    Lhs0Lhs0 += beta_q1*beta_q2*un.dot(m);
                } // m2 loop lhs0
            } // q2 loop lhs0

            // SUBLOOP ON LHS1
            for ( int q2=0; q2<QLhs1; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                {
                    double beta_q2 = betaLhs[q2][m2];
                    auto m = M_R_Lhs0Lhs1[Row][q1][m1][q2][m2].block(0,0,N0,N1)*pn;
                    Lhs0Lhs1 += beta_q1*beta_q2*un.dot(m);
                } // m2 loop on lhs1
            } // q2 loop on lhs1

        } // m1 loop lhs0
    } // q1 loop lhs0

    // LOOP ON LHS1
    for ( int q1=0; q1<QLhs1; q1++ )
    {
        for ( int m1=0; m1<this->M_model->mMaxA(q1); m1++ )
        {
            double beta_q1 = betaLhs[q1][m1];

            // SUBLOOP ON RHS
            for ( int q2=0; q2<QRhs; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxF( 0, q2 ); m2++ )
                {
                    double beta_q2 = betaRhs[q2][m2];
                    Lhs1Rhs += beta_q1*beta_q2*M_R_Lhs1Rhs[Row][q1][m1][q2][m2].head(N1).dot(pn);
                } // m2 loop rhs
            } // q2 loop rhs

            // SUBLOOP ON LHS1
            for ( int q2=0; q2<QLhs1; q2++ )
            {
                for ( int m2=0; m2<this->M_model->mMaxA(q2); m2++ )
                {
                    double beta_q2 = betaLhs[q2][m2];
                    auto m = M_R_Lhs1Lhs1[Row][q1][m1][q2][m2].block(0,0,N1,N1)*pn;
                    Lhs1Lhs1 += beta_q1*beta_q2*pn.dot(m);
                } // m2 loop on lhs1
            } // q2 loop on lhs1

        } // m1 loop lhs1
    } // q1 loop lhs1


    if ( test )
    {
        cout << "test online Residual : \n RhsRhs="<<RhsRhs << ", Lhs0Rhs="<< Lhs0Rhs
             <<", Lhs0Lhs0="<< Lhs0Lhs0 <<", Lhs0Lhs1="<<Lhs0Lhs1
             <<", Lhs1Lhs1="<< Lhs1Lhs1 <<", Lhs1Rhs="<< Lhs1Rhs <<std::endl;

        auto XhRow = this->M_model->functionSpace()->template functionSpace<Row>();
        auto Xh0 = this->M_model->functionSpace()->template functionSpace<0>();
        auto Xh1 = this->M_model->functionSpace()->template functionSpace<1>();

        auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
        auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();

        vector_ptrtype FRow = backend()->newVector( XhRow );
        vector_ptrtype F0 = backend()->newVector( XhRow );
        vector_ptrtype F1 = backend()->newVector( XhRow );
        vector_ptrtype Vn = backend()->newVector( Xh0 );
        vector_ptrtype Pn = backend()->newVector( Xh1 );
        sparse_matrix_ptrtype A0 = backend()->newMatrix( XhRow, Xh0 );
        sparse_matrix_ptrtype A1 = backend()->newMatrix( XhRow, Xh1 );
        vector_ptrtype Rhs =  backend()->newVector( XhRow );
        vector_ptrtype Lhs0 = backend()->newVector( XhRow );
        vector_ptrtype Lhs1 = backend()->newVector( XhRow );

        auto truc0 = un(0)*XN0->primalBasisElement(0);
        for ( int k=1; k<un.size(); k++ )
            truc0 += un(k)*XN0->primalBasisElement(k);
        *Vn = truc0;

        auto truc1 = pn(0)*XN1->primalBasisElement(0);
        for ( int k=1; k<pn.size(); k++ )
            truc1 += pn(k)*XN1->primalBasisElement(k);
        *Pn = truc1;

        sparse_matrix_ptrtype A;
        std::vector<vector_ptrtype> F;
        boost::tie(boost::tuples::ignore, A, F) = this->M_model->update( mu );

        auto const& i_row = A->mapRow().dofIdToContainerId( Row );
        auto const& i_col0 = A->mapCol().dofIdToContainerId( 0 );
        auto const& i_col1 = A->mapCol().dofIdToContainerId( 1 );

        A0 = A->createSubMatrix( i_row, i_col0 );
        A0->multVector( Vn, F0 );
        F0->scale(-1);
        this->M_model->l2solve( Lhs0, F0, Row );

        A1 = A->createSubMatrix( i_row, i_col1 );
        A1->multVector( Pn, F1 );
        F1->scale(-1);
        this->M_model->l2solve( Lhs1, F1, Row );

        FRow = F[0]->createSubVector( i_row );
        this->M_model->l2solve( Rhs, FRow, Row );

        auto RhsRhs_test = this->M_model->scalarProduct(Rhs,Rhs, Row);
        auto Lhs0Rhs_test = 2.*this->M_model->scalarProduct(Lhs0, Rhs, Row);
        auto Lhs0Lhs0_test = this->M_model->scalarProduct(Lhs0, Lhs0, Row);
        auto Lhs0Lhs1_test = 2.*this->M_model->scalarProduct(Lhs0, Lhs1, Row);
        auto Lhs1Rhs_test = 2.*this->M_model->scalarProduct(Lhs1, Rhs, Row);
        auto Lhs1Lhs1_test = this->M_model->scalarProduct(Lhs1, Lhs1, Row);

        cout << "test block residual : \n RhsRhs="<<RhsRhs_test << ", Lhs0Rhs="<< Lhs0Rhs_test
             <<", Lhs0Lhs0="<< Lhs0Lhs0_test <<", Lhs0Lhs1="<<Lhs0Lhs1_test
             <<", Lhs1Lhs1="<< Lhs1Lhs1_test <<", Lhs1Rhs="<< Lhs1Rhs_test <<std::endl;

        double delta_test = math:: abs( RhsRhs_test + Lhs0Rhs_test + Lhs0Lhs0_test + Lhs0Lhs1_test + Lhs1Rhs_test + Lhs1Lhs1_test );
        cout << "Test Residual delta_test"<<Row <<"="<<delta_test<<std::endl;
    }



    double delta = math:: abs( RhsRhs + Lhs0Rhs + Lhs0Lhs0 + Lhs0Lhs1 + Lhs1Rhs + Lhs1Lhs1 );
    return delta;
}


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::testResidual()
{
    using Feel::cout;
    cout << "\n TEST RESIDUAL \n"
         << "WNmu.size()="<<this->M_WNmu->size()<<std::endl;

    std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
    std::vector< vectorN_type > uNduold;

    for (int k=0; k<this->M_N; k++)
    {
        cout << "====================================================\n";
        parameter_type mu_test = this->M_WNmu->at(k);

        this->lb( this->M_N, mu_test, uN, uNdu, uNold, uNduold );
        cout <<"uN =\n"<<uN[0]<<std::endl;
        cout<< "mu_test = [";
        for ( int i=0; i<mu_test.size(); i++ )
            cout<<mu_test( i )<<",";
        cout<<"]"<<std::endl;
        double rez_test0 = onlineResidualSP<0>( this->M_N, mu_test, uN[0], true );
        double rez_test1 = onlineResidualSP<1>( this->M_N, mu_test, uN[0], true );
        cout << "Test Residual rez0="<< std::sqrt(rez_test0)
             << ", rez1="<< std::sqrt(rez_test1) <<std::endl;
        }
    cout << "====================================================\n";
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::max_error_type
CRBSaddlePoint<TruthModelType>::maxErrorBounds( size_type N ) const
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
        double current_err;
        if( this->M_error_type == CRB_EMPIRICAL )
            current_err = empiricalError( N, current_mu, o.template get<0>() );
        else
            current_err = onlineResidual( N, current_mu, uN[0] );
        if ( current_err > err )
        {
            mu = current_mu;
            err = current_err;
        }
    } // loop on M_WNmu_complement

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
    boost::mpi::broadcast( Environment::worldComm(), tuple, proc_having_good_mu );
    mu = tuple.template get<0>();
    err = tuple.template get<1>();

    Feel::cout << std::setprecision(15) << "[CRBSaddlePoint] max error="<< err << std::endl;

    return boost::make_tuple( err, mu, 0, 0 );
}

template<typename TruthModelType>
double
CRBSaddlePoint<TruthModelType>::empiricalError( int N, parameter_type const& mu, std::vector<double> output_vec ) const
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
CRBSaddlePoint<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    LOG(INFO) <<"[CRBSaddlepoint::save] version : "<<version<<std::endl;

    ar & boost::serialization::base_object<super>( *this );

    ar & BOOST_SERIALIZATION_NVP( M_R_RhsRhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs0 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs1 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Lhs1 );

} //save( ... )


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::load( Archive & ar, const unsigned int version )
{
    LOG(INFO) <<"[CRBSaddlePoint::load] version"<< version <<std::endl;

    ar & boost::serialization::base_object<super>( *this );

    ar & BOOST_SERIALIZATION_NVP( M_R_RhsRhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs0 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs1 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Lhs1 );

} // load( ... )




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
