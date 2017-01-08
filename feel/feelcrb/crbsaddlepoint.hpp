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

// #include <boost/multi_array.hpp>
// #include <boost/tuple/tuple.hpp>
// #include "boost/tuple/tuple_io.hpp"
// #include <boost/format.hpp>
// #include <boost/foreach.hpp>
// #include <boost/bimap.hpp>
// #include <boost/bimap/support/lambda.hpp>
// #include <boost/archive/text_oarchive.hpp>
// #include <boost/archive/text_iarchive.hpp>
// #include <boost/math/special_functions/fpclassify.hpp>
// #include <fstream>

// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/list.hpp>
// #include <boost/serialization/string.hpp>
// #include <boost/serialization/version.hpp>
// #include <boost/serialization/split_member.hpp>

// #include <vector>

// #include <Eigen/Core>
// #include <Eigen/LU>
// #include <Eigen/Dense>

// #include <feel/feelalg/solvereigen.hpp>

// #include <feel/feelcore/environment.hpp>
// #include <feel/feelcore/parameter.hpp>
// #include <feel/feelcrb/parameterspace.hpp>
// #include <feel/feelcrb/crbdb.hpp>
// #include <feel/feelcrb/crbscm.hpp>
// #include <feel/feelcore/serialization.hpp>
// #include <feel/feelfilters/exporter.hpp>

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>

namespace Feel
{

/**
 * \class CRBSaddlePoint
 * \brief Certfified Reduced BAsis for Saddle Point Problems
 *
 * \author JB Wahl
 */
po::options_description crbSaddlePointOptions( std::string const& prefix="" );

template<typename TruthModelType>
class CRBSaddlePoint :
        public CRB<TruthModelType>
{
    typedef CRB<TruthModelType> super;

public:
    //@{ // Truth Model
    typedef TruthModelType model_type;
    typedef boost::shared_ptr<model_type> truth_model_ptrtype;
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

    //@{ /// Database
    typedef CRBElementsDB<model_type> crb_elements_db_type;
    typedef boost::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;
    //@}

    typedef std::vector< std::vector< std::vector< std::vector< matrixN_type >>>> blockmatrixN_type;
    typedef std::vector< std::vector< std::vector< vectorN_type >>> blockvectorN_type;


    CRBSaddlePoint( std::string const& name = "defaultname_crb",
                    WorldComm const& worldComm = Environment::worldComm() ) :
        super( name, worldComm ),
        M_N0(0),
        M_N1(0)
        {}

    CRBSaddlePoint( std::string const& name, truth_model_ptrtype const & model ) :
        super( name, model ),
        M_N0(0),
        M_N1(0)
        {}


    //@{ /// Database
    void saveDB();
    bool loadDB();
    //@}

    element_type expansion( parameter_type const& mu , int N=-1, int time_index=-1);
    element_type expansion( vectorN_type const& u , int const N, bool dual ) const;
    element_type expansionSaddlePoint( vectorN_type const& U_coeff, int const N, bool dual ) const;

private :
    void addBasis( element_type& U, element_type& Udu, parameter_type& mu );
    void orthonormalizeBasis( int number_of_added_elements );
    template <typename WNType>
    double orthonormalize( size_type N, WNType& wn, int Nm, int n_space );
    template <typename WNType>
    double checkOrthonormality( int N, const WNType& wn, int n_space ) const;
    void buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field );
    void saveRB();
    void updateAffineDecompositionSize();
    matrix_info_tuple fixedPointPrimal( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold,
                                        std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput=true ) const;

    void exportBasisFunctions();

    void initBlockMatrix();


    CRBDB M_crbdb;



    int M_N0, M_N1;

    blockmatrixN_type M_blockAqm_pr;
    blockvectorN_type M_blockFqm_pr;
    blockvectorN_type M_blockLqm_pr;



    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()


}; // class CRBSaddlePoint


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::addBasis( element_type& U, element_type& Udu, parameter_type& mu )
{
    auto u = U.template elementPtr<0>();
    auto p = U.template elementPtr<1>();
    auto udu = Udu.template elementPtr<0>();
    auto pdu = Udu.template elementPtr<1>();
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();

     tic();
     XN0->addPrimalBasisElement( u );
     XN0->addDualBasisElement( udu );
     M_N0++;
     toc("Add Basis Function 0");
     tic();
     XN1->addPrimalBasisElement( p );
     XN1->addDualBasisElement( pdu );
     M_N1++;
     toc("Add Basis Function 1");

     if ( boption("crb.saddlepoint.add-supremizer") )
     {
         tic();
         auto us = this->M_model->supremizer( mu, p );
         XN0->addPrimalBasisElement( us );
         XN0->addDualBasisElement( us );
         M_N0++;
         toc("Supremizer computation");
     }

}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::orthonormalizeBasis( int number_of_added_elements )
{
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();


    double norm_max = doption(_name="crb.orthonormality-tol");
    int max_iter = ioption(_name="crb.orthonormality-max-iter");

    if( boption("crb.saddlepoint.orthonormalize0") )
    {
        tic();
        double norm = norm_max+1;
        int iter=0;
        double old = 10;
        int n_added = ( boption("crb.saddlepoint.add-supremizer") ) ? 2:1;
        while( norm >= norm_max && iter < max_iter)
        {
            norm = this->orthonormalize( M_N0, XN0->primalRB(), n_added, 0 );
            iter++;
            //if the norm doesn't change
            if( math::abs(old-norm) < norm_max )
                norm=0;
            old=norm;
        }
        XN0->updatePrimalBasisForUse();
        toc("RB Space Orthnormalization #0");
    }
    if( boption("crb.saddlepoint.orthonormalize1") )
    {
        tic();
        double norm = norm_max+1;
        int iter=0;
        double old = 10;
        while( norm >= norm_max && iter < max_iter )
        {
            norm = this->orthonormalize( M_N1, XN1->primalRB(), 1, 1 );
            iter++;
            if( math::abs(old-norm) < norm_max )
                norm=0;
            old=norm;
        }
        XN1->updatePrimalBasisForUse();
        toc("RB Space Orthnormalization #1");
    }

}

template<typename TruthModelType>
template<typename WNType>
double
CRBSaddlePoint<TruthModelType>::orthonormalize( size_type N, WNType& wn, int Nm, int n_space )
{
    int proc_number = this->worldComm().globalRank();
    Feel::cout << "  -- orthonormalization (Gram-Schmidt)\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for N = " << N << "\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for WN = " << wn.size() << "\n";
    DVLOG(2) << "[CRB::orthonormalize] starting ...\n";

    for ( size_type i =N-Nm; i < N; ++i )
    {
        for ( size_type j = 0; j < i; ++j )
        {
            value_type __rij_pr = this->M_model->scalarProduct(  wn[i], wn[ j ], n_space );
            wn[i].add( -__rij_pr, wn[j] );
        }
    }

    // normalize
    for ( size_type i =N-Nm; i < N; ++i )
    {
        value_type __rii_pr = math::sqrt( this->M_model->scalarProduct(  wn[i], wn[i], n_space ) );
        wn[i].scale( 1./__rii_pr );
    }

    DVLOG(2) << "[CRB::orthonormalize] finished ...\n";
    DVLOG(2) << "[CRB::orthonormalize] copying back results in basis\n";

    return this->checkOrthonormality( N , wn, n_space );
}


template <typename TruthModelType>
template <typename WNType>
double
CRBSaddlePoint<TruthModelType>::checkOrthonormality ( int N, const WNType& wn, int n_space ) const
{
    if ( wn.size()==0 )
    {
        throw std::logic_error( "[CRB::checkOrthonormality] ERROR : size of wn is zero" );
    }


    matrixN_type A, I;
    A.setZero( N, N );
    I.setIdentity( N, N );

    for ( int i = 0; i < N; ++i )
    {
        for ( int j = 0; j < N; ++j )
        {
            A( i, j ) = this->M_model->scalarProduct(  wn[i], wn[j], n_space );
        }
    }

    A -= I;
    DVLOG(2) << "orthonormalization: " << A.norm() << "\n";
    if ( this->worldComm().isMasterRank() )
    {
        LOG( INFO ) << "    o check : " << A.norm() << " (should be 0)";
    }
    //FEELPP_ASSERT( A.norm() < 1e-14 )( A.norm() ).error( "orthonormalization failed.");

    return A.norm();
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field )
{
    tic();
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();

    element_type Ur = this->M_model->functionSpace()->element();
    element_type Uc = this->M_model->functionSpace()->element();
    auto ur = Ur.template element<0>();
    auto uc = Uc.template element<0>();
    auto pr = Ur.template element<1>();
    auto pc = Uc.template element<1>();

    int number_of_elements_to_update = number_of_added_elements;
    // in the case of cobuild, we have to update all since affine decomposition has changed
    if( ioption(_name="ser.rb-frequency") != 0 && !this->M_rebuild)
        number_of_elements_to_update = this->M_N;

    int number_of_elements_to_update0 = boption("crb.saddlepoint.add-supremizer") ? 2*number_of_elements_to_update : number_of_elements_to_update;
    // In case of SER use + error estimation, we compute \hat{A}, \hat{F} (resp. \hat{R}) to compute norm of residual (Riesz)
    int ser_error_estimation = this->M_SER_errorEstimation;

    // update Aqm block matrices
    for ( size_type q=0; q<this->M_model->Qa(); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxA(q); m++ )
        {
            M_blockAqm_pr[0][0][q][m].conservativeResize( M_N0, M_N0 );
            M_blockAqm_pr[0][1][q][m].conservativeResize( M_N0, M_N1 );
            M_blockAqm_pr[1][0][q][m].conservativeResize( M_N1, M_N0 );
            M_blockAqm_pr[1][1][q][m].conservativeResize( M_N1, M_N1 );

            for ( size_type i=M_N0-number_of_elements_to_update0; i<M_N0; i++ )
            {
                //update last row of matrix 00
                for ( size_type j=0; j<M_N0; j++ )
                {
                    ur = XN0->primalBasisElement(i);
                    pr.zero();
                    uc = XN0->primalBasisElement(j);
                    pc.zero();
                    M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

                //update last row of matrix 01
                for ( size_type j=0; j<M_N1; j++ )
                {
                    ur=XN0->primalBasisElement(i);
                    pr.zero();
                    uc.zero();
                    pc=XN1->primalBasisElement(j);
                    M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

            }

            for ( size_type i=M_N1-number_of_elements_to_update; i<M_N1; i++ )
            {
                //update last row of matrix 10
                for ( size_type j=0; j<M_N0; j++ )
                {
                    ur.zero();
                    pr = XN1->primalBasisElement(i);
                    uc = XN0->primalBasisElement(j);
                    pc.zero();
                    M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

                // update last row of matrix 11
                for ( size_type j=0; j<M_N1; j++ )
                {
                    ur.zero();
                    pr = XN1->primalBasisElement(i);
                    uc.zero();
                    pc = XN1->primalBasisElement(j);
                    M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }
            }

            for ( size_type j=M_N0-number_of_elements_to_update0; j<M_N0; j++ )
            {
                //update last column of matrix 00
                for ( size_type i=0; i<M_N0; i++)
                {
                    ur=XN0->primalBasisElement(i);
                    pr.zero();
                    uc=XN0->primalBasisElement(j);
                    pr.zero();
                    M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

                //update last column of matrix 10
                for ( size_type i=0; i<M_N1; i++ )
                {
                    ur.zero();
                    pr=XN1->primalBasisElement(i);
                    uc=XN0->primalBasisElement(j);
                    pc.zero();
                    M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }
            }

            for ( size_type j=M_N1-number_of_elements_to_update; j<M_N1; j++ )
            {
                //update last column of matrix 01
                for ( size_type i=0; i<M_N0; i++ )
                {
                    ur=XN0->primalBasisElement(i);
                    pr.zero();
                    uc.zero();
                    pc=XN1->primalBasisElement(j);
                    M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

                //update last column of matrix 11
                for ( size_type i=0; i<M_N1; i++ )
                {
                    ur.zero();
                    pr=XN1->primalBasisElement(i);
                    uc.zero();
                    pc=XN1->primalBasisElement(j);
                    M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                }

            }
        }
    }


    for ( size_type q=0; q<this->M_model->Ql(0); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxF(0, q); m++ )
        {
            M_blockFqm_pr[0][q][m].conservativeResize(M_N0);
            M_blockFqm_pr[1][q][m].conservativeResize(M_N1);

            // update block 0
            for ( size_type l=M_N0-number_of_elements_to_update0; l<M_N0; l++ )
            {
                ur=XN0->primalBasisElement(l);
                pr.zero();
                M_blockFqm_pr[0][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
            }

            // udpate block 1
            for ( size_type l=M_N1-number_of_elements_to_update; l<M_N1; l++ )
            {
                ur.zero();
                pr=XN1->primalBasisElement(l);
                M_blockFqm_pr[1][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
            }
        }
    }

    for ( size_type q=0; q<this->M_model->Ql(this->M_output_index); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxF(this->M_output_index, q); m++ )
        {
            M_blockLqm_pr[0][q][m].conservativeResize(M_N0);
            M_blockLqm_pr[1][q][m].conservativeResize(M_N1);

            // update block 0
            for ( size_type l=M_N0-number_of_elements_to_update0; l<M_N0; l++ )
            {
                ur=XN0->primalBasisElement(l);
                pr.zero();
                M_blockLqm_pr[0][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
            }

            // udpate block 1
            for ( size_type l=M_N1-number_of_elements_to_update; l<M_N1; l++ )
            {
                ur.zero();
                pr=XN1->primalBasisElement(l);
                M_blockLqm_pr[1][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
            }
        }
    }
    toc("Reduced Matrices Built");
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::initBlockMatrix()
{
    M_blockAqm_pr.resize(2);
    M_blockAqm_pr[0].resize(2);
    M_blockAqm_pr[1].resize(2);

    M_blockFqm_pr.resize(2);
    M_blockLqm_pr.resize(2);
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::updateAffineDecompositionSize()
{
    tic();
    if( this->M_rebuild || this->M_N == 0 )
        initBlockMatrix();
    int output_index = this->M_output_index;

    M_N0 = boption("crb.saddlepoint.add-supremizer") ? 2*this->M_N : this->M_N;
    M_N1 = this->M_N;

    for (int r=0; r<2; r++ )
    {
        for ( int c=0; c<2; c++ )
        {
            M_blockAqm_pr[r][c].resize( this->M_model->Qa() );
            for ( int q=0; q<M_blockAqm_pr[r][c].size(); q++ )
            {
                M_blockAqm_pr[r][c][q].resize( this->M_model->mMaxA(q) );
            }
        }

        M_blockFqm_pr[r].resize( this->M_model->Ql(0) );
        for ( int q=0; q<M_blockFqm_pr[r].size(); q++ )
            M_blockFqm_pr[r][q].resize( this->M_model->mMaxF( 0, q) );

        M_blockLqm_pr[r].resize( this->M_model->Ql( output_index ) );
        for ( int q=0; q<M_blockLqm_pr[r].size(); q++ )
            M_blockLqm_pr[r][q].resize( this->M_model->mMaxF( output_index, q ) );
    }
    toc("Update Affine Decomposition Size (SP)");
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveRB()
{

}





template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::matrix_info_tuple
CRBSaddlePoint<TruthModelType>::fixedPointPrimal(  size_type N, parameter_type const& mu,
                                                   std::vector< vectorN_type > & uN,
                                                   std::vector<vectorN_type> &,
                                                   std::vector< double > & output_vector,
                                                   int K, bool print_rb_matrix,
                                                   bool computeOutput ) const
{
    double output=0;
    int N0 = boption("crb.saddlepoint.add-supremizer") ? 2*N:N;
    int N1 = N;

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

    boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu );

    A.setZero( N0+N1,N0+N1);
    for ( size_type q=0; q<Qa; q++ )
    {
        for ( size_type m=0; m<mMaxA[q]; m++ )
        {
            A.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N0);
            A.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N1);
            A.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N1,N0);
            A.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N1,N1);
        }
    }

    F.setZero(N0+N1);
    for ( size_type q=0; q<Qf; q++ )
    {
        for ( size_type m=0; m<mMaxF[q]; m++ )
        {
            F.segment(0,N0) += betaFqm[0][q][m]*M_blockFqm_pr[0][q][m].head(N0);
            F.segment(N0,N1) += betaFqm[0][q][m]*M_blockFqm_pr[1][q][m].head(N1);
        }
    }

    uN[0] = A.lu().solve( F );

    if ( computeOutput )
    {
        L.setZero(N0+N1);
        for ( size_type q=0; q<Ql; q++ )
        {
            for ( size_type m=0; m<mMaxL[q]; m++ )
            {
                L.segment(0,N0) += betaFqm[this->M_output_index][q][m]*M_blockLqm_pr[0][q][m].head(N0);
                L.segment(N0,N1) += betaFqm[this->M_output_index][q][m]*M_blockLqm_pr[1][q][m].head(N1);
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
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansion( parameter_type const& mu , int N , int time_index )
{
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = this->M_N;

    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold;
    std::vector<vectorN_type> uNduold;

    auto o = this->lb( Nwn, mu, uN, uNdu , uNold, uNduold );
    int size = uN.size();


    return expansionSaddlePoint( uN[size-1], Nwn, false );
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansion( vectorN_type const& u , int const N, bool dual ) const
{
    return expansionSaddlePoint( u, N, dual );
}

template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::expansionSaddlePoint( vectorN_type const& U_coeff, int const N, bool dual ) const
{
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();
    auto WN0 = dual ? XN0->dualRB() : XN0->primalRB();
    auto WN1 = dual ? XN1->dualRB() : XN1->primalRB();

    int Nwn = N>0 ? N:this->M_N;

    int N0 = boption("crb.saddlepoint.add-supremizer") ? 2*Nwn:Nwn;
    int N1 = Nwn;

    CHECK( Nwn <= WN0.size() )<< "invalid expansion size\n";
    CHECK( N0+N1 <= U_coeff.size() )<< "invalid expansion size, Nwn="
                                    << N0+N1
                                    << ", U_coeff.size="<<U_coeff.size()<<std::endl;
    CHECK( U_coeff.size() == N0 + N1 )
        << "invalide size of U_coeff, vector can't be cut\n";

    vectorN_type u_coeff = U_coeff.head( N0 );
    vectorN_type p_coeff = U_coeff.tail( N1 );

    element_type U = this->M_model->functionSpace()->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    u = Feel::expansion( WN0, u_coeff, N0 ).container();
    p = Feel::expansion( WN1, p_coeff, N1 ).container();
    return U;
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::exportBasisFunctions()
{
    tic();
    auto XN0 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<0>();
    auto XN1 = this->M_model->rBFunctionSpace()->template rbFunctionSpace<1>();
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

        std::string basis_name = ( boost::format( "u_pr_%1%.0_param") %index  ).str();
        std::string name = basis_name + mu_str;
        int index0 = boption("crb.saddlepoint.add-supremizer") ? 2*index:index;
        e->step(0)->add( name, u_vec[index0] );

        basis_name = ( boost::format( "p_pr_%1%_param") %index ).str();
        name = basis_name + mu_str;
        e->step(0)->add( name, p_vec[index] );

        if (boption("crb.saddlepoint.add-supremizer"))
        {
            basis_name = ( boost::format( "u_pr_%1%.1_param") %index  ).str();
            name = basis_name + mu_str;
            e->step(0)->add( name, u_vec[index0+1] );
        }
    }
    e->save();
    toc("Export Basis Functions");
}

template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlepoint::save] version : "<<version<<std::endl;


} //save( ... )


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlePoint::load] version"<< version <<std::endl;


} // load( ... )

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveDB()
{

} //saveDB()

template<typename TruthModelType>
bool
CRBSaddlePoint<TruthModelType>::loadDB()
{

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
