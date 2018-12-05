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
 */

#ifndef __CRBSADDLEPOINT_H
#define __CRBSADDLEPOINT_H 1

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
    typedef std::shared_ptr<model_type> truth_model_ptrtype;
    typedef typename model_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //@}

    //@{ /// Parameter Space
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    //@}

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;
    typedef double value_type;
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
        super( name, model, stage, prefixExt ),
        M_N0(0),
        M_N1(0),
        M_addSupremizer(boption("crb.saddlepoint.add-supremizer")),
        M_orthonormalize0(boption("crb.saddlepoint.orthonormalize0")),
        M_orthonormalize1(boption("crb.saddlepoint.orthonormalize1"))
        {
        }

public:
    void init()
    {
        using Feel::cout;
        if ( !this->M_rebuild && this->loadDB() )
        {
            cout << "Database CRB SP " << this->lookForDB() << " available and loaded with M_N0="
                 << M_N0 << ", M_N1="<< M_N1 <<", M_N="<<this->M_N <<std::endl;
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
                    M_N0=0;
                    M_N1=0;
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

    //@{ /// Database
    //! save the CRB SP database
    void saveDB() override;
    //! load the CRB SP Database
    bool loadDB() override;
    //@}


    element_type runWithExpansion( parameter_type const& mu , int N=-1, int time_index=-1) override;
    element_type expansion( vectorN_type const& u, int N = -1,  bool dual=false ) const override;
    element_type expansionSaddlePoint( vectorN_type const& U_coeff, int const N, bool dual ) const;

    element_type solve( parameter_type const& mu )
    {
        return this->M_model->solve( mu );
    }

    void offlineResidual( int Ncur, int number_of_added_elements=1 ) override;

    max_error_type maxErrorBounds( size_type N ) const override;

    matrix_info_tuple fixedPointPrimal( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold, std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput=true ) const override;
    void fixedPointDual(  size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN,
                          std::vector< vectorN_type > & uNdu,  std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K=0) const override;
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
    void addBasis( element_type& U, element_type& Udu, parameter_type& mu ) override;
    void orthonormalizeBasis( int number_of_added_elements ) override;
    template <typename WNType>
    double orthonormalize( size_type N, WNType& wn, int Nm, int n_space );
    template <typename WNType>
    double checkOrthonormality( int N, const WNType& wn, int n_space ) const;
    void buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field ) override;
    void saveRB() override;
    void updateAffineDecompositionSize() override;

    void exportBasisFunctions() override;

    void initBlockMatrix();

    template <int Row>
    void offlineResidualSP( int Ncur , int number_of_added_elements );

    double onlineResidual( int Ncur, parameter_type const& mu, vectorN_type Un ) const;
    template <int Row>
    double onlineResidualSP( int Ncur, parameter_type const& mu, vectorN_type Un, bool test=false ) const;
    void testResidual() override;
    double empiricalError( int N, parameter_type const& mu, std::vector<double> output_vec ) const;

    int M_N0, M_N1;
    bool M_addSupremizer;
    bool M_orthonormalize0;
    bool M_orthonormalize1;

    blockmatrixN_type M_blockAqm_pr;
    blockvectorN_type M_blockFqm_pr;
    blockvectorN_type M_blockLqm_pr;
    blockmatrixN_type M_blockAqm_du;
    blockvectorN_type M_blockFqm_du;
    blockvectorN_type M_blockLqm_du;
    blockmatrixN_type M_blockAqm_pr_du;

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

     if ( M_addSupremizer )
     {
         tic();
         auto us = this->M_model->supremizer( mu, U );
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

    if( M_orthonormalize0 )
    {
        tic();
        double norm = norm_max+1;
        int iter=0;
        double old = 10;
        int n_added = ( M_addSupremizer ) ? 2:1;
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
    if( M_orthonormalize1 )
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
        auto & wni = unwrap_ptr( wn[i] );
        for ( size_type j = 0; j < i; ++j )
        {
            auto & wnj = unwrap_ptr( wn[j] );
            value_type __rij_pr = this->M_model->scalarProduct(  wni, wnj, n_space );
            wni.add( -__rij_pr, wnj );
        }
    }

    // normalize
    for ( size_type i =N-Nm; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        value_type __rii_pr = math::sqrt( this->M_model->scalarProduct(  wni, wni, n_space ) );
        wni.scale( 1./__rii_pr );
    }

    DVLOG(2) << "[CRB::orthonormalize] finished ...\n";
    DVLOG(2) << "[CRB::orthonormalize] copying back results in basis\n";

    return this->checkOrthonormality( N , wn, n_space );
} // orthonormalize()


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
} //checkOrthonormality()

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field )
{
    tic();
    int nbBlock = 2;
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

    int number_of_elements_to_update0 = M_addSupremizer ? 2*number_of_elements_to_update : number_of_elements_to_update;
    // In case of SER use + error estimation, we compute \hat{A}, \hat{F} (resp. \hat{R}) to compute norm of residual (Riesz)
    int ser_error_estimation = this->M_SER_errorEstimation;
    if ( ioption("crb.saddlepoint.version")==2 )
        this->M_model->initBlockMatrix();

    // update Aqm block matrices
    for ( size_type q=0; q<this->M_model->Qa(); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxA(q); m++ )
        {
            for( int k = 0; k < nbBlock*nbBlock; ++k )
            {
                int n0 = k%nbBlock == 0 ? M_N0 : M_N1;
                int n1 = k/nbBlock == 0 ? M_N0 : M_N1;
                M_blockAqm_pr[k%nbBlock][k/nbBlock][q][m].conservativeResize( n0, n1 );
                M_blockAqm_du[k%nbBlock][k/nbBlock][q][m].conservativeResize( n0, n1 );
                M_blockAqm_pr_du[k%nbBlock][k/nbBlock][q][m].conservativeResize( n0, n1 );
            }

            for ( size_type i=M_N0-number_of_elements_to_update0; i<M_N0; i++ )
            {
                //update last row of matrix 00
                for ( size_type j=0; j<M_N0; j++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur = XN0->primalBasisElement(i);
                        pr.zero();
                        uc = XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur = XN0->dualBasisElement(i);
                        pr.zero();
                        uc = XN0->dualBasisElement(j);
                        pc.zero();
                        M_blockAqm_du[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur = XN0->dualBasisElement(i);
                        pr.zero();
                        uc = XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr_du[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN0->primalBasisElement(j), 0,0 );
                        M_blockAqm_du[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->dualBasisElement(i), XN0->dualBasisElement(j), 0,0, true );
                        M_blockAqm_pr_du[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN0->dualBasisElement(j), 0,0 );
                    }
                }

                //update last row of matrix 01
                for ( size_type j=0; j<M_N1; j++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur=XN0->primalBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->dualBasisElement(j);
                        M_blockAqm_du[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr_du[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN1->primalBasisElement(j), 0,1 );
                        M_blockAqm_du[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->dualBasisElement(i), XN1->dualBasisElement(j), 0,1, true );
                        M_blockAqm_pr_du[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN1->dualBasisElement(j), 0,1 );
                    }
                }
            }

            for ( size_type i=M_N1-number_of_elements_to_update; i<M_N1; i++ )
            {
                //update last row of matrix 10
                for ( size_type j=0; j<M_N0; j++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur.zero();
                        pr = XN1->primalBasisElement(i);
                        uc = XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur.zero();
                        pr = XN1->dualBasisElement(i);
                        uc = XN0->dualBasisElement(j);
                        pc.zero();
                        M_blockAqm_du[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur.zero();
                        pr = XN1->dualBasisElement(i);
                        uc = XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr_du[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN0->primalBasisElement(j), 1,0 );
                        M_blockAqm_du[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->dualBasisElement(i), XN0->dualBasisElement(j), 1,0, true );
                        M_blockAqm_pr_du[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN0->dualBasisElement(j), 1,0 );
                    }
                }

                // update last row of matrix 11
                for ( size_type j=0; j<M_N1; j++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur.zero();
                        pr = XN1->primalBasisElement(i);
                        uc.zero();
                        pc = XN1->primalBasisElement(j);
                        M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur.zero();
                        pr = XN1->dualBasisElement(i);
                        uc.zero();
                        pc = XN1->dualBasisElement(j);
                        M_blockAqm_du[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur.zero();
                        pr = XN1->dualBasisElement(i);
                        uc.zero();
                        pc = XN1->primalBasisElement(j);
                        M_blockAqm_pr_du[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN1->primalBasisElement(j), 1,1);
                        M_blockAqm_du[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->dualBasisElement(i), XN1->dualBasisElement(j), 1,1, true);
                        M_blockAqm_pr_du[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN1->dualBasisElement(j), 1,1);
                    }
                }
            }

            for ( size_type j=M_N0-number_of_elements_to_update0; j<M_N0; j++ )
            {
                //update last column of matrix 00
                for ( size_type i=0; i<M_N0-number_of_elements_to_update0; i++)
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur=XN0->primalBasisElement(i);
                        pr.zero();
                        uc=XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc=XN0->dualBasisElement(j);
                        pc.zero();
                        M_blockAqm_du[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc=XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr_du[0][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN0->primalBasisElement(j), 0,0 );
                        M_blockAqm_du[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->dualBasisElement(i), XN0->dualBasisElement(j), 0,0, true );
                        M_blockAqm_pr_du[0][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN0->dualBasisElement(j), 0,0 );
                    }
                }

                //update last column of matrix 10
                for ( size_type i=0; i<M_N1-number_of_elements_to_update; i++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur.zero();
                        pr=XN1->primalBasisElement(i);
                        uc=XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur.zero();
                        pr=XN1->dualBasisElement(i);
                        uc=XN0->dualBasisElement(j);
                        pc.zero();
                        M_blockAqm_du[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur.zero();
                        pr=XN1->dualBasisElement(i);
                        uc=XN0->primalBasisElement(j);
                        pc.zero();
                        M_blockAqm_pr_du[1][0][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN0->primalBasisElement(j), 1,0 );
                        M_blockAqm_du[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->dualBasisElement(i), XN0->dualBasisElement(j), 1,0, true);
                        M_blockAqm_pr_du[1][0][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN0->dualBasisElement(j), 1,0 );
                    }
                }
            }

            for ( size_type j=M_N1-number_of_elements_to_update; j<M_N1; j++ )
            {
                //update last column of matrix 01
                for ( size_type i=0; i<M_N0-number_of_elements_to_update0; i++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur=XN0->primalBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->dualBasisElement(j);
                        M_blockAqm_du[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur=XN0->dualBasisElement(i);
                        pr.zero();
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr_du[0][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN1->primalBasisElement(j), 0,1 );
                        M_blockAqm_du[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->dualBasisElement(i), XN1->dualBasisElement(j), 0,1, true );
                        M_blockAqm_pr_du[0][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN0->primalBasisElement(i), XN1->dualBasisElement(j), 0,1 );
                    }
                }

                //update last column of matrix 11
                for ( size_type i=0; i<M_N1-number_of_elements_to_update; i++ )
                {
                    if ( ioption("crb.saddlepoint.version")==1 )
                    {
                        ur.zero();
                        pr=XN1->primalBasisElement(i);
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                        ur.zero();
                        pr=XN1->dualBasisElement(i);
                        uc.zero();
                        pc=XN1->dualBasisElement(j);
                        M_blockAqm_du[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur, true );
                        ur.zero();
                        pr=XN1->dualBasisElement(i);
                        uc.zero();
                        pc=XN1->primalBasisElement(j);
                        M_blockAqm_pr_du[1][1][q][m](i,j) = this->M_model->Aqm( q, m, Uc, Ur );
                    }
                    else if ( ioption("crb.saddlepoint.version")==2 )
                    {
                        M_blockAqm_pr[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN1->primalBasisElement(j), 1,1 );
                        M_blockAqm_du[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->dualBasisElement(i), XN1->dualBasisElement(j), 1,1, true );
                        M_blockAqm_pr_du[1][1][q][m](i,j) = this->M_model->AqmBlock( q, m, XN1->primalBasisElement(i), XN1->dualBasisElement(j), 1,1 );
                    }
                }

            }
        }
    }

    // update Fqm block vectors
    for ( size_type q=0; q<this->M_model->Ql(0); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxF(0, q); m++ )
        {
            M_blockFqm_pr[0][q][m].conservativeResize(M_N0);
            M_blockFqm_pr[1][q][m].conservativeResize(M_N1);
            M_blockFqm_du[0][q][m].conservativeResize(M_N0);
            M_blockFqm_du[1][q][m].conservativeResize(M_N1);

            // update block 0
            for ( size_type l=M_N0-number_of_elements_to_update0; l<M_N0; l++ )
            {
                if ( ioption("crb.saddlepoint.version")==1 )
                {
                    ur=XN0->primalBasisElement(l);
                    pr.zero();
                    M_blockFqm_pr[0][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
                    ur=XN0->dualBasisElement(l);
                    pr.zero();
                    M_blockFqm_du[0][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
                }
                else if ( ioption("crb.saddlepoint.version")==2 )
                {
                    M_blockFqm_pr[0][q][m](l) = this->M_model->FqmBlock( 0, q, m, XN0->primalBasisElement(l), 0 );
                    M_blockFqm_du[0][q][m](l) = this->M_model->FqmBlock( 0, q, m, XN0->dualBasisElement(l), 0 );
                }
            }

            // udpate block 1
            for ( size_type l=M_N1-number_of_elements_to_update; l<M_N1; l++ )
            {
                if ( ioption("crb.saddlepoint.version")==1 )
                {
                    ur.zero();
                    pr=XN1->primalBasisElement(l);
                    M_blockFqm_pr[1][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
                    ur.zero();
                    pr=XN1->dualBasisElement(l);
                    M_blockFqm_du[1][q][m](l) = this->M_model->Fqm( 0, q, m, Ur );
                }
                else if ( ioption("crb.saddlepoint.version")==2 )
                {
                    M_blockFqm_pr[1][q][m](l) = this->M_model->FqmBlock( 0, q, m, XN1->primalBasisElement(l), 1 );
                    M_blockFqm_du[1][q][m](l) = this->M_model->FqmBlock( 0, q, m, XN1->dualBasisElement(l), 1 );
                }
            }
        }
    }

    // update Lqm block vectors
    for ( size_type q=0; q<this->M_model->Ql(this->M_output_index); q++ )
    {
        for ( size_type m=0; m<this->M_model->mMaxF(this->M_output_index, q); m++ )
        {
            M_blockLqm_pr[0][q][m].conservativeResize(M_N0);
            M_blockLqm_pr[1][q][m].conservativeResize(M_N1);
            M_blockLqm_du[0][q][m].conservativeResize(M_N0);
            M_blockLqm_du[1][q][m].conservativeResize(M_N1);

            // update block 0
            for ( size_type l=M_N0-number_of_elements_to_update0; l<M_N0; l++ )
            {
                if ( ioption("crb.saddlepoint.version")==1 )
                {
                    ur=XN0->primalBasisElement(l);
                    pr.zero();
                    M_blockLqm_pr[0][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
                    ur=XN0->dualBasisElement(l);
                    pr.zero();
                    M_blockLqm_du[0][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
                }
                else if ( ioption("crb.saddlepoint.version")==2 )
                {
                    M_blockLqm_pr[0][q][m](l) = this->M_model->FqmBlock( this->M_output_index, q, m, XN0->primalBasisElement(l), 0 );
                    M_blockLqm_du[0][q][m](l) = this->M_model->FqmBlock( this->M_output_index, q, m, XN0->dualBasisElement(l), 0 );
                }
            }

            // udpate block 1
            for ( size_type l=M_N1-number_of_elements_to_update; l<M_N1; l++ )
            {
                if ( ioption("crb.saddlepoint.version")==1 )
                {
                    ur.zero();
                    pr=XN1->primalBasisElement(l);
                    M_blockLqm_pr[1][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
                    ur.zero();
                    pr=XN1->dualBasisElement(l);
                    M_blockLqm_du[1][q][m](l) = this->M_model->Fqm( this->M_output_index, q, m, Ur );
                }
                else if ( ioption("crb.saddlepoint.version")==2 )
                {
                    M_blockLqm_pr[1][q][m](l) = this->M_model->FqmBlock( this->M_output_index, q, m, XN1->primalBasisElement(l), 1 );
                    M_blockLqm_du[1][q][m](l) = this->M_model->FqmBlock( this->M_output_index, q, m, XN1->dualBasisElement(l), 1 );
                }
            }
        }
    }



    if ( ioption("crb.saddlepoint.version")==2 )
        this->M_model->clearBlockMatix();

    toc("Reduced Matrices Built");
} // buildRbmatrix()

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::initBlockMatrix()
{
    int nbBlock = 2;
    M_blockAqm_pr.resize(nbBlock);
    M_blockAqm_du.resize(nbBlock);
    M_blockAqm_pr_du.resize(nbBlock);
    for( int i = 0; i < nbBlock; i++ )
    {
        M_blockAqm_pr[i].resize(nbBlock);
        M_blockAqm_du[i].resize(nbBlock);
        M_blockAqm_pr_du[i].resize(nbBlock);
    }

    M_blockFqm_pr.resize(nbBlock);
    M_blockLqm_pr.resize(nbBlock);
    M_blockFqm_du.resize(nbBlock);
    M_blockLqm_du.resize(nbBlock);
}

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::updateAffineDecompositionSize()
{
    if( this->M_rebuild || this->M_N == 0 )
        initBlockMatrix();
    int output_index = this->M_output_index;

    M_N0 = M_addSupremizer ? 2*this->M_N : this->M_N;
    M_N1 = this->M_N;

    for (int r=0; r<2; r++ )
    {
        for ( int c=0; c<2; c++ )
        {
            M_blockAqm_pr[r][c].resize( this->M_model->Qa() );
            M_blockAqm_du[r][c].resize( this->M_model->Qa() );
            M_blockAqm_pr_du[r][c].resize( this->M_model->Qa() );
            for ( int q=0; q<M_blockAqm_pr[r][c].size(); q++ )
            {
                M_blockAqm_pr[r][c][q].resize( this->M_model->mMaxA(q) );
                M_blockAqm_du[r][c][q].resize( this->M_model->mMaxA(q) );
                M_blockAqm_pr_du[r][c][q].resize( this->M_model->mMaxA(q) );
            }
        }

        M_blockFqm_pr[r].resize( this->M_model->Ql(0) );
        M_blockFqm_du[r].resize( this->M_model->Ql(0) );
        for ( int q=0; q<M_blockFqm_pr[r].size(); q++ )
        {
            M_blockFqm_pr[r][q].resize( this->M_model->mMaxF( 0, q) );
            M_blockFqm_du[r][q].resize( this->M_model->mMaxF( 0, q) );
        }

        M_blockLqm_pr[r].resize( this->M_model->Ql( output_index ) );
        M_blockLqm_du[r].resize( this->M_model->Ql( output_index ) );
        for ( int q=0; q<M_blockLqm_pr[r].size(); q++ )
        {
            M_blockLqm_pr[r][q].resize( this->M_model->mMaxF( output_index, q ) );
            M_blockLqm_du[r][q].resize( this->M_model->mMaxF( output_index, q ) );
        }
    }

    if ( this->M_error_type == CRB_RESIDUAL || this->M_error_type == CRB_RESIDUAL_SCM )
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
    } // if Residual
} // updateAffinedecompositionsize()

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveRB()
{
    this->M_elements_database.saveDB();
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
    bool is_linear = this->M_model->isLinear();
    double output=0;
    double increment = this->M_fixedpointIncrementTol;
    int N0 = M_addSupremizer ? 2*N:N;
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

    if ( is_linear )
    {
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu );

        A.setZero( N0+N1,N0+N1);
        for ( size_type q=0; q<Qa; q++ )
        {
            for ( size_type m=0; m<mMaxA[q]; m++ )
            {
                A.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N0);
                A.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_pr[0][1][q][m].block(0,0,N0,N1);
                A.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_pr[1][0][q][m].block(0,0,N1,N0);
                A.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_pr[1][1][q][m].block(0,0,N1,N1);
            }
        }

        F.setZero(N0+N1);
        for ( size_type q=0; q<Qf; q++ )
        {
            for ( size_type m=0; m<mMaxF[q]; m++ )
            {
                F.head(N0) += betaFqm[0][q][m]*M_blockFqm_pr[0][q][m].head(N0);
                F.tail(N1) += betaFqm[0][q][m]*M_blockFqm_pr[1][q][m].head(N1);
            }
        }
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
            // if ( this->M_useRbSpaceContextEim && this->M_hasRbSpaceContextEim )
            //     boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
            //         this->M_model->computeBetaQm( uN[0], mu/*, N*/ );
            // else
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                    this->M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
            A.setZero( N0+N1,N0+N1);
            for ( size_type q=0; q<Qa; q++ )
            {
                for ( size_type m=0; m<mMaxA[q]; m++ )
                {
                    A.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_pr[0][0][q][m].block(0,0,N0,N0);
                    A.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_pr[0][1][q][m].block(0,0,N0,N1);
                    A.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_pr[1][0][q][m].block(0,0,N1,N0);
                    A.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_pr[1][1][q][m].block(0,0,N1,N1);
                }
            }

            F.setZero(N0+N1);
            for ( size_type q=0; q<Qf; q++ )
            {
                for ( size_type m=0; m<mMaxF[q]; m++ )
                {
                    F.head(N0) += betaFqm[0][q][m]*M_blockFqm_pr[0][q][m].head(N0);
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
                L.head(N0) += betaFqm[this->M_output_index][q][m]*M_blockLqm_pr[0][q][m].head(N0);
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
    int N0 = M_addSupremizer ? 2*N:N;
    int N1 = N;

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
                Adu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_du[0][0][q][m].block(0,0,N0,N0);
                Adu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_du[0][1][q][m].block(0,0,N0,N1);
                Adu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_du[1][0][q][m].block(0,0,N1,N0);
                Adu.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_du[1][1][q][m].block(0,0,N1,N1);
            }
        }
        Fdu.setZero( N0+N1 );
        for ( size_type q = 0; q < Ql ; ++q )
            {
                for(int m=0; m < mMaxF[q]; m++)
                {
                    Fdu.head(N0) -= betaFqm[this->M_output_index][q][m]*M_blockLqm_du[0][q][m].head( N0 );
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
                    Adu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_du[0][0][q][m].block(0,0,N0,N0);
                    Adu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_du[0][1][q][m].block(0,0,N0,N1);
                    Adu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_du[1][0][q][m].block(0,0,N1,N0);
                    Adu.block( N0, N0, N1, N1 ) += betaAqm[q][m]*M_blockAqm_du[1][1][q][m].block(0,0,N1,N1);
                }
            }

            Fdu.setZero( N0+N1 );
            for ( size_type q=0; q<Ql; q++ )
            {
                for ( size_type m=0; m<mMaxF[q]; m++ )
                {
                    Fdu.head(N0) += betaFqm[0][q][m]*M_blockFqm_du[0][q][m].head(N0);
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
    int N = uN[0].size();
    int Ni = M_addSupremizer ? N/3 : N/2;
    int N0 = M_addSupremizer ? 2*Ni : Ni;
    int N1 = Ni;

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
        // if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
        //     boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = this->M_model->computeBetaQm( uN[0], mu/*, N*/ );
        // else
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                this->M_model->computeBetaQm( this->expansion( uN[0], Ni ), mu );
    }
    for(size_type q = 0;q < this->M_model->Ql(0); ++q)
    {
        for(int m=0; m < this->M_model->mMaxF(0,q); m++)
        {
            Fdu.head(N0) += betaFqm[0][q][m]*M_blockFqm_du[0][q][m].head(N0);
            Fdu.tail(N1) += betaFqm[0][q][m]*M_blockFqm_du[1][q][m].head(N1);
        }
    }
    for(size_type q = 0;q < this->M_model->Qa(); ++q)
    {
        for(int m=0; m < this->M_model->mMaxA(q); m++)
        {
            Aprdu.block( 0, 0, N0, N0 ) += betaAqm[q][m]*M_blockAqm_pr_du[0][0][q][m].block(0,0,N0,N0);
            Aprdu.block( 0, N0, N0, N1 ) += betaAqm[q][m]*M_blockAqm_pr_du[0][1][q][m].block(0,0,N0,N1);
            Aprdu.block( N0, 0, N1, N0 ) += betaAqm[q][m]*M_blockAqm_pr_du[1][0][q][m].block(0,0,N1,N0);
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

    bool optimize = boption(_name="crb.optimize-offline-residual") ;
    int N0 = M_addSupremizer ? 2*Ncur:Ncur;
    int N1 = Ncur;
    int n_added0 = M_addSupremizer ? 2*number_of_added_elements:number_of_added_elements;
    int n_added1 = number_of_added_elements;

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
                this->M_model->l2solveSP( Z1, Fqm[q][m], Row );
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
                this->M_model->l2solveSP( Z1, W, Row );
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
                this->M_model->l2solveSP( Z1, W, Row );
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

    this->M_model->clearBlockMatix();

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
    int N0 = M_addSupremizer ? 2*Ncur:Ncur;
    int N1 = Ncur;

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
        this->M_model->l2solveSP( Lhs0, F0, Row );

        A1 = A->createSubMatrix( i_row, i_col1 );
        A1->multVector( Pn, F1 );
        F1->scale(-1);
        this->M_model->l2solveSP( Lhs1, F1, Row );

        FRow = F[0]->createSubVector( i_row );
        this->M_model->l2solveSP( Rhs, FRow, Row );

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
typename CRBSaddlePoint<TruthModelType>::element_type
CRBSaddlePoint<TruthModelType>::runWithExpansion( parameter_type const& mu , int N , int time_index )
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
CRBSaddlePoint<TruthModelType>::expansion( vectorN_type const& u, int N, bool dual ) const
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

    int N0 = M_addSupremizer ? 2*Nwn:Nwn;
    int N1 = Nwn;

    CHECK( Nwn <= WN0.size() )<< "invalid expansion size\n";
    CHECK( N0+N1 <= U_coeff.size() )<< "invalid expansion size, Nwn="
                                    << N0+N1
                                    << ", U_coeff.size="<<U_coeff.size()<<std::endl;
    CHECK( U_coeff.size() == N0 + N1 )
        << "invalide size of U_coeff, vector can't be cut\n";

    vectorN_type u_coeff = U_coeff.head( N0 );
    vectorN_type p_coeff = U_coeff.tail( N1 );

    //Feel::cout << U_coeff << std::endl;

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
        int index0 = M_addSupremizer ? 2*index:index;
        e->step(0)->add( name, unwrap_ptr( u_vec[index0] ) );

        basis_name = ( boost::format( "p_pr_%1%_param") %index ).str();
        name = basis_name + mu_str;
        e->step(0)->add( name, unwrap_ptr( p_vec[index] ) );

        if (M_addSupremizer)
        {
            basis_name = ( boost::format( "u_pr_%1%.1_param") %index  ).str();
            name = basis_name + mu_str;
            e->step(0)->add( name, unwrap_ptr( u_vec[index0+1] ) );
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

    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_N0 );
    ar & BOOST_SERIALIZATION_NVP( M_N1 );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr_du );

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

    int proc_number = this->worldComm().globalRank();
    LOG(INFO) <<"[CRBSaddlePoint::load] version"<< version <<std::endl;

    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_N0 );
    ar & BOOST_SERIALIZATION_NVP( M_N1 );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr_du );

    ar & BOOST_SERIALIZATION_NVP( M_R_RhsRhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Rhs );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs0 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs0Lhs1 );
    ar & BOOST_SERIALIZATION_NVP( M_R_Lhs1Lhs1 );

} // load( ... )

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveDB()
{
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );

    if ( ofs )
    {
        // boost::archive::text_oarchive oa( ofs );
        boost::archive::binary_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
    this->saveJson();
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
        // boost::archive::text_iarchive ia( ifs );
        boost::archive::binary_iarchive ia( ifs );
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
