/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-09

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file crbmodel.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-09
 */
#ifndef __CRBModel_H
#define __CRBModel_H 1

#include <boost/shared_ptr.hpp>

#include <vector>


#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>


namespace Feel
{
enum class CRBModelMode
{
    PFEM = 0, SCM = 1, CRB = 2, SCM_ONLINE=3, CRB_ONLINE=4
};

/**
 * \class CRBModel
 * \brief Certified Reduced Basis Model class
 *
 * This class implements the requirements over a model to be usable by the
 * certified reduced basis method.
 *
 * \tparam ModelType the type of the finite element method model
 *
 * The FEM model type should derive from this class and fill the vector of
 * matrices M_Aq and vector of vectors M_Fq
 *
 * @author Christophe Prud'homme
 * @see crb
 */
template<typename ModelType>
class CRBModel
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type ParameterSpaceDimension = ModelType::ParameterSpaceDimension;
    static const bool is_time_dependent = ModelType::is_time_dependent;

    //@}

    /** @name Typedefs
     */
    //@{

    //! model type
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    //! value_type
    typedef typename model_type::value_type value_type;
    //! mesh type
    typedef typename ModelType::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::vector_type vector_type;

    typedef typename model_type::eigen_matrix_type eigen_matrix_type;

    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename model_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename model_type::parameter_type parameter_type;
    typedef typename model_type::parameter_ptrtype parameter_ptrtype;
    typedef typename model_type::sampling_type sampling_type;
    typedef typename model_type::sampling_ptrtype sampling_ptrtype;


    typedef typename std::vector< std::vector < element_ptrtype > > initial_guess_type;
    //typedef Eigen::VectorXd theta_vector_type;
    //typedef Eigen::VectorXd vectorN_type;
    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype,
                                  sparse_matrix_ptrtype,
                                  std::vector<vector_ptrtype>,
                                  vector_ptrtype
                                  > offline_merge_type;


    typedef typename boost::tuple<std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector< std::vector<vector_ptrtype> > >,
                                  std::vector< std::vector< vector_ptrtype> >
                                  > affine_decomposition_type;


    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>,
                                   beta_vector_type
                                   > betaqm_type;

    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                       typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                  typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                     fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                     >::type >::type >::type index_vector_type;



    //@}

    /** @name Constructors, destructor
     */
    //@{

    CRBModel()
        :
        M_Aqm(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_mode( CRBModelMode::PFEM ),
        M_model( new model_type() ),
        M_backend( backend_type::build( BACKEND_PETSC ) ),
        M_B()
    {
        this->init();
    }

    CRBModel( po::variables_map const& vm, CRBModelMode mode = CRBModelMode::PFEM  )
        :
        M_Aqm(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm( vm ),
        M_mode( mode ),
        M_model( new model_type( vm ) ),
        M_backend( backend_type::build( vm ) ),
        M_B()
    {
        this->init();
    }

    /**
     * \param model the model to be used
     */
    CRBModel( model_ptrtype & model )
        :
        M_Aqm(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm(),
        M_mode( CRBModelMode::PFEM ),
        M_model( model ),
        M_backend( backend_type::build( model->vm ) ),
        M_B()
    {
        this->init();
    }

    /**
     * copy constructor
     */
    CRBModel( CRBModel const & o )
        :
        M_Aqm( o.M_Aqm ),
        M_Mqm( o.M_Mqm ),
        M_Fqm( o.M_Fqm ),
        M_is_initialized( o.M_is_initialized ),
        M_vm( o.M_vm ),
        M_mode( o.M_mode ),
        M_model(  o.M_model ),
        M_backend( o.M_backend ),
        M_B( o.M_B )
    {
        this->init();
    }

    //! destructor
    virtual ~CRBModel()
    {}

    //! initialize the model (mesh, function space, operators, matrices, ...)
    FEELPP_DONT_INLINE void init()
    {
        if ( M_is_initialized )
            return;

        M_is_initialized=true;

        if ( M_mode != CRBModelMode::CRB_ONLINE &&
                M_mode != CRBModelMode::SCM_ONLINE )
        {
            //the model is already initialized
            //std::cout << "  -- init FEM  model\n";
            //M_model->init();
            this->initB();
        }
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRBModel& operator=( CRBModel const & o )
    {
        if ( this != &o )
        {
            M_Aqm = o.M_Aqm;
            M_Fqm = o.M_Fqm;
            M_Mqm = o.M_Mqm;
            M_model = o.M_model;
            M_backend = o.M_backend;
            M_B = o.M_B;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return  the \p variables_map
     */
    po::variables_map vm() const
    {
        return M_vm;
    }

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    virtual sparse_matrix_ptrtype newMatrix() const
    {
        return M_model->newMatrix();
    }

    /**
     * create a new vector
     * \return the newly created vector
     */
    virtual vector_ptrtype newVector() const
    {
        return M_model->newVector();
    }

    /**
     * \brief Returns the matrix associated with the \f$H_1\f$ inner product
     */
    sparse_matrix_ptrtype const& innerProduct() const
    {
        return M_B;
    }

    /**
     * \brief Returns the matrix associated with the \f$H_1\f$ inner product
     */
    sparse_matrix_ptrtype const& h1() const
    {
        return M_B;
    }

    //!  Returns the function space
    functionspace_ptrtype  functionSpace() const
    {
        return M_model->functionSpace();
    }

    //! return the number of \f$\mu\f$ independent terms for the bilinear form
    size_type Qa() const
    {
        return M_model->Qa();
    }

    //! return the number of \f$\mu\f$ independent terms for the bilinear form ( time dependent )
    //int Qm() const { return 1; }//return M_model->Qm(); }

    size_type Qm() const
    {
        return Qm( mpl::bool_<model_type::is_time_dependent>() );
    }
    size_type Qm( mpl::bool_<true> ) const
    {
        return M_model->Qm();
    }
    size_type Qm( mpl::bool_<false> ) const
    {
        return 1;
    }

    size_type Qmf() const
    {
        return M_model->QInitialGuess();
    }

    int mMaxA(int q )
    {
        return M_model->mMaxA( q );
    }


    int mMaxM( int q )
    {
        return mMaxM(q, mpl::bool_<model_type::is_time_dependent>() );
    }
    int mMaxM( int q, mpl::bool_<true> )
    {
        return M_model->mMaxM( q );
    }
    int mMaxM( int q , mpl::bool_<false> )
    {
        return 1;
    }

    int mMaxMF( int q ) const
    {
        return M_model->mMaxInitialGuess( q );
    }

    int mMaxF(int output_index, int q )
    {
        return M_model->mMaxF( output_index, q );
    }

    //! return the number of outputs
    size_type Nl() const
    {
        return M_model->Nl();
    }

    //! return the number of \f$\mu\f$ independent terms for the right hand side
    size_type Ql( int l ) const
    {
        return M_model->Ql( l );
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_model->parameterSpace();
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the mesh characteristic length to \p s
     */
    void setMeshSize( double s )
    {
        M_model->setMeshSize( s );
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \brief compute the betaqm given \p mu
     */
    betaqm_type computeBetaQm( parameter_type const& mu , double time=0 )
    {
        return computeBetaQm( mu , mpl::bool_<model_type::is_time_dependent>(), time  );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<true>, double time=0 )
    {
        return M_model->computeBetaQm( mu , time );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false>, double time=0 )
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm, betaInitialGuessQm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type>,
            beta_vector_type >
        steady_beta;

        steady_beta = M_model->computeBetaQm( mu , time );
        betaAqm = steady_beta.get<0>();
        betaFqm = steady_beta.get<1>();
        betaInitialGuessQm = steady_beta.get<2>();

        betaMqm.resize( 1 );
        betaMqm[0].resize( 1 );
        betaMqm[0][0] = 1 ;

        return boost::make_tuple( betaMqm, betaAqm, betaFqm, betaInitialGuessQm );
    }


    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , double time=0 )
    {
        return computeBetaQm( T , mu , mpl::bool_<model_type::is_time_dependent>(), time );
    }
    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<true>, double time=0 )
    {
        return M_model->computeBetaQm( T, mu, time );
    }
    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<false>, double time=0 )
    {
        beta_vector_type betaAqm, betaMqm, betaInitialGuessQm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<  beta_vector_type,
                       std::vector<beta_vector_type>,
                       beta_vector_type >
            steady_beta;

        steady_beta = M_model->computeBetaQm(T, mu , time );
        betaAqm = steady_beta.get<0>();
        betaFqm = steady_beta.get<1>();
        betaInitialGuessQm = steady_beta.get<2>();
        betaMqm.resize( 1 );
        betaMqm[0].resize( 1 );
        betaMqm[0][0] = 1 ;

        return boost::make_tuple( betaMqm, betaAqm, betaFqm, betaInitialGuessQm );
    }


    element_ptrtype initialGuess( parameter_type const& mu );
    element_ptrtype initialGuess( parameter_type const& mu , mpl::bool_<true>);
    element_ptrtype initialGuess( parameter_type const& mu , mpl::bool_<false>);

    /**
     * \brief update the model wrt \p mu
     */
    offline_merge_type update( parameter_type const& mu,  double time=0 )
    {
        M_model->computeBetaQm( mu , time );
        return offlineMerge( mu );
    }
    offline_merge_type update( parameter_type const& mu, element_type const& T, double time=0 )
    {
        M_model->computeBetaQm(T, mu , time );
        return offlineMerge( mu );
    }

    /**
     * \brief Compute the affine decomposition of the various forms
     *
     * This function defined in the \p M_model assembles the parameter
     * independant part of the affine decomposition of the bilinear and linear
     * forms.
     */
    affine_decomposition_type computeAffineDecomposition()
    {
        return computeAffineDecomposition( mpl::bool_<model_type::is_time_dependent>() );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<true> )
    {
        initial_guess_type initial_guess;
        boost::tie( M_Mqm, M_Aqm, M_Fqm, initial_guess ) = M_model->computeAffineDecomposition();
        assembleMF( initial_guess );
        return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm, M_MFqm );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        initial_guess_type initial_guess;
        boost::tie( M_Aqm, M_Fqm, initial_guess ) = M_model->computeAffineDecomposition();
        assembleMassMatrix();
        assembleMF( initial_guess );
        return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm, M_MFqm );
    }


    /**
     * \brief the inner product \f$h1(\xi_i, \xi_j) = \xi_j^T H_1 \xi_i\f$
     *
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c A_q
     *
     * \return the inner product \f$h1(\xi_i, \xi_j) = \xi_j^T H_1 \xi_i\f$
     */
    value_type h1( element_type const& xi_i, element_type const& xi_j  )
    {
        return M_B->energy( xi_j, xi_i );
    }
    /**
     * \brief the inner product \f$h1(\xi_i, \xi_j) = \xi_j^T H_1 \xi_i\f$
     *
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c A_q
     *
     * \return the inner product \f$h1(\xi_i, \xi_j) = \xi_j^T H_1 \xi_i\f$
     */
    value_type h1( element_type const& xi_i  )
    {
        return M_B->energy( xi_i, xi_i );
    }



    /**
     * \brief Returns the matrix \c Aq[q][m] of the affine decomposition of the bilinear form
     *
     * \param q and m are index of the component in the affine decomposition
     * \param transpose transpose \c A_q
     *
     * \return the matrix \c Aq[q][m] of the affine decomposition of the bilinear form
     */
    sparse_matrix_ptrtype Aqm( uint16_type q, uint16_type m, bool transpose = false ) const
    {
        if ( transpose )
            return M_Aqm[q][m]->transpose();

        return M_Aqm[q][m];
    }


    /**
     * \brief Returns the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (time dependent)
     *
     * \param q and m are index of the component in the affine decomposition
     * \param transpose transpose \c M_q
     *
     * \return the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (ime dependent)
     */
    sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false ) const
    {
        if ( transpose )
            return M_Mqm[q][m]->transpose();

        return M_Mqm[q][m];
    }


    vector_ptrtype MFqm( uint16_type q, uint16_type m ) const
    {
        return M_MFqm[q][m];
    }


    /**
     * \brief the inner product \f$a_{qm}(\xi_i, \xi_j) = \xi_j^T A_{qm} \xi_i\f$
     *
     * \param q and m index of the component in the affine decomposition
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c A_{qm}
     *
     * \return the inner product \f$a_qm(\xi_i, \xi_j) = \xi_j^T A_{qm} \xi_i\f$
     */
    value_type Aqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false )
    {
        return M_Aqm[q][m]->energy( xi_j, xi_i, transpose );
    }

    /**
     * \brief the inner product \f$m_{qm}(\xi_i, \xi_j) = \xi_j^T M_{qm} \xi_i\f$
     *
     * \param q and m index of the component in the affine decomposition
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c M_{qm}
     *
     * \return the inner product \f$m_{qm}(\xi_i, \xi_j) = \xi_j^T M_{qm} \xi_i\f$
     */
    value_type Mqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false )
    {
        return M_Mqm[q][m]->energy( xi_j, xi_i, transpose );
    }


    /**
     * \brief Returns the vector coefficients
     */
    beta_vector_type const& betaAqm() const
    {
        return M_model->betaAqm();
    }

    /**
     * \brief Returns the vector coefficients
     */
    beta_vector_type const& betaMqm() const
    {
        return mpl::bool_<model_type::is_time_dependent>();
    }
    beta_vector_type const& betaMqm( mpl::bool_<true> ) const
    {
        return M_model->betaMqm();
    }
    beta_vector_type const& betaMqm( mpl::bool_<false> ) const
    {
        beta_vector_type vect;
        return  vect;
    }
    beta_vector_type const& betaInitialGuessQm( mpl::bool_<true> ) const
    {
        return M_model->betaInitialGuessQm();
    }
    /**
     * \brief Returns the value of the \f$A_{qm}\f$ coefficient at \f$\mu\f$
     */
    value_type betaAqm( int q, int m ) const
    {
        return M_model->betaAqm(q,m);
    }

    /**
     * \brief Returns the value of the \f$M_{qm}\f$ coefficient at \f$\mu\f$
     */
    value_type betaMqm( int q, int m ) const
    {
        return betaMqm( q, m, mpl::bool_<model_type::is_time_dependent>() );
    }
    value_type betaMqm( int q, int m, mpl::bool_<false> ) const
    {
        return 0;
    }
    value_type betaMqm( int q, int m, mpl::bool_<true> ) const
    {
        return M_model->betaMqm(q,m);
    }


    value_type betaInitialGuessQm( int q, int m ) const
    {
        return M_model->betaInitialGuessQm(q,m);
    }

    /**
     * \brief Returns the vector coefficients
     */
    std::vector<beta_vector_type> const& betaL() const
    {
        return M_model->betaL();
    }

    /**
     * \brief Returns the value of the \f$L_qm\f$ coefficient at \f$\mu\f$
     */
    value_type betaL( int l, int q, int m ) const
    {
        return M_model->betaL( l, q , m);
    }


    /**
     * \brief the vector \c Fq[q][m] of the affine decomposition of the right hand side
     *
     * \return the vector associated with \f$F_{qm}\f$
     */
    vector_ptrtype Fqm( uint16_type l, uint16_type q, int m ) const
    {
        return M_Fqm[l][q][m];
    }
    element_ptrtype  InitialGuessQm( uint16_type q, int m ) const
    {
        return InitialGuessQm[q][m];
    }

    /**
     * \brief the inner product \f$f_{qm}(\xi) = \xi^T F_{qm} \f$
     *
     * Denote \f$F_{qm}\f$ the algebraic representation of the linear form associated
     * with the right hand side.
     *
     * \param q and m index of the component in the affine decomposition
     * \param xi an element of the function space
     *
     * \return the inner product \f$f_{qm}(\xi) = \xi^T F_{qm} \f$
     */
    value_type Fqm( uint16_type l, uint16_type q,  uint16_type m, element_type const& xi )
    {
        element_ptrtype eltF( new element_type( M_model->functionSpace() ) );
        for(int i=0; i<eltF->localSize();i++)
            eltF->operator()(i)=M_Fqm[l][q][m]->operator()(i);
        return inner_product( *eltF , xi );
        //return inner_product( *M_Fqm[l][q][m], xi );
    }

    /**
     * \brief the inner product \f$f_{qm}(\xi) = \xi^T F_{qm} \f$
     *
     * Denote \f$F_{qm}\f$ the algebraic representation of the linear form associated
     * with the right hand side.
     *
     * \param q and m index of the component in the affine decomposition
     * \param xi an element of the function space
     *
     * \return the inner product \f$f_{qm}(\xi) = \xi^T F_{qm} \f$
     */
    value_type Fqm( uint16_type l, uint16_type q, uint16_type m, element_ptrtype const& xi )
    {
        return inner_product( M_Fqm[l][q][m], xi );
    }

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& X, vector_type const& Y )
    {
        return M_model->scalarProduct( X, Y );
    }
    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        return M_model->scalarProduct( X, Y );
    }


    /**
     * returns the scalar product used to assemble POD matrix of the vector x and vector y
     */
    double scalarProductForPod( vector_type const& X, vector_type const& Y )
    {
        return scalarProductForPod( X, Y ,mpl::bool_<model_type::is_time_dependent>() );
    }
    double scalarProductForPod( vector_type const& X, vector_type const& Y , mpl::bool_<true> )
    {
        return M_model->scalarProductForPod( X, Y );
    }
    double scalarProductForPod( vector_type const& X, vector_type const& Y , mpl::bool_<false> )
    {
        return 0;
    }


    /**
     * returns the scalar product used to assemble POD matrix of the vector x and vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        return scalarProductForPod( X, Y , mpl::bool_<model_type::is_time_dependent>() );
    }
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y , mpl::bool_<true> )
    {
        return M_model->scalarProductForPod( X, Y );
    }
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y , mpl::bool_<false> )
    {
        return 0;
    }


    /**
     * solve the model for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu )
    {
        return M_model->solve( mu );
    }

    /**
     * solve the model for a given parameter \p mu
     */
    element_type solveRB( parameter_type const& mu )
    {
        return M_model->solveRB( mu );
    }

    /**
     * solve the model for a given parameter \p mu
     */
    void solve( parameter_type const& mu, element_ptrtype& u )
    {
        return M_model->solve( mu, u );
    }

    /**
     * solve the model for a given parameter \p mu and \p L as right hand side
     * \param transpose if true solve the transposed(dual) problem
     */
    void solve( parameter_type const& mu, element_ptrtype& u, vector_ptrtype const& L, bool transpose = false )
    {
        return M_model->solve( mu, u, L, transpose );
    }

    /**
     * solve \f$M u = f\f$ where \f$ M \f$ is the matrix associated to the \f$ L_2 \f$
     * norm
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f )
    {
        return M_model->l2solve( u, f );
    }


    /**
     * \brief solve \f$A x = f\f$
     *
     * \note if \p tranpose is true then solve for \f$A^T x = f\f$.
     *
     * \param A matrix
     * \param x solution vector
     * \param f right hand side vector
     * \param transpose if is true solve for \f$A^T x = f\f$, otherwise solve \f$A x = f\f$.
     *
     * \return a tuple with the number of iterations used and the residual
     */
    boost::tuple<int, value_type>
    solve( sparse_matrix_ptrtype const& A,
           vector_ptrtype & x,
           vector_ptrtype const& f,
           bool tranpose = false
         )
    {
        //return M_model->solve( A, x, f );
    }

    /**
     * \brief export a vector of elements
     *
     * \param v a vector of \c shared_ptr<> elements of the functions space
     */
#if 0
    bool
    exportResults( double time, std::vector<element_ptrtype> const& v )
    {
        //return model->export( time, v );
    }
#endif
    /**
     * run the model
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );
    std::map<double, boost::tuple<double,double,double,int, double,double> > run();
    std::map<double, boost::tuple<double,double,double,int, double,double> > run( mpl::bool_<true> );
    std::map<double, boost::tuple<double,double,double,int, double,double> > run( mpl::bool_<false> );
    std::map<double, boost::tuple<double,double,double,int,double,double> > runq();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu )
    {
        return M_model->output( output_index, mu );
    }

    int computeNumberOfSnapshots()
    {
        return computeNumberOfSnapshots( mpl::bool_<model_type::is_time_dependent>() );
    }
    int computeNumberOfSnapshots( mpl::bool_<true> )
    {
        return M_model->computeNumberOfSnapshots();
    }
    int computeNumberOfSnapshots( mpl::bool_<false> )
    {
        return 1;
    }


    double timeStep()
    {
        return timeStep( mpl::bool_<model_type::is_time_dependent>() );
    }
    double timeStep( mpl::bool_<true> )
    {
        double timestep;

        if ( M_model->isSteady() ) timestep=1e30;

        else timestep = M_model->timeStep();

        return timestep;
    }
    double timeStep( mpl::bool_<false> )
    {
        return 1e30;
    }

    double timeInitial()
    {
        return timeInitial( mpl::bool_<model_type::is_time_dependent>() );
    }
    double timeInitial( mpl::bool_<true> )
    {
        return M_model->timeInitial();
    }
    double timeInitial( mpl::bool_<false> )
    {
        return 0;
    }

    double timeFinal()
    {
        return timeFinal( mpl::bool_<model_type::is_time_dependent>() );
    }
    double timeFinal( mpl::bool_<true> )
    {
        double timefinal;

        if ( M_model->isSteady() ) timefinal=1e30;

        else timefinal = M_model->timeFinal();

        return timefinal;
    }
    double timeFinal( mpl::bool_<false> )
    {
        return 1e30;
    }

    int timeOrder()
    {
        return timeOrder( mpl::bool_<model_type::is_time_dependent>() );
    }
    int timeOrder( mpl::bool_<true> )
    {
        return M_model->timeOrder();
    }
    int timeOrder( mpl::bool_<false> )
    {
        return 0;
    }


    bool isSteady()
    {
        return isSteady( mpl::bool_<model_type::is_time_dependent>() );
    }
    bool isSteady( mpl::bool_<true> )
    {
        return M_model->isSteady();
    }
    bool isSteady( mpl::bool_<false> )
    {
        return true;
    }


    void initializationField( element_ptrtype& initial_field,parameter_type const& mu )
    {
        return initializationField( initial_field,mu,mpl::bool_<model_type::is_time_dependent>() );
    }
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu,mpl::bool_<true> )
    {
        return M_model->initializationField( initial_field,mu );
    }
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu,mpl::bool_<false> ) {};


    //@}



protected:


    //! affine decomposition terms for the left hand side
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;

    //! affine decomposition terms ( time dependent )
    mutable std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;

    mutable std::vector< std::vector<vector_ptrtype> > M_MFqm;

    //! affine decomposition terms for the right hand side
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

private:

    bool M_is_initialized;

    //! variables_map
    po::variables_map M_vm;

    //! mode for CRBModel
    CRBModelMode M_mode;

    //! model
    model_ptrtype M_model;

    backend_ptrtype M_backend;

    // ! matrix associated with inner product
    sparse_matrix_ptrtype M_B;
    sparse_matrix_ptrtype M_H1;


    //! initialize the matrix associated with the \f$H_1\f$ inner product
    void initB();

    /**
     * \brief given \p mu merge the Aq and Fq into A and F respectively
     *
     * \param mu the parameter at which the matrix A and vector F are assembled
     *
     */
    offline_merge_type offlineMerge( parameter_type const& mu );

    void assembleMassMatrix( );
    void assembleMassMatrix( mpl::bool_<true> );
    void assembleMassMatrix( mpl::bool_<false> );

    void assembleMF( initial_guess_type & initial_guess );
    void assembleMF( initial_guess_type & initial_guess, mpl::bool_<true> );
    void assembleMF( initial_guess_type & initial_guess, mpl::bool_<false> );

};




template <typename ModelType>
struct AssembleMassMatrixInCompositeCase
{

    //! mesh type
    typedef typename ModelType::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;

    //! function space type
    typedef typename ModelType::functionspace_type functionspace_type;
    typedef typename ModelType::functionspace_ptrtype functionspace_ptrtype;

    //! element of the functionspace type
    typedef typename ModelType::element_type element_type;
    typedef typename ModelType::element_ptrtype element_ptrtype;


    AssembleMassMatrixInCompositeCase( element_type const u ,
                                       element_type const v ,
                                       CRBModel<ModelType> & crb_model)
        :
        M_composite_u ( u ),
        M_composite_v ( v ),
        M_crb_model( crb_model )
    {}

    template< typename T >
    void
    operator()( const T& t ) const
    {

        using namespace Feel::vf;

        auto u = M_composite_u.template element< T::value >();
        auto v = M_composite_v.template element< T::value >();
        auto Xh = M_composite_u.functionSpace();
        mesh_ptrtype mesh = Xh->mesh();

        form2( _test=Xh, _trial=Xh, _matrix=M_crb_model.Mqm(0,0) ) +=
               integrate( _range=elements( mesh ), _expr=idt( u )*id( v ) );

    }

    element_type  M_composite_u;
    element_type  M_composite_v;
    CRBModel<ModelType> & M_crb_model;
};


template <typename ModelType>
struct AssembleMFInCompositeCase
{

    //! mesh type
    typedef typename ModelType::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;

    //! function space type
    typedef typename ModelType::functionspace_type functionspace_type;
    typedef typename ModelType::functionspace_ptrtype functionspace_ptrtype;

    //! element of the functionspace type
    typedef typename ModelType::element_type element_type;
    typedef typename ModelType::element_ptrtype element_ptrtype;

    typedef typename std::vector< std::vector < element_ptrtype > > initial_guess_type;

    AssembleMFInCompositeCase( element_type  const v ,
                               initial_guess_type  const initial_guess ,
                               CRBModel<ModelType> & crb_model)
        :
        M_composite_v ( v ),
        M_composite_initial_guess ( initial_guess ),
        M_crb_model ( crb_model )
    {}

    template< typename T >
    void
    operator()( const T& t ) const
    {

        auto v = M_composite_v.template element< T::value >();
        auto Xh = M_composite_v.functionSpace();
        mesh_ptrtype mesh = Xh->mesh();

        for(int q = 0; q < M_crb_model.Qmf(); q++)
        {
            for( int m = 0; m < M_crb_model.mMaxMF(q); m++)
            {
                auto vectFM = M_crb_model.MFqm(q,m);
                auto ini = M_composite_initial_guess[q][m]->template element< T::value >();
                form1( _test=Xh, _vector=vectFM ) +=
                    integrate ( _range=elements( mesh ), _expr=Feel::vf::idv( ini )*Feel::vf::id( v ) );
            }
        }
    }

    element_type  M_composite_v;
    initial_guess_type  M_composite_initial_guess;
    CRBModel<ModelType> & M_crb_model;
};




template<typename TruthModelType>
void
CRBModel<TruthModelType>::initB()
{

    //the matrix associated with H1 scalar product is now given by the model
    M_B = M_model->innerProduct();

#if 0
    LOG(INFO) << "[CRBModel::initB] initialize scalar product\n";
    M_B = M_backend->newMatrix( M_model->functionSpace(), M_model->functionSpace() );
    using namespace Feel::vf;
    typename functionspace_type::element_type u( M_model->functionSpace() );
    form2( M_model->functionSpace(), M_model->functionSpace(), M_B, _init=true ) =
        integrate( elements( M_model->functionSpace()->mesh() ),
                   gradt( u )*trans( grad( u ) ) );

    M_B->close();

    auto M = M_backend->newMatrix( M_model->functionSpace(), M_model->functionSpace() );
    form2( M_model->functionSpace(), M_model->functionSpace(), M, _init=true ) =
        integrate( elements( M_model->functionSpace()->mesh() ),
                   idt( u )*id( u ) );
    M_B->printMatlab( "ipB.m" );
    M->printMatlab( "ipM.m" );
    M->close();
    LOG(INFO) << "[CRBModel::initB] starting eigen solve\n";
#if 0
    SolverEigen<double>::eigenmodes_type modesmin=
        eigs( _matrixA=M_B,
              _matrixB=M,
              _problem=( EigenProblemType )GHEP,
              _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].as<int>(),
              //_spectrum=LARGEST_MAGNITUDE,
              _spectrum=SMALLEST_MAGNITUDE,
              //_transform=SINVERT,
              _ncv=M_vm["solvereigen-ncv"].as<int>(),
              _nev=M_vm["solvereigen-nev"].as<int>(),
              _tolerance=M_vm["solvereigen-tol"].as<double>(),
              _maxit=M_vm["solvereigen-maxiter"].as<int>()
            );
    double eigmin = 1;

    if ( modesmin.empty() || modesmin.begin()->second.get<0>()<1e-6 )
    {
        LOG(INFO) << "coercivity constant not computed, taking 1\n";
    }

    else
    {
        eigmin = modesmin.begin()->second.get<0>();
    }

    LOG(INFO) << "[CRBModel::initB] coercivity constant (tau) = " << eigmin << "\n";
#else
    double eigmin = 1;
#endif
    M_B->addMatrix( eigmin, M );

#endif
}


template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMassMatrix( void )
{
    static const bool is_composite = functionspace_type::is_composite;
    return assembleMassMatrix( mpl::bool_< is_composite >() );
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMassMatrix( mpl::bool_<false> )
{
    using namespace Feel::vf;
    auto Xh = M_model->functionSpace();
    M_Mqm.resize( 1 );
    M_Mqm[0].resize( 1 );
    M_Mqm[0][0] = M_backend->newMatrix( _test=Xh , _trial=Xh );
    auto mesh = Xh->mesh();
    element_type u ( Xh , "u" );
    element_type v ( Xh , "v" );
    form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=idt( u )*id( v )  );
    M_Mqm[0][0]->close();
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMassMatrix( mpl::bool_<true> )
{
    auto Xh = M_model->functionSpace();

    element_type u ( Xh , "u" );
    element_type v ( Xh , "v" );

    index_vector_type index_vector;
    int size = functionspace_type::nSpaces;
    M_Mqm.resize( 1 );
    M_Mqm[0].resize(1);
    M_Mqm[0][0]=M_backend->newMatrix( _test=Xh , _trial=Xh );

    AssembleMassMatrixInCompositeCase<TruthModelType> assemble_mass_matrix_in_composite_case ( u , v , *this);
    fusion::for_each( index_vector, assemble_mass_matrix_in_composite_case );

    M_Mqm[0][0]->close();
}



template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMF( initial_guess_type & initial_guess )
{
    static const bool is_composite = functionspace_type::is_composite;
    return assembleMF( initial_guess, mpl::bool_< is_composite >() );
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMF( initial_guess_type & initial_guess, mpl::bool_<true> )
{
    auto Xh = M_model->functionSpace();
    auto mesh = Xh->mesh();
    element_type v ( Xh , "v" );

    int q_max= this->Qmf();
    M_MFqm.resize( q_max );
    for(int q = 0; q < q_max; q++ )
    {
        int m_max= this->mMaxMF(q);
        M_MFqm[q].resize( m_max );
        for(int m = 0; m < m_max; m++ )
            M_MFqm[q][m] = M_backend->newVector( Xh );
    }


    index_vector_type index_vector;
    AssembleMFInCompositeCase<TruthModelType> assemble_mf_in_composite_case ( v , initial_guess , *this);
    fusion::for_each( index_vector, assemble_mf_in_composite_case );


    for(int q = 0; q < q_max; q++ )
    {
        for(int m = 0; m < this->mMaxMF(q); m++ )
            M_MFqm[q][m]->close();
    }

}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMF( initial_guess_type & initial_guess, mpl::bool_<false> )
{
    using namespace Feel::vf;
    auto Xh = M_model->functionSpace();
    auto mesh = Xh->mesh();
    element_type v ( Xh , "v" );

    int q_max= this->Qmf();
    M_MFqm.resize( q_max );
    for(int q = 0; q < q_max; q++ )
    {
        int m_max= this->mMaxMF(q);
        M_MFqm[q].resize( m_max );
        for(int m = 0; m < m_max; m++ )
        {
            M_MFqm[q][m] = M_backend->newVector( Xh );
            std::cout<<"M_MFqm["<<q<<"]["<<m<<"]->size() : "<<M_MFqm[q][m]->size()<<std::endl;
            std::cout<<"initial_guess["<<q<<"]["<<m<<"] : "<<initial_guess[q][m]->size()<<std::endl;
            form1( _test=Xh, _vector=M_MFqm[q][m]) =
                integrate( _range=elements( mesh ), _expr=idv( initial_guess[q][m] )*id( v )  );
            std::cout<<"bahhhh"<<std::endl;
            M_MFqm[q][m]->close();
        }
    }
}



template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_ptrtype
CRBModel<TruthModelType>::initialGuess( parameter_type const& mu )
{
    return initialGuess( mu , mpl::bool_<model_type::is_time_dependent>() );
}
template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_ptrtype
CRBModel<TruthModelType>::initialGuess( parameter_type const& mu , mpl::bool_<false> )
{
    auto Xh = M_model->functionSpace();
    element_ptrtype initial_guess = Xh->elementPtr();
    initial_guess_type initial_guess_vector;
    boost::tie( boost::tuples::ignore, boost::tuples::ignore, initial_guess_vector ) = M_model->computeAffineDecomposition();

    for ( size_type q = 0; q < Qmf(); ++q )
    {
        for ( size_type m = 0; m < mMaxMF(q); ++m )
        {
            element_type temp = Xh->element();
            temp = *initial_guess_vector[q][m];
            temp.scale( this->betaInitialGuessQm( q , m ) );
            *initial_guess += temp;
        }
    }
    return initial_guess;
}
template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_ptrtype
CRBModel<TruthModelType>::initialGuess( parameter_type const& mu , mpl::bool_<true> )
{
    auto Xh = M_model->functionSpace();
    element_ptrtype initial_guess = Xh->elementPtr();
    initial_guess_type initial_guess_vector;
    boost::tie( boost::tuples::ignore, boost::tuples::ignore, boost::tuples::ignore, initial_guess_vector ) = M_model->computeAffineDecomposition();

    for ( size_type q = 0; q < Qmf(); ++q )
    {
        for ( size_type m = 0; m < mMaxMF(q); ++m )
        {
            element_type temp = Xh->element();
            temp = *initial_guess_vector[q][m];
            temp.scale( this->betaInitialGuessQm( q , m ) );
            *initial_guess += temp;
        }
    }
    return initial_guess;
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::offline_merge_type
CRBModel<TruthModelType>::offlineMerge( parameter_type const& mu )
{

#if 0
    sparse_matrix_ptrtype A( M_backend->newMatrix(
                                                  _test=M_model->functionSpace(),
                                                  _trial=M_model->functionSpace()
                                                  ) );

    sparse_matrix_ptrtype M( M_backend->newMatrix(
                                                  _test=M_model->functionSpace(),
                                                  _trial=M_model->functionSpace()
                                                  ) );
    vector_ptrtype MF( M_backend->newVector(M_model->functionSpace()) );

#else

    auto A = this->newMatrix();
    auto M = this->newMatrix();
#endif

    std::vector<vector_ptrtype> F( Nl() );


    for ( size_type q = 0; q < Qa(); ++q )
    {
        for(size_type m = 0; m < mMaxA(q); ++m )
            A->addMatrix( this->betaAqm( q , m ), M_Aqm[q][m] );
    }

    if( Qm() > 0 )
    {
        for ( size_type q = 0; q < Qm(); ++q )
        {
            for(size_type m = 0; m < mMaxM(q) ; ++m )
                M->addMatrix( this->betaMqm( q , m ), M_Mqm[q][m] );
        }
    }

    //element_ptrtype InitialGuess = M_model->functionSpace()->elementPtr() ;
    vector_ptrtype MF = M_backend->newVector( M_model->functionSpace() ) ;

    if ( Qmf() > 0 )
    {
        for ( size_type q = 0; q < Qmf(); ++q )
        {
            for ( size_type m = 0; m < mMaxMF(q); ++m )
            {
                MF->add( this->betaInitialGuessQm( q , m ), *M_MFqm[q][m] );
            }
        }
    }

    for ( size_type l = 0; l < Nl(); ++l )
    {
        F[l] = M_backend->newVector( M_model->functionSpace() );
        F[l]->close();
        F[l]->zero();

        for ( size_type q = 0; q < Ql( l ); ++q )
        {
            for ( size_type m = 0; m < mMaxF(l,q); ++m )
                F[l]->add( this->betaL( l, q , m ), M_Fqm[l][q][m] );
        }
    }

    return boost::make_tuple( M, A, F, MF );
}



template<typename TruthModelType>
void
CRBModel<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    M_model->run( X, N, Y, P );
}



template<typename TruthModelType>
std::map<double, boost::tuple<double,double,double,int,double,double> >
CRBModel<TruthModelType>::run()
{
    return run( mpl::bool_<model_type::is_time_dependent>() );
}

template<typename TruthModelType>
std::map<double, boost::tuple<double,double,double,int,double,double> >
CRBModel<TruthModelType>::run( mpl::bool_<true> )
{
    // parameter space
    parameterspace_ptrtype Dmu = M_model->parameterSpace();

    // sampling of parameter space
    sampling_ptrtype Xi( new sampling_type( Dmu ) );
    Xi->equidistribute( 10 );
    parameter_type mu( Dmu );

    this->computeAffineDecomposition();

    std::map<double, boost::tuple<double,double,double,int,double,double> > res;
    BOOST_FOREACH( mu, *Xi )
    {
        // first for a particular parameter
        M_model->solve( mu );

        // solve for an eigenproblem

        sparse_matrix_ptrtype M,A,B,symmA;
        vector_ptrtype F;

        boost::tie( M, A, F ) = this->update( mu );


        A->printMatlab( "run_A.m" );
        symmA = M_model->newMatrix();
        A->symmetricPart( symmA );
        symmA->printMatlab( "run_symmA.m" );
        B = this->innerProduct();
        B->printMatlab( "run_B.m" );
        int nconv;
        boost::timer t_emin;
        double eigenvalue_real, eigenvalue_imag;
        vector_ptrtype eigenvector;
        // solve  for eigenvalue problem at \p mu
        SolverEigen<double>::eigenmodes_type modes=
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=LARGEST_MAGNITUDE,
                  _transform=SINVERT,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );

        double time_eigmin = t_emin.elapsed();
        double eigmin = modes.begin()->second.template get<0>();

        boost::timer t_emax;
        modes=
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=LARGEST_MAGNITUDE,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );
        double time_eigmax = t_emax.elapsed();
        double eigmax = modes.rbegin()->second.template get<0>();

        std::cout << std::setprecision( 4 ) << mu << " "
                  << std::setprecision( 16 ) << eigmin << " " << eigmax
                  << " " << functionSpace()->nDof() << "\n";

        res[mu( 0 )] = boost::make_tuple( mu( 0 ), eigmin, eigmax, functionSpace()->nDof(), time_eigmin, time_eigmax );
    }

    return res;
}

template<typename TruthModelType>
std::map<double, boost::tuple<double,double,double,int,double,double> >
CRBModel<TruthModelType>::run( mpl::bool_<false> )
{
    // parameter space
    parameterspace_ptrtype Dmu = M_model->parameterSpace();

    // sampling of parameter space
    sampling_ptrtype Xi( new sampling_type( Dmu ) );
    Xi->equidistribute( 10 );
    parameter_type mu( Dmu );

    this->computeAffineDecomposition();

    std::map<double, boost::tuple<double,double,double,int,double,double> > res;
    BOOST_FOREACH( mu, *Xi )
    {
        // first for a particular parameter
        M_model->solve( mu );

        // solve for an eigenproblem

        sparse_matrix_ptrtype A,B,symmA;
        vector_ptrtype F;

        boost::tie( A, F ) = this->update( mu );

        A->printMatlab( "run_A.m" );
        symmA = M_model->newMatrix();
        A->symmetricPart( symmA );
        symmA->printMatlab( "run_symmA.m" );
        B = this->innerProduct();
        B->printMatlab( "run_B.m" );
        int nconv;
        boost::timer t_emin;
        double eigenvalue_real, eigenvalue_imag;
        vector_ptrtype eigenvector;
        // solve  for eigenvalue problem at \p mu
        SolverEigen<double>::eigenmodes_type modes=
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=LARGEST_MAGNITUDE,
                  _transform=SINVERT,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );

        double time_eigmin = t_emin.elapsed();
        double eigmin = modes.begin()->second.template get<0>();

        boost::timer t_emax;
        modes=
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=LARGEST_MAGNITUDE,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );
        double time_eigmax = t_emax.elapsed();
        double eigmax = modes.rbegin()->second.template get<0>();

        std::cout << std::setprecision( 4 ) << mu << " "
                  << std::setprecision( 16 ) << eigmin << " " << eigmax
                  << " " << functionSpace()->nDof() << "\n";

        res[mu( 0 )] = boost::make_tuple( mu( 0 ), eigmin, eigmax, functionSpace()->nDof(), time_eigmin, time_eigmax );
    }

    return res;
}

template<typename TruthModelType>
std::map<double, boost::tuple<double,double,double,int,double,double> >
CRBModel<TruthModelType>::runq()
{
    // parameter space
    parameterspace_ptrtype Dmu = M_model->parameterSpace();

    // sampling of parameter space
    sampling_ptrtype Xi( new sampling_type( Dmu ) );
    Xi->equidistribute( 10 );
    parameter_type mu( Dmu );

    mu = Xi->min().template get<0>();
    boost::tie( M_Aqm, M_Fqm ) = this->computeAffineDecomposition();
    M_model->solve( mu );
    std::cout << "done with setup\n";
    std::map<double, boost::tuple<double,double,double,int,double,double> > res;

    for ( int q = 0; q < M_Aqm.size(); ++q )
    {
        std::cout << "================================================================================\n";
        std::cout << "= q = " << q << " / " << M_Aqm.size() << "\n";
        sparse_matrix_ptrtype B,symmA;
        symmA = M_model->newMatrix();
        M_Aqm[q]->symmetricPart( symmA );
        B = this->innerProduct();
        std::ostringstream os;
        os << "symmAq_" << q << ".m";
        symmA->printMatlab( os.str() );
        os.str( "" );
        os << "B_" << q << ".m";
        B->printMatlab( os.str() );
        os.str( "" );
        os << "Aq_" << q << ".m";
        M_Aqm[q]->printMatlab( os.str() );
        int nconv;
        double eigenvalue_real, eigenvalue_imag;
        vector_ptrtype eigenvector;
        SolverEigen<double>::eigenmodes_type modes;
        boost::timer t_emin;
#if 1
        // solve  for eigenvalue problem at \p mu
        modes =
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=SMALLEST_MAGNITUDE,
                  //_transform=SINVERT,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );
#endif
        double time_eigmin = t_emin.elapsed();
        double eigmin = modes.empty()?0:modes.begin()->second.template get<0>();
        boost::timer t_emax;
        modes=
            eigs( _matrixA=symmA,
                  _matrixB=B,
                  //_problem=(EigenProblemType)GHEP,
                  _solver=( EigenSolverType )M_vm["solvereigen-solver-type"].template as<int>(),
                  _spectrum=LARGEST_MAGNITUDE,
                  _ncv=M_vm["solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["solvereigen-maxiter"].template as<int>()
                );
        double time_eigmax = t_emax.elapsed();
        double eigmax = modes.empty()?0:modes.rbegin()->second.template get<0>();

        std::cout << std::setprecision( 4 ) << mu << " "
                  << std::setprecision( 16 ) << eigmin << " " << eigmax
                  << " " << functionSpace()->nDof() << "\n";

        res[double( q )] = boost::make_tuple( double( q ), eigmin, eigmax, functionSpace()->nDof(), time_eigmin, time_eigmax );
    }

    return res;
}
}
#endif /* __CRBModel_H */
