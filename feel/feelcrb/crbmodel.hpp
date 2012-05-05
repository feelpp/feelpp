/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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


    typedef Eigen::VectorXd theta_vector_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype, sparse_matrix_ptrtype, std::vector<vector_ptrtype> > offline_merge_type;


    typedef typename boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<sparse_matrix_ptrtype>, std::vector< std::vector<vector_ptrtype> > > affine_decomposition_type;

    typedef typename boost::tuple<theta_vector_type,theta_vector_type,std::vector<theta_vector_type> > thetaq_type;




    //@}

    /** @name Constructors, destructor
     */
    //@{

    CRBModel()
        :
        M_Aq(),
        M_Mq(),
        M_Fq(),
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
        M_Aq(),
        M_Mq(),
        M_Fq(),
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
        M_Aq(),
        M_Mq(),
        M_Fq(),
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
        M_Aq( o.M_Aq ),
        M_Mq( o.M_Mq ),
        M_Fq( o.M_Fq ),
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
    ~CRBModel()
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
            std::cout << "  -- init FEM  model\n";
            M_model->init();
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
            M_Aq = o.M_Aq;
            M_Fq = o.M_Fq;
            M_Mq = o.M_Mq;
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
        return 0;
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
     * \brief compute the thetaq given \p mu
     */
    thetaq_type computeThetaq( parameter_type const& mu , double time=0 )
    {
        return computeThetaq( mu , mpl::bool_<model_type::is_time_dependent>(), time  );
    }
    thetaq_type computeThetaq( parameter_type const& mu , mpl::bool_<true>, double time=0 )
    {
        return M_model->computeThetaq( mu , time );
    }
    thetaq_type computeThetaq( parameter_type const& mu , mpl::bool_<false>, double time=0 )
    {
        theta_vector_type theta_aq;
        theta_vector_type theta_mq;
        std::vector<theta_vector_type>  theta_fq;
        boost::tuple<theta_vector_type, std::vector<theta_vector_type> > steady_theta;
        steady_theta = M_model->computeThetaq( mu , time );
        theta_aq = steady_theta.get<0>();
        theta_fq = steady_theta.get<1>();
        return boost::make_tuple( theta_mq, theta_aq, theta_fq );
    }



    /**
     * \brief update the model wrt \p mu
     */
    offline_merge_type update( parameter_type const& mu,  double time=0 )
    {
        M_model->computeThetaq( mu , time );
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
        boost::tie( M_Mq, M_Aq, M_Fq ) = M_model->computeAffineDecomposition();
        return M_model->computeAffineDecomposition();
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        boost::tie( M_Aq, M_Fq ) = M_model->computeAffineDecomposition();
        return boost::make_tuple( M_Mq, M_Aq, M_Fq );
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
     * \brief Returns the matrix \c Aq[q] of the affine decomposition of the bilinear form
     *
     * \param q the index of the component in the affine decomposition
     * \param transpose transpose \c A_q
     *
     * \return the matrix \c Aq[q] of the affine decomposition of the bilinear form
     */
    sparse_matrix_ptrtype  Aq( uint16_type q, bool transpose = false ) const
    {
        if ( transpose )
            return M_Aq[q]->transpose();

        return M_Aq[q];
    }


    /**
     * \brief Returns the matrix \c Mq[q] of the affine decomposition of the bilinear form (time dependent)
     *
     * \param q the index of the component in the affine decomposition
     * \param transpose transpose \c M_q
     *
     * \return the matrix \c Mq[q] of the affine decomposition of the bilinear form (ime dependent)
     */
    sparse_matrix_ptrtype  Mq( uint16_type q, bool transpose = false ) const
    {
        if ( transpose )
            return M_Mq[q]->transpose();

        return M_Mq[q];
    }


    /**
     * \brief the inner product \f$a_q(\xi_i, \xi_j) = \xi_j^T A_q \xi_i\f$
     *
     * \param q the index of the component in the affine decomposition
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c A_q
     *
     * \return the inner product \f$a_q(\xi_i, \xi_j) = \xi_j^T A_q \xi_i\f$
     */
    value_type Aq( uint16_type q, element_type const& xi_i, element_type const& xi_j, bool transpose = false )
    {
        return M_Aq[q]->energy( xi_j, xi_i, transpose );
    }

    /**
     * \brief the inner product \f$m_q(\xi_i, \xi_j) = \xi_j^T M_q \xi_i\f$
     *
     * \param q the index of the component in the affine decomposition
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c M_q
     *
     * \return the inner product \f$m_q(\xi_i, \xi_j) = \xi_j^T M_q \xi_i\f$
     */
    value_type Mq( uint16_type q, element_type const& xi_i, element_type const& xi_j, bool transpose = false )
    {
        return M_Mq[q]->energy( xi_j, xi_i, transpose );
    }


    /**
     * \brief Returns the vector coefficients
     */
    theta_vector_type const& thetaAq() const
    {
        return M_model->thetaAq();
    }

    /**
     * \brief Returns the vector coefficients
     */
    theta_vector_type const& thetaMq() const
    {
        return mpl::bool_<model_type::is_time_dependent>();
    }
    theta_vector_type const& thetaMq( mpl::bool_<true> ) const
    {
        return M_model->thetaMq();
    }
    theta_vector_type const& thetaMq( mpl::bool_<false> ) const
    {
        theta_vector_type vect;
        return  vect;
    }

    /**
     * \brief Returns the value of the \f$A_q\f$ coefficient at \f$\mu\f$
     */
    value_type thetaAq( int q ) const
    {
        return M_model->thetaAq( q );
    }

    /**
     * \brief Returns the value of the \f$M_q\f$ coefficient at \f$\mu\f$
     */
    value_type thetaMq( int q ) const
    {
        return thetaMq( q,mpl::bool_<model_type::is_time_dependent>() );
    }
    value_type thetaMq( int q, mpl::bool_<true> ) const
    {
        return M_model->thetaMq( q );
    }
    value_type thetaMq( int q, mpl::bool_<false> ) const
    {
        return 0;
    }

    /**
     * \brief Returns the vector coefficients
     */
    std::vector<theta_vector_type> const& thetaL() const
    {
        return M_model->thetaL();
    }

    /**
     * \brief Returns the value of the \f$L_q\f$ coefficient at \f$\mu\f$
     */
    value_type thetaL( int l, int q ) const
    {
        return M_model->thetaL( l, q );
    }


    /**
     * \brief the vector \c Fq[q] of the affine decomposition of the right hand side
     *
     * \return the vector associated with \f$F_q\f$
     */
    vector_ptrtype  Fq( uint16_type l, uint16_type q ) const
    {
        return M_Fq[l][q];
    }

    /**
     * \brief the inner product \f$f_q(\xi) = \xi^T F_q \f$
     *
     * Denote \f$F_q\f$ the algebraic representation of the linear form associated
     * with the right hand side.
     *
     * \param q the index of the component in the affine decomposition
     * \param xi an element of the function space
     *
     * \return the inner product \f$f_q(\xi) = \xi^T F_q \f$
     */
    value_type Fq( uint16_type l, uint16_type q, element_type const& xi )
    {
        return inner_product( *M_Fq[l][q], xi );
    }

    /**
     * \brief the inner product \f$f_q(\xi) = \xi^T F_q \f$
     *
     * Denote \f$F_q\f$ the algebraic representation of the linear form associated
     * with the right hand side.
     *
     * \param q the index of the component in the affine decomposition
     * \param xi an element of the function space
     *
     * \return the inner product \f$f_q(\xi) = \xi^T F_q \f$
     */
    value_type Fq( uint16_type l, uint16_type q, element_ptrtype const& xi )
    {
        return inner_product( M_Fq[l][q], xi );
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
    void solve( parameter_type const& mu )
    {
        return M_model->solve( mu );
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
    std::vector<sparse_matrix_ptrtype> M_Aq;

    //! affine decomposition terms ( time dependent )
    std::vector<sparse_matrix_ptrtype> M_Mq;

    //! affine decomposition terms for the right hand side
    std::vector<std::vector<vector_ptrtype> > M_Fq;


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


};

template<typename TruthModelType>
void
CRBModel<TruthModelType>::initB()
{
    Log() << "[CRBModel::initB] initialize scalar product\n";
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
    Log() << "[CRBModel::initB] starting eigen solve\n";
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
        Log() << "coercivity constant not computed, taking 1\n";
    }

    else
    {
        eigmin = modesmin.begin()->second.get<0>();
    }

    Log() << "[CRBModel::initB] coercivity constant (tau) = " << eigmin << "\n";
#else
    double eigmin = 1;
#endif
    M_B->addMatrix( eigmin, M );
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::offline_merge_type
CRBModel<TruthModelType>::offlineMerge( parameter_type const& mu )
{

#if 1
    sparse_matrix_ptrtype A( M_backend->newMatrix(
                                                  _test=M_model->functionSpace(),
                                                  _trial=M_model->functionSpace()
                                                  ) );

    sparse_matrix_ptrtype M( M_backend->newMatrix(
                                                  _test=M_model->functionSpace(),
                                                  _trial=M_model->functionSpace()
                                                  ) );

#else

    auto A = this->newMatrix();
    auto M = this->newMatrix();

#endif
    std::vector<vector_ptrtype> F( Nl() );

    *A = *M_Aq[0];
    A->scale( this->thetaAq( 0 ) );
    for ( size_type q = 1; q < Qa(); ++q )
    {
        A->addMatrix( this->thetaAq( q ), M_Aq[q] );
    }


    *M = *M_Mq[0];
    M->scale( this->thetaMq( 0 ) );
    for ( size_type q = 1; q < Qm(); ++q )
    {
        M->addMatrix( this->thetaMq( q ), M_Mq[q] );
    }


    for ( size_type l = 0; l < Nl(); ++l )
    {
        F[l] = M_backend->newVector( M_model->functionSpace() );
        F[l]->close();
        F[l]->zero();

        for ( size_type q = 0; q < Ql( l ); ++q )
        {
            F[l]->add( this->thetaL( l, q ), M_Fq[l][q] );
        }

    }

    return boost::make_tuple( M, A, F );
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
    boost::tie( M_Aq, M_Fq ) = this->computeAffineDecomposition();
    M_model->solve( mu );
    std::cout << "done with setup\n";
    std::map<double, boost::tuple<double,double,double,int,double,double> > res;

    for ( int q = 0; q < M_Aq.size(); ++q )
    {
        std::cout << "================================================================================\n";
        std::cout << "= q = " << q << " / " << M_Aq.size() << "\n";
        sparse_matrix_ptrtype B,symmA;
        symmA = M_model->newMatrix();
        M_Aq[q]->symmetricPart( symmA );
        B = this->innerProduct();
        std::ostringstream os;
        os << "symmAq_" << q << ".m";
        symmA->printMatlab( os.str() );
        os.str( "" );
        os << "B_" << q << ".m";
        B->printMatlab( os.str() );
        os.str( "" );
        os << "Aq_" << q << ".m";
        M_Aq[q]->printMatlab( os.str() );
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
