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
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/operatorlinearfree.hpp>
#include <feel/feeldiscr/operatorlinearcomposite.hpp>
#include <feel/feeldiscr/fsfunctionallinearfree.hpp>
#include <feel/feeldiscr/fsfunctionallinearcomposite.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcore/pslogger.hpp>


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
class CRBModel : public boost::enable_shared_from_this<CRBModel<ModelType> >
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

    //! reduced basis function space type
    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename model_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;


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
    typedef Eigen::VectorXd vectorN_type;
    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype,
                                  sparse_matrix_ptrtype,
                                  std::vector<vector_ptrtype>
                                  > offline_merge_type;


    typedef typename boost::tuple<std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector< std::vector<vector_ptrtype> > >
                                  > affine_decomposition_type;


    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > betaqm_type;

    //! time discretization
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef Preconditioner<double> preconditioner_type;
    typedef boost::shared_ptr<preconditioner_type> preconditioner_ptrtype;

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
        M_InitialGuessV(),
        M_InitialGuessVector(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm( vm ),
        M_mode( mode ),
        M_model( new model_type( vm ) ),
        M_backend( backend_type::build( vm ) ),
        M_backend_primal( backend_type::build( vm , "backend-primal" ) ),
        M_backend_dual( backend_type::build( vm , "backend-dual" ) ),
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
        M_InitialGuessV(),
        M_InitialGuessVector(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm(),
        M_mode( CRBModelMode::PFEM ),
        M_model( model ),
        M_backend( backend_type::build( model->vm ) ),
        M_backend_primal( backend_type::build( model->vm , "backend-primal" ) ),
        M_backend_dual( backend_type::build( model->vm , "backend-dual") ),
        M_B()
    {
        this->init();
    }

    CRBModel( model_ptrtype & model , CRBModelMode mode )
        :
        M_Aqm(),
        M_InitialGuessV(),
        M_InitialGuessVector(),
        M_Mqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm(),
        M_mode( mode ),
        M_model( model ),
        M_backend( backend_type::build( Environment::vm() ) ),
        M_backend_primal( backend_type::build( Environment::vm() , "backend-primal" ) ),
        M_backend_dual( backend_type::build( Environment::vm() , "backend-dual" ) ),
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
        M_InitialGuessV( o.M_InitialGuessV ),
        M_InitialGuessVector( o.M_InitialGuessVector ),
        M_Mqm( o.M_Mqm ),
        M_Fqm( o.M_Fqm ),
        M_is_initialized( o.M_is_initialized ),
        M_vm( o.M_vm ),
        M_mode( o.M_mode ),
        M_model(  o.M_model ),
        M_backend( o.M_backend ),
        M_backend_primal( o.M_backend_primal ),
        M_backend_dual( o.M_backend_dual ),
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

        M_preconditioner_primal = preconditioner(_pc=(PreconditionerType) M_backend_primal->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                                 _backend= M_backend_primal,
                                                 _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_primal->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                                 _worldcomm=M_backend_primal->comm(),
                                                 _prefix=M_backend_primal->prefix() ,
                                                 _rebuild=true);
        M_preconditioner_dual = preconditioner(_pc=(PreconditionerType) M_backend_dual->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                               _backend= M_backend_dual,
                                               _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_dual->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                               _worldcomm=M_backend_dual->comm(),
                                               _prefix=M_backend_dual->prefix() ,
                                               _rebuild=true);
        M_is_initialized=true;

        if( ! M_model->isInitialized() )
        {
            LOG( INFO ) << "CRBModel Model is not initialized";
            M_model->initModel();
            M_model->setInitialized( true );
        }

        initB();

        if ( M_mode != CRBModelMode::CRB_ONLINE &&
                M_mode != CRBModelMode::SCM_ONLINE )
        {
            //the model is already initialized
            //std::cout << "  -- init FEM  model\n";
            //M_model->init();
        }

        auto Xh = M_model->functionSpace();
        u = Xh->element();
        v = Xh->element();

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
            M_backend_primal = o.M_backend_primal;
            M_backend_dual = o.M_backend_dual;
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
    sparse_matrix_ptrtype  innerProduct()
    {
        return M_B;
    }

    /**
     * \brief Returns the matrix associated with the inner product
     */
    sparse_matrix_ptrtype const& innerProductForMassMatrix() const
    {
        return M_model->innerProductForMassMatrix();
    }

    /**
     * \brief Returns the matrix associated with the inner product
     */
    sparse_matrix_ptrtype  innerProductForMassMatrix()
    {
        return M_model->innerProductForMassMatrix();
    }




    /**
     * \brief Returns the matrix associated with the \f$H_1\f$ inner product
     */
    sparse_matrix_ptrtype const& h1() const
    {
        return M_B;
    }

    sparse_matrix_ptrtype h1()
    {
        return M_B;
    }

    //!  Returns the function space
    functionspace_ptrtype  functionSpace() const
    {
        return M_model->functionSpace();
    }

    //!  Returns the reduced basis function space
    rbfunctionspace_ptrtype  rBFunctionSpace() const
    {
        return M_model->rBFunctionSpace();
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
        if( M_model->constructOperatorCompositeM() )
            return functionspace_type::nSpaces;
        else
            return 1;
    }

    int QInitialGuess() const
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

    int mMaxInitialGuess( int q ) const
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


    parameter_type refParameter()
    {
        return M_model->refParameter();
    }
    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the mesh characteristic length to \p s
     */
#if 0
    void setMeshSize( double s )
    {
        M_model->setMeshSize( s );
    }
#endif

    //@}

    /** @name  Methods
     */
    //@{

    beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) const
    {
        return M_model->computeBetaInitialGuess( mu );
    }

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
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            steady_beta;

        steady_beta = M_model->computeBetaQm( mu , time );
        betaAqm = steady_beta.get<0>();
        betaFqm = steady_beta.get<1>();

        int nspace = functionspace_type::nSpaces;
        //if model provides implementation of operator composite M
        if ( M_model->constructOperatorCompositeM() )
        {
            betaMqm.resize( nspace );
            for(int q=0; q<nspace; q++)
            {
                betaMqm[q].resize(1);
                betaMqm[q][0] = 1 ;
            }
        }
        else
        {
            betaMqm.resize( 1 );
            betaMqm[0].resize(1);
            betaMqm[0][0] = 1 ;
        }

        return boost::make_tuple( betaMqm, betaAqm, betaFqm );
    }

    betaqm_type computeBetaQm( vector_ptrtype const& T, parameter_type const& mu , double time=0 )
    {
        auto solution = M_model->functionSpace()->element();
        solution = *T;
        return computeBetaQm( solution , mu , time );
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
        beta_vector_type betaAqm, betaMqm ;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            steady_beta;

        steady_beta = M_model->computeBetaQm(T, mu , time );
        betaAqm = steady_beta.get<0>();
        betaFqm = steady_beta.get<1>();

        int nspace = functionspace_type::nSpaces;
        if ( M_model->constructOperatorCompositeM() )
        {
            betaMqm.resize( nspace );
            for(int q=0; q<nspace; q++)
            {
                betaMqm[q].resize(1);
                betaMqm[q][0] = 1 ;
            }
        }
        else
        {
            betaMqm.resize( 1 );
            betaMqm[0].resize(1);
            betaMqm[0][0] = 1 ;
        }


        return boost::make_tuple( betaMqm, betaAqm, betaFqm );
    }


    element_ptrtype assembleInitialGuess( parameter_type const& mu );

    /**
     * \brief update the model wrt \p mu
     */
    offline_merge_type update( parameter_type const& mu,  double time=0 )
    {
        auto all_beta = this->computeBetaQm( mu , time );
        offline_merge_type offline_merge;

        if( option(_name="crb.stock-matrices").template as<bool>() )
            offline_merge = offlineMerge( all_beta , mu );
        else
            offline_merge = offlineMergeOnFly( all_beta, mu );

        return offline_merge;
    }
    offline_merge_type update( parameter_type const& mu, element_type const& T, double time=0 )
    {
#if !defined(NDEBUG)
        mu.check();
#endif
        auto all_beta = this->computeBetaQm(  T , mu , time );

        offline_merge_type offline_merge;

        if( option(_name="crb.stock-matrices").template as<bool>() )
            offline_merge = offlineMerge( all_beta , mu );
        else
            offline_merge = offlineMergeOnFly( all_beta, mu );

        return offline_merge;

    }

    element_type solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu );
    element_type solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu );
    element_type solveFemUsingOfflineEim( parameter_type const& mu );


    /**
     * initialize the model
     */
    void initModel()
    {
        return M_model->initModel();
    }
    void setInitialized( const bool & b)
    {
        return M_model->setInitialized( b );
    }
    bool isInitialized()
    {
        return M_model->isInitialized();
    }



    /**
     * returns list of eim objects ( scalar continuous)
     */
    typename model_type::funs_type scalarContinuousEim()
    {
        return M_model->scalarContinuousEim();
    }

    /**
     * returns list of eim objects ( scalar discontinuous)
     */
    typename model_type::funsd_type scalarDiscontinuousEim()
    {
        return M_model->scalarDiscontinuousEim();
    }


    struct ComputeNormL2InCompositeCase
    {

        ComputeNormL2InCompositeCase( element_type const composite_u1 ,
                                      element_type const composite_u2 )
        :
            M_composite_u1 ( composite_u1 ),
            M_composite_u2 ( composite_u2 )
        {}

        template< typename T >
        void
        operator()( const T& t ) const
        {
            int i = T::value;
            if( i == 0 )
                M_vec.resize( 1 );
            else
                M_vec.conservativeResize( i+1 );

            auto u1 = M_composite_u1.template element< T::value >();
            auto u2 = M_composite_u2.template element< T::value >();
            mesh_ptrtype mesh = u1.functionSpace()->mesh();
            double norm  = normL2(_range=elements( mesh ),_expr=( idv(u1)-idv(u2) ) );
            M_vec(i)= norm ;
        }

        double norm()
        {
            return M_vec.sum();
        }

        mutable vectorN_type M_vec;
        element_type M_composite_u1;
        element_type M_composite_u2;
    };


    double computeNormL2( element_type u1 , element_type u2 )
    {
        static const bool is_composite = functionspace_type::is_composite;
        return computeNormL2( u1, u2 , mpl::bool_< is_composite >() );
    }
    double computeNormL2( element_type u1 , element_type u2 , mpl::bool_<false>)
    {
        double norm=normL2(_range=elements( u2.mesh() ),_expr=( idv(u1)-idv(u2) ) );
        return norm;
    }
    double computeNormL2( element_type u1 , element_type u2 , mpl::bool_<true>)
    {
        //in that case we accumulate the norm of every components
        ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u1, u2 );
        index_vector_type index_vector;
        fusion::for_each( index_vector, compute_normL2_in_composite_case );
        return compute_normL2_in_composite_case.norm();
    }


    /**
     * Access to a pre-assemble objects given by the model
     * The key of the map is a tuple of two double (one for "q" index and other for "m" index)
     * The q^th term of the affine decomposition can have m terms from the EIM
     * Let MAT the matrix associated to index (q,m)
     * MAT can have contributions from elements of the mesh, or faces, or both
     * that is why we have a vector of pre-assemble objects (and not only one pre-assemble object)
     * associated with a tuple <q,m>
     */
    //bilinear form
    operatorcomposite_ptrtype operatorCompositeA() const
    {
        return M_model->operatorCompositeA();
    }
    //linear form
    std::vector< functionalcomposite_ptrtype > functionalCompositeF() const
    {
        return M_model->functionalCompositeF();
    }
    //mass matrix
    operatorcomposite_ptrtype operatorCompositeM() const
    {
        return operatorCompositeM( mpl::bool_<model_type::is_time_dependent>() );
    }
    operatorcomposite_ptrtype operatorCompositeM( mpl::bool_<true> ) const
    {
        return M_model->operatorCompositeM();
    }
    operatorcomposite_ptrtype operatorCompositeM( mpl::bool_<false> ) const
    {
        bool constructed_by_model = M_model->constructOperatorCompositeM();
        if( constructed_by_model )
            return M_model->operatorCompositeM();
        else
            return preAssembleMassMatrix();
    }


    /**
     * \brief Compute the affine decomposition of the various forms
     *
     * This function assembles the parameter independant part of
     * the affine decomposition of the bilinear and linear forms.
     * This function will assemble and stock all matrices/vector
     * associated to the affine decomposition and must be called
     * ONLY if crb.stock-matrices=true
     */
    affine_decomposition_type computeAffineDecomposition()
    {
        return computeAffineDecomposition( mpl::bool_<model_type::is_time_dependent>() );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<true> )
    {
        boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();

        if( M_Aqm.size() == 0 )
        {
            auto compositeM = operatorCompositeM();
            int q_max = this->Qm();
            M_Mqm.resize( q_max);
            for(int q=0; q<q_max; q++)
            {
                int m_max = this->mMaxM(q);
                M_Mqm[q].resize(m_max);
                for(int m=0; m<m_max;m++)
                {
                    auto operatorfree = compositeM->operatorlinear(q,m);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Mqm[q][m]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Mqm[q][m]);//fill the matrix
                }//m
            }//q

            auto compositeA = operatorCompositeA();
            q_max = this->Qa();
            M_Aqm.resize( q_max);
            for(int q=0; q<q_max; q++)
            {
                int m_max = this->mMaxA(q);
                M_Aqm[q].resize(m_max);
                for(int m=0; m<m_max;m++)
                {
                    auto operatorfree = compositeA->operatorlinear(q,m);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Aqm[q][m]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Aqm[q][m]);//fill the matrix
                }//m
            }//q

            auto vector_compositeF = functionalCompositeF();
            int number_outputs = vector_compositeF.size();
            M_Fqm.resize(number_outputs);
            for(int output=0; output<number_outputs; output++)
            {
                auto composite_f = vector_compositeF[output];
                int q_max = this->Ql(output);
                M_Fqm[output].resize( q_max);
                for(int q=0; q<q_max; q++)
                {
                    int m_max = this->mMaxF(output,q);
                    M_Fqm[output][q].resize(m_max);
                    for(int m=0; m<m_max;m++)
                    {
                        auto operatorfree = composite_f->functionallinear(q,m);
                        auto space = operatorfree->space();
                        M_Fqm[output][q][m]= M_backend->newVector( space );
                        operatorfree->containerPtr(M_Fqm[output][q][m]);//fill the vector
                    }//m
                }//q
            }//output
        }

        return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        initial_guess_type initial_guess;
        boost::tie( M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();

        if ( M_Aqm.size() > 0 )
        {
            assembleMassMatrix();
        }
        else
        {
            auto compositeM = operatorCompositeM();
            int q_max = this->Qm();
            M_Mqm.resize( q_max);
            for(int q=0; q<q_max; q++)
            {
                int m_max = this->mMaxM(q);
                M_Mqm[q].resize(m_max);
                for(int m=0; m<m_max;m++)
                {
                    auto operatorfree = compositeM->operatorlinear(q,m);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Mqm[q][m]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Mqm[q][m]);//fill the matrix
                }//m
            }//q

            auto compositeA = operatorCompositeA();
            q_max = this->Qa();
            M_Aqm.resize( q_max);
            for(int q=0; q<q_max; q++)
            {
                int m_max = this->mMaxA(q);
                M_Aqm[q].resize(m_max);
                for(int m=0; m<m_max;m++)
                {
                    auto operatorfree = compositeA->operatorlinear(q,m);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Aqm[q][m]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Aqm[q][m]);//fill the matrix
                }//m
            }//q

            auto vector_compositeF = functionalCompositeF();
            int number_outputs = M_model->Nl();
            M_Fqm.resize(number_outputs);
            for(int output=0; output<number_outputs; output++)
            {
                auto composite_f = vector_compositeF[output];
                int q_max = this->Ql(output);
                M_Fqm[output].resize(q_max);
                for(int q=0; q<q_max; q++)
                {
                    int m_max = this->mMaxF(output,q);
                    M_Fqm[output][q].resize(m_max);
                    for(int m=0; m<m_max;m++)
                    {
                        auto functionalfree = composite_f->functionallinear(q,m);
                        auto space = functionalfree->space();
                        M_Fqm[output][q][m]= M_backend->newVector( space );
                        functionalfree->containerPtr(M_Fqm[output][q][m]);//fill the vector
                    }//m
                }//q
            }//output
        }

        return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
    }


    std::vector< std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition( )
    {
        return M_model->computeInitialGuessAffineDecomposition( );
    }

    std::vector< std::vector<element_ptrtype> > computeInitialGuessVAffineDecomposition( )
    {
        initial_guess_type initial_guess_v;
        initial_guess_v = M_model->computeInitialGuessAffineDecomposition();
        this->assembleInitialGuessV( initial_guess_v);
        for(int q=0; q<M_InitialGuessVector.size(); q++)
        {
            for(int m=0; m<M_InitialGuessVector[q].size(); m++)
                *M_InitialGuessV[q][m] = *M_InitialGuessVector[q][m];
        }
        return M_InitialGuessV;
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
    value_type h1( element_type const& xi_i, element_type const& xi_j  ) const
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
    value_type h1( element_type const& xi_i  ) const
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
    sparse_matrix_ptrtype Aqm( uint16_type q, uint16_type m, bool transpose = false )
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Aqm[q][m]->transpose();

            return M_Aqm[q][m];
        }
        else
        {
            auto composite = operatorCompositeA();
            auto opfree = composite->operatorlinear(q,m);
            size_type pattern = opfree->pattern();
            auto trial = opfree->domainSpace();
            auto test = opfree->dualImageSpace();
            auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
            opfree->matPtr( matrix ); //assemble the matrix

            if ( transpose )
                return matrix->transpose();

            return matrix;
        }
    }


    /**
     * \brief Returns the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (time dependent)
     *
     * \param q and m are index of the component in the affine decomposition
     * \param transpose transpose \c M_q
     *
     * \return the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (ime dependent)
     */
    const sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false ) const
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Mqm[q][m]->transpose();

            return M_Mqm[q][m];
        }
        else
        {
            auto composite = operatorCompositeM();
            auto opfree = composite->operatorlinear(q,m);
            size_type pattern = opfree->pattern();
            auto trial = opfree->domainSpace();
            auto test = opfree->dualImageSpace();
            auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
            opfree->matPtr( matrix ); //assemble the matrix

            if ( transpose )
                return matrix->transpose();

            return matrix;
        }
    }

    /**
     * \brief Returns the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (time dependent)
     *
     * \param q and m are index of the component in the affine decomposition
     * \param transpose transpose \c M_q
     *
     * \return the matrix \c Mq[q][m] of the affine decomposition of the bilinear form (ime dependent)
     */
    sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false )
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Mqm[q][m]->transpose();

            return M_Mqm[q][m];
        }
        else
        {
            auto composite = operatorCompositeM();
            auto opfree = composite->operatorlinear(q,m);
            size_type pattern = opfree->pattern();
            auto trial = opfree->domainSpace();
            auto test = opfree->dualImageSpace();
            auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
            opfree->matPtr( matrix ); //assemble the matrix

            if ( transpose )
                return matrix->transpose();

            return matrix;
        }
    }

    const vector_ptrtype InitialGuessVector( uint16_type q, uint16_type m ) const
    {
        return M_InitialGuessVector[q][m];
    }

    vector_ptrtype InitialGuessVector( uint16_type q, uint16_type m )
    {
        return M_InitialGuessVector[q][m];
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
    value_type Aqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            //in this case matrices have already been stocked
            return M_Aqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        else
        {
            auto composite = operatorCompositeA();
            auto opfree = composite->operatorlinear(q,m);
            size_type pattern = opfree->pattern();
            auto trial = opfree->domainSpace();
            auto test = opfree->dualImageSpace();
            auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
            opfree->matPtr( matrix ); //assemble the matrix

            return matrix->energy( xi_j, xi_i, transpose );
        }
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
    value_type Mqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            //in this case matrices have already been stocked
            return M_Mqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        else
        {
            auto composite = operatorCompositeM();
            auto opfree = composite->operatorlinear(q,m);
            size_type pattern = opfree->pattern();
            auto trial = opfree->domainSpace();
            auto test = opfree->dualImageSpace();
            auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
            opfree->matPtr( matrix ); //assemble the matrix

            return matrix->energy( xi_j, xi_i, transpose );
        }
    }




    beta_vector_type const& betaInitialGuessQm( mpl::bool_<true> ) const
    {
        return M_model->betaInitialGuessQm();
    }


    /**
     * \brief the vector \c Fq[q][m] of the affine decomposition of the right hand side
     *
     * \return the vector associated with \f$F_{qm}\f$
     */
    vector_ptrtype Fqm( uint16_type l, uint16_type q, int m ) const
    {
        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            return M_Fqm[l][q][m];
        }
        else
        {
            auto vector_composite = functionalCompositeF();
            auto composite = vector_composite[l];
            auto functional = composite->functionallinear(q,m);
            auto space = functional->space();
            auto vector = M_backend->newVector( space );
            functional->containerPtr( vector );
            return vector;
        }
    }

    element_ptrtype  InitialGuessQm( uint16_type q, int m ) const
    {
        return M_InitialGuessV[q][m];
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

        value_type result=0;

        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            result = inner_product( *M_Fqm[l][q][m] , xi );
        }
        else
        {
            auto vector_composite = functionalCompositeF();
            auto composite = vector_composite[l];
            auto functional = composite->functionallinear(q,m);
            auto space = functional->space();
            auto vector = M_backend->newVector( space );
            functional->containerPtr( vector );
            result = inner_product( *vector, xi ) ;
        }

        return result;
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
        value_type result=0;

        if( option(_name="crb.stock-matrices").template as<bool>() )
        {
            result = inner_product( *M_Fqm[l][q][m] , xi );
        }
        else
        {
            auto vector_composite = functionalCompositeF();
            auto composite = vector_composite[l];
            auto functional = composite->functionallinear(q,m);
            auto space = functional->space();
            auto vector = M_backend->newVector( space );
            functional->containerPtr( vector );
            result = inner_product( *vector, xi ) ;
        }

        return result;
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
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_type const& X, vector_type const& Y )
    {
        return M_model->scalarProductForMassMatrix( X, Y );
    }
    /**
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        return M_model->scalarProductForMassMatrix( X, Y );
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
        //return this->solveFemUsingAffineDecompositionFixedPoint( mu );
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

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type &u , bool need_to_solve=false)
    {
        return M_model->output( output_index, mu , u , need_to_solve );
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

    vectorN_type computeStatistics ( Eigen::VectorXd vector , std::string name )
    {
        return M_model->computeStatistics( vector , name );
    }

    void readConvergenceDataFromFile( std::vector< vectorN_type > & vector, std::string filename )
    {
        return M_model->readConvergenceDataFromFile( vector, filename );
    }

    void generateGeoFileForOutputPlot(  vectorN_type outputs, vectorN_type parameter, vectorN_type estimated_error )
    {
        return M_model->generateGeoFileForOutputPlot( outputs, parameter, estimated_error );
    }

    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector, std::string filename )
    {
        return M_model->writeConvergenceStatistics( vector, filename );
    }

    void writeVectorsExtremumsRatio(std::vector< vectorN_type > const& vector1, std::vector< vectorN_type > const& vector2, std::string filename )
    {
        return M_model->writeVectorsExtremumsRatio( vector1, vector2, filename );
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

    mutable std::vector< std::vector<element_ptrtype> > M_InitialGuessV;
    mutable std::vector< std::vector<vector_ptrtype> > M_InitialGuessVector;

    //! affine decomposition terms ( time dependent )
    mutable std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;

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
    backend_ptrtype M_backend_primal;
    backend_ptrtype M_backend_dual;

    // ! matrix associated with inner product
    sparse_matrix_ptrtype M_B;
    sparse_matrix_ptrtype M_H1;

    beta_vector_type M_dummy_betaMqm;

    //! initialize the matrix associated with the \f$H_1\f$ inner product
    void initB();

    /**
     * \brief given \p mu merge the Aq and Fq into A and F respectively
     *
     * \param mu the parameter at which the matrix A and vector F are assembled
     *
     */
    offline_merge_type offlineMerge( betaqm_type const& all_beta, parameter_type const& mu );
    offline_merge_type offlineMergeOnFly( betaqm_type const& all_beta , parameter_type const& mu  );

    void assembleMassMatrix( );
    void assembleMassMatrix( mpl::bool_<true> );
    void assembleMassMatrix( mpl::bool_<false> );

    operatorcomposite_ptrtype preAssembleMassMatrix( ) const ;
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<true> ) const ;
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<false> ) const ;

    void assembleInitialGuessV( initial_guess_type & initial_guess );
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<true> );
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<false> );
    element_type u,v;

    preconditioner_ptrtype M_preconditioner_primal;
    preconditioner_ptrtype M_preconditioner_dual;

};


template <typename ModelType>
struct PreAssembleMassMatrixInCompositeCase
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


    typedef OperatorLinearComposite<space_type, space_type> operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    PreAssembleMassMatrixInCompositeCase( element_type const u ,
                                          element_type const v )
        :
        M_composite_u ( u ),
        M_composite_v ( v )
    {
        auto Xh = M_composite_u.functionSpace();
        op_mass = opLinearComposite( _imageSpace=Xh, _domainSpace=Xh );
    }

    template< typename T >
    void
    operator()( const T& t ) const
    {

        using namespace Feel::vf;

        auto u = M_composite_u.template element< T::value >();
        auto v = M_composite_v.template element< T::value >();
        auto Xh = M_composite_u.functionSpace();
        mesh_ptrtype mesh = Xh->mesh();

        auto expr = integrate( _range=elements( mesh ) , _expr=trans( idt( u ) )*id( v ) ) ;
        auto opfree = opLinearFree( _imageSpace=Xh, _domainSpace=Xh, _expr=expr );

        //each composant of the affine decomposition
        //is the subspace contribution
        int q=T::value;
        int m=0;
        auto tuple = boost::make_tuple( q , m );
        op_mass->addElement( tuple , opfree );
    }


    operatorcomposite_ptrtype opmass()
    {
        return op_mass;
    }

    element_type  M_composite_u;
    element_type  M_composite_v;
    operatorcomposite_ptrtype op_mass;

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
                                       boost::shared_ptr<CRBModel<ModelType> > crb_model)
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

        form2( _test=Xh, _trial=Xh, _matrix=M_crb_model->Mqm(0,0) ) +=
            integrate( _range=elements( mesh ), _expr=trans( idt( u ) )*id( v ) );

    }

    element_type  M_composite_u;
    element_type  M_composite_v;
    mutable boost::shared_ptr<CRBModel<ModelType>  > M_crb_model;
};

template <typename ModelType>
struct AssembleInitialGuessVInCompositeCase
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

    AssembleInitialGuessVInCompositeCase( element_type  const v ,
                                          initial_guess_type  const initial_guess ,
                                          boost::shared_ptr<CRBModel<ModelType> > crb_model)
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
        int q_max = M_crb_model->QInitialGuess();
        for(int q = 0; q < q_max; q++)
        {
            int m_max = M_crb_model->mMaxInitialGuess(q);
            for( int m = 0; m < m_max ; m++)
            {
                auto initial_guess_qm = M_crb_model->InitialGuessVector(q,m);
                auto view = M_composite_initial_guess[q][m]->template element< T::value >();
                form1( _test=Xh, _vector=initial_guess_qm ) +=
                    integrate ( _range=elements( mesh ), _expr=trans( Feel::vf::idv( view ) )*Feel::vf::id( v ) );
            }
        }

    }

    element_type  M_composite_v;
    initial_guess_type  M_composite_initial_guess;
    mutable boost::shared_ptr<CRBModel<ModelType> > M_crb_model;
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
              _solver=( EigenSolverType )M_vm["solvereigen.solver-type"].as<int>(),
              //_spectrum=LARGEST_MAGNITUDE,
              _spectrum=SMALLEST_MAGNITUDE,
              //_transform=SINVERT,
              _ncv=M_vm["solvereigen.ncv"].as<int>(),
              _nev=M_vm["solvereigen.nev"].as<int>(),
              _tolerance=M_vm["solvereigen.tol"].as<double>(),
              _maxit=M_vm["solvereigen.maxiter"].as<int>()
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


//create a vector of preassemble objects
//for the mass matrix

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix() const
{
    static const bool is_composite = functionspace_type::is_composite;
    return preAssembleMassMatrix( mpl::bool_< is_composite >() );
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix( mpl::bool_<false> ) const
{

    auto Xh = M_model->functionSpace();
    auto mesh = Xh->mesh();

    auto expr=integrate( _range=elements( mesh ) , _expr=idt( u )*id( v ) );
    auto op_mass = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh  );
    auto opfree = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr );
    opfree->setName("mass operator (automatically created)");
    //in this case, the affine decompositon has only one element
    op_mass->addElement( boost::make_tuple(0,0) , opfree );
    return op_mass;
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix( mpl::bool_<true> ) const
{
    auto Xh = M_model->functionSpace();

    index_vector_type index_vector;
    PreAssembleMassMatrixInCompositeCase<TruthModelType> preassemble_mass_matrix_in_composite_case ( u , v );
    fusion::for_each( index_vector, preassemble_mass_matrix_in_composite_case );

    auto op_mass = preassemble_mass_matrix_in_composite_case.opmass();
    return op_mass;
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
    form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=idt( u )*id( v )  );
    M_Mqm[0][0]->close();
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleMassMatrix( mpl::bool_<true> )
{
    auto Xh = M_model->functionSpace();

    index_vector_type index_vector;
    int size = functionspace_type::nSpaces;
    M_Mqm.resize( 1 );
    M_Mqm[0].resize(1);
    M_Mqm[0][0]=M_backend->newMatrix( _test=Xh , _trial=Xh );

    AssembleMassMatrixInCompositeCase<TruthModelType> assemble_mass_matrix_in_composite_case ( u , v , this->shared_from_this());
    fusion::for_each( index_vector, assemble_mass_matrix_in_composite_case );

    M_Mqm[0][0]->close();
}



template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleInitialGuessV( initial_guess_type & initial_guess )
{
    static const bool is_composite = functionspace_type::is_composite;
    return assembleInitialGuessV( initial_guess, mpl::bool_< is_composite >() );
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<true> )
{
    auto Xh = M_model->functionSpace();
    auto mesh = Xh->mesh();

    int q_max= this->QInitialGuess();
    M_InitialGuessV.resize( q_max );
    M_InitialGuessVector.resize( q_max );
    for(int q = 0; q < q_max; q++ )
    {
        int m_max= this->mMaxInitialGuess(q);
        M_InitialGuessV[q].resize( m_max );
        M_InitialGuessVector[q].resize( m_max );
        for(int m = 0; m < m_max; m++ )
        {
            M_InitialGuessV[q][m] = Xh->elementPtr();
            M_InitialGuessVector[q][m] = M_model->newVector();
        }
    }

    index_vector_type index_vector;
    AssembleInitialGuessVInCompositeCase<TruthModelType> assemble_initial_guess_v_in_composite_case ( v , initial_guess , this->shared_from_this());
    fusion::for_each( index_vector, assemble_initial_guess_v_in_composite_case );

    for(int q = 0; q < q_max; q++ )
    {
        int m_max = this->mMaxInitialGuess(q) ;
        for(int m = 0; m < m_max; m++ )
            M_InitialGuessVector[q][m]->close();
    }

}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<false> )
{
    using namespace Feel::vf;
    auto Xh = M_model->functionSpace();
    auto mesh = Xh->mesh();

    int q_max= this->QInitialGuess();
    M_InitialGuessV.resize( q_max );
    M_InitialGuessVector.resize( q_max );
    for(int q = 0; q < q_max; q++ )
    {
        int m_max= this->mMaxInitialGuess(q);
        M_InitialGuessV[q].resize( m_max );
        M_InitialGuessVector[q].resize( m_max );
        for(int m = 0; m < m_max; m++ )
        {
            M_InitialGuessV[q][m] = Xh->elementPtr();
            M_InitialGuessVector[q][m] = M_model->newVector();
            form1( _test=Xh, _vector=M_InitialGuessVector[q][m]) =
                integrate( _range=elements( mesh ), _expr=idv( initial_guess[q][m] )*id( v )  );
            M_InitialGuessVector[q][m]->close();
        }
    }
}


template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_ptrtype
CRBModel<TruthModelType>::assembleInitialGuess( parameter_type const& mu )
{
    auto Xh = M_model->functionSpace();
    element_ptrtype initial_guess = Xh->elementPtr();
    initial_guess_type vector_initial_guess;
    beta_vector_type beta;
    vector_initial_guess = M_model->computeInitialGuessAffineDecomposition();
    beta = M_model->computeBetaInitialGuess( mu );

    int q_max = vector_initial_guess.size();
    for ( size_type q = 0; q < q_max; ++q )
    {
        int m_max=vector_initial_guess[q].size();
        for ( size_type m = 0; m < m_max ; ++m )
        {
            element_type temp = Xh->element();
            temp = *vector_initial_guess[q][m];
            temp.scale( beta[q][m] );
            *initial_guess += temp;
        }
    }
    return initial_guess;
}


template<typename TruthModelType>
typename CRBModel<TruthModelType>::offline_merge_type
CRBModel<TruthModelType>::offlineMergeOnFly(betaqm_type const& all_beta, parameter_type const& mu  )
{
    //recovery from the model of free operators Aqm in A = \sum_q \sum_m \beta_qm( u ; mu ) A_qm( u, v )
    auto compositeA = this->operatorCompositeA();
    //idem with other operator / functional
    auto compositeM = this->operatorCompositeM();
    auto vector_compositeF = this->functionalCompositeF();

    //acces to beta coefficients
    auto beta_M = all_beta.template get<0>();
    auto beta_A = all_beta.template get<1>();
    //warning : beta_F is a vector of beta_coefficients
    auto beta_F = all_beta.template get<2>();

    //associate beta coefficients to operators
    compositeA->setScalars( beta_A );
    compositeM->setScalars( beta_M );

    //merge
    auto A = M_model->newMatrix();
    auto M = M_model->newMatrix();
    compositeA->sumAllMatrices( A );
    //auto A = compositeA->sumAllMatrices();
    //auto M = compositeM->sumAllMatrices();
    compositeM->sumAllMatrices( M );

    std::vector<vector_ptrtype> F( Nl() );

    for(int output=0; output<Nl(); output++)
    {
        auto compositeF = vector_compositeF[output];
        compositeF->setScalars( beta_F[output] );
        F[output] = M_model->newVector();
        compositeF->sumAllVectors( F[output] );
    }

    return boost::make_tuple( M, A, F );
}


template<typename TruthModelType>
typename CRBModel<TruthModelType>::offline_merge_type
CRBModel<TruthModelType>::offlineMerge( betaqm_type const& all_beta , parameter_type const& mu )
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

#else

    auto A = this->newMatrix();
    auto M = this->newMatrix();
    //size_type pattern = operatorCompositeA()->operatorlinear(0,0)->pattern();

#endif
    std::vector<vector_ptrtype> F( Nl() );

    //acces to beta coefficients
    auto beta_M = all_beta.template get<0>();
    auto beta_A = all_beta.template get<1>();
    //warning : beta_F is a vector of beta_coefficients
    auto beta_F = all_beta.template get<2>();

    A->zero();
    for ( size_type q = 0; q < Qa(); ++q )
    {
        for(size_type m = 0; m < mMaxA(q); ++m )
            A->addMatrix( beta_A[q][m], M_Aqm[q][m] );
    }

    if( Qm() > 0 )
    {
        for ( size_type q = 0; q < Qm(); ++q )
        {
            for(size_type m = 0; m < mMaxM(q) ; ++m )
                M->addMatrix( beta_M[q][m] , M_Mqm[q][m] );
        }
    }

    for ( size_type l = 0; l < Nl(); ++l )
    {
        F[l] = M_backend->newVector( M_model->functionSpace() );
        F[l]->zero();

        for ( size_type q = 0; q < Ql( l ); ++q )
        {
            for ( size_type m = 0; m < mMaxF(l,q); ++m )
                F[l]->add( beta_F[l][q][m] , M_Fqm[l][q][m] );
        }
        F[l]->close();
    }

    return boost::make_tuple( M, A, F );
}


template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemUsingOfflineEim( parameter_type const& mu )
{
    auto Xh= this->functionSpace();

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype A;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    element_ptrtype InitialGuess = Xh->elementPtr();
    vector_ptrtype Rhs( M_backend->newVector( Xh ) );
    auto u = Xh->element();

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
    }
    else
    {
        time_initial=this->timeInitial();
        time_step=this->timeStep();
        time_final=this->timeFinal();
        this->initializationField( InitialGuess, mu );
    }

    mybdf->setTimeInitial( time_initial );
    mybdf->setTimeStep( time_step );
    mybdf->setTimeFinal( time_final );

    double bdf_coeff ;
    auto vec_bdf_poly = M_backend->newVector( Xh );

    for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
    {
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        boost::tie(M, A, F) = this->update( mu , mybdf->time() );
        *Rhs = *F[0];
        if( !isSteady() )
        {
            A->addMatrix( bdf_coeff, M );
            Rhs->addVector( *vec_bdf_poly, *M );
        }
        M_preconditioner_primal->setMatrix( A );
        M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs , _prec=M_preconditioner_primal);
        mybdf->shiftRight(u);
    }

    return u;

}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu )
{
    auto Xh= this->functionSpace();

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype A;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    element_ptrtype InitialGuess = Xh->elementPtr();
    auto u = Xh->element();
    auto uold = Xh->element();
    vector_ptrtype Rhs( M_backend->newVector( Xh ) );

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
        //we want to have the initial guess given by function update
        InitialGuess = this->assembleInitialGuess( mu ) ;
    }
    else
    {
        time_initial=this->timeInitial();
        time_step=this->timeStep();
        time_final=this->timeFinal();
        this->initializationField( InitialGuess, mu );
    }

    mybdf->setTimeInitial( time_initial );
    mybdf->setTimeStep( time_step );
    mybdf->setTimeFinal( time_final );

    u=*InitialGuess;
    double norm=0;
    int iter=0;

    double bdf_coeff ;
    auto vec_bdf_poly = M_backend->newVector( Xh );

    int max_fixedpoint_iterations  = this->vm()["crb.max-fixedpoint-iterations"].template as<int>();
    double increment_fixedpoint_tol  = this->vm()["crb.increment-fixedpoint-tol"].template as<double>();

    for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
    {
        iter=0;
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        do {
            boost::tie(M, A, F) = this->update( mu , u , mybdf->time() );
            *Rhs = *F[0];

            if( !isSteady() )
            {
                A->addMatrix( bdf_coeff, M );
                Rhs->addVector( *vec_bdf_poly, *M );
            }
            uold = u;
            M_preconditioner_primal->setMatrix( A );
            M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs , _prec=M_preconditioner_primal);
            if( option(_name="crb.use-linear-model").template as<bool>() )
                norm = 0;
            else
                norm = this->computeNormL2( uold , u );
            iter++;
        } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
        mybdf->shiftRight(u);
    }
    return u;
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu )
{
    int output_index = option(_name="crb.output-index").template as<int>();

    auto Xh= this->functionSpace();

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype A,Adu;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    auto udu = Xh->element();
    auto uold = Xh->element();
    vector_ptrtype Rhs( M_backend->newVector( Xh ) );
    auto dual_initial_field = Xh->elementPtr();

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
        //InitialGuess = this->assembleInitialGuess( mu ) ;
    }
    else
    {
        time_initial=this->timeFinal()+this->timeStep();
        time_step=-this->timeStep();
        time_final=this->timeInitial()+this->timeStep();
    }

    mybdf->setTimeInitial( time_initial );
    mybdf->setTimeStep( time_step );
    mybdf->setTimeFinal( time_final );

    double norm=0;
    int iter=0;

    double bdf_coeff ;
    auto vec_bdf_poly = M_backend->newVector( Xh );

    if ( this->isSteady() )
        udu.zero() ;
    else
    {
        boost::tie( M, A, F) = this->update( mu , mybdf->timeInitial() );
        *Rhs=*F[output_index];
        M_preconditioner_dual->setMatrix( M );
        M_backend_dual->solve( _matrix=M, _solution=dual_initial_field, _rhs=Rhs, _prec=M_preconditioner_dual );
        udu=*dual_initial_field;
    }


    int max_fixedpoint_iterations  = option(_name="crb.max-fixedpoint-iterations").template as<int>();
    double increment_fixedpoint_tol  = option(_name="crb.increment-fixedpoint-tol").template as<double>();
    for( mybdf->start(udu); !mybdf->isFinished(); mybdf->next() )
    {
        iter=0;
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        do {
            boost::tie(M, A, F) = this->update( mu , udu , mybdf->time() );

            if( ! isSteady() )
            {
                A->addMatrix( bdf_coeff, M );
                Rhs->zero();
                *vec_bdf_poly = bdf_poly;
                Rhs->addVector( *vec_bdf_poly, *M );
            }
            else
            {
                *Rhs = *F[output_index];
                Rhs->close();
                Rhs->scale( -1 );
            }

            if( option("crb.use-symmetric-matrix").template as<bool>() )
                Adu = A;
            else
                A->transpose( Adu );

            uold = udu;
            M_preconditioner_dual->setMatrix( Adu );
            M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs , _prec=M_preconditioner_dual);

            if( option(_name="crb.use-linear-model").template as<bool>() )
                norm = 0;
            else
                norm = this->computeNormL2( uold , udu );
            iter++;
        } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
        mybdf->shiftRight(udu);
    }
    return udu;
}


template<typename TruthModelType>
void
CRBModel<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    M_model->run( X, N, Y, P );
}



}
#endif /* __CRBModel_H */
