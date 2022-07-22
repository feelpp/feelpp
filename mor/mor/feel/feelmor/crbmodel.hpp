/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-09

  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)

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
#ifndef FEELPP_MOR_CRBModel_H
#define FEELPP_MOR_CRBModel_H

//#include <boost/shared_ptr.hpp>

#include <vector>


#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelts/bdf.hpp>
//#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/operatorlinearfree.hpp>
#include <feel/feeldiscr/operatorlinearcomposite.hpp>
#include <feel/feeldiscr/fsfunctionallinearfree.hpp>
#include <feel/feeldiscr/fsfunctionallinearcomposite.hpp>

#include <feel/feelmor/parameterspace.hpp>
//#include <feel/feelcore/pslogger.hpp>
#include <feel/feelalg/aitken.hpp>

//#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelmor/crbmodeldb.hpp>
#include <feel/feelmor/crbmodelbase.hpp>

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
class CRBModel : public CRBModelBase,
                 public std::enable_shared_from_this<CRBModel<ModelType> >
{
public:


    /** @name Constants
     */
    //@{

    static const bool is_time_dependent = ModelType::is_time_dependent;
    static const bool is_linear = ModelType::is_linear;

    //@}

    /** @name Typedefs
     */
    //@{

    //! model type
    typedef ModelType model_type;
    typedef std::shared_ptr<ModelType> model_ptrtype;

    //! value_type
    typedef typename model_type::value_type value_type;
    //! mesh type
    typedef typename ModelType::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;

    //! function space type
    typedef typename model_type::space_type functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    //! reduced basis function space type
    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename model_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;


    //! element of the functionspace type
    typedef typename model_type::space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename model_type::backend_type backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::vector_type vector_type;

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef std::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename model_type::parameter_type parameter_type;
    typedef std::shared_ptr<parameter_type> parameter_ptrtype;


    typedef typename model_type::parameterspace_type::sampling_type sampling_type;
    typedef typename model_type::parameterspace_type::sampling_ptrtype sampling_ptrtype;


    typedef typename std::vector< std::vector < element_ptrtype > > initial_guess_type;
    //typedef Eigen::VectorXd theta_vector_type;
    typedef Eigen::VectorXd vectorN_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype,
                                  sparse_matrix_ptrtype,
                                  std::vector<vector_ptrtype>
                                  > offline_merge_type;


    typedef typename boost::tuple<std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector< std::vector<vector_ptrtype> > >
                                  > affine_decomposition_type;

    typedef typename boost::tuple< sparse_matrix_ptrtype,
                                   sparse_matrix_ptrtype,
                                   std::vector<vector_ptrtype>
                                   > monolithic_type;

    typedef typename boost::tuple< std::map<int,double>,
                                   std::map<int,double>,
                                   std::vector< std::map<int,double> >
                                   > eim_interpolation_error_type ;



    typedef typename model_type::affine_decomposition_light_type affine_decomposition_light_type;
    typedef typename model_type::betaq_type betaq_type;
    typedef typename model_type::beta_vector_type beta_vector_type;
    typedef typename model_type::beta_vector_light_type beta_vector_light_type;

    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > betaqm_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >,
                               boost::tuple<
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >
                               >::type incomplete_betaqm_type;

    //! time discretization
    typedef Bdf<space_type>  bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef std::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef std::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef Preconditioner<double> preconditioner_type;
    typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;

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

    CRBModel( std::string const& name, crb::stage stage, int level = 0 )
        :
        CRBModel( name, Environment::randomUUID( true ), stage, level )
    {}
    CRBModel( std::string const& name, uuids::uuid const& uid, crb::stage stage, int level = 0 )
        :
        CRBModel( std::make_shared<CRBModelDB>(name,uid),
                  std::make_shared<model_type>(), stage, level )
        {}
    CRBModel( std::string const& name, model_ptrtype const& model, crb::stage stage, int level = 0 )
        :
        CRBModel( std::make_shared<CRBModelDB>(name,Environment::randomUUID( true )), model, stage, level )
        {}

    CRBModel( std::shared_ptr<CRBModelDB> crbModelDb, model_ptrtype const& model, crb::stage stage, int level = 0 )
        :
        M_crbModelDb( crbModelDb ),
        M_level( level ),
        M_Aqm(),
        M_InitialGuessV(),
        M_InitialGuessVector(),
        M_Mqm(),
        M_Fqm(),
        M_model( model ),
        M_prefix( model->prefix() ),
        M_is_initialized( false ),
        M_backend( (stage==crb::stage::offline)?backend():nullptr ),
        M_backend_primal( (stage==crb::stage::offline)?backend( _name="backend-primal"):nullptr ),
        M_backend_dual( (stage==crb::stage::offline)?backend( _name="backend-dual"):nullptr ),
        M_backend_l2( (stage==crb::stage::offline)?backend( _name="backend-l2"):nullptr ),
        M_fixedpointUseAitken( boption(_prefix=M_prefix,_name="crb.fixedpoint.aitken") ),
        M_alreadyCountAffineDecompositionTerms( false ),
        M_isSteadyModel( !model_type::is_time_dependent || boption(_prefix=M_prefix,_name="crb.is-model-executed-in-steady-mode") ),
        M_numberOfTimeStep( 1 ),
        M_has_eim( false ),
        M_useSER( ioption(_prefix=M_prefix,_name="ser.rb-frequency") || ioption(_prefix=M_prefix,_name="ser.eim-frequency") ),
        M_usePrimalPc(boption(_prefix=M_prefix,_name="crb.use-primal-pc")),
        M_useSymmMat(boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix")),
        M_stockMatrices(boption(_prefix=M_prefix,_name="crb.stock-matrices")),
        M_serErrorEstimation(boption(_prefix=M_prefix,_name="ser.error-estimation")),
        M_crbUseNewton(boption(_prefix=M_prefix,_name="crb.use-newton")),
        M_fixedpointMaxIt(ioption(_prefix=M_prefix,_name="crb.fixedpoint.maxit")),
        M_fixedpointIncrTol(doption(_prefix=M_prefix,_name="crb.fixedpoint.increment-tol")),
        M_fixedpointVerbose(boption(_prefix=M_prefix,_name="crb.fixedpoint.verbose")),
        M_outputIndex(ioption(_prefix=M_prefix,_name="crb.output-index")),
        M_useLinearModel(boption(_prefix=M_prefix,_name="crb.use-linear-model"))
        {

            M_model->attach( M_crbModelDb );
            bool M_rebuildDb = boption(_prefix=M_prefix,_name="crb.rebuild-database");
            int M_dbLoad = ioption(_prefix=M_prefix, _name="crb.db.load" );
            std::string M_dbFilename = soption(_prefix=M_prefix, _name="crb.db.filename");
            std::string M_dbId = soption(_prefix=M_prefix,_name="crb.db.id");
            int M_dbUpdate = ioption(_prefix=M_prefix, _name="crb.db.update" );
            if ( !M_rebuildDb )
            {
                switch ( M_dbLoad )
                {
                case 0 :
                    M_crbModelDb->updateIdFromDBFilename( M_dbFilename );
                    break;
                case 1:
                    M_crbModelDb->updateIdFromDBLast( crb::last::created );
                    break;
                case 2:
                    M_crbModelDb->updateIdFromDBLast( crb::last::modified );
                    break;
                case 3:
                    M_crbModelDb->updateIdFromId( M_dbId );
                    break;
                }
            }
            else
            {
                switch ( M_dbUpdate )
                {
                case 0 :
                    M_crbModelDb->updateIdFromDBFilename( M_dbFilename );
                    break;
                case 1:
                    M_crbModelDb->updateIdFromDBLast( crb::last::created );
                    break;
                case 2:
                    M_crbModelDb->updateIdFromDBLast( crb::last::modified );
                    break;
                case 3:
                    M_crbModelDb->updateIdFromId( M_dbId );
                    break;
                default:
                    // don't do anything and let the system pick up a new unique id
                    break;
                }
            }

            if ( stage == crb::stage::offline )
                this->init();
        }


    FEELPP_DEPRECATED CRBModel( bool doInit = true )
        :
        CRBModel( doInit?crb::stage::offline:crb::stage::online, 0 )
        {}

    FEELPP_DEPRECATED CRBModel( int level, bool doInit = true )
        :
        CRBModel( doInit?crb::stage::offline:crb::stage::online, level )
    {
    }

    FEELPP_DEPRECATED CRBModel( CRBModelMode mode, int level=0, bool doInit = true )
        :
        CRBModel( doInit?crb::stage::offline:crb::stage::online, level )
        {}

    FEELPP_DEPRECATED CRBModel( model_ptrtype const& model , bool doInit = true )
        :
        CRBModel( model, doInit?crb::stage::offline:crb::stage::online, 0 )
        {
        }

    FEELPP_DEPRECATED CRBModel( model_ptrtype const& model , CRBModelMode mode, bool doInit = true )
        :
        CRBModel( model, doInit?crb::stage::offline:crb::stage::online, 0 )
        {}

    /**
     * copy constructor
     */
    CRBModel( CRBModel const & o )
        :
        M_level( o.M_level ),
        M_Aqm( o.M_Aqm ),
        M_InitialGuessV( o.M_InitialGuessV ),
        M_InitialGuessVector( o.M_InitialGuessVector ),
        M_Mqm( o.M_Mqm ),
        M_Fqm( o.M_Fqm ),
        M_model(  o.M_model ),
        M_is_initialized( o.M_is_initialized ),
        M_backend( o.M_backend ),
        M_backend_primal( o.M_backend_primal ),
        M_backend_dual( o.M_backend_dual ),
        M_backend_l2( o.M_backend_l2 ),
        M_alreadyCountAffineDecompositionTerms( o.M_alreadyCountAffineDecompositionTerms ),
        M_isSteadyModel( o.M_isSteadyModel ),
        M_numberOfTimeStep( o.M_numberOfTimeStep ),
        M_useSER( o.M_useSER )
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

        M_has_eim=false;

        if( M_usePrimalPc )
        {
            M_preconditioner_primal = preconditioner(_pc=(PreconditionerType) M_backend_primal->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                                     _backend= M_backend_primal,
                                                     _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_primal->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                                     _worldcomm=M_backend_primal->comm(),
                                                     _prefix=M_backend_primal->prefix() ,
                                                     _rebuild=M_useSER);
        }

        M_preconditioner_dual = preconditioner(_pc=(PreconditionerType) M_backend_dual->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                               _backend= M_backend_dual,
                                               _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_dual->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                               _worldcomm=M_backend_dual->comm(),
                                               _prefix=M_backend_dual->prefix() ,
                                               _rebuild=M_useSER);
        M_preconditioner_l2 = preconditioner(_pc=(PreconditionerType) M_backend_l2->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                             _backend= M_backend_l2,
                                             _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_l2->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                             _worldcomm=M_backend_l2->comm(),
                                             _prefix=M_backend_l2->prefix() ,
                                             _rebuild=M_useSER);
        M_is_initialized=true;

        if( ! M_model->isInitialized() )
        {
            LOG( INFO ) << "CRBModel Model is not initialized";
            M_model->initModel();
            M_model->setInitialized( true );
        }

        M_model->buildGinacBetaExpressions( M_model->parameterSpace()->min() );

        auto Xh = M_model->functionSpace();
        M_u = Xh->element();
        M_v = Xh->element();

        bool symmetric = M_useSymmMat;
        bool stock = M_stockMatrices;

        // Access to eim eventually built in initModel
        M_has_eim = this->scalarContinuousEim().size() > 0 || this->scalarDiscontinuousEim().size() > 0 || M_model->hasDeim();

        if( stock )
            this->computeAffineDecomposition();
        else
            this->countAffineDecompositionTerms(); //already called in computeAffineDecomposition() if stock

        if( this->hasEim() || (!symmetric) )
        {
            CHECK( stock )<<"There is some work to do before using operators free when using EIM, for now we compute (and stock matrices) affine decomposition to assemble the inner product \n";

            //in this case, we use linear part of bilinear form a
            //as the inner product
            auto muref = this->refParameter();
            M_inner_product_matrix = this->newMatrix();
            M_inner_product_matrix->zero();
            if ( M_linearAqm.size() )
            {
                auto betaqm = computeBetaLinearDecompositionA( muref );
                for ( size_type q = 0; q < M_QLinearDecompositionA; ++q )
                {
                    for(size_type m = 0; m < mMaxLinearDecompositionA(q); ++m )
                        M_inner_product_matrix->addMatrix( betaqm[q][m], M_linearAqm[q][m] );
                }
            }

            //check that the matrix is filled, else we take energy matrix
            double norm=M_inner_product_matrix->l1Norm();
            if( norm == 0 && symmetric )
            {
                M_inner_product_matrix = M_model->energyMatrix();
            }
        }
        else
        {
            //in this case, we use bilinear form a
            //as the inner product
            M_inner_product_matrix = M_model->energyMatrix();
            CHECK( symmetric )<< "You use energy matrix as inner product but you specified that bilinear form a() is not symmetric !\n";
        }
        M_preconditioner_l2->setMatrix( M_inner_product_matrix );


        if ( this->isSteady() )
            M_numberOfTimeStep=1;
        else
        {
            M_numberOfTimeStep=0;
            for ( double t=timeInitial();t<=(timeFinal()+1e-9);t+=timeStep() )
                ++M_numberOfTimeStep;

            // TODO VINCENT : maybe not the good place and must be rebuild at each rbspace update
            M_onlineTimeStepping = bdf(_space=this->rBFunctionSpace(),
                                       _initial_time=this->timeInitial(),_final_time=this->timeFinal(),
                                       _time_step=this->timeStep(),
                                       _save=false);
        }



        if ( countoption( _name="crb.copy-files-inside-db.path",_prefix=this->prefix() ) > 0 )
        {
            int cpt = 0;
            for ( std::string const& filepathstr : vsoption( _name="crb.copy-files-inside-db.path",_prefix=this->prefix() ) )
            {
                auto filepath = fs::path( Environment::expand(filepathstr) );
                M_model->addModelData( fmt::format("copy_files_inside_db_{}",cpt++), filepath );
            }
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
            M_backend_primal = o.M_backend_primal;
            M_backend_dual = o.M_backend_dual;
            M_backend_l2 = o.M_backend_l2;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    std::string const& prefix() const { return M_prefix; }

    virtual bool useMonolithicRbSpace() { return !model_type::by_block; }

    //!
    //! world communicator
    //!
    WorldComm const& worldComm() const { return M_model->worldComm(); }

    //!
    //! world communicator
    //!
    worldcomm_ptr_t const& worldCommPtr() const { return M_model->worldCommPtr(); }

    /**
     * \return  the \p variables_map
     */
    po::variables_map const& vm() const
    {
        return Environment::vm();
    }

    /**
     * \return model
     */
    model_ptrtype const& model() const { return M_model; }

    /**
     * \return model
     */
    model_ptrtype & model() { return M_model; }

    //!
    //! get the id of the model
    //!
    uuids::uuid uuid() const { return M_model->uuid(); }

    //!
    //! in case of hierarchy of models, return level index.
    //! default value is 0
    //!
    int level() const { return M_level; }

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const
    {
        auto Xh = M_model->functionSpace();
        return M_backend->newMatrix( _test=Xh, _trial=Xh );
    }

    /**
     * create a new vector
     * \return the newly created vector
     */
    vector_ptrtype newVector() const
    {
        auto Xh = M_model->functionSpace();
        return M_backend->newVector( Xh );
    }

    /**
     * \brief Returns the matrix associated with the inner product
     * linked to energy norm
     */
    sparse_matrix_ptrtype const& energyMatrix() const
    {
        //return M_model->innerProduct();
        return M_inner_product_matrix;
    }

    /**
     * \brief Returns the matrix associated with the inner product
     * linked to energy norm
     */
    sparse_matrix_ptrtype  energyMatrix()
    {
        //return M_model->innerProduct();
        return M_inner_product_matrix;
    }

    /**
     * \brief Returns the matrix associated with the inner product
     */
    sparse_matrix_ptrtype const& massMatrix() const
    {
        return M_model->massMatrix();
    }

    /**
     * \brief Returns the matrix associated with the inner product
     */
    sparse_matrix_ptrtype  massMatrix()
    {
        return M_model->massMatrix();
    }

    /**
     * \brief Returns the matrix associated with the inner product
     * used to perform the POD (parabolic case)
     */
    sparse_matrix_ptrtype const& innerProductForPod() const
    {
        //return M_model->innerProductForPod();
        return M_inner_product_matrix;
    }

    /**
     * \brief Returns the matrix associated with the inner product
     * used to perform the POD (parabolic case)
     */
    sparse_matrix_ptrtype  innerProductForPod()
    {
        return M_inner_product_matrix;
        //return M_model->innerProductForPod();
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
    virtual size_type Qa() const
    {
        return M_Qa;
    }

    int QLinearDecompositionA() const
    {
        return M_QLinearDecompositionA;
    }

    //! return the number of \f$\mu\f$ independent terms for the bilinear form ( time dependent )
    //int Qm() const { return 1; }//return M_model->Qm(); }

    size_type Qm() const
    {
        return Qm( mpl::bool_<model_type::is_time_dependent>() );
    }
    size_type Qm( mpl::bool_<true> ) const
    {
        return M_Qm;
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
        return QInitialGuess( mpl::bool_<model_type::is_linear>() );
    }
    int QInitialGuess( mpl::bool_<true> ) const
    {
        return 0;
    }
    int QInitialGuess( mpl::bool_<false> ) const
    {
        return M_model->QInitialGuess();
    }


    int mMaxA(int q )
    {
        int size=M_mMaxA.size();
        CHECK( q < size ) << "mMaxA called with the bad q index "<<q<<" and max is "<<size<<"\n";
        return M_mMaxA[ q ];
    }


    int mMaxM( int q )
    {
        return mMaxM(q, mpl::bool_<model_type::is_time_dependent>() );
    }
    int mMaxM( int q, mpl::bool_<true> )
    {
        int size=M_mMaxM.size();
        CHECK( q < size ) << "mMaxM called with the bad q index "<<q<<" and max is "<<size<<"\n";
        return M_mMaxM[ q ];
    }
    int mMaxM( int q , mpl::bool_<false> )
    {
        return 1;
    }

    int mMaxLinearDecompositionA( int q ) const
    {
        return M_mMaxLinearDecompositionA[q];
    }

    int mMaxInitialGuess( int q ) const
    {
        return M_model->mMaxInitialGuess( q );
    }

    int mMaxF(int output_index, int q )
    {
        bool goodidx=true;
        int size=M_mMaxF.size();
        if( output_index >= size )
            goodidx=false;
        else
        {
            size=M_mMaxF[output_index].size();
            if( q >= size )
                goodidx=false;
        }
        CHECK( goodidx )<<"mMaxF functions called with bad index ! output index : "<<output_index<<" and q : "<<q<<"\n";
        return M_mMaxF[output_index][q];
    }

    //! return the number of outputs
    size_type Nl() const
    {
        return M_Nl;
    }

    //! return the number of \f$\mu\f$ independent terms for the right hand side
    virtual size_type Ql( int l ) const
    {
        return M_Ql[l];
        //return M_model->Ql( l );
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_model->parameterSpace();
    }

    void adaptMesh( parameter_type const& mu )
    {
        return M_model->adaptMesh( mu );
    }

    parameter_type refParameter()
    {
        bool user_specify_parameters = M_model->referenceParametersGivenByUser();
        auto Dmu = M_model->parameterSpace();
        auto muref = Dmu->min();
        if( user_specify_parameters )
            muref = M_model->refParameter();
        return muref;
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
     * \brief compute the beta coefficients
     * mu : parameter
     * time : the time
     * only_time_dependent_terms : return evaluation time-dependent terms only if true, else return evaluation of all terms
     */
    betaqm_type computeBetaQm( parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        return computeBetaQm( mu , mpl::bool_<model_type::is_time_dependent>(), time , only_time_dependent_terms );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<true>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            beta_vector_type,
            std::vector<beta_vector_type> >
            beta_coefficients;

        //check if a light version is available
        auto beta_light = M_model->computeBetaQ( mu, time , only_time_dependent_terms );
        auto betaAlight = beta_light.template get<0>();
        if( betaAlight.size() == 0 )
        {
            beta_coefficients = M_model->computeBetaQm( mu , time , only_time_dependent_terms );
        }
        else
        {
            beta_coefficients = extendBetaCoefficients( beta_light );
        }

        return beta_coefficients;
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            steady_beta;


        //check if a light version is available
        auto beta_light = M_model->computeBetaQ( mu );
        auto betaAlight = beta_light.template get<0>();
        if( betaAlight.size() == 0 )
        {
            steady_beta = M_model->computeBetaQm( mu );
        }
        else
        {
            steady_beta = extendBetaCoefficients( beta_light );
        }

        betaAqm = steady_beta.template get<0>();
        betaFqm = steady_beta.template get<1>();

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

    betaqm_type computeBetaQm( vector_ptrtype const& T, parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        auto solution = M_model->functionSpace()->element();
        solution = *T;
        return computeBetaQm( solution , mu , time , only_time_dependent_terms);
    }
    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        return computeBetaQm( T , mu , mpl::bool_<model_type::is_time_dependent>(), time , only_time_dependent_terms);
    }
    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<true>, double time=0 , bool only_time_dependent_terms=false )
    {
        return M_model->computeBetaQm( T, mu, time , only_time_dependent_terms );
    }
    betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<false>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm, betaMqm ;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            steady_beta;

        steady_beta = M_model->computeBetaQm(T, mu );
        betaAqm = steady_beta.template get<0>();
        betaFqm = steady_beta.template get<1>();

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

    betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        return computeBetaQm( urb , mu , mpl::bool_<model_type::is_time_dependent>(), time , only_time_dependent_terms);
    }
    betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu , mpl::bool_<true>, double time=0 , bool only_time_dependent_terms=false )
    {
        return M_model->computeBetaQm( urb, mu, time , only_time_dependent_terms );
    }
    betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu , mpl::bool_<false>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm, betaMqm ;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            steady_beta;

        steady_beta = M_model->computeBetaQm( urb, mu );
        betaAqm = steady_beta.template get<0>();
        betaFqm = steady_beta.template get<1>();

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


    betaqm_type computePicardBetaQm( parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        return computePicardBetaQm( mu , mpl::bool_<model_type::is_time_dependent>(), time , only_time_dependent_terms );
    }

    betaqm_type computePicardBetaQm( parameter_type const& mu , mpl::bool_<true>, double time=0 , bool only_time_dependent_terms=false )
    {
        boost::tuple<beta_vector_type,beta_vector_type,std::vector<beta_vector_type> >  beta_coefficients;
        beta_coefficients = M_model->computeBetaQm( mu , time , only_time_dependent_terms );
        return beta_coefficients;
    }
    betaqm_type computePicardBetaQm( parameter_type const& mu , mpl::bool_<false>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<beta_vector_type,std::vector<beta_vector_type> > steady_beta;

        steady_beta = M_model->computePicardBetaQm( mu );

        betaAqm = steady_beta.template get<0>();
        betaFqm = steady_beta.template get<1>();

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
    betaqm_type computePicardBetaQm( vector_ptrtype const& T, parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        auto solution = M_model->functionSpace()->element();
        solution = *T;
        return computePicardBetaQm( solution , mu , time , only_time_dependent_terms);
    }
    betaqm_type computePicardBetaQm( element_type const& T, parameter_type const& mu , double time=0 , bool only_time_dependent_terms=false )
    {
        return computePicardBetaQm( T , mu , mpl::bool_<model_type::is_time_dependent>(), time , only_time_dependent_terms);
    }
    betaqm_type computePicardBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<true>, double time=0 , bool only_time_dependent_terms=false )
    {
        return M_model->computePicardBetaQm( T, mu, time , only_time_dependent_terms );
    }
    betaqm_type computePicardBetaQm( element_type const& T, parameter_type const& mu , mpl::bool_<false>, double time=0 , bool only_time_dependent_terms=false )
    {
        beta_vector_type betaAqm, betaMqm ;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<beta_vector_type,std::vector<beta_vector_type> > steady_beta;

        steady_beta = M_model->computePicardBetaQm(T, mu );
        betaAqm = steady_beta.template get<0>();
        betaFqm = steady_beta.template get<1>();

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
    void assemble(){ return M_model->assemble(); }
    /**
     * \brief update the model wrt \p mu
     */
    offline_merge_type update( parameter_type const& mu,  double time=0 , bool only_time_dependent_terms=false )
    {
        auto all_beta = this->computeBetaQm( mu , time , only_time_dependent_terms );
        offline_merge_type offline_merge;
        if( M_stockMatrices )
            offline_merge = offlineMerge( all_beta , only_time_dependent_terms );
        else
            offline_merge = offlineMergeOnFly( all_beta, only_time_dependent_terms );

        return offline_merge;
    }

    offline_merge_type update( parameter_type const& mu, element_type const& T, double time=0 , bool only_time_dependent_terms=false )
    {
#if !defined(NDEBUG)
        mu.check();
#endif
        auto all_beta = this->computeBetaQm(  T , mu , time , only_time_dependent_terms );

        offline_merge_type offline_merge;

        if( M_stockMatrices )
            offline_merge = offlineMerge( all_beta , only_time_dependent_terms );
        else
            offline_merge = offlineMergeOnFly( all_beta, only_time_dependent_terms );

        return offline_merge;

    }
    offline_merge_type update( parameter_type const& mu, vectorN_type const& urb, double time=0 , bool only_time_dependent_terms=false )
    {
#if !defined(NDEBUG)
        mu.check();
#endif
        auto all_beta = this->computeBetaQm( urb, mu , time , only_time_dependent_terms );
        offline_merge_type offline_merge;
        if( M_stockMatrices )
            offline_merge = offlineMerge( all_beta , only_time_dependent_terms );
        else
            offline_merge = offlineMergeOnFly( all_beta, only_time_dependent_terms );

        return offline_merge;

    }

    element_type solveFemMonolithicFormulation( parameter_type const& mu );
    element_type solveFemDualMonolithicFormulation( parameter_type const& mu );
    element_type solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu );
    element_type solveFemUsingAffineDecompositionNewton( parameter_type const& mu );
    void solveFemUpdateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype & J , const parameter_type & mu);
    void solveFemUpdateResidual( const vector_ptrtype& X, vector_ptrtype& R , const parameter_type & mu);
    bool updateJacobian( vector_ptrtype const& X, std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm);
    bool updateJacobian( element_type const& X, std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm);
    bool updateResidual( vector_ptrtype const& X, std::vector< std::vector< std::vector<vector_ptrtype> > >& Rqm);
    bool updateResidual( element_type const& X, std::vector< std::vector< std::vector<vector_ptrtype> > >& Rqm);
    element_type solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu );
    element_type solveFemUsingOfflineEim( parameter_type const& mu );


    /**
     * \brief update model description in property_tree
     * \param ptree to update
     */
    void updatePropertyTree( boost::property_tree::ptree & ptree ) const
    {
        M_model->updatePropertyTree( ptree );

        boost::property_tree::ptree ptree_mMaxA;
        for (int k=0;k<M_mMaxA.size();++k )
        {
            boost::property_tree::ptree ptree_mMaxA_value;
            ptree_mMaxA_value.put( "",M_mMaxA[k] );
            ptree_mMaxA.push_back( std::make_pair("", ptree_mMaxA_value) );
        }
        boost::property_tree::ptree ptree_mMaxM;
        for (int k=0;k<M_mMaxM.size();++k )
        {
            boost::property_tree::ptree ptree_mMaxM_value;
            ptree_mMaxM_value.put( "",M_mMaxM[k] );
            ptree_mMaxM.push_back( std::make_pair("", ptree_mMaxM_value) );
        }
        boost::property_tree::ptree ptree_mMaxF;
        for (int k=0;k<M_mMaxF.size();++k )
        {
            boost::property_tree::ptree ptree_mMaxFsub;
            for (int k2=0;k2<M_mMaxF[k].size();++k2 )
            {
                boost::property_tree::ptree ptree_mMaxFsub_value;
                ptree_mMaxFsub_value.put( "",M_mMaxF[k][k2] );
                ptree_mMaxFsub.push_back( std::make_pair("", ptree_mMaxFsub_value) );
            }
            ptree_mMaxF.push_back( std::make_pair("", ptree_mMaxFsub) );
        }
        boost::property_tree::ptree ptree_mMaxLinearDecompositionA;
        for (int k=0;k<M_mMaxLinearDecompositionA.size();++k )
        {
            boost::property_tree::ptree ptree_mMaxLinearDecompositionA_value;
            ptree_mMaxLinearDecompositionA_value.put( "",M_mMaxLinearDecompositionA[k] );
            ptree_mMaxLinearDecompositionA.push_back( std::make_pair("", ptree_mMaxLinearDecompositionA_value) );
        }


        auto ptreeAffineDecompositionOptional = ptree.get_child_optional( "affine-decomposition" );
        if ( ptreeAffineDecompositionOptional )
        {
            ptreeAffineDecompositionOptional->add_child( "mMaxA", ptree_mMaxA );
            ptreeAffineDecompositionOptional->add_child( "mMaxM", ptree_mMaxM );
            ptreeAffineDecompositionOptional->add_child( "mMaxF", ptree_mMaxF );
            ptreeAffineDecompositionOptional->add_child( "mMaxLinearDecompositionA", ptree_mMaxLinearDecompositionA );
        }
        else
        {
            boost::property_tree::ptree ptreeAffineDecomposition;
            ptreeAffineDecomposition.add_child( "mMaxA", ptree_mMaxA );
            ptreeAffineDecomposition.add_child( "mMaxM", ptree_mMaxM );
            ptreeAffineDecomposition.add_child( "mMaxF", ptree_mMaxF );
            ptreeAffineDecomposition.add_child( "mMaxLinearDecompositionA", ptree_mMaxLinearDecompositionA );

            ptree.add_child( "affine-decomposition", ptreeAffineDecomposition );
        }


        if ( !this->isSteady() )
        {
            nl::json j_timestepping;
            j_timestepping["time_step"] = this->bdfModel()->timeStep();
            j_timestepping["time_initial"] = this->bdfModel()->timeInitial();
            j_timestepping["time_final"] = this->bdfModel()->timeFinal();

            boost::property_tree::ptree pt_timestepping;
            std::istringstream istr( j_timestepping.dump() );
            boost::property_tree::read_json( istr, pt_timestepping );

            ptree.add_child( "time_stepping", pt_timestepping );
        }

    }

    /**
     * \brief copy additional files in a directory (can be usefull in online loading)
     * \param dirDataBase : directory copy
     */
    void copyAdditionalModelFiles( std::string const& dir )
        {
            if ( !M_model )
                return;

            if ( !M_model->additionalModelData().empty() )
            {
                for ( auto & [key,mdata] : M_model->additionalModelData() )
                {
                    // TODO move this code in AdditionalModelData class
                    mdata.prepareSave();
                    std::string const& relPath = mdata.relativeFilePathInDatabase();
                    fs::path newFilePath = fs::path(dir)/relPath;
                    fs::path parentNewFilePath = newFilePath.parent_path();
                    if ( M_model->worldComm().isMasterRank() )
                    {
                        if ( !fs::exists( parentNewFilePath ) )
                            fs::create_directories( parentNewFilePath );
                        if ( mdata.template has<nl::json>() )
                        {
                            auto const& jsonData = mdata.template data<nl::json>();
                            fs::ofstream o(newFilePath);
                            o << jsonData.dump(/*1*/);
                        }
                        else if ( mdata.template has<fs::path>() )
                        {
                            fs::path const& inputPath = mdata.template data<fs::path>();
                            boost::system::error_code ec;
                            fs::copy_file( inputPath, newFilePath, fs::copy_option::overwrite_if_exists, ec );
                        }
                    }
                    mdata.setOnDisk();
                }
                M_model->worldComm().barrier();
            }
        }

    /**
     * \brief load CrbModel from json
     * \param input json filename
     */
    void loadJson( std::string const& filename, std::string const& childname = "" )
        {
            if ( !fs::exists( filename ) )
            {
                LOG(INFO) << "Could not find " << filename << std::endl;
                return;
            }

            // first the underlying model
            M_model->loadJson( filename, "crbmodel" );

#if 0 // VINCENT
            auto json_str_wo_comments = removeComments(readFromFile(filename));
            //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

            boost::property_tree::ptree ptree;
            std::istringstream istr( json_str_wo_comments );
            boost::property_tree::read_json( istr, ptree );
            if ( childname.empty() )
                this->setup( ptree );
            else
            {
                auto const& ptreeChild = ptree.get_child( childname );
                this->setup( ptreeChild );
            }
#else
            std::ifstream ifs( filename );
            nl::json jparsed = nl::json::parse(ifs,nullptr,true,true);
            this->setup( childname.empty()? jparsed : jparsed.at( childname ) );
#endif
        }

    void setup( boost::property_tree::ptree const& ptree )
        {
            //auto const& ptreeCrbModel = ptree.get_child( "crbmodel" );
            auto const& ptreeAffineDecomposition = /*ptreeCrbModel*/ptree.get_child( "affine-decomposition" );

            M_mMaxA.clear();
            for ( auto const& item : ptreeAffineDecomposition.get_child("mMaxA") )
                M_mMaxA.push_back( item.second.template get_value<int>() );
            M_Qa=M_mMaxA.size();

            M_mMaxM.clear();
            for ( auto const& item : ptreeAffineDecomposition.get_child("mMaxM") )
                M_mMaxM.push_back( item.second.template get_value<int>() );
            M_Qm=M_mMaxM.size();

            M_mMaxF.clear();
            M_Ql.clear();
            for ( auto const& item : ptreeAffineDecomposition.get_child("mMaxF") )
            {
                std::vector<int> sizeloaded;
                for ( auto const& item2 : item.second.get_child("") )
                    sizeloaded.push_back( item2.second.template get_value<int>() );
                M_mMaxF.push_back( sizeloaded );
                M_Ql.push_back( sizeloaded.size() );
            }
            M_Nl = M_mMaxF.size();

            M_mMaxLinearDecompositionA.clear();
            for ( auto const& item : ptreeAffineDecomposition.get_child("mMaxLinearDecompositionA") )
                M_mMaxLinearDecompositionA.push_back( item.second.template get_value<int>() );
            M_QLinearDecompositionA = M_mMaxLinearDecompositionA.size();

            M_alreadyCountAffineDecompositionTerms = true;
            M_is_initialized=true;
        }
    void setup( nl::json const& jarg )
        {
            if ( jarg.contains( "affine-decomposition" ) )
            {
                M_mMaxA.clear();
                M_mMaxM.clear();
                M_mMaxF.clear();
                M_Ql.clear();
                M_mMaxLinearDecompositionA.clear();

                auto const& j_affdec = jarg.at( "affine-decomposition" );
                if ( j_affdec.contains( "mMaxA" ) )
                {
                    for ( auto const& [j_mMaxA_key,j_mMaxA_val] : j_affdec.at( "mMaxA" ).items() )
                    {
                        if ( j_mMaxA_val.is_number_integer() )
                            M_mMaxA.push_back( j_mMaxA_val.template get<int>() );
                        else if ( j_mMaxA_val.is_string() )
                            if ( !j_mMaxA_val.template get<std::string>().empty() )
                                M_mMaxA.push_back( std::stoi( j_mMaxA_val.template get<std::string>() ) );
                    }
                }
                M_Qa = M_mMaxA.size();

                if ( j_affdec.contains( "mMaxM" ) )
                {
                    for ( auto const& [j_mMaxM_key,j_mMaxM_val] : j_affdec.at( "mMaxM" ).items() )
                    {
                        if ( j_mMaxM_val.is_number_integer() )
                            M_mMaxM.push_back( j_mMaxM_val.template get<int>() );
                        else if ( j_mMaxM_val.is_string() )
                            if ( !j_mMaxM_val.template get<std::string>().empty() )
                            M_mMaxM.push_back( std::stoi( j_mMaxM_val.template get<std::string>() ) );
                    }
                }
                M_Qm = M_mMaxM.size();

                if ( j_affdec.contains( "mMaxF" ) )
                {
                    for ( auto const& [j_mMaxF_key,j_mMaxF_val] : j_affdec.at( "mMaxF" ).items() )
                    {
                        std::vector<int> sizeloaded;
                        for ( auto const& [j_mMaxF_sub_key,j_mMaxF_sub_val] : j_mMaxF_val.items() )
                        {
                            if ( j_mMaxF_sub_val.is_number_integer() )
                                sizeloaded.push_back( j_mMaxF_sub_val.template get<int>() );
                            else if ( j_mMaxF_sub_val.is_string() )
                                if ( !j_mMaxF_sub_val.template get<std::string>().empty() )
                                    sizeloaded.push_back( std::stoi( j_mMaxF_sub_val.template get<std::string>() ) );
                        }
                        M_mMaxF.push_back( sizeloaded );
                        M_Ql.push_back( sizeloaded.size() );
                    }
                }
                M_Nl = M_mMaxF.size();

                if ( j_affdec.contains( "mMaxLinearDecompositionA" ) )
                {
                    for ( auto const& [j_mMaxLinearDecompositionA_key,j_mMaxLinearDecompositionA_val] : j_affdec.at( "mMaxLinearDecompositionA" ).items() )
                    {
                        if ( j_mMaxLinearDecompositionA_val.is_number_integer() )
                            M_mMaxLinearDecompositionA.push_back( j_mMaxLinearDecompositionA_val.template get<int>() );
                        else if ( j_mMaxLinearDecompositionA_val.is_string() )
                            if ( !j_mMaxLinearDecompositionA_val.template get<std::string>().empty() )
                                M_mMaxLinearDecompositionA.push_back( std::stoi( j_mMaxLinearDecompositionA_val.template get<std::string>() ) );
                    }
                }
                M_QLinearDecompositionA = M_mMaxLinearDecompositionA.size();

                M_alreadyCountAffineDecompositionTerms = true;
            }

            M_is_initialized=true;
        }

    void postOnlineSetup( nl::json const& jarg )
        {
            //std::cout << "postOnlineSetup : " << jarg.dump(1) << std::endl;

            if ( jarg.contains( "time_stepping" ) )
            {
                auto const& j_timestepping = jarg.at( "time_stepping" );

                std::map<std::string,double> mapTimeDataDouble;
                for ( std::string const& dataName : {"time_step","time_initial","time_final"} )
                {
                    auto const& j_timestepping_dn = j_timestepping.at( dataName );
                    if ( j_timestepping_dn.is_string() )
                        mapTimeDataDouble[dataName] = std::stod( j_timestepping_dn.template get<std::string>() );
                    else if ( j_timestepping_dn.is_number() )
                        mapTimeDataDouble[dataName] = j_timestepping_dn.template get<double>();
                }
                M_onlineTimeStepping = bdf(_space=this->rBFunctionSpace(),
                                           _initial_time=mapTimeDataDouble.at("time_initial"),
                                           _final_time=mapTimeDataDouble.at("time_final"),
                                           _time_step=mapTimeDataDouble.at("time_step"),
                                           _save=false);
            }

            //M_model->postOnlineSetup( jarg, *this );

        }

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

    void updateRbSpaceContextEim()
    {
        M_model->updateRbSpaceContextEim();
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

    typename model_type::deim_vector_type deimVector() const
    {
        return M_model->deimVector();
    }
    typename model_type::mdeim_vector_type mdeimVector() const
    {
        return M_model->mdeimVector();
    }
    void updateRbInDeim( std::vector<std::vector<int>> subN={} )
    {
        auto deim_vector = this->deimVector();
        auto mdeim_vector = this->mdeimVector();

        for ( auto deim : deim_vector )
            deim->updateRb( rBFunctionSpace(), subN );
        for ( auto mdeim : mdeim_vector )
            mdeim->updateRb( rBFunctionSpace(), subN );
    }

    struct ComputeNormL2InCompositeCase
    {
        typedef double result_type;

        ComputeNormL2InCompositeCase( element_type const composite_u1 ,
                                      element_type const composite_u2 )
        :
            M_composite_u1 ( composite_u1 ),
            M_composite_u2 ( composite_u2 )
        {}

        template< typename T >
        result_type
        operator()( result_type const& r, const T& t ) const
        {
            auto u1 = M_composite_u1.template element< T::value >();
            auto u2 = M_composite_u2.template element< T::value >();
            double norm  = normL2(_range=u1.functionSpace()->template rangeElements< 0 >(),
                                  _expr=( idv(u1)-idv(u2) ) );
            return r + norm;
        }
    private :
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
        double norm = normL2(_range=u2.functionSpace()->template rangeElements<0>(),
                             _expr=( idv(u1)-idv(u2) ) );
        return norm;
    }
    double computeNormL2( element_type u1 , element_type u2 , mpl::bool_<true>)
    {
        //in that case we accumulate the norm of every components
        index_vector_type index_vector;
        return fusion::fold( index_vector, double(0.), ComputeNormL2InCompositeCase( u1,u2 ) );
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

    operatorcomposite_ptrtype operatorCompositeLightA() const
    {
        return M_model->operatorCompositeLightA();
    }
    //linear form
    std::vector< functionalcomposite_ptrtype > functionalCompositeLightF() const
    {
        return M_model->functionalCompositeLightF();
    }
    //mass matrix
    operatorcomposite_ptrtype operatorCompositeLightM() const
    {
        return operatorCompositeLightM( mpl::bool_<model_type::is_time_dependent>() );
    }
    operatorcomposite_ptrtype operatorCompositeLightM( mpl::bool_<true> ) const
    {
        return M_model->operatorCompositeLightM();
    }
    operatorcomposite_ptrtype operatorCompositeLightM( mpl::bool_<false> ) const
    {
        bool constructed_by_model = M_model->constructOperatorCompositeM();
        bool light_version=true;
        if( constructed_by_model )
            return M_model->operatorCompositeLightM();
        else
            return preAssembleMassMatrix( light_version );
    }


    /**
     * \brief count numer of terms in the affine decomposition
     * either using vectors from affine decomposition if matrices are stored
     * or using operators free
     */
    void countAffineDecompositionTerms()
    {
        if( M_alreadyCountAffineDecompositionTerms && !M_useSER)
            return;
        else
            M_alreadyCountAffineDecompositionTerms=true;

        if ( M_Aqm.size() > 0 )
        {
            if ( this->hasEimError() )
            {

                M_has_eim=true;
                //when using EIM we need to be careful
                //if we want to estimate the error then there
                //is an "extra" basis fuction.
                auto all_errors = eimInterpolationErrorEstimation();
                auto errorsM = all_errors.template get<0>();
                auto errorsA = all_errors.template get<1>();
                auto errorsF = all_errors.template get<2>();
                std::map<int,double>::iterator it;
                auto endA = errorsA.end();
                auto endM = errorsM.end();

                M_Qm=M_Mqm.size();
                M_mMaxM.resize(M_Qm);
                for(int q=0; q<M_Qm; q++)
                {
                    auto itq = errorsM.find(q);
                    if( itq != endM) //found
                    {
                        //we have an eim error estimation
                        M_mMaxM[q]=M_Mqm[q].size()-1;
                    }
                    else //not found
                    {
                        M_mMaxM[q]=M_Mqm[q].size();
                    }
                }

                M_Qa=M_Aqm.size();
                M_mMaxA.resize(M_Qa);
                for(int q=0; q<M_Qa; q++)
                {
                    auto itq = errorsA.find(q);
                    if( itq != endA ) //found
                    {
                        //we have an eim error estimation
                        M_mMaxA[q]=M_Aqm[q].size()-1;
                    }
                    else //not found
                    {
                        M_mMaxA[q]=M_Aqm[q].size();
                    }
                }

                M_Nl=M_Fqm.size();
                M_Ql.resize(M_Nl);
                M_mMaxF.resize(M_Nl);
                for(int output=0; output<M_Nl; output++)
                {
                    M_Ql[output]=M_Fqm[output].size();
                    M_mMaxF[output].resize(M_Ql[output]);
                    auto endF = errorsF[output].end();
                    for(int q=0; q<M_Ql[output]; q++)
                    {
                        auto itq = errorsF[output].find(q);
                        if( itq != endF ) //found
                        {
                            //we have an eim error estimation
                            M_mMaxF[output][q]=M_Fqm[output][q].size()-1;
                        }
                        else //not found
                        {
                            M_mMaxF[output][q]=M_Fqm[output][q].size();
                        }
                    }

                    M_linearAqm = M_model->computeLinearDecompositionA();
                    M_QLinearDecompositionA = M_linearAqm.size();
                    M_mMaxLinearDecompositionA.resize( M_QLinearDecompositionA );
                    for(int q=0; q<M_QLinearDecompositionA; q++)
                    {
                        M_mMaxLinearDecompositionA[q] = M_linearAqm[q].size();
                    }
                }
            }//EIM error
            else
            {
                M_Qm=M_Mqm.size();
                M_mMaxM.resize(M_Qm);
                for(int q=0; q<M_Qm; q++)
                {
                    M_mMaxM[q]=M_Mqm[q].size();
                    // if( M_mMaxM[q] > 1 )
                    //     M_has_eim=true;
                }

                M_Qa=M_Aqm.size();
                M_mMaxA.resize(M_Qa);
                for(int q=0; q<M_Qa; q++)
                {
                    M_mMaxA[q]=M_Aqm[q].size();
                    // if( M_mMaxA[q] > 1 )
                    //     M_has_eim=true;
                }

                M_Nl=M_Fqm.size();
                M_Ql.resize(M_Nl);
                M_mMaxF.resize(M_Nl);
                for(int output=0; output<M_Nl; output++)
                {
                    M_Ql[output]=M_Fqm[output].size();
                    M_mMaxF[output].resize(M_Ql[output]);
                    for(int q=0; q<M_Ql[output]; q++)
                    {
                        M_mMaxF[output][q]=M_Fqm[output][q].size();
                        //std::cout << "countAffineDecomposition terms : M_mMaxF[" << output << "][" << q << "] = " << M_mMaxF[output][q] << std::endl;
                        // if( M_mMaxF[output][q] > 1 )
                        //     M_has_eim=true;
                    }
                }

                if( M_has_eim && !is_linear )
                {
                    M_linearAqm = M_model->computeLinearDecompositionA();
                    M_QLinearDecompositionA = M_linearAqm.size();
                    M_mMaxLinearDecompositionA.resize( M_QLinearDecompositionA );
                    for(int q=0; q<M_QLinearDecompositionA; q++)
                    {
                        M_mMaxLinearDecompositionA[q] = M_linearAqm[q].size();
                    }
                }

            }
        }
        else
        {
            //operators free
            auto compositeAlight = operatorCompositeLightA();
            if( compositeAlight )
            {
                auto compositeMlight = operatorCompositeLightM();
                M_mMaxM = compositeMlight->countAllContributions();
                M_Qm=M_mMaxM.size();
                M_mMaxA = compositeAlight->countAllContributions();
                M_Qa=M_mMaxA.size();
                auto vector_compositeFlight = functionalCompositeLightF();
                int number_outputs = vector_compositeFlight.size();
                M_Nl=number_outputs;
                M_mMaxF.resize(number_outputs);
                M_Ql.resize(number_outputs);
                for(int output=0; output<number_outputs; output++)
                {
                    auto compositeF = vector_compositeFlight[output];
                    M_mMaxF[output] = compositeF->countAllContributions();
                    M_Ql[output] = M_mMaxF[output].size();
                }
            }//light version
            else
            {
                M_has_eim=true;
                auto compositeM = operatorCompositeM();
                M_mMaxM = compositeM->countAllContributions();
                M_Qm=M_mMaxM.size();

                auto compositeA = operatorCompositeA();
                //it is important to check that the user provides operators free
                //else we don't have any affine decomposition
                CHECK( compositeA )<<"Very important ERROR !!! You have not implemented computeAffineDecomposition or operatorCompositeA !\n";
                M_mMaxA = compositeA->countAllContributions();
                M_Qa=M_mMaxA.size();
                auto vector_compositeF = functionalCompositeF();
                int number_outputs = vector_compositeF.size();
                M_Nl=number_outputs;
                M_mMaxF.resize(number_outputs);
                M_Ql.resize(number_outputs);
                for(int output=0; output<number_outputs; output++)
                {
                    auto compositeF = vector_compositeF[output];
                    M_mMaxF[output] = compositeF->countAllContributions();
                    M_Ql[output] = M_mMaxF[output].size();
                }
            }
        }//classic version
    }

    incomplete_betaqm_type extendBetaCoefficients( betaq_type beta_light )
    {
        return extendBetaCoefficients( beta_light , mpl::bool_< is_time_dependent >() );
    }
    incomplete_betaqm_type extendBetaCoefficients( betaq_type beta_light , mpl::bool_<true> )
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;

        auto betaMq = beta_light.template get<0>();
        auto betaAq = beta_light.template get<1>();
        auto betaFq = beta_light.template get<2>();
        int qsize=betaMq.size();
        betaMqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            betaMqm[q].resize(1);
            betaMqm[q][0]=betaMq[q];
        }
        qsize=betaAq.size();
        betaAqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            betaAqm[q].resize(1);
            betaAqm[q][0]=betaAq[q];
        }
        int nboutputs=betaFq.size();
        betaFqm.resize(nboutputs);
        for(int output=0; output<nboutputs;output++)
        {
            qsize=betaFq[output].size();
            betaFqm[output].resize(qsize);

            for(int q=0; q<qsize; q++)
            {
                betaFqm[output][q].resize(1);
                betaFqm[output][q][0]=betaFq[output][q];
            }
        }

        return boost::make_tuple( betaMqm, betaAqm , betaFqm );

    }
    incomplete_betaqm_type extendBetaCoefficients( betaq_type beta_light , mpl::bool_<false> )
    {
        beta_vector_type betaAqm;
        std::vector<beta_vector_type>  betaFqm;

        auto betaAq = beta_light.template get<0>();
        auto betaFq = beta_light.template get<1>();
        int qsize=betaAq.size();
        betaAqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            betaAqm[q].resize(1);
            betaAqm[q][0]=betaAq[q];
        }
        int nboutputs=betaFq.size();
        betaFqm.resize(nboutputs);
        for(int output=0; output<nboutputs;output++)
        {
            qsize=betaFq[output].size();
            betaFqm[output].resize(qsize);

            for(int q=0; q<qsize; q++)
            {
                betaFqm[output][q].resize(1);
                betaFqm[output][q][0]=betaFq[output][q];
            }
        }

        return boost::make_tuple( betaAqm , betaFqm );

    }


    void extendAffineDecomposition( affine_decomposition_light_type tuple )
    {
        return extendAffineDecomposition(  tuple , mpl::bool_< is_time_dependent >() );
    }
    void extendAffineDecomposition( affine_decomposition_light_type tuple , mpl::bool_<true> )
    {
        auto Xh=M_model->functionSpace();
        auto Mq=tuple.template get<0>();
        auto Aq=tuple.template get<1>();
        auto Fq=tuple.template get<2>();

        int qsize=Mq.size();
        M_Mqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            M_Mqm[q].resize(1);
            M_Mqm[q][0]=M_backend->newMatrix( _test=Xh, _trial=Xh );
            M_Mqm[q][0]=Mq[q];
        }
        qsize=Aq.size();
        M_Aqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            M_Aqm[q].resize(1);
            M_Aqm[q][0]=M_backend->newMatrix( _test=Xh, _trial=Xh );
            M_Aqm[q][0]=Aq[q];
        }
        int nb_output=Fq.size();
        M_Fqm.resize(nb_output);
        for(int output=0; output<nb_output; output++)
        {
            qsize=Fq[output].size();
            M_Fqm[output].resize( qsize );
            for(int q=0; q<qsize; q++)
            {
                M_Fqm[output][q].resize(1);
                M_Fqm[output][q][0]=M_backend->newVector( Xh );
                M_Fqm[output][q][0]=Fq[output][q];
            }
        }
    }
    void extendAffineDecomposition( affine_decomposition_light_type tuple , mpl::bool_<false> )
    {
        auto Xh=M_model->functionSpace();
        auto Aq=tuple.template get<0>();
        auto Fq=tuple.template get<1>();

        int qsize=Aq.size();
        M_Aqm.resize( qsize );
        for(int q=0; q<qsize; q++)
        {
            M_Aqm[q].resize(1);
            M_Aqm[q][0]=M_backend->newMatrix( _test=Xh, _trial=Xh );
            M_Aqm[q][0]=Aq[q];
        }
        int nb_output=Fq.size();
        M_Fqm.resize(nb_output);
        for(int output=0; output<nb_output; output++)
        {
            qsize=Fq[output].size();
            M_Fqm[output].resize( qsize );
            for(int q=0; q<qsize; q++)
            {
                M_Fqm[output][q].resize(1);
                M_Fqm[output][q][0]=M_backend->newVector( Xh );
                M_Fqm[output][q][0]=Fq[output][q];
            }
        }
    }



    /**
     * \brief compute the monolithic formulation
     * that is to say with no affine decomposition
     */
    monolithic_type computeMonolithicFormulation( parameter_type const& mu )
    {
        return computeMonolithicFormulation( mu, mpl::bool_<model_type::is_time_dependent>() );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu , mpl::bool_<true> )
    {
        boost::tie( M_monoM, M_monoA, M_monoF ) = M_model->computeMonolithicFormulation( mu );
        return boost::make_tuple( M_monoM, M_monoA, M_monoF );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu , mpl::bool_<false> )
    {
        boost::tie( M_monoA, M_monoF ) = M_model->computeMonolithicFormulation( mu );
        return boost::make_tuple( M_monoM, M_monoA, M_monoF );
    }

    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& u )
    {
        return computeMonolithicFormulationU( mu, u, mpl::bool_<model_type::is_time_dependent>() );
    }
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& u , mpl::bool_<true> )
    {
        boost::tie( M_monoM, M_monoA, M_monoF ) = M_model->computeMonolithicFormulationU( mu, u );
        return boost::make_tuple( M_monoM, M_monoA, M_monoF );
    }
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& u ,  mpl::bool_<false> )
    {
        boost::tie( M_monoA, M_monoF ) = M_model->computeMonolithicFormulationU( mu, u );
        return boost::make_tuple( M_monoM, M_monoA, M_monoF );
    }


    /*
     * return true if the model use EIM and need to comute associated error
     * usefull to the CRB error estimation
     */
    bool hasEimError()
    {
        auto all_errors = eimInterpolationErrorEstimation();
        auto errorsM = all_errors.template get<0>();
        auto errorsA = all_errors.template get<1>();
        auto errorsF = all_errors.template get<2>();
        int sizeM = errorsM.size();
        int sizeA = errorsA.size();
        int sizeF = errorsF.size();
        bool b=false;
        if( sizeM > 0 || sizeA > 0 || sizeF > 0 )
        {
            b=true;
        }
        return b;
    }

    /*
     * return true if the model uses EIM
     */
    bool hasEim()
    {
        return M_has_eim || hasDeim();
    }

    bool hasDeim()
    {
        return M_model->hasDeim();
    }
    /*
     * return true if the model uses SER
     */
    bool useSER()
    {
        return M_useSER;
    }


    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN )
    {
        return eimInterpolationErrorEstimation( mu, uN, mpl::bool_<model_type::is_time_dependent>() );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN, mpl::bool_<true> )
    {
        return M_model->eimInterpolationErrorEstimation( mu , uN );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN, mpl::bool_<false> )
    {
        std::map< int, double > errorsM;
        std::map< int, double > errorsA;
        std::vector< std::map< int, double > > errorsF;
        boost::tie(  errorsA , errorsF ) = M_model->eimInterpolationErrorEstimation( mu ,uN );
        return boost::make_tuple( errorsM, errorsA, errorsF );
    }

    eim_interpolation_error_type eimInterpolationErrorEstimation( )
    {
        return eimInterpolationErrorEstimation( mpl::bool_<model_type::is_time_dependent>() );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( mpl::bool_<true> )
    {
        return M_model->eimInterpolationErrorEstimation( );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( mpl::bool_<false> )
    {
        std::map< int, double > errorsM;
        std::map< int, double > errorsA;
        std::vector< std::map< int, double > > errorsF;
        boost::tie(  errorsA , errorsF ) = M_model->eimInterpolationErrorEstimation( );
        return boost::make_tuple( errorsM, errorsA, errorsF );
    }



    /**
     * \brief compute the linear part of the affine decomposition of a()
     * it is needed to compute the norm of the error for nonlinear models using EIM
     */
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA()
    {
        return M_model->computeLinearDecompositionA();
    }

    /**
     * \brief compute beta coefficients associetd to the linear part of the affine decomposition of a()
     * it is needed to compute the norm of the error for nonlinear models using EIM
     */
    beta_vector_type computeBetaLinearDecompositionA ( parameter_type const& mu ,  double time=0 )
    {
        return M_model->computeBetaLinearDecompositionA( mu , time );
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
        std::vector< sparse_matrix_ptrtype > Aq;
        std::vector< sparse_matrix_ptrtype > Mq;
        std::vector< std::vector< vector_ptrtype > > Fq;

        //first check is model provides a small affine decomposition or not
        boost::tie( Mq, Aq, Fq ) = M_model->computeAffineDecompositionLight();

        if( Aq.size() == 0 )
        {
            boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();
            if( M_serErrorEstimation && M_crbUseNewton )
            {
                auto RF = M_model->computePicardAffineDecomposition();
                M_RF_Aqm = RF.template get<1>();
                M_RF_Fqm = RF.template get<2>();
            }
        }
        else
        {
            auto tuple=boost::make_tuple(Mq,Aq,Fq);
            this->extendAffineDecomposition( tuple );
        }

        this->countAffineDecompositionTerms();

        if( M_Aqm.size() == 0 )
        {
            //in that case the model must provides operators free

            auto compositeAlight = operatorCompositeLightA();
            if( compositeAlight )
            {

                auto compositeMlight = operatorCompositeLightM();
                int q_max = this->Qm();
                M_Mqm.resize( q_max);
                for(int q=0; q<q_max; q++)
                {
                    M_Mqm[q].resize(1);
                    auto operatorfree = compositeMlight->operatorlinear(q);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Mqm[q][0]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Mqm[q][0]);//fill the matrix
                }//q

                q_max = this->Qa();
                M_Aqm.resize( q_max);
                for(int q=0; q<q_max; q++)
                {
                    M_Aqm[q].resize(1);
                    auto operatorfree = compositeAlight->operatorlinear(q);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Aqm[q][0]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Aqm[q][0]);//fill the matrix
                }//q

                auto vector_compositeFlight = functionalCompositeLightF();
                int number_outputs = vector_compositeFlight.size();
                M_Fqm.resize(number_outputs);
                for(int output=0; output<number_outputs; output++)
                {
                    auto composite_f = vector_compositeFlight[output];
                    int q_max = this->Ql(output);
                    M_Fqm[output].resize( q_max);
                    for(int q=0; q<q_max; q++)
                    {
                        M_Fqm[output][q].resize(1);
                        auto operatorfree = composite_f->functionallinear(q);
                        auto space = operatorfree->space();
                        M_Fqm[output][q][0]= M_backend->newVector( space );
                        operatorfree->containerPtr(M_Fqm[output][q][0]);//fill the vector
                    }//q
                }//output
            }//end of light version
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
            }//end of classic version
        }

        return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        initial_guess_type initial_guess;
        std::vector< sparse_matrix_ptrtype > Aq;
        std::vector< std::vector< vector_ptrtype > > Fq;

        //first check is model provides a small affine decomposition or not
        boost::tie( Aq, Fq ) = M_model->computeAffineDecompositionLight();
        if( Aq.size() == 0 )
        {
            boost::tie( M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();

            if( M_serErrorEstimation && M_crbUseNewton )
            {
                auto RF = M_model->computePicardAffineDecomposition();
                M_RF_Aqm = RF.template get<0>();
                M_RF_Fqm = RF.template get<1>();
            }
        }
        else
        {
            auto tuple=boost::make_tuple(Aq,Fq);
            this->extendAffineDecomposition( tuple );
        }
        this->countAffineDecompositionTerms();

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

            auto compositeAlight = operatorCompositeLightA();
            if( compositeAlight )
            {
                int q_max = this->Qa();
                M_Aqm.resize( q_max);
                for(int q=0; q<q_max; q++)
                {
                    M_Aqm[q].resize(1);
                    auto operatorfree = compositeAlight->operatorlinear(q);
                    size_type pattern = operatorfree->pattern();
                    auto trial = operatorfree->domainSpace();
                    auto test=operatorfree->dualImageSpace();
                    M_Aqm[q][0]= M_backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
                    operatorfree->matPtr(M_Aqm[q][0]);//fill the matrix
                }//q

                auto vector_compositeFlight = functionalCompositeLightF();
                int number_outputs = vector_compositeFlight.size();
                M_Fqm.resize(number_outputs);
                for(int output=0; output<number_outputs; output++)
                {
                    auto composite_f = vector_compositeFlight[output];
                    int q_max = this->Ql(output);
                    M_Fqm[output].resize( q_max);
                    for(int q=0; q<q_max; q++)
                    {
                        M_Fqm[output][q].resize(1);
                        auto operatorfree = composite_f->functionallinear(q);
                        auto space = operatorfree->space();
                        M_Fqm[output][q][0]= M_backend->newVector( space );
                        operatorfree->containerPtr(M_Fqm[output][q][0]);//fill the vector
                    }//q
                }//output
            }// light version
            else
            {
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
                int number_outputs = M_Nl;
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
            }//classic version
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

        this->assembleInitialGuessV( initial_guess_v); //Compute \int (q_initialguess) v

        for(int q=0; q<M_InitialGuessVector.size(); q++)
        {
            for(int m=0; m<M_InitialGuessVector[q].size(); m++)
                *M_InitialGuessV[q][m] = *M_InitialGuessVector[q][m];
        }
        return M_InitialGuessV;
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
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Aqm[q][m]->transpose();

            return M_Aqm[q][m];
        }
        else
        {
            auto compositeLight=operatorCompositeLightA();
            if( compositeLight )
            {
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto opfree = compositeLight->operatorlinear(q);
                size_type pattern = opfree->pattern();
                auto trial = opfree->domainSpace();
                auto test = opfree->dualImageSpace();
                auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
                opfree->matPtr( matrix ); //assemble the matrix
                if ( transpose )
                    return matrix->transpose();
                return matrix;
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
    }

    /**
     * Returns the matrix \c Aqm ( if Newton is used, Aqm(q,m) returns the Jacobian)
     */
    std::vector< std::vector<sparse_matrix_ptrtype> > RF_Aqm(){ return M_RF_Aqm; }
    /**
     * Returns the matrix \c Fqm ( if Newton is used, Fqm(0,q,m) returns the Residual)
     */
    std::vector< std::vector< vector_ptrtype > > RF_Fqm(){ return M_RF_Fqm[0]; }

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
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Mqm[q][m]->transpose();

            return M_Mqm[q][m];
        }
        else
        {
            auto compositeLight=operatorCompositeLightM();
            if( compositeLight )
            {
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto opfree = compositeLight->operatorlinear(q);
                size_type pattern = opfree->pattern();
                auto trial = opfree->domainSpace();
                auto test = opfree->dualImageSpace();
                auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
                opfree->matPtr( matrix ); //assemble the matrix
                if ( transpose )
                    return matrix->transpose();
                return matrix;
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
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            if ( transpose )
                return M_Mqm[q][m]->transpose();

            return M_Mqm[q][m];
        }
        else
        {
            auto compositeLight=operatorCompositeLightM();
            if( compositeLight )
            {
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto opfree = compositeLight->operatorlinear(q);
                size_type pattern = opfree->pattern();
                auto trial = opfree->domainSpace();
                auto test = opfree->dualImageSpace();
                auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
                opfree->matPtr( matrix ); //assemble the matrix
                if ( transpose )
                    return matrix->transpose();
                return matrix;
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
     * \brief the inner product \f$linear a_{qm}(\xi_i, \xi_j) = \xi_j^T LinearA_{qm} \xi_i\f$
     *
     * \param q and m index of the component in the affine decomposition
     * \param xi_i an element of the function space
     * \param xi_j an element of the function space
     * \param transpose transpose \c LinearA_{qm}
     *
     * \return the inner product \f$linear a_qm(\xi_i, \xi_j) = \xi_j^T LinearA_{qm} \xi_i\f$
     */
    value_type linearDecompositionAqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
    {
        bool stock = M_stockMatrices;
        if( stock )
        {
            //in this case matrices have already been stocked
            return M_linearAqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        CHECK( stock )<<" This check is here because nothing is done to deal with linear decomposition of a (model using EIM) when using "
                      <<" operators free and option crb.stock-matrices=false.\nFor now we suppose that crb.stock-matrices=true.\n "
                      <<" So make sure that all developments have been done to deal with crb.stock-matrices=false before delete this CHECK.\n";
        return 0;
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
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            return M_Aqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        else
        {
            auto compositeLight=operatorCompositeLightA();
            if( compositeLight )
            {
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto opfree = compositeLight->operatorlinear(q);
                size_type pattern = opfree->pattern();
                auto trial = opfree->domainSpace();
                auto test = opfree->dualImageSpace();
                auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
                opfree->matPtr( matrix ); //assemble the matrix

                return matrix->energy( xi_j, xi_i, transpose );
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
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            //M_Mqm[q][m]->printMatlab("mass_matrix.m");
            return M_Mqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        else
        {
            auto compositeLight=operatorCompositeLightM();
            if( compositeLight )
            {
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto opfree = compositeLight->operatorlinear(q);
                size_type pattern = opfree->pattern();
                auto trial = opfree->domainSpace();
                auto test = opfree->dualImageSpace();
                auto matrix = M_backend->newMatrix( _trial=trial , _test=test , _pattern=pattern );
                opfree->matPtr( matrix ); //assemble the matrix
                return matrix->energy( xi_j, xi_i, transpose );
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
        if( M_stockMatrices )
        {
            return M_Fqm[l][q][m];
        }
        else
        {
            auto compositeLight=operatorCompositeLightA();
            if( compositeLight )
            {
                auto vector_compositeLight=functionalCompositeLightF();
                CHECK( m == 0 )<<"Error, try to use light version of operators free with EIM \n";
                auto composite = vector_compositeLight[l];
                auto functional = composite->functionallinear(q);
                auto space = functional->space();
                auto vector = M_backend->newVector( space );
                functional->containerPtr( vector );
                return vector;
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
    }

    value_type Jqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
    {
        if( M_stockMatrices )
        {
            //in this case matrices have already been stocked
            if ( M_Jqm.empty() )
                return M_Aqm[q][m]->energy( xi_j, xi_i, transpose );
            else
                return M_Jqm[q][m]->energy( xi_j, xi_i, transpose );
        }
        else
        {
            CHECK( false ) << "TODO : Jqm operatorComposite";
            return 0;
        }
    }

    element_ptrtype InitialGuessQm( uint16_type q, int m ) const
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

        if( M_stockMatrices )
        {
            result = inner_product( *M_Fqm[l][q][m] , xi );
        }
        else
        {
            auto compositeLight=operatorCompositeLightA();
            if( compositeLight )
            {
                auto vector_compositeLight=functionalCompositeLightF();
                auto composite = vector_compositeLight[l];
                auto functional = composite->functionallinear(q);
                auto space = functional->space();
                auto vector = M_backend->newVector( space );
                functional->containerPtr( vector );
                result = inner_product( *vector, xi ) ;
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

        if( M_stockMatrices )
        {
            result = inner_product( *M_Fqm[l][q][m] , xi );
        }
        else
        {
            auto compositeLight=operatorCompositeLightA();
            if( compositeLight )
            {
                auto vector_compositeLight=functionalCompositeLightF();
                auto composite = vector_compositeLight[l];
                auto functional = composite->functionallinear(q);
                auto space = functional->space();
                auto vector = M_backend->newVector( space );
                functional->containerPtr( vector );
                result = inner_product( *vector, xi ) ;
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
        }

        return result;
    }

    value_type InitialGuessVqm( uint16_type q, uint16_type m, element_type const& xi)
    {
        value_type result=0;
        result = inner_product( *M_InitialGuessV[q][m] , xi );
        return result;
    }


    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& X, vector_type const& Y )
    {
        return M_inner_product_matrix->energy( X, Y );
    }
    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        return M_inner_product_matrix->energy( X, Y );
    }


    /**
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_type const& X, vector_type const& Y )
    {
        auto M = M_model->massMatrix();
        return M->energy( X, Y );
    }
    /**
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        auto M = M_model->massMatrix();
        return M->energy( X, Y );
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
        return M_inner_product_matrix->energy( X, Y );
        //auto M = M_model->innerProductForPod();
        //return M->energy( X, Y );
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
        return M_inner_product_matrix->energy( X, Y );
        //return M_model->scalarProductForPod( X, Y );
    }
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y , mpl::bool_<false> )
    {
        return 0;
    }

    preconditioner_ptrtype preconditionerPrimal() const { return M_preconditioner_primal; }
    preconditioner_ptrtype preconditionerDual() const { return M_preconditioner_dual; }
    preconditioner_ptrtype preconditionerL2() const { return M_preconditioner_l2; }
    preconditioner_ptrtype preconditionerPrimal()  { return M_preconditioner_primal; }
    preconditioner_ptrtype preconditionerDual()  { return M_preconditioner_dual; }
    preconditioner_ptrtype preconditionerL2()  { return M_preconditioner_l2; }

    /**
     * solve the model for a given parameter \p mu
     */
    virtual element_type solve( parameter_type const& mu )
    {
        element_type solution;// = M_model->functionSpace()->element();
        if( is_linear )
            solution = this->solveFemUsingAffineDecompositionFixedPoint( mu );
        else
            solution = M_model->solve( mu );
        return solution;
    }


    /**
     * solve \f$M u = f\f$ where \f$ M \f$ is the matrix associated to
     * the energy norm
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f )
    {
        M_backend_l2->solve( _matrix=M_inner_product_matrix,  _solution=u, _rhs=f , _prec=M_preconditioner_l2 );
        //return M_model->l2solve( u, f );
    }


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

    gmsh_ptrtype createStructuredGrid( std::vector<int> components_vary, std::vector<parameter_type> extremums,
                                       std::vector<int> cuttings, std::vector<double> time_cuttings, bool time_vary )
    {
        return M_model->createStructuredGrid( components_vary, extremums, cuttings, time_cuttings, time_vary );
    }

    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector, std::string filename , std::string extra="")
    {
        return M_model->writeConvergenceStatistics( vector, filename , extra );
    }

    void writeVectorsExtremumsRatio(std::vector< vectorN_type > const& vector1, std::vector< vectorN_type > const& vector2, std::string filename )
    {
        return M_model->writeVectorsExtremumsRatio( vector1, vector2, filename );
    }

    std::shared_ptr<TSBase> onlineTimeStepping() const { return M_onlineTimeStepping; }

    bdf_ptrtype /*const&*/ bdfModel() const
    {
        return M_model->bdfModel();
    }

    int numberOfTimeStep() const override { return 100+1; } // return M_numberOfTimeStep; }

    double timeStep() const override
    {
        return 5;//VINCENT
        double timestep = 1e30;
        if ( !this->isSteady() )
            timestep = this->bdfModel()->timeStep();
        return timestep;
    }
    double timeInitial() const override
    {
        return 0;//VINCENT
        double timeinitial = 0.;
        if ( !this->isSteady() )
            timeinitial = this->bdfModel()->timeInitial();
        return timeinitial;
    }
    double timeFinal() const override
    {
        return 500;//VINCENT
        double timefinal=1e30;
        if ( !this->isSteady() )
            timefinal = this->bdfModel()->timeFinal();
        return timefinal;
    }
    int timeOrder() const
    {
        return 1;//VINCENT
        int order = 0;
        if ( !this->isSteady() )
            order = this->bdfModel()->timeOrder();
        return order;
    }


    bool isSteady() const override
    {
        return M_isSteadyModel;
    }

    bool isLinear() const
    {
        return is_linear ;
    }

    virtual bool isTrilinear() const
        {return false;}

    bool hasZeroMeanPressure()
        { return this->M_model->hasZeroMeanPressure(); }

    void initializationField( element_ptrtype& initial_field,parameter_type const& mu )
    {
        return initializationField( initial_field,mu,mpl::bool_<model_type::is_time_dependent>() );
    }
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu,mpl::bool_<true> )
    {
        return M_model->initializationField( initial_field,mu );
    }
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu,mpl::bool_<false> ) {}


    typename model_type::displacement_field_ptrtype meshDisplacementField( parameter_type const& mu )
    {
        return M_model->meshDisplacementField(mu);
    }
    bool hasDisplacementField() const { return M_model->hasDisplacementField(); }
    //@}



protected:

    std::shared_ptr<CRBModelDB> M_crbModelDb;
    std::string M_prefix;
    int M_level = 0;

    //! affine decomposition terms for the left hand side
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_linearAqm;

    mutable std::vector< std::vector<element_ptrtype> > M_InitialGuessV;
    mutable std::vector< std::vector<vector_ptrtype> > M_InitialGuessVector;

    //! affine decomposition terms ( time dependent )
    mutable std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;

    //! affine decomposition terms for the right hand side
    std::vector< std::vector <std::vector<vector_ptrtype>> > M_Fqm;

    sparse_matrix_ptrtype M_monoA;
    sparse_matrix_ptrtype M_monoM;
    std::vector<vector_ptrtype> M_monoF;

    //! model
    model_ptrtype M_model;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Jqm;
    std::vector< std::vector< std::vector<vector_ptrtype> > > M_Rqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_RF_Aqm;
    std::vector< std::vector< std::vector<vector_ptrtype> > > M_RF_Fqm;


protected:


    bool M_is_initialized = false;

    //! mode for CRBModel


    backend_ptrtype M_backend;
    backend_ptrtype M_backend_primal;
    backend_ptrtype M_backend_dual;
    backend_ptrtype M_backend_l2;

    beta_vector_type M_dummy_betaMqm;

    /**
     * \brief given \p mu merge the Aq and Fq into A and F respectively
     *
     * \param mu the parameter at which the matrix A and vector F are assembled
     *
     */
    offline_merge_type offlineMerge( betaqm_type const& all_beta,  bool only_time_dependent_terms=false );
    offline_merge_type offlineMergeOnFly( betaqm_type const& all_beta , bool only_time_dependent_terms=false );

    void assembleMassMatrix( );
    void assembleMassMatrix( mpl::bool_<true> );
    void assembleMassMatrix( mpl::bool_<false> );

    operatorcomposite_ptrtype preAssembleMassMatrix( bool light_version=false ) const ;
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<true> , bool light_version ) const ;
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<false>, bool light_version ) const ;

    void assembleInitialGuessV( initial_guess_type & initial_guess );
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<true> );
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<false> );
    element_type M_u, M_v;

    preconditioner_ptrtype M_preconditioner_primal;
    preconditioner_ptrtype M_preconditioner_dual;
    preconditioner_ptrtype M_preconditioner_l2;
    bool M_fixedpointUseAitken;

    //number of terms in affine decomposition
    int M_Qa; //A
    int M_Qm; //M
    int M_Nl; //number of outputs
    int M_QLinearDecompositionA;
    std::vector<int> M_Ql;//F associated to given output
    std::vector<int> M_mMaxA;//number of sub-terms (using EIM)
    std::vector<int> M_mMaxM;//number of sub-terms (using EIM)
    std::vector<int> M_mMaxLinearDecompositionA;
    std::vector< std::vector<int> > M_mMaxF;//number of sub-terms (using EIM)

    bool M_alreadyCountAffineDecompositionTerms;

    sparse_matrix_ptrtype M_inner_product_matrix;

    bool M_isSteadyModel;
    int M_numberOfTimeStep;

    bool M_has_eim;
    bool M_useSER;

    bool M_usePrimalPc;
    bool M_useSymmMat;
    bool M_stockMatrices;
    bool M_serErrorEstimation;
    bool M_crbUseNewton;
    int M_fixedpointMaxIt;
    double M_fixedpointIncrTol;
    bool M_fixedpointVerbose;
    int M_outputIndex;
    bool M_useLinearModel;

    std::shared_ptr<TSBase> M_onlineTimeStepping;
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
    typedef std::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    PreAssembleMassMatrixInCompositeCase( element_type const& u ,
                                          element_type const& v,
                                          operatorcomposite_ptrtype op_mass )
        :
        M_composite_u ( u ),
        M_composite_v ( v ),
        M_op_mass( op_mass )
        {}

    template< typename T >
    void
    operator()( const T& t ) const
    {

        using namespace Feel::vf;

        auto u = M_composite_u.template element< T::value >();
        auto v = M_composite_v.template element< T::value >();
        auto Xh = M_composite_u.functionSpace();

        auto expr = integrate( _range=u.functionSpace()->template rangeElements< 0 >(),
                               _expr=trans( idt( u ) )*vf::id( v ) ) ;
        auto opfree = opLinearFree( _imageSpace=Xh, _domainSpace=Xh, _expr=expr );

        //each composant of the affine decomposition
        //is the subspace contribution
        int q=T::value;
        int m=0;
        auto tuple = boost::make_tuple( q , m );
        M_op_mass->addElement( tuple , opfree );
    }

private :
    element_type const&  M_composite_u;
    element_type const&  M_composite_v;
    operatorcomposite_ptrtype M_op_mass;

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


    AssembleMassMatrixInCompositeCase( element_type const& u ,
                                       element_type const& v ,
                                       CRBModel<ModelType>* crb_model)
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

        auto u = this->M_composite_u.template element< T::value >();
        auto v = this->M_composite_v.template element< T::value >();
        auto Xh = M_composite_u.functionSpace();
        auto Vh = u.functionSpace();

        form2( _test=Xh, _trial=Xh, _matrix=M_crb_model->Mqm(0,0) ) +=
            integrate( _range=Vh->template rangeElements< 0 >(),
                       _expr=trans( idt( u ) )*vf::id( v ) );

    }
private :
    element_type const&  M_composite_u;
    element_type const&  M_composite_v;
    mutable CRBModel<ModelType>* M_crb_model;
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

    AssembleInitialGuessVInCompositeCase( element_type const& v ,
                                          initial_guess_type const& initial_guess ,
                                          std::shared_ptr<CRBModel<ModelType> > crb_model)
        :
        M_composite_v ( v ),
        M_composite_initial_guess ( initial_guess ),
        M_crb_model ( crb_model )
    {}

    template< typename T >
    void
    operator()( const T& t ) const
    {
        using namespace Feel::vf;

        auto v = M_composite_v.template element< T::value >();
        auto Xh = M_composite_v.functionSpace();
        auto Vh = v.functionSpace();
        auto range = Vh->template rangeElements< 0 >();
        int q_max = M_crb_model->QInitialGuess();
        for(int q = 0; q < q_max; q++)
        {
            int m_max = M_crb_model->mMaxInitialGuess(q);
            for( int m = 0; m < m_max ; m++)
            {
                auto view = M_composite_initial_guess[q][m]->template element< T::value >();
                form1( _test=Xh, _vector=M_crb_model->InitialGuessVector(q,m) ) +=
                    integrate ( _range=range, _expr=trans( idv( view ) )*vf::id( v ) );
            }
        }
    }
private :
    element_type const&  M_composite_v;
    initial_guess_type const&  M_composite_initial_guess;
    mutable std::shared_ptr<CRBModel<ModelType> > M_crb_model;
};






//create a vector of preassemble objects
//for the mass matrix

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix( bool light_version ) const
{
    static const bool is_composite = functionspace_type::is_composite;
    return preAssembleMassMatrix( mpl::bool_< is_composite >() , light_version );
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix( mpl::bool_<false> , bool light_version ) const
{

    auto Xh = M_model->functionSpace();

    auto expr = integrate( _range=Xh->template rangeElements<0>() , _expr=inner( idt( M_u ),vf::id( M_v ) ) );
    auto op_mass = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh  );
    auto opfree = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr );
    opfree->setName("mass operator (automatically created)");
    //in this case, the affine decompositon has only one element
    if( light_version )
    {
        //note that only time independent problems need to
        //inform if we use light version or not of the affine decomposition
        //i.e. if we use EIM expansion or not
        //because transient models provide already a decomposition of the mass matrix
        //so it is not automatically created
        op_mass->addElement( 0 , opfree );
    }
    else
    {
        op_mass->addElement( boost::make_tuple(0,0) , opfree );
    }
    return op_mass;
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::operatorcomposite_ptrtype
CRBModel<TruthModelType>::preAssembleMassMatrix( mpl::bool_<true> , bool light_version ) const
{
    auto Xh = M_model->functionSpace();
    operatorcomposite_ptrtype op_mass = opLinearComposite( _imageSpace=Xh, _domainSpace=Xh );
    index_vector_type index_vector;
    fusion::for_each( index_vector, PreAssembleMassMatrixInCompositeCase<TruthModelType>( M_u , M_v, op_mass ) );
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
    form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[0][0] ) =
        integrate( _range=Xh->template rangeElements<0>(), _expr=inner(idt( M_u ),vf::id( M_v ) )  );
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

    //M_Mqm[0][0]->printMatlab("mass_matrix_before.m");
    fusion::for_each( index_vector, AssembleMassMatrixInCompositeCase<TruthModelType>( M_u , M_v , this ) );

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
            M_InitialGuessVector[q][m] = this->newVector();
        }
    }

    index_vector_type index_vector;
    fusion::for_each( index_vector, AssembleInitialGuessVInCompositeCase<TruthModelType>( M_v , initial_guess , this->shared_from_this() ) );

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
    auto range = Xh->template rangeElements<0>();

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
            M_InitialGuessVector[q][m] = this->newVector();
            form1( _test=Xh, _vector=M_InitialGuessVector[q][m]) =
                integrate( _range=range, _expr=inner( idv( initial_guess[q][m] ),vf::id( M_v ) )  );
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
CRBModel<TruthModelType>::offlineMergeOnFly(betaqm_type const& all_beta, bool only_time_dependent_terms )
{

    //recovery from the model of free operators Aqm in A = \sum_q \sum_m \beta_qm( u ; mu ) A_qm( u, v )
    //idem with other operator / functional

    operatorcomposite_ptrtype compositeM;
    operatorcomposite_ptrtype compositeA;
    std::vector< functionalcomposite_ptrtype > vector_compositeF;

    if( ! only_time_dependent_terms )
    {

        auto compositeAlight=this->operatorCompositeLightA();
        if( compositeAlight )
        {
            compositeA=compositeAlight;
            compositeM = this->operatorCompositeLightM();
            vector_compositeF = this->functionalCompositeLightF();
        }
        else
        {
            compositeA = this->operatorCompositeA();
            compositeM = this->operatorCompositeM();
            vector_compositeF = this->functionalCompositeF();
        }
        //acces to beta coefficients
        auto beta_M = all_beta.template get<0>();
        auto beta_A = all_beta.template get<1>();

        //associate beta coefficients to operators
        compositeA->setScalars( beta_A );
        compositeM->setScalars( beta_M );
    }

    //merge
    auto A = this->newMatrix();
    auto M = this->newMatrix();
    if( ! only_time_dependent_terms )
    {
        compositeA->sumAllMatrices( A );
        compositeM->sumAllMatrices( M );
    }

    //warning : beta_F is a vector of beta_coefficients
    auto beta_F = all_beta.template get<2>();
    std::vector<vector_ptrtype> F( Nl() );
    for(int output=0; output<Nl(); output++)
    {
        auto compositeF = vector_compositeF[output];
        compositeF->setScalars( beta_F[output] );
        F[output] = this->newVector();
        compositeF->sumAllVectors( F[output] );
    }

    return boost::make_tuple( M, A, F );

}


template<typename TruthModelType>
typename CRBModel<TruthModelType>::offline_merge_type
CRBModel<TruthModelType>::offlineMerge( betaqm_type const& all_beta , bool only_time_dependent_terms )
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

    if( ! only_time_dependent_terms )
    {
        //acces to beta coefficients
        auto beta_M = all_beta.template get<0>();
        auto beta_A = all_beta.template get<1>();

        A->zero();
        for ( size_type q = 0; q < Qa(); ++q )
        {
            for(size_type m = 0; m < mMaxA(q); ++m )
            {
                A->addMatrix( beta_A[q][m], M_Aqm[q][m] );
            }
        }

        if( Qm() > 0 )
        {
            for ( size_type q = 0; q < Qm(); ++q )
            {
                for(size_type m = 0; m < mMaxM(q) ; ++m )
                    M->addMatrix( beta_M[q][m] , M_Mqm[q][m] );
            }
        }
    }

    std::vector<vector_ptrtype> F( Nl() );

    //warning : beta_F is a vector of beta_coefficients
    auto beta_F = all_beta.template get<2>();

    for ( size_type l = 0; l < Nl(); ++l )
    {
        F[l] = M_backend->newVector( M_model->functionSpace() );
        F[l]->zero();

        for ( size_type q = 0; q < Ql( l ); ++q )
        {
            for ( size_type m = 0; m < mMaxF(l,q); ++m )
            {
                F[l]->add( beta_F[l][q][m] , M_Fqm[l][q][m] );
            }
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
        if( M_usePrimalPc )
        {
            M_preconditioner_primal->setMatrix( A );
            M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs , _prec=M_preconditioner_primal);
        }
        else
        {
            M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs, _rebuild=true);
        }
        mybdf->shiftRight(u);
    }

    return u;

}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemMonolithicFormulation( parameter_type const& mu )
{
    auto Xh= this->functionSpace();

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype A;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    element_ptrtype InitialGuess = Xh->elementPtr();

    auto u = Xh->element();

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
        // !!
        //some stuff needs to be done here
        //to deal with non linear problems
        //and have an initial guess
        // !!
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

    auto uold = Xh->element();
    u=uold;

    int max_fixedpoint_iterations  = this->vm()["crb.fixedpoint.maxit"].template as<int>();
    double increment_fixedpoint_tol  = this->vm()["crb.fixedpoint.increment-tol"].template as<double>();
    int iter=0;
    double norm=0;

    for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
    {
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        iter=0;
        do {

            if( is_linear )
                boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );
            else
                boost::tie(M, A, F) = this->computeMonolithicFormulationU( mu , u );

            if( !isSteady() )
            {
                A->addMatrix( bdf_coeff, M );
                F[0]->addVector( *vec_bdf_poly, *M );
            }
            uold=u;
            if( M_usePrimalPc )
            {
                M_preconditioner_primal->setMatrix( A );
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] , _prec=M_preconditioner_primal);
            }
            else
            {
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0], _rebuild=true);
            }

            mybdf->shiftRight(u);

            if( is_linear )
                norm = 0;
            else
                norm = this->computeNormL2( uold , u );
            iter++;
        } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
    }

    return u;

}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu )
{
    auto Xh = this->functionSpace();
    auto u = Xh->element("u");
    auto uold = Xh->element("u_old");

    sparse_matrix_ptrtype A;
    std::vector<vector_ptrtype> F;

    int max_fixedpoint_iterations = M_fixedpointMaxIt;
    double increment_fixedpoint_tol = M_fixedpointIncrTol;
    bool fixedPointVerbose = M_fixedpointVerbose;
    double norm=0;
    int iter=0;


    if ( this->isSteady() )
    {
        if ( is_linear )
        {
            boost::tie(boost::tuples::ignore, A, F) = this->update( mu );
            if( M_usePrimalPc )
            {
                M_preconditioner_primal->setMatrix( A );
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] , _prec=M_preconditioner_primal);
            }
            else
            {
                auto ret = M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] );
                if  ( !ret.template get<0>() )
                {
                    LOG(INFO)<<"[CRB] WARNING : we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
                    Feel::cout << "Rebuild backend primal\n";
                    M_backend_primal = backend(_name="backend-primal", _rebuild=true );
                    M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] );
                }
            }
        }
        else
        {
            u = *this->assembleInitialGuess( mu );

            bool useAitkenRelaxation = M_fixedpointUseAitken;
            auto residual = Xh->element();
            auto aitkenRelax = aitken( _space=Xh, _tolerance=increment_fixedpoint_tol );
            aitkenRelax.initialize( _residual=residual, _currentElt=u );
            aitkenRelax.restart();
            bool fixPointIsFinished = false;
            do {
                uold = u;

                boost::tie(boost::tuples::ignore, A, F) = this->update( mu , u );

                if( M_usePrimalPc )
                {
                    M_preconditioner_primal->setMatrix( A );
                    M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] , _prec=M_preconditioner_primal);
                }
                else
                {
                    backend( _name="backend-primal")->solve( _matrix=A , _solution=u, _rhs=F[0] );
                }

                iter++;
                if ( !useAitkenRelaxation )
                {
                    norm = this->computeNormL2( uold , u );
                    fixPointIsFinished = norm < increment_fixedpoint_tol || iter >= max_fixedpoint_iterations;
                    if ( this->worldComm().isMasterRank() && fixedPointVerbose )
                        std::cout << "[solveFemUsingAffineDecompositionFixedPoint] iteration " << iter << ", increment_norm = " <<  norm << "\n";
                }
                else
                {
                    residual = u-uold;
                    aitkenRelax.apply2(_newElt=u,_residual=residual, _currentElt=u );
                    if ( this->worldComm().isMasterRank() && fixedPointVerbose )
                        aitkenRelax.printInfo();
                    ++aitkenRelax;
                    if ( aitkenRelax.isFinished() )
                        fixPointIsFinished=true;
                }
            } while( !fixPointIsFinished );//norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );

        }
        return u;
    }



    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype M;
    vector_ptrtype Rhs( M_backend->newVector( Xh ) );
    element_ptrtype initialGuess = Xh->elementPtr();

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
        //we want to have the initial guess given by function update
        initialGuess = this->assembleInitialGuess( mu ) ;
    }
    else
    {
        time_initial=this->timeInitial();
        time_step=this->timeStep();
        time_final=this->timeFinal();
        this->initializationField( initialGuess, mu );
    }

    mybdf->setTimeInitial( time_initial );
    mybdf->setTimeStep( time_step );
    mybdf->setTimeFinal( time_final );

    u=*initialGuess;

    double bdf_coeff ;
    auto vec_bdf_poly = M_backend->newVector( Xh );


    for( mybdf->start(*initialGuess); !mybdf->isFinished(); mybdf->next() )
    {
        iter=0;
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        do {
            if( is_linear )
                boost::tie(M, A, F) = this->update( mu , mybdf->time() );
            else
                boost::tie(M, A, F) = this->update( mu , u , mybdf->time() );
            *Rhs = *F[0];

            if( !isSteady() )
            {
                A->addMatrix( bdf_coeff, M );
                Rhs->addVector( *vec_bdf_poly, *M );
            }
            uold = u;
            if( M_usePrimalPc )
            {
                M_preconditioner_primal->setMatrix( A );
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs , _prec=M_preconditioner_primal);
            }
            else
            {
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs);
            }

            if( is_linear )
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
CRBModel<TruthModelType>::solveFemUsingAffineDecompositionNewton( parameter_type const& mu )
{
    sparse_matrix_ptrtype J = this->newMatrix();
    vector_ptrtype R = this->newVector();

    // auto initialguess = this->functionSpace()->elementPtr();
    // initialguess = this->assembleInitialGuess( mu ) ;
    auto initialguess = this->assembleInitialGuess( mu ) ;


    boost::tie( boost::tuples::ignore , M_Jqm, M_Rqm ) = this->computeAffineDecomposition();

    M_backend_primal->nlSolver()->jacobian = std::bind( &CRBModel<TruthModelType>::solveFemUpdateJacobian,
                                                          std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );
    M_backend_primal->nlSolver()->residual = std::bind( &CRBModel<TruthModelType>::solveFemUpdateResidual,
                                                          std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );
    //M_backend_primal->nlSolver()->setType( TRUST_REGION );

    auto solution = this->functionSpace()->element();
    if( M_usePrimalPc )
    {
        M_preconditioner_primal->setMatrix( J );
        M_backend_primal->nlSolve(_jacobian=J, _solution=solution, _residual=R,_prec=M_preconditioner_primal);
    }
    else
        M_backend_primal->nlSolve(_jacobian=J, _solution=solution, _residual=R);
    return solution;
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::solveFemUpdateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype & J , const parameter_type & mu)
{
    J->zero();

    beta_vector_type betaJqm;
    this->updateJacobian( X, M_Jqm );
    boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = this->computeBetaQm( X , mu , 0 );

    for ( size_type q = 0; q < this->Qa(); ++q )
    {
        for(int m=0; m<this->mMaxA(q); m++)
        {
            J->addMatrix( betaJqm[q][m], M_Jqm[q][m] );
        }
    }
}

template<typename TruthModelType>
void
CRBModel<TruthModelType>::solveFemUpdateResidual( const vector_ptrtype& X, vector_ptrtype& R , const parameter_type & mu)
{
    R->zero();
    std::vector< beta_vector_type > betaRqm;
    this->updateResidual( X, M_Rqm );
    boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = this->computeBetaQm( X , mu , 0 );

    for ( size_type q = 0; q < this->Ql( 0 ); ++q )
    {
        for(int m=0; m<this->mMaxF(0,q); m++)
            R->add( betaRqm[0][q][m] , *M_Rqm[0][q][m] );
    }
}

template<typename TruthModelType>
bool
CRBModel<TruthModelType>::updateJacobian( vector_ptrtype const& X, std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm )
{
    element_type u = this->functionSpace()->element();
    u = *X;
    return this->updateJacobian( u, Jqm );
}
template<typename TruthModelType>
bool
CRBModel<TruthModelType>::updateJacobian( element_type const& X, std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm )
{
    return M_model->updateJacobian( X, Jqm );
}

template<typename TruthModelType>
bool
CRBModel<TruthModelType>::updateResidual( vector_ptrtype const& X, std::vector< std::vector< std::vector<vector_ptrtype> > >& Rqm)
{
    element_type u = this->functionSpace()->element();
    u = *X;
    return this->updateResidual( u, Rqm );
}
template<typename TruthModelType>
bool
CRBModel<TruthModelType>::updateResidual( element_type const& X, std::vector< std::vector< std::vector<vector_ptrtype> > >& Rqm)
{
    return M_model->updateResidual( X, Rqm );
}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemDualMonolithicFormulation( parameter_type const& mu )
{
    int output_index = M_outputIndex;

    auto Xh= this->functionSpace();

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype A,Adu;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    vector_ptrtype Rhs( M_backend->newVector( Xh ) );
    element_ptrtype InitialGuess = Xh->elementPtr();
    auto dual_initial_field = Xh->elementPtr();

    auto udu = Xh->element();

    double time_initial;
    double time_step;
    double time_final;

    if ( this->isSteady() )
    {
        time_initial=0;
        time_step = 1e30;
        time_final = 1e30;
        // !!
        //some stuff needs to be done here
        //to deal with non linear problems
        //and have an initial guess
        // !!
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

    double bdf_coeff ;
    auto vec_bdf_poly = M_backend->newVector( Xh );

    if ( this->isSteady() )
        udu.zero() ;
    else
    {
        boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );
        *Rhs=*F[output_index];
        M_preconditioner_dual->setMatrix( M );
        M_backend_dual->solve( _matrix=M, _solution=dual_initial_field, _rhs=Rhs, _prec=M_preconditioner_dual );
        udu=*dual_initial_field;
    }

    for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
    {
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;

        if ( is_linear )
            boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );
        else
            boost::tie(M, A, F) = this->computeMonolithicFormulationU( mu , udu );

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

        if( M_useSymmMat )
            Adu = A;
        else
            A->transpose( Adu );

        M_preconditioner_dual->setMatrix( Adu );

        M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs , _prec=M_preconditioner_dual);

        mybdf->shiftRight(udu);
    }

    return udu;

}

template<typename TruthModelType>
typename CRBModel<TruthModelType>::element_type
CRBModel<TruthModelType>::solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu )
{
    int output_index = M_outputIndex;

    auto Xh= this->functionSpace();
    auto udu = Xh->element();
    auto uold = Xh->element();
    sparse_matrix_ptrtype A, Adu;
    std::vector<vector_ptrtype> F;
    vector_ptrtype Rhs = M_backend_dual->newVector( Xh );

    int max_fixedpoint_iterations  = M_fixedpointMaxIt;
    double increment_fixedpoint_tol  = M_fixedpointIncrTol;
    bool fixedPointVerbose = M_fixedpointVerbose;

    double norm=0;
    int iter=0;

    if ( this->isSteady() )
    {
        if ( is_linear )
        {
            boost::tie(boost::tuples::ignore, A, F) = this->update( mu , 0 );

            *Rhs = *F[output_index];
            Rhs->scale( -1 );

            if( M_useSymmMat )
                Adu = A;
            else
            {
                Adu = M_backend_dual->newMatrix(_test=Xh,_trial=Xh);
                A->transpose( Adu );
            }


            //uold = udu;
            M_preconditioner_dual->setMatrix( Adu );
            M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs , _prec=M_preconditioner_dual);
        }
        else
        {
            Adu = M_backend_dual->newMatrix(_test=Xh,_trial=Xh);
            bool useAitkenRelaxation = M_fixedpointUseAitken;
            auto residual = Xh->element();
            auto aitkenRelax = aitken( _space=Xh );
            aitkenRelax.initialize( _residual=residual, _currentElt=udu );
            aitkenRelax.restart();
            bool fixPointIsFinished = false;
            do {
                uold = udu;

                boost::tie(boost::tuples::ignore, A, F) = this->update( mu , udu );
                //Adu = A;
                A->transpose( Adu );
                *Rhs = *F[output_index];
                Rhs->scale( -1 );

                M_preconditioner_dual->setMatrix( Adu );
                M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs , _prec=M_preconditioner_dual );

                iter++;
                if ( !useAitkenRelaxation )
                {
                    norm = this->computeNormL2( uold , udu );
                    fixPointIsFinished = norm < increment_fixedpoint_tol || iter >= max_fixedpoint_iterations;
                    if ( this->worldComm().isMasterRank() && fixedPointVerbose )
                        std::cout << "[solveFemDualUsingAffineDecompositionFixedPoint] iteration " << iter << ", increment_norm = " <<  norm << "\n";
                }
                else
                {
                    residual = udu-uold;
                    aitkenRelax.apply2(_newElt=udu,_residual=residual, _currentElt=udu );
                    if ( this->worldComm().isMasterRank() && fixedPointVerbose )
                        aitkenRelax.printInfo();
                    ++aitkenRelax;
                    if ( aitkenRelax.isFinished() )
                        fixPointIsFinished=true;
                }
            } while( !fixPointIsFinished );//norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
        }
        return udu;
    }

    bdf_ptrtype mybdf;
    mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
    sparse_matrix_ptrtype M;
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


    for( mybdf->start(udu); !mybdf->isFinished(); mybdf->next() )
    {
        iter=0;
        bdf_coeff = mybdf->polyDerivCoefficient( 0 );
        auto bdf_poly = mybdf->polyDeriv();
        *vec_bdf_poly = bdf_poly;
        do {
            if( is_linear )
                boost::tie(M, A, F) = this->update( mu , mybdf->time() );
            else
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

            if( M_useSymmMat )
                Adu = A;
            else
                A->transpose( Adu );

            uold = udu;
            M_preconditioner_dual->setMatrix( Adu );
            M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs , _prec=M_preconditioner_dual);

            if( M_useLinearModel )
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
