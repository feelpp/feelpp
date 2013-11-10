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
   \file crbmodeltrilinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-14
 */
#ifndef __CRBModelTrilinear_H
#define __CRBModelTrilinear_H 1

#include <boost/shared_ptr.hpp>

#include <vector>


#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>


namespace Feel
{

/**
 * \class CRBModelTrilinear
 * \brief Certified Reduced Basis Model trilinear class
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
class CRBModelTrilinear
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

    CRBModelTrilinear()
        :
        M_Aqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_mode( CRBModelMode::PFEM ),
        M_model( new model_type() ),
        M_backend( backend_type::build( BACKEND_PETSC ) )
    {
        this->init();
    }

    CRBModelTrilinear( po::variables_map const& vm, CRBModelMode mode = CRBModelMode::PFEM  )
        :
        M_Aqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm( vm ),
        M_mode( mode ),
        M_model( new model_type( vm ) ),
        M_backend( backend_type::build( vm ) )
    {
        this->init();
    }

    /**
     * \param model the model to be used
     */
    CRBModelTrilinear( model_ptrtype & model )
        :
        M_Aqm(),
        M_Fqm(),
        M_is_initialized( false ),
        M_vm(),
        M_mode( CRBModelMode::PFEM ),
        M_model( model ),
        M_backend( backend_type::build( model->vm ) )
    {
        this->init();
    }

    /**
     * copy constructor
     */
    CRBModelTrilinear( CRBModelTrilinear const & o )
        :
        M_Aqm( o.M_Aqm ),
        M_Fqm( o.M_Fqm ),
        M_is_initialized( o.M_is_initialized ),
        M_vm( o.M_vm ),
        M_mode( o.M_mode ),
        M_model(  o.M_model ),
        M_backend( o.M_backend )
    {
        this->init();
    }

    //! destructor
    virtual ~CRBModelTrilinear()
    {}

    //! initialize the model (mesh, function space, operators, matrices, ...)
    FEELPP_DONT_INLINE void init()
    {
        if ( M_is_initialized )
            return;

        M_is_initialized=true;

        if( ! M_model->isInitialized() )
        {
            LOG( INFO ) << "CRBModel Model is not initialized";
            M_model->initModel();
            M_model->setInitialized( true );
        }

        M_is_initialized=true;
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRBModelTrilinear& operator=( CRBModelTrilinear const & o )
    {
        if ( this != &o )
        {
            M_Aqm = o.M_Aqm;
            M_Fqm = o.M_Fqm;
            M_model = o.M_model;
            M_backend = o.M_backend;
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
        return M_model->Qa()-M_model->QaTri();
    }


    //! return the number of \f$\mu\f$ independent terms for the trilinear form
    size_type QaTri() const
    {
        return M_model->QaTri();
    }

    sparse_matrix_ptrtype computeTrilinearForm( element_type const& xi )
    {
        return M_model->computeTrilinearForm( xi );
    }

    sparse_matrix_ptrtype jacobian( element_type const& xi )
    {
        return M_model->jacobian( xi );
    }

    vector_ptrtype residual( element_type const& xi )
    {
        return M_model->residual( xi );
    }


    vectorN_type computeStatistics ( Eigen::VectorXd vector , std::string name )
    {
        return M_model->computeStatistics( vector , name );
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
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false>, double time=0)
    {
        beta_vector_type betaAqm;
        beta_vector_type betaMqm;
        std::vector<beta_vector_type>  betaFqm;
        boost::tuple<
            beta_vector_type,
            std::vector<beta_vector_type> >
            model_beta;

        model_beta = M_model->computeBetaQm( mu );
        betaAqm = model_beta.get<0>();
        betaFqm = model_beta.get<1>();

        return boost::make_tuple( betaMqm, betaAqm, betaFqm );
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
        std::vector< std::vector<sparse_matrix_ptrtype> > mass;
        boost::tie( M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();
        // to have compatibility with SCM, we need to provide the same interface than CRBModel
        return boost::make_tuple( mass, M_Aqm, M_Fqm );
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
     * \brief the vector \c Fq[q][m] of the affine decomposition of the right hand side
     *
     * \return the vector associated with \f$F_{qm}\f$
     */
    vector_ptrtype Fqm( uint16_type l, uint16_type q, int m ) const
    {
        return M_Fqm[l][q][m];
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
     * solve \f$M u = f\f$ where \f$ M \f$ is the matrix associated to the \f$ L_2 \f$
     * norm
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f )
    {
        return M_model->l2solve( u, f );
    }


    /**
     * run the model
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve=false  )
    {
        return M_model->output( output_index, mu , u , need_to_solve);
    }



    double timeStep()
    {
        return timeStep( mpl::bool_<model_type::is_time_dependent>() );
    }
    double timeStep( mpl::bool_<true> )
    {
        double timestep;

        if ( M_model->isSteady() )
            timestep=1e30;
        else
            timestep = M_model->timeStep();
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
        if ( M_model->isSteady() )
            timefinal=1e30;
        else
            timefinal = M_model->timeFinal();
        return timefinal;
    }
    double timeFinal( mpl::bool_<false> )
    {
        return 1e30;
    }


    //only to compile
    element_type solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu )
    {
        auto zero=this->functionSpace()->element();
        return zero;
    }
    element_type solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu )
    {
        auto zero=this->functionSpace()->element();
        return zero;
    }

    element_type solveFemUsingOfflineEim( parameter_type const& mu ){};
    offline_merge_type result_offline_merge_type;
    offline_merge_type update( parameter_type const& mu,  double time=0 )  // for scm
    {
        return result_offline_merge_type ;
    }
    sparse_matrix_ptrtype const& innerProduct() const //  for the scm
    {
        return sparse_matrix ;
    }
    sparse_matrix_ptrtype const& innerProductForMassMatrix() const //  for the scm
    {
        return sparse_matrix ;
    }

    sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false ) const
    {
        return sparse_matrix;
    }
    value_type Mqm( uint16_type q, uint16_type m, element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
    {
        return 0;
    }
    size_type Qm() const
    {
        return 0;
    }
    int mMaxA( int q )
    {
        return 1;
    }
    int mMaxM( int q )
    {
        return 1;
    }
    int mMaxMF( int q )
    {
        return 1;
    }
    int mMaxF(int output_index, int q )
    {
        return 1;
    }

    //@}
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


protected:

    //! affine decomposition terms for the left hand side ( bilinear )
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;

    //! affine decomposition terms for the right hand side
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

private:

    bool M_is_initialized;

    //! variables_map
    po::variables_map M_vm;

    //! mode for CRBModelTrilinear
    CRBModelMode M_mode;

    //! model
    model_ptrtype M_model;

    backend_ptrtype M_backend;

    sparse_matrix_ptrtype sparse_matrix; // only to compile

};




template<typename TruthModelType>
void
CRBModelTrilinear<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    M_model->run( X, N, Y, P );
}



}
#endif /* __CRBModelTrilinear_H */
