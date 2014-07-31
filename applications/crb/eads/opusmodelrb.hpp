/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-12-10

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file opusmodel.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-10
 */
#ifndef __OpusModelRB_H
#define __OpusModelRB_H 1

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

/**/
#include <opusdata.hpp>
#include <opusmodelbase.hpp>
#include <opusdefs.hpp>
#include <eads.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

/**/
namespace Feel
{
Feel::po::options_description opusModelOptions();

//    template<typename SpaceType> class OpusModelThermalRB;
//template<typename SpaceType> class OpusModelFluidPoiseuille;
//template<typename SpaceType> class OpusModelFluidPoiseuilleRB;
//template<typename SpaceType> class OpusModelFluidOseen;

/*!
  \class OpusModelRB
  \brief Opus ModelRB


  @author Christophe Prud'homme
  @see
*/


class ParameterDefinition
{
public:
    static const uint16_type ParameterSpaceDimension = 5;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int OrderU=2, int OrderP=OrderU-1, int OrderT=OrderP>
class FunctionSpaceDefinition
{
public:

    static const uint16_type Dim = 2;

    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;

    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<3>, mpl::int_<11>, mpl::int_<13> >,mpl::vector<mpl::int_<4>, mpl::int_<11>, mpl::int_<14> > > >  discontinuity_type;
    typedef bases<Lagrange<OrderT, Scalar, discontinuity_type> > temp_basis_type;

#if defined( OPUS_WITH_THERMAL_DISCONTINUITY )
    typedef FunctionSpace<mesh_type, temp_basis_type, discontinuity_type,  Periodicity<Periodic<> > > temp_functionspace_type;
#else
    typedef FunctionSpace<mesh_type, temp_basis_type, Periodicity<Periodic<> > > temp_functionspace_type;
#endif

    typedef temp_functionspace_type space_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};

template<int OrderU=2, int OrderP=OrderU-1, int OrderT=OrderP>
class OpusModelRB : public OpusModelBase, public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<OrderU,OrderP,OrderT> >,
                    public boost::enable_shared_from_this< OpusModelRB<OrderU,OrderP,OrderT> >
{
    typedef OpusModelBase super;
public:

    typedef ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition<OrderU,OrderP,OrderT> > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    typedef OpusModelRB< OrderU, OrderP, OrderT > model_type;

    /** @name Typedefs
     */
    //@{
    static const uint16_type Dim = 2;

    // #if FEELPP_GNUC_AT_LEAST(4,6)
#if 0
    static constexpr double kmin ;
    static constexpr double kmax ;

    static constexpr double rmin ;
    static constexpr double rmax ;

    static constexpr double Qmin ;
    static constexpr double Qmax ;
#else
    static const double kmin ;
    static const double kmax ;

    static const double rmin ;
    static const double rmax ;

    static const double Qmin ;
    static const double Qmax ;
#endif
    static const uint16_type ParameterSpaceDimension = 5;

    //@}
    /** @name Typedefs
     */
    //@{


    typedef double value_type;


    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Simplex<1,1,2> entity12_type;
    typedef Mesh<entity12_type> mesh12_type;
    typedef boost::shared_ptr<mesh12_type> mesh12_ptrtype;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    /* temperature */
    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<3>, mpl::int_<11>, mpl::int_<13> >,mpl::vector<mpl::int_<4>, mpl::int_<11>, mpl::int_<14> > > >  discontinuity_type;
    typedef bases<Lagrange<OrderT, Scalar, discontinuity_type> > temp_basis_type;
    //typedef bases<Lagrange<OrderT, Scalar> > temp_basis_type;
    typedef bases<Lagrange<OrderT, Vectorial> > grad_temp_basis_type;
    //typedef Periodic<1,2,value_type> periodic_type;
    typedef Periodic<> periodic_type;


    typedef temp_basis_type basis_type;

#if defined( OPUS_WITH_THERMAL_DISCONTINUITY )
    typedef FunctionSpace<mesh_type, temp_basis_type, discontinuity_type,  Periodicity<Periodic<> > > temp_functionspace_type;
    typedef ReducedBasisSpace< model_type, mesh_type, temp_basis_type, discontinuity_type,  Periodicity<Periodic<> > > temp_rbfunctionspace_type;
#else
    typedef FunctionSpace<mesh_type, temp_basis_type, Periodicity<Periodic<> > > temp_functionspace_type;
    typedef ReducedBasisSpace< model_type, mesh_type, temp_basis_type, Periodicity<Periodic<> > > temp_rbfunctionspace_type;
#endif
    typedef boost::shared_ptr<temp_functionspace_type> temp_functionspace_ptrtype;
    typedef boost::shared_ptr<temp_rbfunctionspace_type> temp_rbfunctionspace_ptrtype;
    typedef typename temp_functionspace_type::element_type temp_element_type;
    typedef boost::shared_ptr<temp_element_type> temp_element_ptrtype;

    typedef temp_functionspace_type functionspace_type;
    typedef temp_functionspace_ptrtype functionspace_ptrtype;
    typedef temp_rbfunctionspace_type rbfunctionspace_type;
    typedef temp_rbfunctionspace_ptrtype rbfunctionspace_ptrtype;
    typedef temp_element_ptrtype element_ptrtype;
    typedef temp_element_type element_type;
    typedef temp_functionspace_type space_type;

    typedef FunctionSpace<mesh_type, grad_temp_basis_type> grad_temp_functionspace_type;
    typedef boost::shared_ptr<grad_temp_functionspace_type> grad_temp_functionspace_ptrtype;
    typedef typename grad_temp_functionspace_type::element_type grad_temp_element_type;
    typedef boost::shared_ptr<grad_temp_element_type> grad_temp_element_ptrtype;

    /* 1D profil */
    typedef bases<Lagrange<1, Scalar> > p1_basis_type;
    typedef FunctionSpace<mesh12_type, p1_basis_type, value_type> p1_functionspace_type;
    typedef boost::shared_ptr<p1_functionspace_type> p1_functionspace_ptrtype;
    typedef typename p1_functionspace_type::element_type p1_element_type;
    typedef boost::shared_ptr<p1_element_type> p1_element_ptrtype;

    /* P0 */
    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar, Discontinuous > > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;
    typedef boost::shared_ptr<p0_element_type> p0_element_ptrtype;

    //typedef Feel::OpusModelThermalRB<temp_functionspace_type> thermal_operator_type;
    //typedef boost::shared_ptr<thermal_operator_type> thermal_operator_ptrtype;

    /* velocity */
    typedef bases<Lagrange<OrderU, Vectorial> > velocity_basis_type;
    typedef FunctionSpace<mesh_type, velocity_basis_type, value_type> velocity_functionspace_type;
    typedef boost::shared_ptr<velocity_functionspace_type> velocity_functionspace_ptrtype;
    typedef typename velocity_functionspace_type::element_type velocity_element_type;
    typedef boost::shared_ptr<velocity_element_type> velocity_element_ptrtype;

    /* pressure */
    typedef bases<Lagrange<OrderP, Scalar> > pressure_basis_type;
    typedef FunctionSpace<mesh_type, pressure_basis_type, value_type> pressure_functionspace_type;
    typedef boost::shared_ptr<pressure_functionspace_type> pressure_functionspace_ptrtype;
    typedef typename pressure_functionspace_type::element_type pressure_element_type;
    typedef boost::shared_ptr<pressure_element_type> pressure_element_ptrtype;
    /* fluid */
    typedef bases<Lagrange<OrderU, Vectorial>,Lagrange<OrderP, Scalar> > fluid_basis_type;

    typedef FunctionSpace<mesh_type, fluid_basis_type, value_type> fluid_functionspace_type;
    typedef boost::shared_ptr<fluid_functionspace_type> fluid_functionspace_ptrtype;
    typedef typename fluid_functionspace_type::element_type fluid_element_type;
    typedef boost::shared_ptr<fluid_element_type> fluid_element_ptrtype;
    typedef typename fluid_element_type::template sub_element<0>::type fluid_element_0_type;
    typedef typename fluid_element_type::template sub_element<1>::type fluid_element_1_type;


    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;


    /* time */
    typedef Bdf<temp_functionspace_type>  temp_bdf_type;
    typedef boost::shared_ptr<temp_bdf_type> temp_bdf_ptrtype;


    typedef std::vector< std::vector< double > > beta_vector_type;
    typedef boost::tuple<beta_vector_type, beta_vector_type, std::vector<beta_vector_type> > beta_vectors_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
        > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    OpusModelRB();
    OpusModelRB( po::variables_map const& vm );
    OpusModelRB( OpusModelRB const & );
    ~OpusModelRB();

    void initModel();
    void initParametrization();
    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const;

    // \return the number of terms in affine decomposition of mass matrix
    int Qm() const;


    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const;

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const;

    /**
     * \param q the qth term of the decomposition
     * \return the number of terms needed by the EIM
     * as this model don't use the EIMit always return 1
     */
    int mMaxA( int q );
    int mMaxM( int q );
    int mMaxF( int output_index, int q );


    /**
     * \brief Returns the function space
     */
    functionspace_ptrtype functionSpace()
    {
        return M_Th;
    }

    /**
     * \brief Returns the function space
     */
    rbfunctionspace_ptrtype rBFunctionSpace()
    {
        return M_RbTh;
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }

    parameter_type refParameter()
    {
        return M_Dmu->min();
    }


    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    beta_vectors_type
    computeBetaQm( parameter_type const& mu , double time=0 );

    beta_vectors_type
    computeBetaQm( element_type const &T,parameter_type const& mu , double time=0 );

    /**
     * \brief return the coefficient vector
     */
    beta_vector_type const& betaAqm() const
    {
        return M_betaAqm;
    }

    /**
     * \brief return the coefficient vector
     */
    beta_vector_type const& betaMqm() const
    {
        return M_betaMqm;
    }


    /**
     * \brief return the coefficient vector
     */
    std::vector<beta_vector_type> const& betaL() const
    {
        return M_betaL;
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type betaAqm( int q, int m ) const
    {
        return M_betaAqm[q][m];
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type betaMqm( int q, int m ) const
    {
        return M_betaMqm[q][m];
    }

    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type betaL( int l, int q, int m ) const
    {
        return M_betaL[l][q][m];
    }

    /**
     * return true if initialized (init() was called), false otherwise
     */
    bool isInitialized() const
    {
        return M_is_initialized;
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
        M_meshSize = s;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;

    /**
     * create a new vector
     * \return the newly created vector
     */
    vector_ptrtype newVector() const;

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition()
    {
        return boost::make_tuple( M_Mqm, M_Aqm, M_L );
    }

    std::vector< std::vector< element_ptrtype > > computeInitialGuessAffineDecomposition()
    {
        std::vector< std::vector<element_ptrtype> > q;
        q.resize(1);
        q[0].resize(1);
        element_ptrtype elt ( new element_type ( M_Th ) );
        q[0][0] = elt;
        return q;
    }


    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& T );

    /**
     * solve the model for a given parameter \p mu and \p L as right hand side
     * \param transpose if true solve the transposed(dual) problem
     */
    void solve( parameter_type const& mu, element_ptrtype& u, vector_ptrtype const& L, bool transpose = false );

    /**
     * solve for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );

    /**
     * H1 scalar product
     */
    sparse_matrix_ptrtype energyMatrix();

    /**
     * update the PDE system with respect to \param mu
     */
    void update( parameter_type const& mu, const double time=0 );

    void run();
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    //void exportResults( double time, temp_element_type& T, fluid_element_type& U, bool force_export = false );
    void exportResults( double time, temp_element_type& T , parameter_type const& mu );


    std::vector<double> sigmaQ( double k,double r, double Q );

    /**
     * returns the scalar product of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );

    /**
     * returns the scalar product used for POD of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product used for POD of the vector x and vector y
     */
    double scalarProductForPod( vector_type const& x, vector_type const& y );


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);


    /**
     * return true if we want to be in a steady state
     */
    bool isSteady()
    {
        return M_is_steady;
    }

    /**
     * return value of the time step
     */
    double timeStep()
    {
        return M_temp_bdf->timeStep();
    }

    /**
     * return value of the time final
     */
    double timeFinal()
    {
        return M_temp_bdf->timeFinal();
    }

    /**
     * return value of the time initial
     */
    double timeInitial()
    {
        return M_temp_bdf->timeInitial();
    }

    /**
     * return value of the time order
     */
    double timeOrder()
    {
        return M_temp_bdf->timeOrder();
    }

    /**
     * return number of snapshots used
     */
    int computeNumberOfSnapshots()
    {
        return M_temp_bdf->timeFinal()/M_temp_bdf->timeStep();
    }

    /**
     * return initialization filed used
     */
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu ) ;

    temp_bdf_ptrtype bdfModel(){ return M_temp_bdf;}
    //@}

protected:

private:

private:


    backend_ptrtype backend;
    backend_ptrtype backendM;

    double M_meshSize;

    bool M_is_steady;

    bool M_is_initialized;


    mesh_ptrtype M_mesh;
    mesh_ptrtype M_mesh_air;
    mesh12_ptrtype M_mesh_line;
    mesh12_ptrtype M_mesh_cross_section_2;

    parameterspace_ptrtype M_Dmu;

    int Nmax;
    int Taille_rb;
    double epsilon1;
    double epsilon2;
    double tau;
    double M_dt;


    p0_space_ptrtype M_P0h;
    p0_element_ptrtype domains;
    p0_element_ptrtype rhoC;
    p0_element_ptrtype k;
    p0_element_ptrtype Q;
    p1_functionspace_ptrtype M_P1h;
    temp_functionspace_ptrtype M_Th;
    temp_rbfunctionspace_ptrtype M_RbTh;
    grad_temp_functionspace_ptrtype M_grad_Th;
    //thermal_operator_ptrtype M_thermal;
    element_ptrtype pT,pV;

    boost::shared_ptr<export_type> M_exporter;

    sparse_matrix_ptrtype D,M,Mass,Mpod;
    std::vector<vector_ptrtype> L;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_L;

    beta_vector_type M_betaAqm;
    beta_vector_type M_betaMqm;
    std::vector<beta_vector_type> M_betaL;


    temp_bdf_ptrtype M_temp_bdf;
    double M_bdf_coeff;
    element_ptrtype M_bdf_poly;
    double M_T0;//initial value of temperature

};


//#if !FEELPP_GNUC_AT_LEAST(4,6)
template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::kmin = 0.2;
template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::kmax = 150.;


template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::rmin=0.1;
template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::rmax=100;


template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::Qmin=0;
template<int OrderU, int OrderP, int OrderT>
const double
OpusModelRB<OrderU,OrderP,OrderT>::Qmax=1000000.;
//#endif

typedef OpusModelRB<2,1,2> opusmodel212_type;
}
#endif /* __OpusModelRB_H */
