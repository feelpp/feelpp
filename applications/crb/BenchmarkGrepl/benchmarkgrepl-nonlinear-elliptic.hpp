/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2014-01-19

  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file benchmarkgrepl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \author Cecile Daversin <daversin@math.unistra.fr>

   date 2014-01-19
 */
#ifndef FEELPP_BENCHMARKGREPLNONLINEARELLIPTIC_HPP
#define FEELPP_BENCHMARKGREPLNONLINEARELLIPTIC_HPP 1

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

#include <applications/crb/BenchmarkGrepl/benchmarkgrepl-options.hpp>

namespace Feel
{

namespace BenchmarkGreplNonlinearElliptic_Definition
{

class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int Order, int Dim>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef bases<Lagrange<Order + 2, Scalar> > basis_type_eimg;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef FunctionSpace<mesh_type, basis_type_eimg, value_type> space_type_eimg;

    typedef typename space_type::element_type element_type;

};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef typename FunctionSpaceDefinition::space_type_eimg space_type_eimg;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type_eimg, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type_eimg, space_type , parameterspace_type> fund_type;

};

} // namespace BenchmarkGreplNonlinearElliptic_Definition

/**
 * \class BenchmarkGreplNonlinearElliptic
 * \brief brief description
 *
 * This is from the paper
 * EFFICIENT REDUCED-BASIS TREATMENT OF NONAFFINE
 * AND NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS
 *  authors :
 * Martin A. Grepl, Yvon Maday, Ngoc C. Nguyen and Anthony T. Patera
 * ESAIM: Mathematical Modelling and Numerical Analysis
 * --
 * 5. Nonaffine monotonic elliptic equations
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
template<int Order, int Dim>
class BenchmarkGreplNonlinearElliptic :
        public ModelCrbBase< BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                             BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>,
                             NonLinear,
                             BenchmarkGreplNonlinearElliptic_Definition::EimDefinition<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                                                                                       BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim> > >
{
public:

    typedef ModelCrbBase<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                         BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>,
                         NonLinear,
                         BenchmarkGreplNonlinearElliptic_Definition::EimDefinition<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                                                                                   BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim> > > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    typedef double value_type;

    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>::space_type space_type;
    typedef typename boost::shared_ptr<space_type> space_ptrtype;
    typedef typename BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>::space_type_eimg space_type_eimg;
    typedef typename boost::shared_ptr<space_type_eimg> space_ptrtype_eimg;

    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::monolithic_type monolithic_type;
    typedef typename super_type::eim_interpolation_error_type eim_interpolation_error_type;

    typedef std::vector< std::vector<sparse_matrix_ptrtype> > vector_sparse_matrix;

    using super_type::computeBetaQm;
    typedef boost::tuple<beta_vector_type,  std::vector<beta_vector_type> > beta_type;

    typedef BenchmarkGreplNonlinearElliptic<Order,Dim> self_type;

    BenchmarkGreplNonlinearElliptic()
        :
        super_type( "BenchMarkGreplNonlinearElliptic" + std::to_string(Order) + "-" + std::to_string(Dim) + "D" ),
        M_use_newton( boption(_name="crb.use-newton") ),
        M_useSerErrorEstimation( boption(_name="ser.error-estimation") )
        {}

    //! initialization of the model
    void initModel();
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir );

    //@}

    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu )
    {

        this->M_betaInitialGuess.resize( 1 );
        this->M_betaInitialGuess[0].resize( 1 );

        this->M_betaInitialGuess[0][0] = 1;

        return this->M_betaInitialGuess;
    }


    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    beta_type
    computeBetaQm( element_type const& T,parameter_type const& mu )
    {
        auto eim_g = this->scalarContinuousEim()[0];
        vectorN_type beta_g = eim_g->beta( mu , T );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
    {
        auto eim_g = this->scalarContinuousEim()[0];
        vectorN_type beta_g = eim_g->beta( mu , urb );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computeBetaQm( parameter_type const& mu )
    {
        auto eim_g = this->scalarContinuousEim()[0];
        vectorN_type beta_g = eim_g->beta( mu );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);

        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computePicardBetaQm( element_type const& T,parameter_type const& mu )
    {
        auto eim_g = this->scalarContinuousEim()[0];
        vectorN_type beta_g = eim_g->beta( mu , T );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    beta_type
    computePicardBetaQm( parameter_type const& mu )
    {
        auto eim_g = this->scalarContinuousEim()[0];
        vectorN_type beta_g = eim_g->beta( mu );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    void fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu)
    {
        auto eim_g = this->scalarContinuousEim()[0];
        int M = eim_g->mMax();

        auto beta_g=*betas[0];

        if( M_use_newton )
        {
            this->M_betaJqm[0][0] = 1;
            this->M_betaJqm[1].resize(M); //needed if cobuild
            for(int m=0; m<M; m++)
            {
                this->M_betaJqm[1][m] = mu(1)*beta_g(m);
            }
            this->M_betaJqm[2][0] = mu(0);

            this->M_betaRqm[0][0][0] = this->computeBetaInitialGuess( mu )[0][0];
            this->M_betaRqm[0][1].resize(M);
            for(int m=0; m<M; m++)
            {
                this->M_betaRqm[0][1][m] = beta_g(m);
            }
            this->M_betaRqm[0][2][0] = -100;
            //output
            this->M_betaRqm[1][0][0] = 1;
        }
        //else
        if( !M_use_newton || M_useSerErrorEstimation )
        {
            this->M_betaAqm[0][0]=1;
            this->M_betaFqm[0][0].resize(M);
            for(int m=0; m<M; m++)
            {
                this->M_betaFqm[0][0][m]=-beta_g(m) ;
            }
            this->M_betaFqm[0][1][0]=100;
            this->M_betaFqm[1][0][0]=1;
        }

    }

    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu , double time=1e30 )
    {
        beta_vector_type beta;
        beta.resize(1);
        beta[0].resize(1);
        beta[0][0]=1;
        return beta;
    }

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();
    affine_decomposition_type computePicardAffineDecomposition();
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA();

    void assemble();

    std::vector< std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition( );


    //@}


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false );

    gmsh_ptrtype createGeo( double hsize );

    void assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void assembleJacobianWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & Jqm);
    void assembleFunctionalWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & RF_Aqm,
                                                   std::vector< std::vector<std::vector<vector_ptrtype> > >& RF_Fqm);
    bool updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void updateResidualMonolithic(vector_ptrtype const& X, vector_ptrtype & R, parameter_type const& mu);
    void updateJacobianMonolithic(vector_ptrtype const& X, sparse_matrix_ptrtype & J, parameter_type const& mu);
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution );

    element_type solve( parameter_type const& mu );

private:

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_ptrtype pT;
    parameter_type M_mu;

    sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    bool M_use_newton;
    bool M_useSerErrorEstimation;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

};


template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    this->M_betaAqm.resize( 1 );
    this->M_betaAqm[0].resize( 1 );
    this->M_betaFqm.resize(2);
    this->M_betaFqm[0].resize( 2 );
    //this->M_betaFqm[0][0].resize( M );
    this->M_betaFqm[0][1].resize( 1 );
    this->M_betaFqm[1].resize( 1 );
    this->M_betaFqm[1][0].resize( 1 );

    boost::shared_ptr<space_type_eimg> Xh_eimg;
    if ( !pT )
        pT.reset( new element_type );

    auto const& ptreeEim = ptree.get_child( "eim" );
    auto const& ptreeEimg = ptreeEim.get_child( "eim_g" );
    std::string dbnameEimg = ptreeEimg.template get<std::string>( "database-filename" );

    auto eim_g = eim( _model=boost::dynamic_pointer_cast< BenchmarkGreplNonlinearElliptic<Order,Dim> >( this->shared_from_this() ),
                      _element=*pT,
                      _space=Xh_eimg,
                      _parameter=M_mu,
                      _expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*_e1 ) - 1 ),
                      //_sampling=Pset,
                      _name="eim_g",
                      _filename=dbnameEimg,
                      _directory=dbDir );
    this->addEim( eim_g );
}
template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::initModel()
{
    std::string mshfile_name = option("mshfile").as<std::string>();

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        double hsize=doption(_name="hsize");
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name = "benchmarkgrepl",
                                             _shape = "hypercube",
                                             _dim = Dim,
                                             _h=hsize,
                                             _xmin=0,_xmax=1,
                                             _ymin=0,_ymax=1 ) );
    }
    else
    {
        mesh = loadGMSHMesh( _mesh=new mesh_type,
                             _filename=option("mshfile").as<std::string>(),
                             _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    }


    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    space_ptrtype_eimg Xh_eimg = space_type_eimg::New( mesh );
    this->setFunctionSpaces( Xh );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
    }

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );
    static auto M_U = *pT;

    auto mu_min = this->Dmu->element();
    mu_min <<  0.01, 0.01;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << 10, 10;
    this->Dmu->setMax( mu_max );

    M_mu = this->Dmu->element();

    auto Pset = this->Dmu->sampling();
    //specify how many elements we take in each direction
    std::vector<size_type> N(2);
    int Ne = ioption(_name="trainset-eim-size");
    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();

    std::ifstream file ( supersamplingname );

    //40 elements in each direction
    N[0]=Ne; N[1]=Ne;

    //interpolation points are located on different proc
    //so we can't distribute parameters on different proc as in crb case
    //else for a given mu we are not able to evaluate g at a node wich
    //is not on the same proc than mu (so it leads to wrong results !)
    bool all_proc_same_sampling=true;

    if( ! file )
    {
        Pset->equidistributeProduct( N , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto eim_g = eim( _model=boost::dynamic_pointer_cast< BenchmarkGreplNonlinearElliptic<Order,Dim> >( this->shared_from_this() ),
                      _element=*pT,
                      _space=Xh_eimg,
                      _parameter=M_mu,
                      //_expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*idv(*pT) ) - 1 ),
                      _expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*_e1 ) - 1 ),
                      _sampling=Pset,
                      _name="eim_g" );

    this->addEim( eim_g );

    int M = this->scalarContinuousEim()[0]->mMax();

    //resize data structure
    if( M_use_newton )
    {
        this->M_betaJqm.resize( 3 );
        this->M_betaJqm[0].resize( 1 );
        this->M_betaJqm[1].resize( M );
        this->M_betaJqm[2].resize( 1 );
        this->M_betaRqm.resize( 2 );
        this->M_betaRqm[0].resize( 3 );
        this->M_betaRqm[0][0].resize( 1 );
        this->M_betaRqm[0][1].resize( M );
        this->M_betaRqm[0][2].resize( 1 );
        this->M_betaRqm[1].resize( 1 );
        this->M_betaRqm[1][0].resize( 1 );

        this->M_Jqm.resize( 3 );
        this->M_Jqm[0].resize( 1 );
        this->M_Jqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        this->M_Jqm[1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        }
        this->M_Jqm[2].resize( 1 );
        this->M_Jqm[2][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Rqm.resize( 2 );
        this->M_Rqm[0].resize( 3 );
        this->M_Rqm[1].resize( 1 );
        this->M_Rqm[0][0].resize(1);
        this->M_Rqm[0][0][0] = backend()->newVector( this->Xh );
        this->M_Rqm[0][1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        }
        this->M_Rqm[0][2].resize(1);
        this->M_Rqm[0][2][0] = backend()->newVector( this->Xh );
        this->M_Rqm[1][0].resize( 1 );
        this->M_Rqm[1][0][0]= backend()->newVector( this->Xh );
    }

    if( !M_use_newton || M_useSerErrorEstimation )
    {
        this->M_betaAqm.resize( 1 );
        this->M_betaAqm[0].resize( 1 );
        this->M_betaFqm.resize(2);
        this->M_betaFqm[0].resize( 2 );
        this->M_betaFqm[0][0].resize( M );
        this->M_betaFqm[0][1].resize( 1 );
        this->M_betaFqm[1].resize( 1 );
        this->M_betaFqm[1][0].resize( 1 );

        this->M_Aqm.resize( 1 );
        this->M_Aqm[0].resize( 1 );
        this->M_Aqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Fqm.resize( 2 );
        this->M_Fqm[0].resize( 2 );
        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
        }
        this->M_Fqm[0][1].resize(1);
        this->M_Fqm[0][1][0]=backend()->newVector( this->Xh );
        this->M_Fqm[1].resize(1);
        this->M_Fqm[1][0].resize(1);
        this->M_Fqm[1][0][0]=backend()->newVector( this->Xh );
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();

    assemble();

} // BenchmarkGreplNonlinearElliptic::init

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::assembleJacobianWithAffineDecomposition( std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto u = Xh->element(); //trial
    auto eim_g = this->scalarContinuousEim()[0];
    int M = eim_g->mMax();
    double gamma = doption(_name="gamma");

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idt(u)*id(v)/hFace()
                   - (gradt(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idt(u) );

    Jqm[0][0]->close();
    Jqm[1].resize(M);
    for(int m=0; m<M; m++)
    {
        this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        form2( _test=Xh, _trial=Xh, _matrix=Jqm[1][m] ) =
            integrate( _range = elements(mesh),
                       _expr = idv(eim_g->q(m))*idt(u)*id(v) );
        Jqm[1][m]->close();
    }

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[2][0] ) =
        integrate( _range= elements( mesh ), _expr = idt(u)*id(v) );
    Jqm[2][0]->close();

}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto eim_g = this->scalarContinuousEim()[0];
    int M = eim_g->mMax();

    auto u = Xh->element();
    u = *this->M_InitialGuess[0][0];
    double gamma = doption(_name="gamma");

    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );
    Rqm[0][0][0]->close();

    Rqm[0][1].resize(M);

    for(int m=0; m<M; m++)
    {
        this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        form1( _test=Xh, _vector=Rqm[0][1][m] ) =
            integrate( _range= elements( mesh ), _expr=( idv(eim_g->q(m))*id(v) ) );
        Rqm[0][1][m]->close();
    }

    if( Dim == 2 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=Rqm[1][0][0] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    Rqm[0][2][0]->close();
    Rqm[1][0][0]->close();
}

template <int Order, int Dim>
bool
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto Xh = this->Xh;
    auto u = Xh->element();
    auto v = Xh->element(); //test
    double gamma = option(_name="gamma").template as<double>();
    //u = *X;
    u = X;
    Rqm[0][0][0] = backend()->newVector( this->Xh );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    Rqm[0][0][0]->close();
    // update associated coefficient
    // this->M_betaRqm[0][0][0] = 1;
    return true;
}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateJacobianMonolithic( vector_ptrtype const& X,
                                                                  sparse_matrix_ptrtype & J,
                                                                  parameter_type const& mu )
{
    //Here we don't use EIM approximation and so no affine decomposition
    auto Xh = this->Xh;
    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    double gamma = doption(_name="gamma");
    auto g = exp( mu(1)*idv(u) );

    J->zero();
    form2( _test=Xh, _trial=Xh, _matrix=J ) =
        integrate( _range= elements( mesh ),_expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = boundaryfaces( mesh ),
                                                          _expr = gamma*idt(u)*id(v)/hFace()
                                                          - (gradt(u)*vf::N())*id(v)
                                                          - (grad(v)*vf::N())*idt(u) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range = elements(mesh), _expr = mu(0)*exp( mu(1)*idv(u) )*idt(u)*id(v) );

    J->close();

}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateResidualMonolithic(vector_ptrtype const& X,
                                                                 vector_ptrtype & R,
                                                                 parameter_type const& mu )
{
    auto Xh = this->Xh;
    double gamma = doption(_name="gamma");

    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    auto g = exp( mu(1)*idv(u) );

    R->zero();
    form1( _test=Xh, _vector=R ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= mu(0)/mu(1)*( exp( mu(1)*idv(u) ) * id(v) ) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= -mu(0)/mu(1)*( id(v)) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    R->close();

}

template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::assemble()
{
    auto Xh = this->Xh;
    auto u = Xh->element();

    if( M_use_newton )
    {
        assembleJacobianWithAffineDecomposition( this->M_Jqm );
        assembleResidualWithAffineDecomposition( this->M_Rqm );
    }

    if( !M_use_newton || M_useSerErrorEstimation )
    {
        auto v = Xh->element();
        auto eim_g = this->scalarContinuousEim()[0];
        int M = eim_g->mMax();
        double gamma = doption(_name="gamma");

        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) =
            integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) += integrate( _range = boundaryfaces( mesh ),
                                                                              _expr = gamma*idt(u)*id(v)/hFace()
                                                                              - (gradt(u)*vf::N())*id(v)
                                                                              - (grad(v)*vf::N())*idt(u) );
        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
            form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) =
                integrate( _range= elements( mesh ), _expr=( idv(eim_g->q(m))*id(v) ) );
            this->M_Fqm[0][0][m]->close();
        }

        if( Dim == 2 )
            form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
                integrate( _range= elements( mesh ),_expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
        else if ( Dim == 3 )
            form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
                integrate( _range= elements( mesh ),_expr=sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );
        this->M_Fqm[0][1][0]->close();

        form1( _test=Xh, _vector=this->M_Fqm[1][0][0] ) =
            integrate( _range= elements( mesh ),
                       _expr=id(v) );
        this->M_Fqm[1][0][0]->close();

    }//no newton

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::monolithic_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution)
{
    auto Xh = this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    M_monoA = backend()->newMatrix( Xh, Xh );
    M_monoF.resize(2);
    M_monoF[0] = backend()->newVector( this->functionSpace() );
    M_monoF[1] = backend()->newVector( Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) += integrate( _range = boundaryfaces( mesh ),
                                                                _expr = gamma*idt(u)*id(v)/hFace()
                                                                - (gradt(u)*vf::N())*id(v)
                                                                - (grad(v)*vf::N())*idt(u) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if ( Dim == 3 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=M_monoF[1] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    M_monoF[0]->close();
    M_monoF[1]->close();
    M_monoA->close();

    return boost::make_tuple( M_monoA, M_monoF );

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::element_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::solve( parameter_type const& mu )
{
    auto Xh=this->Xh;
    sparse_matrix_ptrtype J = backend()->newMatrix( Xh, Xh);
    vector_ptrtype R = backend()->newVector( Xh );
    backend()->nlSolver()->jacobian = boost::bind( &self_type::updateJacobianMonolithic,
                                                   boost::ref( *this ), _1, _2, mu );
    backend()->nlSolver()->residual = boost::bind( &self_type::updateResidualMonolithic,
                                                   boost::ref( *this ), _1, _2, mu );

    auto solution = Xh->element();
    backend()->nlSolve(_jacobian=J, _solution=solution, _residual=R);

    return solution;

}


template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::vector_sparse_matrix
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeLinearDecompositionA()
{
    auto Xh=this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    this->M_linearAqm.resize(1);
    this->M_linearAqm[0].resize(1);
    this->M_linearAqm[0][0] = backend()->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( boundaryfaces( mesh ), gamma*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );

    this->M_linearAqm[0][0]->close();
    return this->M_linearAqm;

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::affine_decomposition_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeAffineDecomposition()
{
    if( M_use_newton )
        return boost::make_tuple( this->M_Jqm, this->M_Rqm );
    else
        return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::affine_decomposition_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computePicardAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order, int Dim>
std::vector< std::vector< typename BenchmarkGreplNonlinearElliptic<Order,Dim>::element_ptrtype > >
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeInitialGuessAffineDecomposition( )
{
    return M_InitialGuess;
}


template<int Order, int Dim>
double BenchmarkGreplNonlinearElliptic<Order,Dim>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve )
{
    //CHECK( need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;
    if ( output_index==0 )
    {
        if( Dim == 2 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * idv( solution ) ).evaluate()(0,0);
        else if ( Dim == 3 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * idv( solution ) ).evaluate()(0,0);
    }
    if( output_index==1 )
    {
        s = integrate( elements( mesh ), idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}

}

#endif /* __BenchmarkGreplNonlinearElliptic_H */
