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


namespace Feel
{

po::options_description
makeBenchmarkGreplNonlinearEllipticOptions()
{
    po::options_description bgoptions( "BenchmarkGreplNonlinearElliptic options" );
    bgoptions.add_options()
        ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
        ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "hsize")
        ( "trainset-eim-size", Feel::po::value<int>()->default_value( 40 ), "EIM trainset is built using a equidistributed grid 40 * 40 by default")
        ( "gamma", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
        ;
    return bgoptions;
}
AboutData
makeBenchmarkGreplNonlinearEllipticAbout( std::string const& str = "benchmarkGrepl" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Benchmark Grepl",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}


class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int Order>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    typedef typename space_type::element_type element_type;

};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;


    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fund_type;

};

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
 * 3. Nonaffine nonlinear coercive elliptic equations
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
template<int Order>
class BenchmarkGreplNonlinearElliptic :
    public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<Order> ,NonLinear, EimDefinition<ParameterDefinition, FunctionSpaceDefinition<Order> > >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition<Order>, NonLinear, EimDefinition<ParameterDefinition,FunctionSpaceDefinition<Order> > > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    typedef double value_type;

    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename FunctionSpaceDefinition<Order>::space_type space_type;

    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::monolithic_type monolithic_type;
    typedef typename super_type::eim_interpolation_error_type eim_interpolation_error_type;

    typedef std::vector< std::vector<sparse_matrix_ptrtype> > vector_sparse_matrix;

    using super_type::computeBetaQm;

    typedef BenchmarkGreplNonlinearElliptic<Order> self_type;

    //! initialization of the model
    void initModel();
    //@}

    std::string modelName()
    {
        std::ostringstream ostr;
        ostr << "BenchMarkGreplNonlinearElliptic" <<  Order;
        return ostr.str();
    }

    //\return the list of EIM objects
    virtual funs_type scalarContinuousEim() const
    {
        return M_funs;
    }

    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu )
    {

        this->M_betaInitialGuess.resize( 1 );
        this->M_betaInitialGuess[0].resize( 1 );

        this->M_betaInitialGuess[0][0] = 0;

        return this->M_betaInitialGuess;
    }


    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu )
    {
        auto eim_g = M_funs[0];
        auto eim_u = M_funs[1];
        int M = eim_g->mMax();
        int Mu = eim_u->mMax();
        vectorN_type beta_g = eim_g->beta( mu , T );
        vectorN_type beta_u = eim_u->beta( mu , T );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);
        betas.push_back(&beta_u);

        fillBetaQm(betas, mu);
        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    boost::tuple<beta_vector_type,  std::vector<beta_vector_type>  >
    computeBetaQm( parameter_type const& mu )
    {
        auto eim_g = M_funs[0];
        auto eim_u = M_funs[1];
        int M = eim_g->mMax();
        int Mu = eim_u->mMax();
        vectorN_type beta_g = eim_g->beta( mu );
        vectorN_type beta_u = eim_u->beta( mu );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);
        betas.push_back(&beta_u);

        fillBetaQm(betas, mu);
        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    void fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu)
    {
        auto eim_g = M_funs[0];
        auto eim_u = M_funs[1];
        int M = eim_g->mMax();
        int Mu = eim_u->mMax();
        auto beta_g=*betas[0];
        auto beta_u=*betas[1];
        if( M_use_newton )
        {
            this->M_betaJqm[0][0] = 1;
            for(int m=0; m<M; m++)
            {
                this->M_betaJqm[1][m] = mu(0)*beta_g(m);
            }
            //rhs
            for(int m=0; m<Mu; m++)
            {
                this->M_betaRqm[0][0][m] =  beta_u(m) ;
            }
            for(int m=0; m<M; m++)
            {
                this->M_betaRqm[0][1][m] = mu(0)/mu(1)*( beta_g(m) );
            }
            this->M_betaRqm[0][2][0] = -mu(0)/mu(1);
            this->M_betaRqm[0][3][0] = -100;
            //output
            this->M_betaRqm[1][0][0] = 1;
        }
        else
        {
            this->M_betaAqm[0][0]=1;
            for(int m=0; m<M; m++)
            {
                //this->M_betaFqm[0][0][m]=mu(0)/mu(1)*( beta(m) - 1);
                this->M_betaFqm[0][0][m]=-mu(0)/mu(1)* beta_g(m) ;
            }
            //this->M_betaFqm[0][0][0]=0 ;
            this->M_betaFqm[0][1][0]=mu(0)/mu(1) ;
            this->M_betaFqm[0][2][0]=100;
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
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA();

    void assemble();

    std::vector< std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition( );


    //@}


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    gmsh_ptrtype createGeo( double hsize );

    // void assembleResidualWithAffineDecomposition(element_type const& X,std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    // void assembleJacobianWithAffineDecomposition(element_type const& X,std::vector<std::vector<sparse_matrix_ptrtype> > & Jqm);
    void assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void assembleJacobianWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & Jqm);
    void updateResidualMonolithic(vector_ptrtype const& X, vector_ptrtype & R, parameter_type const& mu);
    void updateJacobianMonolithic(vector_ptrtype const& X, sparse_matrix_ptrtype & J, parameter_type const& mu);
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution );

    element_type solve( parameter_type const& mu );

private:

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    element_ptrtype pT;
    parameter_type M_mu;

    funs_type M_funs;

   sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    bool M_use_newton;
    int M_Qj;
    int M_Qr0;
    int M_Qr1;
    int M_Qa;
    int M_Qf0;
    int M_Qf1;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

};



template<int Order>
gmsh_ptrtype
BenchmarkGreplNonlinearElliptic<Order>::createGeo( double hsize )
{
    gmsh_ptrtype gmshp( new Gmsh );
    std::ostringstream ostr;
    double H = hsize;
    ostr <<"Point (1) = {0, 0, 0,"<< H <<"};\n"
         <<"Point (2) = {1, 0, 0,"<< H <<"};\n"
         <<"Point (3) = {1, 1, 0,"<< H <<"};\n"
         <<"Point (4) = {0, 1, 0,"<< H <<"};\n"
         <<"Line (11) = {1,2};\n"
         <<"Line (12) = {2,3};\n"
         <<"Line (13) = {3,4};\n"
         <<"Line (14) = {4,1};\n"
         <<"Line Loop (21) = {11, 12, 13, 14};\n"
         <<"Plane Surface (30) = {21};\n"
         <<"Physical Line (\"boundaries\") = {11,12,13,14};\n"
         <<"Physical Surface (\"Omega\") = {30};\n"
         ;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "benchmarkgrepl_geo";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

template<int Order>
void BenchmarkGreplNonlinearElliptic<Order>::initModel()
{
    M_Qj=2;
    M_Qr0=4;
    M_Qr1=1;
    M_Qa=1;
    M_Qf0=3;
    M_Qf1=1;

    M_use_newton = boption(_name="crb.use-newton");

    std::string mshfile_name = option("mshfile").as<std::string>();

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        double hsize=doption(_name="hsize");
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc = createGeo( hsize ) );
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
    auto Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
    }

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    auto mu_min = this->Dmu->element();
    mu_min <<  0.01, 0.01;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << 10, 10;
    this->Dmu->setMax( mu_max );

    M_mu = this->Dmu->element();

    auto Pset = this->Dmu->sampling();
    //specify how many elements we take in each direction
    std::vector<int> N(2);
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

    //auto eim_g = eim( _model=boost::enable_shared_from_this< BenchmarkGreplNonlinearElliptic<Order> >::shared_from_this(),
    auto eim_u = eim( _model=boost::dynamic_pointer_cast< BenchmarkGreplNonlinearElliptic<Order> >( this->shared_from_this() ),
                      _element=*pT,
                      _space=Xh,
                      _parameter=M_mu,
                      _expr= idv(*pT),
                      _sampling=Pset,
                      _name="eim_u" );


    auto eim_g = eim( _model=boost::dynamic_pointer_cast< BenchmarkGreplNonlinearElliptic<Order> >( this->shared_from_this() ),
                      _element=*pT,
                      _space=Xh,
                      _parameter=M_mu,
                      _expr=exp( cst_ref(M_mu(1)) * idv(*pT) ),
                      _sampling=Pset,
                      _name="eim_g" );

    M_funs.push_back( eim_g );
    M_funs.push_back( eim_u );
    int M = eim_g->mMax();
    int Mu = eim_u->mMax();

    if( Environment::worldComm().isMasterRank() )
        std::cout<<"M : "<<M<<" and Mu : "<<Mu<<std::endl;
    //resize data structure

    if( M_use_newton )
    {
        this->M_betaJqm.resize( M_Qj );
        this->M_betaJqm[0].resize( 1 );
        this->M_betaJqm[1].resize( M );
        this->M_betaRqm.resize( 2 );
        this->M_betaRqm[0].resize( M_Qr0 );
        this->M_betaRqm[0][0].resize( Mu );
        this->M_betaRqm[0][1].resize( M );
        this->M_betaRqm[0][2].resize( 1 );
        this->M_betaRqm[0][3].resize( 1 );
        this->M_betaRqm[1].resize( M_Qr1 );
        this->M_betaRqm[1][0].resize( 1 );

        this->M_Jqm.resize( M_Qj );
        this->M_Jqm[0].resize( 1 );
        this->M_Jqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        this->M_Jqm[1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        }
        this->M_Rqm.resize( 2 );
        this->M_Rqm[0].resize( M_Qr0 );
        this->M_Rqm[1].resize( M_Qr1 );
        this->M_Rqm[0][0].resize( Mu );
        for(int m=0; m<Mu; m++)
        {
            this->M_Rqm[0][0][m] = backend()->newVector( this->Xh );
        }
        this->M_Rqm[0][1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        }
        this->M_Rqm[0][2].resize(1);
        this->M_Rqm[0][2][0] = backend()->newVector( this->Xh );
        this->M_Rqm[0][3].resize(1);
        this->M_Rqm[0][3][0] = backend()->newVector( this->Xh );
        this->M_Rqm[1][0].resize( 1 );
        this->M_Rqm[1][0][0]= backend()->newVector( this->Xh );
    }
    else
    {
        this->M_betaAqm.resize( M_Qa );
        this->M_betaAqm[0].resize( 1 );
        this->M_betaFqm.resize(2);
        this->M_betaFqm[0].resize( M_Qf0 );
        this->M_betaFqm[0][0].resize( M );
        this->M_betaFqm[0][1].resize( 1 );
        this->M_betaFqm[0][2].resize( 1 );
        this->M_betaFqm[1].resize( M_Qf1 );
        this->M_betaFqm[1][0].resize( 1 );

        this->M_Aqm.resize( M_Qa );
        this->M_Aqm[0].resize( 1 );
        this->M_Aqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Fqm.resize(2);
        this->M_Fqm[0].resize(M_Qf0);
        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
        }
        this->M_Fqm[0][1].resize(1);
        this->M_Fqm[0][1][0]=backend()->newVector( this->Xh );
        this->M_Fqm[0][2].resize(1);
        this->M_Fqm[0][2][0]=backend()->newVector( this->Xh );
        this->M_Fqm[1].resize(1);
        this->M_Fqm[1][0].resize(1);
        this->M_Fqm[1][0][0]=backend()->newVector( this->Xh );
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();

    assemble();

} // BenchmarkGreplNonlinearElliptic::init

template <int Order>
void
BenchmarkGreplNonlinearElliptic<Order>::assembleJacobianWithAffineDecomposition( std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto u = Xh->element(); //trial
    auto eim_g = M_funs[0];
    int M = eim_g->mMax();
    double gamma = doption(_name="gamma");

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] )
        = integrate( _range= elements( mesh ),
                   _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] )
        += integrate( _range = markedfaces( mesh, "boundaries"),
                      _expr = gamma*idt(u)*id(v)/hFace()
                      - (gradt(u)*vf::N())*id(v) +
                      - (grad(v)*vf::N())*idt(u) );

    Jqm[0][0]->close();
    for(int m=0; m<M; m++)
    {
        form2( _test=Xh, _trial=Xh, _matrix=Jqm[1][m] ) =
            integrate( _range = elements(mesh),
                       _expr = idv(eim_g->q(m))*idt(u)*id(v) );
        Jqm[1][m]->close();
    }

}

//BenchmarkGreplNonlinearElliptic<Order>::assembleResidualWithAffineDecomposition(element_type const& X,
//                                                                                std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
template <int Order>
void
BenchmarkGreplNonlinearElliptic<Order>::assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto eim_g = M_funs[0];
    auto eim_u = M_funs[1];
    int M = eim_g->mMax();
    int Mu = eim_u->mMax();
    double gamma = doption(_name="gamma");


    for(int m=0; m<Mu; m++)
    {
        auto q = eim_u->q(m);
        form1( _test=Xh, _vector=Rqm[0][0][m] )
            = integrate( _range= elements( mesh ),
                         _expr = gradv(q)*trans(grad(v)) );
        form1( _test=Xh, _vector=Rqm[0][0][m] )
            += integrate( _range = markedfaces( mesh, "boundaries"),
                          _expr = gamma*idv(q)*id(v)/hFace()
                          - (gradv(q)*vf::N())*id(v) +
                          - (grad(v)*vf::N())*idv(q) );

        Rqm[0][0][m]->close();
    }
    for(int m=0; m<M; m++)
    {
        form1( _test=Xh, _vector=Rqm[0][1][m] ) =
            integrate( _range= elements( mesh ),
                       _expr=( idv(eim_g->q(m))*id(v) ) );
        Rqm[0][1][m]->close();
    }

    form1( _test=Xh, _vector=Rqm[0][2][0] ) =
        integrate( _range= elements( mesh ),
                   _expr= id(v) );

    form1( _test=Xh, _vector=Rqm[0][3][0] ) =
        integrate( _range= elements( mesh ),
                   _expr=sin(2*M_PI*Px())*cos(2*M_PI*Py()) * id(v) );

    form1( _test=Xh, _vector=Rqm[1][0][0] ) =
        integrate( _range= elements( mesh ),
                   _expr=id(v) );

    Rqm[0][2][0]->close();
    Rqm[0][3][0]->close();
    Rqm[1][0][0]->close();

}

template <int Order>
void
BenchmarkGreplNonlinearElliptic<Order>::updateJacobianMonolithic( vector_ptrtype const& X,
                                                                  sparse_matrix_ptrtype & J,
                                                                  parameter_type const& mu )
{

    //Here we don't use EIM approximation
    //and so no affine decomposition
    auto Xh = this->Xh;
    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    double gamma = doption(_name="gamma");
    auto g = exp( mu(1)*idv(u) );
    auto proj_g = vf::project(_space=Xh, _expr=g);

    form2( _test=Xh, _trial=Xh, _matrix=J )
        = integrate( _range= elements( mesh ),
                   _expr = gradt(u)*trans(grad(v)) ) ;
    form2( _test=Xh, _trial=Xh, _matrix=J )
        += integrate( _range = markedfaces( mesh, "boundaries"),
                     _expr = gamma*idt(u)*id(v)/hFace()
                     - (gradt(u)*vf::N())*id(v) +
                     - (grad(v)*vf::N())*idt(u) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range = elements(mesh),
                   _expr = mu(0) * idv(proj_g)*idt(u)*id(v) );

    J->close();

}

template <int Order>
void
BenchmarkGreplNonlinearElliptic<Order>::updateResidualMonolithic(vector_ptrtype const& X,
                                                                 vector_ptrtype & R,
                                                                 parameter_type const& mu )
{
    auto Xh = this->Xh;
    double gamma = doption(_name="gamma");

    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    auto g = exp( mu(1)*idv(u) );
    auto proj_g = vf::project(_space=Xh, _expr=g);

    form1( _test=Xh, _vector=R ) =
        integrate( _range= elements( mesh ),
                   _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=R )
        += integrate( _range = markedfaces( mesh, "boundaries"),
                      _expr = gamma*idv(u)*id(v)/hFace()
                      - (gradv(u)*vf::N())*id(v) +
                      - (grad(v)*vf::N())*idv(u) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   _expr= mu(0)/mu(1)*( idv(proj_g) * id(v) ) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   _expr= -mu(0)/mu(1)*( id(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   _expr=-100*sin(2*M_PI*Px())*cos(2*M_PI*Py()) * id(v) );

    R->close();

}



template<int Order>
void BenchmarkGreplNonlinearElliptic<Order>::assemble()
{
    auto Xh = this->Xh;
    auto u = Xh->element();

    if( M_use_newton )
    {
        assembleJacobianWithAffineDecomposition( this->M_Jqm );
        assembleResidualWithAffineDecomposition( this->M_Rqm );
    }
    else
    {
        auto v = Xh->element();
        auto eim_g = M_funs[0];
        int M = eim_g->mMax();
        double gamma = doption(_name="gamma");

        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) =
            integrate( _range= elements( mesh ),
                       _expr = gradt(u)*trans(grad(v)) );
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] )
            += integrate( _range = markedfaces( mesh, "boundaries"),
                          _expr = gamma*idt(u)*id(v)/hFace()
                          - (gradt(u)*vf::N())*id(v) +
                          - (grad(v)*vf::N())*idt(u) );

        for(int m=0; m<M; m++)
        {
            form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) =
                integrate( _range= elements( mesh ),
                           _expr=( idv(eim_g->q(m))*id(v) ) );
            this->M_Fqm[0][0][m]->close();
        }

        form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
            integrate( _range= elements( mesh ),
                       _expr= id(v) );
        this->M_Fqm[0][1][0]->close();

        form1( _test=Xh, _vector=this->M_Fqm[0][2][0] ) =
            integrate( _range= elements( mesh ),
                       _expr=sin(2*M_PI*Px())*cos(2*M_PI*Py()) * id(v) );
        this->M_Fqm[0][2][0]->close();

        form1( _test=Xh, _vector=this->M_Fqm[1][0][0] ) =
            integrate( _range= elements( mesh ),
                       _expr=id(v) );
        this->M_Fqm[1][0][0]->close();

    }//no newton

}

template<int Order>
typename BenchmarkGreplNonlinearElliptic<Order>::monolithic_type
BenchmarkGreplNonlinearElliptic<Order>::computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution)
{
    auto Xh = this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    //auto solution=this->solve(mu);

    auto exprg = exp( mu(1)*idv(solution) );
    auto g = vf::project( _space=Xh, _expr=exprg );

    M_monoA = backend()->newMatrix( Xh, Xh );
    M_monoF.resize(2);
    M_monoF[0] = backend()->newVector( this->functionSpace() );
    M_monoF[1] = backend()->newVector( Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) =
        integrate( _range= elements( mesh ),
                   _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA )
        += integrate( _range = markedfaces( mesh, "boundaries"),
                      _expr = gamma*idt(u)*id(v)/hFace()
                      - (gradt(u)*vf::N())*id(v) +
                      - (grad(v)*vf::N())*idt(u) );

    form1( _test=Xh, _vector=M_monoF[0] ) =
        integrate( _range= elements( mesh ),
                   _expr=-mu(0)/mu(1)*( idv(g)*id(v) - id(v) ) );
    form1( _test=Xh, _vector=M_monoF[0] )
        += integrate( _range= elements( mesh ),
                     _expr=100*sin(2*M_PI*Px())*cos(2*M_PI*Py()) * id(v) );

    form1( _test=Xh, _vector=M_monoF[1] ) =
        integrate( _range= elements( mesh ),
                   _expr=id(v) );

    M_monoF[0]->close();
    M_monoF[1]->close();
    M_monoA->close();

    return boost::make_tuple( M_monoA, M_monoF );

}

template<int Order>
typename BenchmarkGreplNonlinearElliptic<Order>::element_type
BenchmarkGreplNonlinearElliptic<Order>::solve( parameter_type const& mu )
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


template<int Order>
typename BenchmarkGreplNonlinearElliptic<Order>::vector_sparse_matrix
BenchmarkGreplNonlinearElliptic<Order>::computeLinearDecompositionA()
{
    auto Xh=this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    this->M_linearAqm.resize(1);
    this->M_linearAqm[0].resize(1);
    this->M_linearAqm[0][0] = backend()->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] ) =
        integrate( _range=elements( mesh ),
                   _expr=gradt( u )*trans( grad( v ) ) );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] )
        = integrate( _range=markedfaces( mesh, "boundaries" ),
                      _expr=gamma*idt( u )*id( v )/h()
                      -gradt( u )*vf::N()*id( v )
                      -grad( v )*vf::N()*idt( u ) );
    this->M_linearAqm[0][0]->close();
    return this->M_linearAqm;

}

template<int Order>
typename BenchmarkGreplNonlinearElliptic<Order>::affine_decomposition_type
BenchmarkGreplNonlinearElliptic<Order>::computeAffineDecomposition()
{
    if( M_use_newton )
        return boost::make_tuple( this->M_Jqm, this->M_Rqm );
    else
        return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order>
std::vector< std::vector< typename BenchmarkGreplNonlinearElliptic<Order>::element_ptrtype > >
BenchmarkGreplNonlinearElliptic<Order>::computeInitialGuessAffineDecomposition( )
{
    return M_InitialGuess;
}


template<int Order>
double BenchmarkGreplNonlinearElliptic<Order>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve , bool export_outputs )
{

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;
    if ( output_index==0 )
    {
        s = integrate( _range= elements( mesh ),
                       _expr=-100*sin(2*M_PI*Px())*cos(2*M_PI*Py()) * idv( solution ) ).evaluate()(0,0);
    }
    if( output_index==1 )
    {
        s = integrate( elements( mesh ), idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}

}

#endif /* __BenchmarkGreplNonlinearElliptic_H */
