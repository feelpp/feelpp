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


namespace Feel
{

po::options_description
makeBenchmarkGreplNonLinearParabolicOptions()
{
    po::options_description bgoptions( "BenchmarkGreplNonLinearParabolic options" );
    bgoptions.add_options()
        ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "hsize")
        ( "trainset-eim-size", Feel::po::value<int>()->default_value( 15 ), "EIM trainset is built using a equidistributed grid 15*15 by default")
        ( "gamma", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
        ;
    return bgoptions.add( bdf_options( "heat" ) );
}
AboutData
makeBenchmarkGreplNonLinearParabolicAbout( std::string const& str = "benchmarkGrepl" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Benchmark Grepl",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Vincent Huber", "developer", "vincent.huber@cemosis.fr", "" );
    return about;
}


class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 3;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int Order>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<3,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef bases<Lagrange<Order, Scalar> > basis_type_eimg;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef FunctionSpace<mesh_type, basis_type_eimg, value_type> space_type_eimg;

    typedef typename space_type::element_type element_type;
    
    static const bool is_time_dependent = true;
    static const bool is_linear = false;

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

/**
 * \class BenchmarkGreplNonLinearParabolic
 * \brief brief description
 *
 * This is from the paper
 * Certified Reduced Basis Methods for Pärametrized 
 * Parabolic Partial Differential Equations with 
 * Non-Affine Source Terms
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @author Vincent Huber
 * @see
 */
template<int Order>
class BenchmarkGreplNonLinearParabolic :
    public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<Order> ,NonLinear&TimeDependent, EimDefinition<ParameterDefinition, FunctionSpaceDefinition<Order> > >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition<Order>, NonLinear&TimeDependent, EimDefinition<ParameterDefinition,FunctionSpaceDefinition<Order> > > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    typedef double value_type;

    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename FunctionSpaceDefinition<Order>::space_type space_type;
    typedef typename boost::shared_ptr<space_type> space_ptrtype;
    typedef typename FunctionSpaceDefinition<Order>::space_type_eimg space_type_eimg;
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
  
    typedef BenchmarkGreplNonLinearParabolic<Order> self_type;
    
    typedef typename super_type::bdf_ptrtype bdf_ptrtype;

    //! initialization of the model
    void initModel();
    //@}

    std::string modelName()
    {
        std::ostringstream ostr;
        ostr << "BenchMarkGreplNonLinearParabolic" <<  Order;
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
        auto eim_g = M_funs[0];
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
    computeBetaQm( parameter_type const& mu )
    {
        auto eim_g = M_funs[0];
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
        auto eim_g = M_funs[0];
        vectorN_type beta_g = eim_g->beta( mu , T );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    beta_type
    computePicardBetaQm( parameter_type const& mu )
    {
        auto eim_g = M_funs[0];
        vectorN_type beta_g = eim_g->beta( mu );

        std::vector<vectorN_type*> betas;
        betas.push_back(&beta_g);

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    void fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu)
    {
        auto eim_g = M_funs[0];
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
        if( !M_use_newton || boption(_name="ser.error-estimation") )
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
    
    bdf_ptrtype bdfModel(){ return M_bdf; }

private:

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_ptrtype pT;
    parameter_type M_mu;

    funs_type M_funs;

   sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    bool M_use_newton;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;
    bdf_ptrtype M_bdf;

};

template<int Order>
void BenchmarkGreplNonLinearParabolic<Order>::initModel()
{

    M_use_newton = boption(_name="crb.use-newton");


    /*
     * First we create the mesh or load it if already exist
     */
    mesh = loadMesh(_mesh = new mesh_type );


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
    
    M_bdf = bdf(_space=Xh, _name="heat", _prefix="heat");

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );
    static auto M_U = *pT;

    auto mu_min = this->Dmu->element();
    mu_min <<  1.08, 1.08, 1.08, 0.5;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << 1.32, 1.32, 1.32, 2;
    this->Dmu->setMax( mu_max );

    M_mu = this->Dmu->element();

    auto Pset = this->Dmu->sampling();
    //specify how many elements we take in each direction
    std::vector<int> N(this->Dmu->dimension());
    int Ne = ioption(_name="trainset-eim-size");
    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();

    std::ifstream file ( supersamplingname );

    //40 elements in each direction
    N[0]=Ne; 
    N[1]=Ne;
    N[2]=Ne;
    N[3]=1;

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

    /*
     * The EIM expansion does not depends on the FE solution !
     * We have to instanciate the eim_no_solve model
     */
    LOG(INFO) << "Starting EIM\n";
    auto eim_g = eim( _model=eim_no_solve(super_type::shared_from_this()),
                      _element=*pT,
                      _space=Xh_eimg,
                      _parameter=M_mu,
                      _expr=
                       exp( -(Px()*Px())/(cst_ref(M_mu(0))*cst_ref(M_mu(0))))
                      *exp( -(Py()*Py())/(cst_ref(M_mu(1))*cst_ref(M_mu(1))))
                      *exp( -(Pz()*Pz())/(cst_ref(M_mu(2))*cst_ref(M_mu(2)))),
                      _sampling=Pset,
                      _name="eim_g" );
    LOG(INFO) << "EIM Finished\n";

    int cpt = 0;
    auto exporter_eim = exporter(_mesh=mesh,_name="eim");
    BOOST_FOREACH( auto mu, *Pset )
    {
      Feel::cout<<"check gM for mu = ["<< mu(0)<<" , "<<mu(1)<<" , "<<mu(2)<<"]"<<std::endl;
      auto g = vf::project(_space=Xh,
                           _range=elements(mesh),
                           _expr= exp( -(Px()*Px())/(mu(0)*mu(0)))
                                 *exp( -(Py()*Py())/(mu(1)*mu(1)))
                                 *exp( -(Pz()*Pz())/(mu(2)*mu(2))));
      auto gm = vf::project(_space=Xh_eimg, _range=elements(mesh), _expr = idv(eim_g->operator()(mu)));
      exporter_eim->step(cpt)->add("g", g);
      exporter_eim->step(cpt)->add("gm", gm);
      exporter_eim->save();
      cpt++;
    }

    M_funs.push_back( eim_g );

    assemble();

} // BenchmarkGreplNonLinearParabolic::init



template<int Order>
void BenchmarkGreplNonLinearParabolic<Order>::assemble()
{
    auto Xh = this->Xh;
    auto u = Xh->element();

    // Mass Matrix
    auto mass = form2(_test=Xh, _trial=Xh);
    mass = integrate(_range=elements(mesh), _expr=inner(idv(u),id(u)));
    this->addMass( {mass, "1"} );

    // Evolution
    ///TODO Mettre le nombre de Peclet
    auto convection = form2(_test=Xh, _trial=Xh);
    convection = integrate(_range=elements(mesh), 
        _expr=inner(vec(-0.1,0,0)*gradt(u),id(u)));
    this->addLhs( { convection, "1"} );

    // Diffusion
    auto diffusion = form2(_test=Xh, _trial=Xh);
    diffusion = integrate(_range=elements(mesh), 
        _expr=inner(gradt(u), grad(u)));
    this->addLhs( { convection, "mu3"} );

    // Rhs (eim)
    // g = sum_{i=1}^{mMax} w_m(mu) g(x)
    for(int i=0; i < M_funs[0]->mMax(); i++)
    {
      auto lhs = form1(_test=Xh);
      lhs = integrate(_range=elements(mesh),
          _expr=inner(idv(M_funs[0]->q(i)),id(u)));
      // je veux que M_funs[0]->beta(mu) * lhs
      this->addLhs( { lhs, "mu0:mu1:mu2"} );
    }

}

template<int Order>
double BenchmarkGreplNonLinearParabolic<Order>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve )
{
    //CHECK( need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;
    if ( output_index==0 )
    {
        s = integrate( _range= markedelements( mesh , "v1"), _expr=idv( solution ) ).evaluate()(0,0);
        s /= (0.5*0.5*0.1); // Volume of v1
    }
    if( output_index==1 )
    {
        s = integrate( elements( mesh ), idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}

}

#endif /* __BenchmarkGreplNonLinearParabolic_H */
