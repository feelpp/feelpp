#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

po::options_description
makeOptions()
{
    po::options_description options( "Thermoelectric" );
    options.add_options()
        ( "thermoelectric.filename", Feel::po::value<std::string>()->default_value("thermoelectric.json"),
          "json file containing application parameters")
        ( "thermoelectric.gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( "thermoelectric.sigma", po::value<double>()->default_value(1), "electric conductivity" )
        ( "thermoelectric.current", po::value<double>()->default_value(1), "current intensity" )
        ( "thermoelectric.k", po::value<double>()->default_value(1), "thermal conductivity" )
        ( "thermoelectric.h", po::value<double>()->default_value(1), "transfer coefficient" )
        ( "thermoelectric.Tw", po::value<double>()->default_value(1), "water's temperatur" )
        ( "thermoelectric.N", po::value<int>()->default_value(100), "number of basis function to use" )
        ( "thermoelectric.trainset-eim-size", po::value<int>()->default_value(40), "size of the eim trainset" )
        ;
    options.add(backend_options("mono") );
    return options;
}

AboutData
makeAbout( std::string const& str = "thermoelectriccrbmodel" )
{
    AboutData about( str.c_str() );
    return about;
}

class ParameterDefinition
{
public:
    using parameterspace_type = ParameterSpace<5>;
};

class FunctionSpaceDefinition
{
public:
    using value_type = double;
    static const uint16_type Order = 1;
    using convex_type =  Simplex<3>;
    using mesh_type = Mesh<convex_type>;
    using fct_base_type = Lagrange<Order,Scalar>;
    using basis_type = bases<fct_base_type, fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    using V_basis_type = bases<fct_base_type>;
    using V_space_type = FunctionSpace<mesh_type, V_basis_type, value_type>;

    using J_basis_type = bases<Lagrange<Order+2, Scalar,Discontinuous> >;
    using J_space_type = FunctionSpace<mesh_type, J_basis_type, Discontinuous, value_type>;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef typename FunctionSpaceDefinition::J_space_type J_space_type;
    using V_space_type = typename FunctionSpaceDefinition::V_space_type;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<V_space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<J_space_type, space_type , parameterspace_type> fund_type;

};

class Thermoelectric : public ModelCrbBase<ParameterDefinition,
                                           FunctionSpaceDefinition,
                                           NonLinear,
                                           EimDefinition<ParameterDefinition,
                                                         FunctionSpaceDefinition> >
{
public:
    using super_type = ModelCrbBase<ParameterDefinition,
                                    FunctionSpaceDefinition,
                                    NonLinear,
                                    EimDefinition<ParameterDefinition,
                                                  FunctionSpaceDefinition> >;
    using value_type = double;
    using element_type = super_type::element_type;
    using element_ptrtype = super_type::element_ptrtype;
    using J_space_type = FunctionSpaceDefinition::J_space_type;
    using V_space_type = FunctionSpaceDefinition::V_space_type;
    using V_element_type = typename V_space_type::element_type;
    using V_view_ptrtype = typename element_type::template sub_element_ptrtype<0>;
    using T_view_ptrtype = typename element_type::template sub_element_ptrtype<1>;
    using parameter_type = super_type::parameter_type;
    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;
    using prop_type = ModelProperties;
    using prop_ptrtype = boost::shared_ptr<prop_type>;
    using vectorN_type = super_type::vectorN_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using beta_type = boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >;
    using affine_decomposition_type = typename super_type::affine_decomposition_type;
    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;

private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;
    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    element_ptrtype pT;
    V_view_ptrtype M_V;
    T_view_ptrtype M_T;
    parameter_type M_mu;

public:
    // Constructors
    Thermoelectric();
    Thermoelectric( mesh_ptrtype mesh );

    // Helpers
    int Qa();
    int Nl();
    int Ql( int l );
    int mQA( int q );
    int mLQF( int l, int q );
    int mCompliantQ( int q );
    void resizeQm();

    // Pure virtual
    void initModel();

    void decomposition();

    beta_vector_type computeBetaInitialGuess( parameter_type const& mu );
    beta_type computeBetaQm( element_type const& T, parameter_type const& mu );
    beta_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu );
    beta_type computeBetaQm( parameter_type const& mu );
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, double time=1e30 );
    void fillBetaQm( parameter_type const& mu, vectorN_type betaEimGradGrad );

    affine_decomposition_type computeAffineDecomposition();
    std::vector<std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA();
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition();

    element_type solve( parameter_type const& mu );
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    int mMaxSigma();
    auto eimSigmaQ(int m);
    vectorN_type eimSigmaBeta( parameter_type const& mu );
    template<typename vec_space_type>
    typename vec_space_type::element_type
    computeTruthCurrentDensity( parameter_type const& mu );
}; // Thermoelectric class

Thermoelectric::Thermoelectric()
    : super_type( "thermoelectric" )
{}

Thermoelectric::Thermoelectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric" )
{
    this->M_mesh = mesh;
}

int Thermoelectric::Qa()
{
    return 4;
}

int Thermoelectric::mQA( int q )
{
    return 1;
}

int Thermoelectric::Nl()
{
    return 1;
}

int Thermoelectric::Ql( int l)
{
    switch(l) {
    case 0:
        return 3;
    default:
        return 0;
    }
}

int Thermoelectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    default:
        return 0;
    }
}

int Thermoelectric::mCompliantQ(int q )
{
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    switch( q )
    {
    case 0:
        return 1;
    case 1:
        return eimGradGrad->mMax();
    case 2:
        return 1;
    default:
        return 0;
    }
}

void Thermoelectric::resizeQm()
{
    M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        M_Aqm[q].resize(mQA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mQA(q));
    }

    M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();
}

void Thermoelectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    M_modelProps = boost::make_shared<prop_type>(Environment::expand( soption("thermoelectric.filename")));
    auto parameters = M_modelProps->parameters();
    if ( parameters.size() != parameterSpace()->dimension() )
        LOG(FATAL) << "number of parameters in thermoelectric.parameters (" << parameters.size()
                   << ") different from dimension of the parameter space ("
                   << parameterSpace()->dimension() << ")";
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        mu_min(i) = parameterPair.second.min();
        mu_max(i) = parameterPair.second.max();
        Dmu->setParameterName(i++, parameterPair.first );
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( M_mesh ) );

    if( !pT )
        pT = element_ptrtype( new element_type( Xh ) );
    M_V = pT->template elementPtr<0>();
    M_T = pT->template elementPtr<1>();

    auto JspaceEim = J_space_type::New( M_mesh );

    auto Pset = this->Dmu->sampling();
    int Ne = ioption(_name="thermoelectric.trainset-eim-size");
    std::vector<size_type> N(parameterSpace()->dimension(), Ne);

    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();
    std::ifstream file ( supersamplingname );
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

    auto eim_gradgrad = eim( _model=boost::dynamic_pointer_cast<Thermoelectric>(this->shared_from_this() ),
                             _element=*M_V,
                             _parameter=M_mu,
                             _expr=M_mu.parameterNamed("sigma")*inner(gradv(*M_V)),
                             _space=JspaceEim,
                             _name="eim_gradgrad",
                             _sampling=Pset );
    this->addEimDiscontinuous( eim_gradgrad );

    this->resizeQm();
    this->decomposition();
}

void Thermoelectric::decomposition()
{
    auto UT = Xh->element();
    auto VP = Xh->element();
    auto u = UT.template element<0>();
    auto v = VP.template element<0>();
    auto t = UT.template element<1>();
    auto p = VP.template element<1>();

    auto gamma = doption("thermoelectric.gamma");

    auto eimGradGrad = this->scalarDiscontinuousEim()[0];

    /***************************** Electro *****************************/
    auto a1 = form2(_test=Xh, _trial=Xh);
    a1 = integrate( elements(M_mesh),
                    inner(gradt(u),grad(v)) );
    M_Aqm[0][0] = a1.matrixPtr();

    // Dirichlet condition
    auto a2 = form2(_test=Xh, _trial=Xh);
    a2 += integrate( markedfaces(M_mesh, {"in", "out"}),
                     gamma/hFace()*inner(idt(u),id(v))
                     -inner(gradt(u)*N(),id(v))
                     -inner(grad(v)*N(),idt(u)) );
    M_Aqm[1][0] = a2.matrixPtr();

    // Dirichlet condition
    auto f1 = form1(_test=Xh);
    f1 = integrate( markedfaces(M_mesh, "in"),
                    gamma/hFace()*id(v) -  grad(v)*N() );
    M_Fqm[0][0][0] = f1.vectorPtr();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(u),grad(v)) );
    M_energy_matrix = m.matrixPtr();

    /***************************** Thermo *****************************/
    auto a3 = form2(_test=Xh, _trial=Xh);
    a3 = integrate( elements(M_mesh),
                    inner(gradt(t), grad(p)) );
    M_Aqm[2][0] = a3.matrixPtr();
    // Robin condition
    auto a4 = form2(_test=Xh, _trial=Xh);
    a4 = integrate( markedfaces(M_mesh, {"Rint","Rext"}),
                    inner(id(t), idt(p)) );
    M_Aqm[3][0] = a4.matrixPtr();

    // right hand side
    for( int m = 0; m < eimGradGrad->mMax(); ++m )
    {
        auto f2 = form1(_test=Xh);
        f2 = integrate(elements(M_mesh),
                       inner(id(p), idv(eimGradGrad->q(m))) );
        M_Fqm[0][1][m] = f2.vectorPtr();
    }

    // Robin condition
    auto f3 = form1(_test=Xh);
    f3 = integrate( markedfaces(M_mesh, {"Rint","Rext"}),
                    id(p) );
    M_Fqm[0][2][0] = f3.vectorPtr();
}

Thermoelectric::beta_vector_type
Thermoelectric::computeBetaInitialGuess( parameter_type const& mu )
{
    Feel::cout << "computeBetaInitalGuess" << std::endl;
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

Thermoelectric::beta_type
Thermoelectric::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    Feel::cout << "computeBetaQm( element_type const& T, parameter_type const& mu )" << std::endl;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu, T );
    this->fillBetaQm(mu, betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Thermoelectric::beta_type
Thermoelectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    Feel::cout << "computeBetaQm( vectorN_type const& urb, parameter_type const& mu )" << std::endl;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu );
    this->fillBetaQm(mu, betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Thermoelectric::beta_type
Thermoelectric::computeBetaQm( parameter_type const& mu )
{
    Feel::cout << "computeBetaQm( parameter_type const& mu )" << std::endl;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu );
    this->fillBetaQm(mu, betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void Thermoelectric::fillBetaQm( parameter_type const& mu, vectorN_type betaEimGradGrad )
{
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    // auto betaEimGradGrad = eimGradGrad->beta( mu );
    Feel::cout << "eimGradGrad M = " << eimGradGrad->mMax() << " and beta =" << std::endl
               << betaEimGradGrad << std::endl;

    M_betaAqm[0][0] = mu.parameterNamed("sigma");
    M_betaAqm[1][0] = mu.parameterNamed("sigma");
    M_betaFqm[0][0][0] = mu.parameterNamed("sigma")*mu.parameterNamed("current");
    M_betaAqm[2][0] = mu.parameterNamed("k");
    M_betaAqm[3][0] = mu.parameterNamed("h");
    for( int m = 0; m < eimGradGrad->mMax(); ++m )
        M_betaFqm[0][1][m] = betaEimGradGrad(m);
    M_betaFqm[0][2][0] = mu.parameterNamed("h")*mu.parameterNamed("Tw");
}

Thermoelectric::beta_vector_type
Thermoelectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    Feel::cout << "computeBetaLinearDecompositionA" << std::endl;
    beta_vector_type beta;
    return beta;
}

Thermoelectric::affine_decomposition_type
Thermoelectric::computeAffineDecomposition()
{
    Feel::cout << "computeAffineDecomposition" << std::endl;
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<Thermoelectric::sparse_matrix_ptrtype> >
Thermoelectric::computeLinearDecompositionA()
{
    Feel::cout << "computeLinearDecompositionA" << std::endl;
    return this->M_linearAqm;
}

std::vector<std::vector<Thermoelectric::element_ptrtype> >
Thermoelectric::computeInitialGuessAffineDecomposition()
{
    Feel::cout << "computeInitialGuessAffineDecomposition" << std::endl;
    return M_InitialGuess;
}

Thermoelectric::element_type
Thermoelectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();
    // auto VT = Xh->element();
    // auto V = VT.template element<0>();
    // auto T = VT.template element<1>();
    // auto phi = Xh->element();
    // auto phiV = phi.template element<0>();
    // auto phiT = phi.template element<1>();

    auto gamma = doption("thermoelectric.gamma");
    auto sigma = mu.parameterNamed("sigma");
    auto current = mu.parameterNamed("current");
    auto k = mu.parameterNamed("k");
    auto h = mu.parameterNamed("h");
    auto Tw = mu.parameterNamed("Tw");

    /***************************** Electro *****************************/
    tic();
    auto a = form2(_test=Vh, _trial=Vh);
    // V
    a = integrate( elements(M_mesh),
                   sigma*inner(gradt(V),grad(phiV)) );

    // V Dirichlet condition
    a += integrate( markedfaces(M_mesh, {"in", "out"}),
                    sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                           -inner(gradt(V)*N(),id(phiV))
                           -inner(grad(phiV)*N(),idt(V)) ) );

    auto f = form1(_test=Vh);
    // V Dirichlet condition
    f = integrate( markedfaces(M_mesh, "in"),
                   sigma*current*(gamma/hFace()*id(phiV) -  grad(phiV)*N() ) );

    a.solve( _solution=V, _rhs=f, _name="mono" );

    auto aT = form2(_test=Th, _trial=Th);
    // T
    aT += integrate( elements(M_mesh),
                     k*inner(gradt(T), grad(phiT)) );

    // T Robin condition
    aT += integrate( markedfaces(M_mesh, {"Rint","Rext"}),
                     h*inner(idt(T), id(phiT)) );

    auto fT = form1(_test=Th);
    // T right hand side
    fT = integrate(elements(M_mesh),
                   id(phiT)*sigma*gradv(V)*trans(gradv(V)) );

    // T Robin condition
    fT += integrate( markedfaces(M_mesh, {"Rint","Rext"}),
                     h*Tw*id(phiT) );

    aT.solve( _solution=T, _rhs=fT, _name="mono" );

    auto solution = Xh->element();
    solution.template element<0>() = V;
    solution.template element<1>() = T;

    auto e = exporter(M_mesh);
    e->add("sol", solution);
    e->save();
    toc("mono");

    return solution;
}

double Thermoelectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == 0 )
    {
        for ( int q = 0; q < Ql(0); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *M_Fqm[output_index][q][0];
            output += M_betaFqm[output_index][q][0]*dot( *eltF, u );
            //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

int Thermoelectric::mMaxSigma()
{
    return 1;
}

auto Thermoelectric::eimSigmaQ(int m)
{
    auto Vh = Xh->template functionSpace<0>();
    auto q = vf::project( _space=Vh, _range=elements(M_mesh), _expr=cst(1.) );
    return q;
}

Thermoelectric::vectorN_type Thermoelectric::eimSigmaBeta( parameter_type const& mu )
{
    vectorN_type beta(1);
    beta(0) = mu.parameterNamed( "sigma" );
    return beta;
}

template<typename vec_space_type>
typename vec_space_type::element_type
Thermoelectric::computeTruthCurrentDensity( parameter_type const& mu )
{
    auto VT = this->solve(mu);
    auto V = VT.template element<0>();
    auto sigma = mu.parameterNamed("sigma");
    auto Vh = Xh->template functionSpace<0>();
    auto j = Vh->element();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
    return j;
}

} // namespace Feel
