#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmesh/complement.hpp>

#include <boost/algorithm/string.hpp>

namespace Feel {

inline
po::options_description
makeMixedPoissonOptions( std::string prefix = "" )
{
    po::options_description mpOptions( "Mixed Poisson HDG options");
    mpOptions.add_options()
        ( prefixvm( prefix, "tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix, "tau_order").c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "picard.itol").c_str(), po::value<double>()->default_value( 1e-4 ), "tolerance" )
        ( prefixvm( prefix, "picard.itmax").c_str(), po::value<int>()->default_value( 10 ), "iterations max" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( prefixvm( prefix, "conductivity_json").c_str(), po::value<std::string>()->default_value( "cond" ), "key for conductivity in json" )
        ( prefixvm( prefix, "conductivityNL_json").c_str(), po::value<std::string>()->default_value( "cond" ), "key for non linear conductivity in json (depends on potential p)" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ;
    return mpOptions;
}

inline po::options_description
makeMixedPoissonLibOptions( std::string prefix = "" )
{
    po::options_description mpLibOptions( "Mixed Poisson HDG Lib options");
    if ( !prefix.empty() )
        mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order>
class MixedPoisson
{
public:

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,1,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;

    // Vh
    using Vh_t = Pdhv_type<mesh_type,Order>;
    using Vh_ptr_t = Pdhv_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
    // Wh
    using Wh_t = Pdh_type<mesh_type,Order>;
    using Wh_ptr_t = Pdh_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
    // Mh
    using Mh_t = Pdh_type<face_mesh_type,Order>;
    using Mh_ptr_t = Pdh_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
    // Ch
    using Ch_t = Pch_type<mesh_type,0>;
    using Ch_ptr_t = Pch_ptrtype<mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    // M0h
    using M0h_t = Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t = Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    using linearAssembly_function_type = boost::function<void ( sparse_matrix_ptrtype& A,vector_ptrtype& F )>;

private:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh; // potential trace
    Ch_ptr_t M_Ch; // Lagrange multiplier
    M0h_ptr_t M_M0h;

    backend_ptrtype M_backend;
    sparse_matrix_ptrtype M_A;
    sparse_matrix_ptrtype M_A_cst;
    vector_ptrtype M_F;
    BlocksBaseVector<double> M_hdg_sol;
    vector_ptrtype M_U;

    Vh_element_ptr_t M_up;
    Wh_element_ptr_t M_pp;

    int M_tau_order;

    bool M_integralCondition;
    bool M_isPicard;

    std::list<std::string> M_dirichletMarkersList;
    std::list<std::string> M_neumannMarkersList;
    std::list<std::string> M_robinMarkersList;
    std::list<std::string> M_integralMarkersList;

    void initGraphs();
    void initGraphsWithIntegralCond();
    void assembleACst();
    void assembleA();

public:
    linearAssembly_function_type M_updateAssembly;

    MixedPoisson( std::string prefix = "" );
    void init( mesh_ptrtype mesh = NULL);
    void solve();
    void assembleF();
    void solveNL();
    template<typename ExprT>
    void updateConductivityTerm( Expr<ExprT> expr, std::string marker = "");
    void updateConductivityTerm( bool isNL = false);
    template<typename ExprT>
    void updatePotentialRHS( Expr<ExprT> expr, std::string marker = "");
    template<typename ExprT>
    void updateFluxRHS( Expr<ExprT> expr, std::string marker = "");
    void exportResults();

    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    Vh_element_ptr_t fluxField() const { return M_up; }
    Wh_element_ptr_t potentialField() const { return M_pp; }
    model_prop_type modelProperties() const { return *M_modelProperties; }
    std::list<std::string> integralMarkersList() const { return M_integralMarkersList; }
};

template<int Dim, int Order>
MixedPoisson<Dim, Order>::MixedPoisson(std::string prefix )
{
    M_prefix = prefix;
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::init( mesh_ptrtype mesh)
{
    tic();

    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;

    // initialize marker lists for each boundary condition type
    M_dirichletMarkersList.clear();
    M_neumannMarkersList.clear();
    M_robinMarkersList.clear();
    M_integralMarkersList.clear();
    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                M_dirichletMarkersList.push_back(marker);
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                M_neumannMarkersList.push_back(marker);
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                M_robinMarkersList.push_back(marker);
            }
        }
    }
    itField = M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                M_integralMarkersList.push_back(marker);
            }
        }
    }
    cout << "Dirichlet : " << M_dirichletMarkersList << std::endl
         << "Neumann : " << M_neumannMarkersList << std::endl
         << "Robin : " << M_robinMarkersList << std::endl;

    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement(faces(M_mesh),[this]( auto const& e ) {
            for( auto marker : this->M_integralMarkersList)
            {
                if ( e.marker().value() == this->M_mesh->markerName( marker ) )
                    return true;
            }
            return false; });
    auto face_mesh = createSubmesh( M_mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );

    toc("mesh",true);

    if ( M_integralMarkersList.empty() )
        M_integralCondition = false;
    else
        M_integralCondition = true;
    if ( boost::icontains(M_modelProperties->model(),"picard") )
        M_isPicard = true;
    else
        M_isPicard = false;

    cout << "Model : " << M_modelProperties->model()
         << " using ";
    if ( M_integralCondition )
        cout << "integral condition on the flux and ";
    if ( M_isPicard )
        cout << "Picard algorithm" << std::endl;
    else
        cout << "linear case" << std::endl;

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    M_Vh = Pdhv<Order>( M_mesh, true );
    M_Wh = Pdh<Order>( M_mesh, true );
    M_Mh = Pdh<Order>( face_mesh, true );
    M_Ch = Pch<0>( M_mesh );
    M_M0h = Pdh<0>( face_mesh, true );

    toc("spaces",true);

    cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;

    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_up = M_Vh->elementPtr( "u" );
    M_pp = M_Wh->elementPtr( "p" );

    tic();
    if ( M_integralCondition )
        this->initGraphsWithIntegralCond();
    else
        this->initGraphs();
    toc("graphs",true);

    tic();
    assembleA();
    toc("assemble A", true);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::initGraphs()
{
    auto phatp = M_Mh->elementPtr( "phat" );

    BlocksBaseGraphCSR hdg_graph(3,3);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();

    M_A = M_backend->newBlockMatrix(_block=hdg_graph);
    M_A_cst = M_backend->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(3);
    hdg_vec(0,0) = M_backend->newVector( M_Vh );
    hdg_vec(1,0) = M_backend->newVector( M_Wh );
    hdg_vec(2,0) = M_backend->newVector( M_Mh );
    M_F = M_backend->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(3);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_U = M_backend->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::initGraphsWithIntegralCond()
{
    auto phatp = M_Mh->elementPtr( "phat" );
    auto mup = M_Ch->elementPtr( "c1" );

    BlocksBaseGraphCSR hdg_graph(4,4);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,0) = stencil( _test=M_Ch,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,1) = stencil( _test=M_Ch,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,2) = stencil( _test=M_Ch,_trial=M_Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,3) = stencil( _test=M_Vh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,3) = stencil( _test=M_Wh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,3) = stencil( _test=M_Mh,_trial=M_Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(3,3) = stencil( _test=M_Ch,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();

    M_A = M_backend->newBlockMatrix(_block=hdg_graph);
    M_A_cst = M_backend->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(4);
    hdg_vec(0,0) = M_backend->newVector( M_Vh );
    hdg_vec(1,0) = M_backend->newVector( M_Wh );
    hdg_vec(2,0) = M_backend->newVector( M_Mh );
    hdg_vec(3,0) = M_backend->newVector( M_Ch );
    M_F = M_backend->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(4);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_hdg_sol(3,0) = mup;
    M_U = M_backend->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::solve()
{
    tic();
    updateConductivityTerm();
    assembleF();

    if ( M_updateAssembly != NULL )
        this->M_updateAssembly( M_A, M_F );
    M_backend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U);
    M_hdg_sol.localize(M_U);
    toc("solve", true);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::solveNL()
{
    tic();
    if ( M_updateAssembly != NULL )
        this->M_updateAssembly( M_A, M_F );
    M_backend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U);
    M_hdg_sol.localize(M_U);
    toc("solve", true);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::assembleA()
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto q = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );
    auto H = M_M0h->element( "H" );
    if ( ioption(prefixvm(M_prefix, "hface") ) == 0 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 1 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 2 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
        H.on( _range=elements(M_M0h->mesh()), _expr=h() );
    // stabilisation parameter
    auto tau_constant = cst(doption(prefixvm(M_prefix, "tau_constant")));

    auto gammaMinusIntegral = complement(boundaryfaces(M_mesh),[this]( auto const& e ) {
            for( auto marker : this->M_integralMarkersList)
            {
                if ( e.marker().value() == this->M_mesh->markerName( marker ) )
                    return true;
            }
            return false; });

    auto a12 = form2( _trial=M_Wh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0, _colstart=1 );
    auto a13 = form2( _trial=M_Mh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0, _colstart=2);
    auto a21 = form2( _trial=M_Vh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=0);
    auto a22 = form2( _trial=M_Wh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=1 );
    auto a23 = form2( _trial=M_Mh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=2);
    auto a31 = form2( _trial=M_Vh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2, _colstart=0);
    auto a32 = form2( _trial=M_Wh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2, _colstart=1);
    auto a33 = form2(_trial=M_Mh, _test=M_Mh,_matrix=M_A_cst,
                     _rowstart=2, _colstart=2);

    // -(p,div(v))_Omega
    a12 += integrate(_range=elements(M_mesh),_expr=-(idt(p)*div(v)));

    // <phat,v.n>_Gamma\Gamma_I
    a13 += integrate(_range=internalfaces(M_mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=gammaMinusIntegral,
                     _expr=idt(phat)*trans(id(v))*N());


    // -(j, grad(w))
    a21 += integrate(_range=elements(M_mesh),_expr=(-grad(w)*idt(u)));
    // <j.n,w>_Gamma
    a21 += integrate(_range=internalfaces(M_mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    a21 += integrate(_range=internalfaces(M_mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    a21 += integrate(_range=boundaryfaces(M_mesh),
                     _expr=(id(w)*trans(idt(u))*N()));


    // <tau p, w>_Gamma
    a22 += integrate(_range=internalfaces(M_mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(idv(H),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(idv(H),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(M_mesh),
                     _expr=(tau_constant * pow(idv(H),M_tau_order)*id(w)*idt(p)));


    // <-tau phat, w>_Gamma\Gamma_I
    a23 += integrate(_range=internalfaces(M_mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(idv(H),M_tau_order)*id(w) )+
                       rightface( pow(idv(H),M_tau_order)*id(w) )));
    a23 += integrate(_range=gammaMinusIntegral,
                     _expr=-tau_constant * idt(phat) * pow(idv(H),M_tau_order)*id(w) );


    // <j.n,mu>_Omega/Gamma
    a31 += integrate(_range=internalfaces(M_mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );


    // <tau p, mu>_Gamma_N
    a32 += integrate(_range=internalfaces(M_mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(idv(H),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(idv(H),M_tau_order)*idt(p) )));


    // <-tau phat, mu>_Omega/Gamma
    a33 += integrate(_range=internalfaces(M_mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(idv(H),M_tau_order) )+
                                                                 rightface( pow(idv(H),M_tau_order) )));

    // BC
    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <phat, mu>_Gamma_D
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=idt(phat) * id(l) );
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <j.n,mu>_Gamma_N
                a31 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));
                // <tau p, mu>_Gamma_N
                a32 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(idv(H),M_tau_order)*idt(p) ) );
                // <-tau phat, mu>_Gamma_N
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(idv(H),M_tau_order) ) );
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                // <j.n,mu>_Gamma_R
                a31 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));
                // <tau p, mu>_Gamma_R
                a32 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(idv(H),M_tau_order)*idt(p) ) );
                // <-tau phat, mu>_Gamma_R
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(idv(H),M_tau_order) ) );
                // <g_R^1 phat, mu>_Gamma_R
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=g*idt(phat) * id(l) );
            }
        }
    }

    if ( M_integralCondition )
    {
        auto a14 = form2(_trial=M_Ch, _test=M_Vh,_matrix=M_A_cst,
                         _rowstart=0,
                         _colstart=3);
        auto a34 = form2(_trial=M_Ch, _test=M_Mh,_matrix=M_A_cst,
                         _rowstart=2,
                         _colstart=3);
        auto a41 = form2(_trial=M_Vh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=0);
        auto a42 = form2(_trial=M_Wh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=1);
        auto a43 = form2(_trial=M_Mh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=2);
        auto a24 = form2(_trial=M_Ch, _test=M_Wh,_matrix=M_A_cst,
                         _rowstart=1,
                         _colstart=3);
        auto a44 = form2(_trial=M_Ch, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=3);

        auto itField = M_modelProperties->boundaryConditions().find( "potential");
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Dirichlet" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    // <lambda, v.n>_Gamma_I
                    a14 += integrate( _range=markedfaces(M_mesh,marker),
                                      _expr=trans(id(u))*N()*idt(nu) );

                    // <lambda, tau w>_Gamma_I
                    a24 += integrate( _range=markedfaces(M_mesh,marker),
                                      _expr=-tau_constant * ( pow(idv(H),M_tau_order)*id(w) ) * idt(nu) );

                    // <j.n, m>_Gamma_I
                    a41 += integrate( _range=markedfaces(M_mesh,marker), _expr=trans(idt(u))*N()*id(nu) );

                    // <tau p, m>_Gamma_I
                    a42 += integrate( _range=markedfaces(M_mesh,marker), _expr=tau_constant *
                          ( pow(idv(H),M_tau_order)*idt(p) ) * id(nu) );

                    // -<lambda2, m>_Gamma_I
                    a44 += integrate( _range=markedfaces(M_mesh,marker), _expr=-pow(idv(H),M_tau_order)*id(nu)*idt(nu) );
                }
            }
        }
    }
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::assembleF()
{
    M_F->zero();

    auto nu = M_Ch->element( "nu" );
    auto l = M_Mh->element( "lambda" );
    auto w = M_Wh->element();
    auto v = M_Vh->element();

    // Building the RHS

    auto rhs1 = form1( _test=M_Wh, _vector=M_F, _rowstart=0);
    auto rhs2 = form1( _test=M_Wh, _vector=M_F, _rowstart=1);
    auto rhs3 = form1( _test=M_Mh, _vector=M_F, _rowstart=2);

    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // (f, w)_Omega
                rhs2 += integrate( _range=markedelements(M_mesh,marker),
                                   _expr=g*id(w));
            }
        }
        itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // <g_D, mu>_Gamma_D
                rhs3 += integrate(_range=markedfaces(M_mesh,marker),
                                  _expr=id(l)*g);
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // <g_N,mu>_Gamma_N
                rhs3 += integrate( _range=markedfaces(M_mesh, marker),
                                   _expr=id(l)*g);
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression2());
                // <g_R^2,mu>_Gamma_R
                rhs3 += integrate( _range=markedfaces(M_mesh, marker),
                                   _expr=id(l)*g);
            }
        }
    }

    itField = M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<3,1>(exAtMarker.expression());
                // (g, v)_Omega
                rhs2 += integrate( _range=markedelements(M_mesh,marker),
                                   _expr=inner(g,id(v)));
            }
        }
    }

    if ( M_integralCondition )
    {
        auto rhs4 = form1( _test=M_Ch, _vector=M_F,
                           _rowstart=3);
        itField = M_modelProperties->boundaryConditions().find( "flux");
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    double meas = integrate( _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
                    auto g = expr(exAtMarker.expression());
                    // <I_target,m>_Gamma_I
                    rhs4 += integrate(_range=markedfaces(M_mesh,marker),
                                      _expr=g*id(nu)/meas);
                }
            }
        }
    }
}

template<int Dim, int Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order>::updateConductivityTerm( Expr<ExprT> expr, std::string marker)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );

    M_A->zero();
    M_A->addMatrix(1.,M_A_cst);

    auto a11 = form2( _trial=M_Vh, _test=M_Vh,_matrix=M_A );
    if ( marker.empty() )
        a11 += integrate( _range=elements(M_mesh), _expr=inner(idt(u),id(v))/expr);
    else
        a11 += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/expr);
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::updateConductivityTerm( bool isNL)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );

    M_A->zero();
    M_A->addMatrix(1.,M_A_cst);

    auto a11 = form2( _trial=M_Vh, _test=M_Vh, _matrix=M_A );
    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( !isNL )
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivity_json")));
            // (sigma^-1 j, v)
            a11 += integrate(_range=markedelements(M_mesh,marker),
                             _expr=(trans(idt(u))*id(v))/cond );
        }
        else
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivityNL_json")), "p", idv(*M_pp));
            // (sigma(p)^-1 j, v)
            a11 += integrate(_range=markedelements(M_mesh,marker),
                             _expr=(trans(idt(u))*id(v))/cond );
        }
    }
}

template<int Dim, int Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order>::updatePotentialRHS( Expr<ExprT> expr, std::string marker)
{
    auto rhs = form1( _test=M_Wh, _vector=M_F, _colstart=1);
    auto w = M_Wh->element();
    if ( marker.empty() )
        rhs += integrate(_range=elements(M_mesh), _expr=expr*id(w) );
    else
        rhs += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(w)) );
}

template<int Dim, int Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order>::updateFluxRHS( Expr<ExprT> expr, std::string marker)
{
    auto rhs = form1( _test=M_Vh, _vector=M_F, _colstart=0);
    auto v = M_Vh->element();
    if ( marker.empty() )
        rhs += integrate(_range=elements(M_mesh), _expr=expr*id(v) );
    else
        rhs += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(v)) );
}

template<int Dim, int Order>
void
MixedPoisson<Dim, Order>::exportResults()
{
    auto e = exporter( M_mesh);

    auto postProcess = M_modelProperties->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "flux" )
                e->add(prefixvm(M_prefix, "flux"), *M_up);
            if ( field == "potential" )
                e->add(prefixvm(M_prefix, "potential"), *M_pp);
        }
    }

    e->save();
}

} // Namespace Feel

