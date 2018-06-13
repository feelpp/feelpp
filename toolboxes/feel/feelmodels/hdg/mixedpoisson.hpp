#ifndef _MIXEDPOISSON2_HPP
#define _MIXEDPOISSON2_HPP

#include <boost/algorithm/string.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/complement.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/hdg/options.hpp>

#define USE_SAME_MAT 1

namespace Feel
{

namespace FeelModels
{

template <int Dim, int Order, int G_Order = 1, int E_Order = 4>
class MixedPoisson : public ModelNumerical
{
  public:
    typedef ModelNumerical super_type;

    static const uint16_type expr_order = Order + E_Order;
    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype;

    typedef MixedPoisson<Dim, Order, G_Order, E_Order> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim, G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim - 1, G_Order, Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;

    // Vh
    using Vh_t = Pdhv_type<mesh_type, Order>;
    using Vh_ptr_t = Pdhv_ptrtype<mesh_type, Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
    // Wh
    using Wh_t = Pdh_type<mesh_type, Order>;
    using Wh_ptr_t = Pdh_ptrtype<mesh_type, Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
    // Mh
    using Mh_t = Pdh_type<face_mesh_type, Order>;
    using Mh_ptr_t = Pdh_ptrtype<face_mesh_type, Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
    // Ch
    using Ch_t = Pch_type<face_mesh_type, 0>;
    using Ch_ptr_t = Pch_ptrtype<face_mesh_type, 0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    using Ch_element_vector_type = std::vector<Ch_element_t>;
    // M0h
    using M0h_t = Pdh_type<face_mesh_type, 0>;
    using M0h_ptr_t = Pdh_ptrtype<face_mesh_type, 0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    //using product_space_t = boost::shared_ptr<ProductSpace<Ch_ptr_t,true>>;
    using product2_space_type = ProductSpaces2<Ch_ptr_t, Vh_ptr_t, Wh_ptr_t, Mh_ptr_t>;
    using product2_space_ptrtype = boost::shared_ptr<product2_space_type>;
    using block_bilinear_type = BlockBilinearForm<product2_space_type>;
    using block_linear_type = BlockLinearForm<product2_space_type>;

    typedef Exporter<mesh_type, G_Order> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

    using op_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Wh_t, Pdh_type<mesh_type, Order>>>;
    using opv_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Vh_t, Pdhv_type<mesh_type, Order>>>;

    using integral_boundary_list_type = std::vector<ExpressionStringAtMarker>;

    // time
    // typedef Lagrange<Order,Scalar,Discontinuous> basis_scalar_type;
    // typedef FunctionSpace<mesh_type, bases<basis_scalar_type>> space_mixedpoisson_type;

    // typedef Bdf<space_mixedpoisson_type>  bdf_type;
    typedef Bdf<Wh_t> bdf_type;
    // typedef Bdf<Pdh_type<mesh_type,Order>> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    //private:
  protected:
    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh; // potential trace
    Ch_ptr_t M_Ch; // Lagrange multiplier
    M0h_ptr_t M_M0h;
    product2_space_type M_ps;

    backend_ptrtype M_backend;
    block_bilinear_type M_a;
    block_linear_type M_rhs;
    vector_ptrtype M_U;

    Vh_element_t M_up;            // flux solution
    Wh_element_t M_pp;            // potential solution
    Ch_element_vector_type M_mup; // potential solution on the integral boundary conditions

    // time discretization
    bdf_ptrtype M_bdf_mixedpoisson;

    int M_tau_order;

    int M_integralCondition;
    int M_useUserIBC;
    integral_boundary_list_type M_IBCList;

    bool M_isPicard;

  public:
    // constructor
    MixedPoisson( std::string const& prefix = "mixedpoisson",
                  WorldComm const& _worldComm = Environment::worldComm(),
                  std::string const& subPrefix = "",
	              ModelBaseRepository const& modelRep = ModelBaseRepository() );

    MixedPoisson( self_type const& MP ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedpoisson",
                             WorldComm const& worldComm = Environment::worldComm(),
                             std::string const& subPrefix = "",
							 ModelBaseRepository const& modelRep = ModelBaseRepository() );

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    M0h_ptr_t traceSpaceOrder0() const { return M_M0h; }
    Ch_ptr_t constantSpace() const { return M_Ch; }

    Vh_element_t fluxField() const { return M_up; }
    Wh_element_t potentialField() const { return M_pp; }
    integral_boundary_list_type integralBoundaryList() const { return M_IBCList; }
    int integralCondition() const { return M_integralCondition; }
    void setIBCList( std::vector<std::string> markersIbc );
    int tau_order() const { return M_tau_order; }
    backend_ptrtype get_backend() { return M_backend; }
    product2_space_type getPS() const { return M_ps; }
	block_bilinear_type get_a() { return M_a; }
    block_linear_type get_rhs() { return M_rhs; }


    // time step scheme
    virtual void createTimeDiscretization();
    bdf_ptrtype timeStepBDF() { return M_bdf_mixedpoisson; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf_mixedpoisson; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepBDF(); }

    // Exporter
    virtual void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr )
    {
        this->exportResults( this->currentTime(), mesh );
        M_exporter->save();
    }
    void exportResults( double Time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr );
    exporter_ptrtype M_exporter;
    exporter_ptrtype exporterMP() { return M_exporter; }

    void init( mesh_ptrtype mesh = nullptr, mesh_ptrtype meshVisu = nullptr );

    virtual void initModel();
    virtual void initSpaces();
    virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void assembleAll();
    virtual void assembleCstPart();
    virtual void assembleNonCstPart();
    void setZero();

    void assembleRHS();
    template <typename ExprT>
    void updateConductivityTerm( Expr<ExprT> expr, std::string marker = "" );
    void updateConductivityTerm( bool isNL = false );

    template <typename ExprT>
    void assembleFluxRHS( Expr<ExprT> expr, std::string marker );
    template <typename ExprT>
    void assemblePotentialRHS( Expr<ExprT> expr, std::string marker );

    void assembleBoundaryCond();
    void assembleRhsBoundaryCond();
    void assembleDirichlet( std::string marker );
    void assembleNeumann( std::string marker );
    template <typename ExprT>
    void assembleRobin( Expr<ExprT> expr1, std::string marker );
    template <typename ExprT>
    void assembleRhsDirichlet( Expr<ExprT> expr, std::string marker );
    template <typename ExprT>
    void assembleRhsNeumann( Expr<ExprT> expr, std::string marker );
    template <typename ExprT>
    void assembleRhsInterfaceCondition( Expr<ExprT> expr, std::string marker );
    template <typename ExprT>
    void assembleRhsRobin( Expr<ExprT> expr2, std::string marker );

    void assembleIBC( int i, std::string marker = "" );
    void assembleRhsIBC( int i, std::string marker = "", double intjn = 0 );

    virtual void solve();
};

#if !defined( FEELPP_MODELS_HDG_NOEXTERN )
extern template class MixedPoisson<2, 0, 1, 4>;
extern template class MixedPoisson<2, 1, 1, 4>;
extern template class MixedPoisson<2, 2, 1, 4>;
extern template class MixedPoisson<2, 3, 1, 4>;

extern template class MixedPoisson<3, 0, 1, 4>;
extern template class MixedPoisson<3, 1, 1, 4>;
extern template class MixedPoisson<3, 2, 1, 4>;
extern template class MixedPoisson<3, 3, 1, 4>;
#endif

template <int Dim, int Order, int G_Order, int E_Order>
MixedPoisson<Dim, Order, G_Order, E_Order>::MixedPoisson( std::string const& prefix,
                                                          WorldComm const& worldComm,
                                                          std::string const& subPrefix,
                                                          ModelBaseRepository const& modelRep )
    : super_type( prefix, worldComm, subPrefix, modelRep ),
      M_useUserIBC( false )
{
    if ( this->verbose() ) Feel::FeelModels::Log( this->prefix() + ".MixedPoisson", "constructor", "start",
                                                  this->worldComm(), this->verboseAllProc() );

    this->setFilenameSaveInfo( prefixvm( this->prefix(), "MixedPoisson.info" ) );

    if ( this->prefix().empty() )
        M_backend = backend( _rebuild = true );
    else
        M_backend = backend( _name = this->prefix(), _rebuild = true );

    M_tau_order = ioption( prefixvm( this->prefix(), "tau_order" ) );

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonSolve.data";
    this->addTimerTool( "Constructor", nameFileConstructor );
    this->addTimerTool( "Solve", nameFileSolve );

    if ( this->verbose() ) Feel::FeelModels::Log( this->prefix() + ".MixedPoisson", "constructor", "finish",
                                                  this->worldComm(), this->verboseAllProc() );
}

template <int Dim, int Order, int G_Order, int E_Order>
typename MixedPoisson<Dim, Order, G_Order, E_Order>::self_ptrtype
MixedPoisson<Dim, Order, G_Order, E_Order>::New( std::string const& prefix,
                                                 WorldComm const& worldComm, std::string const& subPrefix,
                                                 ModelBaseRepository const& modelRep )
{
    return boost::make_shared<self_type>( prefix, worldComm, subPrefix, modelRep );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::setIBCList( std::vector<std::string> markersIBC )
{
    M_useUserIBC = true;
    M_IBCList.clear();
    for ( auto const& marker : markersIBC )
    {
        ExpressionStringAtMarker exAtMark( std::make_tuple( "expression", marker, std::string( "" ), std::string( "" ), std::string( "" ) ) );
        M_IBCList.push_back( exAtMark );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::init( mesh_ptrtype mesh, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type );
    else
        M_mesh = mesh;
    toc( "mesh" );

    tic();
    this->initModel();
    toc( "model" );

    tic();
    this->initSpaces();
    toc( "spaces" );

    if ( !this->isStationary() )
    {
        tic();
        this->createTimeDiscretization();
        this->initTimeStep();
        toc( "timeDiscretization", true );
    }

    tic();
    this->initExporter( meshVisu );
    toc( "exporter" );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::initModel()
{

    // initialize marker lists for each boundary condition type
    auto itField = modelProperties().boundaryConditions().find( "potential" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "Dirichlet" );

        if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker( marker ) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl
                               << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann:";
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker( marker ) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl
                               << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann computed from exact pressure:";
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker( marker ) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl
                               << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Robin:";
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker( marker ) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl
                               << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
    }
    if ( !M_useUserIBC )
    {
        M_IBCList.clear();
        itField = modelProperties().boundaryConditions().find( "flux" );
        if ( itField != modelProperties().boundaryConditions().end() )
        {
            auto mapField = ( *itField ).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                Feel::cout << "Integral:";
                for ( auto const& exAtMarker : ( *itType ).second )
                {
                    std::string marker = exAtMarker.marker();
                    if ( M_mesh->hasFaceMarker( marker ) )
                        Feel::cout << " " << marker;
                    else
                        Feel::cout << std::endl
                                   << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                    M_IBCList.push_back( exAtMarker );
                }
                Feel::cout << std::endl;
            }

            itType = mapField.find( "InterfaceCondition" );
            if ( itType != mapField.end() )
            {
                Feel::cout << "Interface condition:";
                for ( auto const& exAtMarker : ( *itType ).second )
                {
                    std::string marker = exAtMarker.marker();
                    if ( M_mesh->hasFaceMarker( marker ) )
                        Feel::cout << " " << marker;
                    else
                        Feel::cout << std::endl
                                   << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                }
                Feel::cout << std::endl;
            }
        }
    }

    if ( M_IBCList.empty() )
        M_integralCondition = 0;
    else
        M_integralCondition = M_IBCList.size();

	if ( boost::icontains(modelProperties().models().model().equations(),"picard") )
        M_isPicard = true;
    else
        M_isPicard = false;
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::initSpaces()
{
    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement( faces( M_mesh ), [this]( auto const& ewrap ) {
            auto const& e = unwrap_ref( ewrap );
            for( auto exAtMarker : this->M_IBCList)
            {
                if ( e.hasMarker() && e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                    return true;
            }
            return false; } );

    auto face_mesh = createSubmesh( M_mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );

    M_Vh = Pdhv<Order>( M_mesh, true );
    M_Wh = Pdh<Order>( M_mesh, true );
    M_Mh = Pdh<Order>( face_mesh, true );
    // M_Ch = Pch<0>( M_mesh, true );
    M_M0h = Pdh<0>( face_mesh );

    std::vector<std::string> ibc_markers( M_integralCondition );
    for ( int i = 0; i < M_integralCondition; i++ )
    {
        ibc_markers.push_back( M_IBCList[i].marker() );
    }

    auto ibc_mesh = createSubmesh( M_mesh, markedfaces( M_mesh, ibc_markers ), EXTRACTION_KEEP_MESH_RELATION, 0 );
    M_Ch = Pch<0>( ibc_mesh, true );

    Feel::cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
               << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
               << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        Feel::cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;

    auto ibcSpaces = boost::make_shared<ProductSpace<Ch_ptr_t, true>>( M_integralCondition, M_Ch );
    M_ps = product2( ibcSpaces, M_Vh, M_Wh, M_Mh );

    M_up = M_Vh->element( "u" );
    M_pp = M_Wh->element( "p" );

    for ( int i = 0; i < M_integralCondition; i++ )
        M_mup.push_back( M_Ch->element( "mup" ) );

    solve::strategy s = boption( prefixvm( prefix(), "use-sc" ) ) ? solve::strategy::static_condensation : solve::strategy::monolithic;

    M_a.setFunctionSpace( M_ps );
    M_a.setStrategy( s, M_backend );
    M_rhs.setFunctionSpace( M_ps );
    M_rhs.setStrategy( s, M_backend );
    //M_a = blockform2( M_ps, s, M_backend );
    //M_rhs = blockform1( M_ps, s, M_backend );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::initExporter( mesh_ptrtype meshVisu )
{
    std::string geoExportType = "static"; //change_coords_only, change, static
    M_exporter = exporter( _mesh = meshVisu ? meshVisu : this->mesh(),
                           _name = "Export",
                           _geo = geoExportType,
                           _path = this->exporterPath() );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::solve()
{
    tic();

    auto U = M_ps.element();

    tic();
    Feel::cout << "Start solving" << std::endl;
    M_a.solve( _solution = U, _rhs = M_rhs, _condense = boption( prefixvm( prefix(), "use-sc" ) ), _name=prefix() );
    toc( "MixedPoisson : static condensation" );

    toc( "solve" );

    M_up = U( 0_c );
    M_pp = U( 1_c );

    for ( int i = 0; i < M_integralCondition; i++ )
        M_mup[i] = U( 3_c, i );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::setZero()
{
    M_a.zero();
    M_rhs.zero();
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleAll()
{
    this->setZero();

	this->assembleCstPart();
    this->assembleNonCstPart();
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleCstPart()
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto q = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto uI = M_Ch->element( "uI" );

    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );

    // stabilisation parameter
    auto tau_constant = cst( doption( prefixvm( prefix(), "tau_constant" ) ) );

    auto sc_param = boption( prefixvm( prefix(), "use-sc" ) ) ? 0.5 : 1.0;

    auto gammaMinusIntegral = complement( boundaryfaces( M_mesh ),
                                          [this]( auto const& ewrap ) {
                                              auto const& e = unwrap_ref( ewrap );
                                              for ( auto exAtMarker : this->M_IBCList )
                                              {
                                                  if ( e.hasMarker() && e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                                                      return true;
                                              }
                                              return false;
                                          } );

    // -(p,div(v))_Omega
    M_a( 0_c, 1_c ) += integrate( _range = elements( M_mesh ), _expr = -( idt( p ) * div( v ) ) );

    // <phat,v.n>_Gamma\Gamma_I
    M_a( 0_c, 2_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = idt( phat ) * ( leftface( normal( v ) ) + rightface( normal( v ) ) ) );
    M_a( 0_c, 2_c ) += integrate( _range = gammaMinusIntegral,
                                  _expr = idt( phat ) * normal( v ) );

    // (div(j),q)_Omega
    M_a( 1_c, 0_c ) += integrate( _range = elements( M_mesh ), _expr = -( id( w ) * divt( u ) ) );

    // <tau p, w>_Gamma
    M_a( 1_c, 1_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = -tau_constant *
                                          ( leftfacet( idt( p ) ) * leftface( id( w ) ) +
                                            rightfacet( idt( p ) ) * rightface( id( w ) ) ) );
    M_a( 1_c, 1_c ) += integrate( _range = boundaryfaces( M_mesh ),
                                  _expr = -( tau_constant * id( w ) * idt( p ) ) );

    // <-tau phat, w>_Gamma\Gamma_I
    M_a( 1_c, 2_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = tau_constant * idt( phat ) *
                                          ( leftface( id( w ) ) +
                                            rightface( id( w ) ) ) );
    M_a( 1_c, 2_c ) += integrate( _range = gammaMinusIntegral,
                                  _expr = tau_constant * idt( phat ) * id( w ) );

    // <j.n,mu>_Omega/Gamma
    M_a( 2_c, 0_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = id( l ) * ( leftfacet( normalt( u ) ) + rightfacet( normalt( u ) ) ) );

    // <tau p, mu>_Omega/Gamma
    M_a( 2_c, 1_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = tau_constant * id( l ) * ( leftfacet( idt( p ) ) + rightfacet( idt( p ) ) ) );

    // <-tau phat, mu>_Omega/Gamma
    M_a( 2_c, 2_c ) += integrate( _range = internalfaces( M_mesh ),
                                  _expr = -sc_param * tau_constant * idt( phat ) * id( l ) );

    this->assembleBoundaryCond();
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleNonCstPart()
{
    modelProperties().parameters().updateParameterValues();

    this->updateConductivityTerm();
    this->assembleRHS();
    this->assembleRhsBoundaryCond();
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::updateConductivityTerm( Expr<ExprT> expr, std::string marker )
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );

    if ( marker.empty() )
        M_a( 0_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = elements( M_mesh ), _expr = inner( idt( u ), id( v ) ) / expr );
    else
        M_a( 0_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_mesh, marker ), _expr = inner( idt( u ), id( v ) ) / expr );

    // (1/delta_t p, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() )
    {
        M_a( 1_c, 1_c ) += integrate( _range = elements( M_mesh ),
                                      _expr = -( this->timeStepBDF()->polyDerivCoefficient( 0 ) * idt( p ) * id( w ) ) );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::updateConductivityTerm( bool isNL )
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );

    for ( auto const& pairMat : modelProperties().materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( !isNL )
        {
            auto cond = material.getScalar( soption( prefixvm( prefix(), "conductivity_json" ) ) );
            // (sigma^-1 j, v)
            M_a( 0_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_mesh, marker ),
                                          _expr = ( trans( idt( u ) ) * id( v ) ) / cond );
        }
        else
        {
            auto cond = material.getScalar( soption( prefixvm( prefix(), "conductivityNL_json" ) ), "p", idv( M_pp ) );
            // (sigma(p)^-1 j, v)
            M_a( 0_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_mesh, marker ),
                                          _expr = ( trans( idt( u ) ) * id( v ) ) / cond );
        }
    }

    // (1/delta_t p, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() )
    {
        M_a( 1_c, 1_c ) += integrate( _range = elements( M_mesh ),
                                      _expr = -( this->timeStepBDF()->polyDerivCoefficient( 0 ) * idt( p ) * id( w ) ) );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRHS()
{
    // (p_old,w)_Omega
    if ( !this->isStationary() )
    {
        auto bdf_poly = M_bdf_mixedpoisson->polyDeriv();
        this->assemblePotentialRHS( idv( bdf_poly ), "" );
    }

    auto itField = modelProperties().boundaryConditions().find( "potential" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<expr_order>( exAtMarker.expression() );
                if ( !this->isStationary() )
                    g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                this->assemblePotentialRHS( g, marker );
            }
        }
    }

    itField = modelProperties().boundaryConditions().find( "flux" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<Dim, 1, expr_order>( exAtMarker.expression() );
                if ( !this->isStationary() )
                    g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                this->assembleFluxRHS( g, marker );
            }
        }
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleFluxRHS( Expr<ExprT> expr, std::string marker )
{
    auto v = M_Vh->element();

    if ( marker.empty() )
        M_rhs( 0_c ) += integrate( _quad = _Q<expr_order>(), _range = elements( M_mesh ),
                                   _expr = inner( expr, id( v ) ) );
    else
        M_rhs( 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_mesh, marker ),
                                   _expr = inner( expr, id( v ) ) );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assemblePotentialRHS( Expr<ExprT> expr, std::string marker )
{
    auto w = M_Wh->element();

    if ( marker.empty() )
    {
        M_rhs( 1_c ) += integrate( _quad = _Q<expr_order>(), _range = elements( M_mesh ),
                                   _expr = -inner( expr, id( w ) ) );
    }
    else
    {
        M_rhs( 1_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_mesh, marker ),
                                   _expr = -inner( expr, id( w ) ) );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleBoundaryCond()
{
    auto itField = modelProperties().boundaryConditions().find( "potential" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                this->assembleDirichlet( marker );
            }
        }

        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                this->assembleNeumann( marker );
            }
        }

        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                this->assembleNeumann( marker );
            }
        }

        // Robin depends on expression, so it is in Rhs
    }

    for ( int i = 0; i < M_integralCondition; i++ )
        this->assembleIBC( i );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsBoundaryCond()
{

    auto itField = modelProperties().boundaryConditions().find( "potential" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( exAtMarker.isExpression() )
                {
                    auto g = expr<expr_order>( exAtMarker.expression() );
                    if ( !this->isStationary() )
                        g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                    this->assembleRhsDirichlet( g, marker );
                }
                else if ( exAtMarker.isFile() )
                {
                    double g = 0;
                    if ( !this->isStationary() )
                    {
                        Feel::cout << "use data file to set rhs for Dirichlet BC at time " << M_bdf_mixedpoisson->time() << std::endl;
                        LOG( INFO ) << "use data file to set rhs for Dirichlet BC at time " << M_bdf_mixedpoisson->time() << std::endl;

                        // data may depend on time
                        g = exAtMarker.data( M_bdf_mixedpoisson->time() );
                    }
                    else
                        g = exAtMarker.data( 0.1 );

                    LOG( INFO ) << "use g=" << g << std::endl;
                    Feel::cout << "g=" << g << std::endl;

                    this->assembleRhsDirichlet( cst( g ), marker );
                }
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                std::string exprString = exAtMarker.expression();
                int nComp = nbComp( exprString );
                if ( nComp == 1 )
                {
                    auto g = expr<expr_order>( exAtMarker.expression() );
                    if ( !this->isStationary() )
                        g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                    this->assembleRhsNeumann( g, marker );
                }
                else if ( nComp == Dim )
                {
                    auto g = expr<Dim, 1, expr_order>( exAtMarker.expression() );
                    if ( !this->isStationary() )
                        g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                    /*
	                auto M_rhs = blockform1( *M_ps, M_F );
                    auto l = M_Mh->element( "lambda" );

	                // <g_N,mu>_Gamma_N
                    M_rhs(2_c) += integrate(_range=markedfaces(M_mesh, marker),
                                          _expr=trans(g)*N() * id(l));
                    */
                    auto gn = inner( g, N() );
                    this->assembleRhsNeumann( gn, marker );
                }
            }
        }

        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                auto p_ex = expr<expr_order>( exAtMarker.expression() );
                auto gradp_ex = expr( grad<Dim>( p_ex ) );
                if ( !this->isStationary() )
                    gradp_ex.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );

                for ( auto const& pairMat : modelProperties().materials() )
                {
                    auto material = pairMat.second;
                    auto K = material.getDouble( "k" );

                    auto g = expr( -K * trans( gradp_ex ) );
                    auto gn = inner( g, N() );
                    this->assembleRhsNeumann( gn, marker );
                }
            }
        }

        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr<expr_order>( exAtMarker.expression1() );
                auto g2 = expr<expr_order>( exAtMarker.expression2() );

                if ( !this->isStationary() )
                {
                    g1.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                    g2.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                }

                // both depend on expression
                this->assembleRobin( g1, marker );
                this->assembleRhsRobin( g2, marker );
            }
        }
    }

    itField = modelProperties().boundaryConditions().find( "flux" );
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "InterfaceCondition" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<1, 1, expr_order>( exAtMarker.expression() );
                if ( !this->isStationary() )
                    g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
                Feel::cout << "Interface condition on " << marker << ":\t" << g << std::endl;
                this->assembleRhsInterfaceCondition( g, marker );
            }
        }
    }

    for ( int i = 0; i < M_integralCondition; i++ )
        this->assembleRhsIBC( i );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsDirichlet( Expr<ExprT> expr, std::string marker )
{
    auto l = M_Mh->element( "lambda" );

    // <g_D, mu>_Gamma_D
    M_rhs( 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                               _expr = id( l ) * expr );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsNeumann( Expr<ExprT> expr, std::string marker )
{

    auto l = M_Mh->element( "lambda" );

    // <g_N,mu>_Gamma_N
    M_rhs( 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                               _expr = id( l ) * expr );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsInterfaceCondition( Expr<ExprT> expr, std::string marker )
{
    auto l = M_Mh->element( "lambda" );

    // <g_interface,mu>_Gamma_N
    M_rhs( 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedelements( M_Mh->mesh(), marker ),
                               _expr = id( l ) * expr );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsRobin( Expr<ExprT> expr2, std::string marker )
{
    auto l = M_Mh->element( "lambda" );
    // <g_R2,mu>_Gamma_R
    M_rhs( 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                               _expr = id( l ) * expr2 );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsIBC( int i, std::string markerOpt, double intjn )
{
    auto u = M_Vh->element( "u" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto uI = M_Ch->element( "uI" );

    std::string marker;
    Expr<GinacEx<expr_order>> g;

    if ( !markerOpt.empty() )
    {
        marker = markerOpt;
        std::ostringstream f;
        f << std::setprecision( 14 ) << intjn;
        g = expr<expr_order>( f.str() );
    }
    else
    {
        auto exAtMarker = M_IBCList[i];
        marker = exAtMarker.marker();
        if ( exAtMarker.isExpression() )
        {
            g = expr<expr_order>( exAtMarker.expression() );
            if ( !this->isStationary() )
                g.setParameterValues( {{"t", M_bdf_mixedpoisson->time()}} );
        }
        else if ( exAtMarker.isFile() )
        {
            double d = 0;
            if ( !this->isStationary() )
            {
                Feel::cout << "use data file to set rhs for IBC at time " << M_bdf_mixedpoisson->time() << std::endl;
                LOG( INFO ) << "use data file to set rhs for IBC at time " << M_bdf_mixedpoisson->time() << std::endl;

                // data may depend on time
                d = exAtMarker.data( M_bdf_mixedpoisson->time() );
            }
            else
                d = exAtMarker.data( 0.1 );

            // Scale entries if necessary
    		{
                for ( auto const& field : modelProperties().postProcess().exports().fields() )
        		{
            		if ( field == "scaled_flux" )
					{
            			for( auto const& pairMat : modelProperties().materials() )
            			{
                			auto material = pairMat.second;
                			double kk = material.getDouble( "scale_integral_file" );
							d = kk*d;
						}
					}
				}
			}
            LOG( INFO ) << "use g=" << d << std::endl;
            Feel::cout << "g=" << d << std::endl;

            std::ostringstream f;
            f << std::setprecision( 14 ) << d;
            g = expr<expr_order>( f.str() );
        }
    }

    double meas = integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );

    // <I_target,m>_Gamma_I
    M_rhs( 3_c, i ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ), _expr = g * id( nu ) / meas );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleDirichlet( std::string marker )
{
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );

    // <phat, mu>_Gamma_D
    M_a( 2_c, 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = idt( phat ) * id( l ) );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleNeumann( std::string marker )
{
    auto u = M_Vh->element( "u" );
    auto p = M_Wh->element( "p" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );

    // stabilisation parameter
    auto tau_constant = cst( doption( prefixvm( prefix(), "tau_constant" ) ) );

    // <j.n,mu>_Gamma_N
    M_a( 2_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = ( id( l ) * ( normalt( u ) ) ) );
    // <tau p, mu>_Gamma_N
    M_a( 2_c, 1_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = tau_constant * id( l ) * idt( p ) );
    // <-tau phat, mu>_Gamma_N
    M_a( 2_c, 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = -tau_constant * idt( phat ) * id( l ) );
}

template <int Dim, int Order, int G_Order, int E_Order>
template <typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRobin( Expr<ExprT> expr1, std::string marker )
{
    auto u = M_Vh->element( "u" );
    auto p = M_Wh->element( "p" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );
    // stabilisation parameter
    auto tau_constant = cst( doption( prefixvm( prefix(), "tau_constant" ) ) );

    // <j.n,mu>_Gamma_R
    M_a( 2_c, 0_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = ( id( l ) * ( normalt( u ) ) ) );
    // <tau p, mu>_Gamma_R
    M_a( 2_c, 1_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = tau_constant * id( l ) * idt( p ) );
    // <-tau phat, mu>_Gamma_R
    M_a( 2_c, 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = -tau_constant * idt( phat ) * id( l ) );
    // <g_R^1 phat, mu>_Gamma_R
    M_a( 2_c, 2_c ) += integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                  _expr = expr1 * idt( phat ) * id( l ) );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleIBC( int i, std::string markerOpt )
{
    auto u = M_Vh->element( "u" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto uI = M_Ch->element( "uI" );

    // stabilisation parameter
    auto tau_constant = cst( doption( prefixvm( this->prefix(), "tau_constant" ) ) );

    std::string marker;
    auto exAtMarker = M_IBCList[i];

    if ( !markerOpt.empty() )
    {
        marker = markerOpt;
    }
    else
    {
        marker = exAtMarker.marker();
    }

    // Feel::cout << "Matrix marker integral: " << marker << " with line " << i << std::endl;

    // <lambda, v.n>_Gamma_I
    M_a( 0_c, 3_c, 0, i ) += integrate( _range = markedfaces( M_mesh, marker ),
                                        _expr = idt( uI ) * ( normal( u ) ) );

    // <lambda, tau w>_Gamma_I
    M_a( 1_c, 3_c, 1, i ) += integrate( _range = markedfaces( M_mesh, marker ),
                                        _expr = tau_constant * idt( uI ) * id( w ) );

    // <j.n, m>_Gamma_I
    M_a( 3_c, 0_c, i, 0 ) += integrate( _range = markedfaces( M_mesh, marker ), _expr = ( normalt( u ) ) * id( nu ) );

    // <tau p, m>_Gamma_I
    M_a( 3_c, 1_c, i, 1 ) += integrate( _range = markedfaces( M_mesh, marker ),
                                        _expr = tau_constant * idt( p ) * id( nu ) );

    // -<lambda2, m>_Gamma_I
    M_a( 3_c, 3_c, i, i ) += integrate( _range = markedfaces( M_mesh, marker ),
                                        _expr = -tau_constant * id( nu ) * idt( uI ) );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::initTimeStep()
{
    // start or restart time step scheme
    if ( !this->doRestart() )
    {
        // start time step
        auto itField = modelProperties().boundaryConditions().find( "potential" );
        if ( itField != modelProperties().boundaryConditions().end() )
        {
            auto mapField = ( *itField ).second;
            auto itType = mapField.find( "InitialSolution" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : ( *itType ).second )
                {
                    if ( exAtMarker.isExpression() )
                    {
                        auto p_init = expr( exAtMarker.expression() );
                        auto marker = exAtMarker.marker();

                        if ( !this->isStationary() )
                            p_init.setParameterValues( {{"t", this->time()}} );
                        M_pp = project( _space = M_Wh, _range = markedelements( M_mesh, marker ), _expr = p_init );
                        if ( M_integralCondition )
                        {
                            auto mup = integrate( _range = markedfaces( M_mesh, M_IBCList[0].marker() ), _expr = idv( M_pp ) ).evaluate()( 0, 0 );
                            auto meas = integrate( _range = markedfaces( M_mesh, M_IBCList[0].marker() ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );

                            Feel::cout << "Initial integral value of potential on "
                                       << M_IBCList[0].marker() << " : \t " << mup / meas << std::endl;
                        }
                    }
                }
            }
        }

        M_bdf_mixedpoisson->start( M_pp );
        // up current time
        this->updateTime( M_bdf_mixedpoisson->time() );
    }
    else
    {
        // start time step
        M_bdf_mixedpoisson->restart();
        // load a previous solution as current solution
        M_pp = M_bdf_mixedpoisson->unknown( 0 );
        // up initial time
        this->setTimeInitial( M_bdf_mixedpoisson->timeInitial() );
        // restart exporter
        //this->restartPostProcess();
        // up current time
        this->updateTime( M_bdf_mixedpoisson->time() );

        this->log( "MixedPoisson", "initTimeStep", "restart bdf/exporter done" );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::updateTimeStepBDF()
{
    this->log( "MixedPoisson", "updateTimeStepBDF", "start" );
    this->timerTool( "TimeStepping" ).setAdditionalParameter( "time", this->currentTime() );
    this->timerTool( "TimeStepping" ).start();

    int previousTimeOrder = this->timeStepBDF()->timeOrder();

    M_bdf_mixedpoisson->next( M_pp );

    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf_mixedpoisson->time() );

    this->timerTool( "TimeStepping" ).stop( "updateBdf" );
    if ( this->scalabilitySave() ) this->timerTool( "TimeStepping" ).save();
    this->log( "MixedPoisson", "updateTimeStepBDF", "finish" );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::createTimeDiscretization()
{
    this->log( "MixedPoisson", "createTimeDiscretization", "start" );
    this->timerTool( "Constructor" ).start();

    std::string myFileFormat = soption( _name = "ts.file-format" ); // without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = ( boost::format( "_rank%1%_%2%" ) % this->worldComm().rank() % this->worldComm().size() ).str();
    M_bdf_mixedpoisson = bdf( _vm = Environment::vm(), _space = M_Wh,
                              _name = prefixvm( this->prefix(), prefixvm( this->subPrefix(), "p" + suffixName ) ),
                              _prefix = "",
                              _initial_time = this->timeInitial(),
                              _final_time = this->timeFinal(),
                              _time_step = this->timeStep(),
                              _restart = this->doRestart(),
                              _restart_path = this->restartPath(),
                              _restart_at_last_save = this->restartAtLastSave(),
                              _save = this->tsSaveInFile(), _freq = this->tsSaveFreq() );
    M_bdf_mixedpoisson->setfileFormat( myFileFormat );
    M_bdf_mixedpoisson->setPathSave( ( fs::path( this->rootRepository() ) /
                                       fs::path( prefixvm( this->prefix(), ( boost::format( "bdf_o_%1%_dt_%2%" ) % M_bdf_mixedpoisson->bdfOrder() % this->timeStep() ).str() ) ) )
                                         .string() );

    double tElapsed = this->timerTool( "Constructor" ).stop( "createTimeDiscr" );
    this->log( "MixedPoisson", "createTimeDiscretization", ( boost::format( "finish in %1% s" ) % tElapsed ).str() );
}

template <int Dim, int Order, int G_Order, int E_Order>
void MixedPoisson<Dim, Order, G_Order, E_Order>::exportResults( double time, mesh_ptrtype mesh, op_interp_ptrtype Idh, opv_interp_ptrtype Idhv )
{
    this->log( "MixedPoisson", "exportResults", "start" );
    this->timerTool( "PostProcessing" ).start();

    if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC && mesh )
    {
        LOG( INFO ) << "exporting on visualisation mesh at time " << time;
        M_exporter->step( time )->setMesh( mesh );
    }
    else if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
    {
        LOG( INFO ) << "exporting on computational mesh at time " << time;
        M_exporter->step( time )->setMesh( M_mesh );
    }

    // Export computed solutions
    {
        for ( auto const& field : modelProperties().postProcess().exports().fields() )  
        {
            if ( field == "flux" )
            {
                LOG(INFO) << "exporting flux at time " << time;
            
                M_exporter->step( time )->add(prefixvm(prefix(), "flux"), Idhv?(*Idhv)( M_up):M_up );
                if (M_integralCondition)
                {
                    double meas = 0.0;
                    double j_integral = 0;

                    for ( auto exAtMarker : this->M_IBCList )
                    {
                        auto marker = exAtMarker.marker();
                        LOG( INFO ) << "exporting integral flux at time "
                                    << time << " on marker " << marker;
                        j_integral = integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                                _expr = trans( idv( M_up ) ) * N() )
                                         .evaluate()( 0, 0 );
                        meas = integrate( _quad = _Q<expr_order>(), _range = markedfaces( M_mesh, marker ),
                                          _expr = cst( 1.0 ) )
                                   .evaluate()( 0, 0 );
                        Feel::cout << "Integral flux on " << marker << ": " << j_integral << std::endl;
                    }
                    M_exporter->step( time )->add( prefixvm( prefix(), "integralFlux" ), j_integral );
                    M_exporter->step( time )->add( prefixvm( prefix(), "integralVelocity" ), j_integral / meas );
                }
            }
            else if ( field == "scaled_flux" )
            {
                auto scaled_flux = M_Vh->element( "scaled_flux" );
                for ( auto const& pairMat : modelProperties().materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.getScalar( "scale_flux" );
                    scaled_flux.on( _range = markedelements( M_mesh, marker ), _expr = kk * idv( M_up ) );
                }

                M_exporter->step( time )->add( prefixvm( prefix(), "scaled_flux" ), Idhv ? ( *Idhv )( scaled_flux ) : scaled_flux );
            }
            else if ( field == "potential" )
            {
                LOG( INFO ) << "exporting potential at time " << time;
                M_exporter->step( time )->add( prefixvm( prefix(), "potential" ),
                                               Idh ? ( *Idh )( M_pp ) : M_pp );

                for ( int i = 0; i < M_integralCondition; i++ )
                {
                    double export_mup = M_mup[i].max();

                    LOG( INFO ) << "exporting IBC potential " << i << " at time "
                                << time << " value " << export_mup;

                    M_exporter->step( time )->add( prefixvm( prefix(), "cstPotential_1" ), export_mup );

                    Feel::cout << "Integral value of potential(mup) on "
                               << M_IBCList[i].marker() << " : \t " << export_mup << std::endl;
                }
                auto itField = modelProperties().boundaryConditions().find( "Exact solution" );
                if ( itField != modelProperties().boundaryConditions().end() )
                {
                    auto mapField = ( *itField ).second;
                    auto itType = mapField.find( "p_exact" );
                    if ( itType != mapField.end() )
                    {
                        for ( auto const& exAtMarker : ( *itType ).second )
                        {
                            if ( exAtMarker.isExpression() )
                            {
                                auto p_exact = expr( exAtMarker.expression() );
                                if ( !this->isStationary() )
                                    p_exact.setParameterValues( {{"t", time}} );
                                double K = 1;
                                for ( auto const& pairMat : modelProperties().materials() )
                                {
                                    auto material = pairMat.second;
                                    K = material.getDouble( "k" );
                                }
                                auto gradp_exact = grad<Dim>( p_exact );
                                if ( !this->isStationary() )
                                    gradp_exact.setParameterValues( { {"t", time } } );
                                auto u_exact = cst(-K)*trans(gradp_exact);//expr(-K* trans(gradp_exact)) ;

								auto p_exactExport = project( _space=M_Wh, _range=elements(M_mesh), _expr=p_exact );
								auto u_exactExport = project( _space=M_Vh, _range=elements(M_mesh), _expr=u_exact );

                                M_exporter->step( time )->add(prefixvm(prefix(), "p_exact"), p_exactExport );
								M_exporter->step( time )->add(prefixvm(prefix(), "u_exact"), u_exactExport );

                                // auto l2err_u = normL2( _range=elements(M_mesh), _expr= idv(M_up) - u_exact );
                                auto l2err_u = normL2( _range=elements(M_mesh), _expr= idv(M_up) - idv(u_exactExport) );
								auto l2norm_uex = normL2( _range=elements(M_mesh), _expr= u_exact );

                                if (l2norm_uex < 1)
                                    l2norm_uex = 1.0;

                                auto l2err_p = normL2( _range = elements( M_mesh ), _expr = p_exact - idv( M_pp ) );
                                auto l2norm_pex = normL2( _range = elements( M_mesh ), _expr = p_exact );
                                if ( l2norm_pex < 1 )
                                    l2norm_pex = 1.0;

                                Feel::cout << "----- Computed Errors -----" << std::endl;
                                Feel::cout << "||p-p_ex||_L2=\t" << l2err_p / l2norm_pex << std::endl;
                                Feel::cout << "||u-u_ex||_L2=\t" << l2err_u / l2norm_uex << std::endl;
                                Feel::cout << "---------------------------" << std::endl;

                                // Export the errors
                                M_exporter->step( time )->add( prefixvm( prefix(), "p_error_L2" ), l2err_p / l2norm_pex );
                                M_exporter->step( time )->add( prefixvm( prefix(), "u_error_L2" ), l2err_u / l2norm_uex );
                            }
                        }
                    }
                }
            }
            else if ( field == "scaled_potential" )
            {
                auto scaled_potential = M_Wh->element( "scaled_potential" );
                for ( auto const& pairMat : modelProperties().materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.getScalar( "scale_potential" );

                    scaled_potential.on( _range = markedelements( M_mesh, marker ), _expr = kk * idv( M_pp ) );
                }

                M_exporter->step( time )->add( prefixvm( prefix(), "scaled_potential" ), Idh ? ( *Idh )( scaled_potential ) : scaled_potential );

                for ( int i = 0; i < M_integralCondition; i++ )
                {
                    auto scaled_ibc = M_mup[i].max();
                    for( auto const& pairMat : modelProperties().materials() )
                    {
                        auto material = pairMat.second;
                        auto kk_ibc = material.getScalar( "scale_potential" ).evaluate();
                        scaled_ibc = scaled_ibc * kk_ibc;
                    }                    
                    
                    LOG(INFO) << "exporting IBC scaled potential " << i << " at time "
                              << time << " value " << scaled_ibc;
                    M_exporter->step( time )->add(prefixvm(prefix(), "scaled_cstPotential_1"),
                                                  scaled_ibc );
                    Feel::cout << "Integral scaled value of potential(mup) on "
                               << M_IBCList[i].marker() << " : \t " << scaled_ibc << std::endl;
                }
            }
            else if ( field != "state variable" )
            {
                // Import data
                LOG( INFO ) << "importing " << field << " at time " << time;
                double extra_export = 0.0;
                auto itField = modelProperties().boundaryConditions().find( "Other quantities" );
                if ( itField != modelProperties().boundaryConditions().end() )
                {
                    auto mapField = ( *itField ).second;
                    auto itType = mapField.find( field );
                    if ( itType != mapField.end() )
                    {
                        for ( auto const& exAtMarker : ( *itType ).second )
                        {
                            if ( exAtMarker.isExpression() )
                            {
                                LOG( INFO ) << "WARNING: you are trying to export a single expression";
                            }
                            else if ( exAtMarker.isFile() )
                            {
                                if ( !this->isStationary() )
                                {
                                    extra_export = exAtMarker.data( M_bdf_mixedpoisson->time() );
                                }
                                else
                                    extra_export = exAtMarker.data( 0.1 );
                            }
                        }
                    }
                }
                // Transform data if necessary
                LOG( INFO ) << "transforming " << field << "at time " << time;
                std::string field_k = field;
                field_k += "_k";
                double kk = 0.0;
                for ( auto const& pairMat : modelProperties().materials() )
                {
                    auto material = pairMat.second;
                    kk = material.getDouble( field_k );
                }
                if ( std::abs( kk ) > 1e-10 )
                    extra_export *= kk;

                // Export data
                LOG( INFO ) << "exporting " << field << " at time " << time;
                M_exporter->step( time )->add( prefixvm( prefix(), field ), extra_export );
            }
        }
    }

    this->timerTool( "PostProcessing" ).stop( "exportResults" );
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool( "PostProcessing" ).setAdditionalParameter( "time", this->currentTime() );
        this->timerTool( "PostProcessing" ).save();
    }
    this->log( "MixedPoisson", "exportResults", "finish" );
}

} // Namespace FeelModels

} // Namespace Feel

#endif
