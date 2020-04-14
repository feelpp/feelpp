#ifndef _MIXEDPOISSON2_HPP
#define _MIXEDPOISSON2_HPP

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <boost/hana/tuple.hpp>

#include <feel/feelmodels/hdg/enums.hpp>

#define USE_SAME_MAT 1

namespace Feel {

namespace FeelModels {

inline po::options_description
makeMixedPoissonOptions( std::string const&  _prefix = "", std::string const&  _toolbox_prefix = "hdg.poisson" )
{
    std::string prefix = _toolbox_prefix.empty()?"hdg.poisson":_toolbox_prefix;
    if ( !_prefix.empty() )
        prefix = prefixvm( prefix, _prefix );
    po::options_description mpOptions( "Mixed Poisson HDG options" );
    mpOptions.add_options()( "gmsh.submesh", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "tau_constant" ).c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix, "tau_order" ).c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges" ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "hface" ).c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( prefixvm( prefix, "conductivity_json" ).c_str(), po::value<std::string>()->default_value( "cond" ), "key for conductivity in json" )
        ( prefixvm( prefix, "conductivityNL_json" ).c_str(), po::value<std::string>()->default_value( "condNL" ), "key for non linear conductivity in json (depends on potential p)" )
        ( prefixvm( prefix, "use-sc" ).c_str(), po::value<bool>()->default_value( true ), "use static condensation" )
        ( prefixvm( prefix, "error-quadrature").c_str(), po::value<int>()->default_value(10), "quadrature to compute errors" )
        ( prefixvm( prefix, "set-zero-by-init").c_str(), po::value<bool>()->default_value(true), "reinit matrix and vector when setting to zero" )
        ;
    mpOptions.add( modelnumerical_options( prefix ) );
    mpOptions.add( backend_options( prefix + ".sc" ) );
    return mpOptions;
}

inline po::options_description
makeMixedPoissonLibOptions( std::string const&  prefix = "", std::string const&  _toolbox_prefix = "hdg.poisson" )
{
    po::options_description mpLibOptions( "Mixed Poisson HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order, int G_Order = 1, int E_Order = 4>
class MixedPoisson    : public ModelNumerical
{
public:
    typedef ModelNumerical super_type;

    static const uint16_type expr_order = (Order+E_Order)*G_Order;
    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedPoisson<Dim,Order,G_Order,E_Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,G_Order,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;

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
    // Whp
    using Whp_t = Pdh_type<mesh_type,Order+1>;
    using Whp_ptr_t = Pdh_ptrtype<mesh_type,Order+1>;
    using Whp_element_t = typename Whp_t::element_type;
    using Whp_element_ptr_t = typename Whp_t::element_ptrtype;
    // Mh
    using Mh_t = Pdh_type<face_mesh_type,Order>;
    using Mh_ptr_t = Pdh_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
    // Ch
    using Ch_t = Pch_type<face_mesh_type,0>;
    using Ch_ptr_t = Pch_ptrtype<face_mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    using Ch_element_vector_type = std::vector<Ch_element_t>;
    // M0h
    using M0h_t = Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t = Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    using product2_space_type = ProductSpaces2<Ch_ptr_t,Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using product2_space_ptrtype = std::shared_ptr<product2_space_type>;

    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef std::shared_ptr <exporter_type> exporter_ptrtype;

    using op_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Wh_t, Pdh_type<mesh_type,Order>>>;
    using opv_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Vh_t, Pdhv_type<mesh_type,Order>>>;

    using integral_boundary_list_type = std::vector<ExpressionStringAtMarker>;

    // time
    // typedef Lagrange<Order,Scalar,Discontinuous> basis_scalar_type;
    // typedef FunctionSpace<mesh_type, bases<basis_scalar_type>> space_mixedpoisson_type;

    // typedef Bdf<space_mixedpoisson_type>  bdf_type;
    typedef Bdf <Wh_t> bdf_type;
    // typedef Bdf<Pdh_type<mesh_type,Order>> bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;

    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<Wh_t,Vh_t> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

    struct FieldTag
    {
        static auto potential( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
        static auto flux( self_type const* t ) { return ModelFieldTag<self_type,1>( t ); }
    };
    //private:
protected:
    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Whp_ptr_t M_Whp; // postprocess potential
    Mh_ptr_t M_Mh; // potential trace
    Ch_ptr_t M_Ch; // Lagrange multiplier
    M0h_ptr_t M_M0h;
    product2_space_ptrtype M_ps;

    backend_ptrtype M_backend;
    condensed_matrix_ptr_t<value_type> M_A_cst;
#ifndef USE_SAME_MAT
    condensed_matrix_ptr_t<value_type> M_A;
#endif
    condensed_vector_ptr_t<value_type> M_F;
    condensed_matrix_ptr_t<value_type> M_App;
    condensed_vector_ptr_t<value_type> M_Fpp;
    vector_ptrtype M_U;

    Vh_element_t M_up; // flux solution
    Wh_element_t M_pp; // potential solution
    Whp_element_t M_ppp; // postprocess potential solution
    Ch_element_vector_type M_mup; // potential solution on the integral boundary conditions

    // time discretization
    bdf_ptrtype M_bdf_mixedpoisson;

    MixedPoissonPhysics M_physic;
    std::string M_potentialKey;
    std::string M_fluxKey;

    int M_tauOrder;
    double M_tauCst;
    int M_hFace;
    std::string M_conductivityKey;
    std::string M_nlConductivityKey;
    bool M_useSC;

    int M_integralCondition;
    int M_useUserIBC;
    integral_boundary_list_type M_IBCList;

    bool M_isPicard;
    std::map<std::string,value_type> M_paramValues;

    measure_points_evaluation_ptrtype M_measurePointsEvaluation;

    int M_quadError;
    bool M_setZeroByInit;
public:

    // constructor
    MixedPoisson( std::string const& prefix = "hdg.poisson",
                  MixedPoissonPhysics const& physic = MixedPoissonPhysics::None,
                  worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                  std::string const& subPrefix = "",
                  ModelBaseRepository const& modelRep = ModelBaseRepository() );

    MixedPoisson( self_type const& MP ) = default;
    static self_ptrtype New( std::string const& prefix = "hdg.poisson",
                             MixedPoissonPhysics const& physic = MixedPoissonPhysics::None,
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    // Get Methods
    std::string potentialKey() const { return M_potentialKey; }
    std::string fluxKey() const { return M_fluxKey; }
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Whp_ptr_t postPotentialSpace() const { return M_Whp; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    M0h_ptr_t traceSpaceOrder0() const { return M_M0h; }
    Ch_ptr_t constantSpace() const {return M_Ch;}

    Vh_element_t const& fluxField() const { return M_up; }
    Wh_element_t const& potentialField() const { return M_pp; }
    Whp_element_t const& postPotentialField() const { return M_ppp; }
    Vh_element_t & fluxField() { return M_up; }
    Wh_element_t & potentialField() { return M_pp; }
    Whp_element_t & postPotentialField() { return M_ppp; }
    integral_boundary_list_type integralBoundaryList() const { return M_IBCList; }
    int integralCondition() const { return M_integralCondition; }
    void setIBCList(std::vector<std::string> markersIbc);
    backend_ptrtype get_backend() { return M_backend; }
    product2_space_ptrtype getPS() const { return M_ps; }
    condensed_vector_ptr_t<value_type> getF() { return M_F; }

    int tauOrder() const { return M_tauOrder; }
    void setTauOrder(int order) { M_tauOrder = order; }
    double tauCst() const { return M_tauCst; }
    void setTauCst(double cst) { M_tauCst = cst; }
    int hFace() const { return M_hFace; }
    void setHFace(int h) { M_hFace = h; }
    std::string conductivityKey() const { return M_conductivityKey; }
    void setConductivityKey(std::string key) { M_conductivityKey = key; }
    std::string nlConductivityKey() const { return M_nlConductivityKey; }
    void setNlConductivityKey(std::string key) { M_nlConductivityKey = key; }
    bool useSC() const { return M_useSC; }
    void setUseSC(bool sc) { M_useSC = sc; }

    // time step scheme
    virtual void createTimeDiscretization() ;
    bdf_ptrtype timeStepBDF() { return M_bdf_mixedpoisson; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf_mixedpoisson; }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepBDF(); }

    // Exporter
    virtual void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr )
        {
            this->exportResults (this->currentTime(), mesh );
            M_exporter -> save();
        }
    void exportResults ( double Time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  ) ;
    exporter_ptrtype M_exporter;
    exporter_ptrtype exporterMP() { return M_exporter; }

    void init( mesh_ptrtype mesh = nullptr, mesh_ptrtype meshVisu = nullptr);

    virtual void initModel();
    virtual void initSpaces();
    virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void initMatricesAndVector();
    virtual void assembleAll();
    virtual void assembleCstPart();
    virtual void assembleNonCstPart();
    void copyCstPart();
    void setMatricesAndVectorToZero();

    void assembleRHS();
    template<typename ExprT> void updateConductivityTerm( Expr<ExprT> expr, std::string marker = "");
    void updateConductivityTerm( bool isNL = false);

    template<typename ExprT> void assembleFluxRHS( Expr<ExprT> expr, std::string marker);
    template<typename ExprT> void assemblePotentialRHS( Expr<ExprT> expr, std::string marker);

    void assembleBoundaryCond();
    void assembleRhsBoundaryCond();
    void assembleDirichlet( std::string marker);
    void assembleNeumann( std::string marker);
    template<typename ExprT> void assembleRhsDirichlet( Expr<ExprT> expr, std::string marker);
    template<typename ExprT> void assembleRhsNeumann( Expr<ExprT> expr, std::string marker);
    template<typename ExprT> void assembleRhsInterfaceCondition( Expr<ExprT> expr, std::string marker);
    // u.n + g1.p = g2
    template<typename ExprT1, typename ExprT2> void assembleRobin( Expr<ExprT1> const& expr1, Expr<ExprT2> const& expr2, std::string const& marker, bool update_only = false );

    void assembleIBC(int i, std::string marker = "");
    virtual void assembleRhsIBC(int i, std::string marker = "", double intjn = 0);

    virtual void solve();

    void assemblePostProcessCstPart();
    void assemblePostProcessNonCstPart( bool isNL = false );
    template<typename ExprT> void assemblePostProcessRhs( Expr<ExprT> expr, std::string marker = "");
    void solvePostProcess();
    virtual void postProcess( bool isNL = false );

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields(
                modelField<FieldCtx::ID|FieldCtx::GRAD|FieldCtx::GRAD_NORMAL>( FieldTag::potential(this), prefix, MixedPoissonPhysicsMap[M_physic]["potentialK"], this->potentialField(), MixedPoissonPhysicsMap[M_physic]["potentialSymbol"], this->keyword() ),
                modelField<FieldCtx::ID>( FieldTag::flux(this), prefix, MixedPoissonPhysicsMap[M_physic]["fluxK"], this->fluxField(), MixedPoissonPhysicsMap[M_physic]["fluxSymbol"], this->keyword() )
                                                 );
            // return this->modelFields( /*this->potentialField(), */prefix );
        }
    // auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
    //     {
    //         auto field_p = this->spaceElectricPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "potential-electric" ) );
    //         return this->modelFields( field_p, prefix );
    //     }
    // template <typename PotentialFieldType>
    // auto modelFields( /*PotentialFieldType const& field_p, */std::string const& prefix = "" ) const
    //     {
    //     }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            // auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            // auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols electric_P, electric_grad_P(_x,_y,_z), electric_dn_P
            //Feel::cout << seFields.names() << std::endl;
            //std::cout << "Info field potential = " << mfields.field( FieldTag::potential(this), MixedPoissonPhysicsMap[M_physic]["potentialK"] ).functionSpace()->nDof() << std::endl;
            //std::cout << "Info field flux = " << mfields.field( FieldTag::flux(this), MixedPoissonPhysicsMap[M_physic]["fluxK"] ).functionSpace()->nDof() << std::endl;
            return Feel::vf::symbolsExpr( /*seToolbox,*/ seParam, /*seMat,*/ seFields );
        }
#if 0 // NOT USE this one because field not use shared ptr : (object modelfields is temporary here and the fields are stored inside)
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }
#endif
#if 0
    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            auto const& v = mfields.field( FieldTag::potential(this), MixedPoissonPhysicsMap[M_physic]["potentialK"] );

            // generate symbol electric_matName_current_density
            typedef decltype( this->currentDensityExpr(v,"") ) _expr_currentdensity_type;
            std::vector<std::tuple<std::string,_expr_currentdensity_type,SymbolExprComponentSuffix>> currentDensitySymbs;
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physic() ) )
            {
                std::string symbolcurrentDensityStr = prefixvm( this->keyword(), (boost::format("%1%_current_density") %matName).str(), "_");
                auto _currentDensityExpr = this->currentDensityExpr( v, matName );
                currentDensitySymbs.push_back( std::make_tuple( symbolcurrentDensityStr, _currentDensityExpr, SymbolExprComponentSuffix( nDim,1 ) ) );
            }

            return Feel::vf::symbolsExpr( symbolExpr( currentDensitySymbs ) );
        }
#endif
    template <typename ModelFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr )
        {
            bool hasMeasure = false;
            bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
            bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
            bool hasMeasurePoint = this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, mfields );
            if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
                hasMeasure = true;

            if ( hasMeasure )
            {
                if ( !this->isStationary() )
                    this->postProcessMeasuresIO().setMeasure( "time", time );
                this->postProcessMeasuresIO().exportMeasures();
                this->upload( this->postProcessMeasuresIO().pathFile() );
            }
        }
};

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::updateConductivityTerm(Expr<ExprT> expr, std::string marker)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );

#ifdef USE_SAME_MAT
    auto bbf = blockform2( *M_ps, M_A_cst);
#else
    auto bbf = blockform2( *M_ps, M_A);
#endif

    if ( marker.empty() )
        bbf(0_c,0_c) += integrate( _range=elements(support(M_Wh)), _expr=inner(idt(u),id(v))/expr);
    else
        bbf(0_c,0_c) += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/expr);

    // (1/delta_t p, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() ) {
        bbf( 1_c, 1_c ) += integrate(_range=elements(support(M_Wh)),
                                     _expr = -(this->timeStepBDF()->polyDerivCoefficient(0)*idt(p)*id(w)) );
    }

}

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleFluxRHS( Expr<ExprT> expr, std::string marker)
{
    tic();
    auto blf = blockform1( *M_ps, M_F );
    auto v = M_Vh->element();

    if ( marker.empty() )
        blf(0_c) += integrate( _range=elements(support(M_Wh)),
                              _expr=inner(expr,id(v)) );
    else
        blf(0_c) += integrate( _range=markedelements(M_mesh,marker),
                              _expr=inner(expr,id(v)) );
    toc("assembleFluxRhs", this->verbose() || FLAGS_v > 0);
}

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void MixedPoisson<Dim, Order, G_Order, E_Order>::assemblePotentialRHS( Expr<ExprT> expr, std::string marker)
{
    tic();
    auto blf = blockform1( *M_ps, M_F );
    auto w = M_Wh->element();

    if ( marker.empty() )
    {
        blf(1_c) += integrate( _range=elements(support(M_Wh)),
                              _expr=inner(expr,id(w)) );
    }
    else
    {
        blf(1_c) += integrate( _range=markedelements(M_mesh,marker),
                              _expr=inner(expr,id(w)) );
    }
    toc("assemblePotentialRhs", this->verbose() || FLAGS_v > 0);
}

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsDirichlet( Expr<ExprT> expr, std::string marker)
{
    tic();
    auto blf = blockform1( *M_ps, M_F );
    auto l = M_Mh->element( "lambda" );

    // <g_D, mu>_Gamma_D
    blf(2_c) += integrate(_range=markedfaces(support(M_Wh),marker),
                          _expr=id(l)*expr);
    toc("assembleRhsDirichlet", this->verbose() || FLAGS_v > 0);
}

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsNeumann( Expr<ExprT> expr, std::string marker)
{
    tic();
    auto blf = blockform1( *M_ps, M_F );
    auto l = M_Mh->element( "lambda" );

    // <g_N,mu>_Gamma_N
    blf(2_c) += integrate( _range=markedfaces(support(M_Wh), marker),
                          _expr=id(l)*expr);
    toc("assembleRhsNeumann", this->verbose() || FLAGS_v > 0);
}

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsInterfaceCondition( Expr<ExprT> expr, std::string marker)
{
    tic();
    auto blf = blockform1( *M_ps, M_F );
    auto l = M_Mh->element( "lambda" );

    // <g_interface,mu>_Gamma_N
    blf(2_c) += integrate( _range=markedelements(M_Mh->mesh(), marker),
                          _expr=id(l)*expr);
    toc("assembleRhsInterface", this->verbose() || FLAGS_v > 0);
}

/*
 template<int Dim, int Order, int G_Order, int E_Order>
 template<typename ExprT>
 void MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhsRobin( Expr<ExprT> expr1, Expr<ExprT> expr2, std::string marker)
 {
 #ifdef USE_SAME_MAT
 auto bbf = blockform2( *M_ps, M_A_cst);
 #else
 auto bbf = blockform2( *M_ps, M_A);
 #endif
 auto blf = blockform1( *M_ps, M_F );
 auto u = M_Vh->element( "u" );
 auto p = M_Wh->element( "p" );
 auto phat = M_Mh->element( "phat" );
 auto l = M_Mh->element( "lambda" );
 auto H = M_M0h->element( "H" );
 if ( ioption(prefixvm(prefix(), "hface") ) == 0 )
 H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
 else if ( ioption(prefixvm(prefix(), "hface") ) == 1 )
 H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
 else if ( ioption(prefixvm(prefix(), "hface") ) == 2 )
 H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
 else
 H.on( _range=elements(M_M0h->mesh()), _expr=h() );
 // stabilisation parameter
 auto tau_constant = cst(doption(prefixvm(prefix(), "tau_constant")));

 // <j.n,mu>_Gamma_R
 bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
 _expr=( id(l)*(trans(idt(u))*N()) ));
 // <tau p, mu>_Gamma_R
 bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
 _expr=tau_constant * id(l) * ( pow(idv(H),M_tau_order)*idt(p) ) );
 // <-tau phat, mu>_Gamma_R
 bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
 _expr=-tau_constant * idt(phat) * id(l) * ( pow(idv(H),M_tau_order) ) );
 // <g_R^1 phat, mu>_Gamma_R
 bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
 _expr=expr1*idt(phat) * id(l) );
 // <g_R^2,mu>_Gamma_R
 blf(2_c) += integrate( _range=markedfaces(support(M_Wh), marker),
 _expr=id(l)*expr2);
 }
 */

template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT1, typename ExprT2>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::assembleRobin( Expr<ExprT1> const& expr1, Expr<ExprT2> const& expr2, std::string const& marker, bool update_only )
{
    tic();
    auto bbf = blockform2( *M_ps, M_A_cst);

    auto blf = blockform1( *M_ps, M_F );
    auto u = M_Vh->element( "u" );
    auto p = M_Wh->element( "p" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );
    // auto H = M_M0h->element( "H" );
    // if ( ioption(prefixvm(prefix(), "hface") ) == 0 )
    //     H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    // else if ( ioption(prefixvm(prefix(), "hface") ) == 1 )
    //     H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    // else if ( ioption(prefixvm(prefix(), "hface") ) == 2 )
    //     H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    // else
    //     H.on( _range=elements(M_M0h->mesh()), _expr=h() );
    // stabilisation parameter
    auto tau_constant = cst(M_tauCst);

    if ( !update_only )
    {
        // <j.n,mu>_Gamma_R
        bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
                                     _expr=id(l)*normalt(u) );
        // <tau p, mu>_Gamma_R
        bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
                                     _expr=tau_constant * id(l) * idt(p)  );
        // <-tau phat, mu>_Gamma_R
        bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
                                     _expr=-tau_constant * idt(phat) * id(l) );
    }
    // <g_R^1 phat, mu>_Gamma_R
    bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh),marker),
                                 _expr=expr1*idt(phat) * id(l) );
    // <g_R^2,mu>_Gamma_R
    blf(2_c) += integrate( _range=markedfaces(support(M_Wh), marker),
                           _expr=id(l)*expr2);
    toc("assembleRobin", this->verbose() || FLAGS_v > 0);
}


template<int Dim, int Order, int G_Order, int E_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order, E_Order>::assemblePostProcessRhs(Expr<ExprT> expr, std::string marker)
{
    auto pps = product( M_Whp );
    auto ell = blockform1( pps, M_Fpp);
    ell(0_c) += integrate( _range=markedelements(support(M_Wh),marker),
                          _expr=-grad(M_ppp)*idv(M_up)/expr);
}

} // Namespace FeelModels

} // Namespace Feel

#endif
