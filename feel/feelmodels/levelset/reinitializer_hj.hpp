#ifndef _REINITIALIZER_HJ_HPP
#define _REINITIALIZER_HJ_HPP 1

namespace Feel
{

template<typename FunctionSpaceType>
class ReinitializerHJ
: public Reinitializer<FunctionSpaceType>
{
public:
    //--------------------------------------------------------------------//
    // Typedefs
    typedef Reinitializer<FunctionSpaceType> super_type;
    typedef ReinitializerHJ<FunctionSpaceType> self_type;
    typedef boost:shared_ptr<self_type> self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    //--------------------------------------------------------------------//
    // Hamilton-Jacobi advection
    template<typename SpaceType>
    class AdvectionHJ
        : public AdvectionBase<
            typename SpaceType::convex_type, 
            typename SpaceType::basis_type, 
            typename SpaceType::periodicity_0_type
          >
        , public boost::enable_shared_from_this< AdvectionHJ<SpaceType> >
    {
    public:
        typedef AdvectionBase<
            typename SpaceType::convex_type, 
            typename SpaceType::basis_type, 
            typename SpaceType::periodicity_0_type
          > super_type;

        typedef AdvectionHJ<SpaceType> self_type;
        typedef boost::shared_ptr<self_type> self_ptrtype;

        typedef SpaceType functionspace_type;
        typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

        // Constructor
        AdvectionHJ(
                functionspace_ptrtype const& space,
                std::string const& prefix = ""
                );
        
        static self_ptrtype New(
                functionspace_ptrtype const& space,
                std::string const& prefix = ""
                );

        // Initialization
        void init( bool buildModelAlgebraicFactory = true );

        // BC and source term assembly
        void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        bool hasSourceTerm() const { return false; }
    };

    typedef AdvectionHJ<FunctionSpaceType> advectionhj_type;
    typedef boost::shared_ptr<advectionhj_type> advectionhj_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    ReinitializerHJ( functionspace_ptrtype const& space );
    //--------------------------------------------------------------------//
    // Parameters
    double tolerance() const { return M_tolerance; }
    void setTolerance( double tol ) { M_tolerance = tol; }
    
    double timeStep() const { return M_advectionHJ->timeStep(); }
    void setTimeStep( double dt ) { M_advectionHJ->setTimeStep( dt ); }

    int maxIterations() const { return M_maxIterations; }
    void setMaxIterations( int max ) { M_maxIterations = max; }

    double heavisideThickness const { return M_heavisideThickness; }
    void setHeavisideThickness ( double eps ) { M_heavisideThickness = eps; }
    //--------------------------------------------------------------------//
    // Run reinitialization
    element_type run( element_type const& phi );

private:
    advectionhj_ptrtype M_advectionHJ;

    double M_tolerance;
    int M_maxIterations;
    double M_heavisideThickness;
};

#define REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS \
   template<typename FunctionSpaceType>
#define REINITIALIZERHJ_CLASS_TEMPLATE_TYPE \
    ReinitializerHJ<FunctionSpaceType>

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ReinitializerHJ( functionspace_ptrtype const& space )
    : super_type( space )
{
    M_advectionHJ = advectionhj_type::New( space );
    M_advectionHJ->init();
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
typename REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::element_type
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::run( element_type const& phi )
{
    auto mesh = M_advectionHJ->mesh();
    auto space = M_advectionHJ->functionSpace();
    // Set initial value
    M_advectionHJ->fieldSolution() = phi;
    M_advectionHJ->init( false );
    // Init convergence test utilities
    double err_dist;
    double err_disto = integrate(
            elements(this->mesh()),
            vf::abs( sqrt(gradv(phi)*trans(gradv(phi))) - 1.0 )
            ).evaluate()(0,0);

    double rateChangePhiL2;
    double rateChangePhiL2o = 1.;
    double relativeRateChangePhiL2;

    for(int i = 0; i < this->maxIterations(); ++i)
    {
        LOG(INFO) << "iter reinit : " << i << std::endl;

        auto phi_reinit = M_advectionHJ->fieldSolutionPtr();
        auto phi_reinito = M_advectionHJ->timeStepBDF()->unknowns()[1];

        auto phi_sign = idv(phi_reinit) / vf::max( vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit))), 0.92 );

        auto H_reinit = vf::project(
                space, 
                elements(mesh),
                vf::abs(
                    ( phi_sign < -M_heavisideThickness )*vf::constant(0.0)
                    +
                    ( phi_sign >= -M_heavisideThickness && phi_sign <= M_heavisideThickness )*
                    0.5*(1 + phi_sign/M_heavisideThickness + 1/M_PI*vf::sin( M_PI*phi_sign/M_heavisideThickness ) )
                    +
                    ( phi_sign > M_heavisideThickness )*vf::constant(1.0) )
                );
        auto Sign = 2*( idv(H_reinit)-0.5 );

        // Update advection source
        M_advectionHJ->updateSourceAdded( Sign );
        // Update advection velocity
        auto beta = Sign * trans(gradv(phi_reinit)) / vf::max( vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit))), 0.92 );
        M_advection->updateAdvectionVelocity( beta );

        // Solve
        M_advectionHJ->solve();

        // Test convergence
        err_dist = integrate(
                elements(this->mesh()),
                vf::abs( sqrt(gradv(phi)*trans(gradv(phi_reinit))) - 1.0 )
                ).evaluate()(0,0);
        double delta_err = std::abs(err_disto - err_dist) / err_disto;
        LOG(INFO)<<"delta err = "<<delta_err<<"\n";

        rateChangePhiL2 = std::sqrt( integrate(
                elements(mesh),
                (idv(phi_reinit) - idv(phi_reinito))*(idv(phi_reinit) - idv(phi_reinito))
                ).evaluate()(0,0) );
        relativeRateChangePhiL2 = std::abs( rateChangePhiL2 - rateChangePhiL2o) / rateChangePhiL2o;

        if ( (i!=0) && (relativeRateChangePhiL2 <= M_tolerance ) )
        {
            LOG(INFO) << "stop Hamilton Jacobi iterations because tolerence threshold has been reached\n";
            LOG(INFO) << "relative rate of change = " << relativeRateChangePhiL2 << std::endl;
            break;
        }
        
        // Update "old" variables
        err_disto = err_dist;
        rateChangePhiL2o = rateChangePhiL2;
    }

    Feel::cout << "reinit done in " << __iter - start_iter << " iter\n";
    Feel::cout << "final relative rate of change = " << relativeRateChangePhiL2 << std::endl;

    return M_advectionHJ->fieldSolution();
}

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
// Hamilton-Jacobi advection
#define ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS \
    template<typename ConvexType, typename BasisAdvectionType, typename PeriodicityType>
#define ADVECTIONHJ_CLASS_TEMPLATE_TYPE \
    AdvectionHJ<ConvexType, BasisAdvectionType, PeriodicityType>

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::AdvectionHJ(
        functionspace_ptrtype const& space,
        std::string const& prefix )
    : super_type( prefix )
{
    super_type::build( space );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::self_ptrtype
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::New(
        functionspace_ptrtype const& space,
        std::string const& prefix )
{
    return boost::make_shared<self_type>( space, prefix );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::init(
        bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F,
        bool buildCstPart) const
{
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F) const
{}

#undef ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
#undef ADVECTIONHJ_CLASS_TEMPLATE_TYPE

#undef REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
#undef REINITIALIZERHJ_CLASS_TEMPLATE_TYPE

} // namespace Feel

#endif
