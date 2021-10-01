#ifndef _HELFRICHFORCEMODEL_HPP
#define _HELFRICHFORCEMODEL_HPP 1

//#define DEBUG_HELFRICHFORCEMODEL

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>
#include <feel/feelmodels/multifluid/helfrichforceexpr.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType, class FluidMechanicsType>
class HelfrichForceModel
: public virtual InterfaceForcesModel<LevelSetType, FluidMechanicsType>
{
    typedef HelfrichForceModel<LevelSetType, FluidMechanicsType> self_type;
    typedef InterfaceForcesModel<LevelSetType, FluidMechanicsType> super_type;
public:
    typedef typename super_type::levelset_type levelset_type;
    typedef typename super_type::levelset_ptrtype levelset_ptrtype;

    typedef typename super_type::fluidmechanics_type fluidmechanics_type;
    typedef typename super_type::fluidmechanics_ptrtype fluidmechanics_ptrtype;

    typedef typename super_type::mesh_type mesh_type;

    typedef typename levelset_type::space_vectorial_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    typedef typename levelset_type::element_levelset_type element_scalar_type;
    typedef typename levelset_type::element_levelset_ptrtype element_scalar_ptrtype;

    typedef typename fluidmechanics_type::DataUpdateJacobian DataUpdateJacobian;
    typedef typename fluidmechanics_type::DataUpdateResidual DataUpdateResidual;
    typedef typename fluidmechanics_type::DataUpdateLinear DataUpdateLinear;

    static const uint16_type Dim = levelset_type::nDim;

    enum class HelfrichMethod {
        NODAL_PROJECTION, L2_PROJECTION, SMOOTH_PROJECTION, DIFFUSION
    };

    //--------------------------------------------------------------------//
    // Construction
    HelfrichForceModel() = default;
    HelfrichForceModel( HelfrichForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm = fluidmechanics_ptrtype() ) override;

    std::shared_ptr<std::ostringstream> getInfo() const override;

    void loadParametersFromOptionsVm();
    //--------------------------------------------------------------------//
    double bendingModulus() const { return M_helfrichBendingModulus; }
    void setBendingModulus( double k ) { M_helfrichBendingModulus = k; }

    //--------------------------------------------------------------------//
    virtual void updateFluidInterfaceForcesLinearPDE( DataUpdateLinear & data ) const override;
    virtual void updateFluidInterfaceForcesJacobian( DataUpdateJacobian & data ) const override;
    virtual void updateFluidInterfaceForcesResidual( DataUpdateResidual & data ) const override;
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F ) const override;
    void addHelfrichForce( element_ptrtype & F, int impl ) const;

    element_scalar_type curvatureDiffusionWillmore( mpl::int_<1> /*Dim*/ ) const {}
    element_scalar_type curvatureDiffusionWillmore( mpl::int_<2> /*Dim*/ ) const;
    element_scalar_type curvatureDiffusionWillmore( mpl::int_<3> /*Dim*/ ) const;

    //--------------------------------------------------------------------//
    double M_helfrichBendingModulus;
    int M_forceImpl;
    HelfrichInnerDivImplementation M_helfrichInnerDivExprImpl;
    bool M_useIntegrationByParts;
    HelfrichMethod M_helfrichForceMethod;

#ifdef DEBUG_HELFRICHFORCEMODEL
    typedef std::shared_ptr<Exporter<mesh_type, 1>> exporter_ptrtype;
    exporter_ptrtype M_exporter;
    bool M_exporterInitDone;
#endif
};

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm )
{
    super_type::build( prefix, ls, fm );
    this->loadParametersFromOptionsVm();
    // Ensures levelset curvature diffusion is built if needed
    if( M_helfrichForceMethod == HelfrichMethod::DIFFUSION )
    {
        this->levelset()->toolManager()->initCurvatureDiffusion();
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
std::shared_ptr<std::ostringstream> 
HelfrichForceModel<LevelSetType, FluidMechanicsType>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "Helfrich force ("
           << "bending modulus = " << this->M_helfrichBendingModulus
           << ")";

    return _ostr;
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::loadParametersFromOptionsVm()
{
    M_helfrichBendingModulus = doption( _name="helfrich-bending-modulus", _prefix=this->prefix() ); 
    M_forceImpl = ioption( _name="helfrich-force-impl", _prefix=this->prefix() );
    std::string helfrichInnerDivExprImpl = soption( _name="helfrich-inner-div-impl", _prefix=this->prefix() );
    if( helfrichInnerDivExprImpl == "generic" )
        M_helfrichInnerDivExprImpl = HelfrichInnerDivImplementation::GENERIC;
    else if( helfrichInnerDivExprImpl == "distance" )
        M_helfrichInnerDivExprImpl = HelfrichInnerDivImplementation::DISTANCE;
    else
        CHECK( false ) << helfrichInnerDivExprImpl << " is not a valid Helfrich inner div implementation\n";

    M_useIntegrationByParts = boption( _name="helfrich-force.use-integration-by-parts", _prefix=this->prefix() );

    std::string const helfrichProjectionMethod = soption( _name="helfrich-force.projection-method", _prefix=this->prefix() );
    if( helfrichProjectionMethod == "nodal-projection" )
        M_helfrichForceMethod = HelfrichMethod::NODAL_PROJECTION;
    else if( helfrichProjectionMethod == "l2-projection" )
        M_helfrichForceMethod = HelfrichMethod::L2_PROJECTION;
    else if( helfrichProjectionMethod == "smooth-projection" )
        M_helfrichForceMethod = HelfrichMethod::SMOOTH_PROJECTION;
    else if( helfrichProjectionMethod == "diffusion" )
        M_helfrichForceMethod = HelfrichMethod::DIFFUSION;
    else
        CHECK( false ) << helfrichProjectionMethod << " is not in the list of possible projection methods\n";

#ifdef DEBUG_HELFRICHFORCEMODEL
    M_exporter = Feel::exporter(
            _mesh=this->levelset()->mesh(),
            _name="HelfrichForce",
            _geo="static",
            _path=this->levelset()->exporterPath()
            );
    M_exporterInitDone = false;
#endif
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesLinearPDE( DataUpdateLinear & data ) const
{
    CHECK( this->fluid() ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesLinearPDE\n";
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    vector_ptrtype& F = data.rhs();
    auto Xh = this->fluid()->functionSpaceVelocity();
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->fluid()->rowStartInVector() );

    auto const& u = this->fluid()->fieldVelocity();
    auto const& v = u;

    if( BuildNonCstPart )
    {
        if( M_useIntegrationByParts )
        {
            auto phi = this->levelset()->phi();
            auto N = this->levelset()->N();
            auto K = this->levelset()->K();
            auto gradPhi = this->levelset()->gradPhi();
            auto modGradPhi = this->levelset()->modGradPhi();
            auto D = this->levelset()->D();

            auto helfrichInnerDiv = helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl );
            myLinearForm +=
                integrate( _range=this->fluid()->rangeMeshElements(),
                        _expr= -this->bendingModulus() * trans(helfrichInnerDiv) * (
                            trans(gradv(D))*(trans(idv(N))*id(v))
                            + idv(D)*( trans(gradv(N))*id(v)+trans(grad(v))*idv(N) )
                            ),
                        _geomap=this->fluid()->geomap() 
                        );
        }
        else
        {
            this->updateInterfaceForcesImpl( this->M_interfaceForce );
            myLinearForm +=
                integrate( _range=this->fluid()->rangeMeshElements(),
                        _expr= trans(idv(*this->M_interfaceForce))*id(v),
                        _geomap=this->fluid()->geomap() );
        }
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesJacobian( DataUpdateJacobian & data ) const
{
    CHECK( this->fluid() ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesJacobian\n";
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesResidual( DataUpdateResidual & data ) const
{
    CHECK( this->fluid() ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesResidual\n";
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    auto Xh = this->fluid()->functionSpaceVelocity();

    auto u = Xh->element(XVec, this->fluid()->rowStartInVector());
    auto const& v = this->fluid()->fieldVelocity();

    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=this->fluid()->rowStartInVector() );
    if( BuildCstPart )
    {
        if( M_useIntegrationByParts )
        {
            auto phi = this->levelset()->phi();
            auto N = this->levelset()->N();
            auto K = this->levelset()->K();
            auto gradPhi = this->levelset()->gradPhi();
            auto modGradPhi = this->levelset()->modGradPhi();
            auto D = this->levelset()->D();

            auto helfrichInnerDiv = helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl );
            linearForm_PatternCoupled +=
                integrate( _range=this->fluid()->rangeMeshElements(),
                        _expr= this->bendingModulus() * trans(helfrichInnerDiv) * (
                            trans(gradv(D))*(trans(idv(N))*id(v))
                            + idv(D)*( trans(gradv(N))*id(v)+trans(grad(v))*idv(N) )
                            ),
                        _geomap=this->fluid()->geomap() 
                        );
        }
        else
        {
            this->updateInterfaceForcesImpl( this->M_interfaceForce );
            linearForm_PatternCoupled +=
                integrate( _range=this->fluid()->rangeMeshElements(),
                        _expr= -trans(idv(*this->M_interfaceForce))*id(v),
                        _geomap=this->fluid()->geomap() );
        }
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::updateInterfaceForcesImpl( element_ptrtype & F ) const
{
    this->addHelfrichForce( F, M_forceImpl );
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HelfrichForceModel<LevelSetType, FluidMechanicsType>::addHelfrichForce( element_ptrtype & F, int impl ) const
{
#ifdef DEBUG_HELFRICHFORCEMODEL
    if( !M_exporterInitDone )
    {
        if (this->levelset()->doRestart() && this->levelset()->restartPath().empty() )
        {
            Feel::cout << "Restarting inextensibility-force exporter...\n";
            if ( M_exporter->doExport() ) M_exporter->restart(this->levelset()->timeInitial());
        }
        M_exporterInitDone = true;
    }
#endif

    switch (impl)
    {
        case 0:
        {
            //auto phi = this->levelset()->phi();
            //auto N = this->levelset()->N();
            //auto K = this->levelset()->K();
            //auto gradPhi = this->levelset()->gradPhi();
            //auto modGradPhi = this->levelset()->modGradPhi();
            auto phi = this->levelset()->distance();
            auto N = this->levelset()->distanceNormal();
            auto K = this->levelset()->distanceCurvature();
            auto gradPhi = this->levelset()->distanceNormal();
            auto modGradPhi = this->levelset()->modGradPhi();

            switch ( M_helfrichForceMethod )
            {
                case HelfrichMethod::NODAL_PROJECTION:
                {
                    auto helfrichInnerDiv = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * divv(helfrichInnerDiv) * idv(this->levelset()->D()) * idv(gradPhi)
                            );
                }
                break;
                case HelfrichMethod::L2_PROJECTION:
                {
                    auto helfrichInnerDiv = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    auto Fb_div = this->levelset()->projectorL2()->project(
                            _expr= divv(helfrichInnerDiv) 
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * idv(Fb_div) * idv(gradPhi) * idv(this->levelset()->D())
                            );
                }
                break;
                case HelfrichMethod::SMOOTH_PROJECTION:
                {
                    auto helfrichInnerDiv = this->levelset()->smootherVectorial()->project(
                            _expr=helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    auto Fb_div = this->levelset()->smoother()->project(
                            _expr=divv(helfrichInnerDiv) 
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * idv(Fb_div) * idv(gradPhi) * idv(this->levelset()->D())
                            );
                }
                break;
                case HelfrichMethod::DIFFUSION:
                {
                    /*
                     * Rk: Helfrich's energy is $int kappa^2$ while Willmore's one is $int H^2 = int kappa^2 / (dim-1)^2$
                     * so that Helfrich force is (dim-1)^2*Kb*W*N*delta
                     */
                    auto W = this->curvatureDiffusionWillmore( mpl::int_<Dim>() );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=(Dim-1)*(Dim-1)*this->bendingModulus()*idv(W)*idv(N)*idv(this->levelset()->D())
                            );
                }
                break;
            }
        }
        break;
        case 1:
        {
            auto phi = this->levelset()->phi();
            auto N = this->levelset()->N();
            auto K = this->levelset()->K();
            auto Id = vf::Id<Dim, Dim>();
            auto NxN = idv(N)*trans(idv(N));

            //auto gradPhi = this->levelset()->smootherVectorial()->project(
                    //trans(gradv(phi))
                    //);

            auto modGradPhi = this->levelset()->modGradPhi();

            auto modGradPhiK = this->levelset()->smoother()->project( 
                    idv(modGradPhi) * idv(K)
                    );

            auto Fb_par = this->levelset()->projectorL2Vectorial()->project(
                    - idv(K)*idv(K)/2. * idv(N)
                    );
            auto Fb_par_div = this->levelset()->smoother()->project(
                    divv(Fb_par)
                    );
            //auto Fb_ortho = this->levelset()->smootherVectorial()->project(
                    //(Id-NxN)*trans(gradv(modGradPhiK)) / sqrt( trans(idv(gradPhi)) * idv(gradPhi) )
                    //);
            auto Fb_ortho = this->levelset()->smootherVectorial()->project(
                    (trans(gradv(modGradPhiK)) - (gradv(modGradPhiK)*idv(N))*idv(N)) / idv(modGradPhi)
                    );
            auto Fb_ortho_div = this->levelset()->smoother()->project(
                    divv(Fb_ortho)
                    );
            //auto Fb_global = this->levelset()->smootherVectorial()->project(
                    //this->bendingModulus() * (divv(Fb_par) + divv(Fb_ortho)) * idv(N)
                    //);
            auto Fb_global = vf::project(
                    _space=this->levelset()->functionSpaceVectorial(),
                    _range=elements(this->levelset()->mesh()),
                    _expr=this->bendingModulus() * (idv(Fb_par_div) + idv(Fb_ortho_div)) * idv(N)
                    );
            //auto Fb_global = this->levelset()->projectorL2Vectorial()->project(
                    //this->bendingModulus() * (divv(Fb_par) + divv(Fb_ortho)) * idv(N)
                    //);

            auto Fb = vf::project(
                    _space=this->levelset()->functionSpaceVectorial(),
                    _range=elements(this->levelset()->mesh()),
                    _expr=idv(Fb_global) * idv(this->levelset()->D())
                    );

            *F = Fb;

#ifdef DEBUG_HELFRICHFORCEMODEL
            M_exporter->step(this->levelset()->time())->add(
                    "dirac", "dirac", this->levelset()->D() );
            auto Fb_par_D = this->levelset()->projectorL2Vectorial()->project(
                    this->bendingModulus() * divv(Fb_par) * idv(N) * idv(this->levelset()->D())
                    );
            M_exporter->step(this->levelset()->time())->add(
                    "helfrich-parallel", "helfrich-parallel", Fb_par_D );
            auto Fb_ortho_D = this->levelset()->projectorL2Vectorial()->project(
                    this->bendingModulus() * divv(Fb_ortho) * idv(N) * idv(this->levelset()->D())
                    );
            M_exporter->step(this->levelset()->time())->add(
                    "helfrich-ortho", "helfrich-ortho", Fb_ortho_D );
            M_exporter->step(this->levelset()->time())->add(
                    "helfrich-force", "helfrich-force", Fb );
            auto curvD = vf::project(
                    _space=this->levelset()->functionSpace(),
                    _range=elements(this->levelset()->mesh()),
                    _expr=idv(K)*idv(this->levelset()->D())
                    );
            M_exporter->step(this->levelset()->time())->add(
                    "curvature_dirac", "curvature_dirac", curvD );
            M_exporter->save();
#endif
        }
        break;
        default:
            CHECK(false) << "Wrong force implementation\n";
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
typename HelfrichForceModel<LevelSetType, FluidMechanicsType>::element_scalar_type
HelfrichForceModel<LevelSetType, FluidMechanicsType>::curvatureDiffusionWillmore( mpl::int_<2> /*Dim*/ ) const
{
    auto dist = this->levelset()->distance();
    double dt = this->levelset()->toolManager()->curvatureDiffusion()->timeStep();
    // Solve Galpha ( alpha=sqrt(2) )
    auto Galpha = this->levelset()->toolManager()->curvatureDiffusion()->solveGalpha( dist );
    // Solve Gbeta ( beta=1/sqrt(2) )
    auto Gbeta = this->levelset()->toolManager()->curvatureDiffusion()->solveGbeta( dist );

    // W = Lap_s H + H^3/2 = 2/dt^2 * (d - aGalpha - bGbeta) - H^3/2
    // with a = -1 and b = 2
    // Rk: in 2D, H=K
    auto H_expr = idv(this->levelset()->distanceCurvature());
    return vf::project( 
            _space=this->levelset()->functionSpace(), _range=this->levelset()->rangeMeshElements(),
            _expr=( idv(dist) + idv(Galpha) - 2*idv(Gbeta) ) * 2 / (dt*dt) - H_expr*H_expr*H_expr/2.
            );
}

template<typename LevelSetType, typename FluidMechanicsType>
typename HelfrichForceModel<LevelSetType, FluidMechanicsType>::element_scalar_type
HelfrichForceModel<LevelSetType, FluidMechanicsType>::curvatureDiffusionWillmore( mpl::int_<3> /*Dim*/ ) const
{
    auto dist = this->levelset()->distance();
    double dt = this->levelset()->toolManager()->curvatureDiffusion()->timeStep();
    // Solve Galpha ( alpha=sqrt(2) )
    auto Galpha = this->levelset()->toolManager()->curvatureDiffusion()->solveGalpha( dist );
    // Solve Gbeta ( beta=1/sqrt(2) )
    auto Gbeta = this->levelset()->toolManager()->curvatureDiffusion()->solveGbeta( dist );
    // Solve Galpha*d^2
    auto dist2 = vf::project( 
            _space=this->levelset()->functionSpace(),
            _range=this->levelset()->rangeMeshElements(),
            _expr=idv(dist)*idv(dist)
            );
    auto G2alpha = this->levelset()->toolManager()->curvatureDiffusion()->solveGalpha( dist2 );
    // Solve Gbeta*d^2
    auto G2beta = this->levelset()->toolManager()->curvatureDiffusion()->solveGbeta( dist2 );

    // W = Lap_s H + 2*H*(H^2-G) = ( (d-aGalpha-bGbeta)(1-Hd) + H/2 (d^2-aG2alpha-bG2beta) )/dt^2
    // with a = -1 and b = 2
    // Rk: in 3D, H=K/2
    auto H_expr = idv(this->levelset()->distanceCurvature())/2.;
    auto d_expr = idv(dist);
    auto d2_expr = idv(dist2);
    return vf::project( 
            _space=this->levelset()->functionSpace(), _range=this->levelset()->rangeMeshElements(),
            _expr=( (d_expr+idv(Galpha)-2*idv(Gbeta))*(1-H_expr*d_expr) + 0.5*H_expr*(d2_expr+idv(G2alpha)-2*idv(G2beta)) )/(dt*dt)
            );
}

} // namespace FeelModels
} // namespace Feel

#endif
