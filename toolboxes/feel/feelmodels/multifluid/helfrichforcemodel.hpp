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

    typedef typename fluidmechanics_type::DataUpdateJacobian DataUpdateJacobian;
    typedef typename fluidmechanics_type::DataUpdateResidual DataUpdateResidual;
    typedef typename fluidmechanics_type::DataUpdateLinear DataUpdateLinear;

    static const uint16_type Dim = levelset_type::nDim;

    enum class ProjectionMethod {
        NODAL_PROJECTION, L2_PROJECTION, SMOOTH_PROJECTION
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

    //--------------------------------------------------------------------//
    double M_helfrichBendingModulus;
    int M_forceImpl;
    vf::FeelModels::HelfrichInnerDivImplementation M_helfrichInnerDivExprImpl;
    bool M_useIntegrationByParts;
    ProjectionMethod M_helfrichForceProjectionMethod;

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
        M_helfrichInnerDivExprImpl = vf::FeelModels::HelfrichInnerDivImplementation::GENERIC;
    else if( helfrichInnerDivExprImpl == "distance" )
        M_helfrichInnerDivExprImpl = vf::FeelModels::HelfrichInnerDivImplementation::DISTANCE;
    else
        CHECK( false ) << helfrichInnerDivExprImpl << " is not a valid Helfrich inner div implementation\n";

    M_useIntegrationByParts = boption( _name="helfrich-force.use-integration-by-parts", _prefix=this->prefix() );

    std::string const helfrichProjectionMethod = soption( _name="helfrich-force.projection-method", _prefix=this->prefix() );
    if( helfrichProjectionMethod == "nodal-projection" )
        M_helfrichForceProjectionMethod = ProjectionMethod::NODAL_PROJECTION;
    else if( helfrichProjectionMethod == "l2-projection" )
        M_helfrichForceProjectionMethod = ProjectionMethod::L2_PROJECTION;
    else if( helfrichProjectionMethod == "smooth-projection" )
        M_helfrichForceProjectionMethod = ProjectionMethod::SMOOTH_PROJECTION;
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
    auto Xh = this->fluid()->functionSpace();
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->fluid()->rowStartInVector() );

    auto const& U = this->fluid()->fieldVelocityPressure();
    auto u = U.template element<0>();
    auto v = U.template element<0>();

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

            auto helfrichInnerDiv = Feel::vf::FeelModels::helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl );
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
    auto Xh = this->fluid()->functionSpace();

    auto U = Xh->element(XVec, this->fluid()->rowStartInVector());
    auto u = U.template element<0>();
    auto v = U.template element<0>();

    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=this->fluid()->rowStartInVector() );
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

            auto helfrichInnerDiv = Feel::vf::FeelModels::helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl );
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
            auto phi = this->levelset()->phi();
            auto N = this->levelset()->N();
            auto K = this->levelset()->K();
            auto gradPhi = this->levelset()->gradPhi();
            auto modGradPhi = this->levelset()->modGradPhi();

            switch ( M_helfrichForceProjectionMethod )
            {
                case ProjectionMethod::NODAL_PROJECTION:
                {
                    auto helfrichInnerDiv = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=Feel::vf::FeelModels::helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * divv(helfrichInnerDiv) * idv(this->levelset()->D()) * idv(gradPhi)
                            );
                }
                break;
                case ProjectionMethod::L2_PROJECTION:
                {
                    auto helfrichInnerDiv = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=Feel::vf::FeelModels::helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    auto Fb = this->levelset()->projectorL2Vectorial()->project(
                            _expr= divv(helfrichInnerDiv) * idv(gradPhi)
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * idv(Fb) * idv(this->levelset()->D())
                            );
                }
                break;
                case ProjectionMethod::SMOOTH_PROJECTION:
                {
                    auto helfrichInnerDiv = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=Feel::vf::FeelModels::helfrichInnerDivExpr( *N, *K, *modGradPhi, M_helfrichInnerDivExprImpl )
                            );
                    auto Fb = this->levelset()->smootherVectorial()->project(
                            _expr=divv(helfrichInnerDiv) * idv(gradPhi)
                            );
                    *F = vf::project(
                            _space=this->levelset()->functionSpaceVectorial(),
                            _range=this->levelset()->rangeMeshElements(),
                            _expr=this->bendingModulus() * idv(Fb) * idv(this->levelset()->D())
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

} // namespace FeelModels
} // namespace Feel

#endif
