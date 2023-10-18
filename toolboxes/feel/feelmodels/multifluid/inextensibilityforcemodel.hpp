#ifndef _INEXTENSIBILITYFORCEMODEL_HPP
#define _INEXTENSIBILITYFORCEMODEL_HPP 1

//#define DEBUG_INEXTENSIBILITYFORCEMODEL

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType, class FluidMechanicsType>
class InextensibilityForceModel
: public virtual InterfaceForcesModel<LevelSetType, FluidMechanicsType>
{
    typedef InextensibilityForceModel<LevelSetType, FluidMechanicsType> self_type;
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

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype;

    static inline const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    InextensibilityForceModel() = default;
    InextensibilityForceModel( InextensibilityForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm = fluidmechanics_ptrtype() ) override;

    std::shared_ptr<std::ostringstream> getInfo() const override;

    void loadParametersFromOptionsVm();

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F ) const override;

    void updateInterfaceRectangularFunction();

    //--------------------------------------------------------------------//
    element_levelset_ptrtype M_levelsetModGradPhi;
    element_levelset_ptrtype M_interfaceRectangularFunction;

    double M_inextensibilityForceCoefficient;
    double M_inextensibilityForceEpsilon;

#ifdef DEBUG_INEXTENSIBILITYFORCEMODEL
    typedef std::shared_ptr<Exporter<mesh_type, 1>> exporter_ptrtype;
    exporter_ptrtype M_exporter;
    bool M_exporterInitDone;
#endif
};

template<typename LevelSetType, typename FluidMechanicsType>
void
InextensibilityForceModel<LevelSetType, FluidMechanicsType>::build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm )
{
    super_type::build( prefix, ls, fm );
    this->loadParametersFromOptionsVm();

    M_levelsetModGradPhi.reset( new element_levelset_type(this->levelset()->functionSpace(), "ModGradPhi") );
}

template<typename LevelSetType, typename FluidMechanicsType>
std::shared_ptr<std::ostringstream> 
InextensibilityForceModel<LevelSetType, FluidMechanicsType>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "Inextensibility force ("
           << "coeff = " << this->M_inextensibilityForceCoefficient
           << ")";

    return _ostr;
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InextensibilityForceModel<LevelSetType, FluidMechanicsType>::loadParametersFromOptionsVm()
{
    M_inextensibilityForceCoefficient = doption( _name="inextensibility-force-coeff", _prefix=this->prefix() );

    if( Environment::vm().count( prefixvm(this->prefix(),"inextensibility-force-epsilon").c_str() ) )
        M_inextensibilityForceEpsilon = doption( _name="inextensibility-force-epsilon", _prefix=this->prefix() );
    else
        M_inextensibilityForceEpsilon = this->levelset()->thicknessInterface();
    
#ifdef DEBUG_INEXTENSIBILITYFORCEMODEL
    M_exporter = Feel::exporter(
            _mesh=this->levelset()->mesh(),
            _name="InextensibilityForce",
            _geo="static",
            _path=this->levelset()->exporterPath()
            );
    M_exporterInitDone = false;
#endif
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InextensibilityForceModel<LevelSetType, FluidMechanicsType>::updateInterfaceRectangularFunction()
{
    if( !M_interfaceRectangularFunction )
        M_interfaceRectangularFunction.reset( new element_levelset_type(this->levelset()->functionSpace(), "InterfaceRectangularFunction") );

    auto phi = idv(this->levelset()->phi());
    double epsilon = M_inextensibilityForceEpsilon;
    double epsilon_rect = 2.*epsilon;
    double epsilon_delta = (epsilon_rect - epsilon)/2.;
    double epsilon_zero = epsilon + epsilon_delta;

    auto R_expr =
        vf::chi( phi<-epsilon_rect )*vf::constant(0.0)
        +
        vf::chi( phi>=-epsilon_rect )*vf::chi( phi<=-epsilon )*
        1/2*(1 + (phi+epsilon_zero)/epsilon_delta + 1/M_PI*vf::sin( M_PI*(phi+epsilon_zero)/epsilon_delta ) )
        +
        vf::chi( phi>=-epsilon )*vf::chi( phi<=epsilon )*vf::constant(1.0)
        +
        vf::chi( phi>=epsilon )*vf::chi( phi<=epsilon_rect )*
        1/2*(1 - (phi-epsilon_zero)/epsilon_delta - 1/M_PI*vf::sin( M_PI*(phi-epsilon_zero)/epsilon_delta ) )
        +
        vf::chi(phi>epsilon_rect)*vf::constant(0.0)
        ;

    *M_interfaceRectangularFunction = vf::project( 
        _space=this->levelset()->functionSpace(), 
        _range=elements(this->levelset()->mesh()),
        _expr=R_expr
            );
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InextensibilityForceModel<LevelSetType, FluidMechanicsType>::updateInterfaceForcesImpl( element_ptrtype & F ) const
{
#ifdef DEBUG_INEXTENSIBILITYFORCEMODEL
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

    // Update ModGradPhi
    //auto gradPhi = this->levelset()->gradPhi();
    auto phi = this->levelset()->phi();
    //auto gradPhi = this->levelset()->smootherVectorial()->project(
            //trans(gradv(phi))
            //);
    auto N = this->levelset()->N();
    auto K = this->levelset()->K();
    auto Id = vf::Id<Dim, Dim>();
    auto NxN = idv(N)*trans(idv(N));

    //if( this->levelset()->hasReinitialized() )
    //{
        //this->updateInterfaceRectangularFunction();
        ///[>M_levelsetModGradPhi = this->levelset()->projectorL2()->project(
        //*M_levelsetModGradPhi = this->levelset()->smoother()->project(
                ////_expr = sqrt( trans(idv(gradPhi))*idv(gradPhi) ) * (1.-idv(M_interfaceRectangularFunction))
                        ////+ idv(M_levelsetModGradPhi)*idv(M_interfaceRectangularFunction)
                //_expr = idv(this->levelset()->modGradPhi()) * (1.-idv(M_interfaceRectangularFunction))
                        //+ idv(M_levelsetModGradPhi)*idv(M_interfaceRectangularFunction)
                //);
    //}
    //else
    //{
        ///[>M_levelsetModGradPhi = this->levelset()->projectorL2()->project(
        //*M_levelsetModGradPhi = this->levelset()->smoother()->project(
                ////_expr = sqrt( trans(idv(gradPhi))*idv(gradPhi) )
                //_expr = idv(this->levelset()->modGradPhi())
                //);
    //}

    //auto Ep = this->levelset()->projectorL2()->project(
            //_expr=max( idv(M_levelsetModGradPhi)-cst(1.), 0. )
            //);
    //auto Ep = this->levelset()->projectorL2()->project(
            //_expr=idv(M_levelsetModGradPhi)-cst(1.)
            //);
    //auto Ep = *M_levelsetModGradPhi;
    //Ep.add (-1.);
    //auto Ep = this->levelset()->stretch();
    //auto Ep = this->levelset()->smoother()->project(
            //_expr=idv(this->levelset()->stretch())
            //);
    //auto EpRaw = this->levelset()->stretch();
    //auto Ep = this->levelset()->smoother()->project(
            //_expr=max(idv(EpRaw), 0.)
            //);
    auto EpRaw = this->levelset()->stretch();
    auto Ep = vf::project(
            _space=this->levelset()->functionSpace(),
            _range=elements(this->levelset()->mesh()),
            _expr=max(idv(EpRaw), 0)
            );
    //auto EpN = this->levelset()->projectorL2Vectorial()->project(
            //_expr=idv(Ep)*idv(N)
            //);
    //auto GradEp = this->levelset()->smootherVectorial()->project(
            //_expr=trans(gradv(Ep))
            //);
    auto GradEp = this->levelset()->projectorL2Vectorial()->project(
            _expr=trans(gradv(Ep))
            );
    //auto DivEpNtN = this->levelset()->smootherVectorial()->project(
            //_expr=divv(EpN)*idv(N)
            //);

    //auto Fe = vf::project(
            //_space=this->levelset()->functionSpaceVectorial(),
            //_range=elements(this->levelset()->mesh()),
            //_expr=M_inextensibilityForceCoefficient*(idv(GradEp)-idv(DivEpNtN))*idv(this->levelset()->D())
            //);
    //auto Fe = vf::project(
            //_space=this->levelset()->functionSpaceVectorial(),
            //_range=elements(this->levelset()->mesh()),
            //_expr=M_inextensibilityForceCoefficient*( (Id-NxN)*idv(GradEp)-idv(Ep)*idv(K)*idv(N) )*idv(this->levelset()->D())
            //);
    //auto Fe = this->levelset()->smootherVectorial()->project(
            //_expr=M_inextensibilityForceCoefficient*( (Id-NxN)*idv(GradEp)-idv(Ep)*idv(K)*idv(N) )*idv(this->levelset()->D())
            //);
    auto Fe = this->levelset()->projectorL2Vectorial()->project(
            _expr=M_inextensibilityForceCoefficient*( (Id-NxN)*idv(GradEp)-idv(Ep)*idv(K)*idv(N) )*idv(this->levelset()->D())
            );
    *F = Fe;

#ifdef DEBUG_INEXTENSIBILITYFORCEMODEL
    M_exporter->step(this->levelset()->time())->add(
            "gradEp", "gradEp", GradEp );
    M_exporter->step(this->levelset()->time())->add(
            "Ep", "Ep", Ep );
    M_exporter->step(this->levelset()->time())->add(
            "inextensibility-force", "inextensibility-force", Fe );
    auto Fe_ortho = this->levelset()->projectorL2Vectorial()->project(
            _expr=( (Id-NxN)*idv(GradEp) )*idv(this->levelset()->D())
            );
    M_exporter->step(this->levelset()->time())->add(
            "inextensibility-force-ortho", "inextensibility-force-ortho", Fe_ortho );
    auto Fe_par = this->levelset()->projectorL2Vectorial()->project(
            _expr=( -idv(Ep)*idv(K)*idv(N) )*idv(this->levelset()->D())
            );
    M_exporter->step(this->levelset()->time())->add(
            "inextensibility-force-parallel", "inextensibility-force-parallel", Fe_par );
    M_exporter->save();
#endif
}

} // namespace FeelModels
} // namespace Feel

#endif
