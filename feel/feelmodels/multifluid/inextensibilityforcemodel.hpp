#ifndef _INEXTENSIBILITYFORCEMODEL_HPP
#define _INEXTENSIBILITYFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class InextensibilityForceModel
: public virtual InterfaceForcesModel<LevelSetType>
{
    typedef InextensibilityForceModel<LevelSetType> self_type;
    typedef InterfaceForcesModel<LevelSetType> super_type;
public:
    typedef typename super_type::levelset_type levelset_type;
    typedef typename super_type::levelset_ptrtype levelset_ptrtype;

    typedef typename levelset_type::space_levelset_vectorial_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype;

    //--------------------------------------------------------------------//
    // Construction
    InextensibilityForceModel() = default;
    InextensibilityForceModel( InextensibilityForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls );

    void loadParametersFromOptionsVm();

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F );

    void updateInterfaceRectangularFunction();

    //--------------------------------------------------------------------//
    element_levelset_ptrtype M_levelsetModGradPhi;
    element_levelset_ptrtype M_interfaceRectangularFunction;

    double M_inextensibilityForceCoefficient;
};

template<typename LevelSetType>
void
InextensibilityForceModel<LevelSetType>::build( std::string const& prefix, levelset_ptrtype const& ls )
{
    super_type::build( prefix, ls );

    M_levelsetModGradPhi.reset( new element_levelset_type(this->levelset()->functionSpace(), "ModGradPhi") );
}

template<typename LevelSetType>
void
InextensibilityForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_inextensibilityForceCoefficient = doption( _name="inextensibility-force-coeff", _prefix=this->prefix() );
}

template<typename LevelSetType>
void
InextensibilityForceModel<LevelSetType>::updateInterfaceRectangularFunction()
{
    if( !M_interfaceRectangularFunction )
        M_interfaceRectangularFunction.reset( new element_levelset_type(this->levelset()->functionSpace(), "InterfaceRectangularFunction") );

    auto phi = idv(this->levelset()->phi());
    double epsilon = this->levelset()->thicknessInterface();
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
            this->levelset()->functionSpace(), 
            elements(this->levelset()->mesh()),
            R_expr
            );
}

template<typename LevelSetType>
void
InextensibilityForceModel<LevelSetType>::updateInterfaceForcesImpl( element_ptrtype & F )
{
    // Update ModGradPhi
    auto gradPhi = this->levelset()->gradPhi();
    if( this->levelset()->hasReinitialized() )
    {
        this->updateInterfaceRectangularFunction();
        *M_levelsetModGradPhi = this->levelset()->projectorL2()->project(
                _expr = sqrt( trans(idv(gradPhi))*idv(gradPhi) ) * (1.-idv(M_interfaceRectangularFunction))
                        + idv(M_levelsetModGradPhi)*idv(M_interfaceRectangularFunction)
                );
    }
    else
    {
        *M_levelsetModGradPhi = this->levelset()->projectorL2()->project(
                _expr = sqrt( trans(idv(gradPhi))*idv(gradPhi) )
                );
    }

    auto Ep = this->levelset()->projectorL2()->project(
            _expr=max( idv(M_levelsetModGradPhi)-cst(1.), 0. )
            );
    auto EpN = this->levelset()->projectorL2Vectorial()->project(
            _expr=idv(Ep)*idv(this->levelset()->N())
            );
    auto GradEp = this->levelset()->smootherVectorial()->project(
            _expr=trans(gradv(Ep))
            );
    auto DivEpNtN = this->levelset()->smootherVectorial()->project(
            _expr=divv(EpN)*idv(this->levelset()->N())
            );

    *F += vf::project(
            _space=this->levelset()->functionSpaceVectorial(),
            _range=elements(this->levelset()->mesh()),
            _expr=M_inextensibilityForceCoefficient*(idv(GradEp)-idv(DivEpNtN))*idv(this->levelset()->D())
            );
}

} // namespace FeelModels
} // namespace Feel

#endif
