#ifndef _SKALAKFORCEMODEL_HPP
#define _SKALAKFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/hyperelasticforcemodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType, class FluidMechanicsType>
class SkalakForceModel
: public HyperelasticForceModel<LevelSetType, FluidMechanicsType>
{
    typedef SkalakForceModel<LevelSetType, FluidMechanicsType> self_type;
    typedef HyperelasticForceModel<LevelSetType, FluidMechanicsType> super_type;
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

    typedef typename levelset_type::element_cauchygreen_invariant_type element_cauchygreen_invariant_type;
    typedef typename levelset_type::element_cauchygreen_invariant_ptrtype element_cauchygreen_invariant_ptrtype;

    typedef typename super_type::element_cauchygreen_invariant_type element_energyderivative_type;
    typedef typename super_type::element_cauchygreen_invariant_ptrtype element_energyderivative_ptrtype;

    static inline const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    SkalakForceModel() = default;
    SkalakForceModel( SkalakForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm = fluidmechanics_ptrtype() );
    void loadParametersFromOptionsVm();

    std::shared_ptr<std::ostringstream> getInfo() const;

private:
    virtual element_energyderivative_ptrtype const& energyDerivative1Impl() const;
    virtual element_energyderivative_ptrtype const& energyDerivative2Impl() const;

    mutable element_energyderivative_ptrtype M_energyDerivative1;
    mutable element_energyderivative_ptrtype M_energyDerivative2;

    double M_skalakForceStretchModulus;
    double M_skalakForceShearModulus;

};

template<typename LevelSetType, typename FluidMechanicsType>
void
SkalakForceModel<LevelSetType, FluidMechanicsType>::build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm )
{
    super_type::build( prefix, ls, fm, "skalak-force" );
    this->loadParametersFromOptionsVm();
}

template<typename LevelSetType, typename FluidMechanicsType>
void
SkalakForceModel<LevelSetType, FluidMechanicsType>::loadParametersFromOptionsVm()
{
    M_skalakForceStretchModulus = doption( _name="skalak-stretch-modulus", _prefix=this->prefix() );
    M_skalakForceShearModulus = doption( _name="skalak-shear-modulus", _prefix=this->prefix() );
}

template<typename LevelSetType, typename FluidMechanicsType>
std::shared_ptr<std::ostringstream> 
SkalakForceModel<LevelSetType, FluidMechanicsType>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "Skalak force ("
           << "stretch modulus C = " << this->M_skalakForceStretchModulus << ", "
           << "shear modulus B = " << this->M_skalakForceShearModulus
           << ")";

    return _ostr;
}

template<typename LevelSetType, typename FluidMechanicsType>
typename SkalakForceModel<LevelSetType, FluidMechanicsType>::element_energyderivative_ptrtype const&
SkalakForceModel<LevelSetType, FluidMechanicsType>::energyDerivative1Impl() const
{
    if( !M_energyDerivative1 )
        M_energyDerivative1.reset( new element_energyderivative_type(this->levelset()->functionSpace(), "EnergyDerivative1") );

    auto i1 = idv(this->levelset()->cauchyGreenInvariant1());

    *M_energyDerivative1 = vf::project(
        _space=M_energyDerivative1->functionSpace(),
        _expr=0.5*i1*(M_skalakForceStretchModulus * (i1*i1 - 1) - M_skalakForceShearModulus)
        );

    return M_energyDerivative1;
}

template<typename LevelSetType, typename FluidMechanicsType>
typename SkalakForceModel<LevelSetType, FluidMechanicsType>::element_energyderivative_ptrtype const&
SkalakForceModel<LevelSetType, FluidMechanicsType>::energyDerivative2Impl() const
{
    if( !M_energyDerivative2 )
        M_energyDerivative2.reset( new element_energyderivative_type(this->levelset()->functionSpace(), "EnergyDerivative2") );

    auto i2 = idv(this->levelset()->cauchyGreenInvariant2());

    *M_energyDerivative2 = vf::project(
        _space=M_energyDerivative2->functionSpace(),
        _expr=M_skalakForceShearModulus * (i2-0.5)
        );
    return M_energyDerivative2;
}

}
}

#endif //_SKALAKFORCEMODEL_HPP
