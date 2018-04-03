#ifndef _LINEARELASTICFORCEMODEL_HPP
#define _LINEARELASTICFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/hyperelasticforcemodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class LinearElasticForceModel
: public HyperelasticForceModel<LevelSetType>
{
    typedef LinearElasticForceModel<LevelSetType> self_type;
    typedef HyperelasticForceModel<LevelSetType> super_type;
public:
    typedef typename super_type::levelset_type levelset_type;
    typedef typename super_type::levelset_ptrtype levelset_ptrtype;

    typedef typename super_type::mesh_type mesh_type;

    typedef typename levelset_type::space_levelset_vectorial_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype;

    typedef typename levelset_type::element_cauchygreen_invariant_type element_cauchygreen_invariant_type;
    typedef typename levelset_type::element_cauchygreen_invariant_ptrtype element_cauchygreen_invariant_ptrtype;

    typedef typename super_type::element_cauchygreen_invariant_type element_energyderivative_type;
    typedef typename super_type::element_cauchygreen_invariant_ptrtype element_energyderivative_ptrtype;

    static const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    LinearElasticForceModel() = default;
    LinearElasticForceModel( LinearElasticForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls );
    void loadParametersFromOptionsVm();

    boost::shared_ptr<std::ostringstream> getInfo() const;

private:
    virtual element_energyderivative_ptrtype const& energyDerivative1Impl() const;
    virtual element_energyderivative_ptrtype const& energyDerivative2Impl() const;

    mutable element_energyderivative_ptrtype M_energyDerivative1;
    mutable element_energyderivative_ptrtype M_energyDerivative2;

    double M_elasticStretchModulus;
    double M_elasticShearModulus;

};

template<typename LevelSetType>
void
LinearElasticForceModel<LevelSetType>::build( std::string const& prefix, levelset_ptrtype const& ls )
{
    super_type::build( prefix, ls, "linear-elastic-force" );
    this->loadParametersFromOptionsVm();
}

template<typename LevelSetType>
void
LinearElasticForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_elasticStretchModulus = doption( _name="elastic-stretch-modulus", _prefix=this->prefix() );
    M_elasticShearModulus = doption( _name="elastic-shear-modulus", _prefix=this->prefix() );
}

template<typename LevelSetType>
boost::shared_ptr<std::ostringstream> 
LinearElasticForceModel<LevelSetType>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "Linear elastic force ("
           << "stretch modulus = " << this->M_elasticStretchModulus << ", "
           << "shear modulus = " << this->M_elasticShearModulus
           << ")";

    return _ostr;
}

template<typename LevelSetType>
typename LinearElasticForceModel<LevelSetType>::element_energyderivative_ptrtype const&
LinearElasticForceModel<LevelSetType>::energyDerivative1Impl() const
{
    if( !M_energyDerivative1 )
        M_energyDerivative1.reset( new element_energyderivative_type(this->levelset()->functionSpace(), "EnergyDerivative1") );

    auto I1 = idv(this->levelset()->cauchyGreenInvariant1());

#if 0
    *M_energyDerivative1 = vf::project(
        _space=M_energyDerivative1->functionSpace(),
        _expr=M_elasticShearModulus * (I1-2.)
        );
#else // New implementation
    *M_energyDerivative1 = vf::project(
        _space=M_energyDerivative1->functionSpace(),
        _expr=M_elasticStretchModulus * (I1-1.)
        );
#endif
    return M_energyDerivative1;
}

template<typename LevelSetType>
typename LinearElasticForceModel<LevelSetType>::element_energyderivative_ptrtype const&
LinearElasticForceModel<LevelSetType>::energyDerivative2Impl() const
{
    if( !M_energyDerivative2 )
        M_energyDerivative2.reset( new element_energyderivative_type(this->levelset()->functionSpace(), "EnergyDerivative2") );

    auto I2 = idv(this->levelset()->cauchyGreenInvariant2());

#if 0
    *M_energyDerivative2 = vf::project(
        _space=M_energyDerivative2->functionSpace(),
        _expr=M_elasticStretchModulus * (I2-1.)
        );
#else // New implementation
    *M_energyDerivative2 = vf::project(
        _space=M_energyDerivative2->functionSpace(),
        _expr=M_elasticShearModulus * (I2-1.)
        );
#endif
    return M_energyDerivative2;
}

}
}

#endif //_LINEARELASTICFORCEMODEL_HPP
