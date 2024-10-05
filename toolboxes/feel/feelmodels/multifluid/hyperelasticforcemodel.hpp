#ifndef _HYPERELASTICFORCEMODEL_HPP
#define _HYPERELASTICFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType, class FluidMechanicsType>
class HyperelasticForceModel
: public virtual InterfaceForcesModel<LevelSetType, FluidMechanicsType>
{
    typedef HyperelasticForceModel<LevelSetType, FluidMechanicsType> self_type;
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

    typedef typename levelset_type::element_cauchygreen_invariant_type element_cauchygreen_invariant_type;
    typedef typename levelset_type::element_cauchygreen_invariant_ptrtype element_cauchygreen_invariant_ptrtype;

    typedef element_cauchygreen_invariant_type element_energyderivative_type;
    typedef element_cauchygreen_invariant_ptrtype element_energyderivative_ptrtype;

    static inline const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    HyperelasticForceModel() = default;
    HyperelasticForceModel( HyperelasticForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm = fluidmechanics_ptrtype() ) override;
    void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm, std::string const& name );

    void loadParametersFromOptionsVm();

    //--------------------------------------------------------------------//
    element_energyderivative_ptrtype const& energyDerivative1() const;
    element_energyderivative_ptrtype const& energyDerivative2() const;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F ) const override;
    virtual element_energyderivative_ptrtype const& energyDerivative1Impl() const =0;
    virtual element_energyderivative_ptrtype const& energyDerivative2Impl() const =0;

    std::string M_name;

protected:
#ifdef DEBUG_HYPERELASTICFORCEMODEL
    typedef std::shared_ptr<Exporter<mesh_type, 1>> exporter_ptrtype;
    exporter_ptrtype M_exporter;
    bool M_exporterInitDone;
#endif
};

template<typename LevelSetType, typename FluidMechanicsType>
void
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm )
{
    this->build( prefix, ls, fm, "");
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm, std::string const& name )
{
    super_type::build( prefix, ls, fm );
    M_name = name;
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::loadParametersFromOptionsVm()
{
#ifdef DEBUG_HYPERELASTICFORCEMODEL
    M_exporter = Feel::exporter(
            _mesh=this->levelset()->mesh(),
            _name=M_name,
            _geo="static",
            _path=this->levelset()->exporterPath()
            );
    M_exporterInitDone = false;
#endif
}

template<typename LevelSetType, typename FluidMechanicsType>
typename HyperelasticForceModel<LevelSetType, FluidMechanicsType>::element_energyderivative_ptrtype const&
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::energyDerivative1() const
{
    return this->energyDerivative1Impl();
}

template<typename LevelSetType, typename FluidMechanicsType>
typename HyperelasticForceModel<LevelSetType, FluidMechanicsType>::element_energyderivative_ptrtype const&
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::energyDerivative2() const
{
    return this->energyDerivative2Impl();
}

template<typename LevelSetType, typename FluidMechanicsType>
void
HyperelasticForceModel<LevelSetType, FluidMechanicsType>::updateInterfaceForcesImpl( element_ptrtype & F ) const
{
#ifdef DEBUG_HYPERELASTICFORCEMODEL
    if( !M_exporterInitDone )
    {
        if (this->levelset()->doRestart() && this->levelset()->restartPath().empty() )
        {
            Feel::cout << "Restarting " + M_name + " exporter...\n";
            if ( M_exporter->doExport() ) M_exporter->restart(this->levelset()->timeInitial());
        }
        M_exporterInitDone = true;
    }
#endif

    //auto N = this->levelset()->N();
    //auto K = this->levelset()->K();
    //auto Id = vf::Id<Dim, Dim>();
    //auto NxN = idv(N)*trans(idv(N));

    // Left Cauchy-Green strain tensor
    auto A = idv(this->levelset()->leftCauchyGreenTensor());

    // Derivative of the energy wrt the first Cauchy-Green invariant (TrC)
    auto Ep1 = this->energyDerivative1();
    // Derivative of the energy wrt the second Cauchy-Green invariant (TrCofC)
    auto Ep2 = this->energyDerivative2();

    //auto Fe_expr = this->levelset()->projectorL2Tensor2Symm()->project(
            //_expr=2.*(idv(Ep1)*A + idv(Ep2)*(A*trace(A)-A*A) )
            //);
#if 0
    auto Fe_expr = vf::project(
            _space=this->levelset()->functionSpaceTensor2Symm(),
            _expr=2.*(idv(Ep1)*A + idv(Ep2)*(A*trace(A)-A*A) )
            );
#elif 0 // New implementation
    auto N = this->levelset()->N();
    auto Z1 = idv( this->levelset()->cauchyGreenInvariant1() );
    auto C1 = vf::Id<Dim, Dim>() - idv(N)*trans(idv(N));
    auto Z2 = idv( this->levelset()->cauchyGreenInvariant2() );
    auto C2 = 2.*A/trace(A) - C1;
    auto Fe_expr = this->levelset()->functionSpaceTensor2Symm()->element();
    Fe_expr.on(
            _range=this->levelset()->interfaceElements(),
            _expr=idv(Ep1)*Z1*C1 + idv(Ep2)*Z2*C2
            );
#elif 1 // Old implementation
    auto N = this->levelset()->N();
    auto Z1 = idv( this->levelset()->cauchyGreenInvariant1() );
    auto C1 = vf::Id<Dim, Dim>() - idv(N)*trans(idv(N));
    auto M1 = Z1*C1;
    auto M2 = A;
    auto Fe_expr = this->levelset()->functionSpaceTensor2Symm()->element();
    Fe_expr.on(
            //_range=this->levelset()->interfaceElements(),
            _range=elements(this->levelset()->mesh()),
            _expr=idv(Ep1)*M1 + idv(Ep2)*M2
            );
#endif
    //auto Fe = this->levelset()->smootherInterfaceVectorial()->project(
    auto Fe = this->levelset()->smootherVectorial()->project(
    //auto Fe = this->levelset()->projectorL2Vectorial()->project(
            _expr=divv(Fe_expr)
            );
    
    F->zero();
    F->on(
        //_range=elements(this->levelset()->mesh()),
        _range=this->levelset()->rangeDiracElements(),
        _expr=idv(Fe)*idv(this->levelset()->D())
        );

#ifdef DEBUG_HYPERELASTICFORCEMODEL
    M_exporter->step(this->levelset()->time())->add(
            "Ep1", "Ep1", Ep1 );
    M_exporter->step(this->levelset()->time())->add(
            "Ep2", "Ep2", Ep2 );
    M_exporter->step(this->levelset()->time())->add(
            "F", "F", Fe );
    M_exporter->save();
#endif
}

} // namespace FeelModels
} // namespace Feel

#endif //_HYPERELASTICFORCEMODEL_HPP
