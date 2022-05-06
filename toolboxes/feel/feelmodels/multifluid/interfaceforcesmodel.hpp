#ifndef _INTERFACEFORCESMODEL_HPP
#define _INTERFACEFORCESMODEL_HPP 1

namespace Feel {
namespace FeelModels {

//enum class InterfaceForcesModels {
    //NONE, HELFRICH
//};

template<class LevelSetType, class FluidMechanicsType>
class InterfaceForcesModel
{
    typedef InterfaceForcesModel<LevelSetType, FluidMechanicsType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
public:
    typedef LevelSetType levelset_type;
    typedef std::shared_ptr<levelset_type> levelset_ptrtype;
    
    typedef FluidMechanicsType fluidmechanics_type;
    typedef std::shared_ptr<fluidmechanics_type> fluidmechanics_ptrtype;

    typedef typename levelset_type::space_vectorial_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    typedef typename fluidmechanics_type::DataUpdateJacobian DataUpdateJacobian;
    typedef typename fluidmechanics_type::DataUpdateResidual DataUpdateResidual;
    typedef typename fluidmechanics_type::DataUpdateLinear DataUpdateLinear;

    //--------------------------------------------------------------------//
    // Construction
    InterfaceForcesModel() = default;
    InterfaceForcesModel( InterfaceForcesModel const& i ) = default;
    virtual ~InterfaceForcesModel() = default;

    virtual void build( std::string const& prefix, levelset_ptrtype const& ls, fluidmechanics_ptrtype const& fm = fluidmechanics_ptrtype() );

    virtual std::shared_ptr<std::ostringstream> getInfo() const;

    //--------------------------------------------------------------------//
    std::string const& prefix() const { return M_prefix; }
    levelset_ptrtype const& levelset() const { return M_levelset; }
    fluidmechanics_ptrtype const& fluid() const { return M_fluidmechanics; }
    //--------------------------------------------------------------------//
    template<typename InterfaceForceEltPtrType>
    void updateInterfaceForces( InterfaceForceEltPtrType & F, bool overwrite = false) const;
    void updateInterfaceForces( element_ptrtype & F, bool overwrite = false) const;
    //--------------------------------------------------------------------//
    virtual void updateFluidInterfaceForcesLinearPDE( DataUpdateLinear & data ) const;
    virtual void updateFluidInterfaceForcesJacobian( DataUpdateJacobian & data ) const;
    virtual void updateFluidInterfaceForcesResidual( DataUpdateResidual & data ) const;

    element_ptrtype const& lastInterfaceForce() const { return M_interfaceForce; }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    virtual void updateInterfaceForcesImpl( element_ptrtype & F ) const =0;

    //--------------------------------------------------------------------//
    std::string M_prefix;
    levelset_ptrtype M_levelset;
    fluidmechanics_ptrtype M_fluidmechanics;
protected:
    // Cache last force value
    mutable element_ptrtype M_interfaceForce;
};

template<typename LevelSetType, typename FluidMechanicsType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::build(
        std::string const& prefix,
        levelset_ptrtype const& ls,
        fluidmechanics_ptrtype const& fm
        )
{
    M_prefix = prefix;
    M_levelset = ls;
    M_fluidmechanics = fm;
    M_interfaceForce.reset( new element_type(ls->functionSpaceVectorial(), "InterfaceForce") );
}

template<typename LevelSetType, typename FluidMechanicsType>
std::shared_ptr<std::ostringstream>
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}

template<typename LevelSetType, typename FluidMechanicsType>
template<typename InterfaceForceEltPtrType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::updateInterfaceForces( InterfaceForceEltPtrType & F, bool overwrite ) const
{
    this->updateInterfaceForcesImpl( M_interfaceForce );
    if( overwrite )
    {
        F->on( _range=elements(F->mesh()), _expr=idv(M_interfaceForce) );
    }
    else
    {
        *F += vf::project(
                _space=F->functionSpace(),
                _range=elements(F->mesh()),
                _expr=idv(M_interfaceForce)
                );
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::updateInterfaceForces( element_ptrtype & F, bool overwrite ) const
{
    this->updateInterfaceForcesImpl( M_interfaceForce );
    if( overwrite )
        *F = *M_interfaceForce;
    else
        *F += *M_interfaceForce;
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesLinearPDE( DataUpdateLinear & data ) const
{
    CHECK( M_fluidmechanics ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesLinearPDE\n";
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
        this->updateInterfaceForcesImpl( M_interfaceForce );
        myLinearForm +=
            integrate( _range=this->fluid()->rangeMeshElements(),
                       _expr= trans(idv(*M_interfaceForce))*id(v),
                       _geomap=this->fluid()->geomap() );
    }
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesJacobian( DataUpdateJacobian & data ) const
{
    CHECK( M_fluidmechanics ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesJacobian\n";
}

template<typename LevelSetType, typename FluidMechanicsType>
void
InterfaceForcesModel<LevelSetType, FluidMechanicsType>::updateFluidInterfaceForcesResidual( DataUpdateResidual & data ) const
{
    CHECK( M_fluidmechanics ) << "FluidMechanics model must be provided to use updateFluidInterfaceForcesResidual\n";
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
        this->updateInterfaceForcesImpl( M_interfaceForce );
        linearForm_PatternCoupled +=
            integrate( _range=this->fluid()->rangeMeshElements(),
                    _expr= -trans(idv(*M_interfaceForce))*id(v),
                    _geomap=this->fluid()->geomap() );
    }
}

namespace detail {

template<template <typename, typename> class ModelType, typename LevelSetType, typename FluidMechanicsType>
std::unique_ptr<InterfaceForcesModel<LevelSetType, FluidMechanicsType>> createInterfaceForcesModel()
{
    //return new ModelType<LevelSetType>;
    return std::make_unique<ModelType<LevelSetType, FluidMechanicsType>>();
}

} // namespace detail

} // namespace FeelModels
} // namespace Feel

#endif
