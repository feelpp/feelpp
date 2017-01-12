#ifndef _INTERFACEFORCESMODEL_HPP
#define _INTERFACEFORCESMODEL_HPP 1

namespace Feel {
namespace FeelModels {

//enum class InterfaceForcesModels {
    //NONE, HELFRICH
//};

template<class LevelSetType>
class InterfaceForcesModel
{
    typedef InterfaceForcesModel<LevelSetType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
public:
    typedef LevelSetType levelset_type;
    typedef boost::shared_ptr<levelset_type> levelset_ptrtype;

    typedef typename levelset_type::space_levelset_vectorial_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    //--------------------------------------------------------------------//
    // Construction
    InterfaceForcesModel() = default;
    InterfaceForcesModel( InterfaceForcesModel const& i ) = default;
    virtual ~InterfaceForcesModel() = default;

    virtual void build( std::string const& prefix, levelset_ptrtype const& ls );

    virtual boost::shared_ptr<std::ostringstream> getInfo() const;

    //--------------------------------------------------------------------//
    std::string const& prefix() const { return M_prefix; }
    levelset_ptrtype const& levelset() const { return M_levelset; }
    //--------------------------------------------------------------------//
    void updateInterfaceForces( element_ptrtype & F, bool overwrite = false);

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    virtual void loadParametersFromOptionsVm() =0;
    virtual void updateInterfaceForcesImpl( element_ptrtype & F ) =0;

    //--------------------------------------------------------------------//
    std::string M_prefix;

    levelset_ptrtype M_levelset;
};

template<typename LevelSetType>
void
InterfaceForcesModel<LevelSetType>::build(
        std::string const& prefix,
        levelset_ptrtype const& ls 
        )
{
    M_prefix = prefix;
    M_levelset = ls;
    this->loadParametersFromOptionsVm();
}

template<typename LevelSetType>
boost::shared_ptr<std::ostringstream>
InterfaceForcesModel<LevelSetType>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}

template<typename LevelSetType>
void
InterfaceForcesModel<LevelSetType>::updateInterfaceForces( element_ptrtype & F, bool overwrite )
{
    if( overwrite )
        F->zero();
    this->updateInterfaceForcesImpl( F );
}

namespace detail {

template<template <typename> class ModelType, typename LevelSetType>
InterfaceForcesModel<LevelSetType>* createInterfaceForcesModel()
{
    return new ModelType<LevelSetType>;
}

} // namespace detail

} // namespace FeelModels
} // namespace Feel

#endif
