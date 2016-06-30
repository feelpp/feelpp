#ifndef _INTERFACEFORCESMODEL_HPP
#define _INTERFACEFORCESMODEL_HPP 1

namespace Feel {
namespace FeelModels {

template<class SpaceType>
class InterfaceForcesModel
{
    typedef InterfaceForcesModel<SpaceType> self_type;
public:
    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type mesh_ptrtype mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Construction
    InterfaceForcesModel( std::string const& prefix );

};

} // namespace FeelModels
} // namespace Feel

#endif
