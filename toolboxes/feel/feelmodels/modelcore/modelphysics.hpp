/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H 1

#include <type_traits> // TO ADD
#include <feel/feelcore/feel.hpp>
//#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelcore/traits.hpp>
#include <string>
#include <set>
#include <map>
#include <vector>

namespace Feel
{
namespace FeelModels
{

class MaterialPropertyDescription : public std::tuple<std::string,std::vector<std::pair<uint16_type,uint16_type>>>
{
    using super_type = std::tuple<std::string,std::vector<std::pair<uint16_type,uint16_type>>>;
  public :
    using shape_dim_type = std::pair<uint16_type,uint16_type>;
    using shapes_dim_type = std::vector<shape_dim_type>;

    MaterialPropertyDescription() = default;
    MaterialPropertyDescription( std::string const& symbol ) : super_type( symbol, shapes_dim_type{} ) {}
    MaterialPropertyDescription( std::string const& symbol, shape_dim_type const& shape )
        :
        super_type( symbol, shapes_dim_type{} )
    {
        this->add( shape );
    }
    template <typename ContainerShapesType>
    MaterialPropertyDescription( std::string const& symbol, ContainerShapesType const& shapes )
        :
        super_type( symbol, shapes_dim_type{} )
    {
        this->add( shapes );
    }
    MaterialPropertyDescription( MaterialPropertyDescription const& ) = default;
    MaterialPropertyDescription( MaterialPropertyDescription && ) = default;
    MaterialPropertyDescription& operator=( MaterialPropertyDescription const& ) = default;
    MaterialPropertyDescription& operator=( MaterialPropertyDescription && ) = default;

    std::string const& symbol() const { return std::get<0>(*this); }
    shapes_dim_type const& shapes() const { return std::get<1>(*this); }

    static shape_dim_type shape( uint16_type ni, uint16_type nj ) { return std::make_pair(ni,nj); }

    void add( shape_dim_type const& shape )
    {
        auto & shapesContainer = std::get<1>(*this);
        if ( std::find( std::begin( shapesContainer ), std::end( shapesContainer ), shape ) == std::end( shapesContainer ) )
            shapesContainer.push_back( shape );
    }
    template <typename ContainerShapesType>
        void add( ContainerShapesType const& shapes, typename std::enable_if_t<is_iterable_v<ContainerShapesType> >* = nullptr )
    {
        for ( shape_dim_type const & shape : shapes )
            this->add( shape );
    }
};

template <uint16_type Dim>
class ModelPhysics
{
public :
    using material_property_description_type = MaterialPropertyDescription;
    using material_property_shape_dim_type = typename material_property_description_type::shape_dim_type;
    static const uint16_type nDim = Dim;

    ModelPhysics() = default;
    explicit ModelPhysics( std::string const& physic );
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;

    //! return name of physic assiciated to system of pde
    std::string const& physic() const { return M_physic; }
    //! return name of all physic present in system of pde
    std::set<std::string> const& physics() const { return M_physics; }
    //! map between physic to subphysics
    std::map<std::string,std::set<std::string>> const& mapPhysicsToSubphysics() const { return M_mapPhysicsToSubphysics; }
    //! return the material properties description (.i.e. coefficients of pdes)
    std::map<std::string,material_property_description_type> const& materialPropertyDescription() const { return M_materialPropertyDescription; }

protected :
    void addMaterialPropertyDescription( std::string const& propName, std::string const& symbol, std::initializer_list<material_property_shape_dim_type> const& shapes )
        {
            auto itFindProp = M_materialPropertyDescription.find( propName );
            if ( itFindProp == M_materialPropertyDescription.end() )
                M_materialPropertyDescription.emplace( propName, material_property_description_type( symbol, shapes ) );
            else
            {
                CHECK( itFindProp->second.symbol() == symbol ) << "invalid symbol in property " << propName << " : symbol given is " << symbol << "but should be " << itFindProp->second.symbol() ;
                itFindProp->second.add( shapes );
            }
        }

    //private :
    std::string M_physic; // the name of the physic
    std::set<std::string> M_physics; // all physics (example thermo-electric include also heat and electric)
    std::map<std::string,std::set<std::string>> M_mapPhysicsToSubphysics;

    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )

};

} // namespace FeelModels
} // namespace Feel

#endif
