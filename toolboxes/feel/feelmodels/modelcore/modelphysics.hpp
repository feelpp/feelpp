/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

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

#include <feel/feelmodels/modelmodels.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>

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
class ModelPhysics;

template <uint16_type Dim>
class ModelPhysic
{
public :
    using material_property_description_type = MaterialPropertyDescription;
    using material_property_shape_dim_type = typename material_property_description_type::shape_dim_type;
    inline static const uint16_type nDim = Dim;

    //ModelPhysic() = default;
    ModelPhysic( std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysic( ModelPhysic const& ) = default;
    ModelPhysic( ModelPhysic && ) = default;

    static std::shared_ptr<ModelPhysic<Dim>> New( ModelPhysics<Dim> const& mphysics, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );

    std::string const& type() const { return M_type; }
    std::string const& name() const { return M_name; }

    //! return the types of subphysics
    std::set<std::string> const& subphysicsTypes() const { return M_subphysicsTypes; }

    //! return the map of subphysics
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> const& subphysics() const { return M_subphysics; }

    //! return the material properties description (.i.e. coefficients of pdes)
    std::map<std::string,material_property_description_type> const& materialPropertyDescription() const { return M_materialPropertyDescription; }


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

    void physicsShared( std::set<std::string> & res ) const
        {
            res.insert( M_name );
            for ( auto const& [subName,subPhysic] : M_subphysics )
                subPhysic->physicsShared( res );
        }

    //! add a subphysic model
    void addSubphysic( std::shared_ptr<ModelPhysic<nDim>> const& sp ) { M_subphysics[sp->name()] = sp; }

    //! set parameter values in expression
    virtual void setParameterValues( std::map<std::string,double> const& mp ) {}
private :
    std::string M_type, M_name;
    std::set<std::string> M_subphysicsTypes;
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> M_subphysics;
    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )
};

template <uint16_type Dim>
class ModelPhysicFluid : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
public :
    ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicFluid( ModelPhysicFluid const& ) = default;
    ModelPhysicFluid( ModelPhysicFluid && ) = default;

    std::string const& equation() const { return M_equation; }
    void setEquation( std::string const& eq );

    bool gravityForceEnabled() const { return M_gravityForceEnabled; }
    auto const& gravityForceExpr() { return M_gravityForceExpr.template expr<Dim,1>(); }

    //! set parameter values in expression
    void setParameterValues( std::map<std::string,double> const& mp ) override
        {
            if ( M_gravityForceEnabled )
                M_gravityForceExpr.setParameterValues( mp );
        }

private :
    std::string M_equation;
    bool M_gravityForceEnabled;
    ModelExpression M_gravityForceExpr;
};

template <uint16_type Dim>
class ModelPhysicSolid : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
public :
    ModelPhysicSolid( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicSolid( ModelPhysicSolid const& ) = default;
    ModelPhysicSolid( ModelPhysicSolid && ) = default;

    std::string const& equation() const { return M_equation; }
    void setEquation( std::string const& eq );

    std::string const& materialModel() const { return M_materialModel; }
    std::string const& formulation() const { return M_formulation; }
    std::string const& decouplingEnergyVolumicLaw() const { return M_decouplingEnergyVolumicLaw; }
    std::string const& compressibleNeoHookeanVariantName() const { return M_compressibleNeoHookeanVariantName; }

    bool useDisplacementPressureFormulation() const { return M_formulation == "displacement-pressure"; }
private :
    std::string M_equation, M_materialModel, M_formulation;
    std::string M_decouplingEnergyVolumicLaw, M_compressibleNeoHookeanVariantName;
};

template <uint16_type Dim>
class ModelPhysics : virtual public ModelBase
{
public :
    using material_property_description_type = MaterialPropertyDescription;
    using material_property_shape_dim_type = typename material_property_description_type::shape_dim_type;
    inline static const uint16_type nDim = Dim;

    using subphysic_description_type = std::map<std::string,std::tuple<std::string,std::shared_ptr<ModelPhysics<Dim>>>>;

    //ModelPhysics() = default;
    explicit ModelPhysics( std::string const& type ) : ModelBase(""), M_physicType( type ) {}
    ModelPhysics( std::string const& type, ModelBase const& mbase ) : ModelBase(mbase), M_physicType( type ) {}
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;

    //! return all physics registerd
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> const& physics() const { return M_physics; }

    //! return all physics registerd related to a type of physic
    std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>> physics( std::string const& type ) const;

    //! return all physics registerd related to the current type
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> physicsFromCurrentType() const { return this->physics( this->physicType() ); }

    //! return the type of physic at top level
    std::string const& physicType() const { return M_physicType; }

    //! return the name of the physic used by default
    std::string const& physicDefault() const { return M_physicDefault; }

    //! return name of all physics registered
    std::set<std::string> physicsAvailable() const;

    //! return name of all physics registered related to a type of physic
    std::set<std::string> physicsAvailable( std::string const& type ) const;

    //! return name of all physics registered related to the current type
    std::set<std::string> physicsAvailableFromCurrentType() const { return this->physicsAvailable( this->physicType() ); }

    //! return the name of all physic shared from the parent physic \pname
    std::set<std::string> physicsShared( std::string const& pname ) const;


    void initPhysics( std::string const& name, ModelModels const& models, subphysic_description_type const& subPhyicsDesc = subphysic_description_type{} );

    void setPhysics( std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics, std::string const& physicDefault = "" );

protected :

    std::string M_physicType;
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> M_physics;
    std::string M_physicDefault;
};

} // namespace FeelModels
} // namespace Feel

#endif
