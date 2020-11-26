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

    //! return a subphysic from the type of physic (the first one found, and return null shared ptr if not found)
    std::shared_ptr<ModelPhysic<nDim>> subphysicFromType( std::string const& type ) const
        {
            for ( auto const& [name,sp] : this->subphysics() )
                if ( sp->type() == type )
                    return sp;
            return  std::shared_ptr<ModelPhysic<nDim>>{};
        }

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

    //! add a constant (double) parameter called \pname and return the symbol associated
    std::string addParameter( std::string const& pname, double val )
        {
            M_parameterNameToExpr.emplace( pname, ModelExpression(val) );
            return this->symbolFromParameter( pname );
        }
    //! add a parameter called \pname describe by an expression \expr and return the symbol associated
    std::string addParameter( std::string const& pname, std::string const& expr, WorldComm const& worldComm, std::string const& directoryLibExpr )
        {
            M_parameterNameToExpr[pname].setExpr( expr, worldComm, directoryLibExpr );
            return this->symbolFromParameter( pname );
        }

    //! set parameter values in expression
    virtual void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & [ pname, mexpr ] : M_parameterNameToExpr )
                mexpr.setParameterValues( mp );

            // TODO : maybe subphysics??
        }

    std::map<std::string,double>
    toParameterValues() const
        {
            std::map<std::string,double> pv;
            for ( auto const& [pname, mexpr] : M_parameterNameToExpr )
                mexpr.updateParameterValues( this->symbolFromParameter( pname ), pv );
            return pv;
        }

    void updateParameterValues( std::map<std::string,double> & mp )
        {
            this->setParameterValues( mp );

            std::set<std::string> allSymbNames;
            for ( auto const& [pname, mexpr] : M_parameterNameToExpr )
                allSymbNames.insert( pname );

            // erase parameters values from material properties
            for ( auto & [pname, mexpr] : M_parameterNameToExpr )
                mexpr.eraseParameterValues( allSymbNames );

            // get value of symbols evaluables
            int previousParam = 0;
            std::map<std::string,double> curmp;
            while ( true )
            {
                curmp = this->toParameterValues();
                if ( curmp.size() == previousParam )
                    break;
                previousParam = curmp.size();
                this->setParameterValues( curmp );
            }

            // update output parameter values
            for ( auto const& [_name,_value] : curmp )
                mp[_name] = _value;
        }

        auto symbolsExpr() const
        {
            auto tupleSymbolExprs = hana::transform( ModelExpression::expr_shapes, [this](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                    using _expr_type = std::decay_t< decltype( ModelExpression{}.template expr<ni,nj>() ) >;
                    symbol_expression_t<_expr_type> se;
                    for ( auto const& [pname, mexpr] : M_parameterNameToExpr )
                    {
                        if ( !mexpr.template hasExpr<ni,nj>() )
                            continue;
                        auto const& paramExpr = mexpr.template expr<ni,nj>();
                        if ( paramExpr.expression().isEvaluable() )
                            continue;

                        std::string symbolParam = this->symbolFromParameter( pname );
                        se.add( symbolParam, paramExpr, SymbolExprComponentSuffix( ni,nj ) );
                    }
                    return se;
                });
            return Feel::vf::SymbolsExpr( tupleSymbolExprs );
        }

    //! return the full symbol from parameter name
    std::string symbolFromParameter( std::string const& pname ) const { return "physics_" + prefixvm( prefixvm( M_type, M_name, "_" ), pname, "_" ); }

private :
    std::string M_type, M_name;
    std::set<std::string> M_subphysicsTypes;
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> M_subphysics;
    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )
    std::map<std::string,ModelExpression> M_parameterNameToExpr;
};

template <uint16_type Dim>
class ModelPhysicFluid : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
public :
    struct DynamicViscosity
    {
        DynamicViscosity( std::string const& ln )
            :
            M_lawName( ln )
            {
                CHECK( this->checkLawName() ) << "invalid law name " << ln;
            }
        DynamicViscosity( DynamicViscosity const& ) = default;
        DynamicViscosity( DynamicViscosity && ) = default;

        void setLaw( std::string const& ln ) { M_lawName = ln; CHECK( this->checkLawName() ) << "invalid law name " << ln; }
        std::string const& lawName() const { return M_lawName; }

        bool isNewtonianLaw() const { return (this->lawName() == "newtonian"); }
        bool isPowerLaw() const { return (this->lawName() == "power_law"); }
        bool isCarreauLaw() const { return (this->lawName() == "carreau_law"); }
        bool isCarreauYasudaLaw() const { return (this->lawName() == "carreau-yasuda_law"); }

        bool checkLawName() const
            {
                return ( this->isNewtonianLaw() || this->isPowerLaw() || this->isCarreauLaw() || this->isCarreauYasudaLaw() );
            }

    private :
        std::string M_lawName;
    };

    struct Turbulence
    {
        Turbulence() : M_isEnabled( false ) {}
        Turbulence( Turbulence const& ) = default;
        Turbulence( Turbulence && ) = default;

        bool isEnabled() const { return M_isEnabled; }
        std::string const& model() const { return M_model; }

        bool useBoussinesqApproximation() const;
        bool hasTurbulentKineticEnergy() const;

        void setup( pt::ptree const& p );

    private :
        bool M_isEnabled;
        std::string M_model;
    };

    ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicFluid( ModelPhysicFluid const& ) = default;
    ModelPhysicFluid( ModelPhysicFluid && ) = default;

    std::string const& equation() const { return M_equation; }
    void setEquation( std::string const& eq );

    bool gravityForceEnabled() const { return M_gravityForceEnabled; }
    auto const& gravityForceExpr() const { return M_gravityForceExpr.template expr<Dim,1>(); }

    DynamicViscosity const& dynamicViscosity() const { return M_dynamicViscosity; }
    Turbulence const& turbulence() const { return M_turbulence; }

    //! set parameter values in expression
    void setParameterValues( std::map<std::string,double> const& mp ) override
        {
            super_type::setParameterValues( mp );
            if ( M_gravityForceEnabled )
                M_gravityForceExpr.setParameterValues( mp );
        }

private :
    std::string M_equation;
    bool M_gravityForceEnabled;
    ModelExpression M_gravityForceExpr;

    DynamicViscosity M_dynamicViscosity;
    Turbulence M_turbulence;
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

    using model_physic_type = ModelPhysic<nDim>;
    using model_physic_ptrtype = std::shared_ptr<model_physic_type>;
    using subphysic_description_type = std::map<std::string,std::tuple<std::string,std::shared_ptr<ModelPhysics<Dim>>>>;

    //ModelPhysics() = default;
    explicit ModelPhysics( std::string const& type ) : ModelBase(""), M_physicType( type ) {}
    ModelPhysics( std::string const& type, ModelBase const& mbase ) : ModelBase(mbase), M_physicType( type ) {}
    ModelPhysics( std::string const& type, ModelBase && mbase ) : ModelBase(std::move(mbase)), M_physicType( type ) {}
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;
    virtual ~ModelPhysics() = default;

    //void updateInformationObject( pt::ptree & p ) const override {}

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

    auto symbolsExprPhysics( std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> const& mphysics ) const
        {
            using se_type = std::decay_t<decltype(mphysics.begin()->second->symbolsExpr())>;
            se_type res;
            for ( auto const& [name, mphysic] : mphysics )
                res = Feel::vf::symbolsExpr( res, mphysic->symbolsExpr() );
            return res;
        }
    auto symbolsExprPhysicsFromCurrentType() const
        {
            return symbolsExprPhysics( this->physicsFromCurrentType() );
        }
protected :

    std::string M_physicType;
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> M_physics;
    std::string M_physicDefault;
};

} // namespace FeelModels
} // namespace Feel

#endif
