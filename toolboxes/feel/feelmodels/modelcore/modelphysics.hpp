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

/**
 * @brief Description of Model Properties
 * @ingroup ModelCore
 * 
 */
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

/**
 * @brief Physic of the Model
 * @ingroup ModelCore
 * 
 * @tparam Dim real dimension
 */
template <uint16_type Dim>
class ModelPhysic : public JournalWatcher
{
    using super_type = JournalWatcher;
public :
    using material_property_description_type = MaterialPropertyDescription;
    using material_property_shape_dim_type = typename material_property_description_type::shape_dim_type;
    static inline const uint16_type nDim = Dim;
    using physic_id_type = std::pair<std::string,std::string>;

    //ModelPhysic() = default;
    ModelPhysic( std::string const& modeling, std::string const& type, std::string const& name, ModelBase const& mparent, ModelModel const& model = ModelModel{} );
    ModelPhysic( ModelPhysic const& ) = default;
    ModelPhysic( ModelPhysic && ) = default;

    static std::shared_ptr<ModelPhysic<Dim>> New( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );

    void updateInformationObject( nl::json & p ) const override;
    virtual tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const;

    std::string const& modeling() const { return M_modeling; }
    std::string const& type() const { return M_type; }
    std::string const& name() const { return M_name; }

    physic_id_type id() const { return std::make_pair(this->type(),this->name()); }

    //! return the map of subphysics
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> const& subphysics() const { return M_subphysics; }

    //! return a subphysic from the type of physic (the first one found, and return null shared ptr if not found)
    std::shared_ptr<ModelPhysic<nDim>> subphysicFromType( std::string const& type ) const
        {
            for ( auto const& [name,sp] : this->subphysics() )
                if ( sp->type() == type )
                    return sp;
            return  std::shared_ptr<ModelPhysic<nDim>>{};
        }
    //! return all subphysic from the type of physic
    std::vector<std::shared_ptr<ModelPhysic<nDim>>> subphysicsFromType( std::string const& type ) const
        {
            std::vector<std::shared_ptr<ModelPhysic<nDim>>> res;
            for ( auto const& [name,sp] : this->subphysics() )
                if ( sp->type() == type )
                    res.push_back( sp );
            return res;
        }

    //! return material names used by this physic
    std::set<std::string> const& materialNames() const { return M_materialNames; }

    //! set material names
    void setMaterialNames( std::set<std::string> const& matNames ) { M_materialNames = matNames; }

    //! add material names
    void addMaterialNames( std::set<std::string> const& matNames ) { M_materialNames.insert( matNames.begin(), matNames.end() ); }

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


    bool hasSubphysic( physic_id_type const& pId ) const { return M_subphysics.find( pId ) != M_subphysics.end(); }

    //! add a subphysic model
    void addSubphysic( std::shared_ptr<ModelPhysic<nDim>> const& sp, bool force = false )
        {
            if ( force )
                M_subphysics[std::make_pair(sp->type(),sp->name())] = sp;
            else
                M_subphysics.emplace( std::make_pair(sp->type(),sp->name()), sp );
        }

    //! add a constant (double) parameter called \pname and return the symbol associated
    std::string addParameter( std::string const& pname, double val )
        {
            //M_parameterNameToExpr.emplace( pname, ModelExpression(val) );
            M_parameterNameToExpr[pname].reset();
            M_parameterNameToExpr[pname].setExpr( val );
            return this->symbolFromParameter( pname );
        }
    //! add a parameter called \pname describe by an expression \expr and return the symbol associated
    std::string addParameter( std::string const& pname, std::string const& expr, WorldComm const& worldComm, std::string const& directoryLibExpr )
        {
            M_parameterNameToExpr[pname].reset();
            M_parameterNameToExpr[pname].setExpr( expr, worldComm, directoryLibExpr );
            return this->symbolFromParameter( pname );
        }
    std::string addParameter( std::string const& pname, std::string const& expr )
        {
            return this->addParameter( pname,expr, *M_worldComm, M_directoryLibExpr );
        }
    //! add a parameter called \pname describe by a json entry \jarg and return the symbol associated
    std::string addParameter( std::string const& pname, nl::json const& jarg, WorldComm const& worldComm, std::string const& directoryLibExpr )
        {
            M_parameterNameToExpr[pname].reset();
            M_parameterNameToExpr[pname].setExpr( jarg, worldComm, directoryLibExpr );
            return this->symbolFromParameter( pname );
        }
    std::string addParameter( std::string const& pname, nl::json const& jarg )
        {
            return this->addParameter( pname, jarg, *M_worldComm, M_directoryLibExpr );
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

    //! return true if parameter \pname is present, false otherwise
    bool hasParameter( std::string const& pname ) const
        {
            return M_parameterNameToExpr.find( pname ) != M_parameterNameToExpr.end();
        }

    //! return true if parameter \pname of shape MxN is present, false otherwise
    template <int M,int N>
    bool hasParameterExpr( std::string const& pname ) const
        {
            auto itFindParam = M_parameterNameToExpr.find( pname );
            if ( itFindParam == M_parameterNameToExpr.end() )
                return false;
            return itFindParam->second.template hasExpr<M,N>();
        }

    //! return the expression of the parameter \pname
    template <int M,int N,typename SymbolsExprType = symbols_expression_empty_t>
    auto parameterExpr( std::string const& pname, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            if ( !this->hasParameterExpr<M,N>( pname ) )
                CHECK( false ) << "parameter " << pname << " not found or shape incompatible";
            auto itFindParam = M_parameterNameToExpr.find( pname );
            return expr( itFindParam->second.template expr<M,N>(), se );
        }
protected :
    ModelExpression const& parameterModelExpr( std::string const& pname ) const
        {
            auto itFindParam = M_parameterNameToExpr.find( pname );
            CHECK( itFindParam!=M_parameterNameToExpr.end() ) << "no parmeter registered : "<< pname;
            return itFindParam->second;
        }
    void updateTabulateInformationsBasic( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const;
    void updateTabulateInformationsSubphysics( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const;
    void updateTabulateInformationsParameters( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const;

private :
    template <uint16_type TheDim>
    friend class ModelPhysics;

private :
    std::string M_modeling, M_type, M_name;
    worldcomm_ptr_t M_worldComm;
    std::string M_directoryLibExpr;
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> M_subphysics;
    std::set<std::string> M_materialNames;
    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )
    std::map<std::string,ModelExpression> M_parameterNameToExpr;
};

/**
 * @brief Heat Physic Model
 * @ingroup ModelCore
 * 
 * @tparam Dim real dimension
 */
template <uint16_type Dim>
class ModelPhysicHeat : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicHeat<Dim>;
public :

    struct Convection
    {
        Convection( self_type * mparent ) : M_parent( mparent ), M_enabled( false ) {}
        Convection( Convection const& ) = default;
        Convection( Convection && ) = default;

        void setup( nl::json const& jarg );

        bool enabled() const { return M_enabled; }

        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
            {
                return M_parent->template parameterExpr<Dim,1>( "convection", se );
            }

        void updateInformationObject( nl::json & p ) const;
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
    private:
        self_type * M_parent;
        bool M_enabled;
    };
    struct HeatSource
    {
        HeatSource( self_type * mparent, std::string const& name ) : M_parent( mparent ), M_name( name ), M_type("heat-source") {}
        HeatSource( HeatSource const& ) = default;
        HeatSource( HeatSource && ) = default;

        void setup( nl::json const& jarg );

        std::string const& type() const { return M_type; }

        bool givenAsHeatRate() const { return M_type == "heat-rate"; }

        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
            {
                return M_parent->template parameterExpr<1,1>( M_name + "_heatsource", se );
            }

        void updateInformationObject( nl::json & p ) const;
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
    private:
        self_type * M_parent;
        std::string M_name;
        std::string M_type;
    };

    ModelPhysicHeat( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicHeat( ModelPhysicHeat const& ) = default;
    ModelPhysicHeat( ModelPhysicHeat && ) = default;

    std::vector<HeatSource> heatSources() const { return M_heatSources; }

    bool hasConvectionEnabled() const { return M_convection && M_convection->enabled(); }
    Convection const& convection() const { CHECK( M_convection ) << "no convection"; return *M_convection; }

    void setupConvection( nl::json & jarg );

    //! set parameter values in expression
    void setParameterValues( std::map<std::string,double> const& mp ) override
        {
            super_type::setParameterValues( mp );
        }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    template <typename MaterialsPropertiesType, typename MeshType>
    void updateForUse( std::shared_ptr<MaterialsPropertiesType> const& matProp, std::shared_ptr<MeshType> mesh )
        {
            bool requireMeasureMat = false;
            for ( auto const& hs : M_heatSources )
                if ( hs.givenAsHeatRate() )
                {
                    requireMeasureMat = true;
                    break;
                }
            if ( requireMeasureMat )
            {
                double measureMat = 0;
                for ( std::string const& matName : matProp->physicToMaterials( this->id() ) )
                {
                    auto const& rangeMat = matProp->rangeMeshElementsByMaterial( mesh,matName );
                    measureMat += measure(_range=rangeMat);
                }
                this->addParameter( "measure_materials", measureMat );
            }
        }
private:
    std::vector<HeatSource> M_heatSources;
    std::optional<Convection> M_convection;
};

/**
 * @brief Electric Physic Model
 * @ingroup ModelCore
 * 
 * @tparam Dim real dimension of the model
 */
template <uint16_type Dim>
class ModelPhysicElectric : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicElectric<Dim>;
public :
    ModelPhysicElectric( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicElectric( ModelPhysicElectric const& ) = default;
    ModelPhysicElectric( ModelPhysicElectric && ) = default;

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;
private :
};

/**
 * @brief ThermoElectric Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa ThermoElectric
 */
template <uint16_type Dim>
class ModelPhysicThermoElectric : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicHeat<Dim>;
public :
    ModelPhysicThermoElectric( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicThermoElectric( ModelPhysicThermoElectric const& ) = default;
    ModelPhysicThermoElectric( ModelPhysicThermoElectric && ) = default;

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;
private :
};


/**
 * @brief Fluid Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa Fluid
 */
template <uint16_type Dim>
class ModelPhysicFluid : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicFluid<Dim>;
public :

    enum class NavierStokesFormulation { Convective=0, SkewSymmetric, Conservative, Rotational, EMAC };
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
        Turbulence( self_type * mparent ) : M_parent( mparent ), M_isEnabled( false ) {}
        Turbulence( Turbulence const& ) = default;
        Turbulence( Turbulence && ) = default;

        bool isEnabled() const { return M_isEnabled; }
        std::string const& model() const { return M_model; }

        bool useBoussinesqApproximation() const;
        bool hasTurbulentKineticEnergy() const;

        void setup( nl::json const& jarg );

        void setFrictionVelocityWallFunction( std::string const& matName, std::string const& expr )
            {
                this->setParameterOnMaterial( matName, "u_tau_wallfunction", expr );
            }
        std::string frictionVelocityWallFunctionSymbol( std::string const& matName ) const
            {
                return M_parent->symbolFromParameter( this->symbolBaseNameOnMaterial( matName, "u_tau_wallfunction" ) );
            }
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto frictionVelocityWallFunctionExpr( std::string const& matName, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
            {
                return M_parent->template parameterExpr<1,1>( this->symbolBaseNameOnMaterial( matName, "u_tau_wallfunction" ), se );
            }

        std::string vonKarmanConstantSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialNoModelLink("kappa") ); }
        // upper limit on the mixing length (used by k-epsilon)
        std::string mixingLengthLimitSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialNoModelLink("l_mix_lim") ); }

        std::string kEpsilon_c_muSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialLink( "c_mu", "k-epsilon" ) ); }
        std::string kEpsilon_c_1epsilonSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialLink( "c_1epsilon", "k-epsilon" ) ); }
        std::string kEpsilon_c_2epsilonSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialLink( "c_2epsilon", "k-epsilon" ) ); }
        std::string kEpsilon_sigma_kSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialLink( "sigma_k", "k-epsilon" ) ); }
        std::string kEpsilon_sigma_epsilonSymbol() const { return M_parent->symbolFromParameter( this->symbolBaseNameNoMaterialLink( "sigma_epsilon", "k-epsilon" ) ); }


        void setParameterOnMaterial( std::string const& matName, std::string const& propName, std::string const& expr )
            {
                M_parent->addParameter( this->symbolBaseNameOnMaterial( matName,propName ), expr );
            }
        void setParameterOnMaterial( std::string const& matName, std::string const& propName, std::string const& tmodel, std::string const& expr )
            {
                M_parent->addParameter( this->symbolBaseNameOnMaterial( matName,propName,tmodel ), expr );
            }
        void setParameterNoMaterialLink( std::string const& propName, std::string const& expr, std::string const& tmodel )
            {
                M_parent->addParameter( this->symbolBaseNameNoMaterialLink( propName, tmodel ), expr );
            }
        void setParameterNoMaterialLink( std::string const& propName, std::string const& expr )
            {
                M_parent->addParameter( this->symbolBaseNameNoMaterialLink( propName ), expr );
            }
        void setParameterNoMaterialNoModelLink( std::string const& propName, std::string const& expr )
            {
                M_parent->addParameter( this->symbolBaseNameNoMaterialNoModelLink( propName ), expr );
            }

    private :
        std::string prefixSymbolTurbulenceModel() const { return this->prefixSymbolTurbulenceModel( M_model ); }

        std::string prefixSymbolTurbulenceModel( std::string const& tmodel ) const
            {
                std::string prefix = "turbulence";
                if ( tmodel.empty() )
                    return prefix;
                else if ( tmodel == "k-epsilon" )
                    return prefixvm( prefix, "k_epsilon", "_" );
                return prefix;
            }
        std::string symbolBaseNameOnMaterial( std::string const& matName, std::string const& propName, std::string const& tmodel ) const
            {
                std::string prefix = this->prefixSymbolTurbulenceModel( tmodel );
                if ( !matName.empty() )
                    prefix += "_" + matName;
                return prefixvm( prefix,propName,"_");
            }
        std::string symbolBaseNameOnMaterial( std::string const& matName, std::string const& propName ) const
            {
                return this->symbolBaseNameOnMaterial( matName, propName, M_model );
            }
        std::string symbolBaseNameNoMaterialLink( std::string const& propName ) const
            {
                return this->symbolBaseNameNoMaterialLink( propName, M_model );
            }
        std::string symbolBaseNameNoMaterialLink( std::string const& propName, std::string const& tmodel ) const
            {
                return this->symbolBaseNameOnMaterial( "", propName, tmodel );
            }
        std::string symbolBaseNameNoMaterialNoModelLink( std::string const& propName ) const
            {
                return this->symbolBaseNameOnMaterial( "", propName, "" );
            }

    private :
        self_type * M_parent;
        bool M_isEnabled;
        std::string M_model;
    };

    ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicFluid( ModelPhysicFluid const& ) = default;
    ModelPhysicFluid( ModelPhysicFluid && ) = default;

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    std::string const& equation() const { return M_equation; }
    void setEquation( std::string const& eq );

    NavierStokesFormulation navierStokesFormulation() const { return M_navierStokesFormulation; }

    bool gravityForceEnabled() const { return M_gravityForceEnabled; }
    auto const& gravityForceExpr() const { return M_gravityForceExpr.template expr<Dim,1>(); }

    DynamicViscosity const& dynamicViscosity() const { return M_dynamicViscosity; }
    Turbulence const& turbulence() const { return M_turbulence; }
    Turbulence & turbulence() { return M_turbulence; }

    //! set parameter values in expression
    void setParameterValues( std::map<std::string,double> const& mp ) override
        {
            super_type::setParameterValues( mp );
            if ( M_gravityForceEnabled )
                M_gravityForceExpr.setParameterValues( mp );
        }

private :
    std::string M_equation;
    NavierStokesFormulation M_navierStokesFormulation;

    bool M_gravityForceEnabled;
    ModelExpression M_gravityForceExpr;

    DynamicViscosity M_dynamicViscosity;
    Turbulence M_turbulence;
};

/**
 * @brief Solid Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa Solid
 */
template <uint16_type Dim>
class ModelPhysicSolid : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicSolid<Dim>;
public :
    struct BodyForces
    {
        BodyForces( self_type * mparent, std::string const& name ) : M_parent( mparent ), M_name( name ) {}
        BodyForces( BodyForces const& ) = default;
        BodyForces( BodyForces && ) = default;

        void setup( nl::json const& jarg );

        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
            {
                return M_parent->template parameterExpr<Dim,1>( M_name + "_bodyforce", se );
            }

        void updateInformationObject( nl::json & p ) const;
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
    private:
        self_type * M_parent;
        std::string M_name;
    };

    ModelPhysicSolid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicSolid( ModelPhysicSolid const& ) = default;
    ModelPhysicSolid( ModelPhysicSolid && ) = default;

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    std::string const& equation() const { return M_equation; }
    void setEquation( std::string const& eq );

    std::string const& materialModel() const { return M_materialModel; }
    std::string const& formulation() const { return M_formulation; }
    std::string const& decouplingEnergyVolumicLaw() const { return M_decouplingEnergyVolumicLaw; }
    std::string const& compressibleNeoHookeanVariantName() const { return M_compressibleNeoHookeanVariantName; }

    bool useDisplacementPressureFormulation() const { return M_formulation == "displacement-pressure"; }

    std::vector<BodyForces> bodyForces() const { return M_bodyForces; }
private :
    std::string M_equation, M_materialModel, M_formulation;
    std::string M_decouplingEnergyVolumicLaw, M_compressibleNeoHookeanVariantName;
    std::vector<BodyForces> M_bodyForces;
};

/**
 * @brief Fluid Structure Interaction Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa FSI
 */
template <uint16_type Dim>
class ModelPhysicFSI : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicFSI<Dim>;
public :
    ModelPhysicFSI( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicFSI( ModelPhysicFSI const& ) = default;
    ModelPhysicFSI( ModelPhysicFSI && ) = default;

    std::set<std::string> const& interfaceFluid() const { return M_interfaceFluid; }
    std::set<std::string> const& interfaceSolid() const { return M_interfaceSolid; }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;
private :
    std::set<std::string> M_interfaceFluid, M_interfaceSolid;
};

/**
 * @brief CoefficientFormPDE Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa CoefficientFormPDE
 */
template <uint16_type Dim>
class ModelPhysicCoefficientFormPDE : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicCoefficientFormPDE<Dim>;
public :
    struct Infos
    {
        enum class Coefficient { convection=0, diffusion, reaction, firstTimeDerivative, secondTimeDerivative, source, conservativeFluxConvection, conservativeFluxSource, curlCurl };

        using shape_dim_type = std::pair<uint16_type,uint16_type>;
        using shapes_dim_type = std::vector<shape_dim_type>;

        Infos() = default;
        Infos( std::string const& name, nl::json const& jarg );
        Infos( std::string const& name, std::string const& unknownName, std::string const& unknownSymbol, std::string const& unknownBasis )
            :
            M_equationName( name ),
            M_unknownName( unknownName ),
            M_unknownSymbol( unknownSymbol ),
            M_unknownBasis( unknownBasis )
            {}
        Infos( Infos const& ) = default;
        Infos( Infos && ) = default;
        std::string const& equationName() const { return M_equationName; }
        std::string const& unknownName() const { return M_unknownName; }
        std::string const& unknownSymbol() const { return M_unknownSymbol; }
        std::string const& unknownBasis() const { return M_unknownBasis; }

        std::map<Coefficient,std::tuple<std::string,shapes_dim_type>> const& coefficientProperties() const { return M_coefficientProperties; }
        std::string const& coefficientName( Coefficient c ) const { return std::get<0>( M_coefficientProperties.at( c ) ); }
    private :
        std::string M_equationName, M_unknownName, M_unknownSymbol, M_unknownBasis;
        std::map<Coefficient,std::tuple<std::string,shapes_dim_type>> M_coefficientProperties;
    };
    using infos_type = self_type::Infos;
    using infos_ptrtype = std::shared_ptr<infos_type>;

    ModelPhysicCoefficientFormPDE( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicCoefficientFormPDE( ModelPhysicCoefficientFormPDE const& ) = default;
    ModelPhysicCoefficientFormPDE( ModelPhysicCoefficientFormPDE && ) = default;

    std::shared_ptr<infos_type> const& infos() const { return M_infos; }

    using Coefficient = typename infos_type::Coefficient;
    bool hasCoefficient( Coefficient c ) const { return this->hasCoefficient( M_infos->coefficientName( c ) ); }
    ModelExpression const& coefficient( Coefficient c ) const { return this->coefficient( M_infos->coefficientName( c ) ); }

    template <int M,int N,typename SymbolsExprType = symbols_expression_empty_t>
    auto coefficientExpr( Coefficient c, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return this->coefficientExpr<M,N>( M_infos->coefficientName( c ), se );
        }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;
private:
    bool hasCoefficient( std::string const& coeffName ) const { return this->hasParameter( coeffName ); }
    ModelExpression const& coefficient( std::string const& coeffName ) const { return this->parameterModelExpr( coeffName ); }
    template <int M,int N,typename SymbolsExprType = symbols_expression_empty_t>
    auto coefficientExpr( std::string const& coeffName, SymbolsExprType const& se = symbols_expression_empty_t{} ) const { return this->template parameterExpr<M,N>( coeffName, se ); }
private :
    std::shared_ptr<infos_type> M_infos;
};

/**
 * @brief CoefficientFormPDE set Physic Model
 * @ingroup ModelCore
 *
 * @tparam Dim real dimension of the model
 * @sa CoefficientFormPDE
 */
template <uint16_type Dim>
class ModelPhysicCoefficientFormPDEs : public ModelPhysic<Dim>
{
    using super_type = ModelPhysic<Dim>;
    using self_type = ModelPhysicCoefficientFormPDEs<Dim>;
public :
    ModelPhysicCoefficientFormPDEs( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model = ModelModel{} );
    ModelPhysicCoefficientFormPDEs( ModelPhysicCoefficientFormPDEs const& ) = default;
    ModelPhysicCoefficientFormPDEs( ModelPhysicCoefficientFormPDEs && ) = default;

    auto const& pdes() const { return M_pdes; }
    auto & pdes() { return M_pdes; }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

private :
    //std::vector<std::tuple<typename ModelPhysicCoefficientFormPDE<Dim>::infos_type,std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>>> M_pdes;
    using cfpde_infos_type = typename ModelPhysicCoefficientFormPDE<Dim>::infos_type;
    std::vector<std::tuple<std::string,std::shared_ptr<cfpde_infos_type>,std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>>> M_pdes;
};

/**
 * @brief Physic set Model
 * @ingroup ModelCore
 * 
 * @tparam Dim 
 */
template <uint16_type Dim>
class ModelPhysics : virtual public ModelBase
{
public :
    using material_property_description_type = MaterialPropertyDescription;
    using material_property_shape_dim_type = typename material_property_description_type::shape_dim_type;
    static inline const uint16_type nDim = Dim;
    using physic_id_type = std::pair<std::string,std::string>; // (type, name)

    using model_physic_type = ModelPhysic<nDim>;
    using model_physic_ptrtype = std::shared_ptr<model_physic_type>;

protected :

    struct PhysicsTreeNode;
    struct PhysicsTree {
        using physics_ptrtype = std::shared_ptr<ModelPhysics<nDim>>;

        PhysicsTree( physics_ptrtype p )
            :
            M_root( p ) {}

        std::shared_ptr<PhysicsTreeNode> addNode( physics_ptrtype mphysics, model_physic_ptrtype mp )
            {
                auto treeNode = PhysicsTreeNode::New( mphysics,mp );
                M_children.push_back( std::move(treeNode) );
                return M_children.back();
            }

        std::vector<std::shared_ptr<PhysicsTreeNode>> const& children() const { return M_children; }

        void updatePhysics( physics_ptrtype mphysics, ModelModels const& models );

        std::map<physic_id_type,model_physic_ptrtype> collectAllPhysics() const;

    private :
        physics_ptrtype M_root;
        std::vector<std::shared_ptr<PhysicsTreeNode>> M_children;
    };
    struct PhysicsTreeNode {
        using physics_ptrtype = std::shared_ptr<ModelPhysics<nDim>>;

        static std::shared_ptr<PhysicsTreeNode> New( physics_ptrtype p, model_physic_ptrtype mp ) { return std::make_shared<PhysicsTreeNode>(p,mp); }

        PhysicsTreeNode( std::string const& type, physics_ptrtype p, model_physic_ptrtype mp ) : M_data( std::make_tuple(type,p,mp) ) {}
        PhysicsTreeNode( physics_ptrtype p, model_physic_ptrtype mp ) : PhysicsTreeNode( p->keyword(), p, mp ) {}

        void addChild( std::string const& type, physics_ptrtype p, ModelModels const& models );
        void addChild( physics_ptrtype p, ModelModels const& models ) { this->addChild( p->keyword(), p, models ); }

        physics_ptrtype modelPhysics() const { return std::get<1>( M_data ); }
        model_physic_ptrtype physic() const { return std::get<2>( M_data ); }

        std::vector<std::shared_ptr<PhysicsTreeNode>> const& children() const { return M_children; }

        void updateMaterialSupportFromChildren( std::string const& applyType = "intersect" );
    private :
        std::tuple<std::string,physics_ptrtype,model_physic_ptrtype> M_data;
        std::vector<std::shared_ptr<PhysicsTreeNode>> M_children;
        PhysicsTree* M_globalTree;
    };

public:

    //ModelPhysics() = default;
    explicit ModelPhysics( std::string const& modeling ) : ModelBase(""), M_physicModeling( modeling ) {}
    ModelPhysics( std::string const& modeling, ModelBase const& mbase ) : ModelBase(mbase), M_physicModeling( modeling ) {}
    ModelPhysics( std::string const& modeling, ModelBase && mbase ) : ModelBase(std::move(mbase)), M_physicModeling( modeling ) {}
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;
    virtual ~ModelPhysics() = default;

    //void updateInformationObject( nl::json & p ) const override {}
    void updateInformationObjectFromCurrentType( nl::json & p ) const;

    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;


    //! return true if physic id is already registered
    bool hasPhysic( physic_id_type const& pId ) const { return M_physics.find( pId ) != M_physics.end(); }

    //! return physic registered from physic id
    std::shared_ptr<ModelPhysic<Dim>> physic( physic_id_type const& pId ) const { CHECK( this->hasPhysic(pId) ) << "no physic id"; return M_physics.find( pId )->second; }

    //! return all physics registered
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> const& physics() const { return M_physics; }

    //! return all physics registered related to a type of physic
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<Dim>>> physics( std::string const& type ) const;

    //! return all physics registered related to the current type
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> physicsFromCurrentType() const { return this->physics( this->physicType() ); }

    //! return the physic modeling
    std::string const& physicModeling() const { return M_physicModeling; }

    //! return the type of physic at top level
    std::string const& physicType() const { return this->keyword();/*M_physicType;*/ }

    //! return name of all physics registered
    std::set<physic_id_type> physicsAvailable() const;

    //! return name of all physics registered related to a type of physic
    std::set<physic_id_type> physicsAvailable( std::string const& type ) const;

    //! return name of all physics registered related to the current type
    std::set<physic_id_type> physicsAvailableFromCurrentType() const { return this->physicsAvailable( this->physicType() ); }

    // TODO try remove first arg, because it just the shared_ptr of the current object
    void initPhysics( std::shared_ptr<ModelPhysics<nDim>> mphysics, ModelModels const& models );
    void initPhysics( std::shared_ptr<ModelPhysics<nDim>> mphysics, std::function<void(PhysicsTree &)> lambdaInit );


    //! add physics if not already registered
    void addPhysics( std::map<physic_id_type,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics );

    auto symbolsExprPhysics( std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> const& mphysics ) const
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
    std::vector<std::shared_ptr<ModelPhysic<Dim>>>
    initPhysics( std::string const& type, ModelModels const& models, std::set<std::string> const& modelNameRequired, bool isRequired = true );

    virtual void updatePhysics( PhysicsTreeNode & physicsTree, ModelModels const& models ) {}

private :
    template <uint16_type TheDim>
    friend class ModelPhysic;

protected :

    std::string M_physicModeling/*, M_physicType*/;
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> M_physics;
};

} // namespace FeelModels
} // namespace Feel

#endif
