/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelPhysic<Dim>::ModelPhysic( std::string const& type, std::string const& name, ModelBase const& mparent, ModelModel const& model )
    :
    M_type( type ),
    M_name( name ),
    M_worldComm( mparent.worldCommPtr() ),
    M_directoryLibExpr( mparent.repository().expr() )
{

    // save name of submodel but not yet intialized
    for ( std::string const& submodel : model.submodels() )
        M_subphysics[submodel];

    if ( type == "GenericPDE" || type == "GenericPDEs" )
        return;

    if ( M_type == "thermo-electric" )
        M_subphysicsTypes.insert( { "heat","electric" } );
    else if ( M_type == "heat-fluid" )
        M_subphysicsTypes.insert( { "heat","fluid" } );


    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    this->addMaterialPropertyDescription( "density", "rho", { scalarShape } );
    if ( M_type == "heat" ||  M_type == "thermo-electric" || M_type == "heat-fluid" )
    {
        this->addMaterialPropertyDescription( "specific-heat-capacity", "Cp", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-expansion", "beta", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-conductivity", "k", { scalarShape,matrixShape } );
    }
    if ( M_type == "electric" || M_type == "thermo-electric" )
    {
        this->addMaterialPropertyDescription( "electric-conductivity", "sigma", { scalarShape } );
    }
    if ( M_type == "solid" )
    {
        this->addMaterialPropertyDescription( "Young-modulus", "E", { scalarShape } );
        this->addMaterialPropertyDescription( "Poisson-ratio", "nu", { scalarShape } );
        this->addMaterialPropertyDescription( "Lame-first-parameter", "lambda", { scalarShape } );
        this->addMaterialPropertyDescription( "Lame-second-parameter", "mu", { scalarShape } );
        this->addMaterialPropertyDescription( "bulk-modulus", "K", { scalarShape } );
    }
    if ( M_type == "fluid" || M_type == "heat-fluid" )
    {
        this->addMaterialPropertyDescription( "dynamic-viscosity", "mu", { scalarShape } );
        this->addMaterialPropertyDescription( "turbulent-dynamic-viscosity", "mu_t", { scalarShape } );
        this->addMaterialPropertyDescription( "turbulent-kinetic-energy", "tke", { scalarShape } );
        this->addMaterialPropertyDescription( "consistency-index", "mu_k", { scalarShape } );
        this->addMaterialPropertyDescription( "power-law-index", "mu_power_law_n", { scalarShape } );
        this->addMaterialPropertyDescription( "viscosity-min", "mu_min", { scalarShape } );
        this->addMaterialPropertyDescription( "viscosity-max", "mu_max", { scalarShape } );
        this->addMaterialPropertyDescription( "viscosity-zero-shear", "mu_0", { scalarShape } );
        this->addMaterialPropertyDescription( "viscosity-infinite-shear", "mu_inf", { scalarShape } );
        this->addMaterialPropertyDescription( "carreau-law-lambda", "mu_carreau_law_lambda", { scalarShape } );
        this->addMaterialPropertyDescription( "carreau-law-n", "mu_carreau_law_n", { scalarShape } );
        this->addMaterialPropertyDescription( "carreau-yasuda-law-lambda", "mu_carreau_yasuda_law_lambda", { scalarShape } );
        this->addMaterialPropertyDescription( "carreau-yasuda-law-n", "mu_carreau_yasuda_law_n", { scalarShape } );
        this->addMaterialPropertyDescription( "carreau-yasuda-law-a", "mu_carreau_yasuda_law_a", { scalarShape } );
    }

}

template <uint16_type Dim>
std::shared_ptr<ModelPhysic<Dim>>
ModelPhysic<Dim>::New( ModelPhysics<Dim> const& mphysics, std::string const& type, std::string const& name, ModelModel const& model )
{
    if ( type == "heat" )
        return std::make_shared<ModelPhysicHeat<Dim>>( mphysics, name, model );
    if ( type == "fluid" )
        return std::make_shared<ModelPhysicFluid<Dim>>( mphysics, name, model );
    else if ( type == "solid" )
        return std::make_shared<ModelPhysicSolid<Dim>>( mphysics, name, model );
    else
        return std::make_shared<ModelPhysic<Dim>>( type, name, mphysics, model );
}

template <uint16_type Dim>
ModelPhysicHeat<Dim>::ModelPhysicHeat( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model )
    :
    super_type( "heat", name, model )
{

}

template <uint16_type Dim>
ModelPhysicFluid<Dim>::ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model )
    :
    super_type( "fluid", name, mphysics, model ),
    M_equation( "Navier-Stokes" ),
    M_gravityForceEnabled( boption(_name="use-gravity-force",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) ),
    M_dynamicViscosity( soption(_name="viscosity.law",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) ),
    M_turbulence( this )
{

    auto const& j_model = model.jsonProperties();
    for ( std::string const& eqarg : { "equations", "equation" } )
        if ( j_model.contains( eqarg ) )
        {
            auto const& j_model_eq = j_model.at( eqarg );
            if ( j_model_eq.is_string() )
                this->setEquation( j_model_eq.get<std::string>() );
        }


    // setup gravity force
    if ( j_model.contains("gravity") )
    {
        auto const& j_model_gravity = j_model.at("gravity");
        if ( j_model_gravity.is_boolean() )
            M_gravityForceEnabled = j_model_gravity.template get<bool>();
        else if ( j_model_gravity.is_string() )
            M_gravityForceEnabled = boost::lexical_cast<bool>( j_model_gravity.template get<std::string>() );
        else if ( j_model_gravity.is_object() )
        {
            if ( j_model_gravity.contains( "enable" ) )
            {
                auto const& j_model_gravity_enable =  j_model_gravity.at("enable");
                if ( j_model_gravity_enable.is_boolean() )
                    M_gravityForceEnabled = j_model_gravity_enable.template get<bool>();
                else if ( j_model_gravity_enable.is_string() )
                    M_gravityForceEnabled = boost::lexical_cast<bool>( j_model_gravity_enable.template get<std::string>() );
            }
            if ( j_model_gravity.contains( "expr" ) )
                M_gravityForceExpr.setExpr( j_model_gravity.at("expr"),mphysics.worldComm(),mphysics.repository().expr() );
        }
    }
    if ( M_gravityForceEnabled && !M_gravityForceExpr.template hasExpr<Dim,1>() )
    {
        std::string gravityStr;
        if ( mphysics.clovm().count(prefixvm(mphysics.prefix(),"gravity-force").c_str()) )
            gravityStr = soption(_name="gravity-force",_prefix=mphysics.prefix(),_vm=mphysics.clovm());
        else if (Dim == 2 )
            gravityStr = "{0,-9.80665}";
        else if (Dim == 3 )
            gravityStr = "{0,0,-9.80665}";
        M_gravityForceExpr.setExpr( gravityStr,mphysics.worldComm(),mphysics.repository().expr() );
    }

    if ( j_model.contains("viscosity_law") )
    {
        auto const& j_model_viscosity_law = j_model.at("viscosity_law");
        CHECK( j_model_viscosity_law.is_string() ) << "viscosity_law must be a string";
            M_dynamicViscosity.setLaw( j_model_viscosity_law.template get<std::string>() );
    }

    if ( j_model.contains("turbulence") )
    {
        auto const& j_model_turbulence = j_model.at("turbulence");
        CHECK( j_model_turbulence.is_object() ) << "turbulence must be an json object";
        M_turbulence.setup( j_model_turbulence );
    }

}

template <uint16_type Dim>
void
ModelPhysicFluid<Dim>::setEquation( std::string const& eq )
{
    CHECK( eq == "Navier-Stokes" || eq == "Stokes" || eq == "StokesTransient" ) << "invalid equation of fluid : " << eq;
    M_equation = eq;
}

template <uint16_type Dim>
bool
ModelPhysicFluid<Dim>::Turbulence::useBoussinesqApproximation() const
{
    return (M_model == "Spalart-Allmaras") || (M_model == "k-epsilon");
}
template <uint16_type Dim>
bool
ModelPhysicFluid<Dim>::Turbulence::hasTurbulentKineticEnergy() const
{
    return false;// (M_model == "k-epsilon");
}

template <uint16_type Dim>
void
ModelPhysicFluid<Dim>::Turbulence::setup( nl::json const& jarg )
{
    if ( jarg.contains("enable") )
    {
        auto const& j_enable = jarg.at("enable");
        if ( j_enable.is_boolean() )
            M_isEnabled = j_enable.template get<bool>();
        else if ( j_enable.is_string() )
            M_isEnabled = boost::lexical_cast<bool>( j_enable.template get<std::string>() );
    }

    if ( !M_isEnabled )
        return;

    if ( jarg.contains("model") )
    {
        auto const& j_model = jarg.at("model");
        if ( j_model.is_string() )
            M_model = j_model.template get<std::string>();
        CHECK( M_model == "Spalart-Allmaras" || M_model == "k-epsilon" ) << "invalid turubulence model " << M_model;
    }

    this->setParameterNoMaterialNoModelLink( "kappa", "0.41" );
    // upper limit on the mixing length :
    // * default choice (see Comsol doc) : l_mix_lim = 0.5*lbb where lbb is the shortest side of the geometry bounding box.
    // * for complex geometry, l_mix_lim should be given manually
    //this->setParameterNoMaterialNoModelLink( "l_mix_lim", "0.5*0.0635" ); // turbulent pipe case
    //this->setParameterNoMaterialNoModelLink( "l_mix_lim", "0.0508" ); // backward facing step case
    this->setParameterNoMaterialNoModelLink( "l_mix_lim", "0.0381+0.0762" ); // backward facing step case

    this->setParameterNoMaterialLink( "c_mu", "0.09", "k-epsilon" );
    this->setParameterNoMaterialLink( "c_1epsilon", "1.44", "k-epsilon" );
    this->setParameterNoMaterialLink( "c_2epsilon", "1.92", "k-epsilon" );
    this->setParameterNoMaterialLink( "sigma_k", "1.0", "k-epsilon" );
    this->setParameterNoMaterialLink( "sigma_epsilon", "1.3", "k-epsilon" );

}

template <uint16_type Dim>
ModelPhysicSolid<Dim>::ModelPhysicSolid( ModelPhysics<Dim> const& mphysics, std::string const& name, ModelModel const& model )
    :
    super_type( "solid", name, mphysics, model ),
    M_equation( "Hyper-Elasticity" ),
    M_materialModel( soption(_prefix=mphysics.prefix(), _name="material_law",_vm=mphysics.clovm() ) ),
    M_formulation( soption(_prefix=mphysics.prefix(), _name="formulation",_vm=mphysics.clovm() ) ),
    M_decouplingEnergyVolumicLaw( "classic" ),
    M_compressibleNeoHookeanVariantName( "default" )
{
    if ( mphysics.clovm().count(prefixvm(mphysics.prefix(),"model").c_str()) )
        this->setEquation( soption(_name="model",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) );

    auto const& j_model = model.jsonProperties();
    for ( std::string const& eqarg : { "equations", "equation" } )
        if ( j_model.contains( eqarg ) )
        {
            auto const& j_model_eq = j_model.at( eqarg );
            if ( j_model_eq.is_string() )
                this->setEquation( j_model_eq.get<std::string>() );
        }

    if ( j_model.contains( "material-model" ) )
    {
        auto const& j_model_matmodel = j_model.at( "material-model" );
        if ( j_model_matmodel.is_string() )
            M_materialModel = j_model_matmodel.template get<std::string>();
    }
    CHECK( M_materialModel == "StVenantKirchhoff" || M_materialModel == "NeoHookean" ) << "invalid material-model :" << M_materialModel;

    if ( j_model.contains( "formulation" ) )
    {
        auto const& j_model_formulation = j_model.at( "formulation" );
        if ( j_model_formulation.is_string() )
            M_formulation = j_model_formulation.template get<std::string>();
    }
    if ( j_model.contains( "volumetric-strain-energy" ) )
    {
        auto const& j_model_volstrainenergy = j_model.at( "volumetric-strain-energy" );
        if ( j_model_volstrainenergy.is_string() )
            M_decouplingEnergyVolumicLaw = j_model_volstrainenergy.template get<std::string>();
    }
    if ( j_model.contains( "neo-Hookean.variant" ) )
    {
        auto const& j_model_neoHookeanVariant = j_model.at( "neo-Hookean.variant" );
        if ( j_model_neoHookeanVariant.is_string() )
            M_compressibleNeoHookeanVariantName = j_model_neoHookeanVariant.template get<std::string>();
    }

    CHECK( M_decouplingEnergyVolumicLaw == "classic" || M_decouplingEnergyVolumicLaw == "simo1985" ) << "invalid decouplingEnergyVolumicLaw : " << M_decouplingEnergyVolumicLaw;
    CHECK( M_compressibleNeoHookeanVariantName == "default" || M_compressibleNeoHookeanVariantName == "molecular-theory" ||
           M_compressibleNeoHookeanVariantName == "molecular-theory-simo1985" ) << "invalid compressibleNeoHookeanVariantName : " <<  M_compressibleNeoHookeanVariantName;

}

template <uint16_type Dim>
void
ModelPhysicSolid<Dim>::setEquation( std::string const& eq )
{
    CHECK( eq == "Hyper-Elasticity" || eq == "hyperelasticity" || eq == "Elasticity" || eq == "linear-elasticity" || eq == "Generalised-String" ) << "invalid equation of solid : " << eq;
    M_equation = eq;
#if 0
    if ( M_equation == "Elasticity" )
        M_equation = "linear-elasticity";
#endif
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::string const& name, ModelModels const& models, subphysic_description_type const& subPhyicsDesc )
{
    std::string const& type = M_physicType;

    M_physicDefault = name;
    auto const& theGlobalModel = models.model( name );
    if ( std::find_if( theGlobalModel.variants().begin(), theGlobalModel.variants().end(), [&name]( auto const& varMod ) { return name == varMod.first; } ) == theGlobalModel.variants().end() )
        M_physics.emplace( name, ModelPhysic<Dim>::New( *this, type, name, theGlobalModel ) );
    for ( auto const& [variantName,variantModel] : theGlobalModel.variants() )
        M_physics.emplace( variantName, ModelPhysic<Dim>::New( *this, type, variantName, variantModel ) );

    // list of subphysics from description
    std::map<std::string, std::set<std::string>> mapSubPhysicsTypeToDefaultNames;
    std::vector<std::shared_ptr<ModelPhysic<Dim>>> parentPhysicToUpdate;
    for ( auto const& [physicName,physicData] : M_physics )
    {
        for ( std::string const& subPhysicType : physicData->subphysicsTypes() )
        {
            auto itFindSubphysic = std::find_if( physicData->subphysics().begin(), physicData->subphysics().end(), [&subPhysicType]( auto const& e ) {
                    if ( !e.second )
                        return false;
                    return e.second->type() == subPhysicType;
                } );
            if ( itFindSubphysic != physicData->subphysics().end() )
                continue;

            auto itFindTypeInDesc = subPhyicsDesc.find( subPhysicType );
            CHECK( itFindTypeInDesc != subPhyicsDesc.end() ) << "type not given in subPhyicsDesc";
            std::string const& subPhysicDefaultName = std::get<0>( itFindTypeInDesc->second );
            mapSubPhysicsTypeToDefaultNames[subPhysicType].insert( subPhysicDefaultName );
            parentPhysicToUpdate.push_back( physicData );
        }
    }

    // create subphysics
    bool useModelNameInJson = models.useModelName();
    for ( auto const& [subPhysicType,subPhysicDefaultNames] : mapSubPhysicsTypeToDefaultNames )
    {
        auto mSubPhysics = std::get<1>( subPhyicsDesc.find( subPhysicType )->second );
        for ( std::string const& subName : subPhysicDefaultNames )
        {
            auto const& theSubModel = useModelNameInJson && models.hasModel( subName )? models.model( subName ) : ModelModel{};
            if ( std::find_if( theSubModel.variants().begin(), theSubModel.variants().end(), [&subName]( auto const& varMod ) { return subName == varMod.first; } ) == theSubModel.variants().end() )
                M_physics.emplace( subName, ModelPhysic<Dim>::New( *mSubPhysics, subPhysicType, subName, theSubModel ) );
            for ( auto const& [variantName,variantModel] : theSubModel.variants() )
                M_physics.emplace( variantName, ModelPhysic<Dim>::New( *mSubPhysics, subPhysicType, variantName, variantModel ) );
        }
    }

    // update subphysics connection
    for ( auto const& parentPhysic : parentPhysicToUpdate )
    {
        std::string const& parentType = parentPhysic->type();
        std::string const& parentName = parentPhysic->name();

        // if subphysics defined but not init
        auto const& subphysicRegistered = parentPhysic->subphysics();
        for ( auto const& subPhysicPair : subphysicRegistered )
        {
            std::string subPhysicName = subPhysicPair.first;
            auto const& subPhysicData = subPhysicPair.second;
            if ( subPhysicData )
                continue;
            auto itFindSubphysic = std::find_if( M_physics.begin(), M_physics.end(), [&subPhysicName]( auto const& e ) { return /*e.second->type() == subPhysicType &&*/ e.second->name() == subPhysicName; } );
            CHECK( itFindSubphysic != M_physics.end() ) << "subphysic not found";
            parentPhysic->addSubphysic( itFindSubphysic->second );
        }

        // init all subphysics
        for ( std::string const& subPhysicType : parentPhysic->subphysicsTypes() )
        {
            auto itFindType = std::find_if( subphysicRegistered.begin(), subphysicRegistered.end(), [&subPhysicType]( auto const& e ) { return e.second->type() == subPhysicType; } );
            if ( itFindType != subphysicRegistered.end() )
                continue;

            auto itFindTypeInDesc = subPhyicsDesc.find( subPhysicType );
            CHECK( itFindTypeInDesc != subPhyicsDesc.end() ) << "type not given in subPhyicsDesc";
            std::string const& subPhysicDefaultName = std::get<0>( itFindTypeInDesc->second );

            auto itFindSubphysic = std::find_if( M_physics.begin(), M_physics.end(), [&subPhysicType,&subPhysicDefaultName]( auto const& e ) { return e.second->type() == subPhysicType && e.second->name() == subPhysicDefaultName; } );
            CHECK( itFindSubphysic != M_physics.end() ) << "subphysic not found";
            parentPhysic->addSubphysic( itFindSubphysic->second );

        }
    }

}


template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsAvailable() const
{
    std::set<std::string> res;
    for ( auto const& [physicName,physicData] : M_physics )
        res.insert( physicName );
    return res;
}

template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsAvailable( std::string const& type ) const
{
    std::set<std::string> res;
    for ( auto const& [physicName,physicData] : M_physics )
        if ( physicData->type() == type )
            res.insert( physicName );
    return res;
}

template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsShared( std::string const& pname ) const
{
    std::set<std::string> res;
    auto itFindPhysic = M_physics.find( pname );
    if ( itFindPhysic != M_physics.end() )
        itFindPhysic->second->physicsShared( res );
    return res;
}

template <uint16_type Dim>
std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>>
    ModelPhysics<Dim>::physics( std::string const& type ) const
{
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> res;
    for ( auto const& [physicName,physicData] : M_physics )
    {
        if ( physicData->type() == type )
            res[physicName] = physicData;
    }
    return res;
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::setPhysics( std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics, std::string const& physicDefault )
{
    M_physics.insert( thePhysics.begin(), thePhysics.end() );

    if ( !physicDefault.empty() )
        M_physicDefault = physicDefault;
}

template class ModelPhysic<2>;
template class ModelPhysic<3>;
template class ModelPhysics<2>;
template class ModelPhysics<3>;
template class ModelPhysicFluid<2>;
template class ModelPhysicFluid<3>;

} // namespace FeelModels
} // namespace Feel
