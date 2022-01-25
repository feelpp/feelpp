/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelPhysic<Dim>::ModelPhysic( std::string const& modeling, std::string const& type, std::string const& name, ModelBase const& mparent, ModelModel const& model )
    :
    super_type( "ModelPhysic", type+"/"+name ),
    M_modeling( modeling ),
    M_type( type ),
    M_name( name ),
    M_worldComm( mparent.worldCommPtr() ),
    M_directoryLibExpr( mparent.repository().expr() )
{
    M_materialNames = model.materials();

    if ( type == "GenericPDE" || type == "GenericPDEs" )
        return;

    if ( M_type == "thermo-electric" )
    {
        M_mapSubphysicsTypeToNames["heat"];
        M_mapSubphysicsTypeToNames["electric"];
    }
    else if ( M_type == "heat-fluid" )
    {
        M_mapSubphysicsTypeToNames["heat"];
        M_mapSubphysicsTypeToNames["fluid"];
    }

    for ( auto const& [_type,_subm] : model.submodels() )
    {
        CHECK( M_mapSubphysicsTypeToNames.find( _type ) != M_mapSubphysicsTypeToNames.end() ) << "something wrong";
        M_mapSubphysicsTypeToNames[_type].insert(_subm.begin(),_subm.end() );
    }


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
void
ModelPhysic<Dim>::initSubphysics( ModelPhysics<Dim> & mphysics, subphysic_description_type const& subPhysicsDesc, ModelModel const& model, ModelModels const& models  )
{
#if 1 //TODO
    auto submodels = model.submodels();
    for ( auto & [_type,_msubphysics] : subPhysicsDesc )
    {
        std::set<std::string> submodelNames;
        auto itFindSubmodel = submodels.find(_type);
        if ( itFindSubmodel != submodels.end() )
            submodelNames.insert( itFindSubmodel->second.begin(), itFindSubmodel->second.end());

        if ( submodelNames.empty() && models.hasType( _type ) )
        {
            for ( auto const& [_name,_model] : models.models( _type ) )
                submodelNames.insert( _name );
        }

        if ( submodelNames.empty() ) // if empty add all physics with same type
        {
            if ( mphysics.physics( _type ).empty() ) // if no physic, create default physic
            {
                std::string _name = "default";
                auto _mphysic = ModelPhysic<Dim>::New( *_msubphysics, _msubphysics->physicModeling(),  _type, _name, ModelModel(_type,_name) );
                mphysics.addPhysicInternal( _mphysic, models, subPhysicsDesc );
            }
            for ( auto const& [pId,_mphysic] : mphysics.physics( _type ) )
                this->addSubphysic( _mphysic );
        }
        else
        {
            for ( std::string const& _name : submodelNames )
            {
                for ( auto const& [_name2,_model] : models.models( _type ) )
                {
                    if ( _name != _name2 )
                        continue;

                    auto pId = std::make_pair(_type, _name);
                    if ( !mphysics.hasPhysic( pId ) )
                    {
                        auto _mphysic = ModelPhysic<Dim>::New( *_msubphysics, _msubphysics->physicModeling(), _type, _name, _model );
                        mphysics.addPhysicInternal( _mphysic, models, subPhysicsDesc );
                    }
                    if ( !this->hasSubphysic( pId ) )
                    {
                        auto _mphysic = mphysics.physic( pId );
                        this->addSubphysic( _mphysic );
                    }
                }
            }
        }

    }
#endif

#if 0
    for ( auto & [_type,_subm] : M_mapSubphysicsTypeToNames )
    {
        if ( _subm.empty() && models.hasType( _type ) )
        {
            for ( auto const& [_name,_model] : models.models( _type ) )
                _subm.insert( _name );
        }

        CHECK( subPhysicsDesc.find( _type ) != subPhysicsDesc.end() ) << "aiai";
        auto mSubPhysics = std::get<1>( subPhysicsDesc.find( _type )->second );

        if ( _subm.empty() ) // if empty add all physics with same type
        {
            if ( mphysics.physics( _type ).empty() ) // if no physic, create default physic
            {
                std::string _name = "default";
                auto _mphysic = ModelPhysic<Dim>::New( *mSubPhysics, _type, _name, ModelModel(_type,_name) );
                mphysics.addPhysicInternal( _mphysic, models, subPhysicsDesc );
            }
            for ( auto const& [pId,_mphysic] : mphysics.physics( _type ) )
                this->addSubphysic( _mphysic );
        }
        else
        {
            for ( std::string const& _name : _subm )
            {
                for ( auto const& [_name2,_model] : models.models( _type ) )
                {
                    if ( _name != _name2 )
                        continue;

                    auto pId = std::make_pair(_type, _name);
                    if ( !mphysics.hasPhysic( pId ) )
                    {
                        auto _mphysic = ModelPhysic<Dim>::New( *mSubPhysics, _type, _name, _model );
                        mphysics.addPhysicInternal( _mphysic, models, subPhysicsDesc );
                    }
                    if ( !this->hasSubphysic( pId ) )
                    {
                        auto _mphysic = mphysics.physic( pId );
                        this->addSubphysic( _mphysic );
                    }
                 }
            }
         }
    }
#endif

}

template <uint16_type Dim>
std::shared_ptr<ModelPhysic<Dim>>
ModelPhysic<Dim>::New( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
{
    if ( modeling == "heat" )
        return std::make_shared<ModelPhysicHeat<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "fluid" )
        return std::make_shared<ModelPhysicFluid<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "solid" )
        return std::make_shared<ModelPhysicSolid<Dim>>( mphysics, modeling, type, name, model );
    else
        return std::make_shared<ModelPhysic<Dim>>( modeling, type, name, mphysics, model );
}

template <uint16_type Dim>
ModelPhysicHeat<Dim>::ModelPhysicHeat( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    auto const& j_setup = model.setup();
    if ( j_setup.contains( "convection" ) )
    {
        M_convection.emplace( this );
        M_convection->setup( j_setup.at("convection" ) );
    }

    if ( j_setup.contains( "heat-sources" ) )
    {
        auto const& j_setup_heatsources = j_setup.at( "heat-sources" );
        if ( j_setup_heatsources.is_array() )
        {
            for ( auto const& [j_setup_heatsourceskey,j_setup_heatsourcesval] : j_setup_heatsources.items() )
            {
                CHECK( j_setup_heatsourcesval.is_object() ) << "j_setup_heatsourcesval should be an object";
                HeatSource hs( this, fmt::format("heatsource{}",M_heatSources.size()) );
                hs.setup( j_setup_heatsourcesval );
                M_heatSources.push_back( std::move( hs ) );
            }
        }
        else if ( j_setup_heatsources.is_object() )
        {
            HeatSource hs( this, fmt::format("heatsource{}",M_heatSources.size()) );
            hs.setup( j_setup_heatsources );
            M_heatSources.push_back( std::move( hs ) );
        }
    }
}

template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::updateInformationObject( nl::json & p ) const
{
    // TODO move in base class
    p["type"] = this->type();
    p["name"] = this->name();
}


template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::Convection::setup( nl::json const& jarg )
{
    std::string convectionExprName = "convection";
    if ( jarg.is_string() )
    {
        M_parent->addParameter( convectionExprName, jarg );
    }
    else if ( jarg.is_object() )
    {
        CHECK( jarg.contains( "expr" ) ) << "no expr";
        M_parent->addParameter( convectionExprName, jarg.at("expr") );
    }
    if ( !M_parent->template hasParameterExpr<Dim,1>( convectionExprName ) )
        CHECK( false ) << "invalid convection setup";
    M_enabled = true;
}

template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::HeatSource::setup( nl::json const& jarg )
{
    CHECK( jarg.contains( "expr" ) ) << "no expr";
    M_parent->addParameter( M_name/*"heatsource"*/, jarg.at("expr") );
}


template <uint16_type Dim>
ModelPhysicFluid<Dim>::ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model ),
    M_equation( "Navier-Stokes" ),
    M_gravityForceEnabled( boption(_name="use-gravity-force",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) ),
    M_dynamicViscosity( soption(_name="viscosity.law",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) ),
    M_turbulence( this )
{

    auto const& j_setup = model.setup();
    if ( j_setup.is_string() )
    {
        this->setEquation( j_setup.get<std::string>() );
    }
    else if ( j_setup.is_object() )
    {
        // equation
        for ( std::string const& eqarg : { "equations", "equation" } )
            if ( j_setup.contains( eqarg ) )
            {
                auto const& j_setup_eq = j_setup.at( eqarg );
                if ( j_setup_eq.is_string() )
                    this->setEquation( j_setup_eq.get<std::string>() );
            }

        // gravity force
        if ( j_setup.contains("gravity") )
        {
            auto const& j_setup_gravity = j_setup.at("gravity");
            if ( j_setup_gravity.is_boolean() )
                M_gravityForceEnabled = j_setup_gravity.template get<bool>();
            else if ( j_setup_gravity.is_string() )
                M_gravityForceEnabled = boost::lexical_cast<bool>( j_setup_gravity.template get<std::string>() );
            else if ( j_setup_gravity.is_object() )
            {
                if ( j_setup_gravity.contains( "enable" ) )
                {
                    auto const& j_setup_gravity_enable =  j_setup_gravity.at("enable");
                    if ( j_setup_gravity_enable.is_boolean() )
                        M_gravityForceEnabled = j_setup_gravity_enable.template get<bool>();
                    else if ( j_setup_gravity_enable.is_string() )
                        M_gravityForceEnabled = boost::lexical_cast<bool>( j_setup_gravity_enable.template get<std::string>() );
                }
                if ( j_setup_gravity.contains( "expr" ) )
                    M_gravityForceExpr.setExpr( j_setup_gravity.at("expr"),mphysics.worldComm(),mphysics.repository().expr() );
            }
        }

        if ( j_setup.contains("viscosity_law") )
        {
            auto const& j_setup_viscosity_law = j_setup.at("viscosity_law");
            CHECK( j_setup_viscosity_law.is_string() ) << "viscosity_law must be a string";
            M_dynamicViscosity.setLaw( j_setup_viscosity_law.template get<std::string>() );
        }

        if ( j_setup.contains("turbulence") )
        {
            auto const& j_setup_turbulence = j_setup.at("turbulence");
            CHECK( j_setup_turbulence.is_object() ) << "turbulence must be an json object";
            M_turbulence.setup( j_setup_turbulence );
        }
    } // is_object

    // update for use
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
ModelPhysicSolid<Dim>::ModelPhysicSolid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model ),
    M_equation( "Hyper-Elasticity" ),
    M_materialModel( soption(_prefix=mphysics.prefix(), _name="material_law",_vm=mphysics.clovm() ) ),
    M_formulation( soption(_prefix=mphysics.prefix(), _name="formulation",_vm=mphysics.clovm() ) ),
    M_decouplingEnergyVolumicLaw( "classic" ),
    M_compressibleNeoHookeanVariantName( "default" )
{
    if ( mphysics.clovm().count(prefixvm(mphysics.prefix(),"model").c_str()) )
        this->setEquation( soption(_name="model",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) );

    auto const& j_setup = model.setup();
    for ( std::string const& eqarg : { "equations", "equation" } )
        if ( j_setup.contains( eqarg ) )
        {
            auto const& j_setup_eq = j_setup.at( eqarg );
            if ( j_setup_eq.is_string() )
                this->setEquation( j_setup_eq.get<std::string>() );
        }

    if ( j_setup.contains( "material-model" ) )
    {
        auto const& j_setup_matmodel = j_setup.at( "material-model" );
        if ( j_setup_matmodel.is_string() )
            M_materialModel = j_setup_matmodel.template get<std::string>();
    }
    CHECK( M_materialModel == "StVenantKirchhoff" || M_materialModel == "NeoHookean" ) << "invalid material-model :" << M_materialModel;

    if ( j_setup.contains( "formulation" ) )
    {
        auto const& j_setup_formulation = j_setup.at( "formulation" );
        if ( j_setup_formulation.is_string() )
            M_formulation = j_setup_formulation.template get<std::string>();
    }
    if ( j_setup.contains( "volumetric-strain-energy" ) )
    {
        auto const& j_setup_volstrainenergy = j_setup.at( "volumetric-strain-energy" );
        if ( j_setup_volstrainenergy.is_string() )
            M_decouplingEnergyVolumicLaw = j_setup_volstrainenergy.template get<std::string>();
    }
    if ( j_setup.contains( "neo-Hookean.variant" ) )
    {
        auto const& j_setup_neoHookeanVariant = j_setup.at( "neo-Hookean.variant" );
        if ( j_setup_neoHookeanVariant.is_string() )
            M_compressibleNeoHookeanVariantName = j_setup_neoHookeanVariant.template get<std::string>();
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
ModelPhysics<Dim>::addPhysicInternal( model_physic_ptrtype mphysic, ModelModels const& models, subphysic_description_type const& subPhysicsDesc )
{
    M_physics.emplace( std::make_pair(mphysic->type(), mphysic->name()), mphysic );
#if 0
    mphysic->initSubphysics( *this, subPhysicsDesc, models );
#endif
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::string const& type, ModelModels const& models, subphysic_description_type const& subPhyicsDesc )
{
    M_physicType = type;
    //std::string const& type = M_physicType;

    if ( models.hasType( type ) )
    {
        for ( auto const& [_name,_model] : models.models( type ) )
        {
            if ( this->hasPhysic( std::make_pair(type, _name) ) )
                continue;
            auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
            M_physics.emplace( std::make_pair(type, _name), _mphysic );
            _mphysic->initSubphysics( *this, subPhyicsDesc, _model, models );
        }
    }
    else if ( this->physics( type ).empty() )
    {
        // create default physic
        std::string _name = "default";
        auto _model = ModelModel(type,_name);
        auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
        M_physics.emplace( std::make_pair(type, _name), _mphysic );
        _mphysic->initSubphysics( *this, subPhyicsDesc, _model, models );
    }
}


template <uint16_type Dim>
std::set<typename ModelPhysics<Dim>::physic_id_type>
ModelPhysics<Dim>::physicsAvailable() const
{
    std::set<physic_id_type> res;
    for ( auto const& [physicName,physicData] : M_physics )
        res.insert( physicName );
    return res;
}

template <uint16_type Dim>
std::set<typename ModelPhysics<Dim>::physic_id_type>
ModelPhysics<Dim>::physicsAvailable( std::string const& type ) const
{
    std::set<physic_id_type> res;
    for ( auto const& [physicName,physicData] : M_physics )
        if ( physicData->type() == type )
            res.insert( physicName );
    return res;
}
#if 0
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
#endif
template <uint16_type Dim>
std::map<typename ModelPhysics<Dim>::physic_id_type,std::shared_ptr<ModelPhysic<Dim>>>
ModelPhysics<Dim>::physics( std::string const& type ) const
{
    std::map<physic_id_type,std::shared_ptr<ModelPhysic<nDim>>> res;
    for ( auto const& [physicName,physicData] : M_physics )
    {
        if ( physicData->type() == type )
            res[physicName] = physicData;
    }
    return res;
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::setPhysics( std::map<physic_id_type,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics )
{
    M_physics.insert( thePhysics.begin(), thePhysics.end() );
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::updateInformationObjectFromCurrentType( nl::json & p ) const
{
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        p[physicName.first][physicName.second] = physicData->journalSection().to_string();
}


template class ModelPhysic<2>;
template class ModelPhysic<3>;
template class ModelPhysics<2>;
template class ModelPhysics<3>;
template class ModelPhysicHeat<2>;
template class ModelPhysicHeat<3>;
template class ModelPhysicFluid<2>;
template class ModelPhysicFluid<3>;

} // namespace FeelModels
} // namespace Feel
