/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelmarkers.hpp>


namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelPhysic<Dim>::ModelPhysic( std::string const& modeling, std::string const& type, std::string const& name, ModelBase const& mparent, ModelModel const& model )
    :
    super_type( "ModelPhysic"/*, type+"/"+name*/ ), // not put the name or need to use also object counter
    M_modeling( modeling ),
    M_type( type ),
    M_name( name ),
    M_worldComm( mparent.worldCommPtr() ),
    M_directoryLibExpr( mparent.repository().expr() )
{
    M_materialNames = model.materials();

    if ( M_modeling == "GenericPDE" || M_modeling == "GenericPDEs" )
        return;

    //std::cout << "CREATE physics :  " << "M_modeling=" << M_modeling << " M_type=" << M_type << " M_name=" << M_name << std::endl;

    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    this->addMaterialPropertyDescription( "density", "rho", { scalarShape } );
    if ( M_modeling == "heat" ||  M_modeling == "thermo-electric" || M_modeling == "heat-fluid" )
    {
        this->addMaterialPropertyDescription( "specific-heat-capacity", "Cp", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-expansion", "beta", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-conductivity", "k", { scalarShape,matrixShape } );
    }
    if ( M_modeling == "electric" || M_modeling == "thermo-electric" )
    {
        this->addMaterialPropertyDescription( "electric-conductivity", "sigma", { scalarShape } );
    }
    if ( M_modeling == "solid" || M_modeling == "fsi" )
    {
        this->addMaterialPropertyDescription( "Young-modulus", "E", { scalarShape } );
        this->addMaterialPropertyDescription( "Poisson-ratio", "nu", { scalarShape } );
        this->addMaterialPropertyDescription( "Lame-first-parameter", "lambdaLame", { scalarShape } );
        this->addMaterialPropertyDescription( "Lame-second-parameter", "muLame", { scalarShape } );
        this->addMaterialPropertyDescription( "bulk-modulus", "K", { scalarShape } );
    }
    if ( M_modeling == "fluid" || M_modeling == "heat-fluid" || M_modeling == "fsi" )
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
ModelPhysic<Dim>::New( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
{
    if ( modeling == "heat" )
        return std::make_shared<ModelPhysicHeat<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "electric" )
        return std::make_shared<ModelPhysicElectric<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "thermo-electric" )
        return std::make_shared<ModelPhysicThermoElectric<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "fluid" )
        return std::make_shared<ModelPhysicFluid<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "solid" )
        return std::make_shared<ModelPhysicSolid<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "fsi" )
        return std::make_shared<ModelPhysicFSI<Dim>>( mphysics, modeling, type, name, model );
    else
        return std::make_shared<ModelPhysic<Dim>>( modeling, type, name, mphysics, model );
}

template <uint16_type Dim>
void
ModelPhysic<Dim>::updateInformationObject( nl::json & p ) const
{
    p["modeling"] = this->modeling();
    p["type"] = this->type();
    p["name"] = this->name();
    if ( !M_materialNames.empty() )
        p["materials"] = M_materialNames;

    if ( !this->subphysics().empty() )
    {
        std::set<std::string> subtypes;
        for ( auto const& [pId,subphysic] : this->subphysics() )
            subtypes.insert( subphysic->type() );

        for ( std::string const& subtype : subtypes )
        {
            auto subphysics = this->subphysicsFromType( subtype );
            nl::json::array_t jaSubphysicSameType;
            for ( auto const& subphysic : subphysics )
                jaSubphysicSameType.push_back( nl::json({ { subphysic->name(), subphysic->journalSection().to_string() } }) );
            p["subphysics"][subtype] = jaSubphysicSameType;
        }
    }

    if ( !M_parameterNameToExpr.empty() )
    {
        nl::json::array_t jaParams;
        for ( auto const& [pname,pmexpr] : M_parameterNameToExpr )
        {
            std::string symbolParam = this->symbolFromParameter( pname );
            auto [exprStr,compInfo] = pmexpr.exprInformations();
            jaParams.push_back( symbolExprInformations( symbolParam, exprStr, compInfo, pname ) );
        }
        p["parameters"] = jaParams;
    }
}

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysic<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    this->updateTabulateInformationsBasic( jsonInfo, tabInfo, tabInfoProp );
    this->updateTabulateInformationsSubphysics( jsonInfo, tabInfo, tabInfoProp );
    this->updateTabulateInformationsParameters( jsonInfo, tabInfo, tabInfoProp );
    return tabInfo;
}

template <uint16_type Dim>
void
ModelPhysic<Dim>::updateTabulateInformationsBasic( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    Feel::Table tabInfoGeneric;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfoGeneric, jsonInfo, tabInfoProp, { "modeling","type","name" } );
    if ( jsonInfo.contains( "materials" ) )
    {
        auto j_mats = jsonInfo.at( "materials" );
        Feel::Table tabInfoMaterials;
        if ( j_mats.is_array() )
            tabInfoMaterials = TabulateInformationTools::FromJSON::createTableFromArray( j_mats, true );
        if ( tabInfoMaterials.nRow() > 0 )
            tabInfoGeneric.add_row( { "materials", tabInfoMaterials } );
    }
    tabInfoGeneric.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    if ( tabInfoGeneric.nRow() > 0 )
        tabInfo->add( "", TabulateInformations::New( tabInfoGeneric, tabInfoProp ) );
}
template <uint16_type Dim>
void
ModelPhysic<Dim>::updateTabulateInformationsSubphysics( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    if ( !jsonInfo.contains("subphysics") )
        return;

    auto tabInfoByType = TabulateInformationsSections::New( tabInfoProp );
    for ( auto const& [j_subphysicskey,subphysicsval] : jsonInfo.at("subphysics").items() )
    {
        if ( !subphysicsval.is_array() )
            continue;

        Feel::Table tabInfoNameToJournalRef;
        for ( auto const& [j_subphysics_namekey,j_subphysics_nameval] : subphysicsval.items() )
            TabulateInformationTools::FromJSON::addAllKeyToValues(tabInfoNameToJournalRef,j_subphysics_nameval,tabInfoProp);
        tabInfoNameToJournalRef.format()
            .setShowAllBorders( false )
            .setColumnSeparator(":")
            .setHasRowSeparator( false );
        tabInfoByType->add( j_subphysicskey, TabulateInformations::New( tabInfoNameToJournalRef, tabInfoProp ) );
    }
    tabInfo->add( "Subphysics", tabInfoByType );
}
template <uint16_type Dim>
void
ModelPhysic<Dim>::updateTabulateInformationsParameters( nl::json const& jsonInfo, tabulate_informations_sections_ptr_t & tabInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    if ( jsonInfo.contains("parameters") )
        tabInfo->add( "Parameters", TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( jsonInfo.at("parameters"),tabInfoProp,true) );
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
ModelPhysicHeat<Dim>::setupConvection( nl::json & jarg )
{
    if ( !M_convection )
        M_convection.emplace( this );
    M_convection->setup( jarg );
}

template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );

    nl::json & pHeat = p["Heat"];
    std::string eqTermConvection, eqTermSource = "0", eqTermTimeDerivative;
    if ( this->hasConvectionEnabled() )
    {
        eqTermConvection = "u dot nabla T";
        M_convection->updateInformationObject( pHeat["Convection"] );
    }
    if ( !M_heatSources.empty() )
    {
        eqTermSource = "f";
        nl::json & pHeatSources = pHeat["HeatSources"];
        for ( auto const& hs : M_heatSources )
        {
            nl::json phs;
            hs.updateInformationObject( phs );
            pHeatSources.push_back( std::move( phs ) );
        }
    }
    std::string heateq;
    if ( !eqTermTimeDerivative.empty() )
        heateq += eqTermTimeDerivative;
    if ( !heateq.empty() )
        heateq += " + ";
    if ( !eqTermConvection.empty() )
        heateq += eqTermConvection;
    heateq += " - div( k grad T ) = " + eqTermSource;
    pHeat["Equation"] = heateq;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicHeat<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        //tabInfo->add( "", super_type::tabulateInformations( jsonInfo.at("Generic"), tabInfoProp ) );
    }
    if ( jsonInfo.contains("Heat") )
    {
        auto const& jsonInfoHeat = jsonInfo.at("Heat");
        Feel::Table tabInfoHeatEquation;
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfoHeatEquation, jsonInfoHeat, tabInfoProp, { "Equation" } );
        tabInfo->add( "", TabulateInformations::New( tabInfoHeatEquation, tabInfoProp ) );

        if ( jsonInfoHeat.contains("Convection") )
            tabInfo->add( "Convection", Convection::tabulateInformations( jsonInfoHeat.at("Convection"), tabInfoProp ) );
        if ( jsonInfoHeat.contains("HeatSources") )
        {
            auto tabInfoHeatSources = TabulateInformationsSections::New( tabInfoProp );
            for ( auto const& [hskey,hsval] : jsonInfoHeat.at("HeatSources").items() )
                tabInfoHeatSources->add( "", HeatSource::tabulateInformations( hsval, tabInfoProp ) );
            tabInfo->add( "Heat Sources", tabInfoHeatSources );
        }
    }

    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }

    return tabInfo;
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
ModelPhysicHeat<Dim>::Convection::updateInformationObject( nl::json & p ) const
{
    if ( !this->enabled() )
        return;
    auto [exprStr,compInfo] = M_parent->parameterModelExpr("convection").exprInformations();
    p["expr"] = exprStr;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicHeat<Dim>::Convection::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    if ( jsonInfo.contains("expr") )
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "expr" } );
    return TabulateInformations::New( tabInfo, tabInfoProp );
}

template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::HeatSource::setup( nl::json const& jarg )
{
    if ( jarg.contains("type") )
    {
        auto const& j_type = jarg.at("type");
        if ( j_type.is_string() )
            M_type = j_type.template get<std::string>();
    }
    CHECK( M_type == "heat-source" || M_type == "heat-rate" ) << "invalid heat source type " << M_type;
    CHECK( jarg.contains( "expr" ) ) << "no expr";

    if ( M_type == "heat-source" )
    {
        M_parent->addParameter( M_name + "_heatsource", jarg.at("expr") );
    }
    else if ( M_type == "heat-rate" )
    {
        std::string symbHeatRate = M_parent->addParameter( M_name + "_heatrate", jarg.at("expr") );
        std::string symbMeasureMat = M_parent->symbolFromParameter( "measure_materials" );
        std::string exprHeatSource = fmt::format("{0}/{1}:{0}:{1}",symbHeatRate,symbMeasureMat);
        M_parent->addParameter( M_name + "_heatsource", exprHeatSource );
    }
}

template <uint16_type Dim>
void
ModelPhysicHeat<Dim>::HeatSource::updateInformationObject( nl::json & p ) const
{
    p["name"] = M_name;
    p["type"] = M_type;
    if ( this->givenAsHeatRate() )
    {
        auto [exprStr,compInfo] = M_parent->parameterModelExpr(M_name+"_heatrate").exprInformations();
        p["expr"] = exprStr;
    }
    else
    {
        auto [exprStr,compInfo] = M_parent->parameterModelExpr(M_name+"_heatsource").exprInformations();
        p["expr"] = exprStr;
    }
}

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicHeat<Dim>::HeatSource::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "name","type","expr" } );
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    return TabulateInformations::New( tabInfo, tabInfoProp );
}



template <uint16_type Dim>
ModelPhysicElectric<Dim>::ModelPhysicElectric( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    auto const& j_setup = model.setup();
}

template <uint16_type Dim>
void
ModelPhysicElectric<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicElectric<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    return tabInfo;
}

template <uint16_type Dim>
ModelPhysicThermoElectric<Dim>::ModelPhysicThermoElectric( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    auto const& j_setup = model.setup();
}

template <uint16_type Dim>
void
ModelPhysicThermoElectric<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicThermoElectric<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    return tabInfo;
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
ModelPhysicFluid<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicFluid<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    return tabInfo;
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
    if ( j_setup.is_string() )
    {
        this->setEquation( j_setup.get<std::string>() );
    }
    else if ( j_setup.is_object() )
    {
        for ( std::string const& eqarg : { "equations", "equation" } )
            if ( j_setup.contains( eqarg ) )
            {
                auto const& j_setup_eq = j_setup.at( eqarg );
                if ( j_setup_eq.is_string() )
                    this->setEquation( j_setup_eq.get<std::string>() );
            }

        for ( std::string const& eqarg : { "body-forces" } )
            if ( j_setup.contains( eqarg ) )
            {
                auto const& j_setup_bodyforces = j_setup.at( eqarg );
                if ( j_setup_bodyforces.is_array() )
                {
                    for ( auto const& [j_setup_bodyforceskey,j_setup_bodyforcesval] : j_setup_bodyforces.items() )
                    {
                        CHECK( j_setup_bodyforcesval.is_object() ) << "j_setup_bodyforcesval should be an object";
                        BodyForces bf( this, fmt::format("bodyforce{}",M_bodyForces.size()) );
                        bf.setup( j_setup_bodyforcesval );
                        M_bodyForces.push_back( std::move( bf ) );
                    }
                }
                else
                {
                    BodyForces bf( this, fmt::format("bodyforce{}",M_bodyForces.size()) );
                    bf.setup( j_setup_bodyforces );
                    M_bodyForces.push_back( std::move( bf ) );
                }
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
    }

    CHECK( M_decouplingEnergyVolumicLaw == "classic" || M_decouplingEnergyVolumicLaw == "simo1985" ) << "invalid decouplingEnergyVolumicLaw : " << M_decouplingEnergyVolumicLaw;
    CHECK( M_compressibleNeoHookeanVariantName == "default" || M_compressibleNeoHookeanVariantName == "molecular-theory" ||
           M_compressibleNeoHookeanVariantName == "molecular-theory-simo1985" ) << "invalid compressibleNeoHookeanVariantName : " <<  M_compressibleNeoHookeanVariantName;

}

template <uint16_type Dim>
void
ModelPhysicSolid<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );

    nl::json & pSolid = p["Solid"];
    pSolid["Equation"] = M_equation;

    if ( !M_bodyForces.empty() )
    {
        //eqTermSource = "f";
        nl::json & pBodyForces = pSolid["BodyForces"];
        for ( auto const& bf : M_bodyForces )
        {
            nl::json pbf;
            bf.updateInformationObject( pbf );
            pBodyForces.push_back( std::move( pbf ) );
        }
    }

}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicSolid<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
#if 0
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    if ( jsonInfo.contains("Solid") )
    {
        Feel::Table tabInfoSolidEquation;
        auto const& jsonInfoSolid = jsonInfo.at("Solid");
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfoSolidEquation, jsonInfoSolid, tabInfoProp, { "Equation" } );
        tabInfo->add( "", TabulateInformations::New( tabInfoSolidEquation, tabInfoProp ) );
    }
    return tabInfo;
#endif



    
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    if ( jsonInfo.contains("Solid") )
    {
        auto const& jsonInfoSolid = jsonInfo.at("Solid");
        Feel::Table tabInfoSolidEquation;
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfoSolidEquation, jsonInfoSolid, tabInfoProp, { "Equation" } );
        tabInfo->add( "", TabulateInformations::New( tabInfoSolidEquation, tabInfoProp ) );

        if ( jsonInfoSolid.contains("BodyForces") )
        {
            auto tabInfoBodyForces = TabulateInformationsSections::New( tabInfoProp );
            for ( auto const& [bfkey,bfval] : jsonInfoSolid.at("BodyForces").items() )
                tabInfoBodyForces->add( "", BodyForces::tabulateInformations( bfval, tabInfoProp ) );
            tabInfo->add( "Body Forces", tabInfoBodyForces );
        }

    }

    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    return tabInfo;
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
ModelPhysicSolid<Dim>::BodyForces::setup( nl::json const& jarg )
{
    if ( jarg.is_string() )
    {
        M_parent->addParameter( M_name + "_bodyforce", jarg );
    }
    else if ( jarg.is_object() )
    {
        CHECK( jarg.contains( "expr" ) ) << "no expr";
        M_parent->addParameter( M_name + "_bodyforce", jarg.at("expr") );
    }
}
template <uint16_type Dim>
void
ModelPhysicSolid<Dim>::BodyForces::updateInformationObject( nl::json & p ) const
{
    p["name"] = M_name;
    auto [exprStr,compInfo] = M_parent->parameterModelExpr(M_name+"_bodyforce").exprInformations();
    p["expr"] = exprStr;
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicSolid<Dim>::BodyForces::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "name","expr" } );
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    return TabulateInformations::New( tabInfo, tabInfoProp );
}

template <uint16_type Dim>
ModelPhysicFSI<Dim>::ModelPhysicFSI( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    auto const& j_setup = model.setup();
    if ( j_setup.contains( "interface" ) )
    {
        ModelMarkers markers;
        markers.setup( j_setup.at("interface")/*, indexes*/ );
        M_interfaceFluid = markers;
        M_interfaceSolid = markers;
    }
}

template <uint16_type Dim>
void
ModelPhysicFSI<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );

    nl::json & pFSI = p["FSI"];
    pFSI["interface_fluid"] = M_interfaceFluid;
    pFSI["interface_solid"] = M_interfaceSolid;
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicFSI<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Generic") )
    {
        super_type::updateTabulateInformationsBasic( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsSubphysics( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
        super_type::updateTabulateInformationsParameters( jsonInfo.at("Generic"), tabInfo, tabInfoProp );
    }
    return tabInfo;
}


template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::string const& type, ModelModels const& models, bool isRequired )
{
    M_physicType = type;

    if ( models.hasType( type ) )
    {
        for ( auto const& [_nameFromModels,_model] : models.models( type ) )
        {
            std::string _name = _nameFromModels.empty()? type : _nameFromModels; // use name=type if not given
            if ( this->hasPhysic( std::make_pair(type, _name) ) )
                continue;
            auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
            M_physics.emplace( std::make_pair(type, _name), _mphysic );
        }
    }
    else if ( this->physics( type ).empty() && isRequired )
    {
        // create default physic
        std::string _name = "default";
        auto _model = ModelModel(type,_name);
        auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
        M_physics.emplace( std::make_pair(type, _name), _mphysic );
    }
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( PhysicsTree const& physicTree, ModelModels const& models )
{
    auto allMPhysics = physicTree.allPhysics();
    for ( auto const& [type,mphysics,isRequired] : allMPhysics )
    {
        if ( mphysics )
            mphysics->initPhysics( type, models, isRequired );
    }

    // up subphysics by using recursive lambda function
    auto upSubphysics = [&models]( PhysicsTree const& _physicTree )
                            {
                                auto upSubphysicsImpl = [&models]( PhysicsTree const& _physicTreeBis, const auto& ff ) -> void
                                                            {
                                                                auto const& [rootType,rootPhysics,rootIsRequiredxs] = _physicTreeBis.root();
                                                                std::map<std::string,std::map<std::string,std::set<std::string>>> submodelsFromRoot; // subType -> ( rootName -> ( subNames ) )
                                                                if ( models.hasType( rootType ) )
                                                                    for ( auto const& [_rootnameFromModels,_rootmodel] : models.models( rootType ) )
                                                                    {
                                                                        std::string _rootname = _rootnameFromModels.empty()? rootType : _rootnameFromModels; // use name=type if not given
                                                                        for ( auto const& [_subtype,_subnames] : _rootmodel.submodels() )
                                                                            if ( !_subnames.empty() )
                                                                                submodelsFromRoot[_subtype][_rootname] = _subnames;
                                                                    }

                                                                for ( auto st : _physicTreeBis.subtrees() )
                                                                {
                                                                    ff( st, ff );

                                                                    auto const& [subType,subPhysics,subIsRequired] = st.root();
                                                                    rootPhysics->setPhysics( std::get<1>(st.root())->physics() );

                                                                    auto itFind = submodelsFromRoot.find( subType );
                                                                    if ( itFind != submodelsFromRoot.end() ) // link defined in models (i.e. from json)
                                                                    {
                                                                        for ( auto const& [_rootname,_subnames] : itFind->second )
                                                                        {
                                                                            auto rootPhysicId = std::make_pair( rootType, _rootname );
                                                                            auto rootPhysic = rootPhysics->physic( rootPhysicId );

                                                                            for ( auto const& _subname : _subnames )
                                                                            {
                                                                                auto subPhysic = rootPhysics->physic( std::make_pair( subType, _subname ) );
                                                                                rootPhysic->addSubphysic( subPhysic );
                                                                            }
                                                                        }
                                                                    }
                                                                    else // not defined, add all physic with same type
                                                                    {
                                                                        for ( auto const& [_rootPhysicId,_rootPhysic] : rootPhysics->physics( rootType ) )
                                                                            for ( auto const& [_subPhysicId,_subPhysic] : rootPhysics->physics( subType ) )
                                                                                _rootPhysic->addSubphysic( _subPhysic );
                                                                    }
                                                                }

                                                                // update maybe materials names
                                                                for ( auto const& [_rootPhysicId,_rootPhysic] : rootPhysics->physics( rootType ) )
                                                                {
                                                                    if ( _rootPhysic->subphysics().empty() )
                                                                        continue;
                                                                    auto const& rootMatNames = _rootPhysic->materialNames();
                                                                    if ( !rootMatNames.empty() )
                                                                        continue;

                                                                    std::map<std::string,std::set<std::string>> mapSubphysicTypeToMatNames;
                                                                    for ( auto const& [pId,subphysic] : _rootPhysic->subphysics() )
                                                                    {
                                                                        auto const& subMatNames = subphysic->materialNames();
                                                                        if ( !subMatNames.empty() ) // warning, only non empty
                                                                            mapSubphysicTypeToMatNames[subphysic->type()].insert( subMatNames.begin(), subMatNames.end() );
                                                                    }

                                                                    if ( mapSubphysicTypeToMatNames.empty() ) // ok, do nothing
                                                                        continue;

                                                                    // take intersection with subphysics
                                                                    std::set<std::string> intersecMatNames; bool isFirst = true;
                                                                    for ( auto const& [_t,_matnames] : mapSubphysicTypeToMatNames )
                                                                    {
                                                                        if ( isFirst )
                                                                        {
                                                                            intersecMatNames = _matnames;
                                                                            isFirst = false;
                                                                        }
                                                                        else
                                                                        {
                                                                            std::set<std::string> res;
                                                                            std::set_intersection(intersecMatNames.begin(), intersecMatNames.end(), _matnames.begin(), _matnames.end(),
                                                                                                  std::inserter(res, res.begin()));
                                                                            intersecMatNames = res;
                                                                        }
                                                                    }
                                                                    CHECK( !intersecMatNames.empty() ) << "incompatible materials names";
                                                                    _rootPhysic->addMaterialNames( intersecMatNames );
                                                                }

                                                                // check materials names compatibility
                                                            };
                                upSubphysicsImpl( _physicTree, upSubphysicsImpl );
                            };

    upSubphysics( physicTree );

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

template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysics<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New();

    for ( auto const& [jsonInfokey,jsonInfoval] : jsonInfo.items() )
    {
        auto tabInfoByType = TabulateInformationsSections::New();
        std::string const& type = jsonInfokey;
        for ( auto const& [jsonInfovalkey,jsonInfovalval] : jsonInfoval.items() )
        {
            std::string const& name = jsonInfovalkey;
            auto physicId = std::make_pair(type,name);
            if ( !this->hasPhysic( physicId ) )
                continue;
            auto thephysic = this->physic( physicId );
            nl::json::json_pointer jsonPointerInfoPhysic( jsonInfovalval.template get<std::string>() );
            if ( JournalManager::journalData().contains( jsonPointerInfoPhysic ) )
            {
                auto const& jsonInfoPhysic = JournalManager::journalData().at( jsonPointerInfoPhysic );
                tabInfoByType->add( name, thephysic->tabulateInformations( jsonInfoPhysic, tabInfoProp ) );
            }
        }
        tabInfo->add( type, tabInfoByType );
    }

    return tabInfo;
}


template class ModelPhysic<2>;
template class ModelPhysic<3>;
template class ModelPhysics<2>;
template class ModelPhysics<3>;
template class ModelPhysicHeat<2>;
template class ModelPhysicHeat<3>;
template class ModelPhysicElectric<2>;
template class ModelPhysicElectric<3>;
template class ModelPhysicThermoElectric<2>;
template class ModelPhysicThermoElectric<3>;
template class ModelPhysicFluid<2>;
template class ModelPhysicFluid<3>;
template class ModelPhysicSolid<2>;
template class ModelPhysicSolid<3>;
template class ModelPhysicFSI<2>;
template class ModelPhysicFSI<3>;

} // namespace FeelModels
} // namespace Feel
