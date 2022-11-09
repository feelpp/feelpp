/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

#include <feel/feelmodels/modelcore/modelgenericpde.hpp>

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
    else if ( modeling == "multifluid" )
        return std::make_shared<ModelPhysicMultifluid<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "solid" )
        return std::make_shared<ModelPhysicSolid<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "fsi" )
        return std::make_shared<ModelPhysicFSI<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "GenericPDE" )
        return std::make_shared<ModelPhysicCoefficientFormPDE<Dim>>( mphysics, modeling, type, name, model );
    else if ( modeling == "GenericPDEs" )
        return std::make_shared<ModelPhysicCoefficientFormPDEs<Dim>>( mphysics, modeling, type, name, model );
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
void
ModelPhysicFluid<Dim>::DynamicViscosity::setup( nl::json const& jarg )
{
    bool isMultifluid = false;

    if ( jarg.is_string() )
    {
        this->setLaw( jarg.template get<std::string>() );
    }
    else if ( jarg.is_object() )
    {
        CHECK( jarg.contains( "law" ) && jarg.at( "law" ).is_string() ) << "invalid law";
        this->setLaw( jarg.at( "law" ).template get<std::string>() );

        if ( jarg.contains( "multifluid" ) )
        {
            auto const& jarg_multifluid =  jarg.at("multifluid");
            if ( jarg_multifluid.is_boolean() )
                isMultifluid = jarg_multifluid.template get<bool>();
            else if ( jarg_multifluid.is_number_unsigned() )
                isMultifluid = jarg_multifluid.template get<int>() > 0;
            else if ( jarg_multifluid.is_string() )
                isMultifluid = boost::lexical_cast<bool>( jarg_multifluid.template get<std::string>() );
        }
    }

    this->setMultifluid( isMultifluid );
}

template <uint16_type Dim>
void
ModelPhysicFluid<Dim>::DynamicViscosity::setMultifluid( bool b )
{
    M_isMultifluid = b; 
    if( M_isMultifluid )
    {
        M_parent->addMaterialPropertyDescription( "heaviside", "heaviside", { std::make_pair(1,1) } );
    }
}

template <uint16_type Dim>
ModelPhysicFluid<Dim>::ModelPhysicFluid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model ),
    M_equation( "Navier-Stokes" ),
    M_navierStokesFormulation( NavierStokesFormulation::Convective ),
    M_gravityForceEnabled( boption(_name="use-gravity-force",_prefix=mphysics.prefix(),_vm=mphysics.clovm()) ),
    M_dynamicViscosity( soption(_name="viscosity.law",_prefix=mphysics.prefix(),_vm=mphysics.clovm()), this ),
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

        // navier stokes formulation
        for ( std::string const& formulationarg : { "formulation" } )
            if ( j_setup.contains( formulationarg ) )
            {
                 auto const& formulationStr = j_setup.at( formulationarg ).get<std::string>();
                 if ( formulationStr == "convective" )
                     M_navierStokesFormulation = NavierStokesFormulation::Convective;
                 else if ( formulationStr == "skew-symmetric" )
                     M_navierStokesFormulation = NavierStokesFormulation::SkewSymmetric;
                 else if ( formulationStr == "conservative" )
                     M_navierStokesFormulation = NavierStokesFormulation::Conservative;
                 else if ( formulationStr == "rotational" )
                     M_navierStokesFormulation = NavierStokesFormulation::Rotational;
                 else if ( formulationStr == "emac" || formulationStr == "EMAC" )
                     M_navierStokesFormulation = NavierStokesFormulation::EMAC;
                 else
                     CHECK( false) << "wrong formulation value : " << formulationStr;
            }

        // gravity force
        if ( j_setup.contains("gravity") )
        {
            auto const& j_setup_gravity = j_setup.at("gravity");
            if ( j_setup_gravity.is_boolean() )
                M_gravityForceEnabled = j_setup_gravity.template get<bool>();
            else if ( j_setup_gravity.is_number_unsigned() )
                M_gravityForceEnabled = j_setup_gravity.template get<int>() > 0;
            else if ( j_setup_gravity.is_string() )
                M_gravityForceEnabled = boost::lexical_cast<bool>( j_setup_gravity.template get<std::string>() );
            else if ( j_setup_gravity.is_object() )
            {
                if ( j_setup_gravity.contains( "enable" ) )
                {
                    auto const& j_setup_gravity_enable =  j_setup_gravity.at("enable");
                    if ( j_setup_gravity_enable.is_boolean() )
                        M_gravityForceEnabled = j_setup_gravity_enable.template get<bool>();
                    else if ( j_setup_gravity_enable.is_number_unsigned() )
                        M_gravityForceEnabled = j_setup_gravity_enable.template get<int>() > 0;
                    else if ( j_setup_gravity_enable.is_string() )
                        M_gravityForceEnabled = boost::lexical_cast<bool>( j_setup_gravity_enable.template get<std::string>() );
                }
                if ( j_setup_gravity.contains( "expr" ) )
                    M_gravityForceExpr.setExpr( j_setup_gravity.at("expr"),mphysics.worldComm(),mphysics.repository().expr() );
            }
        }

        if ( j_setup.contains("viscosity_law") )
        {
            M_dynamicViscosity.setup( j_setup.at("viscosity_law") );
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

    nl::json & pFluid = p["Fluid"];
    pFluid["equation"] = M_equation;
    if ( M_equation == "Navier-Stokes" )
    {
        std::map<NavierStokesFormulation,std::string> mapFormulation = { { NavierStokesFormulation::Convective, "Convective" },
                                                                         { NavierStokesFormulation::SkewSymmetric, "Skew-Symmetric" },
                                                                         { NavierStokesFormulation::Conservative, "Conservative" },
                                                                         { NavierStokesFormulation::Rotational, "Rotational" },
                                                                         { NavierStokesFormulation::EMAC, "EMAC" } };
        pFluid["formulation"] = mapFormulation.at( this->navierStokesFormulation() );
    }
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
    }
    if ( jsonInfo.contains("Fluid") )
    {
        auto const& jsonInfoFluid = jsonInfo.at("Fluid");
        Feel::Table tabInfoFluid;
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfoFluid, jsonInfoFluid, tabInfoProp, { "equation","formulation" } );
        tabInfo->add( "", TabulateInformations::New( tabInfoFluid, tabInfoProp ) );
    }
    if ( jsonInfo.contains("Generic") )
    {
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
        else if ( j_enable.is_number_unsigned() )
            M_isEnabled = j_enable.template get<int>() > 0;
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
ModelPhysicMultifluid<Dim>::ModelPhysicMultifluid( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model ):
    super_type( modeling, type, name, mphysics, model ),
    M_nFluids( ioption( _name="nfluids", _prefix=mphysics.prefix(), _vm=mphysics.clovm() ) )
{
    auto const& j_setup = model.setup();
    if( j_setup.is_object() )
    {
        if( j_setup.contains( "nfluids" ) )
        {
            auto const& j_setup_nfluids = j_setup.at( "nfluids" );
            if( j_setup_nfluids.is_number_unsigned() )
                M_nFluids = j_setup_nfluids.template get<unsigned int>();
            else if( j_setup_nfluids.is_string() )
                M_nFluids = std::stoi( j_setup_nfluids.template get<std::string>() );
        }
        CHECK( M_nFluids >= 2 ) << "invalid number of fluids (" << M_nFluids << ")";

        if ( j_setup.contains( "outer_fluid" ) )
        {
            auto const& j_outerfluid = j_setup.at( "outer_fluid" );
            if ( j_outerfluid.is_string() )
                M_outerFluid = j_outerfluid.get<std::string>();
            else
                CHECK( false ) << "outer_fluid must be a string";
        }
        else
        {
            // default: fluid
            M_outerFluid = "fluid";
        }

        if ( j_setup.contains( "inner_fluids" ) )
        {
            auto const& j_innerfluids = j_setup.at( "inner_fluids" );
            if ( j_innerfluids.is_string() )
                M_innerFluids.insert( j_innerfluids.get<std::string>() );
            else if ( j_innerfluids.is_array() )
            {
                for ( auto const& [j_innerfluidskey,j_innerfluidsval]: j_innerfluids.items() )
                    if ( j_innerfluidsval.is_string() )
                        M_innerFluids.insert( j_innerfluidsval.get<std::string>() );
            }
        }
        else
        {
            // default: fluid%i
            for( uint32_t n = 0; n < M_nFluids - 1; ++n)
            {
                M_innerFluids.insert( fmt::format( "fluid{}", n ) );
            }
        }
    }

    // Set outer and inner fludis as materials used by MultifluidPhysic
    std::set<std::string> materialNames = M_innerFluids;
    materialNames.insert( M_outerFluid );
    this->addMaterialNames( materialNames );

    // Update properties description
    this->updateMaterialPropertyDescription();
}

template <uint16_type Dim>
void
ModelPhysicMultifluid<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicMultifluid<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
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
ModelPhysicMultifluid<Dim>::updateMaterialPropertyDescription()
{
    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    auto const addMaterialPropertyDescriptionFluid = [this, scalarShape, matrixShape]( const std::string & prefix ) {
        this->addMaterialPropertyDescription( prefixvm( prefix, "dynamic-viscosity", "_" ), prefixvm( prefix, "mu", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "turbulent-dynamic-viscosity", "_" ), prefixvm( prefix, "mu_t", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "turbulent-kinetic-energy", "_" ), prefixvm( prefix, "tke", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "consistency-index", "_" ), prefixvm( prefix, "mu_k", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "power-law-index", "_" ), prefixvm( prefix, "mu_power_law_n", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "viscosity-min", "_" ), prefixvm( prefix, "mu_min", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "viscosity-max", "_" ), prefixvm( prefix, "mu_max", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "viscosity-zero-shear", "_" ), prefixvm( prefix, "mu_0", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "viscosity-infinite-shear", "_" ), prefixvm( prefix, "mu_inf", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "carreau-law-lambda", "_" ), prefixvm( prefix, "mu_carreau_law_lambda", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "carreau-law-n", "_" ), prefixvm( prefix, "mu_carreau_law_n", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "carreau-yasuda-law-lambda", "_" ), prefixvm( prefix, "mu_carreau_yasuda_law_lambda", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "carreau-yasuda-law-n", "_" ), prefixvm( prefix, "mu_carreau_yasuda_law_n", "_" ), { scalarShape } );
        this->addMaterialPropertyDescription( prefixvm( prefix, "carreau-yasuda-law-a", "_" ), prefixvm( prefix, "mu_carreau_yasuda_law_a", "_" ), { scalarShape } );
    };

    // Fluids
    for ( uint32_t n = 0; n < this->nFluids(); ++n )
    {
        std::string const fluidKeyword = fmt::format( "fluid{}", n );
        addMaterialPropertyDescriptionFluid( fluidKeyword );
    }
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
ModelPhysicCoefficientFormPDE<Dim>::Infos::Infos( std::string const& name, nl::json const& jarg )
{
    M_equationName = name;

    if ( jarg.contains("name") )
    {
        auto const& j_name = jarg.at("name");
        if ( j_name.is_string() )
            M_equationName = j_name.template get<std::string>();
    }

    if ( jarg.contains("unknown") )
    {
        auto const& j_unknown = jarg.at("unknown");
        CHECK( j_unknown.contains("name") ) << "require to define unknown.name";
        auto const& j_unknown_name = j_unknown.at("name");
        if ( j_unknown_name.is_string() )
            M_unknownName = j_unknown_name.template get<std::string>();
        CHECK( !M_unknownName.empty() ) << "require to define a non empty unknown.name";

        if ( j_unknown.contains( "symbol" ) )
        {
            auto const& j_unknown_symbol = j_unknown.at( "symbol" );
            if ( j_unknown_symbol.is_string() )
                M_unknownSymbol = j_unknown_symbol.template get<std::string>();
        }
        else
             M_unknownSymbol = M_unknownName;

        if ( j_unknown.contains( "basis" ) )
        {
            auto const& j_unknown_basis = j_unknown.at( "basis" );
            if ( j_unknown_basis.is_string() )
                M_unknownBasis = j_unknown_basis.template get<std::string>();
        }
        else
            M_unknownBasis = "Pch1";
    }

    this->initCoefficientProperties();
}

template <uint16_type Dim>
ModelPhysicCoefficientFormPDE<Dim>::Infos::Infos( std::string const& name, std::string const& unknownName, std::string const& unknownSymbol, std::string const& unknownBasis ):
    M_equationName( name ),
    M_unknownName( unknownName ),
    M_unknownSymbol( unknownSymbol ),
    M_unknownBasis( unknownBasis )
{
    this->initCoefficientProperties();
}

template <uint16_type Dim>
void
ModelPhysicCoefficientFormPDE<Dim>::Infos::initCoefficientProperties()
{
    std::string unknownShape;
    if ( M_unknownBasis == "Pch1" || M_unknownBasis == "Pch2" || M_unknownBasis == "Pdh1" )
        unknownShape = "scalar";
    else if ( M_unknownBasis == "Pchv1" || M_unknownBasis == "Pchv2" || M_unknownBasis == "Ned1h0" )
        unknownShape = "vectorial";
    else
        CHECK( false ) << "invalid unknown.basis : " << M_unknownBasis;


    static constexpr uint16_type nDim = Dim;
    shape_dim_type scalarShape = std::make_pair(1,1);
    shape_dim_type vectorialShape = std::make_pair(nDim,1);
    shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    M_coefficientProperties.emplace( Coefficient::convection, std::make_tuple( "beta", shapes_dim_type{ vectorialShape } ) );
    M_coefficientProperties.emplace( Coefficient::diffusion, std::make_tuple( "c", shapes_dim_type{ scalarShape, matrixShape } ) );
    M_coefficientProperties.emplace( Coefficient::reaction, std::make_tuple( "a", shapes_dim_type{ scalarShape } ) );
    M_coefficientProperties.emplace( Coefficient::firstTimeDerivative, std::make_tuple( "d", shapes_dim_type{ scalarShape } ) );
    M_coefficientProperties.emplace( Coefficient::secondTimeDerivative, std::make_tuple( "m", shapes_dim_type{ scalarShape } ) );

    M_coefficientProperties.emplace( Coefficient::conservativeFluxConvection, std::make_tuple( "alpha", shapes_dim_type{ vectorialShape } ) );

    if ( unknownShape == "scalar" )
    {
        M_coefficientProperties.emplace( Coefficient::source, std::make_tuple( "f", shapes_dim_type{ scalarShape } ) );
        M_coefficientProperties.emplace( Coefficient::conservativeFluxSource, std::make_tuple( "gamma", shapes_dim_type{ vectorialShape } ) );
    }
    else if ( unknownShape == "vectorial" )
    {
        M_coefficientProperties.emplace( Coefficient::source, std::make_tuple( "f", shapes_dim_type{ vectorialShape } ) );
        M_coefficientProperties.emplace( Coefficient::conservativeFluxSource, std::make_tuple( "gamma", shapes_dim_type{ matrixShape } ) );
        M_coefficientProperties.emplace( Coefficient::curlCurl, std::make_tuple( "zeta", shapes_dim_type{ scalarShape } ) );
    }
}


template <uint16_type Dim>
ModelPhysicCoefficientFormPDE<Dim>::ModelPhysicCoefficientFormPDE( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    auto const& mgpde = dynamic_cast<ModelGenericPDE<Dim> const&>( mphysics );
    M_infos = mgpde.M_infos;

    auto const& j_setup = model.setup();

    if ( j_setup.contains("coefficients") )
    {
        auto const& j_setup_coeff = j_setup.at("coefficients");
        for ( auto const& [j_setup_coeffkey,j_setup_coeffval] : j_setup_coeff.items() )
        {
            this->addParameter( j_setup_coeffkey, j_setup_coeffval );
        }

        // check coefficients shape
        for ( auto const& [c,prop] : M_infos->coefficientProperties() )
        {
            if ( !this->hasCoefficient( c ) )
                continue;

            auto const& coeff = this->coefficient( c );
            bool shapeOk = false;
            for ( auto const& [dim1,dim2] : std::get<1>( prop ) )
            {
                if ( coeff.hasExpr( dim1,dim2 ) )
                {
                    shapeOk = true;
                    break;
                }
            }
            CHECK( shapeOk ) << "invalid coefficient expr shape : " << std::get<0>( prop );
        }
    }

}

template <uint16_type Dim>
void
ModelPhysicCoefficientFormPDE<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicCoefficientFormPDE<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
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
ModelPhysicCoefficientFormPDEs<Dim>::ModelPhysicCoefficientFormPDEs( ModelPhysics<Dim> const& mphysics, std::string const& modeling, std::string const& type, std::string const& name, ModelModel const& model )
    :
    super_type( modeling, type, name, mphysics, model )
{
    using model_physic_cfpde_type = ModelPhysicCoefficientFormPDE<Dim>;
    auto const& j_setup = model.setup();

    if ( j_setup.is_string() )
    {
        std::string const& eqName = j_setup.template get<std::string>();
        M_pdes.push_back( std::make_tuple( eqName, typename model_physic_cfpde_type::infos_ptrtype{}, std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>{} ) );
    }
    else if ( j_setup.is_object() )
    {
#if 0
        std::string nameEqDefault = fmt::format("equation{}",M_pdes.size());
        //typename ModelPhysicCoefficientFormPDE<Dim>::infos_type infos( nameEqDefault, j_setup );
        auto infos = std::make_shared<cfpde_infos_type>( nameEqDefault, j_setup );
        M_pdes.push_back( std::make_tuple( infos->equationName(), std::move( infos ), std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>{} ) );
#endif
    }
    else if ( j_setup.is_array() )
    {
        for ( auto const& [j_setupkey,j_setupval] : j_setup.items() )
        {
            if ( j_setupval.is_string() )
            {
                std::string const& eqName = j_setupval.template get<std::string>();
                M_pdes.push_back( std::make_tuple( eqName, typename model_physic_cfpde_type::infos_ptrtype{}, std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>{} ) );
            }
            else
            {
#if 0
                std::string nameEqDefault = fmt::format("equation{}",M_pdes.size());
                //typename ModelPhysicCoefficientFormPDE<Dim>::infos_type infos( nameEqDefault, j_setupval );
                auto infos = std::make_shared<cfpde_infos_type>( nameEqDefault, j_setupval );
                M_pdes.push_back( std::make_tuple( infos->equationName(), std::move( infos ), std::shared_ptr<ModelPhysicCoefficientFormPDE<Dim>>{} ) );
#endif
            }
        }
    }

}

template <uint16_type Dim>
void
ModelPhysicCoefficientFormPDEs<Dim>::updateInformationObject( nl::json & p ) const
{
    super_type::updateInformationObject( p["Generic"] );
}
template <uint16_type Dim>
tabulate_informations_ptr_t
ModelPhysicCoefficientFormPDEs<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
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
ModelPhysics<Dim>::PhysicsTreeNode::addChild( std::string const& type, physics_ptrtype p, ModelModels const& models )
{
    std::set<std::string> submodelsName;
    if ( models.hasType( this->physic()->type() ) )
    {
        for ( auto const& [_nameFromModels,_model] : models.models( this->physic()->type() ) )
            if ( _nameFromModels.empty() || _nameFromModels == this->physic()->name() )
            {
                auto itFindSubmodel = _model.submodels().find( type );
                if ( itFindSubmodel != _model.submodels().end() )
                    submodelsName =  itFindSubmodel->second;
            }
    }
    //std::cout << "hahah submodelsName = " << submodelsName << std::endl;


    auto phycisUsedCollection = p->initPhysics( p->keyword(), models, submodelsName );

    // first : add subtree
    std::vector<std::shared_ptr<PhysicsTreeNode>> treeNodeAdded;
    for ( auto phycisUsed : phycisUsedCollection )
    {
        M_children.push_back( PhysicsTreeNode::New(p,phycisUsed) );
        treeNodeAdded.push_back( M_children.back() );

        std::get<2>(M_data)->addSubphysic( phycisUsed );
    }
    // second : up subphysics
    for ( auto treeNode : treeNodeAdded )
    {
        p->updatePhysics( *treeNode, models );
    }

}


template <uint16_type Dim>
void
ModelPhysics<Dim>::PhysicsTreeNode::updateMaterialSupportFromChildren( std::string const& applyType )
{
    if ( M_children.empty() )
        return;

    std::vector<std::set<std::string>> matNamesForEachSubphysic;
    for ( auto & child : M_children )
    {
        //child->updateSubphysicsMaterial();
        matNamesForEachSubphysic.push_back( child->physic()->materialNames() );
    }


    bool isInterfacePhysics = applyType != "intersect";
    //bool isInterfacePhysics = false;
    if ( isInterfacePhysics )
    {

    }
    else // take intersection
    {
        std::set<std::string> intersecMatNames; bool isFirst = true;
        for ( auto const& _matnames : matNamesForEachSubphysic )
        {
            if ( _matnames.empty() )
                continue;
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
        //CHECK( !intersecMatNames.empty() ) << "incompatible materials names";
        if ( !intersecMatNames.empty() )
        {
            auto currentPhysic = std::get<2>( M_data );
            if ( !currentPhysic->materialNames().empty() )
            {
                // TODO
            }
            currentPhysic->setMaterialNames( intersecMatNames );
        }

    }
}



template <uint16_type Dim>
void
ModelPhysics<Dim>::PhysicsTree::updatePhysics( physics_ptrtype mphysics, ModelModels const& models )
{
    auto physicsUsedCollection = mphysics->initPhysics( mphysics->keyword(), models, {} );

    for ( auto physicsUsed : physicsUsedCollection )
    {
        auto treeNode = this->addNode( mphysics, physicsUsed );
        mphysics->updatePhysics( *treeNode, models );
        //treeNode->updateSubphysicsMaterial();
    }
}

template <uint16_type Dim>
std::map<typename ModelPhysics<Dim>::physic_id_type,typename ModelPhysics<Dim>::model_physic_ptrtype>
ModelPhysics<Dim>::PhysicsTree::collectAllPhysics() const
{
    std::map<physic_id_type,model_physic_ptrtype> res;

    auto upCollect = [&res]( PhysicsTreeNode const& nodeTree, const auto& ff ) -> void
                         {
                             auto physicData = nodeTree.physic();
                             auto [it,isInserted] = res.emplace( physicData->id(), physicData );
                             if ( !isInserted )
                                 return;
                             for ( auto const& child : nodeTree.children() )
                                 ff( *child, ff );
                         };
    for ( auto const& nodeTree : this->children() )
        upCollect( *nodeTree, upCollect );

    return res;
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::shared_ptr<ModelPhysics<nDim>> mphysics, ModelModels const& models )
{
    PhysicsTree physicsTree( mphysics );
    physicsTree.updatePhysics( mphysics, models );
    this->addPhysics( physicsTree.collectAllPhysics() );
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::shared_ptr<ModelPhysics<nDim>> mphysics, std::function<void(PhysicsTree &)> lambdaInit )
{
    PhysicsTree physicsTree( mphysics );
    lambdaInit( physicsTree );
    this->addPhysics( physicsTree.collectAllPhysics() );
}


template <uint16_type Dim>
std::vector<std::shared_ptr<ModelPhysic<Dim>>>
ModelPhysics<Dim>::initPhysics( std::string const& type, ModelModels const& models, std::set<std::string> const& modelNameRequired, bool isRequired )
{
    //M_physicType = type;

    std::vector<std::shared_ptr<ModelPhysic<Dim>>> ret;
    if ( models.hasType( type ) )
    {
        for ( auto const& [_nameFromModels,_model] : models.models( type ) )
        {
            std::string _name = _nameFromModels.empty()? type : _nameFromModels; // use name=type if not given
            if ( !modelNameRequired.empty() &&  modelNameRequired.find( _name ) == modelNameRequired.end() )
                continue;
            auto pId = std::make_pair(type, _name);
            if ( !this->hasPhysic( pId ) )
            {
                auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
                M_physics.emplace( pId, _mphysic );
            }
            ret.push_back( this->physic( pId ) );
        }
    }
    else// if ( modelNameRequired.empty() )
    {
        if ( this->physics( type ).empty()  )
        {
            // create default physic
            std::string _name = "default";
            auto _model = ModelModel(type,_name);
            auto _mphysic = ModelPhysic<Dim>::New( *this, M_physicModeling, type, _name, _model );
            auto pId = std::make_pair(type, _name);
            M_physics.emplace( pId, _mphysic );
            ret.push_back( this->physic( pId ) );
        }
        else
        {
            //CHECK( false ) << "TODO take all";
            for ( auto const& [pId,mphysic] : this->physics( type ) )
                ret.push_back( mphysic );
        }
    }
    return ret;
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
ModelPhysics<Dim>::addPhysics( std::map<physic_id_type,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics )
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
template class ModelPhysicCoefficientFormPDE<2>;
template class ModelPhysicCoefficientFormPDE<3>;
template class ModelPhysicCoefficientFormPDEs<2>;
template class ModelPhysicCoefficientFormPDEs<3>;

} // namespace FeelModels
} // namespace Feel
