#include <feel/feelmodels/hdg/mixedpoisson.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelmesh/complement.hpp>
#include <boost/algorithm/string.hpp>

namespace Feel
{
namespace FeelModels
{

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::MixedPoisson( std::string const& prefix,
                                                MixedPoissonPhysics const& physic,
                                                worldcomm_ptr_t const& worldComm,
                                                std::string const& subPrefix,
                                                ModelBaseRepository const& modelRep )
    : super_type( prefix, MixedPoissonPhysicsMap[physic]["keyword"], worldComm, subPrefix, modelRep ),
      ModelBase( prefix, MixedPoissonPhysicsMap[physic]["keyword"], worldComm, subPrefix, modelRep ),
      M_physic(physic),
      M_potentialKey(MixedPoissonPhysicsMap[physic]["potentialK"]),
      M_fluxKey(MixedPoissonPhysicsMap[physic]["fluxK"]),
      M_tauOrder(ioption( prefixvm(this->prefix(), "tau_order") )),
      M_tauCst(doption( prefixvm(this->prefix(), "tau_constant") )),
      M_hFace(ioption( prefixvm(this->prefix(), "hface") )),
      M_conductivityKey(soption( prefixvm(this->prefix(), "conductivity_json")) ),
      M_nlConductivityKey(soption( prefixvm(this->prefix(),"conductivityNL_json")) ),
      M_useSC(boption( prefixvm(this->prefix(), "use-sc")) ),
      M_useUserIBC(false),
      M_quadError(ioption(prefixvm(this->prefix(), "error-quadrature")) ),
      M_setZeroByInit(boption(prefixvm(this->prefix(), "set-zero-by-init")) )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedPoisson","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());


    if ( this->prefix().empty())
        M_backend = Feel::backend( _rebuild=true);
    else
        M_backend = Feel::backend( _name=this->prefix(), _rebuild=true);

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonSolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedPoisson","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
typename MixedPoisson<Dim,Order, G_Order, E_Order>::self_ptrtype
MixedPoisson<Dim,Order,G_Order, E_Order>::New( std::string const& prefix,
                                               MixedPoissonPhysics const& physic,
                                               worldcomm_ptr_t const& worldComm, std::string const& subPrefix,
                                               ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type> ( prefix,physic,worldComm,subPrefix,modelRep );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void MixedPoisson<Dim,Order,G_Order, E_Order>::setIBCList( std::vector<std::string> markersIBC )
{
    M_useUserIBC = true;
    M_IBCList.clear();
    for( auto const& marker : markersIBC )
    {
        ExpressionStringAtMarker exAtMark(std::make_tuple("expression", marker, std::string(""), std::string(""), std::string("")));
        M_IBCList.push_back(exAtMark);
    }
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::init( mesh_ptrtype mesh, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;
    toc("mesh", FLAGS_v > 0);

    tic();
    this->initModel();
    toc("model", FLAGS_v > 0);

    tic();
    this->initSpaces();
    toc("spaces", FLAGS_v > 0);

    if(!this->isStationary()){
        tic();
        this->createTimeDiscretization();
        this->initTimeStep();
        toc("timeDiscretization", FLAGS_v > 0);
    }

    tic();
    this->initExporter( meshVisu );
    toc("exporter", FLAGS_v > 0);
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initModel()
{
    this->modelProperties().parameters().updateParameterValues();
    M_paramValues = this->modelProperties().parameters().toParameterValues();
    for( auto const& [k,v] : M_paramValues )
    {
        Feel::cout << " - parameter " << k << " : " << v << std::endl;
    }
    this->modelProperties().materials().setParameterValues( M_paramValues );
    //this->modelProperties().boundaryConditions().setParameterValues( paramValues );
    this->modelProperties().postProcess().setParameterValues( M_paramValues );

    // initialize marker lists for each boundary condition type
    auto itField = modelProperties().boundaryConditions().find( M_potentialKey);
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );

        if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann computed from exact pressure:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Robin:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
    }
    if ( !M_useUserIBC )
    {
        M_IBCList.clear();
        itField = modelProperties().boundaryConditions().find( M_fluxKey);
        if ( itField != modelProperties().boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                Feel::cout << "Integral:";
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    if ( M_mesh->hasFaceMarker(marker) )
                        Feel::cout << " " << marker;
                    else
                        Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                    M_IBCList.push_back(exAtMarker);
                }
                Feel::cout << std::endl;
            }

            itType = mapField.find( "InterfaceCondition" );
            if ( itType != mapField.end() )
            {
                Feel::cout << "Interface condition:";
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    if ( M_mesh->hasFaceMarker(marker) )
                        Feel::cout << " " << marker;
                    else
                        Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                }
                Feel::cout << std::endl;
            }
        }
    }

    if ( M_IBCList.empty() )
        M_integralCondition = 0;
    else
        M_integralCondition = M_IBCList.size();

    if ( boost::icontains(modelProperties().models().model().equations(),"picard") )
        M_isPicard = true;
    else
        M_isPicard = false;
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initSpaces()
{
    std::set<std::string> markers;
    for( auto const& m : this->modelProperties().materials() )
    {
        auto const& mat = m.second;
        if( mat.hasPhysics() )
        {
            switch(M_physic)
            {
            case MixedPoissonPhysics::None:
                break;
            case MixedPoissonPhysics::Electric:
                if( !mat.hasPhysics( { "electric","thermo-electric"} ) )
                    continue;
                break;
            case MixedPoissonPhysics::Heat:
                if( !mat.hasPhysics( { "heat","thermo-electric"} ) )
                    continue;
                break;
            }
        }
        for ( std::string const& matmarker : mat.meshMarkers() )
            markers.insert( matmarker );
    }
    M_rangeMeshElements = markedelements(M_mesh, markers);
    M_Vh = Pdhv<Order>( M_mesh, markedelements(M_mesh, markers) );
    M_Wh = Pdh<Order>( M_mesh, markedelements(M_mesh, markers) );
    M_Whp = Pdh<Order+1>( M_mesh, markedelements(M_mesh, markers) );

    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement(faces(support(M_Wh)),[this]( auto const& ewrap ) {
                                                                auto const& e = unwrap_ref( ewrap );
                                                                for( auto exAtMarker : this->M_IBCList)
                                                                {
                                                                    if ( e.hasMarker() && e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                                                                        return true;
                                                                }
                                                                return false; });

    auto face_mesh = createSubmesh( _mesh=M_mesh, _range=complement_integral_bdy, _update=0 );


    M_Mh = Pdh<Order>( face_mesh, true );
    // M_Ch = Pch<0>( M_mesh, true );
    M_M0h = Pdh<0>( face_mesh );

    // we need one space per ibc
    std::vector<std::string> ibc_markers(M_integralCondition);
    for( int i = 0; i < M_integralCondition; i++)
    {
        ibc_markers.push_back(M_IBCList[i].marker());
    }

    auto ibc_mesh = createSubmesh( _mesh=M_mesh, _range=markedfaces(M_mesh, ibc_markers), _update=0 );
    M_Ch = Pch<0>( ibc_mesh, true );

    Feel::cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
               << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
               << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        Feel::cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;

    auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >( M_integralCondition, M_Ch);
    M_ps = std::make_shared<product2_space_type>(product2(ibcSpaces,M_Vh,M_Wh,M_Mh));

    M_up = M_Vh->element( "u" );
    M_pp = M_Wh->element( "p" );
    M_ppp = M_Whp->element( "pp" );

    for( int i = 0; i < M_integralCondition; i++ )
        M_mup.push_back(M_Ch->element("mup"));

    tic();
    this->initMatricesAndVector();
    toc("matrixCondensed", FLAGS_v > 0);
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initMatricesAndVector()
{
    solve::strategy s = M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;
    solve::strategy spp = solve::strategy::local;
    auto pps = product( M_Whp );

    M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *M_backend );
#ifndef USE_SAME_MAT
    M_A = makeSharedMatrixCondensed<value_type>(s,  csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *M_backend );
#endif
    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *M_backend, false);
    M_App = makeSharedMatrixCondensed<value_type>(spp,  csrGraphBlocks(pps, (spp>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), Feel::backend(), true );
    M_Fpp = makeSharedVectorCondensed<value_type>(solve::strategy::local, blockVector(pps), Feel::backend(), false);
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initExporter( mesh_ptrtype meshVisu )
{
    std::string geoExportType="static"; //change_coords_only, change, static
    M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh() ,
                            _name="Export",
                            _geo=geoExportType,
                            _path=this->exporterPath() );

    // point measures
    auto fieldNamesWithSpacePotential = std::make_pair( std::set<std::string>({this->potentialKey()}), this->potentialSpace() );
    auto fieldNamesWithSpaceFlux = std::make_pair( std::set<std::string>({this->fluxKey()}), this->fluxSpace() );
    auto fieldNamesWithSpaces = boost::hana::make_tuple( fieldNamesWithSpacePotential, fieldNamesWithSpaceFlux );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
       M_measurePointsEvaluation->init( evalPoints );
    }

}

} // namespace FeelModels
} // namespace Feel
