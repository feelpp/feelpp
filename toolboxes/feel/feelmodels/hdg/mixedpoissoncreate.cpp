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
                                                worldcomm_ptr_t const& worldComm,
                                                std::string const& subPrefix,
                                                ModelBaseRepository const& modelRep )
    : super_type( prefix, worldComm, subPrefix, modelRep ),
      M_tauOrder(ioption( prefixvm(this->prefix(), "tau_order") )),
      M_tauCst(doption( prefixvm(this->prefix(), "tau_constant") )),
      M_hFace(ioption( prefixvm(this->prefix(), "hface") )),
      M_conductivityKey(soption( prefixvm(this->prefix(), "conductivity_json")) ),
      M_nlConductivityKey(soption( prefixvm(this->prefix(),"conductivityNL_json")) ),
      M_useSC(boption( prefixvm(this->prefix(), "use-sc")) ),
      M_useUserIBC(false)
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedPoisson","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());


    this->setFilenameSaveInfo( prefixvm(this->prefix(),"MixedPoisson.info") );


    if ( this->prefix().empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=this->prefix(), _rebuild=true);

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
                                               worldcomm_ptr_t const& worldComm, std::string const& subPrefix,
                                               ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type> ( prefix,worldComm,subPrefix,modelRep );
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
    toc("mesh");

    tic();
    this->initModel();
    toc("model");

    tic();
    this->initSpaces();
    toc("spaces");

    if(!this->isStationary()){
        tic();
        this->createTimeDiscretization();
        this->initTimeStep();
        toc("timeDiscretization",true);
    }

    tic();
    this->initExporter( meshVisu );
    toc("exporter");
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initModel()
{

    // initialize marker lists for each boundary condition type
    auto itField = modelProperties().boundaryConditions().find( "potential");
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
        itField = modelProperties().boundaryConditions().find( "flux");
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
    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement(faces(M_mesh),[this]( auto const& ewrap ) {
                                                                auto const& e = unwrap_ref( ewrap );
                                                                for( auto exAtMarker : this->M_IBCList)
                                                                {
                                                                    if ( e.hasMarker() && e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                                                                        return true;
                                                                }
                                                                return false; });

    auto face_mesh = createSubmesh( M_mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );


    M_Vh = Pdhv<Order>( M_mesh, true);
    M_Wh = Pdh<Order>( M_mesh, true );
    M_Mh = Pdh<Order>( face_mesh, true );
    // M_Ch = Pch<0>( M_mesh, true );
    M_M0h = Pdh<0>( face_mesh );

    // we need one space per ibc
    std::vector<std::string> ibc_markers(M_integralCondition);
    for( int i = 0; i < M_integralCondition; i++)
    {
        ibc_markers.push_back(M_IBCList[i].marker());
    }

    auto ibc_mesh = createSubmesh( M_mesh, markedfaces(M_mesh, ibc_markers), EXTRACTION_KEEP_MESH_RELATION, 0 );
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

    for( int i = 0; i < M_integralCondition; i++ )
        M_mup.push_back(M_Ch->element("mup"));

    solve::strategy s = M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;

#if 0
    M_A_cst = M_backend->newBlockMatrix(_block=csrGraphBlocks(*M_ps));
    M_A = M_backend->newBlockMatrix(_block=csrGraphBlocks(*M_ps));
    M_F = M_backend->newBlockVector(_block=blockVector(*M_ps), _copy_values=false);
#else
    tic();
    M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *M_backend ); //M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    M_A = makeSharedMatrixCondensed<value_type>(s,  csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *M_backend ); //M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *M_backend, false);//M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    toc("matrixCondensed");
#endif
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
}

} // namespace FeelModels
} // namespace Feel
