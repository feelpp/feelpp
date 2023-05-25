/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes_registered_type.hpp>

#include <feel/feelfilters/loadmesh.hpp>

template <int nDim,int nOrderGeo>
void
runApplicationCoefficientFormPDEs()
{
    using namespace Feel;

    using model_type = FeelModels::coefficient_form_PDEs_t< Simplex<nDim,nOrderGeo> >;
    using mesh_type = typename model_type::mesh_type;

    auto mesh = loadMesh(_mesh=new mesh_type); //, _filename="");
    auto Vh = Pch<2>( mesh );
    auto v = Vh->element();
    v.on(_range=elements(mesh), _expr=2*Px()*Py());
    auto Wh = Pchv<1>( mesh );
    auto w = Wh->element();
    w.on(_range=elements(mesh), _expr=5*P());

    auto cfpdes = std::make_shared<model_type>( "cfpdes" );
    nl::json jsonModelsSetup = nl::json::parse(R"(
    {
        "cfpdes":{
            "equations":"heat"
        },
        "heat":{
            "setup":{
                "unknown":{
                    "basis":"Pch1",
                    "name":"temperature",
                    "symbol":"T"
                },
                "coefficients":{
                    "c":"1"
                }
            }
        }
    }
)");

    nl::json jsonPostProcessSetup = nl::json::parse(R"(
    {
        "cfpdes":
        {
            "Exports":
            {
                "fields":["heat.temperature"],
                "expr":
                {
                    "customA":"meshes_cfpdes_fields_customA:meshes_cfpdes_fields_customA",
                    "customB":"{meshes_cfpdes_fields_customB_0,meshes_cfpdes_fields_customB_1}:meshes_cfpdes_fields_customB_0:meshes_cfpdes_fields_customB_1",
                    "customC":"meshes_cfpdes_fields_customC:meshes_cfpdes_fields_customC",
                    "customD":{ "expr":"meshes_cfpdes_fields_customD:meshes_cfpdes_fields_customD", "representation":"element" }
                }
            }
        }
    }
)");

    nl::json jsonToolboxSetup= {
        { "Models", jsonModelsSetup },
        { "PostProcess", jsonPostProcessSetup }
    };

    cfpdes->setModelProperties( jsonToolboxSetup );
    cfpdes->setMesh( mesh );
    cfpdes->modelMesh().template updateField<mesh_type>( "customA", v, "Pch2" );
    cfpdes->modelMesh().template updateField<mesh_type>( "customB", w, "Pchv1" );
    cfpdes->modelMesh().template updateField<mesh_type>( "customC", cst(3.14)*Px(), "Pch1" );
    cfpdes->modelMesh().template updateField<mesh_type>( "customD", cst(2*3.14)*idv(v), "Pdh0" );
    cfpdes->modelMesh().template updateField<mesh_type>( "customD", Px()*Py(), boundaryelements(mesh), "Pdh0" );
    cfpdes->init();
    cfpdes->printAndSaveInfo();
    cfpdes->exportResults();


    // TODO compute measures and check results!
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description cfpdesoptions( "coefficient-form-pdes options" );
    cfpdesoptions.add( toolboxes_options( "coefficient-form-pdes", "cfpdes" ) );

	Environment env( _argc=argc, _argv=argv,
                     _desc=cfpdesoptions,
                     _about=about(_name="toolboxes_cfpdes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runApplicationCoefficientFormPDEs<2,1>();
    return 0;
}
