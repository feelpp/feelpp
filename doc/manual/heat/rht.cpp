/**
 * @file rht.cpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief Radiative heat transfer 
 * @version 0.1
 * @date 2022-07-22
 * 
 * @copyright Copyright (c) 2022 Université de Strasbourg
 * 
 */
#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel {


    
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "rht options" );
    options.add_options()
    // mesh parameters
    ( "specs", Feel::po::value<std::string>(),
      "json spec file for rht" );

    return options.add( Feel::feel_options() );
}

}  // namespace Feel


int main(int argc, char**argv)
{

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="rht",
                                  _author="Feel++ Consortium",
                                  _email="feelpp@cemosis.fr"));
    auto jsonfile = removeComments(readFromFile(Environment::expand(soption("specs"))));
    std::istringstream istr( jsonfile );
    json specs = json::parse( istr );
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>, _filename=specs["/Meshes/heat/Import/filename"_json_pointer].get<std::string>());
    auto Xh = Pch<2>(mesh);
    auto v = Xh->element();
    auto a = form2( _test=Xh, _trial=Xh );
    auto l = form1( _test=Xh );

    for ( auto [key,material] : specs["/Models/heat/materials"_json_pointer].items() )
    { 
        LOG(INFO) << fmt::format("material {}", material );
        std::string mat = fmt::format("/Materials/{}/k",material.get<std::string>());
        auto k = specs[nl::json::json_pointer(mat)].get<std::string>();
        a += integrate( _range = markedelements( mesh, material.get<std::string>() ), _expr = expr( k ) * gradt( v ) * trans( grad( v ) ) );
    }

    // BC Neumann
    for ( auto& [bc,value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items() )
    { 
        LOG(INFO) << fmt::format("flux {}: {}", bc, value.dump() );
        auto flux = value["expr"].get<std::string>();
        l += integrate( _range = markedfaces( mesh, bc ), _expr = expr(flux) * id( v ) );
    }

    // BC Robin
    for ( auto& [bc,value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
        auto h = value["h"].get<std::string>();
        auto Text = value["Text"].get<std::string>();
        a += integrate( _range = markedfaces( mesh, bc ), _expr = expr(h)*id(v)*idt(v));
        l += integrate( _range = markedfaces( mesh, bc ), _expr = expr( h ) * expr( Text ) * id( v ) );
    }


    // Solve
    a.solve( _rhs=l, _solution=v );


    // compute outputs
    auto m = mean(_range=elements(mesh),_expr=idv(v));
    auto m_root = mean( _range = markedfaces( mesh, "Gamma_root" ), _expr = idv( v ) );
    std::cout << fmt::format( "- mean value: {}", m  ) << std::endl;
    std::cout << fmt::format( "-  min value: {}", v.min()  ) << std::endl;
    std::cout << fmt::format( "-  max value: {}", v.max() ) << std::endl;
    std::cout << fmt::format( "-  max deviation: {}", v.max() - v.min() ) << std::endl;
    std::cout << fmt::format( "-  mean root: {}", m_root ) << std::endl;
    // Export
    auto e = exporter( _mesh = mesh );
    e->addRegions();
    e->add( "T", v );
    e->save();

    return 0;
}
                                