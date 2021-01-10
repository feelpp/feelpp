/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Feb 2020

 Copyright (C) 2020 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#define BOOST_TEST_MODULE test mmg

#include <feel/feelcore/testsuite.hpp>
#include <boost/hana/integral_constant.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/remesh.hpp>

using namespace Feel;
namespace hana =  boost::hana;
using namespace hana::literals;
inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_Mmg options" );
    desc_options.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "remesh.verbose", po::value<int>()->default_value( 1 ), "mesh size" )
        ( "remesh.hmin", po::value<double>(), "Minimal mesh size " )
        ( "remesh.hmax", po::value<double>(), "Maximal mesh size " )
        ( "remesh.hsiz", po::value<double>(), "Constant mesh size " )
        
        ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "test_remesh" ,
                     "test_remesh" ,
                     "0.1",
                     "test Mmg",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2020 Universite de Strasbourg" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )


BOOST_AUTO_TEST_SUITE( mmg )

//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<
    std::pair<boost::mpl::int_<2>,boost::mpl::int_<2>>,
    //std::pair<boost::mpl::int_<2>,boost::mpl::int_<3>>,
    std::pair<boost::mpl::int_<3>,boost::mpl::int_<3>> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( testnD, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;

    typedef Mesh<Simplex<T::first_type::value,1,T::second_type::value>> mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::first_type::value == 2 && T::second_type::value == 2 ) 
                           {
                               return unitSquare();
                           }
                           else if constexpr ( T::first_type::value == 2 && T::second_type::value == 3 )
                           {
                               auto m = unitCube();
                               auto mS = createSubmesh( _mesh=m, _range=boundaryfaces( m ) );
                               return mS;
                           }
                           else if constexpr ( T::first_type::value == 3 )
                           {
                               return unitCube();
                           }
                       };
    auto m = create_mesh();
    auto Xh = Pch<1>( m );
    auto met = Xh->element();
    met.on( _range=elements(m), _expr=expr(soption("functions.s")) );
    
    auto r =  remesher( m );
    
    r.setMetric( met );
    auto out = r.execute();
    out->updateForUse();
    auto Vh = Pch<1>( out );
    auto interp = I( _domainSpace=Xh, _imageSpace=Vh );
    auto eout = exporter( _mesh=out, _name="remeshed" );
    auto v= interp.operator()( met );
    auto w= Vh->element( expr(soption("functions.s") ));
    double err = normL2(_range=elements(out), _expr=idv(v)-idv(w));
    BOOST_MESSAGE( "L2 error norm: " << err );
    eout->addRegions();
    eout->add( "metricout", v );
    eout->save();

    auto ein = exporter( _mesh=m, _name="initial" );
    ein->addRegions();
    ein->add( "metric", met );
    ein->save();

}


// testcase with several successive remeshing steps on a square/cube
BOOST_AUTO_TEST_CASE_TEMPLATE( testMultipleRemesh, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;

    typedef Mesh<Simplex<T::first_type::value,1,T::second_type::value>> mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::first_type::value == 2 && T::second_type::value == 2 ) 
                           {
                               return unitSquare();
                           }
                           else if constexpr ( T::first_type::value == 2 && T::second_type::value == 3 )
                           {
                               auto m = unitCube();
                               auto mS = createSubmesh( _mesh=m, _range=boundaryfaces( m ) );
                               return mS;
                           }
                           else if constexpr ( T::first_type::value == 3 )
                           {
                               return unitCube();
                           }
                       };
    auto m = create_mesh();
    auto out = create_mesh();
    auto ein = exporter( _mesh=m, _name="initial",_geo="change");
    auto eout = exporter( _mesh=out, _name="remeshed",_geo="change" );
    std::string metric[4] = { "(abs(x)+1):x",
                              "x*x:x",
                              "exp(x):x",
                              "(x*x+y*y+1):x:y"};
    int ind =0;
    for(int i=0;i<20;i++)
    {
        
        auto Xh = Pch<1>( m );
        auto met = Xh->element();
        met.on( _range=elements(m), _expr=cst(meshSize)*expr(metric[ind]) );
        auto r =  remesher( m );
    
        r.setMetric( met );
        auto out = r.execute();
        out->updateForUse();
        auto Vh = Pch<1>( out );
        auto interp = I( _domainSpace=Xh, _imageSpace=Vh );
        auto v= interp.operator()( met );
        auto w= Vh->element( cst(meshSize)*expr(metric[ind] ));
        double err = normL2(_range=elements(out), _expr=idv(v)-idv(w));
        BOOST_MESSAGE( "L2 error norm: " << err );
        std::cout << "L2 error norm: " << err << std::endl;
        eout->step(i)->setMesh(out);
        eout->step(i)->add( "metricout", v );
        eout->save();
        
        ein->step(i)->setMesh(m);
        ein->step(i)->add( "metric", met );
        ein->save();
        m = out;
        ind++;
        ind = ind%4; std::cout << "index = " << ind <<std::endl;
    }
}
// test checking vector interpolation after multiple remeshing
BOOST_AUTO_TEST_CASE_TEMPLATE( testVecInterp, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;

    typedef Mesh<Simplex<T::first_type::value,1,T::second_type::value>> mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::first_type::value == 2 && T::second_type::value == 2 ) 
                           {
                               return unitSquare();
                           }
                           else if constexpr ( T::first_type::value == 2 && T::second_type::value == 3 )
                           {
                               auto m = unitCube();
                               auto mS = createSubmesh( _mesh=m, _range=boundaryfaces( m ) );
                               return mS;
                           }
                           else if constexpr ( T::first_type::value == 3 )
                           {
                              // auto mesh = loadMesh(_mesh=new Mesh<Simplex<3,1>>( "mesh" ),_filename="cube_cubeHole.geo");
                               //return mesh;
                               return unitCube();
                           }
                       };
    auto m = create_mesh();
    auto out = create_mesh();
    auto ein = exporter( _mesh=m, _name="initial",_geo="change");
    auto eout = exporter( _mesh=out, _name="remeshed",_geo="change" );
    std::string metric[4] = { "(abs(x)+1):x",
                              "x*x:x",
                              "exp(x):x",
                              "(x*x+y*y+1):x:y"};
    int ind =0;
    for(int i=0;i<20;i++)
    {
        
        auto Xh = Pch<1>( m );
        auto met = Xh->element();
        met.on( _range=elements(m), _expr=cst(meshSize)*expr(metric[ind]) );
        auto Th = Pchv<1>(m);
        auto fieldToInterp = Th->element();
        fieldToInterp.on(_range=elements(m), _expr=P()*cst(meshSize)*expr(metric[ind]) );
        auto r =  remesher( m);
    
        r.setMetric( met );
        out = r.execute();
        out->updateForUse();
        auto Vh = Pchv<1>( out );
        auto interp = I( _domainSpace=Th, _imageSpace=Vh );
        auto v= interp.operator()( fieldToInterp );
        auto w= Vh->element( P()*cst(meshSize)*expr(metric[ind] ));
        double err = normL2(_range=elements(out), _expr=idv(v)-idv(w));
        BOOST_MESSAGE( "L2 error norm: " << err );
        std::cout << "L2 error norm: " << err << std::endl;
        eout->step(i)->setMesh(out);
        eout->step(i)->add( "interpField", v );
        eout->save();
        
        ein->step(i)->setMesh(m);
        ein->step(i)->add( "fieldToInterp", fieldToInterp );
        ein->save();
        
        m = out;
        ind++;
        ind = ind%4; std::cout << ind<<std::endl;
    }
     
    
} 

// testcase with hole and fixed discretization on hole boundary
BOOST_AUTO_TEST_CASE_TEMPLATE( testRemeshMarker, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;

    typedef Mesh<Simplex<T::first_type::value,1,T::second_type::value>> mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = 0.05;//doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::first_type::value == 2 && T::second_type::value == 2 ) 
                           {
                               auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>,_filename="$cfgdir/square_squareHole.geo");
                               std::cout << mesh->hasMarkers({"BdryFixedDiscr"}) << std::endl;
                               return mesh;
                           }
                           else if constexpr ( T::first_type::value == 2 && T::second_type::value == 3 )
                           {
                               auto m = unitCube();
                               auto mS = createSubmesh( _mesh=m, _range=boundaryfaces( m ) );
                               return mS;
                           }
                           else if constexpr ( T::first_type::value == 3 )
                           {
                               auto mesh = loadMesh(_mesh=new Mesh<Simplex<3,1>>( "mesh" ),_filename="cube_cubeHole.geo");
                               return mesh;
                           }
                       };
    auto m = create_mesh();
    auto out = create_mesh();
    auto ein = exporter( _mesh=m, _name="initial",_geo="change");
    auto eout = exporter( _mesh=out, _name="remeshed",_geo="change" );
    std::string metric[4] = { "(abs(x)+1):x",
                              "x*x:x",
                              "exp(x):x",
                              "(x*x+y*y+1):x:y"};
    int ind =0;
    for(int i=0;i<20;i++)
    {
        
        auto Xh = Pch<1>( m );
        auto met = Xh->element();
        met.on( _range=elements(m), _expr=cst(meshSize)*expr(metric[ind]) );
        auto r =  remesher( m,"",m->markerName("BdryFixedDiscr"));
    
        r.setMetric( met );
        out = r.execute();
        out->updateForUse();
        auto Vh = Pch<1>( out );
        auto interp = I( _domainSpace=Xh, _imageSpace=Vh );
        auto v= interp.operator()( met );
        auto w= Vh->element( cst(meshSize)*expr(metric[ind] ));
        double err = normL2(_range=elements(out), _expr=idv(v)-idv(w));
        BOOST_MESSAGE( "L2 error norm: " << err );
        std::cout << "L2 error norm: " << err << std::endl;
        eout->step(i)->setMesh(out);
        eout->step(i)->add( "metricout", v );
        eout->save();
        
        ein->step(i)->setMesh(m);
        ein->step(i)->add( "metric", met );
        ein->save();
        
        m = out;
        ind++;
        ind = ind%4; std::cout << ind<<std::endl;
    }
     
    
}  

// testcase with two materials and fixed discretization on inner surface and elements
BOOST_AUTO_TEST_CASE_TEMPLATE( testRemeshTwoMat, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;

    typedef Mesh<Simplex<T::first_type::value,1,T::second_type::value>> mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = 0.05;//doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::first_type::value == 2 && T::second_type::value == 2 ) 
                           {
                               auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>,_filename="$cfgdir/square_twoMaterials.geo");
                               std::cout << mesh->hasMarkers({"BdryFixedDiscr"}) << std::endl;
                               return mesh;
                           }
                           else if constexpr ( T::first_type::value == 2 && T::second_type::value == 3 )
                           {
                               auto m = unitCube();
                               auto mS = createSubmesh( _mesh=m, _range=boundaryfaces( m ) );
                               return mS;
                           }
                           else if constexpr ( T::first_type::value == 3 )
                           {
                               auto mesh = loadMesh(_mesh=new Mesh<Simplex<3,1>>( "mesh" ),_filename="cube_twoMaterials.geo");
                               return mesh;
                           }
                       };
    auto m = create_mesh();
    auto out = create_mesh();
    auto ein = exporter( _mesh=m, _name="initial",_geo="change");
    auto eout = exporter( _mesh=out, _name="remeshed",_geo="change" );
    std::string metric[4] = { "(abs(x)+1):x",
                              "x*x:x",
                              "exp(x):x",
                              "(x*x+y*y+1):x:y"};
    int ind =0;
    for(int i=0;i<20;i++)
    {
        
        auto Xh = Pch<1>( m );
        auto met = Xh->element();
        met.on( _range=elements(m), _expr=cst(meshSize)*expr(metric[ind]) );
        auto r =  remesher( m,m->markerName("MatTwo"),m->markerName("BdryFixedDiscr"));
    
        r.setMetric( met );
        out = r.execute();
        out->updateForUse();
        auto Vh = Pch<1>( out );
        auto interp = I( _domainSpace=Xh, _imageSpace=Vh );
        auto v= interp.operator()( met );
        auto w= Vh->element( cst(meshSize)*expr(metric[ind] ));
        double err = normL2(_range=elements(out), _expr=idv(v)-idv(w));
        BOOST_MESSAGE( "L2 error norm: " << err );
        std::cout << "L2 error norm: " << err << std::endl;
        eout->step(i)->setMesh(out);
        eout->step(i)->add( "metricout", v );
        eout->save();
        
        ein->step(i)->setMesh(m);
        ein->step(i)->add( "metric", met );
        ein->save();
        
        m = out;
        ind++;
        ind = ind%4; std::cout << ind<<std::endl;
    }
     
    
}  

BOOST_AUTO_TEST_SUITE_END()
