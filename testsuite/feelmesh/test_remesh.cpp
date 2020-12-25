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
BOOST_AUTO_TEST_SUITE_END()
