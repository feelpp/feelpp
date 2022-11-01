/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-07

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_geomap.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-07
 */
#define BOOST_TEST_MODULE geomap testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelmesh/regiontree.hpp>


using namespace Feel;



double f( node<double>::type const& __n )
{
    return ublas::sum( __n );
}
double fx( node<double>::type const& __n )
{
    return __n[0];
}
template<uint16_type Dim, template<uint16_type,int,uint16_type> class Entity>
struct TestInterp
{
    typedef Entity<Dim, 1,Dim> entity_type;
    typedef Reference<entity_type,Dim,1,Dim> ref_entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptr_type;

    typedef typename mesh_type::gm_type gm_type;
    static const size_type gmc_context_v = vm::POINT|vm::JACOBIAN|vm::HESSIAN;
    typedef typename gm_type::template Context<typename mesh_type::element_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gm_type::Inverse gic_type;
public:
    TestInterp()
        :
        M_mesh()
    {}
    void test( double hsize, std::string version  = FEELPP_GMSH_FORMAT_VERSION )
    {
        M_mesh = mesh_ptr_type( new mesh_type );
        VLOG(1) << "testing Interp with file format version " << version << "\n";
        std::string fname;
        //GmshHypercubeDomain<entity_type::nDim,entity_type::nOrder,Entity> td;
        GmshSimplexDomain td( entity_type::nDim,entity_type::nOrder );
        td.setVersion( version );
        td.setCharacteristicLength( hsize );
        fname = td.generate( entity_type::name().c_str() );
        ImporterGmsh<mesh_type> import( fname );
        import.setVersion( version );
        M_mesh->accept( import );

        auto rangeElement = elements( *M_mesh );
        auto el_it = boost::get<1>( rangeElement );
        auto el_en = boost::get<2>( rangeElement );

        ref_entity_type refelem;
        typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( M_mesh->gm(),
                refelem.points() ) );

        typename mesh_type::Inverse meshinv( M_mesh );
        using size_type = typename mesh_type::size_type;
        /* initialisation of the mesh::inverse data structure */
        meshinv.addPoints( M_mesh->points() );
        meshinv.distribute();

        std::vector<boost::tuple<size_type, uint16_type > > itab;

        //boost::tie( boost::tuples::ignore, el_it, el_en ) = elements( *M_mesh );
        std::cout << "refelem = " << refelem.points() << "\n";
        gmc_ptrtype gmc;

        for ( ; el_it != el_en; ++el_it )
        {
            auto const& meshElt = unwrap_ref( *el_it );
            if ( !gmc )
                gmc = M_mesh->gm()->template context<gmc_context_v>( meshElt, __geopc );
            else
                gmc->template update<gmc_context_v>( meshElt );
            gic_type gic( M_mesh->gm(), meshElt );

            meshinv.pointsInConvex( meshElt.id(), itab );

            for ( int q = 0; q < itab.size(); ++q )
            {
                std::cout << "xref = " << meshinv.referenceCoords().find(boost::get<0>( itab[q] ))->second << "\n";
            }

            for ( int q = 0; q < refelem.points().size2(); ++q )
            {
                std::cout << "gmc xref " << q << " = " << gmc->xRef( q ) << "\n";
                std::cout << "is in gmc? = " << gmc->geometricMapping()->isIn( gmc->xRef( q ) ) << "\n";

                gic.setXReal( gmc->xReal( q ) );

                typename ref_entity_type::points_type pts( Dim, 1 );
                ublas::column( pts, 0 ) = gic.xRef();
                std::cout << "gic xref " << q << " = " << gic.xRef() << "\n";
                std::cout << "is in gic? = " << gic.geometricMapping()->isIn( gic.xRef() ) << "\n";
            }

            FEELPP_ASSERT( gic.isIn() )
            ( refelem.points() )( gmc->xReal() )
            ( meshElt.id() ).error( "invalid geometric transformation inversion" );
        }

        VLOG(1) << "testing Interp with file format version " << version << " done\n";

    }
private:
    mesh_ptr_type M_mesh;
};


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_geomap )

BOOST_AUTO_TEST_CASE( test_geomap )
{
    TestInterp<2,Simplex> test_interp;
    test_interp.test( doption(_name="gmsh.hsize"), FEELPP_GMSH_FORMAT_VERSION );
}
BOOST_AUTO_TEST_SUITE_END()
