/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

   Author(s): Thomas Lantz
   Date: 2015-04-27

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
   \file test_multiscaleimage.cpp
   \author Thomas Lantz
   \date 2015-04-27
 */

#define BOOST_TEST_MODULE test_meshstructured
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/meshstructured.hpp>
//#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
//#include <feel/feeldiscr/pdh.hpp>
//#include <feel/feelvf/msi.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
//#include <feel/feelvf/projectors.hpp>
//#include <feel/feelpoly/multiscalequadrature.hpp>
//#include <feel/feelvf/ginac.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelfilters/hbf.hpp>
// #if defined(FEELPP_HAS_FFTW)
// #include <fftw3.h>
// #endif


using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS


#if 0
void runImageHbf  (std::string hbf1, std::string hbf2, std::string solution, double pixelsize) 
{
    auto x = readHBF( hbf1 );
    auto y = readHBF( hbf2 );
    auto z = readHBF( solution );

    auto mesh = std::make_shared<MeshStructured>( x.rows(),
                                                  x.cols(),
                                                  pixelsize,
                                                  cx,cy,
                                                  /*NULL,*/
                                                  Environment::worldCommPtr());
    mesh->components().reset();
    mesh->components().set( size_type(MESH_NO_UPDATE_MEASURES) );
    mesh->updateForUse();

    auto Vh = Pch<1> ( mesh ) ;
    auto u = Vh->element () ;
    auto v = Vh->element () ;
    auto px = Vh->element () ;
    auto py = Vh->element () ;
    Hbf2FeelppStruc h2f( nx, ny, Vh );
    Hbf2FeelppStruc h2fS( z.cols(), z.rows(), Vh );
    px = h2f( x);
    py = h2f( y );
    auto s = h2fS ( z);

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=grad(v)*vec(idv(px),idv(py)));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements(mesh),
                   //_expr=gradt(u)*trans(grad(v)));
                   _expr=gradt(u)*trans(grad(v)));

    a.solve(_rhs=l,_solution=u) ;

    std::cout << "L2 error norm : " << normL2( _range=elements(mesh), _expr=idv(u)-idv(s) ) << "\n";
}
#endif


template <typename GeoShapeType>
auto createStructuredMesh( nl::json const& jsonDesc )
{
    auto mesh = std::make_shared<MeshStructured<GeoShapeType>>( jsonDesc );
    mesh->components().reset();
    mesh->components().set( size_type( MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();
    Feel::cout <<"  mesh->numGlobalElements() " << mesh->numGlobalElements() << std::endl;
    Feel::cout <<"  mesh->numGlobalPoints()   " << mesh->numGlobalPoints()   << std::endl;
    return mesh;
}


template </*int Dim,int PolyOrder=2,*/typename MeshStructuredType >
void runLaplacian( std::shared_ptr<MeshStructuredType> mesh, std::string const& exportName )
{
    static const int Dim = MeshStructuredType::nDim;
    std::string u_exact_str = Dim==1? "3*x^2:x" : (Dim==2? "3*x^2+5*y^2:x:y" : "3*x^2+5*y^2+9*z^2:x:y:z");
    auto u_exact = expr( u_exact_str );
    auto f_str = Dim==1? "-6" : (Dim==2? "-6-10" : "-6-10-18");
    auto f = expr( f_str );

    auto Vh = Pch<2>(mesh/*,true*/);
    // std::cout << "Rank - nDof \t" << Environment::worldComm().localRank() << "\t" <<  Vh->nLocalDof() << std::endl;
    // std::cout << "Rank - Ghost Local Dof :" << Environment::worldComm().localRank() << "\t"<< Vh->dof()->nLocalGhosts() << std::endl;
    auto u = Vh->element();
    auto v = Vh->element();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=f*id(v));
    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements(mesh),
                   _expr=inner(gradt(u),grad(v)) );

    a += on( _range = boundaryfaces(mesh), _rhs = l, _element = u, _expr = u_exact/*cst(0)*/ );

    a.solve(_rhs=l,_solution=u,_rebuild=true) ;

    auto e = exporter( _mesh=mesh,_name= exportName+"_"+MeshStructuredType::shape_type::name() );
    e->addRegions();
    e->add( "u", u );
    e->add( "u_exact", u_exact );
    e->save();

    double err = normL2( _range=elements(mesh), _expr=idv(u)-u_exact );
    BOOST_CHECK_SMALL( err,1e-5 );
}




BOOST_AUTO_TEST_SUITE( meshstructured )

BOOST_AUTO_TEST_CASE( test_json )
{
    nl::json jsonDesc1d = {
        { "Geometry", {
                { "shape", "hyperrectangle" },
                { "domain", {
                        { "origin", { -1.0 } },
                        { "length", { 12 } }
                    } }
            }
        },
        { "Discretisation", {
                { "n_points", 42 }
            }
        }
    };
    nl::json jsonDesc2d = {
        { "Geometry", {
                { "shape", "hyperrectangle" },
                { "domain", {
                        { "origin", { -1.0,3.2 } },
                        { "length", { 12,45.5 } }
                    } }
            }
        },
        { "Discretisation", {
                { "n_points", { 42,18 } }
            }
        }
    };
    nl::json jsonDesc3d = {
        { "Geometry", {
                { "shape", "hyperrectangle" },
                { "domain", {
                        { "origin", { -1.0,3.2,5.5 } },
                        { "length", { 12,45.5,26.2 } }
                    } }
            }
        },
        { "Discretisation", {
                { "n_points", { 23,12,7 } }
            }
        }
    };

    auto mesh1d = createStructuredMesh<Hypercube<1>>(jsonDesc1d);
    BOOST_CHECK_CLOSE( mesh1d->measure(), 12, 1e-7 );
    runLaplacian( mesh1d, "json" );

    auto mesh2d = createStructuredMesh<Hypercube<2>>(jsonDesc2d);
    BOOST_CHECK_CLOSE( mesh2d->measure(), 12*45.5, 1e-7 );
    runLaplacian( mesh2d, "json" );

    auto mesh3d = createStructuredMesh<Hypercube<3>>(jsonDesc3d);
    BOOST_CHECK_CLOSE( mesh3d->measure(), 12*45.5*26.2, 1e-7 );
    runLaplacian( mesh3d, "json" );
}



BOOST_AUTO_TEST_CASE( test_usercoords )
{
    eigen_vector_type<2,double> p0,p1,p2,p3;
    p0[0] = 0; p0[1] = 0;
    p1[0] = 5; p1[1] = 1;
    p2[0] = 4; p2[1] = 4;
    p3[0] = 1; p3[1] = 2;

    double l0 = (p0-p1).norm();
    double l1 = (p1-p2).norm();
    double l2 = (p2-p3).norm();
    double l3 = (p3-p0).norm();

    eigen_vector_type<2,double> dir0,dir1,dir2,dir3;
    dir0 = (p1-p0).normalized();
    dir1 = (p2-p1).normalized();
    dir2 = (p2-p3).normalized();
    dir3 = (p3-p0).normalized();

    int nx = 15;
    int ny = 10;

    double h0 = l0/(nx-1);
    double h1 = l1/(ny-1);
    double h2 = l2/(nx-1);
    double h3 = l3/(ny-1);

    holo3_image<float> coords_x(nx,ny);
    holo3_image<float> coords_y(nx,ny);
    eigen_vector_type<2,double> q0,q1,q2,q3;
    for (int i=0;i<nx;++i)
    {
        q0 = p0 + i*h0*dir0;
        q2 = p3 + i*h2*dir2;
        for (int j=0;j<ny;++j)
        {
            q1 = p1 + j*h1*dir1;
            q3 = p0 + j*h3*dir3;
            // compute intersection of 2 segment
            // see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
            double t = ( (q0(0)-q1(0))*(q1(1)-q3(1)) - (q0(1)-q1(1))*(q1(0)-q3(0)) )/( (q0(0)-q2(0))*(q1(1)-q3(1))-(q0(1)-q2(1))*(q1(0)-q3(0)) );
            auto qij = q0+t*(q2-q0);
            coords_x(i,j) = qij(0);
            coords_y(i,j) = qij(1);
        }
    }


    using geoshape_type = Hypercube<2>;
    auto mesh = std::make_shared<MeshStructured<geoshape_type>>( coords_x, coords_y );
    mesh->components().reset();
    mesh->components().set( size_type( MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    BOOST_CHECK_CLOSE( mesh->measure(), 10, 1e-5 );
    runLaplacian( mesh, "usercoords" );
}

BOOST_AUTO_TEST_SUITE_END()
