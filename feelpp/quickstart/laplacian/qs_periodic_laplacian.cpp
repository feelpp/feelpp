// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

#include <feel/feel.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>

#include <sstream>
#include <stdexcept>

using namespace Feel;

namespace
{
std::string periodicSquareGeo( double h )
{
    std::ostringstream geo;
    geo << "Mesh.MshFileVersion = " << FEELPP_GMSH_FORMAT_VERSION << ";\n"
        << "h=" << h << ";\n"
        << "Point(1) = {0,0,0,h};\n"
        << "Point(2) = {1,0,0,h};\n"
        << "Point(3) = {1,1,0,h};\n"
        << "Point(4) = {0,1,0,h};\n"
        << "Line(1) = {1,2};\n"
        << "Line(2) = {2,3};\n"
        << "Line(3) = {3,4};\n"
        << "Line(4) = {4,1};\n"
        << "Line Loop(5) = {1,2,3,4};\n"
        << "Plane Surface(6) = {5};\n"
        << "Physical Surface(\"Omega\") = {6};\n"
        << "Physical Line(\"Bottom\") = {1};\n"
        << "Physical Line(\"Right\") = {2};\n"
        << "Physical Line(\"Top\") = {3};\n"
        << "Physical Line(\"Left\") = {4};\n"
        << "Periodic Curve {2} = {4} Translate {1,0,0};\n"
        << "Periodic Curve {3} = {1} Translate {0,1,0};\n";
    return geo.str();
}

std::string periodicCubeGeo( double h )
{
    std::ostringstream geo;
    geo << "Mesh.MshFileVersion = " << FEELPP_GMSH_FORMAT_VERSION << ";\n"
        << "SetFactory(\"OpenCASCADE\");\n"
        << "h=" << h << ";\n"
        << "Box(1) = {0,0,0,1,1,1};\n"
        << "Mesh.CharacteristicLengthMin = h;\n"
        << "Mesh.CharacteristicLengthMax = h;\n"
        << "eps=1e-6;\n"
        << "sxm[] = Surface In BoundingBox{-eps,-eps,-eps,eps,1+eps,1+eps};\n"
        << "sxp[] = Surface In BoundingBox{1-eps,-eps,-eps,1+eps,1+eps,1+eps};\n"
        << "sym[] = Surface In BoundingBox{-eps,-eps,-eps,1+eps,eps,1+eps};\n"
        << "syp[] = Surface In BoundingBox{-eps,1-eps,-eps,1+eps,1+eps,1+eps};\n"
        << "szm[] = Surface In BoundingBox{-eps,-eps,-eps,1+eps,1+eps,eps};\n"
        << "szp[] = Surface In BoundingBox{-eps,-eps,1-eps,1+eps,1+eps,1+eps};\n"
        << "Physical Volume(\"Omega\") = {1};\n"
        << "Physical Surface(\"XMin\") = {sxm[]};\n"
        << "Physical Surface(\"XMax\") = {sxp[]};\n"
        << "Physical Surface(\"YMin\") = {sym[]};\n"
        << "Physical Surface(\"YMax\") = {syp[]};\n"
        << "Physical Surface(\"ZMin\") = {szm[]};\n"
        << "Physical Surface(\"ZMax\") = {szp[]};\n"
        << "Periodic Surface {sxp[]} = {sxm[]} Translate {1,0,0};\n"
        << "Periodic Surface {syp[]} = {sym[]} Translate {0,1,0};\n"
        << "Periodic Surface {szp[]} = {szm[]} Translate {0,0,1};\n";
    return geo.str();
}

template<int Dim>
std::string periodicGeo( double h )
{
    if constexpr ( Dim == 2 )
        return periodicSquareGeo( h );
    else
        return periodicCubeGeo( h );
}

template<int Dim>
int runPeriodicLaplacian( double h, int order, bool doExport, double assertL2 )
{
    using mesh_type = Mesh<Simplex<Dim, 1>>;
    RuntimeOrder const runtimeOrder( order );

    auto mesh = [&]() {
        Gmsh gmsh;
        gmsh.setDimension( Dim );
        gmsh.setOrder( 1 );
        gmsh.setVersion( FEELPP_GMSH_FORMAT_VERSION, GMSH_FORMAT_ASCII );
        gmsh.setPrefix( Dim == 2 ? "qs-periodic-laplacian-2d" : "qs-periodic-laplacian-3d" );

        std::string filename;
        bool generated_or_modified = false;
        boost::tie( filename, generated_or_modified ) = gmsh.generate( gmsh.prefix(),
                                                                        periodicGeo<Dim>( h ),
                                                                        true );
        Feel::detail::ignore_unused_variable_warning( generated_or_modified );

        return loadGMSHMesh( _mesh = new mesh_type,
                             _filename = filename,
                             _update = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    }();

    CHECK( mesh && mesh->isPeriodic() ) << "Periodic mesh generation/import failed";

    auto Xh = Pch<Dynamic>( mesh, runtimeOrder );
    auto u = Xh->element( "u" );
    auto v = Xh->element( "v" );

    auto uExact = [&]() {
        if constexpr ( Dim == 2 )
            return expr( "cos(2*pi*x)*cos(2*pi*y):x:y" );
        else
            return expr( "cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z):x:y:z" );
    }();
    auto rhs = [&]() {
        if constexpr ( Dim == 2 )
            return expr( "(1+8*pi*pi)*cos(2*pi*x)*cos(2*pi*y):x:y" );
        else
            return expr( "(1+12*pi*pi)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z):x:y:z" );
    }();

    auto a = form2( _trial = Xh, _test = Xh );
    a = integrate( _range = elements( mesh ),
                   _expr = gradt( u ) * trans( grad( v ) ) + idt( u ) * id( v ) );

    auto l = form1( _test = Xh );
    l = integrate( _range = elements( mesh ), _expr = rhs * id( v ) );

    a.solve( _rhs = l, _solution = u );

    auto errL2 = normL2( _range = elements( mesh ), _expr = idv( u ) - uExact );
    if ( assertL2 > 0 )
    {
        CHECK( errL2 <= assertL2 )
            << "L2 error exceeds threshold: ||u-uex||_L2=" << errL2
            << " > qs.assert-l2=" << assertL2;
    }

    if ( Environment::isMasterRank() )
    {
        std::cout << "Periodic Laplacian on periodic mesh\n"
                  << "  dim   = " << Dim << "\n"
                  << "  order = " << order << "\n"
                  << "  nDof  = " << Xh->nDof() << "\n"
                  << "  ||u-uex||_L2 = " << errL2 << "\n";
    }

    if ( doExport )
    {
        auto e = exporter( _mesh = mesh );
        e->addRegions();
        e->add( "u", u );
        auto uex = vf::project( _space = Xh, _range = elements( mesh ), _expr = uExact );
        e->add( "uexact", uex );
        e->save();
    }

    return 0;
}
}

int main( int argc, char** argv )
{
    po::options_description opts( "Periodic Laplacian options" );
    opts.add( feel_options() );
    opts.add_options()
        ( "qs.dim", po::value<int>()->default_value( 2 ), "problem dimension: 2 or 3" )
        ( "qs.order", po::value<int>()->default_value( 1 ), "polynomial order Pk (k>=1)" )
        ( "qs.hsize", po::value<double>()->default_value( 0.15 ), "mesh size" )
        ( "qs.export", po::value<bool>()->default_value( true ), "export results" )
        ( "qs.assert-l2", po::value<double>()->default_value( 5e-2 ),
          "maximum L2 error threshold (<=0 disables assertion)" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = opts,
                     _about = about( _name = "qs_periodic_laplacian",
                                     _author = "Feel++ Consortium",
                                     _email = "feelpp-devel@feelpp.org" ) );

    auto const dim = ioption( "qs.dim" );
    auto const order = ioption( "qs.order" );
    auto const hsize = doption( "qs.hsize" );
    auto const doExport = boption( "qs.export" );
    auto const assertL2 = doption( "qs.assert-l2" );
    CHECK( order >= 1 ) << "Invalid --qs.order value, expected >= 1";

    if ( dim == 2 )
        return runPeriodicLaplacian<2>( hsize, order, doExport, assertL2 );
    if ( dim == 3 )
        return runPeriodicLaplacian<3>( hsize, order, doExport, assertL2 );

    throw std::invalid_argument( "Invalid --qs.dim value, expected 2 or 3" );
}
