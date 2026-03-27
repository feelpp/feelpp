//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include <iomanip>
#include <sstream>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

namespace
{

po::options_description
makeOptions()
{
    auto opts = feel_options();
    opts.add_options()
        ( "sigma", po::value<double>()->default_value( 58e3 ), "electrical conductivity" )
        ( "current", po::value<double>()->default_value( 800.0 ), "prescribed total outward current through the IBC boundary" )
        ( "potential", po::value<double>()->default_value( 0.03125 ), "Dirichlet electric potential" )
        ( "radius-inner", po::value<double>()->default_value( 1.0 ), "inner radius of the quarter torus/annular sector" )
        ( "radius-outer", po::value<double>()->default_value( 2.0 ), "outer radius of the quarter torus/annular sector" )
        ( "section-height", po::value<double>()->default_value( 1.0 ), "cross-section height in 3D and effective thickness in 2D" )
        ( "angle", po::value<double>()->default_value( M_PI/2.0 ), "sector angle in radians" )
        ( "marker.dirichlet", po::value<std::string>()->default_value( "Dirichlet" ), "Dirichlet marker name" )
        ( "marker.ibc", po::value<std::string>()->default_value( "Ibc" ), "integral boundary condition marker name" )
        ( "marker.neumann", po::value<std::string>()->default_value( "Neumann" ), "homogeneous Neumann marker name" )
        ( "no-solve", po::value<bool>()->default_value( false ), "assemble only, do not solve" );
    return opts;
}

AboutData
makeAbout()
{
    AboutData about( "qs_laplacian_ibc",
                     "qs_laplacian_ibc",
                     "0.1",
                     "Quickstart cG Laplacian with an integral current constraint",
                     AboutData::License_GPL,
                     "Copyright (c) 2026 Feel++ Consortium" );
    about.addAuthor( "Feel++ Consortium", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}

template <int Dim>
std::string
makeExactCurrentDensityString( double coeff )
{
    std::ostringstream os;
    os << std::setprecision( 17 );
    if constexpr ( Dim == 2 )
    {
        os << "{"
           << -coeff << "*y/(x*x+y*y),"
           << coeff << "*x/(x*x+y*y)"
           << "}";
    }
    else
    {
        os << "{"
           << -coeff << "*y/(x*x+y*y),"
           << coeff << "*x/(x*x+y*y),0"
           << "}";
    }
    return os.str();
}

template <int Dim>
std::string
makeCoordinateMapping()
{
    if constexpr ( Dim == 2 )
        return ":x:y";
    else
        return ":x:y:z";
}

template <int Dim, int Order>
int
cg_laplacian_ibc()
{
    using mesh_type = Mesh<Simplex<Dim,1>>;

    auto mesh = loadMesh( _mesh = new mesh_type );

    auto const dirichletMarker = soption( "marker.dirichlet" );
    auto const ibcMarker = soption( "marker.ibc" );
    auto const neumannMarker = soption( "marker.neumann" );

    if ( nelements( markedfaces( mesh, dirichletMarker ), true ) == 0 )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "Missing Dirichlet marker '" << dirichletMarker << "'\n";
        return 1;
    }
    if ( nelements( markedfaces( mesh, ibcMarker ), true ) == 0 )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "Missing IBC marker '" << ibcMarker << "'\n";
        return 1;
    }

    auto Vh = Pch<Order>( mesh );
    auto ibcRange = markedfaces( mesh, ibcMarker );
    auto ibcMesh = mesh->trace( ibcRange );
    auto Jbh = Pch<Order>( ibcMesh );
    auto Mh = Pch<0>( ibcMesh );
    auto ibcSupport = support( Mh );
    auto Xh = product( Vh, Jbh, Mh );

    if ( Mh->nDof() != 1 )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "The IBC boundary submesh must be connected for Pch<0>; found "
                      << Mh->nDof() << " dofs on marker '" << ibcMarker << "'\n";
        return 1;
    }

    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto j = Jbh->element( "j" );
    auto eta = Jbh->element( "eta" );
    auto lambda = Mh->element( "lambda" );
    auto mu = Mh->element( "mu" );

    double const sigma = doption( "sigma" );
    double const current = doption( "current" );
    double const potential = doption( "potential" );
    double const radiusInner = doption( "radius-inner" );
    double const radiusOuter = doption( "radius-outer" );
    double const sectionHeight = doption( "section-height" );
    double const angle = doption( "angle" );

    if ( sigma <= 0.0 || radiusInner <= 0.0 || radiusOuter <= radiusInner || sectionHeight <= 0.0 )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "Invalid electrostatic or geometric parameters\n";
        return 1;
    }

    double const ibcMeasure = integrate( _range = elements( ibcSupport ),
                                         _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    if ( ibcMeasure <= 0.0 )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "IBC marker '" << ibcMarker << "' has zero measure\n";
        return 1;
    }

    auto rhs = blockform1( Xh, solve::strategy::monolithic, backend() );
    auto a = blockform2( Xh, solve::strategy::monolithic, backend() );

    a( 0_c, 0_c ) = integrate( _range = elements( mesh ),
                               _expr = cst( sigma ) * inner( gradt( u ), grad( v ) ) );
    a( 0_c, 1_c ) += integrate( _range = elements( ibcSupport ),
                                _expr = idt( j ) * id( v ) );
    a( 1_c, 0_c ) += integrate( _range = elements( ibcSupport ),
                                _expr = idt( u ) * id( eta ) );
    a( 1_c, 2_c ) += integrate( _range = elements( ibcSupport ),
                                _expr = -idt( lambda ) * id( eta ) );
    a( 2_c, 1_c ) += integrate( _range = elements( ibcSupport ),
                                _expr = idt( j ) * id( mu ) );

    rhs( 2_c ) += integrate( _range = elements( ibcSupport ),
                             _expr = cst( current / ibcMeasure ) * id( mu ) );

    auto U = Xh.element();
    auto uh = U( 0_c );
    rhs.close();
    a.close();
    a.row( 0_c ) += on( _range = markedfaces( mesh, dirichletMarker ),
                        _rhs = rhs( 0_c ),
                        _element = uh,
                        _expr = cst( potential ),
                        _type = "elimination" );
    if ( !boption( "no-solve" ) )
    {
        a.solve( _rhs = rhs, _solution = U );
    }

    u = U( 0_c );
    j = U( 1_c );
    lambda = U( 2_c );

    auto Jh = Pdhv<0>( mesh );
    auto J = Jh->element( "J" );
    J.on( _range = elements( mesh ),
          _expr = -cst( sigma ) * trans( gradv( u ) ) );

    double const currentDirichlet =
        integrate( _range = markedfaces( mesh, dirichletMarker ),
                   _expr = -cst( sigma ) * ( gradv( u ) * N() ) ).evaluate()( 0, 0 );
    double const currentIbc =
        integrate( _range = markedfaces( mesh, ibcMarker ),
                   _expr = -cst( sigma ) * ( gradv( u ) * N() ) ).evaluate()( 0, 0 );
    double const currentIbcLM =
        integrate( _range = elements( ibcSupport ), _expr = idv( j ) ).evaluate()( 0, 0 );
    double currentNeumann = 0.0;
    if ( nelements( markedfaces( mesh, neumannMarker ), true ) > 0 )
    {
        currentNeumann = integrate( _range = markedfaces( mesh, neumannMarker ),
                                    _expr = -cst( sigma ) * ( gradv( u ) * N() ) ).evaluate()( 0, 0 );
    }
    double const lambdaMean =
        integrate( _range = elements( ibcSupport ), _expr = idv( lambda ) ).evaluate()( 0, 0 ) / ibcMeasure;

    double const logRatio = std::log( radiusOuter / radiusInner );
    double const tangentialSlope = current / ( sigma * sectionHeight * logRatio );
    double const exactIbcPotential = potential - tangentialSlope * angle;
    double const exactCurrentCoeff = current / ( sectionHeight * logRatio );

    std::ostringstream pExactStream;
    pExactStream << std::setprecision( 17 )
                 << potential << "-" << tangentialSlope << "*atan2(y,x)"
                 << makeCoordinateMapping<Dim>();
    auto pExact = expr( pExactStream.str() );
    auto JExact = expr<Dim,1>( makeExactCurrentDensityString<Dim>( exactCurrentCoeff ) + makeCoordinateMapping<Dim>() );

    double const potentialL2Error =
        normL2( _range = elements( mesh ),
                _expr = idv( u ) - pExact );
    double const potentialL2Exact =
        normL2( _range = elements( mesh ),
                _expr = pExact );
    double const currentL2Error =
        normL2( _range = elements( mesh ),
                _expr = idv( J ) - JExact );
    double const currentL2Exact =
        normL2( _range = elements( mesh ),
                _expr = JExact );
    double const lambdaAbsError = math::abs( lambdaMean - exactIbcPotential );
    double const lambdaRelError = ( math::abs( exactIbcPotential ) > 1e-16 )
                                      ? lambdaAbsError / math::abs( exactIbcPotential )
                                      : lambdaAbsError;
    double const potentialL2Rel = ( potentialL2Exact > 1e-16 )
                                      ? potentialL2Error / potentialL2Exact
                                      : potentialL2Error;
    double const currentL2Rel = ( currentL2Exact > 1e-16 )
                                    ? currentL2Error / currentL2Exact
                                    : currentL2Error;

    if ( Environment::isMasterRank() )
    {
        std::cout << "sigma = " << sigma
                  << ", current target = " << current
                  << ", IBC measure = " << ibcMeasure << "\n"
                  << "int_" << dirichletMarker << "(J.n) = " << currentDirichlet
                  << ", int_" << ibcMarker << "(J.n) = " << currentIbc
                  << ", int_" << ibcMarker << "(j) = " << currentIbcLM
                  << ", int_" << neumannMarker << "(J.n) = " << currentNeumann
                  << ", balance = " << ( currentDirichlet + currentIbc + currentNeumann ) << "\n"
                  << "lambda mean = " << lambdaMean
                  << ", exact IBC potential = " << exactIbcPotential
                  << ", relative error = " << lambdaRelError << "\n"
                  << "relative L2(p - p_exact) = " << potentialL2Rel
                  << ", relative L2(J - J_exact) = " << currentL2Rel << std::endl;
    }

    auto e = exporter( _mesh = mesh, _name = "qs_laplacian_ibc" );
    e->addRegions();
    e->add( "potential", u );
    e->add( "current_density", J, std::set<std::string>{ "element", "nodal" } );
    e->add( "potential_exact", pExact );
    e->add( "current_density_exact", JExact );
    e->save();

    auto JbhExport = Pch<Order>( ibcMesh );
    auto jExport = JbhExport->element( "current_ibc" );
    jExport = vf::project( _space = JbhExport,
                           _range = elements( ibcMesh ),
                           _expr = idv( j ) );
    auto MhExport = Pch<0>( ibcMesh );
    auto lambdaExport = MhExport->element( "lambda" );
    lambdaExport = vf::project( _space = MhExport,
                                _range = elements( ibcMesh ),
                                _expr = idv( lambda ) );

    auto eb = exporter( _mesh = ibcMesh, _name = "qs_laplacian_ibc_boundary" );
    eb->add( "current_ibc", jExport );
    eb->add( "lambda", lambdaExport );
    eb->save();

    return 0;
}

} // namespace

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc = argc,
                     _argv = argv,
                     _desc = makeOptions(),
                     _about = makeAbout() );

    return cg_laplacian_ibc<FEELPP_DIM,1>();
}
