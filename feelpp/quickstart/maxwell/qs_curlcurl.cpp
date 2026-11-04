/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include "../discretization.hpp"
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace Feel;

namespace
{
// curl curl E + mass E = f, with the exact tangential trace on the boundary.
// The smooth sine fields give a manufactured solution; the polynomial field is a patch test.
template<typename SpacePtr>
int solveCurlCurl( SpacePtr const& Xh )
{
    auto mesh = Xh->mesh();
    bool const polynomial = soption( "solution" ) == "polynomial";
    double const mass = doption( "mass" );
    double const pi = std::acos( -1.0 );
    int const requestedQuad = ioption( "quad" );
    int const quad = requestedQuad > 0 ? requestedQuad :
        std::max( 10, 2 * ( ioption( "order" ) + 2 ) * mesh->order() );
#if FEELPP_DIM == 2
    auto exact = expr<2, 1>( polynomial ? "{-y,x}:x:y" :
                             "{sin(pi*y),sin(pi*x)}:x:y" );
    auto exactCurl = expr( polynomial ? "2" :
                           "pi*(cos(pi*x)-cos(pi*y)):x:y" );
#else
    auto exact = expr<3, 1>( polynomial ? "{-y,x,0}:x:y:z" :
                             "{sin(pi*y),sin(pi*z),sin(pi*x)}:x:y:z" );
    auto exactCurl = expr<3, 1>( polynomial ? "{0,0,2}" :
                                 "{-pi*cos(pi*z),-pi*cos(pi*x),-pi*cos(pi*y)}:x:y:z" );
#endif
    auto forcing = ( polynomial ? mass : mass + pi * pi ) * exact;
    auto u = Xh->element( "E" );
    auto v = Xh->element();

    tic();
    auto l = form1( _test=Xh );
    l = integrate( _range=elements( mesh ), _expr=inner( forcing, id( v ) ), _quad=quad );
    auto a = form2( _trial=Xh, _test=Xh );
    a = integrate( _range=elements( mesh ),
                   _expr=inner( curlt( u ), curl( v ) ) + mass * inner( idt( u ), id( v ) ),
                   _quad=quad );
    // H(curl) boundary DoFs constrain the tangential trace only.
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u, _expr=exact,
             _type="elimination_symmetric" );
    double const assemblyTime = toc( "curl-curl assembly" );

    tic();
    auto const status = a.solve( _rhs=l, _solution=u );
    double const solveTime = toc( "curl-curl solve" );
    if ( !status.isConverged() )
        throw std::runtime_error( "curl-curl linear solve did not converge" );

    double const errorL2 = normL2( _range=elements( mesh ), _expr=idv( u ) - exact, _quad=quad );
    double const errorCurl = normL2( _range=elements( mesh ), _expr=curlv( u ) - exactCurl, _quad=quad );
    double const errorHcurl = std::hypot( errorL2, errorCurl );
    Feel::cout << std::setprecision( 12 )
               << "curlcurl: dim=" << FEELPP_DIM << " family=" << soption( "family" )
               << " order=" << ioption( "order" ) << " geometry.order=" << mesh->order()
               << " ndof=" << Xh->nDof() << " h=" << mesh->hMax()
               << " L2=" << errorL2 << " curl=" << errorCurl << " Hcurl=" << errorHcurl
               << " assembly=" << assemblyTime << " solve=" << solveTime << std::endl;

    if ( !std::isfinite( errorHcurl ) )
        throw std::runtime_error( "curl-curl error is not finite" );
    if ( doption( "error-max" ) > 0 && errorHcurl > doption( "error-max" ) )
        throw std::runtime_error( "curl-curl H(curl) error exceeds error-max" );

    if ( !boption( "no-export" ) )
    {
        auto e = exporter( _mesh=mesh );
        e->add( "E", u );
        e->add( "E_exact", exact );
        e->save();
    }
    return 0;
}
}

int main( int argc, char** argv )
{
    try
    {
        po::options_description options( "Nedelec curl-curl quickstart" );
        options.add_options()
            ( "family", po::value<std::string>()->default_value( "ned1" ), "Nedelec family: ned1 or ned2" )
            ( "order", po::value<int>()->default_value( 0 ), "public Nedelec family order k >= 0" )
            ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric order" )
            ( "solution", po::value<std::string>()->default_value( "sin" ), "exact solution: sin or polynomial" )
            ( "mass", po::value<double>()->default_value( 1 ), "positive mass coefficient" )
            ( "quad", po::value<int>()->default_value( 0 ), "quadrature order (0: automatic)" )
            ( "error-max", po::value<double>()->default_value( -1 ), "maximum H(curl) error (negative: disabled)" )
            ( "no-export", po::value<bool>()->default_value( false )->implicit_value( true ), "disable field export" );
        Environment env( _argc=argc, _argv=argv, _desc=options,
                         _about=about( _name="qs_curlcurl", _author="Feel++ Consortium",
                                       _email="feelpp-devel@feelpp.org" ) );

        auto const order = quickstart::checkedOrder( "order", 0 );
        auto const geometryOrder = quickstart::geometryOrder();
        auto const family = soption( "family" );
        if ( family != "ned1" && family != "ned2" )
            throw std::invalid_argument( "family must be ned1 or ned2" );
        if ( soption( "solution" ) != "sin" && soption( "solution" ) != "polynomial" )
            throw std::invalid_argument( "solution must be sin or polynomial" );
        if ( !std::isfinite( doption( "mass" ) ) || doption( "mass" ) <= 0 )
            throw std::invalid_argument( "mass must be finite and positive" );
        if ( ioption( "quad" ) < 0 )
            throw std::invalid_argument( "quad must be non-negative" );

        auto mesh = loadMesh( _mesh=new Mesh<Simplex<FEELPP_DIM, Dynamic>>( geometryOrder ) );
        quickstart::checkMeshOrder( mesh );
        if ( family == "ned1" )
            return solveCurlCurl( Ned1h<Dynamic>( mesh, order ) );
        return solveCurlCurl( Ned2h<Dynamic>( mesh, order ) );
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
