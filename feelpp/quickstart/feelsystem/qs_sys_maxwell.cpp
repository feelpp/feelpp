/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include "../discretization.hpp"
#include "../maxwell/manufactured.hpp"
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace Feel;

namespace
{
template <typename SpacePtr>
int solveMaxwell( SpacePtr const& Xh )
{
    auto mesh = Xh->mesh();
    bool const polynomial = soption( "solution" ) == "polynomial";
    double const mass = doption( "mass" );
    double const pi = std::acos( -1.0 );
    int const requestedQuad = ioption( "quad" );
    int const quad = requestedQuad > 0 ? requestedQuad :
        std::max( 10, 2 * ( ioption( "order" ) + 2 ) * mesh->order() );
    auto exact = quickstart::maxwell::exactField<FEELPP_DIM>( polynomial );
    auto exactCurl = quickstart::maxwell::exactCurl<FEELPP_DIM>( polynomial );
    auto forcing = ( polynomial ? mass : mass + pi * pi ) * exact;
    auto Eh = Xh->element( "E" );
    auto u = trial( Xh, "E" );
    auto v = test( Xh, "v" );

    // Store the variational problem independently of its assembled system.
    auto problem = ellipticProblem( _space = Xh, _name = "maxwell" )
        .lhs( integrate( _range = elements( mesh ),
                         _expr = inner( curl( u ), curl( v ) ) + mass * inner( u, v ),
                         _quad = quad ) )
        .rhs( integrate( _range = elements( mesh ),
                         _expr = inner( forcing, v ), _quad = quad ) )
        // H(curl) essential data constrain the tangential trace only.
        .dirichlet( _range = boundaryfaces( mesh ), _element = Eh, _expr = exact,
                    _type = "elimination_symmetric" );

    tic();
    // Use the default backend so solver options match qs_curlcurl unchanged.
    auto system = problem.assemble( _backend = backend() );
    double const assemblyTime = toc( "Maxwell Problem API assembly" );

    tic();
    auto const result = system.solve( _solution = Eh );
    double const solveTime = toc( "Maxwell Problem API solve" );
    if ( !result.converged )
        throw std::runtime_error( "Maxwell linear solve did not converge" );

    double const errorL2 = normL2( _range = elements( mesh ), _expr = idv( Eh ) - exact, _quad = quad );
    double const errorCurl = normL2( _range = elements( mesh ), _expr = curlv( Eh ) - exactCurl, _quad = quad );
    double const errorHcurl = std::hypot( errorL2, errorCurl );
    Feel::cout << std::setprecision( 12 )
               << "sys_maxwell: dim=" << FEELPP_DIM << " family=" << soption( "family" )
               << " order=" << ioption( "order" ) << " geometry.order=" << mesh->order()
               << " ndof=" << Xh->nDof() << " h=" << mesh->hMax()
               << " L2=" << errorL2 << " curl=" << errorCurl << " Hcurl=" << errorHcurl
               << " assembly=" << assemblyTime << " solve=" << solveTime << std::endl;

    if ( !std::isfinite( errorHcurl ) )
        throw std::runtime_error( "Maxwell error is not finite" );
    if ( doption( "error-max" ) > 0 && errorHcurl > doption( "error-max" ) )
        throw std::runtime_error( "Maxwell H(curl) error exceeds error-max" );

    if ( !boption( "no-export" ) )
    {
        auto e = exporter( _mesh = mesh );
        e->add( "E", Eh );
        e->add( "E_exact", exact );
        e->save();
    }
    return 0;
}
} // namespace

int main( int argc, char** argv )
{
    try
    {
        po::options_description options( "Problem API Maxwell quickstart" );
        options.add_options()
            ( "family", po::value<std::string>()->default_value( "ned1" ), "Nedelec family: ned1 or ned2" )
            ( "order", po::value<int>()->default_value( 0 ), "public Nedelec family order k >= 0" )
            ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric order" )
            ( "solution", po::value<std::string>()->default_value( "sin" ), "exact solution: sin or polynomial" )
            ( "mass", po::value<double>()->default_value( 1 ), "positive mass coefficient" )
            ( "quad", po::value<int>()->default_value( 0 ), "quadrature order (0: automatic)" )
            ( "error-max", po::value<double>()->default_value( -1 ), "maximum H(curl) error (negative: disabled)" )
            ( "no-export", po::value<bool>()->default_value( false )->implicit_value( true ), "disable field export" );
        Environment env( _argc = argc, _argv = argv, _desc = options,
                         _about = about( _name = "qs_sys_maxwell", _author = "University of Strasbourg",
                                         _email = "feelpp-devel@feelpp.org" ) );

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

        auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, Dynamic>>( geometryOrder ) );
        quickstart::checkMeshOrder( mesh );
        if ( family == "ned1" )
            return solveMaxwell( Ned1h<Dynamic>( mesh, order ) );
        return solveMaxwell( Ned2h<Dynamic>( mesh, order ) );
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
