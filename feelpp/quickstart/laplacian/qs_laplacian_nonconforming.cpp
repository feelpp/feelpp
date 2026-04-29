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
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @date 29 Apr 2026
//! @copyright 2026 Feel++ Consortium

#include <feel/feel.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        po::options_description options( "Nonconforming Laplacian options" );
        options.add_options()
            ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
            ( "no-solve", po::value<bool>()->default_value( false ), "assemble without solving" );
        options.add( feel_options() );

        Environment env( _argc = argc, _argv = argv,
                         _desc = options,
                         _about = about( _name = "qs_laplacian_nonconforming",
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp-devel@feelpp.org" ) );

        Environment::changeRepository( _directory = boost::format( "qs_laplacian/nonconforming/h_%1%/" )
                                                     % doption( _name = "hsize" ) );

        using mesh_type = Mesh<Simplex<2>>;
        using basis_type = bases<CrouzeixRaviart<1, Scalar>>;
        using space_type = FunctionSpace<mesh_type, basis_type>;

        auto mesh = unitSquare( doption( _name = "hsize" ) );
        auto Xh = space_type::New( mesh );

        auto u = Xh->element( "u" );
        auto v = Xh->element( "v" );

        auto u_exact = expr( "x*(1-x)*y*(1-y):x:y" );
        auto dx_exact = expr( "(1-2*x)*y*(1-y):x:y" );
        auto dy_exact = expr( "x*(1-x)*(1-2*y):x:y" );
        auto f = expr( "2*(x*(1-x)+y*(1-y)):x:y" );

        auto l = form1( _test = Xh );
        l = integrate( _range = elements( mesh ),
                       _expr = f * id( v ),
                       _quad = _Q<6>() );

        auto a = form2( _trial = Xh, _test = Xh );
        a = integrate( _range = elements( mesh ),
                       _expr = gradt( u ) * trans( grad( v ) ),
                       _quad = _Q<6>() );
        a += on( _range = boundaryfaces( mesh ),
                 _rhs = l,
                 _element = u,
                 _expr = cst( 0. ) );

        if ( !boption( _name = "no-solve" ) )
            a.solve( _rhs = l, _solution = u );

        if ( Environment::isMasterRank() )
        {
            std::cout << "CR1 scalar dofs: " << Xh->nDof() << "\n";
            if ( !boption( _name = "no-solve" ) )
            {
                double const l2 = normL2( _range = elements( mesh ),
                                          _expr = idv( u ) - u_exact,
                                          _quad = _Q<6>() );
                double const h1semi = math::sqrt( integrate( _range = elements( mesh ),
                                                             _expr = ( dxv( u ) - dx_exact ) * ( dxv( u ) - dx_exact ) +
                                                                     ( dyv( u ) - dy_exact ) * ( dyv( u ) - dy_exact ),
                                                             _quad = _Q<6>() )
                                                   .evaluate()( 0, 0 ) );
                std::cout << "||u-u_h||_L2 = " << l2 << "\n";
                std::cout << "|u-u_h|_H1,broken = " << h1semi << "\n";
            }
        }

        auto e = exporter( _mesh = mesh );
        e->addRegions();
        e->add( "u", u );
        e->add( "u_exact", u_exact );
        e->save();

        return 0;
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
