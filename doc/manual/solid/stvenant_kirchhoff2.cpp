/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-02-24

  Copyright (C) 2008-2012 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2012-2014 Feel++ Consortium

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
#if 0
#include <feel/feel.hpp>
#else
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelts/newmark.hpp>
#endif
#include <indicators/progress_spinner.hpp>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stvenantkirchhoffoptions( "StVenantKirchhoff problem options" );
    stvenantkirchhoffoptions.add_options()
    ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "mesh element shape (simplex or hypercube)" )
    ( "fe-order", Feel::po::value<int>()->default_value( 1 ), "finite element polynomial order (1,2,3)" )
    ( "young-modulus", Feel::po::value<double>()->default_value( 1.4e6 ), "young-modulus" )
    ( "poisson-coeff", Feel::po::value<double>()->default_value( 0.4 ), "poisson-coeff" )
    ( "rho", Feel::po::value<double>()->default_value( 1000 ), "density [kg/m^3]" )
    ( "gravity-cst", Feel::po::value<double>()->default_value( 2 ), "gravity-cst" )
    ;
    return stvenantkirchhoffoptions;
}

template<int Order, typename ConvexType>
int
runStVenantKirchhoff()
{
    using namespace Feel;
    static const uint16_type nDim = ConvexType::nRealDim;
    static_assert( ConvexType::nDim == 2 && ConvexType::nRealDim == 2, "stvenant_kirchhoff2 supports 2D convexes only" );
    static_assert( Order >= 1 && Order <= 3, "stvenant_kirchhoff2 supports FE order in {1,2,3}" );

    double meshSize = doption(_name="gmsh.hsize");
    double youngmodulus= doption(_name="young-modulus");
    double coeffpoisson= doption(_name="poisson-coeff");
    double coefflame2 = youngmodulus/(2*(1+coeffpoisson));// mu
    double coefflame1 = youngmodulus*coeffpoisson/((1+coeffpoisson)*(1-2*coeffpoisson));// lambda
    double rho= doption(_name="rho");
    double gravityCst=doption(_name="gravity-cst");

    typedef Mesh<ConvexType> mesh_type;
    GeoTool::Node x1( (0.4+math::sqrt(0.0096))/2.,0.19 );
    GeoTool::Node x2( 0.6,0.21 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="fixe",_marker4=true);
    R.setMarker(_type="line",_name="free",_marker1=true,_marker2=true,_marker3=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name="domainRectangle" );

    auto Vh = Pchv<Order>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto Res = backend()->newVector( Vh );
    auto Jac = backend()->newMatrix( _test=Vh, _trial=Vh );

    auto e = exporter( _mesh=mesh );

    auto ts = newmark( _space=Vh, _name="structure",_rank_proc_in_files_name=true );
    auto Id = eye<nDim,nDim>();
    auto gravityForce = -rho*gravityCst*oneY();
    std::unique_ptr<indicators::ProgressSpinner> progress;
    if ( Environment::isMasterRank() )
    {
        progress = std::make_unique<indicators::ProgressSpinner>(
            indicators::option::PrefixText{ "time " },
            indicators::option::ShowSpinner{ true },
            indicators::option::ShowElapsedTime{ true },
            indicators::option::ShowRemainingTime{ false },
            indicators::option::ShowPercentage{ false },
            indicators::option::MaxProgress{ std::numeric_limits<std::size_t>::max() }
        );
    }

    // start or restart
    if ( !ts->isRestart() )
    {
      ts->start();
    }
    else
    {
        double ti = ts->restart();
        u = ts->previousUnknown();
        if ( e->doExport() ) e->restart(ti);
    }

    for ( ; !ts->isFinished(); ts->next(u) )
    {
        if ( progress )
        {
            std::ostringstream oss;
            oss << "t=" << std::fixed << std::setprecision(4) << ts->time() << "s";
            progress->set_option( indicators::option::PostfixText{ oss.str() } );
            progress->tick();
        }

        auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
            {
                auto u = Vh->element();
                u = *X;
                auto Fv = Id + gradv(u);
                auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
                auto Sv = coefflame1*trace(Ev)*Id + 2*coefflame2*Ev;
                auto dF = gradt(u);
                auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
                auto dS = coefflame1*trace(dE)*Id + 2*coefflame2*dE;

                auto a = form2( _test=Vh, _trial=Vh, _matrix=J );

                a = integrate( _range=elements( mesh ),
                               _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(v) ) );
                a += integrate( _range=elements( mesh ),
                                _expr= rho*inner( ts->polyDerivCoefficient()*idt(u),id( v ) ) );

                auto RR = backend()->newVector( Vh );
                a += on( _range=markedfaces(mesh,"fixe"),
                         _element=u, _rhs=RR,
                         _expr=zero<nDim,1>() );
            };
        auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
            {
                auto u = Vh->element();
                u = *X;

                auto Fv = Id + gradv(u);
                auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
                auto Sv = coefflame1*trace(Ev)*Id + 2*coefflame2*Ev;

                auto r = form1( _test=Vh, _vector=R );
                r = integrate( _range=elements( mesh ),
                               _expr= inner( val(Fv*Sv) , grad(v) ) );
                r += integrate( _range=elements( mesh ),
                                _expr= -inner(gravityForce,id( v ) ) );
                r += integrate( _range=elements( mesh ),
                                _expr= rho*inner( ts->polyDerivCoefficient()*idv(u) -idv(ts->polyDeriv()),id( v ) ) );

                R->close();
                auto temp = Vh->element();
                temp = *R;
                temp.on( _range=markedfaces(mesh,"fixe"),_expr=zero<nDim,1>() );
                *R = temp;
            };

        u.on( _range=markedfaces(mesh,"fixe"),_expr=zero<nDim,1>() );
        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u,_jacobian=Jac,_residual=Res );

        ts->updateFromDisp(u);
        e->step(ts->time())->add( "displacement", u );
        e->step(ts->time())->add( "velocity", ts->currentVelocity() );
        e->step(ts->time())->add( "acceleration", ts->currentAcceleration() );
        e->save();
    }
    if ( progress )
        progress->mark_as_completed();

    return 0;
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    try
    {
        Environment env( _argc=argc, _argv=argv,
                         _desc=makeOptions(),
                         _about=about(_name="stvenant_kirchhoff2",
                                      _author="Vincent Chabannes",
                                      _email="vincent.chabannes@feelpp.org"));

        std::string shape = soption(_name="shape");
        int order = ioption(_name="fe-order");
        if ( shape == "simplex" || shape == "Simplex" )
        {
            switch ( order )
            {
            case 1: return runStVenantKirchhoff<1,Simplex<2,1,2>>();
            case 2: return runStVenantKirchhoff<2,Simplex<2,1,2>>();
            case 3: return runStVenantKirchhoff<3,Simplex<2,1,2>>();
            default:
                throw std::invalid_argument( "Unsupported fe-order '" + std::to_string(order) + "' (expected 1, 2, or 3)" );
            }
        }
        if ( shape == "hypercube" || shape == "Hypercube" )
        {
            switch ( order )
            {
            case 1: return runStVenantKirchhoff<1,Hypercube<2,1,2>>();
            case 2: return runStVenantKirchhoff<2,Hypercube<2,1,2>>();
            case 3: return runStVenantKirchhoff<3,Hypercube<2,1,2>>();
            default:
                throw std::invalid_argument( "Unsupported fe-order '" + std::to_string(order) + "' (expected 1, 2, or 3)" );
            }
        }

        throw std::invalid_argument( "Unsupported shape '" + shape + "' (expected simplex or hypercube)" );
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}
