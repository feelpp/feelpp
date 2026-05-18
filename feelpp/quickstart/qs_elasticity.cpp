/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
#include "nullspace-rigidbody.hpp"

#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>


int main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        std::string const zeroVectorExpr = FEELPP_DIM == 2 ? "{0,0}" : "{0,0,0}";
        std::string const zeroMomentExpr = FEELPP_DIM == 2 ? "0" : "{0,0,0}";

        po::options_description laplacianoptions( "Elasticity options" );
        laplacianoptions.add_options()
            ( "E", po::value<double>()->default_value( 1.0e6 ), "Young modulus" )
            ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
            ( "weakdir", po::value<bool>()->default_value( false ), "use weak dirichlet" )
            ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
            ( "nullspace", po::value<bool>()->default_value( false ), "add null space" )
            ( "point-force.markers", po::value<std::vector<std::string>>()->multitoken(), "point markers receiving a point force" )
            ( "point-force.expr", po::value<std::string>()->default_value( zeroVectorExpr ), "point force vector expression" )
            ( "point-force.quantity", po::value<std::string>()->default_value( "per-point" ), "point force quantity: per-point or total" )
            ( "point-moment.markers", po::value<std::vector<std::string>>()->multitoken(), "point markers receiving a point moment" )
            ( "point-moment.expr", po::value<std::string>()->default_value( zeroMomentExpr ), "point moment expression: scalar in 2D, vector in 3D" )
            ( "point-moment.quantity", po::value<std::string>()->default_value( "per-point" ), "point moment quantity: per-point or total" )
            ( "cantilever.check", po::value<bool>()->default_value( false ), "check tip displacement against an Euler-Bernoulli cantilever reference" )
            ( "cantilever.tip-marker", po::value<std::string>()->default_value( "tip" ), "point marker used for the cantilever tip displacement check" )
            ( "cantilever.component", po::value<int>()->default_value( FEELPP_DIM-1 ), "displacement component used by the cantilever check" )
            ( "cantilever.length", po::value<double>()->default_value( 1. ), "cantilever beam length" )
            ( "cantilever.height", po::value<double>()->default_value( 1. ), "cantilever bending height" )
            ( "cantilever.thickness", po::value<double>()->default_value( 1. ), "cantilever out-of-plane thickness" )
            ( "cantilever.tip-force", po::value<double>()->default_value( 0. ), "signed transverse tip force used in the cantilever reference" )
            ( "cantilever.tip-moment", po::value<double>()->default_value( 0. ), "signed end moment used in the cantilever reference" )
            ( "cantilever.tolerance.relative", po::value<double>()->default_value( 0.5 ), "relative tolerance for the cantilever check" )
            ( "cantilever.tolerance.absolute", po::value<double>()->default_value( 1e-12 ), "absolute tolerance for the cantilever check" )
            ;

        Environment env( _argc=argc, _argv=argv,
                    _desc=laplacianoptions,
                    _about=about(_name="qs_elasticity",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        tic();
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
        toc("loadMesh");

        tic();
        auto Vh = Pchv<1>( mesh );
        toc("Vh");

        auto u = Vh->element("u");
        auto v = Vh->element("v");
        auto nu = doption(_name="nu");
        auto E = doption(_name="E");
        auto lambda = E*nu/( (1+nu)*(1-2*nu) );
        auto mu = E/(2*(1+nu));
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
        auto sigmat = lambda*trace(deft)*Id + 2*mu*deft;
        auto sigma = lambda*trace(def)*Id + 2*mu*def;
        std::map<std::string,std::string> checkerInputs{
            {"dim", std::to_string( FEELPP_DIM )},
            {"lam1", std::to_string( mu )},
            {"lam2", std::to_string( lambda )},
            {"exact", "1"},
            {"grad_displ", ""},
            {"strain", ""},
            {"stress", ""},
            {"stressn", ""},
            {"f", ""},
            {"c1", ""},
            {"c2", ""}
        };
        auto thechecker = checker( _name="L2/H1 displacement norms",
                                   _solution_key="displ",
                                   _gradient_key="grad_displ",
                                   _inputs=checkerInputs );
        auto locals = thechecker.check()? thechecker.runScript() : Checker::variables_t{};
        int status = 0;

        tic();
        auto l = form1( _test=Vh );
        if ( thechecker.check() )
        {
            auto f = expr<FEELPP_DIM,1>( locals.at( "f" ), "f" );
            l = integrate(_range=elements(mesh),
                          _expr=-inner(f,id(v)));
        }
        else
        {
            auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f"), "f" );
            l = integrate(_range=elements(mesh),
                          _expr=inner(f,id(v)));
        }
        auto markersFromOption = []( std::string const& opt )
        {
            std::set<std::string> markers;
            if ( Environment::vm().count( opt ) )
            {
                auto markerList = Environment::vm()[opt].as<std::vector<std::string>>();
                markers.insert( markerList.begin(), markerList.end() );
            }
            return markers;
        };
        auto pointLoadScale = []( auto const& pointRange, std::string const& quantity, std::string const& label )
        {
            if ( quantity == "per-point" || quantity == "per_point" || quantity == "point" )
                return 1.;
            if ( quantity != "total" )
                throw std::invalid_argument( label + " quantity must be 'per-point' or 'total'" );
            size_type nMarkedPoints = nelements( pointRange, true );
            if ( nMarkedPoints == 0 )
                throw std::invalid_argument( label + " uses quantity=total but selects no marked point" );
            return 1./nMarkedPoints;
        };

        auto pointForceMarkers = markersFromOption( "point-force.markers" );
        if ( !pointForceMarkers.empty() )
        {
            auto pointForceRange = markedpoints( mesh, pointForceMarkers );
            auto pointForceExpr = expr<FEELPP_DIM,1>( soption(_name="point-force.expr"), "point_force" );
            double pointForceScale = pointLoadScale( pointForceRange, soption(_name="point-force.quantity"), "point-force" );
            l += integrate( _range=pointForceRange,
                            _expr=pointForceScale*inner( pointForceExpr, id(v) ) );
        }

        auto pointMomentMarkers = markersFromOption( "point-moment.markers" );
        if ( !pointMomentMarkers.empty() )
        {
            auto pointMomentRange = markedpoints( mesh, pointMomentMarkers );
            double pointMomentScale = pointLoadScale( pointMomentRange, soption(_name="point-moment.quantity"), "point-moment" );
#if FEELPP_DIM == 2
            auto pointMomentExpr = expr( soption(_name="point-moment.expr"), "point_moment" );
            l += integrate( _range=pointMomentRange,
                            _expr=pointMomentScale*pointMomentExpr*omegaz( v ) );
#else
            auto pointMomentExpr = expr<FEELPP_DIM,1>( soption(_name="point-moment.expr"), "point_moment" );
            l += integrate( _range=pointMomentRange,
                            _expr=pointMomentScale*inner( pointMomentExpr, omega( v ) ) );
#endif
        }
        toc("l");

        tic();
        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                    _expr=inner( sigmat, grad(v) ) );

        if ( boption(_name="weakdir") )
        {
            double penaldir = doption(_name="gamma");
            a += integrate(_range=markedfaces(mesh,"Dirichlet"),
                        _expr=-inner(sigmat*N(),id(u)) + inner(-sigma*N()+std::max(2*mu,lambda)*penaldir*id(u)/hFace(),idt(u)) );

        }
        else
        {
            auto g = expr<FEELPP_DIM,1>( thechecker.check()? thechecker.solution() : soption(_name="functions.g"), "g" );
            a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
        }
        toc("a");

        //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
        if ( !boption( "no-solve" ) )
        {
            tic();
            auto b = backend();
            auto rigidBodyModes = std::make_shared<NullSpace<double>>( b, qsNullSpace( Vh ) );
            b->attachNearNullSpace( rigidBodyModes );
            if ( boption(_name="nullspace") )
                b->attachNullSpace( rigidBodyModes );

            a.solve(_rhs=l,_solution=u);
            toc("a.solve");
        }

        if ( boption(_name="cantilever.check") && !boption( "no-solve" ) )
        {
            int component = ioption(_name="cantilever.component");
            if ( component < 0 || component >= FEELPP_DIM )
                throw std::invalid_argument( "cantilever.component must be in [0," + std::to_string( FEELPP_DIM-1 ) + "]" );

            std::string tipMarker = soption(_name="cantilever.tip-marker");
            auto tipRange = markedpoints( mesh, tipMarker );
            size_type nTipPoints = nelements( tipRange, true );
            if ( nTipPoints == 0 )
                throw std::invalid_argument( "cantilever tip marker '" + tipMarker + "' selects no point" );

            auto tipComponentDofs = Vh->dofs( tipRange, static_cast<ComponentType>( component ), false );
            double tipDisplacementSum = 0.;
            size_type nTipComponentDofs = 0;
            for ( auto dofId : tipComponentDofs )
            {
                if ( Vh->dof()->dofGlobalProcessIsGhost( dofId ) )
                    continue;
                tipDisplacementSum += u( dofId );
                ++nTipComponentDofs;
            }
            mpi::all_reduce( Environment::worldComm(), mpi::inplace( tipDisplacementSum ), std::plus<double>() );
            mpi::all_reduce( Environment::worldComm(), mpi::inplace( nTipComponentDofs ), std::plus<size_type>() );
            if ( nTipComponentDofs != nTipPoints )
                throw std::runtime_error( "cantilever tip marker '" + tipMarker + "' resolved to " +
                                          std::to_string( nTipComponentDofs ) + " active component dofs but " +
                                          std::to_string( nTipPoints ) + " marked points" );
            double computed = tipDisplacementSum/nTipComponentDofs;

            double length = doption(_name="cantilever.length");
            double height = doption(_name="cantilever.height");
            double thickness = doption(_name="cantilever.thickness");
            double tipForce = doption(_name="cantilever.tip-force");
            double tipMoment = doption(_name="cantilever.tip-moment");
            if ( length <= 0 )
                throw std::invalid_argument( "cantilever.length must be positive" );
            if ( height <= 0 )
                throw std::invalid_argument( "cantilever.height must be positive" );
            if ( thickness <= 0 )
                throw std::invalid_argument( "cantilever.thickness must be positive" );

            double inertia = thickness*std::pow( height, 3 )/12.;
            double expected = tipForce*std::pow( length, 3 )/( 3.*E*inertia )
                              + tipMoment*std::pow( length, 2 )/( 2.*E*inertia );
            double error = std::abs( computed - expected );
            double tolerance = std::max( doption(_name="cantilever.tolerance.absolute"),
                                         doption(_name="cantilever.tolerance.relative")*std::max( std::abs( expected ), 1e-30 ) );

            if ( Environment::isMasterRank() )
            {
                std::cout << "cantilever Euler-Bernoulli reference check\n"
                          << "  component: " << component << "\n"
                          << "  computed tip displacement: " << computed << "\n"
                          << "  reference tip displacement: " << expected << "\n"
                          << "  absolute error: " << error << "\n"
                          << "  tolerance: " << tolerance << std::endl;
            }

            auto cantileverChecker = checker( _name="Euler-Bernoulli cantilever reference",
                                              _solution="tip displacement" );
            cantileverChecker.setCheck( true );
            auto cantileverError = [computed, expected]( std::string const& )
            {
                return std::map<std::string,double>{
                    { "tip-displacement", std::abs( computed - expected ) }
                };
            };
            auto cantileverTolerance = [tolerance]( std::string const&, std::pair<std::string,double> const& r, double, double )
            {
                if ( r.second <= tolerance )
                    return Checks::EXACT;
                throw CheckerExactFailed( r.second, tolerance );
            };
            status = !cantileverChecker.runOnce( cantileverError, cantileverTolerance, "Euler-Bernoulli " );
        }

        tic();
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "u", u );
        if ( thechecker.check() )
            e->add( "u_exact", expr<FEELPP_DIM,1>( thechecker.solution() ) );
        e->save();
        toc("Exporter");
        int checkerStatus = check( thechecker, u );
        return status || checkerStatus;
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}
