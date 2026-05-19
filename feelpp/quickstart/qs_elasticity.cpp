/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
#include "nullspace-rigidbody.hpp"
#include "qs_elasticity_case.hpp"
#include "qs_elasticity_checks.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef FEELPP_ORDER
#define FEELPP_ORDER 1
#endif

namespace
{
namespace qsec = Feel::Quickstart::ElasticityChecks;
namespace qsecase = Feel::Quickstart::ElasticityCase;

using PointConstraintConfig = qsecase::PointConstraintConfig;
} // namespace

int main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description laplacianoptions( "Elasticity options" );
        qsecase::addOptions( laplacianoptions );

        Environment env( _argc=argc, _argv=argv,
                    _desc=laplacianoptions,
                    _about=about(_name="qs_elasticity",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        auto caseConfig = qsecase::fromEnvironment();

        tic();
#if FEELPP_HYPERCUBE == 1
        using mesh_type = Mesh<Hypercube<FEELPP_DIM,FEELPP_ORDER>>;
#else
        using mesh_type = Mesh<Simplex<FEELPP_DIM,FEELPP_ORDER>>;
#endif
        auto mesh = caseConfig.meshFilename.empty()
            ? loadMesh(_mesh=new mesh_type)
            : loadMesh(_mesh=new mesh_type, _filename=caseConfig.meshFilename);
        if ( caseConfig.expectedElements )
        {
            auto const nElements = nelements( elements( mesh ), true );
            if ( nElements != *caseConfig.expectedElements )
                throw std::runtime_error( "mesh element count check failed: expected " +
                                          std::to_string( *caseConfig.expectedElements ) +
                                          ", got " + std::to_string( nElements ) );
        }
        toc("loadMesh");

        tic();
        auto Vh = Pchv<FEELPP_ORDER>( mesh );
        toc("Vh");

        auto u = Vh->element("u");
        auto v = Vh->element("v");
        auto nu = caseConfig.nu;
        auto E = caseConfig.E;
        auto lambda = E*nu/( (1+nu)*(1-2*nu) );
        auto mu = E/(2*(1+nu));

        // Linear isotropic elasticity: sigma(u) = lambda tr(epsilon(u)) I + 2 mu epsilon(u).
        // sigmat carries the trial-field stress for the stiffness term, sigma the test-field stress
        // used below in the symmetric weak Dirichlet contribution.
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
        std::string const checkerSolution = caseConfig.checker.solution.value_or( std::string{} );
        std::string const checkerScript = caseConfig.checker.script.value_or( std::string{} );
        auto thechecker = checker( _name="L2/H1 displacement norms",
                                   _solution_key="displ",
                                   _gradient_key="grad_displ",
                                   _solution=checkerSolution,
                                   _script=checkerScript,
                                   _compute_pde_coefficients=!checkerScript.empty(),
                                   _inputs=checkerInputs );
        if ( caseConfig.checker.check )
            thechecker.setCheck( *caseConfig.checker.check );
        if ( caseConfig.checker.exact )
            thechecker.setExact( *caseConfig.checker.exact );
        if ( caseConfig.checker.exactTolerance )
            thechecker.setExactTolerance( *caseConfig.checker.exactTolerance );
        if ( caseConfig.checker.orderTolerance )
            thechecker.setOrderTolerance( *caseConfig.checker.orderTolerance );

        bool const useManufacturedScript = thechecker.check() && thechecker.useScript() && !checkerScript.empty();
        auto locals = useManufacturedScript ? thechecker.runScript() : Checker::variables_t{};
        qsec::ElasticityReferenceChecker<FEELPP_DIM> checks( "qs_elasticity" );
        checks.setConfig( caseConfig.referenceChecks );

        tic();
        auto l = form1( _test=Vh );
        // Right-hand side: body force. In scripted checker mode the manufactured script provides
        // the strong-form residual with the opposite sign convention.
        if ( useManufacturedScript )
        {
            auto f = expr<FEELPP_DIM,1>( locals.at( "f" ), "f" );
            l = integrate(_range=elements(mesh),
                          _expr=-inner(f,id(v)));
        }
        else
        {
            auto f = expr<FEELPP_DIM,1>( caseConfig.bodyForceExpression, "f" );
            l = integrate(_range=elements(mesh),
                          _expr=inner(f,id(v)));
        }

        for ( auto const& load : caseConfig.faceTractions )
        {
            auto tractionExpr = expr<FEELPP_DIM,1>( load.expression, "traction" );
            // Neumann virtual work: prescribed surface traction t . v on marked faces.
            l += integrate( _range=markedfaces( mesh, load.markers ),
                            _expr=inner( tractionExpr, id( v ) ) );
        }

        for ( auto const& load : caseConfig.facePressures )
        {
            auto pressureExpr = expr( load.expression, "pressure" );
            // Signed normal pressure: p n . v on the selected face marker.
            l += integrate( _range=markedfaces( mesh, load.marker ),
                            _expr=pressureExpr*inner( N(), id( v ) ) );
        }

        for ( auto const& load : caseConfig.faceTotalForces )
        {
            auto forceRange = markedfaces( mesh, load.markers );
            double area = integrate( _range=forceRange, _expr=cst( 1. ) ).evaluate()( 0,0 );
            if ( area <= 0. )
                throw std::invalid_argument( "face total force selects no marked face" );
            auto totalForceExpr = expr<FEELPP_DIM,1>( load.expression, "face_total_force" );
            // Total force distributed uniformly over its marked faces: (F/|Gamma|) . v.
            l += integrate( _range=forceRange,
                            _expr=inner( ( 1./area )*totalForceExpr, id( v ) ) );
        }

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

        for ( auto const& load : caseConfig.pointForces )
        {
            auto pointForceRange = markedpoints( mesh, load.markers );
            auto pointForceExpr = expr<FEELPP_DIM,1>( load.expression, "point_force" );
            double pointForceScale = pointLoadScale( pointForceRange, load.quantity, "point-force" );
            // Concentrated force virtual work: F . v evaluated at the selected physical points.
            l += integrate( _range=pointForceRange,
                            _expr=pointForceScale*inner( pointForceExpr, id(v) ) );
        }

        for ( auto const& load : caseConfig.pointMoments )
        {
            auto pointMomentRange = markedpoints( mesh, load.markers );
            double pointMomentScale = pointLoadScale( pointMomentRange, load.quantity, "point-moment" );
#if FEELPP_DIM == 2
            auto pointMomentExpr = expr( load.expression, "point_moment" );
            // Concentrated moment virtual work in 2D: M_z omega_z(v).
            l += integrate( _range=pointMomentRange,
                            _expr=pointMomentScale*pointMomentExpr*omegaz( v ) );
#else
            auto pointMomentExpr = expr<FEELPP_DIM,1>( load.expression, "point_moment" );
            // Concentrated moment virtual work in 3D: M . omega(v).
            l += integrate( _range=pointMomentRange,
                            _expr=pointMomentScale*inner( pointMomentExpr, omega( v ) ) );
#endif
        }
        toc("l");

        tic();
        auto a = form2( _trial=Vh, _test=Vh);
        // Internal virtual work: integral_Omega sigma(u) : grad(v).
        a = integrate(_range=elements(mesh),
                    _expr=inner( sigmat, grad(v) ) );

        if ( !caseConfig.dirichletMarkers.empty() && boption(_name="weakdir") )
        {
            double penaldir = doption(_name="gamma");
            // Symmetric Nitsche Dirichlet condition on marked faces:
            // consistency, adjoint-consistency, and penalty terms.
            a += integrate(_range=markedfaces(mesh,caseConfig.dirichletMarkers),
                        _expr=-inner(sigmat*N(),id(u)) + inner(-sigma*N()+std::max(2*mu,lambda)*penaldir*id(u)/hFace(),idt(u)) );

        }
        else if ( !caseConfig.dirichletMarkers.empty() )
        {
            auto g = expr<FEELPP_DIM,1>( thechecker.check()? thechecker.solution() : caseConfig.dirichletExpression, "g" );
            a+=on(_range=markedfaces(mesh,caseConfig.dirichletMarkers), _rhs=l, _element=u, _expr=g );
        }

        auto applyPointConstraintComponent = [&]( PointConstraintConfig const& constraint, int component )
        {
            if ( !constraint.components[component] )
                return;

            auto pointValue = expr( constraint.values[component], "point_constraint" );
            // Point displacement constraints prescribe selected scalar components at physical points.
            switch ( component )
            {
            case 0:
            {
                auto ux = u[ComponentType::X];
                a += on( _range=markedpoints( mesh, constraint.marker ),
                         _rhs=l,
                         _element=ux,
                         _expr=pointValue );
                break;
            }
            case 1:
            {
                auto uy = u[ComponentType::Y];
                a += on( _range=markedpoints( mesh, constraint.marker ),
                         _rhs=l,
                         _element=uy,
                         _expr=pointValue );
                break;
            }
            case 2:
            {
                if constexpr ( FEELPP_DIM == 3 )
                {
                    auto uz = u[ComponentType::Z];
                    a += on( _range=markedpoints( mesh, constraint.marker ),
                             _rhs=l,
                             _element=uz,
                             _expr=pointValue );
                }
                break;
            }
            }
        };

        for ( auto const& constraint : caseConfig.pointConstraints )
            for ( int component = 0; component < FEELPP_DIM; ++component )
                applyPointConstraintComponent( constraint, component );
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

        checks.add( "manufactured displacement", thechecker.check(), [&]() {
            return check( thechecker, u );
        } );

        checks.addCantileverChecks( mesh, Vh, u, E, !boption( "no-solve" ) )
              .addDisplacementProbe( Vh, u, !boption( "no-solve" ) )
              .addSmallStrainStressChecks( mesh, u, lambda, mu, !boption( "no-solve" ) );

        tic();
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "u", u );
        if ( thechecker.check() )
            e->add( "u_exact", expr<FEELPP_DIM,1>( thechecker.solution() ) );
        e->save();
        toc("Exporter");
        return checks.run();
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}
