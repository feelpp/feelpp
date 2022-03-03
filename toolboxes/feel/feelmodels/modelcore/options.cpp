/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-08-01

  Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.hpp>
#include <feel/feelcore/feel.hpp>

namespace Feel
{

Feel::po::options_description modelbase_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"verbose").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose")
        (prefixvm(prefix,"verbose_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        (prefixvm(prefix,"timers.activated").c_str(), Feel::po::value<bool>()->default_value( true ), "timers.activated")
        (prefixvm(prefix,"timers.save-master-rank").c_str(), Feel::po::value<bool>()->default_value( true ), "timers.save-master-rank")
        (prefixvm(prefix,"timers.save-max").c_str(), Feel::po::value<bool>()->default_value( false ), "timers.save-max")
        (prefixvm(prefix,"timers.save-min").c_str(), Feel::po::value<bool>()->default_value( false ), "timers.save-min")
        (prefixvm(prefix,"timers.save-mean").c_str(), Feel::po::value<bool>()->default_value( false ), "timers.save-mean")
        (prefixvm(prefix,"timers.save-all").c_str(), Feel::po::value<bool>()->default_value( false ), "timers.save-all")
        (prefixvm(prefix,"scalability-save").c_str(), Feel::po::value<bool>()->default_value( false ), "save-scalability")
        (prefixvm(prefix,"scalability-reinit-savefile").c_str(), Feel::po::value<bool>()->default_value( false ), "reinit-savefile")
        (prefixvm(prefix,"scalability-path").c_str(), Feel::po::value< std::string >(), "scalability-path")
        (prefixvm(prefix,"scalability-filename").c_str(), Feel::po::value< std::string >(), "scalability-filename")
        (prefixvm(prefix,"upload").c_str(), Feel::po::value< std::string >()->default_value(""), "upload decription")
        ;
    return appliBaseOptions;
}

Feel::po::options_description modelalgebraic_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base Methods Num options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"verbose_solvertimer").c_str(), Feel::po::value<bool>()->default_value( false/*true*/ ), "true or false to view verbose")
        (prefixvm(prefix,"verbose_solvertimer_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        (prefixvm(prefix,"linearsystem-cst-update").c_str(), Feel::po::value< bool >()->default_value( true ), "update matrix and rhs cst part")
        (prefixvm(prefix,"jacobian-linear-update").c_str(), Feel::po::value< bool >()->default_value( true ), "update linear part of the jacobian")
        (prefixvm(prefix,"residual-uselinearjac").c_str(), Feel::po::value< bool >()->default_value( true ), "update linear part of the residual with linear jacobian")
        (prefixvm(prefix,"use-cst-matrix").c_str(), Feel::po::value< bool >()->default_value( true ), "use-cst-matrix")
        (prefixvm(prefix,"use-cst-vector").c_str(), Feel::po::value< bool >()->default_value( true ), "use-cst-vector")
        (prefixvm(prefix,"error-if-solver-not-converged").c_str(), Feel::po::value< bool >()->default_value( true ), "error-if-solver-not-converged")
        (prefixvm(prefix,"clear-preconditioner-after-use").c_str(), Feel::po::value< bool >()->default_value( false ), "clear-preconditioner-after-use")
        (prefixvm(prefix,"graph-print-python").c_str(), Feel::po::value<bool>()->default_value( false ), "print graph in python script")
        (prefixvm(prefix,"graph-print-python-filename").c_str(), Feel::po::value< std::string >(), "filename python graph")

        (prefixvm(prefix,"pseudo-transient-continuation").c_str(), Feel::po::value<bool>()->default_value( false ), "use or not the pseudo-transient-continuation")
        (prefixvm(prefix,"pseudo-transient-continuation.evolution").c_str(), Feel::po::value<std::string>()->default_value( "SER" ), "evolution method : SER, EXPur")
        (prefixvm(prefix,"pseudo-transient-continuation.delta0").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "pseudo-transient-continuation parameter : delta0")
        (prefixvm(prefix,"pseudo-transient-continuation.delta-max").c_str(), Feel::po::value<double>()->default_value( 1.0e12 ), "pseudo-transient-continuation parameter : deltaMax")
        (prefixvm(prefix,"pseudo-transient-continuation.ser-variant").c_str(), Feel::po::value<std::string>()->default_value( "residual" ), "pseudo-transient-continuation ser variant : residual, solution")
        (prefixvm(prefix,"pseudo-transient-continuation.expur.threshold-high").c_str(), Feel::po::value<double>()->default_value( 1. ), "pseudo-transient-continuation parameter : threshold-high")
        (prefixvm(prefix,"pseudo-transient-continuation.expur.threshold-low").c_str(), Feel::po::value<double>()->default_value( 0.1 ), "pseudo-transient-continuation parameter : threshold-low")
        (prefixvm(prefix,"pseudo-transient-continuation.expur.beta-high").c_str(), Feel::po::value<double>()->default_value( 1.5 ), "pseudo-transient-continuation parameter : beta-high")
        (prefixvm(prefix,"pseudo-transient-continuation.expur.beta-low").c_str(), Feel::po::value<double>()->default_value( 0.1 ), "pseudo-transient-continuation parameter : beta-low")

        (prefixvm(prefix,"solver.picard.relaxation-parameter").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "solver.picard.relaxation-parameter")
        ;
    return appliBaseOptions.add( modelbase_options(prefix ) ).add( on_options( prefix ) );//.add( backend_options( prefix ) );
}


Feel::po::options_description modelnumerical_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"filename").c_str(), Feel::po::value<std::string>()->default_value( "" ), "json file describing model properties" )
        //(prefixvm(prefix,"mesh.filename").c_str(), Feel::po::value< std::string >(), "input mesh or geo file")
        (prefixvm(prefix,"geomap").c_str(), Feel::po::value< std::string >()->default_value("opt"), "geomap strategy : ho, opt ")

        ( prefixvm( prefix, "ts.order" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "time order" )

        ( prefixvm( prefix, "checker" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use checker" )

        ( prefixvm( prefix, "exporter.use-static-mesh" ).c_str(), Feel::po::value<bool>()->default_value( true ), "exporter.use-static-mesh" )
        ;

    return appliBaseOptions
        .add( mesh_options( prefix ) )
        .add( gmsh_options( prefix ) )
        .add( modelalgebraic_options( prefix ))
        .add( backend_options( prefix ) )
        .add( ptree_options( prefix ) );
}

Feel::po::options_description
coefficientformpde_options(std::string const& prefix)
{
    Feel::po::options_description cfpdeOptions("coefficient-form-pde options");
    cfpdeOptions.add_options()
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "numeric solver : automatic, Newton, Picard")
        (prefixvm(prefix,"time-stepping").c_str(), Feel::po::value< std::string >()->default_value("BDF"), "time integration schema : BDF, Theta")
        (prefixvm(prefix,"time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")

        (prefixvm(prefix,"stabilization").c_str(), Feel::po::value<bool>()->default_value( false ), "apply stabilization method")
        (prefixvm(prefix,"stabilization.type").c_str(), Feel::po::value<std::string>()->default_value( "gls" ), "supg,gls,unusual-gls")
        (prefixvm(prefix,"stabilization.gls.parameter.method").c_str(), Feel::po::value<std::string>()->default_value( "eigenvalue" ), "method used for compute tau : eigenvalue, doubly-asymptotic-approximation")
        (prefixvm(prefix,"stabilization.gls.parameter.hsize.method").c_str(), Feel::po::value<std::string>()->default_value( "hmin" ), "hmin,h,meas")
        (prefixvm(prefix,"stabilization.gls.parameter.eigenvalue.penal-lambdaK").c_str(), Feel::po::value<double>()->default_value( 0. ), "apply stabilization method")

        (prefixvm(prefix,"stabilization.do-assembly-with-grad-diffusion-coeff").c_str(), Feel::po::value<bool>()->default_value( true ), "apply stabilization method")

        (prefixvm(prefix,"stabilization.gls.shock-capturing").c_str(), Feel::po::value<bool>()->default_value( false ), "apply shock capturing in gls stabilization method")
        (prefixvm(prefix,"stabilization.gls.shock-capturing.quad").c_str(), Feel::po::value<int>()->default_value( -1 ), "apply shock capturing in gls stabilization method")
        ;

    return cfpdeOptions.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}

Feel::po::options_description
coefficientformpdes_options(std::string const& prefix)
{
    Feel::po::options_description cfpdesOptions("coefficient-form-pdes options");
    cfpdesOptions.add_options()
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "numeric solver : automatic, Newton, Picard")

        (prefixvm(prefix,"time-stepping").c_str(), Feel::po::value< std::string >()->default_value("BDF"), "time integration schema : BDF, Theta")
        (prefixvm(prefix,"time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")

        (prefixvm(prefix,"stabilization").c_str(), Feel::po::value<bool>()->default_value( false ), "apply stabilization method")
        (prefixvm(prefix,"stabilization.type").c_str(), Feel::po::value<std::string>()->default_value( "gls" ), "supg,gls,unusual-gls")
        (prefixvm(prefix,"stabilization.gls.parameter.method").c_str(), Feel::po::value<std::string>()->default_value( "eigenvalue" ), "method used for compute tau : eigenvalue, doubly-asymptotic-approximation")
        (prefixvm(prefix,"stabilization.gls.parameter.hsize.method").c_str(), Feel::po::value<std::string>()->default_value( "hmin" ), "hmin,h,meas")
        (prefixvm(prefix,"stabilization.gls.parameter.eigenvalue.penal-lambdaK").c_str(), Feel::po::value<double>()->default_value( 0. ), "apply stabilization method")

        (prefixvm(prefix,"stabilization.gls.shock-capturing").c_str(), Feel::po::value<bool>()->default_value( false ), "apply shock capturing in gls stabilization method")
        ;
    return cfpdesOptions.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}

/**
 * generate options for the fluid solver
 */
Feel::po::options_description
densityviscosity_options(std::string const& prefix)
{
    Feel::po::options_description densityviscosityOptions("DensityViscosity options");
    densityviscosityOptions.add_options()
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1050 ), "density [ kg.m^3]")
        (prefixvm(prefix,"mu").c_str(), Feel::po::value<double>()->default_value( 0.00345 ), "dynamic viscosity [ Pa.s = kg/(m.s^2) ]")
        (prefixvm(prefix,"viscosity.law").c_str(), Feel::po::value< std::string >()->default_value("newtonian"), "newtonian, power_law, walburn-schneck_law, carreau_law, carreau-yasuda_law ")
        (prefixvm(prefix,"viscosity.zero_shear").c_str(), Feel::po::value< double >()->default_value( 0.056 ), "parameter mu_0 for generalized Newtonian [ Pa.s ]  ")
        (prefixvm(prefix,"viscosity.infinite_shear").c_str(), Feel::po::value< double >()->default_value( 0.00345 ), "parameter mu_inf for generalized Newtonian [ Pa.s ] ")
        (prefixvm(prefix,"viscosity.min").c_str(), Feel::po::value< double >()->default_value( 0.00345 ), "min viscosity with power law [ Pa.s ]  ")
        (prefixvm(prefix,"viscosity.max").c_str(), Feel::po::value< double >()->default_value( 0.056 ), "max viscosity with power law [ Pa.s ] ")
        // (prefixvm(prefix,"hematocrit").c_str(), Feel::po::value< double >()->default_value( 40 ), "hematocrit : RBC volume fraction (Vrbc/Vtotal) [ percentage ] ")
        (prefixvm(prefix,"power_law.n").c_str(), Feel::po::value< double >()->default_value( 0.6 ), "parameter n in power_law ")
        (prefixvm(prefix,"power_law.k").c_str(), Feel::po::value< double >()->default_value( 0.035 ), "parameter k in power_law [ Pa.s^n ] ")
        (prefixvm(prefix,"carreau_law.lambda").c_str(), Feel::po::value< double >()->default_value( 3.313 ), "parameter lambda in carreau_law [ s ] ")
        (prefixvm(prefix,"carreau_law.n").c_str(), Feel::po::value< double >()->default_value( 0.3568 ), "parameter n in carreau_law ")
        (prefixvm(prefix,"carreau-yasuda_law.lambda").c_str(), Feel::po::value< double >()->default_value( 1.902 ), "parameter lambda in carreau-yasuda_law [ s ] ")
        (prefixvm(prefix,"carreau-yasuda_law.n").c_str(), Feel::po::value< double >()->default_value( 0.22 ), "parameter n in carreau-yasuda_law ")
        (prefixvm(prefix,"carreau-yasuda_law.a").c_str(), Feel::po::value< double >()->default_value( 1.25 ), "parameter a in carreau-yasuda_law ")
        // (prefixvm(prefix,"walburn-schneck_law.C1").c_str(), Feel::po::value< double >()->default_value( 0.00797 ), "parameter C1 in walburn-schneck_law ")
        // (prefixvm(prefix,"walburn-schneck_law.C2").c_str(), Feel::po::value< double >()->default_value( 0.0608 ), "parameter C2 in walburn-schneck_law ")
        // (prefixvm(prefix,"walburn-schneck_law.C3").c_str(), Feel::po::value< double >()->default_value( 0.00499 ), "parameter C3 in walburn-schneck_law ")
        // (prefixvm(prefix,"walburn-schneck_law.C4").c_str(), Feel::po::value< double >()->default_value( 14.585 ), "parameter C4 in walburn-schneck_law [l/g] ")
        // (prefixvm(prefix,"TPMA").c_str(), Feel::po::value< double >()->default_value( 25.9 ), "parameter TPMA (Total Proteins Minus Albumin) [ g/l ] ")
        ;

    return densityviscosityOptions;
}

Feel::po::options_description
fluidMechanics_options(std::string const& prefix)
{
    Feel::po::options_description fluidOptions("Fluid Mechanics options");
    fluidOptions.add_options()
        (prefixvm(prefix,"model").c_str(), Feel::po::value< std::string >(), "fluid model : Navier-Stokes,Stokes")
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "fluid solver")

        (prefixvm(prefix,"time-stepping").c_str(), Feel::po::value< std::string >()->default_value("BDF"), "time integration schema : BDF, Theta")
        (prefixvm(prefix,"time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")

        (prefixvm(prefix,"use-semi-implicit-time-scheme").c_str(), Feel::po::value<bool>()->default_value( false ), "use-semi-implicit-time-scheme")

        ( prefixvm(prefix,"start-by-solve-newtonian").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-newtonian")
        ( prefixvm(prefix,"start-by-solve-stokes-stationary").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-stokes-stationary")
        ( prefixvm(prefix,"start-by-solve-stokes-stationary.do-export").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-stokes-stationary.do-export")
        //(prefixvm(prefix,"strain_tensor.use-sym-tensor").c_str(), Feel::po::value< bool >()->default_value(true), "sym tensor or not ")

        (prefixvm(prefix,"stabilization-gls").c_str(), Feel::po::value<bool>()->default_value( false ), "apply stabilization method")
        (prefixvm(prefix,"stabilization-gls.type").c_str(), Feel::po::value<std::string>()->default_value( "gls" ), "supg,gls,gls-no-pspg, supg-pspg, pspg")
        (prefixvm(prefix,"stabilization-gls.parameter.method").c_str(), Feel::po::value<std::string>()->default_value( "eigenvalue" ), "method used for compute tau : eigenvalue, doubly-asymptotic-approximation")
        (prefixvm(prefix,"stabilization-gls.parameter.hsize.method").c_str(), Feel::po::value<std::string>()->default_value( "hmin" ), "hmin,h,meas")
        (prefixvm(prefix,"stabilization-gls.parameter.eigenvalue.penal-lambdaK").c_str(), Feel::po::value<double>()->default_value( 0. ), "apply stabilization method")
        (prefixvm(prefix,"stabilization-gls.convection-diffusion.location.expressions").c_str(), Feel::po::value<std::string>(), "expression which determinate the location where stab is applied")

        (prefixvm(prefix,"stabilization-gls.check-viscosity-dependency-on-coordinates").c_str(), Feel::po::value<bool>()->default_value( true ), "Debug option") // allow to not compute the derivativn of viscosity

        (prefixvm(prefix,"stabilisation-pspg").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-gls").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-convection-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-pressure").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-pressure-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-divergence").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-divergence-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        // (prefixvm(prefix,"stabilisation-divergence").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-div-div").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-div-div-beta").c_str(), Feel::po::value< double >()->default_value( 1.0 ), "parameter beta in div-div stab ")
        (prefixvm(prefix,"marked-zones.markedfaces").c_str(), po::value<std::vector<std::string> >()->multitoken(), "marked-zone.markedfaces" )
        (prefixvm(prefix,"marked-zones.elements-from-markedfaces").c_str(), po::value<std::vector<std::string> >()->multitoken(), "marked-zone.elements-from-markedfaces" )
        (prefixvm(prefix,"marked-zones.expressions").c_str(), po::value<std::vector<std::string> >()->multitoken(), "marked-zone.expressions" )
        (prefixvm(prefix,"marked-zones.internal-faces").c_str(), po::value<bool>()->default_value(false), "marked-zone.internal-faces" )
        (prefixvm(prefix,"define-pressure-cst").c_str(), Feel::po::value<bool>()->default_value( false ), "maybe need to define the pressure constant if only dirichlet bc")
        (prefixvm(prefix,"define-pressure-cst.markers").c_str(), po::value<std::vector<std::string> >()->multitoken(), "zone to applied pressure cst : marker1:marker2,marker3 -> pressure cst in marker1:marker2 and another pressure cst in marker3" )
        (prefixvm(prefix,"define-pressure-cst.method").c_str(), Feel::po::value<std::string>()->default_value( "algebraic"), " lagrange-multiplier or penalisation or algebraic")
        (prefixvm(prefix,"define-pressure-cst.penalisation-beta").c_str(), Feel::po::value< double >()->default_value( -10e-10 ), "parameter beta in cstpressure penalisation ")
        (prefixvm(prefix,"stabilisation-convection-energy").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"dirichletbc.type").c_str(), Feel::po::value<std::string>()->default_value( "elimination" ), "elimination, nitsche, lm")

        (prefixvm(prefix,"dirichletbc.lm.savemesh").c_str(), Feel::po::value<bool>()->default_value( false ), "export Lagrange multipliers mesh")
        (prefixvm(prefix,"dirichletbc.lm.use-submesh-relation").c_str(), Feel::po::value<bool>()->default_value( true ), "use submesh relation")
        (prefixvm(prefix,"dirichletbc.nitsche.gamma").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"bc-slip-form").c_str(), Feel::po::value<int>()->default_value( 1 ), "formulation for slip condition (1 or 2)")
        (prefixvm(prefix,"bc-slip-gammaN").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"bc-slip-gammaTau").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"hovisu").c_str(), Feel::po::value<bool>(), "true or false for high order visualisation")
        (prefixvm(prefix,"hovisu.space-used").c_str(), Feel::po::value<std::string>()->default_value( "velocity" ), "velocity , pressure, p1 ")

#if 0
        (prefixvm(prefix,"periodicity.translate-x").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-x")
        (prefixvm(prefix,"periodicity.translate-y").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-y")
        (prefixvm(prefix,"periodicity.translate-z").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-z")
        (prefixvm(prefix,"periodicity.marker1").c_str(), Feel::po::value<std::string>(), "periodicity.marker1 ")
        (prefixvm(prefix,"periodicity.marker2").c_str(), Feel::po::value<std::string>(), "periodicity.marker2 ")
        (prefixvm(prefix,"periodicity.pressure-jump").c_str(), Feel::po::value<double>()->default_value(1.0), "periodicity.pressure-jump ")
#endif
        (prefixvm(prefix,"blockns.type").c_str(), Feel::po::value<std::string>()->default_value("PCD"), "type : PCD,PMM")
        (prefixvm(prefix,"preconditioner.attach-mass-matrix").c_str(), Feel::po::value<bool>()->default_value(false), "attach mass matrix")
        (prefixvm(prefix,"preconditioner.attach-pmm").c_str(), Feel::po::value<bool>()->default_value(false), "attach pressure mass matrix")
        (prefixvm(prefix,"preconditioner.attach-pcd").c_str(), Feel::po::value<bool>()->default_value(false), "attach operator pcd")

        (prefixvm(prefix,"use-velocity-near-null-space").c_str(), Feel::po::value<bool>()->default_value( false ), "use-velocity-near-null-space")
        (prefixvm(prefix,"use-velocity-near-null-space.prefix").c_str(), Feel::po::value<std::string>(), "use-velocity-near-null-space.prefix")

        (prefixvm(prefix,"fluid-outlet.type").c_str(), Feel::po::value<std::string>()->default_value( "free" ), "type : free, windkessel ")
        (prefixvm(prefix,"fluid-outlet.windkessel.coupling").c_str(), Feel::po::value<std::string>()->default_value( "implicit" ), "explicit, implicit ")

        (prefixvm(prefix,"use-gravity-force").c_str(), Feel::po::value<bool>()->default_value(false), "use-gravity-force")
        (prefixvm(prefix,"gravity-force").c_str(), Feel::po::value<std::string>(), "gravity-force : (default is {0,-9.80665} or {0,0,-9.80665}")

        (prefixvm(prefix,"distance-to-wall.enabled").c_str(), Feel::po::value<bool>()->default_value(false), "enable distance-to-wall computation")
        (prefixvm(prefix,"distance-to-wall.markers").c_str(), po::value<std::vector<std::string> >()->multitoken(), "markers used for compute distance-to-wall" )

        (prefixvm(prefix,"use-semi-implicit-turbulence-coupling").c_str(), Feel::po::value<bool>()->default_value( false ), "use-semi-implicit-turbulence-coupling")

        (prefixvm(prefix,"pcd.apply-homogeneous-dirichlet-in-newton").c_str(), Feel::po::value<bool>()->default_value(false), "use-gravity-force")

        (prefixvm(prefix,"body.articulation.method").c_str(), Feel::po::value<std::string>()->default_value("lm"), "lm or p-matrix")

        (prefixvm(prefix,"body.elastic-behavior.use-old-version").c_str(), Feel::po::value<bool>()->default_value(false), "use-gravity-force")
        ;

    fluidOptions
        .add( modelnumerical_options( prefix ) )
        .add( bdf_options( prefix ) )
        .add( ts_options( prefix ) )
        .add( alemesh_options( prefix ) )
        .add( backend_options( prefixvm(prefix,"fluidinlet") ) )
        .add( densityviscosity_options( prefix ) )
        .add( pcd_options( prefix ) )
        .add( coefficientformpdes_options( prefixvm(prefix,"turbulence") ) )

        .add( modelnumerical_options( prefixvm(prefix,"body") ) )
        ;

    return fluidOptions;
}

Feel::po::options_description
solidMechanics_options(std::string const& prefix)
{
    Feel::po::options_description solidOptions("Solid Mechanics options");
    solidOptions.add_options()
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "density")
        (prefixvm(prefix,"youngmodulus").c_str(), Feel::po::value<double>()->default_value( 3.e6 ), "young modulus")
        (prefixvm(prefix,"coeffpoisson").c_str(), Feel::po::value<double>()->default_value( 0.3 ), "poisson coefficient")
        (prefixvm(prefix,"model").c_str(), Feel::po::value< std::string >()/*->default_value("Elasticity")*/, "struct model")
        (prefixvm(prefix,"material_law").c_str(), Feel::po::value< std::string >()->default_value("StVenantKirchhoff"), "StVenantKirchhoff, NeoHookean")
        //(prefixvm(prefix,"mechanicalproperties.compressible.volumic-law").c_str(), Feel::po::value< std::string >()->default_value("classic"), "classic, simo1985")
        //(prefixvm(prefix,"mechanicalproperties.compressible.neohookean.variant").c_str(),
        //Feel::po::value< std::string >()->default_value("default"), "default, molecular-theory, molecular-theory-simo1985")
        (prefixvm(prefix,"formulation").c_str(), Feel::po::value<std::string>()->default_value( "displacement" ), "displacement,displacement-pressure")
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "struct solver")
        (prefixvm(prefix,"time-stepping").c_str(), Feel::po::value< std::string >()->default_value("Newmark"), "time integration schema : Newmark, BDF, Theta")
        (prefixvm(prefix,"time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")

        //(prefixvm(prefix,"time-rho").c_str(), Feel::po::value< double >()->default_value(0.8), " Generalized-Alpha parameter")

        (prefixvm(prefix,"time-initial.displacement.files.directory").c_str(), Feel::po::value<std::string>(), "initial displacemen")
        (prefixvm(prefix,"time-initial.displacement.files.format").c_str(), Feel::po::value<std::string>()->default_value( "hdf5" ), "intial displacement file format")

        (prefixvm(prefix,"1dreduced-geofile").c_str(), Feel::po::value< std::string >(), "input geo file 1dreduced")
        (prefixvm(prefix,"1dreduced-thickness").c_str(), Feel::po::value< double >()->default_value(0.1), "1dreduced-thickness")
        (prefixvm(prefix,"1dreduced-radius").c_str(), Feel::po::value< double >()->default_value(0.5), "1dreduced-radius")

        (prefixvm(prefix,"hovisu").c_str(), Feel::po::value<bool>(), "true or false for high order visualisation")
        (prefixvm(prefix,"hovisu.space-used").c_str(), Feel::po::value<std::string>()->default_value( "displacement" ), "displacement, pressure, p1 ")

        (prefixvm(prefix,"use-null-space").c_str(), Feel::po::value<bool>()->default_value( false ), "use-null-space")
        (prefixvm(prefix,"use-near-null-space").c_str(), Feel::po::value<bool>()->default_value( true ), "use-near-null-space")
        ;
    solidOptions.add( gmsh_options( prefixvm(prefix,"1dreduced") ) );
    return solidOptions.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}

Feel::po::options_description
fluidStructInteraction_options( std::string const& prefix )
{
    Feel::po::options_description FSIoptions("FSI options");
    FSIoptions.add_options()
        (prefixvm(prefix,"hsize").c_str(), Feel::po::value<double>()->default_value( 0.1 ), "characteristic mesh size")
        (prefixvm(prefix,"fluid-mesh.markers").c_str(), po::value<std::vector<std::string> >()->multitoken(), "solid-mesh.markers" )
        (prefixvm(prefix,"solid-mesh.markers").c_str(), po::value<std::vector<std::string> >()->multitoken(), "solid-mesh.markers" )
        (prefixvm(prefix,"solid-mesh.extract-1d-from-fluid-mesh").c_str(), Feel::po::value<bool>()->default_value( false ), "solid-mesh.extract-1d-from-fluid-mesh")
        (prefixvm(prefix,"mesh-save.tag").c_str(),Feel::po::value< std::string >()->default_value(""), "mesh-tag")
        (prefixvm(prefix,"mesh-save.directory").c_str(),Feel::po::value< std::string >(), "mesh-directory")
        (prefixvm(prefix,"mesh-save.force-rebuild").c_str(), Feel::po::value<bool>()->default_value( false ), "mesh-save.force-rebuild")

        (prefixvm(prefix,"conforming-interface").c_str(), Feel::po::value<bool>()->default_value( false ), " fsi interface is conforme?")
        (prefixvm(prefix,"coupling-type").c_str(),Feel::po::value< std::string >()->default_value("Implicit"), " Implicit or Semi-Implicit")
        (prefixvm(prefix,"coupling-bc").c_str(),Feel::po::value< std::string >()->default_value("dirichlet-neumann"), " dirichlet-neumann, robin-robin,robin-neumann")

        (prefixvm(prefix,"fixpoint.tol").c_str(), Feel::po::value<double>()->default_value( 1.e-6 ), "tolerance pt fixe")
        (prefixvm(prefix,"fixpoint.initialtheta").c_str(), Feel::po::value<double>()->default_value( 1. ), "relax aitken parameter")
        (prefixvm(prefix,"fixpoint.min_theta").c_str(), Feel::po::value<double>()->default_value( 1.e-4 ), "min theta parameter")
        (prefixvm(prefix,"fixpoint.maxit").c_str(), Feel::po::value<int>()->default_value( 1000 ), "max iteration")
        (prefixvm(prefix,"fixpoint.minit-convergence").c_str(), Feel::po::value<int>()->default_value( 3 ), "max iteration")
        // additionals options if reuse-prec or reuse-jac are actived
        (prefixvm(prefix,"fluid.reuse-prec.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi fluid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"solid.reuse-prec.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi solid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"fluid.reuse-jac.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi fluid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"solid.reuse-jac.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi solid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"fluid.reuse-prec.activated-after-n-fsi-it").c_str(), Feel::po::value<int>()->default_value( 0 ), "activated-after-n-fsi-it")
        (prefixvm(prefix,"solid.reuse-prec.activated-after-n-fsi-it").c_str(), Feel::po::value<int>()->default_value( 0 ), "activated-after-n-fsi-it")
        (prefixvm(prefix,"fluid.reuse-prec.activated-only-if-greater-than-tol").c_str(), Feel::po::value<double>()->default_value( 1e8 ), "activated-only-if-greater-than-tol")
        (prefixvm(prefix,"solid.reuse-prec.activated-only-if-greater-than-tol").c_str(), Feel::po::value<double>()->default_value( 1e8 ), "activated-only-if-greater-than-tol")
        // coupling-nitsche-family : nitsche,robin-robin,robin-robin-genuine
        (prefixvm(prefix,"coupling-nitsche-family.gamma").c_str(), Feel::po::value<double>()->default_value( 2500 ), "nitsche parameters")
        (prefixvm(prefix,"coupling-nitsche-family.gamma0").c_str(), Feel::po::value<double>()->default_value( 1 ), "nitsche parameters")
        (prefixvm(prefix,"coupling-nitsche-family.alpha").c_str(), Feel::po::value<double>()->default_value( 1 ), "nitsche parameters")
        (prefixvm(prefix,"coupling-nitsche-family.use-aitken").c_str(), Feel::po::value<bool>()->default_value( false ), "use-aitken")
        // coupling-robin-neumann-generalized
        (prefixvm(prefix,"coupling-robin-neumann-generalized.use-mass-matrix-lumped-in-solid").c_str(), Feel::po::value<bool>()->default_value( true ),"use-mass-matrix-lumped-in-solid") 
        (prefixvm(prefix,"coupling-robin-neumann-generalized.use-operator-constant").c_str(), Feel::po::value<double>(),"use-operator-constant")
        (prefixvm(prefix,"coupling-robin-neumann-generalized.use-operator-proportional-to-identity").c_str(), Feel::po::value<bool>()->default_value( false ),"use-operator-proportional-to-identity")
        (prefixvm(prefix,"coupling-robin-neumann-generalized.use-aitken").c_str(), Feel::po::value<bool>()->default_value( false ), "use-aitken")
        (prefixvm(prefix,"coupling-robin-neumann-generalized.strategy-time-step-compatibility").c_str(), Feel::po::value<std::string>()->default_value( "default" ),"strategy-time-step-compatibility")

        (prefixvm(prefix,"transfert-velocity-F2S.use-extrapolation").c_str(), Feel::po::value<bool>()->default_value( true ), "transfert-velocity-F2S.use-extrapolation")
        ;
    return FSIoptions.add( modelnumerical_options(prefix) );
}


Feel::po::options_description
heat_options(std::string const& prefix)
{
    Feel::po::options_description heatOptions("Heat options");
    heatOptions.add_options()
        (prefixvm(prefix,"thermal-conductivity").c_str(), Feel::po::value<double>()->default_value( 1 ), "thermal-conductivity [ W/(m*K) ]")
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1 ), "density [ kg/(m^3) ]")
        (prefixvm(prefix,"heat-capacity").c_str(), Feel::po::value<double>()->default_value( 1 ), "heat-capacity [ J/(kg*K) ]")
        (prefixvm(prefix,"thermal-expansion").c_str(), Feel::po::value<double>()->default_value( 1e-4 ), "thermal-conductivity [ 1/K) ]")
        // (prefixvm(prefix,"use_velocity-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "use-velocity-convection")
        // (prefixvm(prefix,"velocity-convection_is_incompressible").c_str(), Feel::po::value<bool>()->default_value( false ), "velocity-convection-is-incompressible")
        // (prefixvm(prefix,"velocity-convection").c_str(), Feel::po::value<std::string>(), "math expression")
        (prefixvm(prefix,"initial-solution.temperature").c_str(), Feel::po::value<std::string>(), "math expression")

        (prefixvm(prefix,"use-extended-doftable").c_str(), Feel::po::value<bool>()->default_value( false ), "use-extended-doftable")

        (prefixvm(prefix,"stabilization-gls").c_str(), Feel::po::value<bool>()->default_value( false ), "apply stabilization method")
        (prefixvm(prefix,"stabilization-gls.type").c_str(), Feel::po::value<std::string>()->default_value( "gls" ), "supg,gls,unusual-gls")
        (prefixvm(prefix,"stabilization-gls.parameter.method").c_str(), Feel::po::value<std::string>()->default_value( "eigenvalue" ), "method used for compute tau : eigenvalue, doubly-asymptotic-approximation")
        (prefixvm(prefix,"stabilization-gls.parameter.hsize.method").c_str(), Feel::po::value<std::string>()->default_value( "hmin" ), "hmin,h,meas")
        (prefixvm(prefix,"stabilization-gls.parameter.eigenvalue.penal-lambdaK").c_str(), Feel::po::value<double>()->default_value( 0. ), "apply stabilization method")

        (prefixvm(prefix,"stabilization-gls.check-conductivity-dependency-on-coordinates").c_str(), Feel::po::value<bool>()->default_value( true ), "Debug option") // allow to not compute the derivativn of thermal conductivity

        (prefixvm(prefix,"time-stepping").c_str(), Feel::po::value< std::string >()->default_value("BDF"), "time integration schema : BDF, Theta")
        (prefixvm(prefix,"time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")

        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "numeric solver : automatic, Newton, Picard, Linear")
        ;
    return heatOptions.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}
Feel::po::options_description
electricity_options(std::string const& prefix)
{
    Feel::po::options_description electricityOptions("Electricity options");
    electricityOptions.add_options()
        (prefixvm(prefix,"electric-conductivity").c_str(), Feel::po::value<double>()->default_value( 1 ), "electric-conductivity")
        ;
    return electricityOptions.add( modelnumerical_options( prefix ) );
}
Feel::po::options_description
maxwell_options(std::string const& prefix)
{
    Feel::po::options_description maxwellOptions("Maxwell options");
    maxwellOptions.add_options()
        (prefixvm(prefix,"magnetic-permeability").c_str(), Feel::po::value<double>()->default_value( 1 ), "magnetic-permeability")
        (prefixvm(prefix,"regularization-epsilon").c_str(), Feel::po::value<double>()->default_value( 1e-5), "regularization parameter")
        ;
    return maxwellOptions.add( modelnumerical_options( prefix ) );
}
Feel::po::options_description
thermoElectric_options(std::string const& prefix)
{
    Feel::po::options_description thermoElectricOptions("ThermoElectric options");

    thermoElectricOptions.add_options()
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "thermoelectric solver : automatic, Newton, Picard")
        (prefixvm(prefix,"solver-newton.initial-guess.use-linear-thermo-electric").c_str(), Feel::po::value<bool>()->default_value( false ), "solver-newton.initial-guess.use-linear-thermo-electric")
        (prefixvm(prefix,"solver-newton.initial-guess.use-linear-heat").c_str(), Feel::po::value<bool>()->default_value( false ), "solver-newton.initial-guess.use-linear-heat")
        (prefixvm(prefix,"solver-newton.initial-guess.use-linear-electric").c_str(), Feel::po::value<bool>()->default_value( false ), "solver-newton.initial-guess.use-linear-electric")
        ;
    thermoElectricOptions.add( heat_options( prefixvm(prefix,"heat") ) );
    thermoElectricOptions.add( electricity_options( prefixvm(prefix,"electric") ) );
    return thermoElectricOptions.add( modelnumerical_options( prefix ) );
}

Feel::po::options_description
heatFluid_options(std::string const& prefix)
{
    Feel::po::options_description heatFluidOptions("HeatFluid options");
    heatFluidOptions.add_options()
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "numeric solver : automatic, Newton, Picard")
        (prefixvm(prefix,"use-natural-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "use natural convection")
        (prefixvm(prefix,"Boussinesq.ref-temperature").c_str(), Feel::po::value<double>()->default_value( 300. ), "Boussinesq ref-temperature T0")
        (prefixvm(prefix,"gravity-force").c_str(), Feel::po::value<std::string>(), "gravity-force : (default is {0,-9.80665} or {0,0,-9.80665}")
        (prefixvm(prefix,"use-semi-implicit-time-scheme").c_str(), Feel::po::value<bool>()->default_value( false ), "use-semi-implicit-time-scheme")
    ;
    heatFluidOptions.add( heat_options( prefixvm(prefix,"heat") ) );
    heatFluidOptions.add( fluidMechanics_options( prefixvm(prefix,"fluid") ) );
    return heatFluidOptions.add( modelnumerical_options( prefix ) );
}

Feel::po::options_description
advection_options(std::string const& prefix)
{
    Feel::po::options_description advectionOptions("Advection options");
    advectionOptions.add_options()
        (prefixvm(prefix,"model").c_str(), Feel::po::value< std::string >(), "advection model : Advection, Advection-Diffusion, Advection-Diffusion-Reaction")
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >()->default_value( "automatic" ), "advection solver")
        (prefixvm(prefix,"D").c_str(), Feel::po::value<double>()->default_value( 0 ), "diffusion coefficient [ m^2.s^-1 ]")
        (prefixvm(prefix,"R").c_str(), Feel::po::value<double>()->default_value( 0 ), "reaction coefficient [ s^-1 ]")
        (prefixvm(prefix,"advection-velocity").c_str(), Feel::po::value<std::string>(), "math expression")
        (prefixvm(prefix,"stabilization.method").c_str(), Feel::po::value<std::string>()->default_value( "GALS" ), "stabilization method")
        (prefixvm(prefix,"stabilization-cip.coeff").c_str(), Feel::po::value<double>()->default_value( 1 ), "CIP stabilization coefficient")
        (prefixvm(prefix,"export-all").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_all")
        (prefixvm(prefix,"export-advection-velocity").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_advection_velocity")
        (prefixvm(prefix,"export-diffusion").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_diffusion_coeff")
        (prefixvm(prefix,"export-reaction").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_reaction_coeff")
        (prefixvm(prefix,"export-source").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_source_field")

        (prefixvm(prefix,"stabilization-gls.parameter.method").c_str(), Feel::po::value<std::string>()->default_value( "eigenvalue" ), "method used for compute tau : eigenvalue, doubly-asymptotic-approximation")
        (prefixvm(prefix,"stabilization-gls.parameter.hsize.method").c_str(), Feel::po::value<std::string>()->default_value( "hmin" ), "hmin,h,meas")
        (prefixvm(prefix,"stabilization-gls.parameter.eigenvalue.penal-lambdaK").c_str(), Feel::po::value<double>()->default_value( 0. ), "apply stabilization method")

        ;
    return advectionOptions.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}

Feel::po::options_description
redistanciation_fm_options(std::string const& prefix, bool addProjectorsOpts )
{
    Feel::po::options_description redistanciationFMOptions("RedistanciationFM options");
    redistanciationFMOptions.add_options()
        (prefixvm(prefix,"fm-init-method").c_str(), Feel::po::value<std::string>()->default_value("ilp-nodal"), "strategy to initialise the first elements before the fast marching:\nnone = do nothing\nilp-[nodal/l2/smooth] = interface local projection by nodal, L2 or smooth (|grad phi|) projections\nhj = Hamilton Jacoby equation (with parameters given in options)")
        ;

    if( addProjectorsOpts )
    {
        redistanciationFMOptions
            .add( backend_options( prefixvm(prefix, "projector-l2") ) )
            .add( backend_options( prefixvm(prefix, "projector-sm") ) )
            ;
        redistanciationFMOptions.add_options()
            (prefixvm(prefix,"projector-sm.smooth-coeff").c_str(), Feel::po::value<double>()->default_value(0.1), "smoothing coefficient for projector-sm")
            ;
    }

    return redistanciationFMOptions;
}

Feel::po::options_description
redistanciation_hj_options(std::string const& prefix)
{
    Feel::po::options_description redistanciationHJOptions("RedistanciationHJ options");
    redistanciationHJOptions.add_options()
        (prefixvm(prefix,"tol").c_str(), Feel::po::value<double>()->default_value( 0.03 ), "tolerance on residual to \"distance function\" of HJ redistanciated level set")
        (prefixvm(prefix,"time-step").c_str(), Feel::po::value<double>()->default_value( 0.1 ), "time step used in HJ equation")
        (prefixvm(prefix,"max-iter").c_str(), Feel::po::value<int>()->default_value( 15 ), "maximum number of iterations for Hamilton-Jacobi redistanciation")
        (prefixvm(prefix,"thickness-heaviside").c_str(), Feel::po::value<double>()->default_value( 0.1 ), "thickness of the interface (support for Heaviside used to compute sign function)")
        (prefixvm(prefix,"keep-volume").c_str(), Feel::po::value<bool>()->default_value( true ), "use constraint to conserve levelset volume")
        ;

    redistanciationHJOptions.add( advection_options( prefix ) );

    return redistanciationHJOptions;
}

Feel::po::options_description
levelset_options(std::string const& prefix)
{
    Feel::po::options_description levelsetOptions("Levelset options");
    levelsetOptions.add_options()
        //(prefixvm(prefix,"fm-use-markerdirac").c_str(), Feel::po::value<bool>()->default_value( false ), "use markerDirac to mark initially done elements in fast-marching")
        (prefixvm(prefix,"use-regularized-phi").c_str(), Feel::po::value<bool>()->default_value( false ), "use grad(phi)/|grad(phi)| to evaluate Dirac and Heaviside functions")
        (prefixvm(prefix,"h-d-nodal-proj").c_str(), Feel::po::value<bool>()->default_value( true ), "use nodal projection to compute dirac and heaviside functions (if false, use L2 projection)")
        (prefixvm(prefix,"thickness-interface").c_str(), Feel::po::value<double>(), "thickness of the interface (support for Dirac and Heaviside functions)")
        (prefixvm(prefix,"use-adaptive-thickness").c_str(), Feel::po::value<bool>()->default_value( false ), "automatically adapt the thickness of Dirac and Heaviside")
        (prefixvm(prefix,"thickness-interface-rectangular-function").c_str(), Feel::po::value<double>(), "thickness of the interface rectangular function")
        (prefixvm(prefix,"distance-method").c_str(), Feel::po::value<std::string>()->default_value( "fm" ), "levelset distance computation method (none, fm: fast-marching, hj: hamilton-jacobi, renormalisation)")
        (prefixvm(prefix,"redist-method").c_str(), Feel::po::value<std::string>()->default_value( "fm" ), "levelset redistanciation method (none, fm: fast-marching, hj: hamilton-jacobi, renormalisation)")
        (prefixvm(prefix,"use-order1-after-redist").c_str(), Feel::po::value<bool>()->default_value( false ), "Use order 1 time-stepper after redistanciation.")

        (prefixvm(prefix,"redist-initial-value").c_str(), Feel::po::value<bool>()->default_value( false ), "redistanciate levelset after setting initial value")

        (prefixvm(prefix,"gradphi-method").c_str(), Feel::po::value<std::string>()->default_value( "nodal-projection" ), "method to compute gradphi (nodal-projection, l2-projection, smooth-projection, pn-nodal-projection)")
        (prefixvm(prefix,"modgradphi-method").c_str(), Feel::po::value<std::string>()->default_value( "nodal-projection" ), "method to compute gradphi (nodal-projection, l2-projection, smooth-projection, pn-nodal-projection)")
        (prefixvm(prefix,"curvature-method").c_str(), Feel::po::value<std::string>()->default_value( "smooth-projection" ), "method to compute curvature (nodal-projection, l2-projection, smooth-projection, pn-nodal-projection)")
        (prefixvm(prefix,"curvature-diffusion.time-step").c_str(), Feel::po::value<double>()->default_value( 0.01 ), "time step used for the heat equations in diffusion-order1 and diffusion-order2 curvature methods")
        (prefixvm(prefix,"curvature-diffusion.time-discretisation").c_str(), Feel::po::value<std::string>()->default_value( "crank-nicolson" ), "time discretisation scheme used for the heat equations in diffusion-order1 and diffusion-order2 curvature methods")

        (prefixvm(prefix,"projector-sm-scalar.smooth-coeff").c_str(), Feel::po::value<double>()->default_value(0.1), "smoothing coefficient for projector-sm-scalar")
        (prefixvm(prefix,"projector-sm-vectorial.smooth-coeff").c_str(), Feel::po::value<double>()->default_value(0.1), "smoothing coefficient for projector-sm-vectorial")
        (prefixvm(prefix,"projector-sm-scalar-isopn.smooth-coeff").c_str(), Feel::po::value<double>()->default_value(0.1), "smoothing coefficient for projector-sm-scalar-isopn")
        (prefixvm(prefix,"projector-sm-vectorial-isopn.smooth-coeff").c_str(), Feel::po::value<double>()->default_value(0.1), "smoothing coefficient for projector-sm-vectorial-isopn")

        (prefixvm(prefix,"use-gradient-augmented").c_str(), Feel::po::value<bool>()->default_value(false), "Advect modGradPhi independently")
        (prefixvm(prefix,"reinit-gradient-augmented").c_str(), Feel::po::value<bool>()->default_value(false), "Reinit modGradPhi when phi is redistanciated")

        (prefixvm(prefix,"use-stretch-augmented").c_str(), Feel::po::value<bool>()->default_value(false), "Advect stretch independently")
        (prefixvm(prefix,"reinit-stretch-augmented").c_str(), Feel::po::value<bool>()->default_value(false), "Reinit stretch when phi is redistanciated")

        (prefixvm(prefix,"use-cauchy-augmented").c_str(), Feel::po::value<bool>()->default_value(false), "Advect additional backward characteristics to compute Cauchy tensor")
        (prefixvm(prefix,"initial-backward-characteristics").c_str(), Feel::po::value<std::string>(), "Initial  backward characteristics value (default for material at rest is {x,y(,z)})")

        (prefixvm(prefix,"use-space-iso-pn").c_str(), Feel::po::value<bool>()->default_value(false), "use isoPN spaces for the levelset")

        (prefixvm(prefix,"fix-volume").c_str(), Feel::po::value<bool>()->default_value(false), "correct levelset volume after each iteration (using phi->phi+(V-Vi)/L)")
        (prefixvm(prefix,"fix-area").c_str(), Feel::po::value<bool>()->default_value(false), "correct levelset area after each iteration (using phi->phi+mu*(K-Kbar) with Kbar=int(K*delta)/L and mu=(L-L0)/int(K*(K-Kbar)))")

        (prefixvm(prefix,"use-extension-velocity").c_str(), Feel::po::value<bool>()->default_value(false), "use extension velocity which preserves the distance property when advecting the level set")
        (prefixvm(prefix,"extension-velocity.gamma").c_str(), Feel::po::value<double>()->default_value(10), "value for the gamma premultying Nitsche's term to impose the weak BC at the interface")

        (prefixvm(prefix,"do_export_advection").c_str(), Feel::po::value<bool>()->default_value(false), "doExportAdvection")
        (prefixvm(prefix,"do_export_dirac").c_str(), Feel::po::value<bool>(), "doExportDirac")
        (prefixvm(prefix,"do_export_heaviside").c_str(), Feel::po::value<bool>(), "doExportHeaviside")
        (prefixvm(prefix,"do_export_normal").c_str(), Feel::po::value<bool>(), "doExportNormal")
        (prefixvm(prefix,"do_export_curvature").c_str(), Feel::po::value<bool>(), "doExportCurvature")
        (prefixvm(prefix,"do_export_gradphi").c_str(), Feel::po::value<bool>(), "doExportGradPhi")
        (prefixvm(prefix,"do_export_modgradphi").c_str(), Feel::po::value<bool>(), "doExportModGradPhi")
        (prefixvm(prefix,"do_export_dirac").c_str(), Feel::po::value<bool>(), "doExportDirac")
        (prefixvm(prefix,"do_export_heaviside").c_str(), Feel::po::value<bool>(), "doExportHeaviside")
        (prefixvm(prefix,"do_export_normal").c_str(), Feel::po::value<bool>(), "doExportNormal")
        (prefixvm(prefix,"do_export_curvature").c_str(), Feel::po::value<bool>(), "doExportCurvature")
        (prefixvm(prefix,"do_export_advectionvelocity").c_str(), Feel::po::value<bool>(), "doExportAdvectionVelocity")
        (prefixvm(prefix,"do_export_modgradphi-advection").c_str(), Feel::po::value<bool>()->default_value(false), "doExportModGradPhi-Advection")
        (prefixvm(prefix,"do_export_stretch-advection").c_str(), Feel::po::value<bool>()->default_value(false), "doExportStretch-Advection")
        (prefixvm(prefix,"do_export_backward-characteristics-advection").c_str(), Feel::po::value<bool>()->default_value(false), "doExportBackwardCharacteristics-Advection")
        (prefixvm(prefix,"do_export_backwardcharacteristics").c_str(), Feel::po::value<bool>(), "doExportCauchyGreenInvariant1")
        (prefixvm(prefix,"do_export_cauchygreeninvariant1").c_str(), Feel::po::value<bool>(), "doExportCauchyGreenInvariant1")
        (prefixvm(prefix,"do_export_cauchygreeninvariant2").c_str(), Feel::po::value<bool>(), "doExportCauchyGreenInvariant2")
        ;

    levelsetOptions
        .add( advection_options( prefix ) )
        .add( advection_options( prefixvm(prefix, "modgradphi-advection") ) )
        .add( advection_options( prefixvm(prefix, "stretch-advection") ) )
        .add( advection_options( prefixvm(prefix, "backward-characteristics-advection") ) )
        .add( backend_options( prefixvm(prefix, "curvature-diffusion") ) )
        .add( backend_options( prefixvm(prefix, "projector-l2-scalar") ) )
        .add( backend_options( prefixvm(prefix, "projector-l2-vectorial") ) )
        .add( backend_options( prefixvm(prefix, "projector-l2-tensor2symm") ) )
        .add( backend_options( prefixvm(prefix, "projector-l2-scalar-isopn") ) )
        .add( backend_options( prefixvm(prefix, "projector-l2-vectorial-isopn") ) )
        .add( backend_options( prefixvm(prefix, "projector-sm-scalar") ) )
        .add( backend_options( prefixvm(prefix, "projector-sm-vectorial") ) )
        .add( backend_options( prefixvm(prefix, "projector-sm-tensor2symm") ) )
        .add( backend_options( prefixvm(prefix, "projector-sm-scalar-isopn") ) )
        .add( backend_options( prefixvm(prefix, "projector-sm-vectorial-isopn") ) )
        .add( backend_options( prefixvm(prefix, "extension-velocity") ) )
        .add( redistanciation_fm_options( prefixvm(prefix, "redist-fm"), false ) )
        .add( redistanciation_hj_options( prefixvm(prefix, "redist-hj") ) )
        ;

    return levelsetOptions;
}

Feel::po::options_description
interfaceforces_options(std::string const& prefix)
{
    Feel::po::options_description interfaceForcesOptions("InterfaceForces options");
    interfaceForcesOptions.add_options()
        (prefixvm(prefix,"surface-tension.coeff").c_str(), Feel::po::value<double>()->default_value(0.), "Surface tension coefficient" )

        (prefixvm(prefix,"helfrich-bending-modulus").c_str(), Feel::po::value<double>()->default_value(0.), "Helfrich bending modulus k_B" )
        (prefixvm(prefix,"helfrich-force-impl").c_str(), Feel::po::value<int>()->default_value(0), "Implementation of Helfrich force" )
        (prefixvm(prefix,"helfrich-inner-div-impl").c_str(), Feel::po::value<std::string>()->default_value("generic"), "Implementation of Helfrich inner div (generic, distance)" )
        (prefixvm(prefix,"helfrich-force.use-integration-by-parts").c_str(), Feel::po::value<bool>()->default_value(false), "Integrate by parts the div term in Helfrich force" )
        (prefixvm(prefix,"helfrich-force.projection-method").c_str(), Feel::po::value<std::string>()->default_value("nodal-projection"), "Projection method for Helfrich force (nodal-projection, l2-projection, smooth-projection)" )

        (prefixvm(prefix,"inextensibility-force-coeff").c_str(), Feel::po::value<double>()->default_value(1), "Inextensibility force coefficient (Lambda)" )
        (prefixvm(prefix,"inextensibility-force-epsilon").c_str(), Feel::po::value<double>(), "Inextensibility force epsilon" )

        (prefixvm(prefix,"elastic-stretch-modulus").c_str(), Feel::po::value<double>()->default_value(1), "Linear elastic force stretch modulus" )
        (prefixvm(prefix,"elastic-shear-modulus").c_str(), Feel::po::value<double>()->default_value(1), "Linear elastic force shear modulus" )

        (prefixvm(prefix,"skalak-stretch-modulus").c_str(), Feel::po::value<double>()->default_value(1), "Skalak force stretch modulus" )
        (prefixvm(prefix,"skalak-shear-modulus").c_str(), Feel::po::value<double>()->default_value(1), "Skalak force shear modulus" )
        ;

    return interfaceForcesOptions;
}

Feel::po::options_description
multifluid_options(std::string const& prefix, unsigned int nls)
{
    Feel::po::options_description multifluidOptions("MultiFluid options");
    multifluidOptions.add_options()
        (prefixvm(prefix,"nlevelsets").c_str(), Feel::po::value<int>()->default_value( 1 ), "number of levelsets")

        (prefixvm(prefix, "use-picard-iterations").c_str(), Feel::po::value<bool>()->default_value( false ), "solve NS-LS coupling with non-linear Picard iterations")

        (prefixvm(prefix, "use-ls-P1iso-mesh").c_str(), Feel::po::value<bool>()->default_value( false ), "use P1 iso-Pn(u) mesh for the levelset")
        ;

    multifluidOptions
        .add( modelnumerical_options( prefix ) )
        .add( fluidMechanics_options( prefixvm(prefix, "fluid") ) )
        .add( levelset_options( prefixvm(prefix, "levelset") ) )
        ;
    for( uint16_type n = 0; n < nls; ++n )
    {
        std::string levelset_prefix = prefixvm(prefix, (boost::format( "levelset%1%" ) %(n)).str());
        multifluidOptions.add( levelset_options( levelset_prefix ) );
        multifluidOptions.add( densityviscosity_options( levelset_prefix ) );
        multifluidOptions.add( interfaceforces_options( levelset_prefix ) );
        multifluidOptions.add_options()
            // Reinitialization
            (prefixvm(levelset_prefix,"redist-every").c_str(), Feel::po::value<int>()->default_value( 10 ), "redistanciate levelset every n iterations" )
            // Interface forces model
            (prefixvm(levelset_prefix,"interface-forces-model").c_str(), Feel::po::value<std::vector<std::string>>()->multitoken()->composing(), "models for interface forces (helfrich, ...)" )
            // Inextensibility
            (prefixvm(levelset_prefix, "enable-inextensibility").c_str(), Feel::po::value<bool>()->default_value( false ), "enable inextensibility of level set")
            (prefixvm(levelset_prefix, "inextensibility-method").c_str(), Feel::po::value<std::string>()->default_value( "penalty" ), "method to impose level set inextensibility (penalty or lagrange-multiplier)")
            (prefixvm(levelset_prefix, "inextensibility-gamma").c_str(), Feel::po::value<double>()->default_value( 5. ), "coeff for inextensibility by penalty method ")
            ;
    }

    return multifluidOptions;
}

Feel::po::options_description
alemesh_options(std::string const& prefix)
{
    po::options_description desc_options("alemesh options");
    desc_options.add_options()
        (prefixvm(prefix,"alemesh.type").c_str(), Feel::po::value<std::string>()->default_value( "harmonic" ), "type : harmonic, winslow")
        //(prefixvm(prefix,"alemesh.verbose").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose")
        //(prefixvm(prefix,"alemesh.verbose_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        (prefixvm(prefix,"alemesh.verbose_solvertimer").c_str(), Feel::po::value<bool>()->default_value( true ), "true or false to view verbose")
        (prefixvm(prefix,"alemesh.verbose_solvertimer_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        (prefixvm(prefix,"alemesh.winslow.solver").c_str(), Feel::po::value<std::string>()->default_value( "Picard" ), "solver : Picard, Newton, Picard-Newton")
        (prefixvm(prefix,"alemesh.winslow.mesh-adaptation").c_str(), Feel::po::value<bool>()->default_value( false ), "use mesh adaptation")
        (prefixvm(prefix,"alemesh.winslow.mesh-adaptation.scalar-weight").c_str(), Feel::po::value<bool>()->default_value( true ), "use scalar weight with mesh adaptation")
        (prefixvm(prefix,"alemesh.winslow.mesh-adaptation.scalar-weight.expr").c_str(), Feel::po::value<std::string>(), "scalar weight expression with mesh adaptation")
        (prefixvm(prefix,"alemesh.winslow.Picard-Newton.maxit-Picard").c_str(), Feel::po::value<int>()->default_value( 3 ), "nb max iteration of point fixe")
        (prefixvm(prefix,"alemesh.harmonic.use_adaptive_penalisation").c_str(), Feel::po::value<bool>()->default_value( true ), "use tau penalisation in low formulation")
        (prefixvm(prefix,"alemesh.apply-ho-correction").c_str(), Feel::po::value<bool>()->default_value( true ), "do apply ho correction")
        //(prefixvm(prefix,"alemesh.export").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to export")
        ;
    return desc_options
        .add( bdf_options( prefixvm(prefix,"alemesh") ) )
        .add( ts_options( prefixvm(prefix,"alemesh") ) )
        .add( backend_options( prefixvm(prefix,"alemesh") ) )
        .add( modelbase_options( prefixvm(prefix,"alemesh") ) )
        .add( modelalgebraic_options( prefixvm(prefix,"alemesh.harmonic") ) )
        .add( modelalgebraic_options( prefixvm(prefix,"alemesh.winslow") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.winslow.l2proj") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.winslow.metric-derivative") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.winslow") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.ho") ) );
}

Feel::po::options_description
mixedpoisson_options( std::string const& prefix )
{
    po::options_description desc_options("mixedpoisson options");
    desc_options.add_options()
        ( prefixvm( prefix, "tau_constant" ).c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix, "use-sc" ).c_str(), po::value<bool>()->default_value( true ), "use static condensation" )
        ( prefixvm( prefix, "time-stepping").c_str(), Feel::po::value< std::string >()->default_value("BDF"), "time integration schema : BDF, Theta")
        ( prefixvm( prefix, "time-stepping.theta.value").c_str(), Feel::po::value< double >()->default_value(0.5), " Theta value")
        ;
    return desc_options.add( modelnumerical_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) ).add( backend_options( prefix+".sc" ) );
}

Feel::po::options_description
toolboxes_options( std::string const& type, std::string const& prefix )
{
    Feel::po::options_description toolboxesOptions("toolboxes options");

    if (type == "fluid")
        toolboxesOptions.add(fluidMechanics_options(prefix));
    else if (type == "solid")
        toolboxesOptions.add(solidMechanics_options(prefix));
    else if ( type == "heat" )
        toolboxesOptions.add( heat_options(prefix) );
    else if (type == "fsi")
        toolboxesOptions
            .add(fluidMechanics_options("fluid"))
            .add(solidMechanics_options("solid"))
            .add(fluidStructInteraction_options("fsi"));
    else if (type == "advection")
        toolboxesOptions.add(advection_options(prefix));
    else if (type == "levelset")
        toolboxesOptions.add(levelset_options(prefix));
    else if (type == "multifluid")
        toolboxesOptions.add(multifluid_options(prefix));
    else if (type == "electric")
        toolboxesOptions.add(electricity_options(prefix));
    else if (type == "thermo-electric")
        toolboxesOptions.add(thermoElectric_options(prefix));
    else if (type == "heat-fluid")
        toolboxesOptions.add(heatFluid_options(prefix));
    else if (type == "maxwell")
        toolboxesOptions.add(maxwell_options(prefix));
    else if (type == "coefficient-form-pdes")
        toolboxesOptions.add(coefficientformpdes_options(prefix));
    else if (type == "mixedpoisson")
        toolboxesOptions.add(mixedpoisson_options(prefix));
    else
        CHECK( false ) << "invalid type : " << type << " -> must be : fluid, solid, heat, fsi, advection, levelset, multifluid, thermo-electric, heat-fluid, heat-fluid, coefficient-form-pdes, mixedpoisson";

    return toolboxesOptions;
}


} // namespace Feel
