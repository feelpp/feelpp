/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
#include <feel/feelmodels2/modelcore/options.hpp>
#include <feel/feelmodels2/modelcore/feelmodelscoreconstconfig.h>
#include <feel/feelcore/feel.hpp>

namespace Feel
{

Feel::po::options_description envfeelmodels_options(std::string const& prefix)
{
    Feel::po::options_description envfsiOptions("Env FSI options");
    envfsiOptions.add_options()
        (prefixvm(prefix,"verbose").c_str(), Feel::po::value<bool>()->default_value( true ), "true or false to view verbose")
        (prefixvm(prefix,"evalpt-path").c_str(), Feel::po::value< std::string >(), "path of file containing set of point to evaluate")
        (prefixvm(prefix,"comm.useNonStandartCompositeSpace").c_str(), Feel::po::value< bool >()->default_value( false ), "try old build of composite space")
        (prefixvm(prefix,"rebuildOnlyAMeshPartion").c_str(), Feel::po::value< bool >()->default_value( false ), "only rebuild a mesh partition")
        ;
    envfsiOptions.add( Feel::gmsh_options("meshpartitioning") );
    envfsiOptions.add_options()
        (prefixvm(prefix,"meshpartitioning.gmsh.outputfilename").c_str(), Feel::po::value< std::string >(), "filename of partitioned mesh")
        (prefixvm(prefix,"meshpartitioning.gmsh.partitions").c_str(), Feel::po::value< int >(), "number of partition")
        ;

    return envfsiOptions;
}

Feel::po::options_description applibase_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"verbose").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose")
        (prefixvm(prefix,"verbose_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        //(prefixvm(prefix,"verbose_solvertimer").c_str(), Feel::po::value<bool>()->default_value( false/*true*/ ), "true or false to view verbose")
        //(prefixvm(prefix,"verbose_solvertimer_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
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
        ;
    return appliBaseOptions;
}

Feel::po::options_description applibaseMethodsNum_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base Methods Num options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"verbose_solvertimer").c_str(), Feel::po::value<bool>()->default_value( false/*true*/ ), "true or false to view verbose")
        (prefixvm(prefix,"verbose_solvertimer_allproc").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to view verbose for all proc")
        (prefixvm(prefix,"linearsystem-cst-update").c_str(), Feel::po::value< bool >()->default_value( true ), "update matrix and rhs cst part")
        (prefixvm(prefix,"jacobian-linear-update").c_str(), Feel::po::value< bool >()->default_value( true ), "update linear part of the jacobian")
        (prefixvm(prefix,"residual-uselinearjac").c_str(), Feel::po::value< bool >()->default_value( true ), "update linear part of the residual with linear jacobian")
        (prefixvm(prefix,"preconditioner.contribution").c_str(), Feel::po::value< std::string >()->default_value("same_matrix"),
         "contribution in preconditioner : same_matrix, standart, extended ")
        (prefixvm(prefix,"use-cst-matrix").c_str(), Feel::po::value< bool >()->default_value( true ), "use-cst-matrix")
        (prefixvm(prefix,"use-cst-vector").c_str(), Feel::po::value< bool >()->default_value( true ), "use-cst-vector")
        (prefixvm(prefix,"error-if-solver-not-converged").c_str(), Feel::po::value< bool >()->default_value( false ), "error-if-solver-not-converged")
        (prefixvm(prefix,"clear-preconditioner-after-use").c_str(), Feel::po::value< bool >()->default_value( false ), "clear-preconditioner-after-use")
        (prefixvm(prefix,"graph-print-python").c_str(), Feel::po::value<bool>()->default_value( false ), "print graph in python script")
        (prefixvm(prefix,"graph-print-python-filename").c_str(), Feel::po::value< std::string >(), "filename python graph")
        ;
    return appliBaseOptions.add( applibase_options(prefix ) );//.add( backend_options( prefix ) );
}


Feel::po::options_description applibaseNumericalSimulationTransitory_options(std::string const& prefix)
{
    Feel::po::options_description appliBaseOptions("Application Base options");
    appliBaseOptions.add_options()
        (prefixvm(prefix,"filename").c_str(), Feel::po::value<std::string>()->default_value( "" ), "json file describing model properties" )
        // mesh
        (prefixvm(prefix,"geofile").c_str(), Feel::po::value< std::string >(), "input geo file")
        (prefixvm(prefix,"mshfile").c_str(), Feel::po::value< std::string >(), "input msh file")
        (prefixvm(prefix,"geotool-mesh-index").c_str(), Feel::po::value< int >()->default_value( 0 ), "choice mesh")
        (prefixvm(prefix,"geotool-save-directory").c_str(), Feel::po::value< std::string >()->default_value(""), "path for save geo and msh file") // work only with mesh in .mesh file
        (prefixvm(prefix,"geotool-save-name").c_str(), Feel::po::value< std::string >()->default_value(""), "name for save geo and msh file") // work only with mesh in .mesh file
        // other
        (prefixvm(prefix,"rebuild_mesh_partitions").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to rebuild mesh partitions ")
        (prefixvm(prefix,"geomap").c_str(), Feel::po::value< std::string >()->default_value("opt"), "geomap strategy : ho, opt ")
        ;

    for (uint16_type nParam=1;nParam<=FEELMODELS_FSIBASE_NUMBER_OF_PARAMETERS;++nParam)
        appliBaseOptions.add_options()
            ((boost::format("%1%parameter%2%") %prefixvm(prefix,"") %nParam).str().c_str(), Feel::po::value<double>()->default_value( 1.0 ), " a parameter");

    for (uint16_type nParam=1;nParam<=FEELMODELS_FSIBASE_NUMBER_OF_GEOPARAMETERS;++nParam)
        appliBaseOptions.add_options()
            ((boost::format("%1%geo-parameter%2%") %prefixvm(prefix,"") %nParam).str().c_str(), Feel::po::value<double>()->default_value( 1.0 ), " a geo parameter");

    for (uint16_type k=1;k<=FEELMODELS_FSIBASE_NUMBER_OF_GINACEXPR;++k)
        {
            appliBaseOptions.add_options()
                ((boost::format("%1%ginac-expr%2%") %prefixvm(prefix,"") %k).str().c_str(), Feel::po::value<std::string>()->default_value( "" ), " a ginac expr");
            appliBaseOptions.add_options()
                ((boost::format("%1%ginac-name%2%") %prefixvm(prefix,"") %k).str().c_str(),
                 Feel::po::value<std::string>()->default_value( (boost::format("defaultNameGinacExpr%1%")%k).str() ), " name of ginac expr");
        }
    appliBaseOptions.add_options()
        (prefixvm(prefix,"ginac-expr-directory").c_str(), Feel::po::value< std::string >(), "ginac-expr-directory");

    return appliBaseOptions
        .add( applibaseMethodsNum_options( prefix ))
        .add( backend_options( prefix ) );
    //return appliBaseOptions.add( backend_options( prefix ) );
}

/**
 * generate options for the fluid solver
 */
Feel::po::options_description
fluidMechanics_options(std::string const& prefix)
{
    Feel::po::options_description fluidOptions("Fluid Mechanics options");
    fluidOptions.add_options()
        (prefixvm(prefix,"hsize").c_str(), Feel::po::value<double>()->default_value( 0.06 ), "h fluid")
        (prefixvm(prefix,"model").c_str(), Feel::po::value< std::string >()->default_value("Navier-Stokes"), "fluid model : Navier-Stokes,Oseen,Stokes")
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >(), "fluid solver")
        (prefixvm(prefix,"stress_tensor_law").c_str(), Feel::po::value< std::string >()->default_value("newtonian"), "newtonian, power_law, walburn-schneck_law, carreau_law, carreau-yasuda_law ")
        ( prefixvm(prefix,"start-by-solve-newtonian").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-newtonian")
        ( prefixvm(prefix,"start-by-solve-stokes-stationary").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-stokes-stationary")
        ( prefixvm(prefix,"start-by-solve-stokes-stationary.do-export").c_str(), Feel::po::value<bool>()->default_value( false ), "start-by-solve-stokes-stationary.do-export")
        ( prefixvm(prefix,"start-by-solve-stokes-stationary.time-value-used-in-bc").c_str(), Feel::po::value<double>()->default_value( 0. ), "time-value-used-in-bc")
        //(prefixvm(prefix,"strain_tensor.use-sym-tensor").c_str(), Feel::po::value< bool >()->default_value(true), "sym tensor or not ")
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1050 ), "density [ kg.m^3]")
        (prefixvm(prefix,"mu").c_str(), Feel::po::value<double>()->default_value( 0.00345 ), "dynamic viscosity [ Pa.s = kg/(m.s^2) ]")
        (prefixvm(prefix,"viscosity.zero_shear").c_str(), Feel::po::value< double >()->default_value( 0.056 ), "parameter mu_0 for generalized Newtonian [ Pa.s ]  ")
        (prefixvm(prefix,"viscosity.infinite_shear").c_str(), Feel::po::value< double >()->default_value( 0.00345 ), "parameter mu_inf for generalized Newtonian [ Pa.s ] ")
        (prefixvm(prefix,"hematocrit").c_str(), Feel::po::value< double >()->default_value( 40 ), "hematocrit : RBC volume fraction (Vrbc/Vtotal) [ percentage ] ")
        (prefixvm(prefix,"power_law.n").c_str(), Feel::po::value< double >()->default_value( 0.6 ), "parameter n in power_law ")
        (prefixvm(prefix,"power_law.k").c_str(), Feel::po::value< double >()->default_value( 0.035 ), "parameter k in power_law [ Pa.s^n ] ")
        (prefixvm(prefix,"carreau_law.lambda").c_str(), Feel::po::value< double >()->default_value( 3.313 ), "parameter lambda in carreau_law [ s ] ")
        (prefixvm(prefix,"carreau_law.n").c_str(), Feel::po::value< double >()->default_value( 0.3568 ), "parameter n in carreau_law ")
        (prefixvm(prefix,"carreau-yasuda_law.lambda").c_str(), Feel::po::value< double >()->default_value( 1.902 ), "parameter lambda in carreau-yasuda_law [ s ] ")
        (prefixvm(prefix,"carreau-yasuda_law.n").c_str(), Feel::po::value< double >()->default_value( 0.22 ), "parameter n in carreau-yasuda_law ")
        (prefixvm(prefix,"carreau-yasuda_law.a").c_str(), Feel::po::value< double >()->default_value( 1.25 ), "parameter a in carreau-yasuda_law ")
        (prefixvm(prefix,"walburn-schneck_law.C1").c_str(), Feel::po::value< double >()->default_value( 0.00797 ), "parameter C1 in walburn-schneck_law ")
        (prefixvm(prefix,"walburn-schneck_law.C2").c_str(), Feel::po::value< double >()->default_value( 0.0608 ), "parameter C2 in walburn-schneck_law ")
        (prefixvm(prefix,"walburn-schneck_law.C3").c_str(), Feel::po::value< double >()->default_value( 0.00499 ), "parameter C3 in walburn-schneck_law ")
        (prefixvm(prefix,"walburn-schneck_law.C4").c_str(), Feel::po::value< double >()->default_value( 14.585 ), "parameter C4 in walburn-schneck_law [l/g] ")
        (prefixvm(prefix,"TPMA").c_str(), Feel::po::value< double >()->default_value( 25.9 ), "parameter TPMA (Total Proteins Minus Albumin) [ g/l ] ")
        (prefixvm(prefix,"stabilisation-pspg").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-gls").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-convection-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-pressure").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-pressure-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-divergence").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-cip-divergence-gamma").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-divergence").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-div-div").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"stabilisation-div-div-beta").c_str(), Feel::po::value< double >()->default_value( 1.0 ), "parameter beta in div-div stab ")
        (prefixvm(prefix,"define-pressure-cst").c_str(), Feel::po::value<bool>()->default_value( false ), "maybe need to define the pressure constant if only dirichlet bc")
        (prefixvm(prefix,"define-pressure-cst.method").c_str(), Feel::po::value<std::string>()->default_value( "lagrange-multiplier"), " lagrange-multiplier or penalisation")
        (prefixvm(prefix,"define-pressure-cst.penalisation-beta").c_str(), Feel::po::value< double >()->default_value( -10e-10 ), "parameter beta in cstpressure penalisation ")
        //(prefixvm(prefix,"stabilisation-cstpressure").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        //(prefixvm(prefix,"stabilisation-cstpressure-beta").c_str(), Feel::po::value< double >()->default_value( -10e-10 ), "parameter beta in cstpressure stab ")
        (prefixvm(prefix,"stabilisation-convection-energy").c_str(), Feel::po::value<bool>()->default_value( false ), "use stabilisation method")
        (prefixvm(prefix,"dirichletbc.type").c_str(), Feel::po::value<std::string>()->default_value( "elimination" ), "elimination, nitsche, lm")
        (prefixvm(prefix,"dirichletbc.lm.savemesh").c_str(), Feel::po::value<bool>()->default_value( false ), "export Lagrange multipliers mesh")
        (prefixvm(prefix,"dirichletbc.lm.use-submesh-relation").c_str(), Feel::po::value<bool>()->default_value( true ), "use submesh relation")
        (prefixvm(prefix,"dirichletbc.nitsche.gamma").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        /*ASUP*/(prefixvm(prefix,"useweakbc").c_str(), Feel::po::value<bool>()->default_value( false ), "use weak BC?")
        /*ASUP*/(prefixvm(prefix,"weakbccoeff").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"bc-slip-form").c_str(), Feel::po::value<int>()->default_value( 1 ), "formulation for slip condition (1 or 2)")
        (prefixvm(prefix,"bc-slip-gammaN").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"bc-slip-gammaTau").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"hovisu").c_str(), Feel::po::value<bool>(), "true or false for high order visualisation")
        (prefixvm(prefix,"hovisu.space-used").c_str(), Feel::po::value<std::string>()->default_value( "velocity" ), "velocity , pressure, p1 ")
        (prefixvm(prefix,"do_export_meshale").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportMeshALE")
        (prefixvm(prefix,"do_export_vorticity").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportVorticity")
        (prefixvm(prefix,"do_export_normalstress").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportNormalStress")
        (prefixvm(prefix,"do_export_wallshearstress").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportWallShearStress")
        (prefixvm(prefix,"do_export_meshdisplacementoninterface").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportMeshDisplacementOnInterface")
        (prefixvm(prefix,"do_export_viscosity").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportViscosity")
        (prefixvm(prefix,"do_export_all").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportAll")
        (prefixvm(prefix,"fluid-outlet.number").c_str(), Feel::po::value<int>()->default_value( 1 ), "number of fluid outlet")
        (prefixvm(prefix,"fluid-outlet.type").c_str(), Feel::po::value<std::string>()->default_value( "free" ), "type : free, windkessel ")
        (prefixvm(prefix,"fluid-outlet.windkessel.coupling").c_str(), Feel::po::value<std::string>()->default_value( "implicit" ), "explicit, implicit ")

        (prefixvm(prefix,"periodicity.translate-x").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-x")
        (prefixvm(prefix,"periodicity.translate-y").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-y")
        (prefixvm(prefix,"periodicity.translate-z").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "periodicity.translate-z")
        (prefixvm(prefix,"periodicity.marker1").c_str(), Feel::po::value<std::string>(), "periodicity.marker1 ")
        (prefixvm(prefix,"periodicity.marker2").c_str(), Feel::po::value<std::string>(), "periodicity.marker2 ")
        (prefixvm(prefix,"periodicity.pressure-jump").c_str(), Feel::po::value<double>()->default_value(1.0), "periodicity.pressure-jump ")
        ;

        for (uint16_type nBFO=0;nBFO<=5;++nBFO)
            fluidOptions.add_options()
                ((boost::format("%1%fluid-outlet.windkessel.Rd%2%") %prefixvm(prefix,"") %nBFO).str().c_str(), Feel::po::value<double>()->default_value( 1.0 ), " a parameter")
                ((boost::format("%1%fluid-outlet.windkessel.Rp%2%") %prefixvm(prefix,"") %nBFO).str().c_str(), Feel::po::value<double>()->default_value( 1.0 ), " a parameter")
                ((boost::format("%1%fluid-outlet.windkessel.Cd%2%") %prefixvm(prefix,"") %nBFO).str().c_str(), Feel::po::value<double>()->default_value( 1.0 ), " a parameter")
                ;

    return fluidOptions.add( applibaseNumericalSimulationTransitory_options( prefix ) ).add( bdf_options( prefix ) ).add( alemesh_options( prefix ) );
}

Feel::po::options_description
solidMechanics_options(std::string const& prefix)
{
    Feel::po::options_description solidOptions("Solid Mechanics options");
    solidOptions.add_options()
        (prefixvm(prefix,"hsize").c_str(), Feel::po::value<double>()->default_value( 0.02 ), "h struct")
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "density")
        (prefixvm(prefix,"youngmodulus").c_str(), Feel::po::value<double>()->default_value( 3.e6 ), "young modulus")
        (prefixvm(prefix,"coeffpoisson").c_str(), Feel::po::value<double>()->default_value( 0.3 ), "poisson coefficient")
        (prefixvm(prefix,"model").c_str(), Feel::po::value< std::string >()->default_value("Elasticity"), "struct model")
        (prefixvm(prefix,"material_law").c_str(), Feel::po::value< std::string >()->default_value("StVenantKirchhoff"), "StVenantKirchhoff, NeoHookean")
        (prefixvm(prefix,"mechanicalproperties.compressible.volumic-law").c_str(), Feel::po::value< std::string >()->default_value("classic"), "classic, simo1985")
        (prefixvm(prefix,"mechanicalproperties.compressible.neohookean.variant").c_str(),
         Feel::po::value< std::string >()->default_value("default"), "default, molecular-theory, molecular-theory-simo1985")
        (prefixvm(prefix,"use-incompressibility-constraint").c_str(), Feel::po::value<bool>()->default_value( false ), "use incompressibility constraint?")
        (prefixvm(prefix,"solver").c_str(), Feel::po::value< std::string >(), "struct solver")
        (prefixvm(prefix,"time-schema").c_str(), Feel::po::value< std::string >()->default_value("Newmark"), "time integration schema : Newmark, Generalized-Alpha")
        (prefixvm(prefix,"time-rho").c_str(), Feel::po::value< double >()->default_value(0.8), " Generalized-Alpha parameter")
        //(prefixvm(prefix,"useweakbc").c_str(), Feel::po::value<bool>()->default_value( false ), "use weak BC?")
        //(prefixvm(prefix,"weakbccoeff").c_str(), Feel::po::value<double>()->default_value( 10.0 ), "coeff for weak Dirichlet conditions")
        (prefixvm(prefix,"1dreduced-geofile").c_str(), Feel::po::value< std::string >(), "input geo file 1dreduced")
        (prefixvm(prefix,"hovisu").c_str(), Feel::po::value<bool>(), "true or false for high order visualisation")
        (prefixvm(prefix,"hovisu.space-used").c_str(), Feel::po::value<std::string>()->default_value( "displacement" ), "displacement, pressure, p1 ")
        (prefixvm(prefix,"do_export_velocity").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportVelocity")
        (prefixvm(prefix,"do_export_acceleration").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportAcceleration")
        (prefixvm(prefix,"do_export_normalstress").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportNormalStress")
        (prefixvm(prefix,"do_export_velocityinterfacefromfluid").c_str(), Feel::po::value<bool>()->default_value( false ), "doExportVelocityInterfaceFromFluid")
        (prefixvm(prefix,"predictor-disp-type").c_str(), Feel::po::value<std::string>()->default_value( "extrap-with-vel" ), "predictor-disp-type (bdf,extrap-with-vel,identity)")
        (prefixvm(prefix,"use-null-space").c_str(), Feel::po::value<bool>()->default_value( false ), "use-null-space")
        (prefixvm(prefix,"use-near-null-space").c_str(), Feel::po::value<bool>()->default_value( true ), "use-near-null-space")
        ;

    return solidOptions.add( applibaseNumericalSimulationTransitory_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
}

Feel::po::options_description
fluidStructInteraction_options( std::string const& prefix )
{
    Feel::po::options_description FSIoptions("FSI options");
    FSIoptions.add_options()
        (prefixvm(prefix,"solid-mesh.extract-1d-from-fluid-mesh").c_str(), Feel::po::value<bool>()->default_value( false ), "solid-mesh.extract-1d-from-fluid-mesh")
        //(prefixvm(prefix,"solid-mesh.extract-1d-from-fluid-mesh.boundary-markers").c_str(), po::value<std::vector<std::string> >()->multitoken(), "boundary-markers")

        (prefixvm(prefix,"coupling-type").c_str(),Feel::po::value< std::string >()->default_value("Implicit"), " Implicit or Semi-Implicit")
        (prefixvm(prefix,"coupling-bc").c_str(),Feel::po::value< std::string >()->default_value("dirichlet-neumann"), " dirichlet-neumann, robin-robin,robin-neumann")
        (prefixvm(prefix,"fixpoint.tol").c_str(), Feel::po::value<double>()->default_value( 1.e-6 ), "tolerance pt fixe")
        (prefixvm(prefix,"fixpoint.initialtheta").c_str(), Feel::po::value<double>()->default_value( 1. ), "relax aitken parameter")
        (prefixvm(prefix,"fixpoint.min_theta").c_str(), Feel::po::value<double>()->default_value( 1.e-4 ), "min theta parameter")
        (prefixvm(prefix,"fixpoint.maxit").c_str(), Feel::po::value<int>()->default_value( 1000 ), "max iteration")
        (prefixvm(prefix,"conforming-interface").c_str(), Feel::po::value<bool>()->default_value( false ), " fsi interface is conforme?")
        (prefixvm(prefix,"fluid.reuse-prec.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi fluid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"solid.reuse-prec.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi solid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"fluid.reuse-jac.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi fluid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"solid.reuse-jac.rebuild-at-first-fsi-step").c_str(), Feel::po::value<bool>()->default_value( true ), " fsi solid reuse-prec.rebuild-at-first-fsi-step")
        (prefixvm(prefix,"coupling-robin-robin.gamma").c_str(), Feel::po::value<double>()->default_value( 2500 ), "nitsche parameters")
        (prefixvm(prefix,"coupling-robin-robin.gamma0").c_str(), Feel::po::value<double>()->default_value( 1 ), "nitsche parameters")
        ;
    return FSIoptions.add( applibaseMethodsNum_options(prefix) );/*FSIoptions.add(alemesh_options());*/
}


Feel::po::options_description
thermoDynamics_options(std::string const& prefix)
{
    Feel::po::options_description thermoDynamicsOptions("Thermo Dynamics options");
    thermoDynamicsOptions.add_options()
        (prefixvm(prefix,"hsize").c_str(), Feel::po::value<double>()->default_value( 0.02 ), "h struct")
        (prefixvm(prefix,"thermal-conductivity").c_str(), Feel::po::value<double>()->default_value( 1 ), "thermal-conductivity [ W/(m*K) ]")
        (prefixvm(prefix,"rho").c_str(), Feel::po::value<double>()->default_value( 1 ), "density [ kg/(m^3) ]")
        (prefixvm(prefix,"heat-capacity").c_str(), Feel::po::value<double>()->default_value( 1 ), "heat-capacity [ J/(kg*K) ]")
        (prefixvm(prefix,"use_velocity-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "use-velocity-convection")
        (prefixvm(prefix,"velocity-convection_is_incompressible").c_str(), Feel::po::value<bool>()->default_value( false ), "velocity-convection-is-incompressible")
        (prefixvm(prefix,"do_export_all").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_all")
        (prefixvm(prefix,"do_export_velocity-convection").c_str(), Feel::po::value<bool>()->default_value( false ), "do_export_velocity-convection")
        ;
    return thermoDynamicsOptions.add( applibaseNumericalSimulationTransitory_options( prefix ) ).add( bdf_options( prefix ) ).add( ts_options( prefix ) );
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
        (prefixvm(prefix,"alemesh.winslow.solver").c_str(), Feel::po::value<std::string>()->default_value( "newton" ), "solver : linear, ptfixe, newton")
        (prefixvm(prefix,"alemesh.winslow.ptfixe.niter").c_str(), Feel::po::value<int>()->default_value( 10 ), "nb max iteration of point fixe")
        (prefixvm(prefix,"alemesh.winslow.ptfixe.tol").c_str(), Feel::po::value<double>()->default_value( 1e-5 ), "tol point fixe")
        (prefixvm(prefix,"alemesh.winslow.aitken.initialtheta").c_str(), Feel::po::value<double>()->default_value( 0.3 ), "initial theta of aitken relax")
        (prefixvm(prefix,"alemesh.harmonic.use_adaptive_penalisation").c_str(), Feel::po::value<bool>()->default_value( true ), "use tau penalisation in low formulation")
        (prefixvm(prefix,"alemesh.apply-ho-correction").c_str(), Feel::po::value<bool>()->default_value( true ), "do apply ho correction")
        //(prefixvm(prefix,"alemesh.export").c_str(), Feel::po::value<bool>()->default_value( false ), "true or false to export")
        ;
    return desc_options
        .add( bdf_options( prefixvm(prefix,"alemesh") ) )
        .add( backend_options( prefixvm(prefix,"alemesh") ) )
        .add( applibase_options( prefixvm(prefix,"alemesh") ) )
        .add( applibaseMethodsNum_options( prefixvm(prefix,"alemesh.harmonic") ) )
        .add( applibaseMethodsNum_options( prefixvm(prefix,"alemesh.winslow") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.winslow.l2proj") ) )
        .add( backend_options( prefixvm(prefix,"alemesh.ho") ) );
}



Feel::po::options_description
feelmodels_options(std::string type)
{
    Feel::po::options_description FSIoptions("FSI options");

    FSIoptions.add( envfeelmodels_options("master") );

    if (type == "fluid")
        FSIoptions.add(fluidMechanics_options("fluid"));
    else if (type == "solid")
        FSIoptions.add(solidMechanics_options("struct"));
    else if ( type == "thermo-dynamics" )
        FSIoptions.add( thermoDynamics_options("thermo") );
    else if (type == "fsi")
        FSIoptions
            .add(fluidMechanics_options("fluid"))
            .add(solidMechanics_options("struct"))
            .add(fluidStructInteraction_options("fsi"));

    return FSIoptions
        .add( Feel::feel_options() );
}


} // namespace Feel
