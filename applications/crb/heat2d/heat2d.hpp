/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2009-11-13

Copyright (C) 2009-2014 Feel++ Consortium

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
  \file heat2d.hpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2009-11-13
  */
#ifndef FEELPP_HEAT2D_HPP
#define FEELPP_HEAT2D_HPP 1

#include <boost/timer.hpp>

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>


namespace Feel
{
  /** \page CRB Certified Reduce Basis
    \author Vincent HUBER
    \date 2015-15-01
    \tableofcontents
    The goal of this project is to show how to use on a simple case the CRB FrameWork.
    \section PS Problem Setting
    Imagine a flat plate composed with two different materials.
    We want to evaluate the temperature in the whole plane for a various range of materials, that is a range of diffusivity parameters.
    \subsection SE System of Equation
    Let us name that two parameters \f$ k_s \f$ whith \f$ s \in \left\lbrace 1,2\right\rbrace\f$.

    \f{equation*}
    \begin{split}
    -k_s \Delta u &= 0 \\
    \left.u\right|_{\Gamma_b} &= \left.f\right|_{b} \text{ b}\in \left\lbrace 1,3\right\rbrace \\
    \left.\partial_n u\right|_{\Gamma_b} &= 0  \text{ b}\in \left\lbrace 2,4\right\rbrace \\
    \end{split}
    \f}
    \section Geo Geometry
    \image html heat2d_geo.png "Geometry of the model" width=2cm
    Our domain is decomposed into two surfaces (1 and 2).
    We impose a Dirichlet boundary condition on the left and right of the plane and an homogeneous Neumann one on the top and the bottom part.
    The whole geometry is the described by:
    \snippet heat2d.geo geo 
    \section AD Affine decomposition
    We need to indicate what is the affine decomposition to the RB FrameWork.
    We impose weak boundary conditions via \f \gamma \f that we will set at a big value.
    \f{equation*}
    \begin{split}
    k_1 \left(\int_{\Omega_1} \nabla u \cdot \nabla v  - \int_{\Gamma 1} \nabla u \cdot v \vec{n} + \nabla v \cdot u \vec{n} + \gamma u v \right) & \\
    + k_2 \left(\int_{\Omega_2} \nabla u \cdot \nabla v  - \int_{\Gamma 3} \nabla u \cdot v \vec{n} + \nabla v \cdot u \vec{n} + \gamma u v \right) & \\
    = k_1 f_1 \int_{\Gamma_1} \nabla v \cdot \vec{n} + \gamma v \\
    + k_2 f_2 \int_{\Gamma_2} \nabla v \cdot \vec{n} + \gamma v 
    \end{split}
    \f}
    \section IMP Implementation
    We have to define a model, that is a C++ class that provide two functions:
    - initmodel()
    - output()
    \subsection TM The Model
    We have two left hand side :
    \snippet heat2d.hpp lhs1
    That first left hand side depends on the first parameter - nammed mu0.
    \snippet heat2d.hpp lhs2
    That second left hand side depends on the second parameter - nammed mu1.
    and two right and side, depending on the first and second parameters :
    \snippet heat2d.hpp rhs1
    \snippet heat2d.hpp rhs2
    We have to indicate what is the energy matrix:
    \snippet heat2d.hpp energy
    and thus the whole code reads:
    \snippet heat2d.hpp initmodel
    \subsection RUN How does the code run ?
    In the CmakeLists.txt, one has to call the crb_add_model macro with:
    ~~~~~~~~~~~~
    crb_add_model(heat2d Heat2D HDRS heat2d.hpp
    LINK_LIBRARIES ${FEELPP_LIBRARIES}
    CFG heat2d.cfg 
    GEO heat2d.geo )
    ~~~~~~~~~~~~
    That will generate a code to encapsulate the model on the CRB FrameWork.
    Thus, the user has to use the `make -j 4 crb_heat2dapp` command to compile the code.
    \subsection OPT Options
    The CRB FrameWork offers a huge variety of options. We will present the more important :
    \snippet heat2d.cfg geo
    What is the used geometry.
    As long as the basis is constructed over that geometry, you must rebuild the database if you change it.
    \snippet heat2d.cfg fct
    The Boundary Conditions
    As long as the basis is constructed over theses BC, you must rebuild the database if you change it.
    \snippet heat2d.cfg run_mode
    \snippet heat2d.cfg model_opt
    \snippet heat2d.cfg crbopt 
    \snippet heat2d.cfg crbscmopt 
    Off course, one can use the command line to indicate theses options.
    ~~~~~~~~~~~~
    ./crb_heat1dapp --heat2d.run.sampling.size 22 --crb.compute-fem-during-online=false --heat1d.export-solution true
    ~~~~~~~~~~~~
  \subsection EXEC Execution
  The CRB method is divided into two parts:
  - offline,
  - online.

  That is actually one of the biggest interest of that method.
  The code you have compiled is able to deal with that two parts.
  For the first execution you have to create the database, and thus the online part is executed with a random set of parameters in the range we have defined in the code:
  \snippet heat2d.hpp parameters
  At first, the database is builded :
  \verbatim
  \endverbatim
  And after that the evaluations are done on the sampling :
  \verbatim
  \endverbatim
  Then, a summary is provided.
  \verbatim
  \endverbatim
  You will find the outputs in `$FEELPP_WORKDIR/feel/heat2d/Heat2D/np_1`. Open the file `heat2d.case` (or `heat2d-paraview-1.sos`) with paraview to visualize your results.
  */

  po::options_description
makeHeat2DOptions()
{
  po::options_description heat2doptions( "Heat2D options" );
  heat2doptions.add_options()
    ( "gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" )
    //( "k1", po::value<double>()->default_value( 0.1 ), "k1" )
    //( "k2", po::value<double>()->default_value( 0.1 ), "k2" )
    ;
  return heat2doptions;
}
AboutData
makeHeat2DAbout( std::string const& str = "heat2d" )
{
  Feel::AboutData about( /*AppName  */ str.c_str(),
      /*ProgName */ str.c_str(),
      /*Version  */ "0.1",
      /*ShortDesc*/ "2D Heat Benchmark",
      /*Licence  */ Feel::AboutData::License_GPL,
      /*Copyright*/ "Copyright (c) 2009-2014 Feel++ Consortium" );

  about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );
  return about;
}

/**
 * \class Heat2D
 * \brief brief description
 *
 * @author Vincent HUBER
 * @see
 */
class Heat2D : public ModelCrbBase<ParameterSpace<2>, decltype(Pch<3>(Mesh<Simplex<2>>::New()))>
{
  public:

    //! initialisation of the model
    void initModel();

    beta_vector_light_type beta;
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
};


/// [initmodel]
  void
Heat2D::initModel()
{

  CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
  LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
  LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
  /*
   * First we create the mesh
   */
  this->setFunctionSpaces( Pch<3>( loadMesh( _mesh=new Mesh<Simplex<2>> ) ) );

  /// [parameters]
  auto mu_min = Dmu->element();
  mu_min << 1, 1;
  Dmu->setMin( mu_min );
  auto mu_max = Dmu->element();
  mu_max << 20, 20;
  Dmu->setMax( mu_max );
  /// [parameters]

  auto u = Xh->element();
  auto v = Xh->element();
  auto mesh = Xh->mesh();
  //lhs
  /// [lhs1]
  auto a1 = form2( _trial=Xh, _test=Xh);
  a1 = integrate( markedelements( mesh,"SR" ),  gradt( u )*trans( grad( v ) )  )
    - integrate( markedfaces( mesh, "BR"), (gradt(u)*id(v)+grad(v)*idt(u))*N() + doption("gamma")*idt(u)*id(v) );
  this->addLhs( { a1 , "mu0" } );
  /// [lhs1]

  /// [lhs2]
  auto a2 = form2( _trial=Xh, _test=Xh);
  a2 = integrate( markedelements( mesh,"SL" ),  gradt( u )*trans( grad( v ) )  )
    - integrate( markedfaces( mesh, "BL"), (gradt(u)*id(v)+grad(v)*idt(u))*N() + doption("gamma")*idt(u)*id(v) );
  this->addLhs( { a2 , "mu1" } );
  /// [lhs2]

  //rhs
  /// [rhs1]
  auto f0 = form1( _test=Xh );
  f0 = integrate( markedfaces( mesh,"BR" ), -expr(soption("functions.f"))*(grad(v)*N()+doption("gamma")*id(v)) );
  this->addRhs( { f0, "mu0" } );
  /// [rhs1]

  /// [rhs2]
  auto f1 = form1( _test=Xh );
  f1 = integrate( markedfaces( mesh,"BL" ), -expr(soption("functions.g"))*(grad(v)*N()+doption("gamma")*id(v)) );
  this->addRhs( { f1, "mu1" } );
  /// [rhs2]

  /// [output]
  auto out1 = form1( _test=Xh );
  double meas = integrate(elements(mesh),cst(1.)).evaluate()(0,0);
  out1 = integrate( elements( mesh ), id( u )/cst(meas)) ;
  this->addOutput( { out1, "1" } );
  /// [output]

  // auto out2 = form1( _test=Xh );
  // out2 = integrate( elements( mesh ), idv( u )) / integrate( elements( mesh ), cst(1.)).evaluate()(0,0) ;
  // this->addOutput( { out2, "2" } );

  /// [energy]
  auto energy = form2( _trial=Xh, _test=Xh);
  energy = integrate( markedelements( mesh, "SL" ), gradt( u )*trans( grad( v ) )  )
    - integrate(    markedfaces( mesh, "BR" ), ( (id(v)*gradt(u)+idt(u)*grad(v))*N() + doption("gamma")*idt(u)*id(v) ))
    - integrate(    markedfaces( mesh, "BL" ), ( (id(v)*gradt(u)+idt(u)*grad(v))*N() + doption("gamma")*idt(u)*id(v) ))
    + integrate( markedelements( mesh, "SR" ), gradt( u )*trans( grad( v ) )  );
  this->addEnergyMatrix( energy );
  /// [energy]
}
/// [initmodel]

/// [output]
  double
Heat2D::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{

  //CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

  auto mesh = Xh->mesh();
  double output=0;
  // right hand side (compliant)
  if ( output_index == 0 )
  {
    output  = integrate( markedfaces( mesh,"BR" ), -mu(0)*expr(soption("functions.f"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0)
      + integrate( markedfaces( mesh,"BL" ), -mu(1)*expr(soption("functions.g"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0);
  }
  else if ( output_index == 1 )
  {
    output = mean(elements(mesh),idv(u))(0,0); 
  }
  // else if ( output_index == 2 )
  // {
  //     output = mean(elements(mesh),idv(u)).evaluat()(0,0); 
  // }
  else
    throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
  return output;

}
/// [output]

}

#endif /* FEELPP_HEAT2D_HPP */
