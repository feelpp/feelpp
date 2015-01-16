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
 We consider a flat plate composed with two different materials.
 We want to evaluate the temperature in the whole plane for a various range of materials, that is a range of diffusivity parameters.

 \subsection SF Strong Formulation
 Let us name that two parameters \f$ \mu_1 \f$  and \f$\mu_2\f$,
 we then define \f$\mu=(\mu_1,\mu_2)\in\mathcal D\subset\mathbb R\f$.
 The governing equations of our problems are
 \f{equation*}
 \left\{
 \begin{split}
 -\mu_s \Delta u &= 0, \text{ in } \Omega_s \\
 \left.u\right|_{\Gamma_b} &= \left.f\right|_{b}, \text{ on } \Gamma_b, \text{ b}\in \left\lbrace 1,3\right\rbrace \\
 \left.\partial_n u\right|_{\Gamma_b} &= 0, \text{ on } \Gamma_b, \text{ b}\in \left\lbrace 2,4\right\rbrace \\
 \end{split}\right.
 \f}

 \subsection geometry Geometry
 \image html heat2d_geo.png  "Geometry of the model"
 Our domain (\f$\Omega\f$) is decomposed into two surfaces ( \f$\Omega_1\f$ and \f$\Omega_2\f$ ).
 We impose a Dirichlet boundary condition on the left (\f$\Gamma_1\f$) and right (\f$\Gamma_3\f$) of the plane and an homogeneous Neumann one on the top and the bottom part (\f$\Gamma_2\cup\Gamma_4\f$).
 The whole geometry is described in \c "applications/crb/heat2d/heat2d.geo"

 \section Method Certified Reduced Basis Method

 \subsection WF Weak Formulation

 An usual integration by parts of the strong formulation leads to weak the formulation of our problem
 \f{equation*}
 a(u,v;\mu)=l(v;\mu),\ \forall v \in H^1(\Omega)
 \f}
 with
\f{equation*}
 \begin{split}
 a(u,v;\mu)&=k_1 \left(\int_{\Omega_1} \nabla u \cdot \nabla v\ d\Omega  - \int_{\Gamma 1} \nabla u \cdot v \vec{n} + \nabla v \cdot u \vec{n} + \gamma u v \ d\Omega\right)
 + k_2 \left(\int_{\Omega_2} \nabla u \cdot \nabla v\ d\Gamma  - \int_{\Gamma 3} \nabla u \cdot v \vec{n} + \nabla v \cdot u \vec{n} + \gamma u v \ d\Omega\right) \\
 l(v;\mu)&=-k_1 \int_{\Gamma_1} f_1(\nabla v \cdot \vec{n} + \gamma v - k_2)\ d\Gamma -k_2\int_{\Gamma_3} f_2(\nabla v \cdot \vec{n} + \gamma v)\ d\Gamma
 \end{split}
\f}
 Note : We chose to weakly impose the Dirichlet condition using a penalty (\f$\gamma\f$) method.

 \subsection AD Affine Decomposition
 We can now write the affine decompostion, which is an essential component of the Reduce Basis Method :
 \f{equation*}
 \begin{split}
 a(u,v;\mu)&=\mu_0a_0(u,v)+\mu_1a_1(u,v),\\
 l(v;\mu)&=\mu_0f_0(v)+\mu_1f_1(v)
 \end{split}
 \f}

 \section IMP Implementation
 We have to define a model, that is a C++ class that provide two functions:
 - initmodel()
 - output()
 \subsection TM The Model
 We have two terms in the affine decomposition of the left hand side :

 \f$ \mu_0a_0(u,v)\f$
 \snippet heat2d.hpp lhs1
 and \f$ \mu_1a_1(u,v)\f$
 \snippet heat2d.hpp lhs2

 We also have two terms in the affine decompostion of the right hand side :

 \f$ \mu_0f_0(v) \f$
 \snippet heat2d.hpp rhs1
 and \f$ \mu_1f_1(v)\f$
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
 [CRB::offline] strategy 1
 -- primal problem solved in 0.0928769 s
 -- dual problem solved in 0.095052 s
 -- complement of M_WNmu built in 0.000298023 s
 -- time to add the primal basis : 2.19345e-05 s
 -- time to add the dual basis : 9.05991e-06 s
 -- time to add primal and dual basis : 3.00407e-05 s
 -- orthonormalization (Gram-Schmidt)
 -- orthonormalization (Gram-Schmidt)
 -- primal orthonormalization : 0.00181198 s
 -- dual orthonormalization : 0.001683 s
 -- projection on reduced basis space : 0.0170138 s
 -- offlineResidual update starts
 o N=1 QLhs=2 QRhs=2 Qoutput=2
 o initialize offlineResidual in 0.060207s
 o M_C0_pr updated in 0.154062s
 o Lambda_pr updated in 0.085693s
 o Gamma_pr updated in 0.119294s
 o C0_du updated in 0.05943s
 o Lambda_du updated in 0.059346s
 o Gamma_du updated in 0.117562s
 -- offlineResidual updated in 0.676627s
 [CRB maxerror] proc 0 delta_pr : 1.26973456445329 -- delta_du : 1.26973455446957 -- output error : 1.69949977698333
 -- max error bounds computed in 1.47827506065369s
 ============================================================
 \endverbatim
 And after that the evaluations are done on the sampling :
 \verbatim
 CRB mode -- 1/10
 output=38336.4609373806 with 2 basis functions  (error estimation on this output : 4.94095341847669e-09)
 \endverbatim
 Then, a summary is provided.
 \verbatim
 mu0        mu1     FEM Output       FEM Time      RB Output   Error Bounds       CRB Time   output error  Conditionning       l2_error       h1_error
 9.8395e+00 3.3342e+00     3.8336e+04     5.6988e-02     3.8336e+04     1.2888e-13     6.2704e-04     1.2716e-14     0.0000e+00     1.3640e-10     2.3204e-09
 3.7843e+00 1.7388e+00     1.5503e+04     5.1528e-02     1.5503e+04     1.2259e-13     7.1502e-04     4.3765e-14     0.0000e+00     8.8669e-11     1.4988e-09
 1.4905e+00 1.8228e+01     3.5259e+04     4.9509e-02     3.5259e+04     2.7409e-13     7.4887e-04     2.0636e-16     0.0000e+00     1.2138e-11     1.8794e-10
 1.1129e+00 3.7174e+00     9.8869e+03     5.0749e-02     9.8869e+03     1.8923e-13     7.2312e-04     1.8398e-15     0.0000e+00     5.7805e-11     9.1918e-10
 6.6080e+00 2.5546e+00     2.6270e+04     5.1569e-02     2.6270e+04     1.2758e-13     7.1120e-04     7.9628e-14     0.0000e+00     1.1580e-10     1.9645e-09
 1.1150e+01 1.9748e+01     6.9981e+04     5.0854e-02     6.9981e+04     1.4065e-13     7.6890e-04     1.2476e-15     0.0000e+00     5.5962e-11     9.0763e-10
 1.4922e+00 5.5829e+00     1.4251e+04     5.0819e-02     1.4251e+04     1.9755e-13     7.1406e-04     2.2975e-15     0.0000e+00     5.5007e-11     8.7190e-10
 1.0216e+00 1.7240e+01     3.2054e+04     4.7903e-02     3.2054e+04     2.9595e-13     7.3981e-04     1.8159e-15     0.0000e+00     2.3779e-12     3.6511e-11
 2.9945e+00 5.8006e+00     1.9620e+04     5.0702e-02     1.9620e+04     1.4777e-13     7.2622e-04     4.6356e-15     0.0000e+00     5.8542e-11     9.4676e-10
 7.4591e+00 1.8562e+01     5.5708e+04     5.0537e-02     5.5708e+04     1.6726e-13     7.3314e-04     5.6162e-15     0.0000e+00     6.1326e-11     9.8386e-10
 \endverbatim
 You will find the outputs in `$FEELPP_WORKDIR/feel/heat2d/Heat2D/np_1`. Open the file `heat2d.case` (or `heat2d-paraview-1.sos`) with paraview to visualize your results.
 You will find - depending on the boundary conditions and mesh you have used something like:
 \image html heat2d_res.png "Result of the CRB evaluation
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
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( markedelements( mesh,"SR" ),  gradt( u )*trans( grad( v ) )  )
        - integrate( markedfaces( mesh, "BR"), (gradt(u)*id(v)+grad(v)*idt(u))*N() + doption("gamma")*idt(u)*id(v) );
    this->addLhs( { a0 , "mu0" } );
    /// [lhs1]

    /// [lhs2]
    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( markedelements( mesh,"SL" ),  gradt( u )*trans( grad( v ) )  )
        - integrate( markedfaces( mesh, "BL"), (gradt(u)*id(v)+grad(v)*idt(u))*N() + doption("gamma")*idt(u)*id(v) );
    this->addLhs( { a1 , "mu1" } );
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
