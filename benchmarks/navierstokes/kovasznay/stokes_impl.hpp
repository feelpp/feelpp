/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-06-18

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file stokes_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-06-18
 */

namespace Feel
{
template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
Stokes<Dim, _OrderU, _OrderP, Entity>::Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    M_stabP( this->vm()["penalisation"].template as<double>() ),
    M_stabD( this->vm()["penalisation"].template as<double>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{
    M_weak_dirichlet = this->vm().count( "weak" );
    mu = this->vm()["mu"].template as<value_type>();
    Parameter h;
    switch( OrderU )
    {
    case 1:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.05" );
        break;
    case 2:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.1" );
        break;

    case 3:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.035:0.025:0.2" );
        break;
    }
    this->
        //addParameter( Parameter(_name="mu",_type=CONT_ATTR,_latex="\\mu", _values=boost::lexical_cast<std::string>( mu ).c_str()))
        addParameter( Parameter(_name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str()) )
        .addParameter( Parameter(_name="orderU",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderU  ).c_str()) )
        .addParameter( Parameter(_name="orderP",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderP  ).c_str()) )
        .addParameter( h );

    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    depend.push_back(h);
    std::ostringstream oss;
    if ( OrderP == OrderU )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP  );
    else if ( OrderP == OrderU-1 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderU  );
    else if ( OrderP == OrderU-2 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP+1  );
    funcs.push_back(oss.str());
    oss.str("");
    std::vector<std::string> funcs2;
    oss << "h**" << boost::lexical_cast<std::string>( OrderP+1 ) ;
    funcs2.push_back(oss.str());
    oss.str("");
    std::vector<std::string> funcs3;
    if ( OrderP == OrderU-2 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderU ) ;
    else
        oss << "h**" << boost::lexical_cast<std::string>( OrderU+1 ) ;
    funcs3.push_back(oss.str());

    this->
        addOutput( Output(_name="norm_L2_u",_latex="\\left\\| u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs3) )
        .addOutput( Output(_name="norm_H1_u",_latex="\\left\\| u \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs) )
        .addOutput( Output(_name="norm_L2_divu",_latex="\\left\\| \\nabla \\cdot u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs) )
        .addOutput( Output(_name="norm_L2_p",_latex="\\left\\| p \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs2) );


    M_lambda = 1./(2.*mu) - math::sqrt( 1./(4.*mu*mu) + 4.*M_PI*M_PI);
    penalbc = this->vm()["bccoeff"].template as<value_type>();
    M_beta = this->vm()["beta"].template as<value_type>();
}


template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename PressureStabExpr>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::addPressureStabilisation( element_1_type& p,
                                                                 element_1_type& q,
                                                                 PressureStabExpr& p_stabexpr )
{

    if ( is_equal_order && this->vm()["stab-p"].template as<bool>() )
    {
        boost::timer t;
        Log() << "[assembly] add stabilisation terms for equal order approximation ( orderU="
              << OrderU << ", orderP=" << OrderP << " )\n";
        size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
        form2( Xh, Xh, D, _pattern=pattern )  +=
            integrate( internalfaces(mesh), _Q<2*OrderP+2>(),
                       (p_stabexpr)*(trans(jumpt(gradt(p)))*jump(grad(q))) );
        Log() << "[assembly] form2 D equal order stabilisation terms in " << t.elapsed() << "s\n"; t.restart();
    }

}

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename DivStabExpr>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::addDivergenceStabilisation( element_0_type& u,
                                                                   element_0_type& v,
                                                                   DivStabExpr& d_stabexpr )
{

    if ( this->vm()["stab-div"].template as<bool>() )
    {
        boost::timer t;
        Log() << "[assembly] add stabilisation terms for divergence ( orderU="
              << OrderU << ", orderP=" << OrderP << " )\n";
        size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
        form2( Xh, Xh, D, _pattern=pattern )  +=
            integrate( internalfaces(mesh),
                       (d_stabexpr)*(trans(jumpt(divt(u)))*jump(div(v))) );
        Log() << "[assembly] form2 D divergence stabilisation terms in " << t.elapsed() << "s\n"; t.restart();
    }

}

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::run()
{
    this->//addParameterValue( mu )
        addParameterValue( Dim )
        .addParameterValue( OrderU )
        .addParameterValue( OrderP )
        .addParameterValue( this->vm()["hsize"].template as<double>() );

    if (this->preProcessing() == RUN_EXIT) return;

    using namespace Feel::vf;

    boost::timer t;

    /*
     * First we create the mesh : a square [0,1]x[0,1] with characteristic
     * length = meshSize
     */
    Log() << "creating mesh with hsize=" << meshSize << "\n";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name="square",
                                         _shape="hypercube",
                                         _dim=Dim,
                                         _h=meshSize,
                                         _xmin=-0.5,_xmax=1.,
                                         _ymin=-0.5,_ymax=1.5 ) );

    Log() << "mesh created in " << t.elapsed() << "s\n"; t.restart();

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type E( Xh, "u" );
    element_type V( Xh, "v" );
    element_0_type ue = E.template element<0>();
    element_1_type pe = E.template element<1>();
    element_2_type le = E.template element<2>();
    element_0_type u = U.template element<0>();
    element_0_type v = V.template element<0>();
    element_1_type p = U.template element<1>();
    element_1_type q = V.template element<1>();
    element_2_type lambda = U.template element<2>();
    element_2_type nu = V.template element<2>();

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count("export") << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";
    Log() << "functionspace and elements created in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    vector_ptrtype F( M_backend->newVector( Xh ) );

    auto deft = .5*(gradt(u)+trans(gradt(u)));
    auto def = .5*(grad(v)+trans(grad(v)));

    // total stress tensor (trial)
    auto SigmaNt = (-idt(p)*N()+2*mu*deft*N());
    auto SigmaN = (-id(p)*N()+2*mu*def*N());

    //
    // the Kovasznay flow (2D)
    //
    // total stress tensor (test)
    double pi = M_PI;
    auto u1 = val(1. - exp( M_lambda * Px() ) * cos(2.*pi*Py()));
    auto u2 = val((M_lambda/(2.*pi)) * exp( M_lambda * Px() ) * sin(2.*pi*Py()));
    auto u_exact = vec(u1,u2);

    auto du_dx = val(-M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py()));
    auto du_dy = val(2*pi*exp( M_lambda * Px() )*sin(2.*pi*Py()));
    auto dv_dx = val((M_lambda*M_lambda/(2*pi))*exp( M_lambda * Px() )*sin(2.*pi*Py()));
    auto dv_dy = val(M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py()));
    auto grad_exact = (mat<2,2>(du_dx, du_dy, dv_dx, dv_dy));
    auto div_exact = du_dx + dv_dy;

    auto beta = vec(cst(M_beta),cst(M_beta));
    auto convection = grad_exact*beta;

    auto p_exact = val((1-exp(2.*M_lambda*Px()))/2.0);

    auto f1 = val(exp( M_lambda * Px() )*((M_lambda*M_lambda - 4.*pi*pi)*mu*cos(2.*pi*Py()) - M_lambda*exp( M_lambda * Px() )));
    auto f2 = val(exp( M_lambda * Px() )*mu*(M_lambda/(2.*pi))*sin(2.*pi*Py())*(-M_lambda*M_lambda +4*pi*pi));

    auto f = vec(f1,f2)+ convection;

    double meas = integrate( elements(mesh), constant(1.) ).evaluate()( 0, 0 );
    double pmean = integrate( elements(mesh), p_exact ).evaluate()( 0, 0 )/meas;

    v = vf::project( v.functionSpace(), elements(mesh), beta );
    ue = vf::project( ue.functionSpace(), elements(mesh), u_exact );
    Log() << "convection terms projectd in " << t.elapsed() << "s\n"; t.restart();

    // right hand side
    form1( Xh, F, _init=true ) =integrate( elements(mesh), trans(f)*id(v) );
    // impose pmean mean value for the solution through the lagrange multipliers
    form1( Xh, F ) += integrate( elements(mesh), pmean*id(nu) );
    if ( M_weak_dirichlet )
    {
        form1( Xh, F )+= integrate( boundaryfaces(mesh),
                                    trans(u_exact)*(-SigmaN+
                                                    penalbc*id(v)/hFace() +
                                                    penalbc*(trans(id(v))*N())*N()*
                                                    max(sqrt(trans(idv(v))*idv(v)),mu/hFace()) ) );
    }
    Log() << "[stokes] vector local assembly done\n";
    Log() << "form1 F created in " << t.elapsed() << "s\n"; t.restart();
    /*
     * Construction of the left hand side
     */
    D = sparse_matrix_ptrtype(  M_backend->newMatrix( Xh, Xh ) );
    size_type pattern = DOF_PATTERN_COUPLED;
    if ( (is_equal_order &&
          this->vm()["stab-p"].template as<bool>()) ||
         this->vm()["stab-div"].template as<bool>())
        pattern |= DOF_PATTERN_NEIGHBOR;
    Feel::Context graph( pattern );
    Log() << "[stokes] test : " << ( graph.test ( DOF_PATTERN_DEFAULT ) || graph.test ( DOF_PATTERN_NEIGHBOR ) ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_DEFAULT )=" <<  graph.test ( DOF_PATTERN_DEFAULT ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_COUPLED )=" <<  graph.test ( DOF_PATTERN_COUPLED ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_NEIGHBOR)=" <<  graph.test ( DOF_PATTERN_NEIGHBOR ) << "\n";
    Log() << "[assembly] add diffusion terms\n";
    form2( Xh, Xh, D, _init=true, _pattern=pattern );
    Log() << "[assembly] form2 D init in " << t.elapsed() << "s\n"; t.restart();
    form2( Xh, Xh, D )+= integrate( elements(mesh), 2*mu*trace(deft*trans(def)) + trans(gradt(u)*idv(v))*id(v) );
    Log() << "[assembly] form2 D convection and viscous terms in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[assembly] add velocity/pressure terms\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh),- div(v)*idt(p) + divt(u)*id(q) );
    Log() << "[assembly] form2 D velocity/pressure terms in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[assembly] add lagrange multipliers terms for zero mean pressure\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh), id(q)*idt(lambda) + idt(p)*id(nu) );
    Log() << "[assembly] form2 D pressure/multipliers terms in " << t.elapsed() << "s\n"; t.restart();
    if ( M_weak_dirichlet )
    {
        Log() << "[assembly] add terms for weak Dirichlet condition handling\n";
        form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), -trans(SigmaNt)*id(v) );
        form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), -trans(SigmaN)*idt(u) );
        form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), +penalbc*trans(idt(u))*id(v)/hFace() );
        form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), +penalbc*(trans(idt(u))*N())*(trans(id(v))*N())*max(sqrt(trans(idv(v))*idv(v)),mu/hFace()) );
        Log() << "[assembly] form2 D boundary terms in " << t.elapsed() << "s\n"; t.restart();
    }

    if ( math::abs( M_beta ) < 1e-10 )
    {
        double pterm = math::pow(double(OrderU), 7./2.);
        auto p_stabexpr = M_stabP*hFace()*hFace()/(pterm);
        this->addPressureStabilisation( p, q, p_stabexpr  );
    }
    else
    {
        double pterm = math::pow(double(OrderU), 4.);
        auto p_stabexpr = M_stabP*hFace()*hFace()*hFace()/max(mu*pterm,hFace()*sqrt(trans(idv(v))*idv(v)));
        this->addPressureStabilisation( p, q, p_stabexpr  );
    }
    auto d_stabexpr = M_stabD*hFace()*hFace()/math::pow(double(OrderU), 7./2.);
    this->addDivergenceStabilisation( u, v, d_stabexpr  );

    Log() << "[stokes] matrix local assembly done\n";

    D->close();
    F->close();
    if ( !M_weak_dirichlet )
    {
        form2( Xh, Xh, D ) += on( boundaryfaces( mesh ), u, F, u_exact );
    }
    Log() << "[stokes] vector/matrix global assembly done\n";
    Log() << "form2 D created in " << t.elapsed() << "s\n"; t.restart();

    if( this->vm().count( "export-matlab" ) )
    {
        D->printMatlab( "S.m" );
        F->printMatlab( "F.m" );
    }

    vector_ptrtype X( M_backend->newVector( U.functionSpace() ) );
    M_backend->solve( _matrix=D, _solution=X, _rhs=F, _rtolerance=1e-14 );
    U = *X;

    Log() << "system solved in " << t.elapsed() << "s\n"; t.restart();
    Log() << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );
    Log() << "postprocessing done in " << t.elapsed() << "s\n"; t.restart();
    this->exportResults( U, V );
    Log() << "exporting done in " << t.elapsed() << "s\n"; t.restart();



    this->postProcessing();

} // Stokes::run

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::exportResults( element_type& U, element_type& V )
{
    const double pi = M_PI;
    auto u1 = val(1. - exp( M_lambda * Px() ) * cos(2.*pi*Py()));
    auto u2 = val((M_lambda/(2.*pi)) * exp( M_lambda * Px() ) * sin(2.*pi*Py()));
    auto u_exact = vec(u1,u2);

    auto du_dx = val(-M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py()));
    auto du_dy = val(2*pi*exp( M_lambda * Px() )*sin(2.*pi*Py()));
    auto dv_dx = val((M_lambda*M_lambda/(2*pi))*exp( M_lambda * Px() )*sin(2.*pi*Py()));
    auto dv_dy = val(M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py()));
    auto grad_exact = (mat<2,2>(du_dx, du_dy, dv_dx, dv_dy));
    auto div_exact = du_dx + dv_dy;

    auto p_exact = val((1-exp(2.*M_lambda*Px()))/2.0);

    element_0_type u = U.template element<0>();
    element_1_type p = U.template element<1>();
    double u_error_L2_2 = integrate( elements(mesh),
                                     trans(idv(u)-u_exact)*(idv(u)-u_exact) ).evaluate()( 0, 0 );
    double uex_L2_2 = integrate( elements(mesh), trans(u_exact)*(u_exact) ).evaluate()( 0, 0 );
    double u_errorL2 = math::sqrt( u_error_L2_2/uex_L2_2 );
    std::cout << "||u_error||_0/||uex||_0 = " << u_errorL2 << "\n";;


    double u_errorsemiH1 = integrate( elements(mesh),
                                      trace((gradv(u)-grad_exact)*trans(gradv(u)-grad_exact))).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( u_error_L2_2+u_errorsemiH1 );
    double uex_semiH1_2 = integrate( elements(mesh), trace((grad_exact)*trans(grad_exact))).evaluate()( 0, 0 );
    double uex_H1 = math::sqrt( uex_L2_2+uex_semiH1_2 );
    double u_errorH1 = u_error_H1/uex_H1;
    std::cout << "||u_error||_1/||uex||_1 = " << u_errorH1 << "\n";


    double p_errorL2_2 = integrate( elements(mesh), (idv(p)-p_exact)*(idv(p)-p_exact) ).evaluate()( 0, 0 );
    double pex_L2 = integrate( elements(mesh), (p_exact)*(p_exact) ).evaluate()( 0, 0 );
    double p_errorL2 = math::sqrt( p_errorL2_2/pex_L2 );
    std::cout << "||p_error||_0/||pex||_0 = " <<  p_errorL2 << "\n";;

    Log() << "[stokes] solve for D done\n";

    double meas = integrate( elements(mesh), constant(1.) ).evaluate()( 0, 0 );
    double meanpexact = integrate( elements(mesh), p_exact ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 3)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 3)\n";

    double mean_p = integrate( elements(mesh), idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    Log() << "[stokes] mean(p_exact)=" << meanpexact << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p_exact)=" << meanpexact << "\n";

    double mean_div_u = integrate( elements(mesh), divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_errorL2_2 = integrate( elements(mesh), divv(u)*divv(u) ).evaluate()( 0, 0 );
    double uex_div_L2 = integrate( elements(mesh), div_exact ).evaluate()( 0, 0 );
    double uex_n_L2 = integrate( boundaryfaces(mesh), trans(u_exact)*N() ).evaluate()( 0, 0 );
    std::cout << "[stokes] ||div(uexact)||=" << uex_div_L2 << "\n";
    std::cout << "[stokes] ||uexact,n||=" << uex_n_L2 << "\n";
    double div_u_errorL2 = math::sqrt( div_u_errorL2_2 );
    Log() << "[stokes] ||div(u)||_2=" << div_u_errorL2 << "\n";
    std::cout << "[stokes] ||div(u)||=" << div_u_errorL2 << "\n";

    this->addOutputValue( u_errorL2 )
        .addOutputValue( u_errorH1 )
        .addOutputValue( div_u_errorL2 )
        .addOutputValue( p_errorL2 );

    if ( exporter->doExport() )
    {
        exporter->step(0)->setMesh( U.functionSpace()->mesh() );
        exporter->step(0)->add( "u", U.template element<0>() );
        exporter->step(0)->add( "p", U.template element<1>() );
        exporter->step(0)->add( "u_exact", V.template element<0>() );
        exporter->step(0)->add( "p_exact", V.template element<1>() );
        exporter->save();
    }
} // Stokes::export

}
