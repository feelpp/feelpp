/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Feb 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_STVENANTKIRCHHOFF_HPP
#define FEELPP_STVENANTKIRCHHOFF_HPP 1

#include <functional>
#include <tuple>
#include <feel/feel.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/newmark.hpp>


namespace Feel
{

enum ModelProperty { AXISYMM=1 };

po::options_description
stvenantkirchhoff_options(std::string prefix)
{
    po::options_description stvenantkirchhoffoptions( "StVenantKirchhoff problem options" );
    stvenantkirchhoffoptions.add_options()
        ( prefixvm( prefix, "young-modulus").c_str(), Feel::po::value<double>()->default_value( 1.4e6 ), "young-modulus" )
        ( prefixvm( prefix, "poisson-coeff").c_str(), Feel::po::value<double>()->default_value( 0.4 ), "poisson-coeff" )
        ( prefixvm( prefix, "rho").c_str(), Feel::po::value<double>()->default_value( 1000 ), "density [kg/m^3]" )
        ( prefixvm( prefix, "gravity").c_str(), Feel::po::value<std::string>()->default_value( "{0,0}" ), "gravity force expression" )
        ( prefixvm( prefix, "gravity-cst").c_str(), Feel::po::value<double>()->default_value( 2 ), "gravity-cst" )
        ( prefixvm( prefix, "verbose").c_str(), Feel::po::value<bool>()->default_value( true ), "verbose" )
        ;
    return stvenantkirchhoffoptions.add( backend_options(prefix) ).add( ts_options(prefix) );;
}

//po::options_description
//stvenantkirchhoff_options( std::string prefix );

template<typename DisplSpaceType, int prop=0>
class StVenantKirchhoff
{
public:
    using self_type = StVenantKirchhoff<DisplSpaceType,prop>;
    using displacement_space_type = typename mpl::if_<is_shared_ptr<DisplSpaceType>,
                                                      mpl::identity<typename DisplSpaceType::element_type>,
                                                      mpl::identity<DisplSpaceType>>::type::type;
    using mesh_type = typename displacement_space_type::mesh_type;
    using mesh_ptrtype = typename displacement_space_type::mesh_ptrtype;
    using displacement_space_ptrtype = boost::shared_ptr<displacement_space_type>;
    using displacement_type = typename displacement_space_type::element_type;
    using property_space_ptrtype = Pdh_ptrtype<mesh_type,0>;
    using properties = typename Pdh_type<mesh_type,0>::element_type;
    using exporter_ptrtype = boost::shared_ptr<Exporter<mesh_type>>;
    using ts_ptrtype = boost::shared_ptr<Newmark<displacement_space_type>>;
    StVenantKirchhoff() = delete;
    StVenantKirchhoff( std::string name, displacement_space_ptrtype Xh );

    static constexpr int dim = mesh_type::nDim;

    // @return the time stepping strategy
    ts_ptrtype ts() { return nm; }
    
    void setDensity( double r ) { rho = r; }
    void setGravityConstant( double g ) { gravity = g; }

    void setGravityForce( vector_field_expression<dim,1,2> e ) { gravityForce = e; }
    
    template < typename ExprT >
    void updateRho( Expr<ExprT> const& e )
        {
            P0Rho.on( _range=elements(mesh), _expr=e );
        }
    template < typename ExprT >
    void updateCoefflame1(Expr<ExprT> const& e )
        {
            P0Coefflame1.on(_range=elements(mesh), _expr=e );
        }
    template < typename ExprT >
    void updateCoefflame2(Expr<ExprT> const& e)
        {
            P0Coefflame2.on(_range=elements(mesh), _expr=e );
        }
    void init();

    using SolveData = SolverNonLinear<double>::SolveData;
    SolveData solve();

    displacement_type const& displacement() const { return u; }
    
    void exportResults();
    

    void updateResidualAxiSymm(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobianAxiSymm(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    
    void updateResidual(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);

private:
    std::string name;
    bool verbose;
    displacement_space_ptrtype Dh;
    property_space_ptrtype P0h;
    mesh_ptrtype mesh;
    displacement_type u;
    properties P0Rho, P0Coefflame1, P0Coefflame2;
    map_vector_field<dim,1,2> dirichlet_conditions;
    map_vector_field<dim,1,2> neumann_conditions;
    vector_field_expression<dim,1,2> gravityForce;
    double rho;
    double gravity;
    double youngmodulus, coeffpoisson;
    double coefflame1, coefflame2;

    ts_ptrtype nm;
    
    backend_ptrtype M_backend;
    vector_ptrtype Res;
    sparse_matrix_ptrtype Jac;


    exporter_ptrtype e;
};

template<typename DisplSpaceType, int props>
StVenantKirchhoff<DisplSpaceType,props>::StVenantKirchhoff( std::string n, displacement_space_ptrtype Xh )
    :
    name( n ),
    verbose( boption(_name=prefixvm(name,"verbose")) ),
    Dh( Xh ),
    P0h( Pdh<0>( Dh->mesh() ) ),
    mesh( Dh->mesh() ),
    u( Xh->element() ),
    P0Rho( P0h->element() ),
    P0Coefflame1( P0h->element() ),
    P0Coefflame2( P0h->element() ),
    nm( newmark( _space=Dh, _name=name,_rank_proc_in_files_name=true ) ),
    M_backend( backend( _name=name ) ),
    Res( M_backend->newVector( Dh ) ),
    Jac( M_backend->newMatrix( Dh, Dh ) ),
    e ( exporter( _mesh=mesh ) )
{
    tic();
    dirichlet_conditions = BoundaryConditionFactory::instance().getVectorFields<dim> ( "displacement", "Dirichlet" );
    neumann_conditions = BoundaryConditionFactory::instance().getVectorFields<dim> ( "displacement", "Neumann" );
    gravityForce = expr<dim,1,2>(soption(prefixvm(name,"gravity")));
    youngmodulus=doption(_name=prefixvm(name,"young-modulus"));
    coeffpoisson=doption(_name=prefixvm(name,"poisson-coeff"));
    coefflame2 = youngmodulus/(2*(1+coeffpoisson));// mu
    coefflame1 = youngmodulus*coeffpoisson/((1+coeffpoisson)*(1-2*coeffpoisson));// lambda
    rho=doption(_name=prefixvm(name,"rho"));
    gravity=doption(_name=prefixvm(name,"gravity-cst"));

    this->updateRho( cst(rho) );
    this->updateCoefflame1( cst(coefflame1 ) );
    this->updateCoefflame2( cst(coefflame2 ) );

    toc("StVenantKirchhoff constructor", verbose || FLAGS_v > 0 );
    
    
}
template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::init()
{
    tic();
    // start or restart
    if ( !nm->isRestart() )
    {
        nm->start();
        //svk.setInitialGuess();
    }
    else
    {
        double ti = nm->restart();
        u = nm->previousUnknown();
        if ( e->doExport() )
            e->restart(ti);
    }
    toc("StVenantKirchhoff::init", verbose || FLAGS_v > 0);
}
template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::updateResidualAxiSymm( const vector_ptrtype& X, vector_ptrtype& R )
{
    u = *X;

    auto Id = eye<dim,dim>();
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;

    auto r = form1( _test=Dh, _vector=R );
    r = integrate( _range=elements( mesh ),
                   _expr= inner( val(Fv*Sv) , grad(u) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= -inner(idv(P0Rho)*gravityForce,id( u ) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= rho*inner( nm->polyDerivCoefficient()*idv(u) -idv(nm->polyDeriv()),id( u ) ) );
    for( auto & n : neumann_conditions )
    {
        // update n with respect to the current time in case it depends on time
        expression(n).setParameterValues({{"t",nm->time()}});
        r += integrate( _range=markedfaces( mesh, marker(n) ),
                        _expr=trans(expression(n))*id(u) );
    }

    R->close();
    auto temp = Dh->element();
    temp = *R;
    for( auto const& d : dirichlet_conditions )
    {
        temp.on( _range=markedfaces( mesh, marker(d)), _expr=zero<dim,1>() );
    }
    *R = temp;

}
template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::updateJacobianAxiSymm(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    u = *X;
    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;
    auto dF = gradt(u);
    auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
    auto dS = idv(P0Coefflame1)*trace(dE)*Id + 2*idv(P0Coefflame2)*dE;
    
    auto a = form2( _test=Dh, _trial=Dh, _matrix=J );
    
    a = integrate( _range=elements( mesh ),
                   _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
    
    a += integrate( _range=elements( mesh ),
                    _expr= rho*inner( nm->polyDerivCoefficient()*idt(u),id( u ) ) );
    
    auto RR = backend()->newVector( Dh );
    for( auto const& d : dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u, _rhs=RR,
                 _expr=zero<dim,1>() );
    }
}

template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    u = *X;

    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;

    auto r = form1( _test=Dh, _vector=R );
    r = integrate( _range=elements( mesh ),
                   _expr= inner( val(Fv*Sv) , grad(u) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= -inner(idv(P0Rho)*gravityForce,id( u ) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= rho*inner( nm->polyDerivCoefficient()*idv(u) -idv(nm->polyDeriv()),id( u ) ) );
    for( auto & n : neumann_conditions )
    {
        // update n with respect to the current time in case it depends on time
        expression(n).setParameterValues({{"t",nm->time()}});
        r += integrate( _range=markedfaces( mesh, marker(n) ),
                        _expr=trans(expression(n))*id(u) );
    }

    R->close();
    auto temp = Dh->element();
    temp = *R;
    for( auto const& d : dirichlet_conditions )
    {
        temp.on( _range=markedfaces( mesh, marker(d)), _expr=zero<dim,1>() );
    }
    *R = temp;

}
template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    u = *X;
    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;
    auto dF = gradt(u);
    auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
    auto dS = idv(P0Coefflame1)*trace(dE)*Id + 2*idv(P0Coefflame2)*dE;
    
    auto a = form2( _test=Dh, _trial=Dh, _matrix=J );
    
    a = integrate( _range=elements( mesh ),
                   _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
    a += integrate( _range=elements( mesh ),
                    _expr= rho*inner( nm->polyDerivCoefficient()*idt(u),id( u ) ) );
    auto RR = backend()->newVector( Dh );
    for( auto const& d : dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u, _rhs=RR,
                 _expr=zero<dim,1>() );
    }
}

template<typename DisplSpaceType, int props>
typename StVenantKirchhoff<DisplSpaceType,props>::SolveData
StVenantKirchhoff<DisplSpaceType,props>::solve()
{
    tic();
    // make sure that the initial guess satisfies the boundary conditions
    for( auto const& d : dirichlet_conditions )
    {
        u.on( _range=markedfaces( mesh, marker(d)), _expr=expression(d));
    }

    using namespace std;
    using std::get;
    auto r = std::bind( &self_type::updateResidual, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);
    auto j = std::bind( &self_type::updateJacobian, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
    M_backend->nlSolver()->residual = r;
    M_backend->nlSolver()->jacobian = j;
        

    auto d = M_backend->nlSolve( _solution=u,_jacobian=Jac,_residual=Res );
    toc("StVenantKirchhoff::solve", verbose || FLAGS_v > 0 );
    return d;
}
template<typename DisplSpaceType, int props>
void
StVenantKirchhoff<DisplSpaceType,props>::exportResults()
{
    tic();
    if ( nm->isSteady() )
    {
        e->add( "displacement", u );
        e->save();
    }
    else
    {
        e->step(nm->time())->add( "displacement", u );
        e->step(nm->time())->add( "velocity", nm->currentVelocity() );
        e->step(nm->time())->add( "acceleration", nm->currentAcceleration() );
        e->save();
    }
    toc("StVenantKirchhoff::exportResults", verbose);
}
}

#endif
