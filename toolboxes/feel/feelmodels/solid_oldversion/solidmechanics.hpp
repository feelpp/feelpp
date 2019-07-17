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
#ifndef FEELPP_SOLIDMECHANICS_HPP
#define FEELPP_SOLIDMECHANICS_HPP 1

#include <functional>
#include <tuple>
#include <feel/feel.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feelmodels/modelproperties.hpp>


namespace Feel
{

enum ModelProperty { AXISYMM=1 };

po::options_description
solidmechanics_options(std::string prefix)
{
    po::options_description solidmechanicsoptions( "SolidMechanics problem options" );
    solidmechanicsoptions.add_options()
        ( prefixvm( prefix, "filename").c_str(), Feel::po::value<std::string>()->default_value( "" ), "json file describing model properties" )
        ( prefixvm( prefix, "model").c_str(), Feel::po::value<std::string>()->default_value( "linear" ), "linear, hyperelastic" )
        ( prefixvm( prefix, "gravity").c_str(), Feel::po::value<std::string>()->default_value( "{0,0}" ), "gravity force expression" )
        ( prefixvm( prefix, "gravity-cst").c_str(), Feel::po::value<double>()->default_value( 2 ), "gravity-cst" )
        ( prefixvm( prefix, "verbose").c_str(), Feel::po::value<bool>()->default_value( true ), "verbose" )
        ;
    return solidmechanicsoptions.add( backend_options(prefix) ).add( backend_options(prefix+".l2p") ).add( ts_options(prefix) );;
}

//po::options_description
//solidmechanics_options( std::string prefix );

template<typename DisplSpaceType>
class SolidMechanics
{
public:
    using self_type = SolidMechanics<DisplSpaceType>;
    using displacement_space_type = typename mpl::if_<is_shared_ptr<DisplSpaceType>,
                                                      mpl::identity<typename DisplSpaceType::element_type>,
                                                      mpl::identity<DisplSpaceType>>::type::type;
    using mesh_type = typename displacement_space_type::mesh_type;
    using mesh_ptrtype = typename displacement_space_type::mesh_ptrtype;
    using displacement_space_ptrtype = std::shared_ptr<displacement_space_type>;
    using displacement_type = typename displacement_space_type::element_type;
    static constexpr int dim = mesh_type::nDim;
    static constexpr int order = displacement_space_type::basis_type::nOrder;

    // stresses
    using stress_space_type = Pchm_type<mesh_type,order>;
    using stress_space_ptrtype = Pchm_ptrtype<mesh_type,order>;
    using stress_type = typename Pchm_type<mesh_type,order>::element_type;

    using equivalent_stress_space_type = Pch_type<mesh_type,order>;
    using equivalent_stress_space_ptrtype = Pch_ptrtype<mesh_type,order>;
    using equivalent_stress_type = typename Pch_type<mesh_type,order>::element_type;

    // material properties
    using property_space_ptrtype = Pdh_ptrtype<mesh_type,0>;
    using property_type = typename Pdh_type<mesh_type,0>::element_type;
    using properties = typename Pdh_type<mesh_type,0>::element_type;

    // exporter
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_type>>;

    // time stepping
    using ts_ptrtype = std::shared_ptr<Newmark<displacement_space_type>>;

    // Projection operator
    using l2_projector_ptrtype = projector_ptrtype<stress_space_type,stress_space_type>;
    
    SolidMechanics() = delete;
    SolidMechanics( std::string name, displacement_space_ptrtype Xh );

    // @return the time stepping strategy
    ts_ptrtype ts() { return nm; }
    
    void setGravityConstant( double g ) { gravity = g; }

    void setGravityForce( vector_field_expression<dim,1,2> e ) { gravityForce = e; }

    bool isLinear() const { return model=="linear"; }

    void setMaterialProperties()
        {
            for( auto const& m : props.materials() )
            {
                auto const& mat = m.second;
                auto const& matmarker = m.first;
                LOG(INFO) << "set material " << mat.name() << " associated to marker : " << matmarker<< "\n";
                P0Rho.on( _range=markedelements(mesh,matmarker), _expr=cst(mat.rho()) );
                //youngmodulus*coeffpoisson/((1+coeffpoisson)*(1-2*coeffpoisson));// lambda
                P0Coefflame1.on( _range=markedelements(mesh,matmarker), _expr=cst(mat.E()*mat.nu()/((1+mat.nu())*(1-2*mat.nu()) ) ) );
                // youngmodulus/(2*(1+coeffpoisson));// mu
                P0Coefflame2.on( _range=markedelements(mesh,matmarker), _expr=cst(mat.E()/(2*mat.nu()) ) );
            }
        }
    
    template < typename ExprT >
    void updateRho( Expr<ExprT> const& e )
        {
            P0Rho.on( _range=elements(mesh), _expr=e );
        }
    /**
     * @return the density of the material
     */
    property_type const&  density() const { return P0Rho; }
    
    /**
     * update the first Lame coefficient 
     * \f$Lambda = \frac{\nu * E}{(1+\nu)*(1-2*\nu)}\f$ 
     * where \f$E\f$ and \f$\nu\f$ are the Young
     * modulus and Poission's ratio respectively
     */
    template < typename ExprT >
    void updateCoefflame1(Expr<ExprT> const& e )
        {
            P0Coefflame1.on(_range=elements(mesh), _expr=e );
        }
    /**
     * @return the first Lame coefficient of the material \f$\lambda\f$
     */
    property_type const&  coeffLame1() const { return P0Coefflame1; }
    /**
     * update the second Lame coefficient 
     * \f$mu = \frac{E}{(2*(1+\nu)}\f$ 
     * where \f$E\f$ and \f$\nu\f$ are the Young
     * modulus and Poission's ratio respectively
     */
    template < typename ExprT >
    void updateCoefflame2(Expr<ExprT> const& e)
        {
            P0Coefflame2.on(_range=elements(mesh), _expr=e );
        }
    /**
     * @return the second Lame coefficient of the material \f$\mu\f$
     */
    property_type const&  coeffLame2() const { return P0Coefflame2; }

    /**
     * initialize the model
     */
    void init();

    /**
     * solve for the elasticity model
     */
    using SolveData = SolverNonLinear<double>::SolveData;
    SolveData solve();

    /**
     * @return the displacement
     */
    displacement_type const& displacement() const { return u; }

    /**
     * export the results from the model as defined in the json data file
     */
    void exportResults();
    void exportResults( double time );

    /**
     * @return compute the equivalent stress
     */

    void updateResidualAxiSymm(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobianAxiSymm(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    
    void updateResidual(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);

    /**
     * @return materials set
     */
    ModelMaterials const& materials() const { return props.materials(); }
    
    /**
     * @return material associated to mesh marker \c m
     */
    ModelMaterial const& material( std::string const& m ) const { return props.materials().material(m); }

    /**
     * @return the equivalent stress or von mises stress
     * \f$ \sqrt{\frac{3}{2} \sigma':\sigma'} \f$ where 
     * \f$\sigma'\f$ is the deviatoric stress tensor equal to
     * \f$\sigma'=\sigma - \frac{1}{3}tr(\sigma)I\f$ where \f$I\f$ 
     * is the identity matrix
     */
    equivalent_stress_type computeEquivalentStress() const;

    equivalent_stress_type const& equivalentStress() const { if (!Seh) this->computeEquivalentStress(); return sigma_vm; }
    /**
     * @return the equivalent stress space
     */
    equivalent_stress_space_ptrtype equivalentStressSpacePtr() const { if (!Seh) Seh = Pch<order>(mesh); return Seh; }
    
private:
    template<typename SpacePtrType>
    void initNullSpace( SpacePtrType Xh_vec, mpl::int_<2> )
        {
            K = nullspace_ptr( Xh_vec, oneX(), oneY(), vec(Py(),-Px()) );
        }
    template<typename SpacePtrType>
    void initNullSpace( SpacePtrType Xh_vec, mpl::int_<3> )
        {
            K = nullspace_ptr( Xh_vec, oneX(), oneY(), oneZ(),
                               vec(Py(),-Px(),cst(0.)),
                               vec(-Pz(),cst(0.),Px()),
                               vec(cst(0.),Pz(),-Py()) );
        }
    void addDirichlet( form2_type<displacement_space_type, displacement_space_type>& a,
                       form1_type<displacement_space_type>& l );
private:
    std::string name;
    ModelProperties props;

    std::string model;
    bool verbose;
    displacement_space_ptrtype Dh;
    property_space_ptrtype P0h;
    mesh_ptrtype mesh;
    displacement_type u;
    properties P0Rho, P0Coefflame1, P0Coefflame2;
    map_vector_field<dim,1,2> dirichlet_conditions;
    map_scalar_field<2> dx_dirichlet_conditions, dy_dirichlet_conditions, dz_dirichlet_conditions;
    map_vector_field<dim,1,2> neumann_conditions;
    vector_field_expression<dim,1,2> gravityForce;
    double gravity;
    double coefflame1, coefflame2;
    mutable stress_space_ptrtype Sh;
    mutable stress_type S;
    mutable equivalent_stress_space_ptrtype Seh;
    mutable equivalent_stress_type sigma_vm;
    mutable l2_projector_ptrtype l2p;
    ts_ptrtype nm;
    
    backend_ptrtype M_backend;
    vector_ptrtype Res;
    sparse_matrix_ptrtype Jac;
    std::shared_ptr<NullSpace<double>> K;

    
    exporter_ptrtype e,et;
};

template<typename DisplSpaceType>
SolidMechanics<DisplSpaceType>::SolidMechanics( std::string n, displacement_space_ptrtype Xh )
    :
    name( n ),
    props( Environment::expand( soption( _name=prefixvm(name,"filename")) ) ),
    model( soption( _name=prefixvm(name,"model")) ),
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
    e ( exporter( _mesh=mesh, _prefix=props.shortName(), _geo="static" ) ),
    et ( exporter( _mesh=mesh, _prefix=props.shortName()+"+t", _geo="static" ) )
    
{
    tic();
    dirichlet_conditions = props.boundaryConditions().template getVectorFields<dim> ( "displacement", "Dirichlet" );
    dx_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "displacement_x", "Dirichlet" );
    dy_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "displacement_y", "Dirichlet" );
    dz_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "displacement_z", "Dirichlet" );
    
    neumann_conditions = props.boundaryConditions().template getVectorFields<dim> ( "displacement", "Neumann" );
    gravityForce = expr<dim,1,2>(soption(prefixvm(name,"gravity")));
    gravity=doption(_name=prefixvm(name,"gravity-cst"));
    
    this->setMaterialProperties();

    initNullSpace( Dh, mpl::int_<dim>() );
    if ( !dirichlet_conditions.size() &&
         !dx_dirichlet_conditions.size() &&
         !dy_dirichlet_conditions.size() &&
         !dz_dirichlet_conditions.size() && nm->isSteady() )
        M_backend->attachNullSpace( toBackend( M_backend, K ) );
    if ( nm->isSteady() )
        M_backend->attachNearNullSpace( toBackend( M_backend, K ) );
    toc("SolidMechanics constructor", verbose || FLAGS_v > 0 );
    
    
}
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::init()
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
    toc("SolidMechanics::init", verbose || FLAGS_v > 0);
}
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::updateResidualAxiSymm( const vector_ptrtype& X, vector_ptrtype& R )
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
                    _expr= idv(P0Rho)*inner( nm->polyDerivCoefficient()*idv(u) -idv(nm->polyDeriv()),id( u ) ) );
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
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::updateJacobianAxiSymm(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
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
                    _expr= idv(P0Rho)*inner( nm->polyDerivCoefficient()*idt(u),id( u ) ) );
    
    auto RR = backend()->newVector( Dh );
    for( auto const& d : dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u, _rhs=RR,
                 _expr=zero<dim,1>() );
    }
}

template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    u = *X;
    

    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + (!isLinear())*(0.5*trans(gradv(u))*gradv(u));
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;

    auto r = form1( _test=Dh, _vector=R );
    r = integrate( _range=elements( mesh ),
                   _expr= inner( val(Fv*Sv) , grad(u) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= -inner(idv(P0Rho)*gravityForce,id( u ) ) );
    r += integrate( _range=elements( mesh ),
                    _expr= idv(P0Rho)*inner( nm->polyDerivCoefficient()*idv(u) -idv(nm->polyDeriv()),id( u ) ) );
    neumann_conditions.setParameterValues( props.parameters().toParameterValues() );
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
        temp.on( _range=markededges( mesh, marker(d)), _expr=zero<dim,1>() );
    }
    
    for( auto const& d : dx_dirichlet_conditions )
    {
        temp[Component::X].on( _range=markedfaces( mesh, marker(d)), _expr=cst(0.) );
        temp[Component::X].on( _range=markededges( mesh, marker(d)), _expr=cst(0.) );
    }
    
    for( auto const& d : dy_dirichlet_conditions )
    {
        temp[Component::Y].on( _range=markedfaces( mesh, marker(d)), _expr=cst(0.) );
        temp[Component::Y].on( _range=markededges( mesh, marker(d)), _expr=cst(0.) );
    }
    
    for( auto const& d : dz_dirichlet_conditions )
    {
        temp[Component::Z].on( _range=markedfaces( mesh, marker(d)), _expr=cst(0.) );
        temp[Component::Z].on( _range=markededges( mesh, marker(d)), _expr=cst(0.) );
    }
    *R = temp;

}
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::addDirichlet( form2_type<displacement_space_type, displacement_space_type>& a,
                                              form1_type<displacement_space_type>& l )
{
    for( auto const& d : dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u, _rhs=l,
                 _expr=zero<dim,1>() );
        a += on( _range=markededges( mesh, marker(d)), _element=u, _rhs=l,
                 _expr=zero<dim,1>() );
    }
    for( auto const& d : dx_dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u[Component::X], _rhs=l, _expr=cst(0.) );
        a += on( _range=markededges( mesh, marker(d)), _element=u[Component::X], _rhs=l, _expr=cst(0.) );
    }
    for( auto const& d : dy_dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u[Component::Y], _rhs=l, _expr=cst(0.) );
        a += on( _range=markededges( mesh, marker(d)), _element=u[Component::Y], _rhs=l, _expr=cst(0.) );
    }
    for( auto const& d : dz_dirichlet_conditions )
    {
        a += on( _range=markedfaces( mesh, marker(d)), _element=u[Component::Z], _rhs=l, _expr=cst(0.) );
        a += on( _range=markededges( mesh, marker(d)), _element=u[Component::Z], _rhs=l, _expr=cst(0.) );
    }

}
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    u = *X;

    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + (!isLinear())*(0.5*trans(gradv(u))*gradv(u));
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;
    auto dF = gradt(u);
    auto dE = sym(gradt(u)) + (!isLinear())*0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
    auto dS = idv(P0Coefflame1)*trace(dE)*Id + 2*idv(P0Coefflame2)*dE;
    
    auto a = form2( _test=Dh, _trial=Dh, _matrix=J );
    
    a = integrate( _range=elements( mesh ),
                   _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
    a += integrate( _range=elements( mesh ),
                    _expr= idv(P0Rho)*inner( nm->polyDerivCoefficient()*idt(u),id( u ) ) );
    auto RR = form1( _test=Dh );
    //auto RR = backend()->newVector( Dh );
    addDirichlet( a, RR );
    
}

template<typename DisplSpaceType>
typename SolidMechanics<DisplSpaceType>::SolveData
SolidMechanics<DisplSpaceType>::solve()
{
    tic();
    dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    //dirichlet_conditions.setParameterValues( { {"t",t->time()}}  );
    dx_dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    //dx_dirichlet_conditions.setParameterValues( { {"t",t->time()}}  );
    dy_dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    //dy_dirichlet_conditions.setParameterValues( { {"t",t->time()}}  );
    dz_dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    //dz_dirichlet_conditions.setParameterValues( { {"t",t->time()}}  );
    
    // make sure that the initial guess satisfies the boundary conditions
    for( auto const& d : dirichlet_conditions )
    {
        u.on( _range=markedfaces( mesh, marker(d)), _expr=expression(d));
        u.on( _range=markededges( mesh, marker(d)), _expr=expression(d));
    }
    for( auto const& d : dx_dirichlet_conditions )
    {
        u[Component::X].on( _range=markedfaces( mesh, marker(d)), _expr=expression(d) );
        u[Component::X].on( _range=markededges( mesh, marker(d)), _expr=expression(d) );
    }
    for( auto const& d : dy_dirichlet_conditions )
    {
        u[Component::Y].on( _range=markedfaces( mesh, marker(d)), _expr=expression(d) );
        u[Component::Y].on( _range=markededges( mesh, marker(d)), _expr=expression(d) );
    }
    for( auto const& d : dz_dirichlet_conditions )
    {
        u[Component::Z].on( _range=markedfaces( mesh, marker(d)), _expr=expression(d) );
        u[Component::Z].on( _range=markededges( mesh, marker(d)), _expr=expression(d) );
    }

    SolveData d;
    if ( isLinear() )
    {
        
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto a = form2( _test=Dh, _trial=Dh);
        a = integrate( _range=elements( mesh ),
                       _expr=idv(P0Coefflame1)*divt( u )*div( u )  +
                       2.*idv(P0Coefflame2)*trace( trans( deft )*def ) );
        auto l = form1(_test=Dh);
        neumann_conditions.setParameterValues( props.parameters().toParameterValues() );
        for( auto const& n : neumann_conditions )
        {
            l += integrate( _range=markedfaces( mesh, marker(n) ),
                            _expr=trans(expression(n))*id(u) );
        }
        addDirichlet( a, l );
        d = a.solve( _solution=u, _rhs=l );
    }
    else
    {
        using namespace std;
        using std::get;
        auto r = std::bind( &self_type::updateResidual, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);
        auto j = std::bind( &self_type::updateJacobian, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
        M_backend->nlSolver()->residual = r;
        M_backend->nlSolver()->jacobian = j;
        

        d = M_backend->nlSolve( _solution=u,_jacobian=Jac,_residual=Res );
    }
    toc("SolidMechanics::solve", verbose || FLAGS_v > 0 );
    return d;
}
template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::exportResults()
{
    tic();
    if ( nm->isSteady() )
    {
#if defined(FEELPP_HAS_HDF5)
        mesh->saveHDF5(props.shortName()+"_mesh.h5");
#endif
        
        for ( auto const& o : props.postProcess()["Fields"] )
        {
            if ( o == "displacement" )
            {
                e->add( "displacement", u );
#if defined(FEELPP_HAS_HDF5)
                u.saveHDF5(props.shortName()+"_displacement.h5");
#endif
            }
            if ( o == "lame" )
            {
                e->add( "lambda", P0Coefflame1 );
                e->add( "mu", P0Coefflame2 );
            }
            if ( o == "stress" )
            {
                tic();
                computeEquivalentStress();
                e->add( "S", S );
                e->add( "sigma_vm", sigma_vm );
#if defined(FEELPP_HAS_HDF5)
                //S.saveHDF5(props.shortName()+"_s.h5");
                sigma_vm.saveHDF5(props.shortName()+"_se.h5");
#endif
                
                toc("compute stress tensor");
            }
        }
        e->save();
    }
    else
    {
        et->step(nm->time())->add( "displacement", u );
        et->step(nm->time())->add( "velocity", nm->currentVelocity() );
        et->step(nm->time())->add( "acceleration", nm->currentAcceleration() );
        et->save();
    }
    toc("SolidMechanics::exportResults", verbose);
}

template<typename DisplSpaceType>
void
SolidMechanics<DisplSpaceType>::exportResults( double t )
{
    tic();
    for ( auto const& o : props.postProcess()["Fields"] )
    {
        if ( o == "displacement" )
        {
            et->step(t)->add( "displacement", u );
        }
        if ( o == "lame" )
        {
            et->add( "lambda", P0Coefflame1 );
            et->add( "mu", P0Coefflame2 );
        }
        if ( o == "stress" )
        {
            tic();
            computeEquivalentStress();
            et->add( "S", S );
            et->add( "sigma_vm", sigma_vm );
#if defined(FEELPP_HAS_HDF5)
            //S.saveHDF5(props.shortName()+"_s.h5");
            //sigma_vm.saveHDF5(props.shortName()+"_se.h5");
#endif
            
            toc("compute stress tensor");
        }
    }
    et->save();

    toc("SolidMechanics::exportResults", verbose);
}

template<typename DisplSpaceType>
typename SolidMechanics<DisplSpaceType>::equivalent_stress_type 
SolidMechanics<DisplSpaceType>::computeEquivalentStress() const
{
#if 1
    if ( !Sh )
    {
        Sh = Pchm<order>(mesh);
        S = Sh->element();
    }
    if ( !l2p )
    {
        //auto l2p = opProjection(_domainSpace=Sh,_imageSpace=Sh,_backend=backend(_prefix=name+".l2p"));
        l2p = opProjection(_domainSpace=Sh,_imageSpace=Sh);
    }
#endif
    if ( !Seh )
    {
        Seh = Pch<order>(mesh);
        sigma_vm = Seh->element();
        std::cout << "Seh done" << std::endl;
    }
    
    auto Id = eye<dim,dim>();    
    auto Fv = Id + gradv(u);
    auto Ev = sym(gradv(u)) + (!isLinear())*(0.5*trans(gradv(u))*gradv(u));
    auto Sv = idv(P0Coefflame1)*trace(Ev)*Id + 2*idv(P0Coefflame2)*Ev;
    std::cout << "expr done" << std::endl;
    S = l2p->projectL2( elements(mesh), Sv );
    //S.on( _range=elements(mesh), _expr=Sv );
    
    //auto sigma_dev = idv(S)-trace(idv(S))*Id/3;
    auto sigma_dev = Sv-trace(Sv)*Id/3;
    std::cout << "sigma done" << std::endl;
    sigma_vm.on( _range=elements(mesh), _expr=sqrt(3.*trace(sigma_dev*trans(sigma_dev))/2) );
    std::cout << "sigma_vm done" << std::endl;
    return sigma_vm;
}

} // Feel

#endif
