/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
            Goncalo Pena  <gpena@mat.uc.pt>
 Date: 02 Oct 2014

 Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_PRECONDITIONERBTPCD_HPP
#define FEELPP_PRECONDITIONERBTPCD_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/operatorpcd.hpp>

namespace Feel
{
template< typename space_type >
class PreconditionerBTPCD : public Preconditioner<typename space_type::value_type>
{
public:

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type velocity_space_type;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype pressure_space_ptrtype;

    typedef typename velocity_space_type::element_type velocity_element_type;
    typedef typename pressure_space_type::element_type pressure_element_type;

    typedef typename space_type::value_type value_type;

    static const uint16_type Dim = space_type::nDim;
    static const uint16_type uOrder = velocity_space_type::basis_type::nOrder;
    static const uint16_type pOrder = pressure_space_type::basis_type::nOrder;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef typename OperatorPCD<space_type>::type op_pcd_type;
    typedef typename OperatorPCD<space_type>::ptrtype op_pcd_ptrtype;
    
    /**
     * \param Xh velocity/pressure space type
     * \param bcFlags the boundary conditions flags
     * \param A the full matrix \f$(F B^T;B 0)\f$
     * \param nu viscosity
     * \param alpha mass term
     */
    PreconditionerBTPCD( space_ptrtype Xh, 
                         std::map<std::string, std::set<flag_type> > bcFlags, 
                         sparse_matrix_ptrtype A,
                         double nu, 
                         double alpha = 0 );

    PreconditionerBTPCD( const PreconditionerBTPCD& tc );

    void initialize();

    void assembleHelmholtz( double nu, double alpha = 0 );

    void assembleDivergence();

    void assembleGradient();

    void assembleSchurApp( double nu, double alpha = 0 );

    template< typename Expr_convection, typename Expr_bc >
    void update( Expr_convection const& expr_b, Expr_bc const& g );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        applyInverse( X, Y );
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;

    // other function
    int SetUseTranspose( bool UseTranspose)
    {
        return(false);
    }

    double NormInf() const
    {
        return 0;
    }

    const char * Label () const
    {
        return("Triangular Block Preconditioner");
    }

    bool UseTranspose() const
    {
        return(false);
    }

    bool HasNormInf () const
    {
        return(false);
    }

    virtual ~PreconditionerBTPCD(){};


private:

    backend_ptrtype M_b;

    space_ptrtype M_Xh;
    
    velocity_space_ptrtype M_Vh;
    pressure_space_ptrtype M_Qh;

    sparse_matrix_ptrtype M_A, M_helm, G, M_div;
    op_mat_ptrtype divOp, helmOp;

    mutable vector_ptrtype M_rhs, M_aux, M_vin,M_pin, M_vout, M_pout;

    mutable element_type U;
    velocity_element_type u, v;
    pressure_element_type p, q;

    op_pcd_ptrtype pcdOp;

    std::map<std::string, std::set<flag_type> > M_bcFlags;

    op_mat_ptrtype precHelm;
};






template < typename space_type >
PreconditionerBTPCD<space_type>::PreconditionerBTPCD( space_ptrtype Xh, 
                                                      std::map<std::string, std::set<flag_type> > bcFlags,
                                                      sparse_matrix_ptrtype A,
                                                      double nu, double alpha )
    :
    M_b(backend()),
    M_Xh( Xh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_A( A ),
    M_helm ( M_b->newMatrix( M_Vh, M_Vh ) ),
    G( M_b->newMatrix( M_Vh, M_Vh ) ),
    M_div ( M_b->newMatrix( _trial=M_Vh, _test=M_Qh ) ),
    M_rhs( M_b->newVector( M_Vh )  ),
    M_aux( M_b->newVector( M_Qh )  ),
    M_vin( M_b->newVector( M_Vh )  ),
    M_pin( M_b->newVector( M_Qh )  ),
    M_vout( M_b->newVector( M_Vh )  ),
    M_pout( M_b->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_bcFlags( bcFlags )
{
    LOG(INFO) << "Alpha: " << alpha << "\n";

    mpi::timer ti;
    initialize();

    LOG(INFO) << "Assemble Helmholtz Operator...\n";
    assembleHelmholtz( nu, alpha );
    LOG(INFO) << "[PreconditionerBTPCD] Constructor: Helmholtz assembly\n" ;

    LOG(INFO) << "Assemble Divergence Operator...\n";
    assembleDivergence();
    LOG(INFO) << "[PreconditionerBTPCD] Constructor: Divergence assembly time\n";

    LOG(INFO) << "Assemble Schur Operator...\n";
    assembleSchurApp(nu, alpha);
    LOG(INFO) << "[PreconditionerBTPCD] Constructor: Schur approximation assembly time:\n";
}



template < typename space_type >
PreconditionerBTPCD<space_type>::PreconditionerBTPCD( const PreconditionerBTPCD& tc )
    :
    M_Xh( tc.M_Xh ),
    M_Vh( tc.M_Vh ),
    M_Qh( tc.M_Qh ),
    M_helm( tc.M_helm ),
    G( tc.G ),
    M_div( tc.M_div ),
    divOp( tc.divOp ),
    helmOp( tc.helmOp ),
    M_rhs( tc.M_rhs ),
    M_aux( tc.M_aux ),
    M_vin( tc.M_vin ),
    M_pin( tc.M_pin ),
    M_vout( tc.M_vout ),
    M_pout( tc.M_pout ),
    u( tc.u ),
    v( tc.v ),
    q( tc.q ),
    pcdOp( tc.pcdOp ),
    M_bcFlags( tc.M_bcFlags ),
    precHelm( tc.precHelm )
{
    //Log() << "Call to Preconditioner BTPCD copy constructor...\n";
}

template < typename space_type >
void
PreconditionerBTPCD<space_type>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template < typename space_type >
void
PreconditionerBTPCD<space_type>::assembleHelmholtz( double nu, double alpha  )
{
    auto a = form2( _trial=M_Vh, _test=M_Vh, _matrix=M_helm );

    if ( alpha != 0 )
        a = integrate( _range=elements(M_Xh->mesh()),  _expr=alpha*trans(idt(u))*id(v) + nu*(trace(trans(gradt(u))*grad(v))) );
    else
        a = integrate( _range=elements(M_Xh->mesh()),  _expr=nu*(trace(trans(gradt(u))*grad(v))) );

    /*
    for( auto dir : M_bcFlags["Dirichlet"] )
        {
            a += integrate( markedfaces(M_Xh->mesh(), dir), - nu*trans(gradt(u)*N())*id(v) );
        }
     */
    M_helm->close();
}

template < typename space_type >
void
PreconditionerBTPCD<space_type>::assembleDivergence()
{
    auto a = form2( _trial=M_Vh, _test=M_Qh, _matrix=M_div );
    a = integrate( _range=elements(M_Xh->mesh()), _expr=divt( u ) * id( q ) );
    M_div->close();

    divOp = op( M_div, "DN");
}

template < typename space_type >
void
PreconditionerBTPCD<space_type>::assembleSchurApp( double nu, double alpha )
{
    LOG(INFO) << "Assembling schur complement";
    pcdOp = boost::make_shared<op_pcd_type>( M_Xh, M_bcFlags, nu, alpha );
    LOG(INFO) << "Assembling schur complement done";
}



template < typename space_type >
template< typename Expr_convection, typename Expr_bc >
void
PreconditionerBTPCD<space_type>::update( Expr_convection const& expr_b,
                                         Expr_bc const& g )
{
    static bool init_G = true;

    boost::timer ti;

    if ( !init_G )
        G->zero();

    auto lg = form2( _trial=M_Vh, _test=M_Vh, _matrix=G );

    lg = integrate( _range=elements(M_Xh->mesh()),
                    _expr=trans( gradt(u)*val(expr_b) )*id(v) );

    G->close();

    G->addMatrix( 1.0, M_helm );

    std::set<flag_type>::const_iterator diriIter;
    for( auto dir : M_bcFlags["Dirichlet"] )
    {
        lg += on( _range=markedfaces(M_Xh->mesh(), dir ), _element=u, _rhs=M_rhs, _expr=g.find(M_Xh->mesh()->markerName(dir))->second );
    }
    G->close();

    helmOp = op( G, "Fu" );

    ti.restart();
    pcdOp->update( expr_b, g );

    init_G = false;
}



template < typename space_type >
int
PreconditionerBTPCD<space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
    U = X;
    U.close();
#if 0
    Y.setZero();
    Y.add( 1., X );
    Y.close();
    return 1;
#endif

    LOG(INFO) << "Create velocity/pressure component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    if ( boption("btpcd.cd") )
    {
        LOG(INFO) << "velocity block : apply inverse convection diffusion...\n";
        helmOp->applyInverse(*M_vin, *M_vout);
    }
    else
    {
        *M_vout = *M_vin;
        M_vout->close();
    }

    LOG(INFO) << "pressure/velocity block : apply divergence...\n";
    divOp->apply( *M_vout, *M_pout );

    *M_aux = *M_pin;
    M_aux->close();
    M_aux->add( -1.0, *M_pout );
    M_aux->close();

    if ( boption("btpcd.pcd") )
    {
        LOG(INFO) << "pressure block: Solve for the pressure convection diffusion...\n";
        CHECK(pcdOp) << "Invalid PCD oeprator\n";
        CHECK(M_aux) << "Invalid aux vector\n";
        CHECK(M_pout) << "Invalid aux vector\n";

        pcdOp->applyInverse( *M_aux, *M_pout );
        LOG(INFO) << "pressure block: Solve for the pressure convection diffusion done\n";
    }
    else
    {
        *M_pout = *M_aux;
        M_pout->close();
    }
    LOG(INFO) << "Update output velocity/pressure...\n";

    U.template element<0>() = *M_vout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();

    return 0;
}
namespace meta
{
template< typename space_type >
struct btpcd
{
    typedef PreconditionerBTPCD<space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::btpcd<typename parameter::value_type<Args, tag::space>::type::element_type>::ptrtype ),
                                 btpcd,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( bc, *)
                                   ( matrix, *)
                                   )
                                 ( optional
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( nu,  *, doption("mu") )
                                   ( alpha, *, 0. )
                                   )
                                 )
{
    typedef typename meta::btpcd<typename parameter::value_type<Args, tag::space>::type::element_type>::ptrtype pbtpcd_t;
    typedef typename meta::btpcd<typename parameter::value_type<Args, tag::space>::type::element_type>::type btpcd_t;
    pbtpcd_t p( new btpcd_t( space, bc, matrix, nu, alpha ) );
    return p;
} // btcpd
} // Feel
#endif
