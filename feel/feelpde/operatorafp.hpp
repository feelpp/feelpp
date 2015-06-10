/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
            Goncalo Pena  <gpena@mat.uc.pt>
 Date: 02 Oct 2014

 Copyright (C) 2014-2015 Feel++ Consortium

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
#ifndef FEELPP_OPERATORAFP_HPP
#define FEELPP_OPERATORAFP_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>

namespace Feel
{

template<typename space_type, typename mu_space_type>
class OperatorAFP : public OperatorBase<typename space_type::value_type>
{
    typedef OperatorBase<typename space_type::value_type> super;
public:

    typedef OperatorAFP<space_type, mu_space_type> type;
    typedef boost::shared_ptr<type> ptrtype;

    typedef typename space_type::value_type value_type;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<mu_space_type> mu_space_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type potential_space_type;
    typedef typename space_type::template sub_functionspace<1>::type lagrange_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype potential_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype lagrange_space_ptrtype;
    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;

    typedef typename mu_space_type::element_type mu_element_type;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef boost::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorCompose<op_inv_type, op_mat_type> comp1_type;
    typedef boost::shared_ptr<comp1_type> comp1_ptrtype;

    typedef OperatorCompose<op_mat_type, comp1_type> comp2_type;
    typedef boost::shared_ptr<comp2_type> comp2_ptrtype;

    typedef super op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;


    static const uint16_type Dim = space_type::nDim;
    static const uint16_type pOrder = lagrange_space_type::basis_type::nOrder;

    OperatorAFP( space_ptrtype Qh,
                 mu_space_ptrtype Mh, 
                 sparse_matrix_ptrtype A,
                 backend_ptrtype b1,
                 backend_ptrtype b2,
                 backend_ptrtype b3,
                 BoundaryConditions const& bcFlags,
                 std::string const& p
                 );

    OperatorAFP( const OperatorAFP& tc ) = default;
    OperatorAFP( OperatorAFP&& tc ) = default;
    OperatorAFP& operator=( const OperatorAFP& tc ) = default;
    OperatorAFP& operator=( OperatorAFP&& tc ) = default;

    void initialize();

    template < typename ExprConvection, typename ExprBC >
    void update( ExprConvection const& expr_b, ExprBC const& ebc );

    void setProblemType( std::string prob_type )
        {
            M_prob_type = prob_type;
        }

    std::string problemType() const
        {
            return M_prob_type;
        }

    ~OperatorAFP() {};

    int apply(const vector_type& X, vector_type& Y) const;
    int applyInverse(const vector_type& X, vector_type& Y) const;

private:
    backend_ptrtype M_amg1; // 14.a
    backend_ptrtype M_amg2; // 14.b
    backend_ptrtype M_amg3; // 5
    
    space_ptrtype M_Xh;
    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    mu_space_type M_Mh;

    potential_element_type u, v;
    lagrange_element_type p, q;
    mu_element_type mu;

    sparse_matrix_ptrtype M_mass, M_diff, M_conv, G, M_B, M_massv_inv, M_A;
    vector_ptrtype rhs;

    op_mat_ptrtype massOp, diffOp, convOp;

    BoundaryConditions M_bcFlags;
    std::string M_prefix;
    
    op_ptrtype precOp;

    std::string M_prob_type;

    void applyBC( sparse_matrix_ptrtype& A );
};




template < typename space_type, typename mu_space_type>
OperatorAFP<space_type, mu_space_type>::OperatorAFP( 
                                      space_ptrtype Qh,
                                      mu_space_ptrtype Mh,
                                      sparse_matrix_ptrtype A,
                                      backend_ptrtype b1,
                                      backend_ptrtype b2,
                                      backend_ptrtype b3,
                                      BoundaryConditions const& bcFlags,
                                      std::string const& p)
    :
    super( Qh->template functionSpace<1>()->mapPtr(), "AFP", false, false ),
    M_amg1( b1 ),
    M_amg2( b2 ),
    M_amg3( b3 ),
    M_Xh( Qh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_Mh( Mh ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    mu( M_Mh, "m" ),
    M_mass(      M_amg1->newMatrix(M_Qh, M_Qh) ),
    M_diff(      M_amg1->newMatrix(M_Qh, M_Qh) ),
    M_conv(      M_amg1->newMatrix(M_Qh, M_Qh) ),
    G(           M_amg1->newMatrix(M_Qh, M_Qh) ),
    M_B(         M_amg1->newMatrix(_trial=M_Vh, _test=M_Qh) ),
    M_massv_inv( M_amg1->newMatrix(_trial=M_Vh, _test=M_Vh) ),
    rhs(         M_amg1->newVector( M_Qh ) ),
    M_A( A ),
    M_bcFlags( bcFlags ),
    M_prefix( p )
{
    initialize();
#if 0
    this->assembleMass();
    this->assembleDiffusion();
#endif
}

template < typename space_type, typename mu_space_type>
void
OperatorAFP<space_type, mu_space_type>::initialize()
{
    rhs->zero();
    rhs->close();
}

template < typename space_type, typename mu_space_type>
template < typename ExprConvection, typename ExprBC >
void
OperatorAFP<space_type, mu_space_type>::update( ExprConvection const& expr_b,
                                 ExprBC const& ebc )
{
    tic();
    auto conv  = form2( _test=M_Qh, _trial=M_Qh, _matrix=G );
    G->zero();

    conv = integrate( _range=elements(M_Qh->mesh()), _expr=(trans(expr_b)*trans(gradt(p)))*id(q));
    conv += integrate( _range=elements(M_Qh->mesh()), _expr=idv(mu)*gradt(p)*trans(grad(q)));

    if ( soption("blockns.pcd.inflow") == "Robin" )
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
            LOG(INFO) << "Setting Robin condition on " << dir.marker();
            if ( ebc.find( dir.marker() ) != ebc.end() ){
                //conv += integrate( _range=markedfaces(M_Qh->mesh(), dir.marker()), _expr=-M_rho*trans(ebc.find(dir.marker())->second)*N()*idt(p)*id(q));
            }
        }

    G->close();

    this->applyBC(G);

    static bool init_G = false;

    //if ( !init_G )
    {
        // S = F G^-1 M
        LOG(INFO) << "[OperatorAFP] setting pcd operator...\n";
        if ( ioption("blockns.pcd.order") == 1 )
            precOp = compose( massOp, compose(inv(op(G,"Fp")),diffOp) );
        else
            precOp = compose( diffOp, compose(inv(op(G,"Fp")),massOp) );
        LOG(INFO) << "[OperatorAFP] setting pcd operator done.\n";
        init_G = true;
    }
    toc("Operator::AFP update",FLAGS_v>0);
}

#if 0



template < typename space_type, typename mu_space_type >
void
OperatorAFP<space_type,mu_space_type>::assembleMass()
{
    tic();
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(q) );
    M_mass->close();
    massOp = op( M_mass, "Mp" );
    toc("OperatorAFP::mass assembly",FLAGS_v>0);
}

template < typename space_type, typename mu_space_type >
void
OperatorAFP<space_type,mu_space_type>::assembleDiffusion()
{
    tic();
    if ( soption("blockns.pcd.diffusion") == "Laplacian" )
    {
        auto d = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_diff );
        d = integrate( _range=elements(M_Qh->mesh()), _expr=gradt(p)*trans(grad(q)));
        
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
            LOG(INFO) << "Diffusion Setting Dirichlet condition on lagrange on " << cond.marker();
            if ( boption("blockns.weakdir" ) )
                d+= integrate( markedfaces(M_Qh->mesh(),cond.marker()), _expr=-gradt(p)*N()*id(p)-grad(p)*N()*idt(p)+doption("penaldir")*idt(p)*id(p)/hFace() );
            else
                d += on( markedfaces(M_Qh->mesh(),cond.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
        //this->applyBC(M_diff);
    }

    diffOp = op( M_diff, "Ap" );
    toc("OperatorAFP::diffusion assembly",FLAGS_v>0);
}
#endif

template < typename space_type, typename mu_space_type >
void
OperatorAFP<space_type,mu_space_type>::applyBC( sparse_matrix_ptrtype& A )
{
#if 0
    tic();
    auto a = form2( _test=M_Qh, _trial=M_Qh, _matrix=A );

    if ( soption("blockns.pcd.inflow") != "Robin" )
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
            a += on( markedfaces(M_Qh->mesh(),dir.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }

    // on neumann boundary on potential, apply Dirichlet condition on lagrange
    if ( soption("blockns.pcd.outflow") == "Dirichlet" )
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
            if ( boption("blockns.weakdir" ) )
                a+= integrate( markedfaces(M_Qh->mesh(),cond.marker()), _expr=-idv(mu)*gradt(p)*N()*id(p)-idv(mu)*grad(p)*N()*idt(p)+doption("penaldir")*idt(p)*id(p)/hFace() );
            else
                a += on( markedfaces(M_Qh->mesh(),cond.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
    rhs->close();
    A->close();
    toc("OperatorAFP::BC apply",FLAGS_v>0);
#endif
}

template < typename space_type, typename mu_space_type >
int
OperatorAFP<space_type,mu_space_type>::apply(const vector_type& X, vector_type& Y) const
{
    return precOp->apply( X, Y );
}
template < typename space_type, typename mu_space_type >
int
OperatorAFP<space_type,mu_space_type>::applyInverse(const vector_type& X, vector_type& Y) const
{
    return precOp->applyInverse( X, Y );
}


} // Feel

#endif
