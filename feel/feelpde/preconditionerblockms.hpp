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
#ifndef FEELPP_PRECONDITIONERBlockMS_HPP
#define FEELPP_PRECONDITIONERBlockMS_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/operatorafp.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{
template< typename space_type, typename mu_space_type >
class PreconditionerBlockMS : public Preconditioner<typename space_type::value_type>
{
    typedef Preconditioner<typename space_type::value_type> super;
public:

    enum Type
    {
        AFP    = 0, // augmentation free preconditioner
        SIMPLE = 2 // 
    };
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<mu_space_type> mu_space_ptrtype;
    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type potential_space_type;
    typedef typename space_type::template sub_functionspace<1>::type lagrange_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype potential_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype lagrange_space_ptrtype;
    
    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;
    typedef typename mu_space_type::element_type element_mu_type;

    typedef typename space_type::value_type value_type;

    static const uint16_type Dim = space_type::nDim;
    static const uint16_type aOrder = potential_space_type::basis_type::nOrder;
    static const uint16_type pOrder = lagrange_space_type::basis_type::nOrder;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef typename OperatorAFP<space_type,mu_space_type>::type op_afp_type;
    typedef typename OperatorAFP<space_type,mu_space_type>::ptrtype op_afp_ptrtype;
    typedef OperatorBase<value_type> op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;
    
    /**
     * \param t Kind of prec
     * \param Xh potential/lagrange space type
     * \param Mh Permeability space type
     * \param bcFlags the boundary conditions flags
     * \param s name of backend
     * \param A the full matrix 
     */
    PreconditionerBlockMS( std::string t,
                           space_ptrtype Xh,
                           mu_space_ptrtype Mh, 
                           BoundaryConditions bcFlags,
                           std::string const& s,
                           sparse_matrix_ptrtype A);
    
    Type type() const { return M_type; }
    void setType( std::string t );
    
    void initialize();
    
    template< typename Expr_convection, typename Expr_bc >
    void update( sparse_matrix_ptrtype A, element_mu_type mu );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        // résoudre équation 11
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;
    int guess( vector_type& U ) const;

    virtual ~PreconditionerBlockMS(){};

private:
    void createSubMatrices();
private:

    Type M_type;
    
    backend_ptrtype M_cg1; // 12
    backend_ptrtype M_cg2; // 16

    space_ptrtype M_Xh;
    mu_space_ptrtype M_Mh;
    
    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;
    
    sparse_matrix_ptrtype M_Pm, P_p, M_l, M_q, M_c; // 10
    op_mat_ptrtype divOp, helmOp;

    mutable vector_ptrtype M_rhs, M_aux, M_vin,M_pin, M_vout, M_pout;

    mutable element_type U;
    potential_element_type u, v;
    lagrange_element_type p, q;
    element_mu_type mu;

    op_afp_ptrtype afpOp;
    op_ptrtype pm;
    
    BoundaryConditions M_cg1cFlags;
    std::string M_prefix;
    op_mat_ptrtype precHelm;
};






template < typename space_type, typename mu_space_type >
PreconditionerBlockMS<space_type,mu_space_type>::PreconditionerBlockMS( std::string t,
                                                          space_ptrtype Xh, 
                                                          mu_space_ptrtype Mh, 
                                                          BoundaryConditions bcFlags,
                                                          std::string const& p,
                                                          sparse_matrix_ptrtype A )
    :
    M_type( AFP ),
    mu( mu ),
    M_cg1(backend()),
    M_cg2(backend()),
    M_Xh( Xh ),
    M_Mh( Mh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_Vh_indices( M_Vh->nLocalDofWithGhost() ),
    M_Qh_indices( M_Qh->nLocalDofWithGhost() ),
    M_rhs( M_cg1->newVector( M_Vh )  ),
    M_aux( M_cg1->newVector( M_Vh )  ),
    M_vin( M_cg1->newVector( M_Vh )  ),
    M_pin( M_cg1->newVector( M_Qh )  ),
    M_vout( M_cg1->newVector( M_Vh )  ),
    M_pout( M_cg1->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_cg1cFlags( bcFlags ),
    M_prefix( p )
{
    tic();
    LOG(INFO) << "[PreconditionerBlockMS] setup starts";
    this->setMatrix( A );
    std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
    std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

    this->createSubMatrices();

    
    initialize();

    tic();    
    this->setType ( t );
    toc( "[PreconditionerBlockMS] setup done ", FLAGS_v > 0 );
}

template < typename space_type, typename mu_space_type >
void
PreconditionerBlockMS<space_type,mu_space_type>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template < typename space_type, typename mu_space_type >
void
PreconditionerBlockMS<space_type,mu_space_type>::createSubMatrices()
{
#if 0
    tic();
    M_F = this->matrix()->createSubMatrix( M_Vh_indices, M_Vh_indices, true );
    M_B = this->matrix()->createSubMatrix( M_Qh_indices, M_Vh_indices );
    M_Bt = this->matrix()->createSubMatrix( M_Vh_indices, M_Qh_indices );
    helmOp = op( M_F, "Fu" );
    divOp = op( M_Bt, "Bt");
    toc( "PreconditionerBlockMS::createSubMatrix(Fu,B^T)", FLAGS_v > 0 );
#endif

}
template < typename space_type, typename mu_space_type >
void
PreconditionerBlockMS<space_type,mu_space_type>::setType( std::string t )
{
    if ( t == "AFP") M_type = AFP;
    if ( t == "SIMPLE") M_type = SIMPLE;
#if 0
    LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
    switch( M_type )
    {
    case AFP:
        tic();
        afpOp = boost::make_shared<op_afp_type>( M_Xh, this->matrix(), M_cg1, M_cg1cFlags, M_prefix, M_mu, M_rho, M_alpha );
        this->setSide( super::LEFT );
        toc( "Preconditioner::setType AFP", FLAGS_v > 0 );
        break;
    break;
    case SIMPLE:
        break;
    }
#endif
}

template < typename space_type, typename mu_space_type >
template< typename Expr_convection, typename Expr_bc >
void
PreconditionerBlockMS<space_type,mu_space_type>::update( sparse_matrix_ptrtype A,element_mu_type mu )
{
#if 0
    tic();
    this->setMatrix( A );
    this->createSubMatrices();
    
    if ( type() == AFP )
    {
        tic();
        afpOp->update( expr_b, g );
        toc( "Preconditioner::update AFP", FLAGS_v > 0 );
        
    }
    toc( "Preconditioner::update", FLAGS_v > 0 );
#endif
}



template < typename space_type, typename mu_space_type >
int
PreconditionerBlockMS<space_type,mu_space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
#if 0
    tic();
    U = X;
    U.close();
    LOG(INFO) << "Create potential/lagrange component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();
    *M_aux = *M_vin;
    M_aux->close();
    
    if ( this->type() == AFP )
    {
        if ( boption("blockms.afp") )
        {
            LOG(INFO) << "lagrange blockms: Solve for the lagrange convection diffusion...\n";
            CHECK(afpOp) << "Invalid AFP oeprator\n";
            CHECK(M_aux) << "Invalid aux vector\n";
            CHECK(M_pout) << "Invalid aux vector\n";
            tic();
            afpOp->applyInverse( *M_pin, *M_pout );
            M_pout->scale(-1);
            M_pout->close();
            toc("PreconditionerBlockMS::applyInverse AFP::S^-1",FLAGS_v>0);
            
            LOG(INFO) << "lagrange blockms: Solve for the lagrange convection diffusion done\n";
        }
        else
        {
            *M_pout = *M_pin;
            M_pout->close();
        }
    }
    LOG(INFO) << "lagrange/potential blockms : apply divergence...\n";
    tic();
    divOp->apply( *M_pout, *M_vout );
    
    
    M_aux->add( -1.0, *M_vout );
    M_aux->close();
    toc("PreconditionerBlockMS::applyInverse apply B^T",FLAGS_v>0);
    
    if ( boption("blockms.cd") )
    {
        
        LOG(INFO) << "potential blockms : apply inverse convection diffusion...\n";
        tic();
        helmOp->applyInverse(*M_aux, *M_vout);
        toc("PreconditionerBlockMS::applyInverse Fu^-1",FLAGS_v>0);
    }
    else
    {
        *M_vout = *M_vin;
        M_vout->close();
    }


    LOG(INFO) << "Update output potential/lagrange...\n";
    tic();
    U.template element<0>() = *M_vout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();
    toc("PreconditionerBlockMS::applyInverse update solution",FLAGS_v>0);
    toc("PreconditionerBlockMS::applyInverse" );
#endif
    return 0;
}

template < typename space_type, typename mu_space_type >
int
PreconditionerBlockMS<space_type,mu_space_type>::guess ( vector_type& Y ) const
{
#if 0
  U = Y;
    U.close();

    LOG(INFO) << "Create potential/lagrange component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    LOG(INFO) << "lagrange/potential blockms : apply divergence...\n";
    divOp->apply( *M_pout, *M_vin );

    M_aux->zero();
    M_aux->add( -1.0, *M_vin );
    M_aux->close();

    LOG(INFO) << "potential blockms : apply inverse convection diffusion...\n";
    helmOp->applyInverse(*M_aux, *M_vin);
    LOG(INFO) << "Update output potential/lagrange...\n";

    U.template element<0>() = *M_vin;
    U.template element<1>() = *M_pin;
    U.close();
    Y=U;
    Y.close();
#endif
    return 0;
}
namespace meta
{
template< typename space_type , typename mu_space_type >
struct blockms
{
    typedef PreconditionerBlockMS<space_type, mu_space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::blockms<
      typename parameter::value_type<Args, tag::space>::type::element_type,
      typename parameter::value_type<Args, tag::space2>::type::element_type>::ptrtype ),
                                 blockms,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( space2, *)
                                   ( matrix, *)
                                   )
                                 ( optional
                                   ( type, *, "AFP")
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                   )
                                 )
{
#if 0
    typedef typename meta::blockms<typename parameter::value_type<Args, tag::space>::type::element_type>::ptrtype pblockms_t;
    typedef typename meta::blockms<typename parameter::value_type<Args, tag::space>::type::element_type>::type blockms_t;
    pblockms_t p( new blockms_t( type, space, space2, bc, prefix, matrix ) );
    return p;
#endif
} // btcpd
} // Feel
#endif
