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
#ifndef FEELPP_PRECONDITIONERBlockNS_HPP
#define FEELPP_PRECONDITIONERBlockNS_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/operatorpcd.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{
template< typename space_type >
class PreconditionerBlockNS : public Preconditioner<typename space_type::value_type>
{
    typedef Preconditioner<typename space_type::value_type> super;
public:

    enum Type
    {
        PCD = 0, // pressure convection diffusion
        PMM=1, // pressure mass matrix
        SIMPLE=2 // 
    };
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
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
    typedef OperatorBase<value_type> op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;
    
    /**
     * \param Xh velocity/pressure space type
     * \param bcFlags the boundary conditions flags
     * \param A the full matrix \f$(F B^T;B 0)\f$
     * \param nu viscosity
     * \param alpha mass term
     */
    PreconditionerBlockNS( std::string t,
                           space_ptrtype Xh, 
                           BoundaryConditions bcFlags,
                           std::string const& s,
                           sparse_matrix_ptrtype A,
                           double nu, 
                           double alpha = 0 );
    
    Type type() const { return M_type; }
    void setType( std::string t );
    
    void initialize();

    void assembleHelmholtz( double nu, double alpha = 0 );

    void assembleDivergence();

    void assembleGradient();

    void assembleSchurApp( double nu, double alpha = 0 );

    template< typename Expr_convection, typename Expr_bc >
    void update( sparse_matrix_ptrtype A, Expr_convection const& expr_b, Expr_bc const& g );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        applyInverse( X, Y );
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;
    int guess( vector_type& U ) const;

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
        return("Triangular Blockns Preconditioner");
    }

    bool UseTranspose() const
    {
        return(false);
    }

    bool HasNormInf () const
    {
        return(false);
    }

    virtual ~PreconditionerBlockNS(){};


private:

    Type M_type;
    value_type M_nu, M_alpha;
    
    backend_ptrtype M_b;

    space_ptrtype M_Xh;
    
    velocity_space_ptrtype M_Vh;
    pressure_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;
    
    sparse_matrix_ptrtype M_helm, G, M_div, M_F, M_B, M_Bt, M_mass, M_massv_inv;
    op_mat_ptrtype divOp, helmOp;

    mutable vector_ptrtype M_rhs, M_aux, M_vin,M_pin, M_vout, M_pout;

    mutable element_type U;
    velocity_element_type u, v;
    pressure_element_type p, q;

    op_pcd_ptrtype pcdOp;
    op_ptrtype pm;
    
    BoundaryConditions M_bcFlags;
    std::string M_prefix;
    op_mat_ptrtype precHelm;
};






template < typename space_type >
PreconditionerBlockNS<space_type>::PreconditionerBlockNS( std::string t,
                                                          space_ptrtype Xh, 
                                                          BoundaryConditions bcFlags,
                                                          std::string const& p,
                                                          sparse_matrix_ptrtype A,
                                                          double nu, double alpha )
    :
    M_type( PCD ),
    M_nu( nu ),
    M_alpha( alpha ),
    M_b(backend()),
    M_Xh( Xh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_Vh_indices( M_Vh->nLocalDofWithGhost() ),
    M_Qh_indices( M_Qh->nLocalDofWithGhost() ),
    M_helm ( M_b->newMatrix( M_Vh, M_Vh ) ),
    G( M_b->newMatrix( M_Vh, M_Vh ) ),
    M_div ( M_b->newMatrix( _trial=M_Vh, _test=M_Qh ) ),
    M_mass( M_b->newMatrix( _trial=M_Qh, _test=M_Qh) ),
    M_massv_inv( M_b->newMatrix( _trial=M_Vh, _test=M_Vh) ),
    M_rhs( M_b->newVector( M_Vh )  ),
    M_aux( M_b->newVector( M_Vh )  ),
    M_vin( M_b->newVector( M_Vh )  ),
    M_pin( M_b->newVector( M_Qh )  ),
    M_vout( M_b->newVector( M_Vh )  ),
    M_pout( M_b->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_bcFlags( bcFlags ),
    M_prefix( p )
{
    tic();
    LOG(INFO) << "[PreconditionerBlockNS] setup starts";
    this->setMatrix( A );
    std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
    std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

    tic();
    M_F = A->createSubMatrix( M_Vh_indices, M_Vh_indices, true );
    M_B = A->createSubMatrix( M_Qh_indices, M_Vh_indices );
    M_Bt = A->createSubMatrix( M_Vh_indices, M_Qh_indices );
    toc( "BlockNS create sub matrix done", FLAGS_v > 0 );
    tic();
    helmOp = op( M_F, "Fu" );
    divOp = op( M_Bt, "Bt");
    toc( "BlockNS convection-diffusion and gradient operators done ", FLAGS_v > 0 );

    
    initialize();

    tic();    
    this->setType ( t );
    toc( "[PreconditionerBlockNS] setup done ", FLAGS_v > 0 );
}

template < typename space_type >
void
PreconditionerBlockNS<space_type>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template < typename space_type >
void
PreconditionerBlockNS<space_type>::setType( std::string t )
{
    if ( t == "PCD") M_type = PCD;
    if ( t == "PMM") M_type = PMM;
    if ( t == "SIMPLE") M_type = SIMPLE;

    LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
    switch( M_type )
    {
    case PCD:
        tic();
        pcdOp = boost::make_shared<op_pcd_type>( M_Xh, this->matrix(), M_b, M_bcFlags, M_prefix, M_nu, M_alpha );
        this->setSide( super::RIGHT );

        toc( "Assembling schur complement done", FLAGS_v > 0 );

        break;
    case PMM:
    {
        tic();
        auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
        m = integrate( elements(M_Qh->mesh()), idt(p)*id(q)/M_nu );
        M_mass->close();
        if ( boption( "blockns.pmm.diag" ) )
        {
            pm = diag( op( M_mass, "Mp" ) );
        }
        else
        {
            pm = op( M_mass, "Mp" );
        }
        this->setSide( super::RIGHT );
        
        toc( "Assembling pressure mass matrix schur complement appoximation done", FLAGS_v > 0 );
    }
    break;
    case SIMPLE:
        break;
    }
}

template < typename space_type >
template< typename Expr_convection, typename Expr_bc >
void
PreconditionerBlockNS<space_type>::update( sparse_matrix_ptrtype A,
                                         Expr_convection const& expr_b,
                                         Expr_bc const& g )
{

    tic();
    M_F = A->createSubMatrix( M_Vh_indices, M_Vh_indices, true );
        
    helmOp = op( M_F, "Fu" );
        
    toc("BlockNS convection-diffusion operator updated", FLAGS_v > 0 );
    if ( type() == PCD )
    {
        tic();
        pcdOp->update( expr_b, g );
        toc("BlockNS pressure convection-diffusion operator updated", FLAGS_v > 0 );
    }
}



template < typename space_type >
int
PreconditionerBlockNS<space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
    U = X;
    U.close();
    LOG(INFO) << "Create velocity/pressure component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();
    *M_aux = *M_vin;
    M_aux->close();
    
    if ( this->type() == PMM )
    {
        VLOG(2) << "applying mass matrix";
        pm->applyInverse( *M_pin, *M_pout );
        M_pout->scale(-1);
        M_pout->close();
    }
    if ( this->type() == PCD )
    {
        if ( boption("blockns.pcd") )
        {
            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion...\n";
            CHECK(pcdOp) << "Invalid PCD oeprator\n";
            CHECK(M_aux) << "Invalid aux vector\n";
            CHECK(M_pout) << "Invalid aux vector\n";
            
            pcdOp->applyInverse( *M_pin, *M_pout );
            M_pout->scale(-1);
            M_pout->close();
            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion done\n";
        }
        else
        {
            *M_pout = *M_pin;
            M_pout->close();
        }
    }
    LOG(INFO) << "pressure/velocity blockns : apply divergence...\n";
    divOp->apply( *M_pout, *M_vout );
    
    M_aux->add( -1.0, *M_vout );
    M_aux->close();

    if ( boption("blockns.cd") )
    {
        LOG(INFO) << "velocity blockns : apply inverse convection diffusion...\n";
        helmOp->applyInverse(*M_aux, *M_vout);
    }
    else
    {
        *M_vout = *M_vin;
        M_vout->close();
    }


    LOG(INFO) << "Update output velocity/pressure...\n";

    U.template element<0>() = *M_vout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();

    return 0;
}

template < typename space_type >
int
PreconditionerBlockNS<space_type>::guess ( vector_type& Y ) const
{
    U = Y;
    U.close();

    LOG(INFO) << "Create velocity/pressure component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    LOG(INFO) << "pressure/velocity blockns : apply divergence...\n";
    divOp->apply( *M_pout, *M_vin );

    M_aux->zero();
    M_aux->add( -1.0, *M_vin );
    M_aux->close();

    LOG(INFO) << "velocity blockns : apply inverse convection diffusion...\n";
    helmOp->applyInverse(*M_aux, *M_vin);
    LOG(INFO) << "Update output velocity/pressure...\n";

    U.template element<0>() = *M_vin;
    U.template element<1>() = *M_pin;
    U.close();
    Y=U;
    Y.close();

    return 0;
}
namespace meta
{
template< typename space_type >
struct blockns
{
    typedef PreconditionerBlockNS<space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type>::ptrtype ),
                                 blockns,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( matrix, *)
                                   ( type, *)
                                   )
                                 ( optional
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                   ( nu,  *, doption("mu") )
                                   ( alpha, *, 0. )
                                   )
                                 )
{
    typedef typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type>::ptrtype pblockns_t;
    typedef typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type>::type blockns_t;
    pblockns_t p( new blockns_t( type, space, bc, prefix, matrix, nu, alpha ) );
    return p;
} // btcpd
} // Feel
#endif
