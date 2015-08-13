/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

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

#if FM_DIM == 2
#define curl_op curlx
#define curlt_op curlxt
#define curlv_op curlxv
#else
#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/preconditioneras.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{
template< typename space_type, typename coef_space_type >
    class PreconditionerBlockMS : public Preconditioner<typename space_type::value_type>
{
    typedef Preconditioner<typename space_type::value_type> super;
public:

    enum Type
    {
        AFP    = 0, // augmentation free preconditioner
        SIMPLE = 2 // 
    };
    typedef typename backend_type::solve_return_type solve_return_type;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<coef_space_type> coef_space_ptrtype;
    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type potential_space_type;
    typedef typename space_type::template sub_functionspace<1>::type lagrange_space_type;
    typedef typename boost::shared_ptr<potential_space_type> potential_space_ptrtype;
    typedef typename boost::shared_ptr<lagrange_space_type> lagrange_space_ptrtype;

    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;
    typedef typename coef_space_type::element_type element_coef_type;

    typedef typename space_type::value_type value_type;

    static const uint16_type Dim = space_type::nDim;

    typedef PreconditionerAS<space_type, coef_space_type> pc_as_type;
    typedef boost::shared_ptr<pc_as_type> pc_as_ptrtype;

    typedef OperatorMatrix<value_type> op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    /**
     * \param t Kind of prec (Simple or AFP)
     * \param Xh potential/lagrange space type
     * \param Mh Permeability space type
     * \param bcFlags the boundary conditions flags
     * \param s name of backend
     * \param A the full matrix
     */
    PreconditionerBlockMS( std::string t,
                           space_ptrtype Xh,
                           coef_space_ptrtype Mh,
                           BoundaryConditions bcFlags,
                           std::string const& s,
                           sparse_matrix_ptrtype A);

    Type type() const { return M_type; }

    void setType( std::string t );

    void update( sparse_matrix_ptrtype A, element_coef_type mu );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        this->applyInverse(X,Y);
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;

    virtual ~PreconditionerBlockMS(){};

    void printMatSize(int i, std::ostream &os)
    {
        if(i == 1)
        os << M_11->graph()->size();
        else if(i==2)
        os << M_L->graph()->size();
    }
    void printIter(int i, std::ostream & os){
        os << "<ul>";
        if(i == 1)
            for(auto i : NbIter1)
                os << "<li>" << i.nIterations() << " -- " << i.residual() << " -- " << i.isConverged() << "</li>";
        else if(i == 2)
            for(auto i : NbIter2)
                os << "<li>" << i.nIterations() << " -- " << i.residual() << " -- " << i.isConverged() << "</li>";
        else
            os << "<li>**error**</li>" << std::endl;
        os << "</ul>";
    }

private:
    Type M_type;

    backend_ptrtype M_backend;
    mutable std::vector<solve_return_type> NbIter1, NbIter2;
    space_ptrtype M_Xh;
    coef_space_ptrtype M_Mh;

    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;

    mutable vector_ptrtype M_uin;
    mutable vector_ptrtype M_uout;
    mutable vector_ptrtype M_pin;
    mutable vector_ptrtype M_pout;

    mutable element_type U;

    sparse_matrix_ptrtype M_11;
    sparse_matrix_ptrtype M_mass;
    sparse_matrix_ptrtype M_L;

    element_coef_type M_er; // permittivity

    op_ptrtype M_22Op;
    op_ptrtype M_11Op; // if not augmented spaces

    value_type M_k; // wave number

    BoundaryConditions M_bcFlags;
    std::string M_prefix;

    potential_element_type u;
    lagrange_element_type phi;

    pc_as_ptrtype M_pcAs;

};

template < typename space_type, typename coef_space_type >
PreconditionerBlockMS<space_type,coef_space_type>::PreconditionerBlockMS(std::string t,                // Type
                                                                         space_ptrtype Xh,             // (u)x(p)
                                                                         coef_space_ptrtype Mh,        // mu
                                                                         BoundaryConditions bcFlags,   // bc
                                                                         std::string const& p,         // prefix
                                                                         sparse_matrix_ptrtype AA )     // The matrix
:
    M_type( AFP ),
    M_backend(backend()),           // the backend associated to the PC
    M_Xh( Xh ),
    M_Mh( Mh ),
    M_Vh( Xh->template functionSpace<0>() ),
    M_Qh( Xh->template functionSpace<1>() ),
    M_Vh_indices( M_Vh->nLocalDofWithGhost() ),
    M_Qh_indices( M_Qh->nLocalDofWithGhost() ),
    M_uin( M_backend->newVector( M_Vh )  ),
    M_uout( M_backend->newVector( M_Vh )  ),
    M_pin( M_backend->newVector( M_Qh )  ),
    M_pout( M_backend->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    M_mass(M_backend->newMatrix(M_Vh,M_Vh)),
    M_L(M_backend->newMatrix(M_Qh,M_Qh)),
    M_er( M_Mh, "er" ),
    M_k(0.),
    M_bcFlags( bcFlags ),
    M_prefix( p ),
    u(M_Vh, "u"),
    phi(M_Qh, "phi")
{
    tic();
    LOG(INFO) << "[PreconditionerBlockMS] setup starts";
    this->setMatrix( AA ); // Needed only if worldComm > 1

    this->setType ( t );

    if(this->type() == AFP){
        //this->setMatrix( AA );
        /* Indices are need to extract sub matrix */
        std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
        std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

        M_11 = AA->createSubMatrix( M_Vh_indices, M_Vh_indices, true);

        map_vector_field<FM_DIM,1,2> m_dirichlet_u { M_bcFlags.getVectorFields<FM_DIM> ( "u", "Dirichlet" ) };
        map_scalar_field<2> m_dirichlet_p { M_bcFlags.getScalarFields<2> ( "phi", "Dirichlet" ) };

        /* Compute the mass matrix */
        auto f2A = form2(_test=M_Vh, _trial=M_Vh, _matrix=M_mass);
        auto f1A = form1(_test=M_Vh);
        f2A = integrate(_range=elements(M_Vh->mesh()), _expr=inner(idt(u),id(u))); // M
        for(auto const & it : m_dirichlet_u )
        {
            LOG(INFO) << "Applying " << it.second << " on " << it.first << " for blockms.11\n";
            f2A += on(_range=markedfaces(M_Vh->mesh(),it.first), _expr=it.second,_rhs=f1A, _element=u, _type=soption("blockms.11.on.type"));
        }

        /* Compute the L matrix */
        auto f2B = form2(_test=M_Qh,_trial=M_Qh, _matrix=M_L);
        auto f1B = form1(_test=M_Qh);
        f2B = integrate(_range=elements(M_Qh->mesh()), _expr=inner(gradt(phi), grad(phi)));
        for(auto const & it : m_dirichlet_p)
        {
            LOG(INFO) << "Applying " << it.second << " on " << it.first << " for blockms.22\n";
            f2B += on(_range=markedfaces(M_Qh->mesh(),it.first),_element=phi, _expr=it.second, _rhs=f1B, _type=soption("blockms.22.on.type"));
        }

        /* Initialize the blockAS prec */
        if(soption("blockms.11.pc-type") == "AS")
        {
            M_pcAs = blockas(_space=M_Xh,
                             _space2=M_Mh,
                             _bc = M_bcFlags);
        }
    }
    toc( "[PreconditionerBlockMS] setup done ", FLAGS_v > 0 );
}

template < typename space_type, typename coef_space_type >
    void
PreconditionerBlockMS<space_type,coef_space_type>::setType( std::string t )
{
    if ( t == "AFP") M_type = AFP;
    else if ( t == "SIMPLE") M_type = SIMPLE;
    LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
}

template < typename space_type, typename coef_space_type >
//template< typename Expr_convection, typename Expr_bc >
    void
PreconditionerBlockMS<space_type,coef_space_type>::update( sparse_matrix_ptrtype A, element_coef_type mu )
{
    tic();
    this->setMatrix( A );
    if(this->type() == AFP){
        M_er.on(_range=elements(M_Mh->mesh()), _expr=cst(1.));;

        LOG(INFO) << "Create sub Matrix\n";

        /*
         * AA = [[ A - k^2 M, B^t],
         *      [ B        , 0  ]]
         * We need to extract A-k^2 M and add it M to form A+(1-k^2) M = A+g M
         */
        // Is the zero() necessary ?
        if(boption("blockms.rebuild_11"))
        {
             // calculer matrice A + g M
            map_vector_field<FM_DIM,1,2> m_dirichlet_u { M_bcFlags.getVectorFields<FM_DIM> ( "u", "Dirichlet" ) };
            auto f2A = form2(_test=M_Vh, _trial=M_Vh,_matrix=M_11);
            auto f1A = form1(_test=M_Vh);
            f2A = integrate(_range=elements(M_Vh->mesh()), _expr=cst(1.)/idv(mu)*trans(curlt_op(u))*curl_op(u) // mu^-1 A
                                                                +cst(1.-M_k*M_k)*idv(M_er)*inner(idt(u),id(u))); // g M
            for(auto const & it : m_dirichlet_u )
            {
                LOG(INFO) << "Applying " << it.second << " on " << it.first << " for blockms.22\n";
                f2A += on(_range=markedfaces(M_Vh->mesh(),it.first), _expr=it.second,_rhs=f1A, _element=u, _type=soption("blockms.11.on.type"));
            }
        }
        else
        {
            M_11->zero();
            A->updateSubMatrix( M_11, M_Vh_indices, M_Vh_indices); // M_11 = A-k^2 %
            M_11->addMatrix(1.0,M_mass);                           // A-k^2 M + M = A+(1-k^2) M
        }

        M_11Op = op(M_11, "blockms.11");

        // Is the zero() necessary ?
        //M_22->zero();
        //M_22->addMatrix(1.0,M_L);
        //auto f2B = form2(_test=M_Qh,_trial=M_Qh, _matrix=M_22);
        //auto f1B = form1(_test=M_Qh);
        //for(auto const & it : m_dirichlet_p)
        //    f2B += on(_range=markedfaces(M_Qh->mesh(),it.first),_element=phi, _expr=it.second, _rhs=f1B, _type=soption("blockms.22.on.type"));
        M_22Op = op(M_L, "blockms.22");

        if(soption("blockms.11.pc-type") == "AS")
        {
            M_pcAs->update(M_11, M_L, mu);
            M_11Op->setPc( M_pcAs );
        }
    }
    toc( "[PreconditionerBlockMS] update", FLAGS_v > 0 );
}

template < typename space_type, typename coef_space_type >
int
PreconditionerBlockMS<space_type,coef_space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
    tic();
    U = X;
    U.close();
    *M_uin = U.template element<0>();
    M_uin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    // Solve eq (12)
    if ( this->type() == AFP )
    {
        // solve here eq 15 : Pm v = c
        M_11Op->applyInverse(*M_uin,*M_uout);
        M_uout->close();
        NbIter1.push_back(M_11Op->solveReturn());

        // solve here eq 16
        M_22Op->applyInverse(*M_pin,*M_pout);
        M_pout->close();
        NbIter2.push_back(M_22Op->solveReturn());
    }
    else if( this->type() == SIMPLE )
    {
        // Nothing is done here
        *M_uout = *M_uin;
        M_uout->close();
        *M_pout = *M_pin;
        M_pout->close();
    }


    LOG(INFO) << "Update output potential/lagrange...\n";
    U.template element<0>() = *M_uout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();
    toc("[PreconditionerBlockMS] applyInverse update solution",FLAGS_v>0);
    return 0;
}

namespace meta
{
template< typename space_type , typename coef_space_type >
    struct blockms
{
    typedef PreconditionerBlockMS<space_type, coef_space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::blockms<
                                   typename parameter::value_type<Args, tag::space >::type::element_type,
                                   typename parameter::value_type<Args, tag::space2 >::type::element_type
                                   >::ptrtype ),
                                 blockms,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( space2, *)
                                   ( matrix, *)
                                 )
                                 ( optional
                                   ( type, *, soption("blockms.type"))
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                 )
                               )
{
    typedef typename meta::blockms<
        typename parameter::value_type<Args, tag::space>::type::element_type,
        typename parameter::value_type<Args, tag::space2>::type::element_type
                     >::ptrtype pblockms_t;

    typedef typename meta::blockms<
        typename parameter::value_type<Args, tag::space>::type::element_type,
        typename parameter::value_type<Args, tag::space2>::type::element_type
                     >::type blockms_t;
    pblockms_t p( new blockms_t( type, space, space2, bc, prefix, matrix ) );
    return p;
} // blockms
} // Feel
#endif
