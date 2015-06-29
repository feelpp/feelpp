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
#ifndef FEELPP_PRECONDITIONERAS_HPP
#define FEELPP_PRECONDITIONERAS_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
//#include <feel/feelpde/operatoras.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{
template< typename space_type, typename coef_space_type >
    class PreconditionerAS : public Preconditioner<typename space_type::value_type>
{
    typedef Preconditioner<typename space_type::value_type> super;
public:

    enum Type
    {
        AS    = 0, // augmentation free preconditioner
        SIMPLE = 2 // 
    };
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
    // Vh
    typedef typename space_type::template sub_functionspace<0>::type potential_space_type;
    typedef typename space_type::template sub_functionspace<1>::type lagrange_space_type;
    // Qh
    typedef typename space_type::template sub_functionspace<0>::ptrtype potential_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype lagrange_space_ptrtype;

    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;
    typedef typename coef_space_type::element_type element_coef_type;

    typedef typename space_type::value_type value_type;

    static const uint16_type Dim = space_type::nDim;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef boost::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorBase<value_type> op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;

    /**
     * \param t Kind of prec
     * \param Xh potential space type
     * \param Mh Permeability space type
     * \param bcFlags the boundary conditions flags
     * \param s name of backend
     * \param k value
     */
    PreconditionerAS( std::string t,
                      space_ptrtype Xh,
                      coef_space_ptrtype Mh, 
                      BoundaryConditions bcFlags,
                      std::string const& s,
                      double k
                      );

    Type type() const { return M_type; }
    void setType( std::string t );

    void initialize();

    void update( sparse_matrix_ptrtype A, sparse_matrix_ptrtype L, element_coef_type mu );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        this->applyInverse(X,Y); 
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;
    int guess( vector_type& U ) const;

    virtual ~PreconditionerAS(){};

private:
    void createMatrices();
private:

    Type M_type;

    space_ptrtype M_Xh;
    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    coef_space_ptrtype M_Mh;

    sparse_matrix_ptrtype 
        M_p, M_p_t,
        M_c, M_c_t; // 10

    op_ptrtype SimpleOp;
    op_mat_ptrtype opMat1,
                   opMat2;
    op_inv_ptrtype opMatInv;

    mutable vector_ptrtype 
        M_r, M_r_t,
        M_uout,
        M_t, M_s,
        M_y, M_y_t,
        M_z, M_z_t;

    mutable potential_element_type U;
    element_coef_type M_mu;

    op_ptrtype M_diagPm, // diag(Pm)
               M_lgqOp,    // \bar{L} + gamma \bar{Q}
               M_lOp      // L
               ;

    BoundaryConditions M_bcFlags;
    std::string M_prefix;

    value_type M_k; // k
    value_type M_g; // gamma
};

template < typename space_type, typename coef_space_type >
PreconditionerAS<space_type,coef_space_type>::PreconditionerAS( std::string t,
                                                              space_ptrtype Xh, 
                                                              coef_space_ptrtype Mh, 
                                                              BoundaryConditions bcFlags,
                                                              std::string const& p,
                                                              double k)
    :
        M_type( AS ),
        M_Xh( Xh ),
        M_Vh(Xh->template functionSpace<0>() ),
        M_Qh(Xh->template functionSpace<1>() ),
        M_Mh( Mh ),
        M_p_t(backend()->newMatrix(M_Qh, M_Vh)),
        M_c_t(backend()->newMatrix(M_Qh, M_Vh)),
        M_t(backend()->newVector(M_Vh)),
        M_s(backend()->newVector(M_Vh)),
        U( M_Vh, "U" ),
        M_mu(M_Mh, "mu"),
        M_bcFlags( bcFlags ),
        M_prefix( p ),
        M_k(k),
        M_g(1.-k*k)
{
    tic();
    LOG(INFO) << "[PreconditionerAS] setup starts";
    //std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
    //std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

    this->createMatrices();

    initialize();

    this->setType ( t );
    toc( "[PreconditionerAS] setup done ", FLAGS_v > 0 );
}

template < typename space_type, typename coef_space_type >
    void
PreconditionerAS<space_type,coef_space_type>::initialize()
{
#if 0
    tic();
    toc( "PreconditionerAS::initialize", FLAGS_v > 0 );
#endif
}

template < typename space_type, typename coef_space_type >
    void
PreconditionerAS<space_type,coef_space_type>::createMatrices()
{
#if 0
    tic();
    toc( "PreconditionerAS::createMatrix", FLAGS_v > 0 );
#endif
}
template < typename space_type, typename coef_space_type >
    void
PreconditionerAS<space_type,coef_space_type>::setType( std::string t )
{
    if ( t == "AS")     M_type = AS;
    if ( t == "SIMPLE") M_type = SIMPLE;
}

template < typename space_type, typename coef_space_type >
//template< typename Expr_convection, typename Expr_bc >
    void
PreconditionerAS<space_type,coef_space_type>::update( sparse_matrix_ptrtype Pm, 
                                                    sparse_matrix_ptrtype L,
                                                    element_coef_type mu )
{
    // Mettre à jour les matrices
    tic();

    if(this->type() == AS){
    // A = Pm
    M_diagPm = diag ( op (Pm, "blockms.11.diag")) ;
    M_lOp    = op( L, "blockms.11.2");

    auto u = M_Qh->element("u");
    auto f21 = form2(M_Qh, M_Qh); // L + g Q - see 14.a
    f21 = integrate(_range=elements(M_Qh->mesh()),
                    _expr=1./idv(mu)*inner(grad(u),gradt(u)) // should be wrong
                    + M_g*inner(id(u),idt(u)));
    M_lgqOp  = op( f21.matrixPtr(), "blockms.11.1");
    //M_pm = f21.matrixPtr();
#if 0
    this->createMatrices();

    if ( type() == AS )
    {
        tic();
        afpOp->update( expr_b, g );
        toc( "Preconditioner::update AS", FLAGS_v > 0 );

    }
#endif
    }
    else if(this->type() == SIMPLE)
    {
        auto uu = M_Vh->element("uu");
        auto f22 = form2(M_Vh, M_Vh); 
        f22 = integrate(_range=elements(M_Vh->mesh()),
                        _expr=inner(id(uu),idt(uu)));
        SimpleOp = op( f22.matrixPtr(),"blockms.11.1");
    }
    toc( "PreconditionerAS::update", FLAGS_v > 0 );
}



template < typename space_type, typename coef_space_type >
int
PreconditionerAS<space_type,coef_space_type>::applyInverse ( const vector_type& X /*R*/, vector_type& Y /*W*/) const
{
    /*
     * We solve Here P_v w = r
     */
#if 1
    // Decompose les éléments
    tic();
    U = X;
    U.close();
    toc("Element created");

    // résout l'équation 12
    if ( this->type() == AS )
    {
        tic();
        // RHS calculation
        auto tmp = backend()->newVector(M_Vh);
        *tmp = Y;
        /*
         * pi_curl = opInterpolation(_domainspace=Qh3, _imagespace=Vh);
         *
         * Project R with P^t in s
         * M_lgqOp->applyInverse(s,y);
         * Project R with C^t in t
         * M_lOp->applyInverse(t,z);
         * M_diagOp->applyInverse(r,A)
         * Project y with P -> B
         * Project z with 1/g C -> C
         * Answer is A+B+C
         */
        tmp->close();
        std::cout 
            << M_p_t->size1() <<"\t"
            << M_p_t->size2() <<"\t"
            << tmp->size() <<"\t"
            << M_s->size() <<std::endl;
        M_p_t->multVector(tmp,M_s); // M_s = P^t r
        M_c_t->multVector(tmp,M_t); // M_t = C^t r
        M_lgqOp->applyInverse(*M_s,*M_y);  // here solve 14.a
        toc("14.a solved");
        tic();
        M_lOp->applyInverse(*M_t,*M_z);  // here solve 14.b
        toc("14.b solved");
        if(FLAGS_v>3)
        {
            M_p->print();
            M_c->print();
            tmp->print();
            M_s->print();
            M_t->print();
            M_y->print();
            M_z->print();
        }
        tic();
        M_diagPm->applyInverse(tmp,M_r_t);// Here compute diag(pm^-1)
        M_p->multVector(M_y,M_y_t);       // M_p*M_y -> M_y_t
        M_c->multVector(M_z,M_z_t);       // M_c*M_z -> M_z_t
        M_z_t->scale(1./M_g);             // 1/g*M_c*M_z -> M_z_t
        M_r_t->add(*M_y_t);
        M_y_t->add(*M_z_t);
        toc("15 assembled");
        *M_uout = *M_y_t; // 15;
    }
    else if( this->type() == SIMPLE ){
        SimpleOp->applyInverse(X, Y);
        *M_uout = Y;
    }
    else{
        Y=X;
        *M_uout = Y;
    }

    tic();
    LOG(INFO) << "coucou\n"; Y=*M_uout;
    LOG(INFO) << "coucou\n"; Y.close();
    toc("PreconditionerAS::applyInverse" );
#endif
    return 0;
}

template < typename space_type, typename coef_space_type >
int
PreconditionerAS<space_type,coef_space_type>::guess ( vector_type& Y ) const
{
#if 0
    U = Y;
    U.close();

    LOG(INFO) << "Create potential/lagrange component...\n";
    *M_ain = U.template element<0>();
    M_ain->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    LOG(INFO) << "lagrange/potential blockas : apply divergence...\n";
    divOp->apply( *M_pout, *M_ain );

    M_aux->zero();
    M_aux->add( -1.0, *M_ain );
    M_aux->close();

    LOG(INFO) << "potential blockas : apply inverse convection diffusion...\n";
    helmOp->applyInverse(*M_aux, *M_ain);
    LOG(INFO) << "Update output potential/lagrange...\n";

    U.template element<0>() = *M_ain;
    U.template element<1>() = *M_pin;
    U.close();
    Y=U;
    Y.close();
#endif
    return 0;
}
namespace meta
{
template< typename space_type , typename coef_space_type >
    struct blockas
{
    typedef PreconditionerAS<space_type, coef_space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::blockas<
                                   typename parameter::value_type<Args, tag::space>::type::element_type,
                                   typename parameter::value_type<Args, tag::space2>::type::element_type
                                   >::ptrtype ),
                                 blockas,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( space2, *)
                                 )
                                 ( optional
                                   ( type, *, "AS")
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                   ( h, (double), 0.0 ) // for k ...
                                 )
                               )
{
#if 1
    typedef typename meta::blockas<
        typename parameter::value_type<Args, tag::space>::type::element_type,
        typename parameter::value_type<Args, tag::space2>::type::element_type
                     >::ptrtype pas_ptrtype;

    typedef typename meta::blockas<
        typename parameter::value_type<Args, tag::space>::type::element_type,
        typename parameter::value_type<Args, tag::space2>::type::element_type
                     >::type pas_type;
    pas_ptrtype p( new pas_type( type, space, space2, bc, prefix, h ) );
    return p;
#endif
} // 
} // Feel
#endif
