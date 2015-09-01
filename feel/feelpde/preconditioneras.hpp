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


#include <feel/feeldiscr/operatorinterpolation.hpp>
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

    static const uint16_type Dim = space_type::nDim;
    typedef typename space_type::value_type value_type;

    enum Type
    {
        AS    = 0, // augmentation free preconditioner
        SIMPLE = 2
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
    // Qh3 - temporary version
    typedef Lagrange<1,Vectorial> lag_v_type;
    typedef FunctionSpace<mesh_type, bases< lag_v_type >> lag_v_space_type;
    typedef boost::shared_ptr<lag_v_space_type> lag_v_space_ptrtype;

    // Elements
    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;
    typedef typename coef_space_type::element_type element_coef_type;
    typedef typename lag_v_space_type::element_type lag_v_element_type;

    // operatrors
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
                      sparse_matrix_ptrtype Pm,
                      double k
                      );

    Type type() const { return M_type; }
    void setType( std::string t );

    //void update( sparse_matrix_ptrtype A, sparse_matrix_ptrtype L, element_coef_type mu );
    void update( sparse_matrix_ptrtype A, element_coef_type mu );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        this->applyInverse(X,Y);
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;

    virtual ~PreconditionerAS(){};

private:

    Type M_type;

    space_ptrtype M_Xh;
    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    lag_v_space_ptrtype M_Qh3;
    coef_space_ptrtype M_Mh;

    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;
    std::vector< std::vector<size_type> > M_Qh3_indices;

    sparse_matrix_ptrtype M_P;
    sparse_matrix_ptrtype M_Pt;
    sparse_matrix_ptrtype M_C;
    sparse_matrix_ptrtype M_Ct;

    op_ptrtype SimpleOp;
    op_mat_ptrtype opMat1;
    op_mat_ptrtype opMat2;
    op_inv_ptrtype opMatInv;

    mutable vector_ptrtype A; // A = inv(diagPm)r
    mutable vector_ptrtype B; // B = Py
    mutable vector_ptrtype C; // C = (1/g) Cz

    mutable vector_ptrtype M_r;
    mutable vector_ptrtype M_r_t;
    mutable vector_ptrtype M_uout;
    mutable vector_ptrtype M_diagPm;
    mutable vector_ptrtype M_t;
    mutable vector_ptrtype M_s;
    mutable vector_ptrtype M_y;
    mutable vector_ptrtype M_y_t;
    mutable vector_ptrtype M_z;
    mutable vector_ptrtype M_z_t;

    mutable potential_element_type U;
    element_coef_type M_mu;
    element_coef_type M_er;  // permittivity

    op_ptrtype M_lgqOp;   // \bar{L} + gamma \bar{Q}
    op_ptrtype M_lOp;   // L

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
                                                                sparse_matrix_ptrtype Pm,
                                                                double k )
    :
        M_type( AS ),
        M_Xh( Xh ),
        M_Vh(Xh->template functionSpace<0>() ),
        M_Qh(Xh->template functionSpace<1>() ),
        M_Mh( Mh ),
        M_Vh_indices( M_Vh->nLocalDofWithGhost() ),
        M_Qh_indices( M_Qh->nLocalDofWithGhost() ),
        M_Qh3_indices( Dim ),
        A(backend()->newVector(M_Vh)),
        B(backend()->newVector(M_Vh)),
        C(backend()->newVector(M_Vh)),
        M_r(backend()->newVector(M_Vh)),
        M_r_t(backend()->newVector(M_Vh)),
        M_uout(backend()->newVector(M_Vh)),
        M_diagPm(backend()->newVector(M_Vh)),
        //M_t(backend()->newVector(M_Vh)),
        U( M_Vh, "U" ),
        M_mu(M_Mh, "mu"),
        M_er(M_Mh, "er"),
        M_bcFlags( bcFlags ),
        M_prefix( p ),
        M_k(k),
        M_g(1.-k*k)
{
    tic();
    LOG(INFO) << "[PreconditionerAS] setup starts";
    this->setMatrix( Pm ); // Needed only if worldComm > 1

    // QH3 : Lagrange vectorial space type
    M_Qh3 = lag_v_space_type::New(Xh->mesh());

    // Block 11.1
    M_s = backend()->newVector(M_Qh3);
    M_y = backend()->newVector(M_Qh3);
    // Block 11.2
    M_z = backend()->newVector(M_Qh);
    M_t = backend()->newVector(M_Qh);
    
    // Create the interpolation and keep only the matrix
    auto pi_curl = I(_domainSpace=M_Qh3, _imageSpace=M_Vh);
    auto Igrad   = Grad( _domainSpace=M_Qh, _imageSpace=M_Vh);

    M_P = pi_curl.matPtr();
    M_C = Igrad.matPtr();
    
    M_Pt = backend()->newMatrix(M_Qh3,M_Vh);
    M_Ct = backend()->newMatrix(M_Qh3,M_Vh);

    M_P->transpose(M_Pt,MATRIX_TRANSPOSE_UNASSEMBLED);
    M_C->transpose(M_Ct,MATRIX_TRANSPOSE_UNASSEMBLED);

    LOG(INFO) << "size of M_C = " << M_C->size1() << ", " << M_C->size2() << std::endl;
    LOG(INFO) << "size of M_P = " << M_P->size1() << ", " << M_P->size2() << std::endl;

    // Create vector of indices to create subvectors/matrices
    std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 ); // Vh indices in Xh
    std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() ); // Qh indices in Xh

    // "Components" of Qh3
    auto Qh3_dof_begin = M_Qh3->dof()->dofPointBegin();
    auto Qh3_dof_end = M_Qh3->dof()->dofPointEnd();
    int idx=0;
    for(int i=0; i<M_Qh3_indices.size(); i++)
        M_Qh3_indices[i].resize( M_Qh3->compSpace()->nLocalDofWithGhost() );

    int dof_comp, dof_idx;
    for( auto it = Qh3_dof_begin; it!= Qh3_dof_end; it++ )
    {
        dof_comp = it->template get<2>(); //Component
        dof_idx = it->template get<1>(); //Global index
        M_Qh3_indices[dof_comp][idx] = dof_idx;
        if( dof_comp == Dim - 1 )
            idx++;
    }

    this->setType ( t );
    toc( "[PreconditionerAS] setup done ", FLAGS_v > 0 );
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
                                                      //sparse_matrix_ptrtype L,
                                                      element_coef_type mu )
{
    tic();
    M_er.on(_range=elements(M_Mh->mesh()), _expr=cst(1.));
    
    map_vector_field<FM_DIM,1,2> m_dirichlet_u { M_bcFlags.getVectorFields<FM_DIM> ( "u", "Dirichlet" ) };

    if(this->type() == AS)
    {
        // A = Pm
        // M_diagPm = compose(diag ( op (Pm, "blockms.11")),op(Pm,"blockms.11.diag")) ;
        backend()->diag(Pm,M_diagPm);

        /*
         * hat(L) = 1/mu L
         * bar(L) = diag( hat(L), hat(L), hat(L) )
         * bar(Q) = diag( er*Q, er*Q, er*Q ) with Q = mass matrix on Qh3
         * blockms.11.1 <=> bar(L) + g*bar(Q) y = s = Pt*r
         * blockms.11.2 <=> bar(L) z = t = trans(C)*r
         *
         */

        auto u = M_Qh->element("u");

        auto f11_1 = form2(M_Qh, M_Qh);
        auto rhs_11_1 = form1(M_Qh);
        f11_1 = integrate(_range=elements(M_Qh->mesh()),
                          _expr=1./idv(mu)*inner(grad(u),gradt(u))
                          + M_g*idv(M_er)*inner(id(u),idt(u)) );
        // TODO : boundary conditions ?
        for(auto const & it : m_dirichlet_u )
        {
            LOG(INFO) << "Applying 0 on " << it.first << " for blockms.11_1\n";
            f11_1 += on(_range=markedfaces(M_Qh->mesh(),it.first), _expr=cst(0.),_rhs=rhs_11_1, _element=u);
        }
        // Operator hat(L) + g Q
        M_lgqOp = op( f11_1.matrixPtr(), "blockms.11.1");

        auto f11_2 = form2(M_Qh, M_Qh); // hat(L)
        // TODO: do not rebuild the matrix - use f11_1
        f11_2 = integrate(_range=elements(M_Qh->mesh()),
                          _expr=1./idv(mu)*inner(grad(u),gradt(u)) );
        // TODO : boundary conditions ?
        for(auto const & it : m_dirichlet_u )
        {
            LOG(INFO) << "Applying 0 on " << it.first << " for blockms.11_2\n";
            f11_2 += on(_range=markedfaces(M_Qh->mesh(),it.first), _expr=cst(0.),_rhs=rhs_11_1, _element=u);
        }
        // Operator hat(L)
        M_lOp = op( f11_2.matrixPtr(), "blockms.11.2" );
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
     * With P_v^-1 = diag(P_m)^-1 (=A)
     *              + P (\bar L + g \bar Q) P^t (=B)
     *              + C (L^-1) C^T (=C)
     */

    tic();
    U = X;
    U.close();
    toc("Element created",FLAGS_v>0);

    // solve equ (12)
    if ( this->type() == AS )
    {
        tic();
        // RHS calculation
        *M_r = U;
        M_r->close();

        // étape A
        // diag(Pm)^-1*r
        //M_diagPm->pointwiseDivide(*M_r,*A);
        A->pointwiseDivide(*M_r,*M_diagPm);

        // s = P^t r
        M_Pt->multVector(M_r,M_s);

        auto y1 = M_y->createSubVector(M_Qh3_indices[0], true);
        auto y2 = M_y->createSubVector(M_Qh3_indices[1], true);
#if FM_DIM == 3
        auto y3 = M_y->createSubVector(M_Qh3_indices[2], true);
#endif
        auto s1 = M_s->createSubVector(M_Qh3_indices[0], true);
        auto s2 = M_s->createSubVector(M_Qh3_indices[1], true);
#if FM_DIM == 3
        auto s3 = M_s->createSubVector(M_Qh3_indices[2], true);
#endif
        /*
         * hat(L) + g Q is a (Qh,Qh) matrix
         * [[ hat(L) + g Q, 0  ,     0   ],    [ y1 ]    [ s1 ]
         * [   0,   hat(L) + g Q,    0   ], *  [ y2 ] =  [ s2 ]
         * [   0,     0   , hat(L) + g Q ]]    [ y3 ]    [ s3 ]
         */
        M_lgqOp->applyInverse(s1,y1);
        M_lgqOp->applyInverse(s2,y2);
#if FM_DIM == 3
        M_lgqOp->applyInverse(s3,y3);
#endif

        // y = [ y1, y2, y3 ]
        M_y->updateSubVector(y1, M_Qh3_indices[0]);
        M_y->updateSubVector(y2, M_Qh3_indices[1]);
#if FM_DIM == 3
        M_y->updateSubVector(y3, M_Qh3_indices[2]);
#endif
        M_y->close();
        // étape B 
        // P*y
        M_P->multVector(M_y,B);

        // t = C^t r
        M_Ct->multVector(M_r,M_t);

        // 14.b : hat(L) z = t
        M_lOp->applyInverse(M_t,M_z);
        M_z->close();

        // étape C 
        // M_C z
        M_C->multVector(M_z,C);
        C->scale(1./M_g);

        A->add(*C);
        A->add(*B);

        C->close();
        B->close();
        A->close();

        toc("15 assembled",FLAGS_v>0);
        *M_uout = *A; // 15 : w = A + B + C
    }
    else if( this->type() == SIMPLE )
    {
        SimpleOp->applyInverse(X, Y);
        *M_uout = Y;
    }
    else
    {
        Y=U;
        *M_uout = Y;
    }
    M_uout->close();

    tic();
    Y=*M_uout;
    Y.close();
    toc("PreconditionerAS::applyInverse", FLAGS_v>0 );

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
                                   ( matrix, *)
                                 )
                                 ( optional
                                   ( type, *, "AS")
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                   ( h, (double), doption("parameters.k") ) // for k ...
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
    pas_ptrtype p( new pas_type( type, space, space2, bc, prefix, matrix, h ) );
    return p;
#endif
} //
} // Feel
#endif
