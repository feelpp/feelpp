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
#include <feel/feelpde/operatoras.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelpde/preconditioneras.hpp>
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
    typedef typename space_type::template sub_functionspace<0>::ptrtype potential_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype lagrange_space_ptrtype;

    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;
    typedef typename coef_space_type::element_type element_coef_type;

    typedef typename space_type::value_type value_type;

    static const uint16_type Dim = space_type::nDim;

    typedef OperatorBase<value_type> op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef typename OperatorAS<potential_space_type,coef_space_type>::type op_as_type;
    typedef typename OperatorAS<potential_space_type,coef_space_type>::ptrtype op_as_ptrtype;

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

    void initialize();

    void update( sparse_matrix_ptrtype A, element_coef_type mu );

    void apply( const vector_type & X, vector_type & Y ) const
    {
      this->applyInverse(X,Y); 
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;
    int guess( vector_type& U ) const;

    virtual ~PreconditionerBlockMS(){};

    private:
    void createSubMatrices();

    Type M_type;

    backend_ptrtype M_backend;
    space_ptrtype M_Xh;
    coef_space_ptrtype M_Mh;

    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;

    mutable vector_ptrtype 
      M_uin,
      M_uout, 
      M_pin, 
      M_pout;

    mutable element_type U;

    sparse_matrix_ptrtype M_11;
    element_coef_type M_mu, // permeability
                      M_er;  // permittivity

    op_as_ptrtype  M_asOp; // Augmented Spaces
    op_ptrtype M_22Op; // 
    op_ptrtype M_11Op;     // if not augmented spaces

    value_type M_k; // wave number

    BoundaryConditions M_bcFlags;
    std::string M_prefix;

    potential_element_type u;
    lagrange_element_type phi;
  };

  template < typename space_type, typename coef_space_type >
    PreconditionerBlockMS<space_type,coef_space_type>::PreconditionerBlockMS( 
        std::string t,                // Type
        space_ptrtype Xh,             // (u)x(p)
        coef_space_ptrtype Mh,        // mu
        BoundaryConditions bcFlags,   // bc
        std::string const& p,         // prefix
        sparse_matrix_ptrtype A )     // The matrix
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
      M_11(M_backend->newMatrix(M_Vh,M_Vh)),
      M_mu( M_Mh, "mu" ),
      M_er( M_Mh, "er" ),
      M_k(0.),
      M_bcFlags( bcFlags ),
      M_prefix( p ),
      u(M_Vh, "u"),
      phi(M_Qh, "phi")
  {
    tic();
    LOG(INFO) << "[PreconditionerBlockMS] setup starts";
    this->setMatrix( A );
    std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
    std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

    this->createSubMatrices();

    initialize();

    this->setType ( t );
    toc( "[PreconditionerBlockMS] setup done ", FLAGS_v > 0 );
  }

  template < typename space_type, typename coef_space_type >
    void
    PreconditionerBlockMS<space_type,coef_space_type>::initialize()
    {
    }

  template < typename space_type, typename coef_space_type >
    void
    PreconditionerBlockMS<space_type,coef_space_type>::createSubMatrices()
    {
#if 0
      tic();
      toc( "PreconditionerBlockMS::createSubMatrix(M_22,M_11)", FLAGS_v > 0 );
#endif

    }
  template < typename space_type, typename coef_space_type >
    void
    PreconditionerBlockMS<space_type,coef_space_type>::setType( std::string t )
    {
      if ( t == "AFP") M_type = AFP;
      if ( t == "SIMPLE") M_type = SIMPLE;
      LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
    }

  template < typename space_type, typename coef_space_type >
    //template< typename Expr_convection, typename Expr_bc >
    void
    PreconditionerBlockMS<space_type,coef_space_type>::update( sparse_matrix_ptrtype A, element_coef_type mu )
    {
      tic();
      this->setMatrix( A );
      M_mu.on(_range=elements(M_Mh->mesh()), _expr=idv(mu));;
      M_er.on(_range=elements(M_Mh->mesh()), _expr=cst(1.));;

      LOG(INFO) << "Create sub Matrix\n";
      // calculer matrice A + g M
      auto f2A = form2(_test=M_Vh, _trial=M_Vh,_matrix=M_11);
      auto f1A = form1(_test=M_Vh);
      f2A = integrate(_range=elements(M_Vh->mesh()), _expr=cst(1.)/idv(M_mu)*trans(curlt_op(u))*curl_op(u) // mu^-1 A
          +cst(1.-M_k*M_k)*inner(idt(u),id(u))); // g M
      f2A += on(_range=boundaryfaces(M_Vh->mesh()), _expr=zero<FM_DIM,1>(),_rhs=f1A, _element=u);

      if(soption("blockms.11.pc-type") == "AS"){
        // donner à M_asOp->update(M_asOp,M_mu);
        M_asOp = boost::make_shared<op_as_type>( M_Vh, M_Mh, M_11, M_backend, M_bcFlags, M_prefix );
        LOG(INFO) << "M_asOp->update()\n";
        M_asOp->update(M_11, M_mu);
      }
      else{
        LOG(INFO) << "M_11Op->update()\n";
        M_11Op = op(M_11, "blockms.11");
      }

      // calcule matrice L
      auto f2B = form2(_trial=M_Qh, _test=M_Qh);
      auto f1B = form1(_test=M_Qh);
      f2B = integrate(_range=elements(M_Qh->mesh()), _expr=idv(M_er)*inner(gradt(phi), grad(phi)));
      f2B += on(_range=boundaryfaces(M_Qh->mesh()),_element=phi, _expr=cst(0.), _rhs=f1B); // rajouter option elimination_keep-diag
      M_22Op = op(f2B.matrixPtr(), "blockms.22");

      toc( "Preconditioner::update", FLAGS_v > 0 );
    }

  template < typename space_type, typename coef_space_type >
    int
    PreconditionerBlockMS<space_type,coef_space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
    {
      // Mettre à jour les opérateurs M_asOp et M_22Op
      // Decompose les éléments
      tic();
      U = X;
      U.close();
      *M_uin = U.template element<0>();
      M_uin->close();
      *M_pin = U.template element<1>();
      M_pin->close();

      // résout l'équation 12
      if ( this->type() == AFP )
      {
        tic();
        // solve here eq 15 : Pm v = c
        // We can use the AS preconditioner or one given thanks to PETSc
        if(soption("blockms.11.pc-type") == "AS"){ // fictious space ?
          M_asOp->applyInverse(*M_uin,*M_uout);
          M_uout->close();
        }
        else{
          M_11Op->applyInverse(*M_uin,*M_uout);
          M_uout->close();
        }
        toc("blockms.11 solved",FLAGS_v>0);

        tic();
        // solve here eq 16
        M_22Op->applyInverse(*M_pin,*M_pout);
        M_pout->close();
        toc("blockms.Lag solved",FLAGS_v>0);
      }
      else if( this->type() == SIMPLE )
      {
        tic();
        // Nothing is done here
        *M_uout = *M_uin;
        M_uout->close();
        *M_pout = *M_pin;
        M_pout->close();
        toc("Dummy solved",FLAGS_v>0);
      }


      LOG(INFO) << "Update output potential/lagrange...\n";
      tic();
      U.template element<0>() = *M_uout;
      U.template element<1>() = *M_pout;
      U.close();
      Y=U;
      Y.close();
      toc("PreconditionerBlockMS::applyInverse update solution",FLAGS_v>0);
      return 0;
    }

  template < typename space_type, typename coef_space_type >
    int
    PreconditionerBlockMS<space_type,coef_space_type>::guess ( vector_type& Y ) const
    {
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
  } // btcpd
} // Feel
#endif
