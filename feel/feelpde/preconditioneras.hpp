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
#ifndef FEELPP_PRECONDITIONERAS_HPP
#define FEELPP_PRECONDITIONERAS_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/operatoras.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{
  template< typename space_type, typename mu_space_type >
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
    typedef boost::shared_ptr<mu_space_type> mu_space_ptrtype;
    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef typename mu_space_type::element_type element_mu_type;

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
     * \param A the full matrix 
     */
    PreconditionerAS( std::string t,
        space_ptrtype Xh,
        mu_space_ptrtype Mh, 
        BoundaryConditions bcFlags,
        std::string const& s,
        sparse_matrix_ptrtype A);

    Type type() const { return M_type; }
    void setType( std::string t );

    void initialize();

    //template< typename Expr_convection, typename Expr_bc >
    void update( sparse_matrix_ptrtype A, element_mu_type mu );

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

    backend_ptrtype M_cg1; // 14.a 
    backend_ptrtype M_cg2; // 14.b

    space_ptrtype M_Xh;
    mu_space_ptrtype M_Mh;

    sparse_matrix_ptrtype 
      M_pm, M_pm_inv,
      M_p, M_id, 
      M_l, 
      M_q, 
      M_c; // 10

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

    mutable element_type U;
    element_mu_type mu;

    op_ptrtype pm;

    BoundaryConditions M_bcFlags;
    std::string M_prefix;
    op_mat_ptrtype precHelm;

    value_type gamma;
  };

  template < typename space_type, typename mu_space_type >
    PreconditionerAS<space_type,mu_space_type>::PreconditionerAS( std::string t,
        space_ptrtype Xh, 
        mu_space_ptrtype Mh, 
        BoundaryConditions bcFlags,
        std::string const& p,
        sparse_matrix_ptrtype A )
    :
      M_type( AS ),
      M_cg1(backend()),
      M_Xh( Xh ),
      M_Mh( Mh ),
      U( M_Xh, "U" ),
      M_bcFlags( bcFlags ),
      M_prefix( p ),
      gamma(1.)
  {
    tic();
    LOG(INFO) << "[PreconditionerAS] setup starts";
    this->setMatrix( A );
    //std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
    //std::iota( M_Qh_indices.begin(), M_Qh_indices.end(), M_Vh->nLocalDofWithGhost() );

    this->createMatrices();

    initialize();

    this->setType ( t );
    toc( "[PreconditionerAS] setup done ", FLAGS_v > 0 );
  }

  template < typename space_type, typename mu_space_type >
    void
    PreconditionerAS<space_type,mu_space_type>::initialize()
    {
    }

  template < typename space_type, typename mu_space_type >
    void
    PreconditionerAS<space_type,mu_space_type>::createMatrices()
    {
#if 1
      tic();
      M_pm = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_pm_inv = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_p  = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_l  = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_q  = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_c  = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_id = M_cg1->newMatrix(M_Xh, M_Xh ); 
      M_r  = M_cg1->newVector(M_Xh); 
      M_r_t= M_cg1->newVector(M_Xh); 
      M_z  = M_cg1->newVector(M_Xh); 
      M_z_t= M_cg1->newVector(M_Xh); 
      M_y  = M_cg1->newVector(M_Xh); 
      M_y_t= M_cg1->newVector(M_Xh); 
      M_uout= M_cg1->newVector(M_Xh); 
      M_t  = M_cg1->newVector(M_Xh); 
      M_s  = M_cg1->newVector(M_Xh); 
      toc( "PreconditionerAS::createMatrix", FLAGS_v > 0 );
#endif

    }
  template < typename space_type, typename mu_space_type >
    void
    PreconditionerAS<space_type,mu_space_type>::setType( std::string t )
    {
      if ( t == "AS") M_type = AS;
      if ( t == "SIMPLE") M_type = SIMPLE;
#if 1
      LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
      switch( M_type )
      {
        case AS:
          tic();
          M_q->scale(gamma);
          M_q->addMatrix(1.,*M_l);
          opMat1 = op( M_q, "blockms.11.1");
          opMat2 = op( M_l, "blockms.11.2");
          this->setSide( super::LEFT );
          toc( "Preconditioner::setType AS", FLAGS_v > 0 );
          break;
        case SIMPLE:
          SimpleOp = op(M_id, "blockms.11.1");
          break;
      }
#endif
    }

  template < typename space_type, typename mu_space_type >
    //template< typename Expr_convection, typename Expr_bc >
    void
    PreconditionerAS<space_type,mu_space_type>::update( sparse_matrix_ptrtype A, element_mu_type mu )
    {
      // Mettre à jour les matrices
      tic();
      this->setMatrix( A );
      auto u = M_Xh->element("u");
      auto f21 = form2(M_Xh, M_Xh, _matrix=M_l);
      f21 = integrate(elements(M_Xh->mesh()), inner(id(u),idt(u)));
      M_l->close();
      M_p=M_l;
      M_c=M_l;
      M_pm_inv=M_l;
      opMat1 = op(M_l, "blockms.11.1");
      opMat2 = op(M_l, "blockms.11.2");
#if 0
      this->createMatrices();

      if ( type() == AS )
      {
        tic();
        afpOp->update( expr_b, g );
        toc( "Preconditioner::update AS", FLAGS_v > 0 );

      }
#endif
      toc( "PreconditionerAS::update", FLAGS_v > 0 );
    }



  template < typename space_type, typename mu_space_type >
    int
    PreconditionerAS<space_type,mu_space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
    {
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
        auto tmp = M_cg1->newVector(Y.mapPtr());
        *tmp = Y;
        tmp->close();
        M_p->multVector(tmp,M_s); // M_s = P^t r
        M_c->multVector(tmp,M_t); // M_t = C^t r
        opMat1->applyInverse(*M_s,*M_y);  // here solve 14.a
        toc("14.a solved");
        tic();
        opMat2->applyInverse(*M_t,*M_z);  // here solve 14.b
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
        //opMatInv->apply(M_pm,M_pm_inv); // Here compute diag(pm^-1)
        M_p->multVector(M_y,M_y_t);       // M_p*M_y -> M_y_t
        M_c->multVector(M_z,M_z_t);       // M_c*M_z -> M_z_t
        M_z_t->scale(1./gamma);           // 1/g*M_c*M_z -> M_z_t
        M_pm_inv->multVector(M_r, M_r_t); // M_pm_inv*M_r -> M_r_t
        M_r_t->add(*M_y_t);
        M_y_t->add(*M_z_t);
        toc("15 assembled");
        LOG(INFO) << "ocuouc\n"; *M_uout = *M_y_t; // 15;
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
      LOG(INFO) << "ocuouc\n"; Y=*M_uout;
      LOG(INFO) << "ocuouc\n"; Y.close();
      toc("PreconditionerAS::applyInverse update solution",FLAGS_v>0);
      toc("PreconditionerAS::applyInverse" );
#endif
      return 0;
    }

  template < typename space_type, typename mu_space_type >
    int
    PreconditionerAS<space_type,mu_space_type>::guess ( vector_type& Y ) const
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
    template< typename space_type , typename mu_space_type >
      struct blockas
      {
        typedef PreconditionerAS<space_type, mu_space_type> type;
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
    pas_ptrtype p( new pas_type( type, space, space2, bc, prefix, matrix ) );
    return p;
#endif
  } // 
} // Feel
#endif
