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
#ifndef FEELPP_OPERATORAS_HPP
#define FEELPP_OPERATORAS_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelpde/preconditioneras.hpp>

namespace Feel
{

  template<typename space_type, typename mu_space_type>
    class OperatorAS : public OperatorBase<typename space_type::value_type>
  {
    typedef OperatorBase<typename space_type::value_type> super;
    public:

    typedef OperatorAS<space_type, mu_space_type> type;
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

    typedef typename mu_space_type::element_type element_mu_type;

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

    typedef PreconditionerAS<space_type, mu_space_type> pc_as_type;
    typedef boost::shared_ptr<pc_as_type> pc_as_ptrtype;


    static const uint16_type Dim = space_type::nDim;

    OperatorAS( space_ptrtype Qh,        // potential x lagrange
        mu_space_ptrtype Mh,     // permeability
        sparse_matrix_ptrtype A, // the full matrix
        backend_ptrtype b1,      // to solve 14.a
        BoundaryConditions const& bcFlags, // what for ?
        std::string const& p     // the prefix
        );

    OperatorAS( const OperatorAS& tc ) = default;
    OperatorAS( OperatorAS&& tc ) = default;
    OperatorAS& operator=( const OperatorAS& tc ) = default;
    OperatorAS& operator=( OperatorAS&& tc ) = default;

    void initialize();

    void update( sparse_matrix_ptrtype, element_mu_type );

    void setProblemType( std::string prob_type )
    {
      M_prob_type = prob_type;
    }

    std::string problemType() const
    {
      return M_prob_type;
    }

    ~OperatorAS() {};

    int apply(const vector_type& X, vector_type& Y) const;
    int applyInverse(const vector_type& X, vector_type& Y) const;

    private:
    backend_ptrtype M_backend;

    space_ptrtype M_Xh;
    mu_space_ptrtype M_Mh;

    element_mu_type M_mu;

    sparse_matrix_ptrtype M_pm,
                          M_p, 
                          M_l, 
                          M_q,
                          M_c, 
                          M_A;
    vector_ptrtype rhs;

    op_mat_ptrtype diapOp, pOp, lOp, qOp, cOp;

    BoundaryConditions M_bcFlags;
    std::string M_prefix;

    op_mat_ptrtype precOp1; 
    op_mat_ptrtype precOp2; 

    pc_as_ptrtype M_pcAs;

    std::string M_prob_type;

    void applyBC( sparse_matrix_ptrtype& A );
  };




  template < typename space_type, typename mu_space_type>
    OperatorAS<space_type, mu_space_type>::OperatorAS( 
        space_ptrtype Qh,
        mu_space_ptrtype Mh,
        sparse_matrix_ptrtype A,
        backend_ptrtype b1,
        BoundaryConditions const& bcFlags,
        std::string const& p)
    :
      super( Qh->mapPtr(), "blockms.11", false, false ),
      M_backend( b1 ),
      M_Xh( Qh ),
      M_Mh( Mh ),
      M_mu( M_Mh, "m" ),
      M_A( A ),
      rhs(         M_backend->newVector( M_Xh ) ),
      M_bcFlags( bcFlags ),
      M_prefix( p )
  {
    initialize();
    M_pcAs = boost::make_shared<pc_as_type>(soption("blockms.11.pc-type"), M_Xh, M_Mh, bcFlags, p, A);
    this->setPc(M_pcAs);
#if 0
    this->assembleMass();
    this->assembleDiffusion();
#endif
  }

  template < typename space_type, typename mu_space_type>
    void
    OperatorAS<space_type, mu_space_type>::initialize()
    {
      rhs->zero();
      rhs->close();
    }

  template < typename space_type, typename mu_space_type>
    void
    OperatorAS<space_type, mu_space_type>::update( sparse_matrix_ptrtype A, element_mu_type mu )
    {
      tic();
      M_A = A;
      M_mu.on(_range=elements(M_Mh->mesh()), _expr=idv(mu));
      M_pcAs->update( M_A, M_mu);
#if 0
      auto conv  = form2( _test=M_Mh, _trial=M_Mh, _matrix=G );
      G->zero();

      conv = integrate( _range=elements(M_Mh->mesh()), _expr=(trans(expr_b)*trans(gradt(p)))*id(q));
      conv += integrate( _range=elements(M_Mh->mesh()), _expr=idv(mu)*gradt(p)*trans(grad(q)));

      if ( soption("blockns.pcd.inflow") == "Robin" )
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
          LOG(INFO) << "Setting Robin condition on " << dir.marker();
          if ( ebc.find( dir.marker() ) != ebc.end() ){
            //conv += integrate( _range=markedfaces(M_Mh->mesh(), dir.marker()), _expr=-M_rho*trans(ebc.find(dir.marker())->second)*N()*idt(p)*id(q));
          }
        }

      G->close();

      this->applyBC(G);

      static bool init_G = false;

      //if ( !init_G )
      {
        // S = F G^-1 M
        LOG(INFO) << "[OperatorAS] setting pcd operator...\n";
        if ( ioption("blockns.pcd.order") == 1 )
          precOp = compose( massOp, compose(inv(op(G,"Fp")),diffOp) );
        else
          precOp = compose( diffOp, compose(inv(op(G,"Fp")),massOp) );
        LOG(INFO) << "[OperatorAS] setting pcd operator done.\n";
        init_G = true;
      }
#endif
      toc("Operator::AS update",FLAGS_v>0);
    }

#if 0



  template < typename space_type, typename mu_space_type >
    void
    OperatorAS<space_type,mu_space_type>::assembleMass()
    {
      tic();
      auto m = form2( _test=M_Mh, _trial=M_Mh, _matrix=M_mass );
      m = integrate( elements(M_Mh->mesh()), idt(p)*id(q) );
      M_mass->close();
      massOp = op( M_mass, "Mp" );
      toc("OperatorAS::mass assembly",FLAGS_v>0);
    }

  template < typename space_type, typename mu_space_type >
    void
    OperatorAS<space_type,mu_space_type>::assembleDiffusion()
    {
      tic();
      if ( soption("blockns.pcd.diffusion") == "Laplacian" )
      {
        auto d = form2( _test=M_Mh, _trial=M_Mh, _matrix=M_diff );
        d = integrate( _range=elements(M_Mh->mesh()), _expr=gradt(p)*trans(grad(q)));

        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
          LOG(INFO) << "Diffusion Setting Dirichlet condition on lagrange on " << cond.marker();
          if ( boption("blockns.weakdir" ) )
            d+= integrate( markedfaces(M_Mh->mesh(),cond.marker()), _expr=-gradt(p)*N()*id(p)-grad(p)*N()*idt(p)+doption("penaldir")*idt(p)*id(p)/hFace() );
          else
            d += on( markedfaces(M_Mh->mesh(),cond.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
        //this->applyBC(M_diff);
      }

      diffOp = op( M_diff, "Ap" );
      toc("OperatorAS::diffusion assembly",FLAGS_v>0);
    }
#endif

  template < typename space_type, typename mu_space_type >
    void
    OperatorAS<space_type,mu_space_type>::applyBC( sparse_matrix_ptrtype& A )
    {
#if 0
      tic();
      auto a = form2( _test=M_Mh, _trial=M_Mh, _matrix=A );

      if ( soption("blockns.pcd.inflow") != "Robin" )
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
          a += on( markedfaces(M_Mh->mesh(),dir.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }

      // on neumann boundary on potential, apply Dirichlet condition on lagrange
      if ( soption("blockns.pcd.outflow") == "Dirichlet" )
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
          if ( boption("blockns.weakdir" ) )
            a+= integrate( markedfaces(M_Mh->mesh(),cond.marker()), _expr=-idv(mu)*gradt(p)*N()*id(p)-idv(mu)*grad(p)*N()*idt(p)+doption("penaldir")*idt(p)*id(p)/hFace() );
          else
            a += on( markedfaces(M_Mh->mesh(),cond.marker()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
      rhs->close();
      A->close();
      toc("OperatorAS::BC apply",FLAGS_v>0);
#endif
    }

  template < typename space_type, typename mu_space_type >
    int
    OperatorAS<space_type,mu_space_type>::apply(const vector_type& X, vector_type& Y) const
    {
      LOG(INFO) << "OperatorAS::apply\n";
      return applyInverse( X, Y );
    }
  template < typename space_type, typename mu_space_type >
    int
    OperatorAS<space_type,mu_space_type>::applyInverse(const vector_type& X, vector_type& Y) const
    {
      LOG(INFO) << "OperatorAS::applyInverse\n";
      M_pcAs->update(M_A, M_mu);
      return M_pcAs->applyInverse(X,Y);
    }


} // Feel

#endif
