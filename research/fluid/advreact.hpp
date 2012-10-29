/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-04

  Copyright (C) 2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file advreact.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-04
 */
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>

#include <feel/feelvf/vf.hpp>
#include <gmm_iter_solvers.h>

namespace Feel
{

// A streamline diffusion stabilized advection-reaction solver.
// It solves the following equation:
// sigma * phi + (betax, betay)^T grad phi = f on the domain
//                                     phi = g on the inflow boundary
// the inflow boundary is detected through betax and betay

template<class Space, uint16_type imOrder = 2*Space::basis_type::nOrder>
class AdvReact
{
public:

    // -- TYPEDEFS --
    typedef Space space_type;

    static const uint16_type order = space_type::basis_type::nOrder;

    static const uint16_type Dim = space_type::nDim;
    typedef typename space_type::value_type value_type;

    /* matrix */
    typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    /* mesh */
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    /* space */
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* quadrature */
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;

    AdvReact( const space_ptrtype& space )
        :
        M_mesh( space->mesh() ),
        M_Xh( space ),
        M_nTot( space->nDof() ),
        M_phi( space, "phi" ),
        M_M(),
        M_F( M_nTot ),
        M_sol( M_nTot ),
        M_updated( false ),
        M_timers(),
        M_stats(),
        M_noisy( 0 ),
        M_maxiter( 1000 ),
        M_fillin( 2 ),
        M_threshold( 1.e-3 ),
        M_tol( 2.e-10 ),
        M_stabcoeff( 0.5 )
    {
        M_stats["nelt"] = M_mesh->elements().size();
        M_stats["ndof"] = M_nTot;
    }

    // setting of options
    void set_noisy( int noisy )
    {
        M_noisy = noisy;
    }
    void set_maxiter( int maxiter )
    {
        M_maxiter = maxiter;
    }
    void set_fillin( int fillin )
    {
        M_fillin = fillin;
    }
    void set_threshold( double threshold )
    {
        M_threshold = threshold;
    }
    void set_tol( double tol )
    {
        M_tol = tol;
    }
    void set_stabcoeff( double stabcoeff )
    {
        M_stabcoeff = stabcoeff;
    }

    // update operator and rhs with given expressions
    template<typename Esigma,
             typename Ebetax, typename Ebetay,
             typename Ef, typename Eg>
    void update( Esigma sigma,
                 Ebetax betax, Ebetay betay,
                 Ef f, Eg g );

    // solve system, call not needed, but possible
    void solve();

    // result
    const element_type phi()
    {
        solve();
        return M_phi;
    }

private:

    /**
     * solve non-symmetric system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solveNonSym( Mat const& D, Vec1& u, Vec2 const& F );

private:

    mesh_ptrtype M_mesh;

    space_ptrtype M_Xh;

    uint16_type M_nTot;

    element_type M_phi;

    sparse_matrix_type M_M;
    std::vector<value_type> M_F;
    std::vector<value_type> M_sol;

    bool M_updated;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;
    std::map<std::string,double> M_stats;

    int M_noisy;
    int M_maxiter;
    int M_fillin;
    double M_threshold;
    double M_tol;
    double M_stabcoeff;

}; // class AdvReact



template<class Space, uint16_type imOrder>
template<typename Esigma,
         typename Ebetax, typename Ebetay,
         typename Ef, typename Eg>
void AdvReact<Space, imOrder>::update( Esigma sigma,
                                       Ebetax betax, Ebetay betay,
                                       Ef f, Eg g )
{
    M_updated = true;

    using namespace Feel::vf;

    M_timers["assembly"].first.restart();

    value_type eps = type_traits<value_type>::epsilon();
    form( M_Xh, M_Xh, M_M ) =
        integrate( elements( *M_mesh ), im_type(),
                   ( ( sigma )*idt( M_phi )
                     + ( ( betax )*dxt( M_phi ) + ( betay )*dyt( M_phi ) ) )
                   * ( id( M_phi ) +
                       M_stabcoeff*h()/( sqrt( ( ( betax )^2.0 )+( ( betay )^2.0 ) )+eps ) *
                       ( ( betax )*dx( M_phi ) + ( betay )*dy( M_phi ) )
                     )
                 )
        + integrate( boundaryfaces( *M_mesh ), im_type(),
                     -chi( ( ( betax )*Nx()+( betay )*Ny() )<0 )* /* inflow */
                     ( ( betax )*Nx()+( betay )*Ny() )*idt( M_phi )*id( M_phi )
                   );
    M_M.close();
    //     M_M.printMatlab("M_advreact.m");

    VectorUblas<value_type> rhs( M_phi.size() );

    form( M_Xh, rhs ) =
        integrate( elements( *M_mesh ), im_type(),
                   ( f )
                   * ( id( M_phi ) +
                       M_stabcoeff*h()/( sqrt( ( ( betax )^2.0 )+( ( betay )^2.0 ) )+eps ) *
                       ( ( betax )*dx( M_phi ) + ( betay )*dy( M_phi ) )
                     )
                 )
        + integrate( boundaryfaces( *M_mesh ), im_type(),
                     -chi( ( ( betax )*Nx()+( betay )*Ny() )<0 )* /* inflow */
                     ( ( betax )*Nx()+( betay )*Ny() )*( g )*id( M_phi )
                   );

    std::copy( rhs.begin(), rhs.end(), M_F.begin() );

    M_timers["assembly"].second += M_timers["assembly"].first.elapsed();

} // update

template<class Space, uint16_type imOrder >
void AdvReact<Space, imOrder>::solve()
{
    // -- make sure solve is needed
    if ( !M_updated )
    {
        return;
    }

    M_updated = false;

    // -- solve
    solveNonSym( M_M, M_sol, M_F );

    // get result...
    std::copy( M_sol.begin(), M_sol.end(), M_phi.container().begin() );

    VLOG(1) << "[timer] run(): assembly: " << M_timers["assembly"].second
            << "\n";

} // AdvReact::solve


template<class Space, uint16_type imOrder>
template<typename Mat, typename Vec1, typename Vec2>
void
AdvReact<Space, imOrder>::solveNonSym( Mat const& D,
                                       Vec1& u,
                                       Vec2 const& F )
{
    M_timers["solver"].first.restart();

    gmm::iteration iter( M_tol );
    iter.set_noisy( M_noisy );
    iter.set_maxiter( M_maxiter );

    gmm::ilutp_precond<typename sparse_matrix_type::matrix_type>
    P( D.mat(), M_fillin, M_threshold );
    //     gmm::diagonal_precond<typename sparse_matrix_type::matrix_type>
    //         P( D.mat() );
    //     gmm::identity_matrix P;
    gmm::bicgstab( D.mat(), u, F, P, iter );

    if ( !iter.converged() )
    {
        std::cerr << "[AdvReact] nonsymmetric linear solver didn't converge\n";
    }

    M_timers["solver"].second += M_timers["solver"].first.elapsed();
    VLOG(1) << "[timer] solveNonSym(): " << M_timers["solver"].second << "\n";
} // AdvReact::solveNonSym

} // Feel
