/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-03

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
   \file oseen.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-03
 */

/* Oseen solver for stable finite element pairs,
   weak Dirichlet boundary conditions [1],
   building u-p-coupled matrix solved by bicgstab, preconditioned by ILUp

   equation solved:
   sigma*u + (betax, betay)*grad u - nu laplace u + grad p = f on the domain
                                                     div u = 0 on the domain
                                                         u = g on boundary

   References:
   [1] E. Burman, M. Fernandez and P. Hansbo: Continuous interior penalty finite
       element method for Oseen's equations, SIAM J. Numer. Anal. Vol. 44,
       No. 3, pp. 1248-1274
*/

#ifndef _OSEEN_HPP_
#define _OSEEN_HPP_

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>

#include <feel/feelvf/vf.hpp>
#include <gmm_iter_solvers.h>

namespace Feel
{
template<class Space_u,
         class Space_p,
         uint16_type imOrder = 3*Space_u::basis_type::nOrder-1>
class Oseen
{
public:

    // -- TYPEDEFS --
    typedef Space_u space_u_type;
    typedef Space_p space_p_type;

    static const uint16_type uOrder = space_u_type::basis_type::nOrder;
    static const uint16_type pOrder = space_p_type::basis_type::nOrder;

    static const uint16_type Dim = space_u_type::nDim;
    typedef typename space_u_type::value_type value_type;

    /* matrix */
    typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    /* mesh */
    typedef typename space_u_type::mesh_type mesh_type;
    typedef typename space_u_type::mesh_ptrtype mesh_ptrtype;

    /* space */
    typedef boost::shared_ptr<space_u_type> space_u_ptrtype;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef typename space_u_type::element_type element_u_type;
    typedef typename space_p_type::element_type element_p_type;

    /* quadrature */
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;

    Oseen( const space_u_ptrtype& spaceU,
           const space_p_ptrtype& spaceP )
        :
        M_mesh( spaceU->mesh() ),
        M_Xh( spaceU ),
        M_Yh( spaceP ),
        M_nTot( Dim * spaceU->nDof() + spaceP->nDof() ),
        ux( spaceU, "ux" ),
        uy( spaceU, "uy" ),
        p( spaceP, "p" ),
        M_M(),
        M_F( M_nTot ),
        M_up( M_nTot ),
        M_bcCoeff( 100 ),
        M_updated( false ),
        M_timers(),
        M_stats(),
        M_noisy( 0 ),
        M_maxiter( 1000 ),
        M_fillin( 2 ),
        M_threshold( 1.e-3 ),
        M_tol( 2.e-10 )
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
    void set_bccoeff( double bccoeff )
    {
        M_bcCoeff = bccoeff;
    }
    void set_tol( double tol )
    {
        M_tol = tol;
    }

    // update operator and rhs with given expressions
    template<typename Esigma, typename Enu,
             typename Ebetax, typename Ebetay,
             typename Efx, typename Efy,
             typename Egx, typename Egy>
    void update( Esigma sigma, Enu nu,
                 Ebetax betax, Ebetay betay,
                 Efx fx, Efy fy,
                 Egx gx, Egy gy );

    // solve system, call not needed, but possible
    void solve();

    // results
    const element_u_type& velocityX()
    {
        solve();
        return ux;
    }
    const element_u_type& velocityY()
    {
        solve();
        return uy;
    }
    const element_p_type& pressure()
    {
        solve();
        return p;
    }

private:

    /**
     * solve non-symmetric system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solveNonSym( Mat const& D, Vec1& u, Vec2 const& F );

private:

    mesh_ptrtype M_mesh;

    space_u_ptrtype M_Xh;
    space_p_ptrtype M_Yh;

    uint16_type M_nTot;

    element_u_type ux;
    element_u_type uy;
    element_p_type p;

    sparse_matrix_type M_M;
    std::vector<value_type> M_F;
    std::vector<value_type> M_up;

    double M_bcCoeff;

    bool M_updated;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;
    std::map<std::string,double> M_stats;

    int M_noisy;
    int M_maxiter;
    int M_fillin;
    double M_threshold;
    double M_tol;

}; // class Oseen



template<class Space_u, class Space_p, uint16_type imOrder>
template<typename Esigma, typename Enu,
         typename Ebetax, typename Ebetay,
         typename Efx, typename Efy,
         typename Egx, typename Egy>
void Oseen<Space_u, Space_p, imOrder>::update( Esigma sigma, Enu nu,
        Ebetax betax, Ebetay betay,
        Efx fx, Efy fy,
        Egx gx, Egy gy )
{
    M_updated = true;

    using namespace Feel::vf;

    M_timers["assembly"].first.restart();

    // --- Construction of velocity convection-diffusion-reaction matrix Lu
    sparse_matrix_type Lu;
    form( M_Xh, M_Xh, Lu ) =
        integrate( elements( *M_mesh ), im_type(),
                   ( nu )*( dxt( ux )*dx( ux ) + dyt( ux )*dy( ux ) )
                   + ( ( sigma )*idt( ux ) +
                       ( betax )*dxt( ux ) +
                       ( betay )*dyt( ux ) ) * id( ux )
                 );
    // use last arg 'false'
    // to tell the form to _not_ initialise its representation
    form( M_Xh, M_Xh, Lu, false ) +=
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   ( nu )*( - ( gradt( ux )*N() )*id ( ux )
                            - ( grad ( ux )*N() )*idt( ux )
                            + M_bcCoeff/hFace() * idt( ux ) * id( ux ) )
                 );

    // --- Construction of derivative matrices Dx and Dy
    sparse_matrix_type Dx;
    form( M_Xh, M_Yh, Dx ) =
        integrate( elements( *M_mesh ), im_type(),
                   -dx( ux )*idt( p )
                 );
    // use last arg 'false'
    // to tell the form to _not_ initialise its representation
    form( M_Xh, M_Yh, Dx, false ) +=
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   id( ux )*Nx()*idt( p )
                 );
    Dx.close();
    //     Dx.printMatlab("Dx.m");
    //     std::cout << "Dx: " << Dx.size1() << "x" << Dx.size2() << std::endl;

    sparse_matrix_type Dy;
    form( M_Xh, M_Yh, Dy ) =
        integrate( elements( *M_mesh ), im_type(),
                   -dy( uy )*idt( p )
                 );
    // use last arg 'false'
    // to tell the form to _not_ initialise its representation
    form( M_Xh, M_Yh, Dy, false ) +=
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   id( uy )*Ny()*idt( p )
                 );
    Dy.close();
    //     Dy.printMatlab("Dy.m");

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    VectorUblas<value_type> rhsx( ux.size() );
    VectorUblas<value_type> rhsy( uy.size() );

    form( M_Xh, rhsx ) = integrate( elements( *M_mesh ), im_type(),
                                    ( fx )*id( ux )
                                  );
    // use last arg 'false' to tell the form to
    // _not_ initialise its representation
    form( M_Xh, rhsx, false ) +=
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   ( nu )*( gx )*( M_bcCoeff/hFace()*id( ux ) - grad( ux )*N() )
                 );

    form( M_Xh, rhsy ) = integrate( elements( *M_mesh ), im_type(),
                                    ( fy )*id( uy )
                                  );
    // use last arg 'false' to tell the form to
    // _not_ initialise its representation
    form( M_Xh, rhsy, false ) +=
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   ( nu )*( gy )*( M_bcCoeff/hFace()*id( uy ) -grad( uy )*N() )
                 );

    // --- Construction of convection-diffusion-reaction
    //     operators on velocity space Auxx, Auxy and Auyy
    sparse_matrix_type Auxx( Lu );
    form( M_Xh, M_Xh, Auxx, false ) +=
        integrate( elements( *M_mesh ), im_type(),
                   ( nu )*dx( ux )*dxt( ux )
                 ) +
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   vf::sqrt( ( ( betax )^2.0 )+( ( betay )^2.0 ) )
                   *M_bcCoeff*idt( ux )*id( ux )*( Nx()^2.0 ) -
                   ( nu )*dxt( ux )*Nx()*id( ux )
                 );
    Auxx.close();

    sparse_matrix_type Auxy;
    form( M_Xh, M_Xh, Auxy ) =
        integrate( elements( *M_mesh ), im_type(),
                   ( nu )*dxt( uy )*dy( ux )
                 ) +
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   -( nu )*dxt( uy )*Ny()*id( ux )
                 );
    Auxy.close();

    sparse_matrix_type Auyy( Lu );
    // use last arg 'false' to tell the form to
    // _not_ initialise its representation
    form( M_Xh, M_Xh, Auyy, false ) +=
        integrate( elements( *M_mesh ), im_type(),
                   ( nu )*dy( uy )*dyt( uy )
                 ) +
        integrate( boundaryfaces( *M_mesh ), im_type(),
                   vf::sqrt( ( ( betax )^2.0 )+( ( betay )^2.0 ) )
                   *M_bcCoeff*idt( uy )*id( uy )*( Ny()^2.0 ) -
                   ( nu )*dyt( uy )*Ny()*id( uy )
                 );
    Auyy.close();

    M_timers["assembly"].second += M_timers["assembly"].first.elapsed();

    // --- Build full matrix
    M_M.init( M_nTot, M_nTot, M_nTot, M_nTot );
    gmm::sub_interval gmmRangeUx( 0, ux.size() );
    gmm::sub_interval gmmRangeUy( ux.size(), uy.size() );
    gmm::sub_interval gmmRangeP ( 2*ux.size(), p.size() );
    gmm::copy( Auxx.mat(),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUx, gmmRangeUx ) );
    gmm::copy( Auxy.mat(),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUx, gmmRangeUy ) );
    gmm::copy( gmm::transposed( Auxy.mat() ),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUy, gmmRangeUx ) );
    gmm::copy( Auyy.mat(),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUy, gmmRangeUy ) );
    gmm::copy( Dx.mat(),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUx, gmmRangeP  ) );
    gmm::copy( Dy.mat(),
               gmm::sub_matrix( M_M.wmat(), gmmRangeUy, gmmRangeP  ) );
    gmm::copy( gmm::transposed( gmm::scaled( Dx.mat(), -1.0 ) ),
               gmm::sub_matrix( M_M.wmat(), gmmRangeP,  gmmRangeUx ) );
    gmm::copy( gmm::transposed( gmm::scaled( Dy.mat(), -1.0 ) ),
               gmm::sub_matrix( M_M.wmat(), gmmRangeP,  gmmRangeUy ) );
    M_M.close();
    //     M_M.printMatlab("M_oseen.m");

    // --- Build rhs vector
    std::copy( rhsx.begin(), rhsx.end(), M_F.begin() );
    std::copy( rhsy.begin(), rhsy.end(), M_F.begin()+ux.size() );

} // update

template<class Space_u, class Space_p, uint16_type imOrder>
void Oseen<Space_u, Space_p, imOrder>::solve()
{
    // -- make sure solve is needed
    if ( !M_updated )
    {
        return;
    }

    M_updated = false;

    // -- solve
    solveNonSym( M_M, M_up, M_F );

    // -- extract solution
    std::copy( M_up.begin(), M_up.begin()+ux.size(), ux.begin() );
    std::copy( M_up.begin()+ux.size(),
               M_up.begin()+2*ux.size(), uy.begin() );
    std::copy( M_up.begin()+2*ux.size(),
               M_up.begin()+2*ux.size()+p.size(), p.begin() );

    VLOG(1) << "[timer] run():     init: " << M_timers["init"].second << "\n";
    VLOG(1) << "[timer] run(): assembly: " << M_timers["assembly"].second
            << "\n";

} // Oseen::solve


template<class Space_u, class Space_p, uint16_type imOrder>
template<typename Mat, typename Vec1, typename Vec2>
void
Oseen<Space_u, Space_p, imOrder>::solveNonSym( Mat const& D,
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
        std::cerr << "[Oseen] nonsymmetric linear solver didn't converge\n";
    }

    M_timers["solver"].second = M_timers["solver"].first.elapsed();
    VLOG(1) << "[timer] solveNonSym(): " << M_timers["solver"].second << "\n";
} // Oseen::solveNonSym

} // Feel

#endif /* _OSEEN_HPP_ */
