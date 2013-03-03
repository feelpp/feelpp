/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-03

  Copyright (C) 2006 EPFL
  Copyright (C) 2011 Université Joseph Fourier

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

/* Oseen solver with interior penalty stabilization, and
   weak Dirichlet and homogeneous Neumann boundary conditions [1],
   building u-p-coupled matrix solved by an adaptive gmm backend

   equation solved:
   sigma u + (grad u)*beta - div(2 nu sgrad u) + grad p = f on domain
                                      div u + epsilon p = c on domain
                                                      u = g on D-bdry
                                   2 nu sgrad u n - p n = 0 on N-bdry

   where sgrad is the symmetric gradient tensor.

   References:
   [1] E. Burman, M. Fernandez and P. Hansbo:
       Continuous interior penalty finite element method for Oseen's equations
       SIAM J. Numer. Anal. Vol. 44, No. 3, pp. 1248-1274
*/
#include <set>

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/oseendata.hpp>

namespace Feel
{


template<class Space,
         uint16_type imOrder = 3*boost::fusion::result_of::value_at_c<typename Space::basis_type, 0>::basis_type::nOrder-1,
         template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
class Oseen
{
public:

    // -- TYPEDEFS --
    typedef Space space_type;

    static const uint16_type Dim = space_type::nDim;
    typedef typename space_type::value_type value_type;

    /* backend */
    typedef Backend<value_type> backend_type;


    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /* mesh */
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    /* P0 basis */
    typedef bases<Lagrange<0,Scalar,Discontinuous> > basis_i_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_i_type, Discontinuous,value_type > space_i_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_u_type;
    typedef typename element_type::template sub_element<1>::type element_p_type;
    typedef typename space_i_type::element_type element_i_type;

    /* quadrature */
    static const uint16_type uOrder = boost::fusion::result_of::value_at_c<typename Space::basis_type, 0>::type::nOrder;
    static const uint16_type pOrder = boost::fusion::result_of::value_at_c<typename Space::basis_type, 1>::type::nOrder;
    template<int Order>
    struct im
    {
        typedef IM<Dim, Order, value_type, Entity> type;
    };

    Oseen( const space_ptrtype& space,
           const backend_ptrtype& backend,
           const std::set<flag_type>& dirichletFlags,
           const std::set<flag_type>& neumannFlags,
           bool weak_dirichlet = true,
           bool export_matlab = false
         );

    Oseen( const space_ptrtype& space,
           const backend_ptrtype& backend,
           const std::set<flag_type>& dirichletFlags,
           const std::set<flag_type>& neumannFlags,
           po::variables_map const& vm,
           std::string const& prefix = ""
         );

    // setting of options
    void setBCCoeffDiff( value_type bccoeffdiff )
    {
        M_bcCoeffDiff = bccoeffdiff;
    }

    void setBCCoeffConv( value_type bccoeffconv )
    {
        M_bcCoeffConv = bccoeffconv;
    }

    void setStabCoeffDiv( value_type stabcoeffdiv )
    {
        M_stabCoeffDiv = stabcoeffdiv;
    }

    void setStabCoeffP( value_type stabcoeffp )
    {
        M_stabCoeffP = stabcoeffp;
    }

    void decouplePstab( element_p_type& pStab, value_type theta )
    {
        M_pStab = &pStab;
        M_thetaPstab = theta;
    }

    void couplePstab()
    {
        M_pStab = 0;
        M_thetaPstab = -1.0;
    }

    void setEpsCompress( value_type epscompress )
    {
        M_epsCompress = epscompress;
    }

    void setDivDivCoeff( value_type divdivcoeff );

    // update operator and rhs with given expressions
    template<typename ItRange, typename EsigmaInc,
             typename EnuInc, typename EnuAbs,
             typename Ebeta, typename Ef, typename Ec, typename Eg,
             typename EnoSlip>
    void update( const ItRange& itRange,
                 const EsigmaInc& sigmaInc,
                 const EnuInc& nuInc,
                 const EnuAbs& nuAbs,
                 const Ebeta& beta,
                 const Ef& f,
                 const Ec& c,
                 const Eg& g,
                 const EnoSlip& noSlip,
                 bool updateStabilization = true );

    // solve system, call not needed, but possible
    void solve();

    // results
    const element_u_type& velocity()
    {
        solve();
        return u;
    }
    const element_p_type& pressure()
    {
        solve();
        return p;
    }

    element_type& solution()
    {
        return M_sol;
    }
    const element_type& solution() const
    {
        return M_sol;
    }

    value_type stabilizationEnergy() const;

private:

    /**
     * solve non-symmetric system
     */
    void solveNonSym( sparse_matrix_ptrtype const& D, element_type& u, vector_ptrtype const& F );

private:

    mesh_ptrtype M_mesh;

    space_ptrtype M_space;

    element_type M_sol;
    element_u_type u;
    element_p_type p;

    sparse_matrix_ptrtype M_matrixFull;
    sparse_matrix_ptrtype M_matrixAu;
    sparse_matrix_ptrtype M_matrixStab;
    vector_ptrtype M_vectorRhsFull;
    vector_ptrtype M_vectorSolFull;

    value_type M_bcCoeffDiff;
    value_type M_bcCoeffConv;

    value_type M_stabCoeffDiv;
    value_type M_stabCoeffP;
    element_p_type* M_pStab;
    value_type M_thetaPstab;

    value_type M_epsCompress;
    value_type M_divDivCoeff;

    bool M_updated;

    backend_ptrtype M_backend;

    std::set<flag_type> M_dirichletFlags;
    std::set<flag_type> M_neumannFlags;

    space_i_ptrtype M_space_i;
    element_i_type  M_stabWeightP;

    bool M_weak_dirichlet;
    bool M_export_matlab;

}; // class Oseen

template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
Oseen<Space, imOrder, Entity>::Oseen( const space_ptrtype& space,
                                      const backend_ptrtype& backend,
                                      const std::set<flag_type>& dirichletFlags,
                                      const std::set<flag_type>& neumannFlags,
                                      bool weak_dirichlet,
                                      bool export_matlab )
    :
    M_mesh( space->mesh() ),
    M_space( space ),
    M_sol( space, "sol" ),
    u( M_sol.template element<0>() ),
    p( M_sol.template element<1>() ),
    M_matrixFull( backend->newMatrix( M_space, M_space ) ),
    M_matrixAu( backend->newMatrix( M_space, M_space ) ),
    M_matrixStab( backend->newMatrix( M_space, M_space ) ),
    M_vectorRhsFull( backend->newVector( M_space ) ),
    M_vectorSolFull( backend->newVector( M_space ) ),
    M_bcCoeffDiff( 100 ),
    M_bcCoeffConv( 100 ),
    M_stabCoeffDiv( 0.0 ),
    M_stabCoeffP( 0.0 ),
    M_pStab( 0 ),
    M_thetaPstab( -1.0 ),
    M_epsCompress( 0.0 ),
    M_divDivCoeff( 0.0 ),
    M_updated( false ),
    M_backend( backend ),
    M_dirichletFlags( dirichletFlags ),
    M_neumannFlags( neumannFlags ),
    M_space_i( new space_i_type( space->mesh() ) ),
    M_stabWeightP( M_space_i, "cs" ),
    M_weak_dirichlet( weak_dirichlet ),
    M_export_matlab( export_matlab )
{
    using namespace Feel::vf;

    element_u_type& v = u;
    element_p_type& q = p;

    DVLOG(2) << "[Oseen::Oseen] -(p,div v) + (q,div u)\n";
    form2( _test=M_space, _trial=M_space, _matrix=M_matrixAu, _init=true ) =
        integrate( elements( M_mesh ),
                   - div( v ) * idt( p )
                 )+
        integrate( elements( M_mesh ),
                   + id( q )  * divt( u )
                 );


    //M_matrixStab.close();

    // --- Set rhs to zero (pressure component is not set otherwise)
    M_vectorRhsFull->zero();

} // Oseen constructor

template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
Oseen<Space, imOrder, Entity>::Oseen( const space_ptrtype& space,
                                      const backend_ptrtype& backend,
                                      const std::set<flag_type>& dirichletFlags,
                                      const std::set<flag_type>& neumannFlags,
                                      po::variables_map const& vm,
                                      std::string const& prefix )
    :
    M_mesh( space->mesh() ),
    M_space( space ),
    M_sol( space, "sol" ),
    u( M_sol.template element<0>() ),
    p( M_sol.template element<1>() ),
    M_matrixFull( backend->newMatrix( M_space, M_space ) ),
    M_matrixAu( backend->newMatrix( M_space, M_space ) ),
    M_matrixStab( backend->newMatrix( M_space, M_space ) ),
    M_vectorRhsFull( backend->newVector( M_space ) ),
    M_vectorSolFull( backend->newVector( M_space ) ),
    M_bcCoeffDiff( 0.0 ),
    M_bcCoeffConv( 0.0 ),
    M_stabCoeffDiv( 0.0 ),
    M_stabCoeffP( 0.0 ),
    M_epsCompress( 0.0 ),
    M_divDivCoeff( 0.0 ),
    M_updated( false ),
    M_backend( backend ),
    M_dirichletFlags( dirichletFlags ),
    M_neumannFlags( neumannFlags ),
    M_space_i( new space_i_type( space->mesh() ) ),
    M_stabWeightP( M_space_i, "cs" ),
    M_weak_dirichlet( true ),
    M_export_matlab( true )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    M_bcCoeffDiff = vm[_prefix+"oseen-bc-coeff-diff"].template as<double>();
    M_bcCoeffConv = vm[_prefix+"oseen-bc-coeff-conv"].template as<double>();
    M_stabCoeffDiv = vm[_prefix+"oseen-stab-coeff-div"].template as<double>();
    M_stabCoeffP = vm[_prefix+"oseen-stab-coeff-p"].template as<double>();
    M_epsCompress = vm[_prefix+"oseen-eps-compress"].template as<double>();
    M_divDivCoeff = vm[_prefix+"oseen-divdiv-coeff"].template as<double>();
    M_weak_dirichlet = vm[_prefix+"oseen-weak-dirichlet"].template as<bool>();
    M_export_matlab = vm[_prefix+"oseen-export-matlab"].template as<bool>();

    using namespace Feel::vf;

    element_u_type& v = u;
    element_p_type& q = p;

    DVLOG(2) << "[Oseen::Oseen] -(p,div v) + (q,div u)\n";
    form2( _test=M_space, _trial=M_space, _matrix=M_matrixAu, _init=true ) =
        integrate( elements( M_mesh ),
                   - div( v ) * idt( p )
                 )+
        integrate( elements( M_mesh ),
                   + id( q )  * divt( u )
                 );

    //M_matrixStab.close();

    // --- Set rhs to zero (pressure component is not set otherwise)
    M_vectorRhsFull->zero();

} // Oseen constructor

template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void
Oseen<Space, imOrder, Entity>::setDivDivCoeff( value_type divdivcoeff )
{
    using namespace Feel::vf;

    element_u_type& v = u;

    value_type delta = divdivcoeff - M_divDivCoeff;

    if ( delta != 0.0 )
    {
        DVLOG(2) << "[Oseen::set_divdivcoeff] (h gamma div u,div v)\n";
        form2( _test=M_space, _trial=M_space, _matrix=M_matrixAu ) +=
            integrate( elements( M_mesh ),
                       h()*delta*div( v )*divt( u )
                     );
        M_divDivCoeff = divdivcoeff;
    }
}
template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
typename Oseen<Space, imOrder, Entity>::value_type
Oseen<Space, imOrder, Entity>::stabilizationEnergy() const
{
    vector_ptrtype temp( M_backend->newVector( M_space ) );
    //backend_type::applyMatrix( M_matrixStab, M_sol.container(), temp );
    M_backend->prod( *M_matrixStab, M_sol.container(), *temp );
    return M_backend->dot( *temp, M_sol.container() );
}


template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename ItRange, typename EsigmaInc,
         typename EnuInc, typename EnuAbs,
         typename Ebeta, typename Ef, typename Ec, typename Eg,
         typename EnoSlip>
void Oseen<Space, imOrder, Entity>::update( const ItRange& itRange,
        const EsigmaInc& sigmaInc,
        const EnuInc& nuInc,
        const EnuAbs& nuAbs,
        const Ebeta& beta,
        const Ef& f,
        const Ec& c,
        const Eg& g,
        const EnoSlip& noSlip,
        bool updateStabilization )
{
    element_u_type& v = u;
    element_p_type& q = p;

    M_updated = true;

    using namespace Feel::vf;

    auto bcCoeffVeloFull = ( nuAbs )*( noSlip )*M_bcCoeffDiff/hFace()-chi( ( trans( beta )*N() ) < 0 ) * trans( beta )*N();
    auto bcCoeffVeloNorm = M_bcCoeffConv * vf::max( vf::sqrt( trans( beta )*( beta ) ), ( nuAbs )/hFace() );

    // --- right hand side rhs

    // rhs volume terms
    DVLOG(2) << "[Oseen::update] (f,v)+(c,q)\n";
    form1( M_space, M_vectorRhsFull, _init=true ) =
        integrate( elements( M_mesh ), val( trans( f ) )*id( v ) + val( c )*id( q ) );
    DVLOG(2) << "[Oseen::update] (f,v)+(c,q) done\n";

    if ( M_weak_dirichlet )
    {
        // rhs boundary terms
        for ( auto diriIter = M_dirichletFlags.begin(),
                diriEnd = M_dirichletFlags.end();
                diriIter != diriEnd; ++diriIter )
        {
            DVLOG(2) << "[Oseen::update] <bc(g),v>_Gamma_"
                           << *diriIter << "\n";
            form1( M_space, M_vectorRhsFull ) +=
                integrate( markedfaces( M_mesh, *diriIter ),
                           ( -val( noSlip )*( val( nuAbs )*trans( N() )
                                              *( grad( v )+trans( grad( v ) ) )
                                              + id( q )*trans( N() ) )
                             + val( bcCoeffVeloFull )*trans( id( v ) )
                             + val( bcCoeffVeloNorm )*( trans( id( v ) )*N() )*trans( N() ) )
                           * val( g )
                         );
            DVLOG(2) << "[Oseen::update] <bc(g),v>_Gamma_"
                           << *diriIter << " done\n";
        }
    } // M_weak_dirichlet

    // --- full matrix

    // incremental volume terms
    DVLOG(2) << "[Oseen::update] adding (nu D(u),D(v)) + (sigma u,v)\n";
    form2( M_space, M_space, _matrix=M_matrixAu ) +=
        integrate( itRange,
                   val( nuInc )*trace( grad( v )*( gradt( u )+trans( gradt( u ) ) ) )
                   + val( sigmaInc )*trans( id( v ) )*idt( u )
                 );
    DVLOG(2) << "[Oseen::update] adding (nu D(u),D(v)) + (sigma u,v) done\n";

    //M_matrixFull = M_matrixAu;
    //M_matrixFull->addMatrix( 1.0, M_matrixAu );

    // convective terms
    DVLOG(2) << "[Oseen::update] ((beta*grad)u,v)\n";
    form2( M_space, M_space, M_matrixFull, _init=true ) =
        integrate( elements( M_mesh ),
                   trans( id( v ) ) * ( gradt( u )*val( beta ) )
                 );
#if 0

    // --- stabilization matrices
    if ( ( M_stabCoeffDiv != 0.0 ) ||
            ( M_stabCoeffP   != 0.0 ) ||
            ( M_epsCompress  != 0.0 ) )
    {
        if ( updateStabilization )
        {
            form2( M_space, M_space, M_matrixStab, _init=true );
#if 1

            if ( M_stabCoeffDiv != 0.0 )
            {
                DVLOG(2) << "[Oseen::update] <gamma_div [div u],[div v]>_Gamma_int\n";
                form( M_space, M_space, M_matrixStab ) +=
                    integrate( internalfaces( M_mesh ), typename im<imOrder-1>::type(),
                               val( M_stabCoeffDiv * vf::pow( hFace(),2.0 ) *
                                    vf::sqrt( trans( beta )*( beta ) ) )
                               * trans( jump( div( v ) ) ) * jumpt( divt( u ) )
                             );
            }

            if ( M_stabCoeffP != 0.0 )
            {
                DVLOG(2) << "[Oseen::update] <gamma_p [grad p],[grad q]>_Gamma_int\n";
                M_stabWeightP =
                    vf::project( M_space_i, elements( M_mesh ),
                                 M_stabCoeffP * vf::pow( h(),3.0 ) /
                                 max( h() * vf::sqrt( trans( beta ) * ( beta ) ),
                                      ( nuAbs ) )
                               );

                if ( M_pStab )
                {
                    form( M_space, M_space, M_matrixStab ) +=
                        integrate( internalfaces( M_mesh ),
                                   typename im<2*pOrder-2>::type(),
                                   val( M_thetaPstab
                                        * maxface( idv( M_stabWeightP ) ) )
                                   * ( leftface  ( grad ( q )*N() ) *
                                       leftfacet ( gradt( p )*N() ) +
                                       rightface ( grad ( q )*N() ) *
                                       rightfacet( gradt( p )*N() ) )
                                 );
                }

                else
                {
                    form( M_space, M_space, M_matrixStab ) +=
                        integrate( internalfaces( M_mesh ),
                                   typename im<2*pOrder-2>::type(),
                                   maxface( idv( M_stabWeightP ) )
                                   * jump( grad( q ) ) * jumpt( gradt( p ) )
                                 );
                }
            }

#endif

            if ( M_epsCompress > 0 )
            {
                VLOG(1) << "[Oseen] adding pseudo compressibility term with coeff= " << M_epsCompress << "\n";
                form( M_space, M_space, M_matrixStab ) +=
                    integrate( elements( M_mesh ), typename im<2*pOrder>::type(),
                               M_epsCompress * id( q ) * idt( p )
                             );
                VLOG(1) << "[Oseen] adding pseudo compressibility term with coeff= " << M_epsCompress << " done\n";
            }

            M_matrixStab->close();
        }

    }

#endif

#if 0

    if ( M_pStab && ( M_stabCoeffP != 0.0 ) )
    {
        DVLOG(2) << "[Oseen::update] added rhs stab\n";
        const int staborder = 2*pOrder-2;
        form( M_space, *M_vectorRhsFull ) +=
            integrate( internalfaces( M_mesh ), typename im<staborder>::type(),
                       val( maxface( idv( M_stabWeightP ) ) )
                       * ( - ( leftface  ( grad ( q )*N() ) *
                               rightfacev( gradv( *M_pStab )*N() ) +
                               rightface ( grad ( q )*N() ) *
                               leftfacev ( gradv( *M_pStab )*N() ) )
                           + ( M_thetaPstab - 1.0 )
                           * ( ( leftface  ( grad ( q )*N() ) *
                                 leftfacev ( gradv( *M_pStab )*N() ) +
                                 rightface ( grad ( q )*N() ) *
                                 rightfacev ( gradv( *M_pStab )*N() ) )
                             )
                         )
                     );
        DVLOG(2) << "[Oseen::update] added rhs stab done\n";
    }

#endif

    // boundary terms
    for ( auto diriIter = M_dirichletFlags.begin(),
            diriEnd = M_dirichletFlags.end();
            diriIter != diriEnd; ++diriIter )
    {
        if ( M_weak_dirichlet )
        {
            DVLOG(2) << "[Oseen::update] <bc(u),v>_Gamma_"
                           << *diriIter << "\n";
            form2( M_space, M_space, M_matrixFull ) +=
                integrate( markedfaces( M_mesh, *diriIter ),
                           -val( noSlip )*( val( ( nuAbs )*trans( N() ) )
                                            *( ( grad( v )+trans( grad( v ) ) ) *idt( u ) +
                                               ( gradt( u )+trans( gradt( u ) ) )*id( v ) )
                                            + id( q )*trans( N() )*idt( u )
                                            - idt( p )*trans( N() )*id( v ) )
                           + val( bcCoeffVeloFull )*( trans( id( v ) )*idt( u ) )
                           + val( bcCoeffVeloNorm )
                           * ( trans( id( v ) )*N() ) * ( trans( N() )*idt( u ) )
                         );
            DVLOG(2) << "[Oseen::update] <bc(u),v>_Gamma_"
                           << *diriIter << " done\n";
        }

        else
        {
            M_matrixFull->close();
            M_vectorRhsFull->close();
            DVLOG(2) << "[Oseen::update] on(u)_Gamma_"
                           << *diriIter << "\n";
            form2( M_space, M_space, M_matrixFull ) +=
                on( markedfaces( *M_mesh, *diriIter ),
                    u, *M_vectorRhsFull, g );
            DVLOG(2) << "[Oseen::update] on(u)_Gamma_"
                           << *diriIter << "done \n";
        }
    }

    //
    // close global matrix and vector to make them usable by the
    // solvers
    //
    M_matrixFull->close();
    M_vectorRhsFull->close();

    // added pressure and viscous terms. \attention addMatrix can be
    // called only on closed matrices
    //M_matrixFull->addMatrix( 1.0, M_matrixStab );
    VLOG(1) << "[Oseen] added stabilisation matrix\n";
    M_matrixFull->addMatrix( 1.0, M_matrixAu );
    VLOG(1) << "[Oseen] added standard oseen matrix\n";

    if ( M_export_matlab )
    {
        M_matrixFull->printMatlab( "M_oseen.m" );
        M_vectorRhsFull->printMatlab( "F_oseen.m" );
    }

} // update

template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void Oseen<Space, imOrder, Entity>::solve()
{
    // -- make sure solve is needed
    if ( !M_updated )
    {
        return;
    }

    M_updated = false;

#if 0
    // -- solve
    solveNonSym( M_matrixFull, M_vectorSolFull, M_vectorRhsFull );

    // -- extract solution
    std::copy( M_vectorSolFull.begin(),
               M_vectorSolFull.end(),
               M_sol.begin() );
#else

    // -- solve
    solveNonSym( M_matrixFull, M_sol, M_vectorRhsFull );

#endif


} // Oseen::solve


template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void
Oseen<Space, imOrder, Entity>::solveNonSym( sparse_matrix_ptrtype const& D,
        element_type  & u,
        vector_ptrtype  const& F )
{
    M_backend->solve( _matrix=D, _solution=u, _rhs=F );
} // Oseen::solveNonSym



} // Feel
