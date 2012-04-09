/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-05-22

  Copyright (C) 2007, 2009 Université Joseph Fourier (Grenoble I)

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
   \file backendgmm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-05-22
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backendgmm.hpp>
#include <feel/feelalg/matrixtriplet.hpp>

namespace Feel
{
// -- CONSTRUCTOR --
template<typename T>
BackendGmm<T>::BackendGmm()
    :
    super(),
    M_iter( 2.e-10, 0, 1000 ), // tol, noisy and maxiter
    M_fillin( 2 ),
    M_threshold( 1.e-3 ),
    M_isSymmetric( false ),
    M_isDirect( true ),
    M_ilut( 0 ),
    M_ilutp( 0 ),
    M_solver_type( "umfpack" ),
    M_precond_type( "ilut" ),
    M_restart( 200 )
#if defined( FEELPP_HAS_UMFPACK )
    , M_umfpack()
#endif // FEELPP_HAS_UMFPACK
{}

template<typename T>
BackendGmm<T>::BackendGmm( po::variables_map const& vm, std::string const& prefix )
    :
    super( vm, prefix ),
    M_iter( 2.e-10, 0, 1000 ), // tol, noisy and maxiter
    M_fillin( 2 ),
    M_threshold( 1.e-3 ),
    M_isSymmetric( false ),
    M_isDirect( true ),
    M_ilut( 0 ),
    M_ilutp( 0 ),
    M_solver_type( "umfpack" ),
    M_precond_type( "ilut" ),
    M_restart( 200 )
#if defined( FEELPP_HAS_UMFPACK )
    , M_umfpack()
#endif // FEELPP_HAS_UMFPACK
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    this->set_solver_type( vm[_prefix+"gmm-solver-type"].template as<std::string>() );
    this->set_preconditioner_type( vm[_prefix+"gmm-pc-type"].template as<std::string>() );
    this->set_noisy( vm[_prefix+"gmm-verbose"].template as<int>() );
    this->set_maxiter( vm[_prefix+"gmm-maxiter"].template as<int>() );
    this->set_fillin( vm[_prefix+"gmm-fillin"].template as<int>() );
    this->set_threshold( vm[_prefix+"gmm-threshold"].template as<double>() );
    this->set_tol( vm[_prefix+"gmm-tolerance"].template as<double>() );
}

template<typename T>
boost::shared_ptr<MatrixTriplet<T> >
BackendGmm<T>:: toTriplet( sparse_matrix_type const& _m )
{
    gmm_sparse_matrix_type const& m( dynamic_cast<gmm_sparse_matrix_type const&>( _m ) );
    int nr = m.mat().nr;
    int nc = m.mat().nc;
    int nz = m.mat().jc[ nr ];
    std::vector<int> Ti( nz );
    std::vector<int> Tj( nz );
    std::vector<double> Tx( nz );
    unsigned int const* jc = m.mat().jc;
    typedef typename gmm_sparse_matrix_type::matrix_type csr_matrix_type;
    typedef typename gmm::linalg_traits<csr_matrix_type>::const_sub_row_type row_type;

    for ( int j = 0; j < nr; ++j )
    {
        row_type row = gmm::mat_const_row( m.mat(), j );
        typename gmm::linalg_traits<row_type>::const_iterator
        it = gmm::vect_const_begin( row ), ite = gmm::vect_const_end( row );

        for ( int k = 0; it != ite; ++it, ++k )
        {
            Ti[jc[j] + k] = j;
            Tj[jc[j] + k] = it.index();
            Tx[jc[j] + k] = *it;
        }
    }

    return boost::shared_ptr<MatrixTriplet<T> >( new MatrixTriplet<T>( nr, nc, Ti, Tj, Tx ) );
}

template<typename T>
typename BackendGmm<T>::solve_return_type
BackendGmm<T>::solve( sparse_matrix_type const& _A,
                      vector_type& _x,
                      const vector_type& _b )
{
    bool reusePC = ( this->precMatrixStructure() == SAME_PRECONDITIONER );

    gmm_sparse_matrix_type const& A( dynamic_cast<gmm_sparse_matrix_type const&>( _A ) );
    gmm_vector_type      & x( dynamic_cast<gmm_vector_type      &>( _x ) );
    gmm_vector_type const& b( dynamic_cast<gmm_vector_type const&>( _b ) );

    if ( M_solver_type == "superlu" )
    {
#if defined(GMM_USES_SUPERLU)
        Debug() << "[BackendGmm] trying direct solve with SuperLU\n";
        double condest;

        try
        {
            gmm::SuperLU_solve( A.mat(), x, b, condest );
        }

        catch ( dal::failure_error& e )
        {
            Debug() << "[BackendGmm] " << e.what() << "\n";
        }

        Debug() << "[BackendGmm] SuperLU condition estimate = "
                << condest << "\n";
        return boost::make_tuple( false, -1, 0 );
#else
        Error() << "Invalid solver type : " << M_solver_type << "\n";
#endif
    }

    else if ( M_solver_type == "umfpack" )
    {
#if defined( FEELPP_HAS_UMFPACK )
        M_umfpack.setMatrix( *BackendGmm<value_type>::toTriplet( A ) );

        if ( M_isSymmetric )
            M_umfpack.setStrategy( UMFPACK_STRATEGY_SYMMETRIC );

        VectorUblas<value_type> xx( x.size() );
        M_umfpack.solve( xx.vec(), b.vec() );
        std::copy( xx.begin(), xx.end(), x.begin() );

        return boost::make_tuple( false, -1, 0 );
#else
        Error() << "Invalid solver type : " << M_solver_type << "\n";
#endif
    }

    std::vector<value_type> sx( x.size() );
    std::copy( x.begin(), x.end(), sx.begin() );
    std::vector<value_type> sb( b.size() );
    std::copy( b.begin(), b.end(), sb.begin() );

    M_iter.init();

    if ( M_isSymmetric )
    {
        //         gmm::ildltt_precond<typename sparse_matrix_type::matrix_type>
        //             P(A.mat(), M_fillin, M_threshold);
        gmm::identity_matrix P;
        gmm::cg( A.mat(), sx, sb, P, M_iter );
    }

    else
    {
        if ( M_precond_type == "ilut" )
        {
            if ( !reusePC || ( M_ilut.get() == 0 ) )
            {
                // delete old preconditioner before creating new
                M_ilut.reset( 0 );
                // delete ilutp preconditioner when using ilut
                M_ilutp.reset( 0 );
                M_ilut.reset( new ilut_type( A.mat(),
                                             M_fillin,
                                             M_threshold ) );
            }

            if ( M_solver_type == "bicgstab" )
                gmm::bicgstab( A.mat(), sx, sb, *M_ilut, M_iter );

            if ( M_solver_type == "gmres" )
                gmm::gmres( A.mat(), sx, sb, *M_ilut, M_restart,
                            M_iter );
        }

        if ( M_precond_type == "ilutp" )
        {
            if ( !reusePC || ( M_ilutp.get() == 0 ) )
            {
                // delete ilut preconditioner when using ilutp
                M_ilut.reset( 0 );
                // delete old preconditioner before creating new
                M_ilutp.reset( 0 );
                M_ilutp.reset( new ilutp_type( A.mat(),
                                               M_fillin,
                                               M_threshold ) );
            }

            if ( M_solver_type == "bicgstab" )
                gmm::bicgstab( A.mat(), sx, sb, *M_ilutp, M_iter );

            if ( M_solver_type == "gmres" )
                gmm::gmres( A.mat(), sx, sb, *M_ilutp, M_restart,
                            M_iter );
        }

        if ( M_precond_type == "diag" )
        {
            diagprec_type P( A.mat() );

            if ( M_solver_type == "bicgstab" )
                gmm::bicgstab( A.mat(), sx, sb, P, M_iter );

            if ( M_solver_type == "gmres" )
                gmm::gmres( A.mat(), sx, sb, P, M_restart, M_iter );
        }

        if ( M_precond_type == "id" )
        {
            gmm::identity_matrix P;

            if ( M_solver_type == "bicgstab" )
                gmm::bicgstab( A.mat(), sx, sb, P, M_iter );

            if ( M_solver_type == "gmres" )
                gmm::gmres( A.mat(), sx, sb, P, M_restart, M_iter );
        }
    }

    std::copy( sx.begin(), sx.end(), x.begin() );

    return boost::make_tuple( M_iter.converged(), M_iter.get_iteration(), M_iter.residual() );

} // BackendGmm::solve

//
// Instantiation
//
template class BackendGmm<double>;

/**
 * \return the command lines options of the gmm backend
 */
po::options_description backendgmm_options( std::string const& prefix,
        BackendGmmDefaults defaults )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "BackendGmm " + prefix + " solver options" );
    _options.add_options()
    // solver options
    ( ( _prefix+"gmm-solver-type" ).c_str(),
      Feel::po::value<std::string>()->default_value( defaults.solver_type ),
      "umfpack, superlu, cg, bicgstab, gmres" )

    // preconditioner options
    ( ( _prefix+"gmm-pc-type" ).c_str(),
      Feel::po::value<std::string>()->default_value( defaults.pc_type ),
      "ilut, ilutp, diag, id" )
    ( ( _prefix+"gmm-threshold" ).c_str(),
      Feel::po::value<double>()->default_value( defaults.threshold ),
      "threshold value for preconditioners" )
    ( ( _prefix+"gmm-fillin" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.fillin ),
      "fill-in level value for preconditioners" )

    // solver control options
    ( ( _prefix+"gmm-restart" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.restart ),
      "number of iterations before solver restarts (gmres)" )
    ( ( _prefix+"gmm-verbose" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.verbose ),
      "(=0,1,2) print solver iterations" )
    ( ( _prefix+"gmm-maxiter" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.maxiter ),
      "set maximum number of iterations" )
    ( ( _prefix+"gmm-tolerance" ).c_str(),
      Feel::po::value<double>()->default_value( defaults.tolerance ),
      "set solver tolerance" )
    ;
    return _options;
}

} // Feel
