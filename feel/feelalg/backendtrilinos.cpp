/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2007-05-30

   Copyright (C) 2007, 2009 Universit� Joseph Fourier (Grenoble I)

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
   \file backendpetsc.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-05-30
*/
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backendtrilinos.hpp>

namespace Feel
{
#if defined ( FEELPP_HAS_TRILINOS_EPETRA )

/**
 * \return the command lines options of the trilinos backend
 */
    po::options_description backendtrilinos_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "BackendTrilinos " + prefix + " solver options" );
    _options.add_options()
    // solver options
    ( ( _prefix+"trilinos-solver-type" ).c_str(), Feel::po::value<std::string>()->default_value( "gmres" ), "cg, bicgstab, gmres, cg-cond, gmres-cond" )

    // preconditioner options
    ( ( _prefix+"trilinos-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "block-jacobi" ), "id, domain-decomposition, block-jacobi, neumann" )
    ( ( _prefix+"trilinos-drop" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "drop tolerance for preconditioners" )
    ( ( _prefix+"trilinos-fillin" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "fill-in level value for preconditioners" )
    ( ( _prefix+"trilinos-overlap" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "domain overlap value for preconditioners" )
    ( ( _prefix+"trilinos-overlap-type" ).c_str(), Feel::po::value<std::string>()->default_value( "symmetric" ), "type of domain decomposition overlap for preconditioners" )
    ( ( _prefix+"trilinos-order" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "order of the Neumann polynomial expansion for preconditioners" )
    ( ( _prefix+"trilinos-local-parts" ).c_str(), Feel::po::value<int>()->default_value( 4 ), "check Trilinos documentation" )
    ( ( _prefix+"trilinos-local-solver-type" ).c_str(), Feel::po::value<std::string>()->default_value( "Amesos" ), "Local domain subsolver for Ifpack" )
    ( ( _prefix+"trilinos-local-solver" ).c_str(), Feel::po::value<std::string>()->default_value( "Amesos_Umfpack" ), "(Amesos_KLU, Amesos_Umfpack, Amesos_Superlu ) Local domain direct subsolver for Ifpack" )


    // solver control options
    ( ( _prefix+"trilinos-restart" ).c_str(), Feel::po::value<int>()->default_value( 20 ), "number of iterations before solver restarts (gmres)" )
    ( ( _prefix+"trilinos-verbose" ).c_str(), Feel::po::value<std::string>()->default_value( "none" ), "(none, summary, last, all) print solver iterations" )
    ( ( _prefix+"trilinos-maxiter" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "set maximum number of iterations" )
    ( ( _prefix+"trilinos-tolerance" ).c_str(), Feel::po::value<double>()->default_value( 2e-10 ), "set solver tolerance" )
    ( ( _prefix+"trilinos-residual" ).c_str(), Feel::po::value<std::string>()->default_value( "rhs" ), "set solver residual calculation" )
    ;
    return _options;
}


    BackendTrilinos::BackendTrilinos( po::variables_map const& vm, std::string const& prefix, WorldComm const& worldComm )
    :
    super( vm, prefix, worldComm ),
    M_options(),
    M_prec_type( "" ),
    M_Prec()
{
    this->M_backend = BackendType::BACKEND_TRILINOS;

    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    set_maxiter( vm[_prefix+"trilinos-maxiter"].as<int>() );
    set_tol( vm[_prefix+"trilinos-tolerance"].as<double>() );
    set_verbose( vm[_prefix+"trilinos-verbose"].as<std::string>() );

    set_drop( vm[_prefix+"trilinos-drop"].as<double>() );
    set_overlap( vm[_prefix+"trilinos-overlap"].as<int>() );
    set_restart( vm[_prefix+"trilinos-restart"].as<int>() );

    set_fillin( vm[_prefix+"trilinos-fillin"].as<int>() );
    set_order( vm[_prefix+"trilinos-order"].as<int>() );
    set_residual( vm[_prefix+"trilinos-residual"].as<std::string>() );

    set_solver( vm[_prefix+"trilinos-solver-type"].as<std::string>() );
    set_prec( vm[_prefix+"trilinos-pc-type"].as<std::string>() );

    if ( M_prec_type == "Ifpack" )
    {
        set_prec_local_solver_type( vm[_prefix+"trilinos-local-solver-type"].as<std::string>() );
        set_overlap_type( vm[_prefix+"trilinos-overlap-type"].as<std::string>() );
        set_localparts( vm[_prefix+"trilinos-local-parts"].as<int>() );
        set_prec_local_solver( vm[_prefix+"trilinos-local-solver"].as<std::string>() );
    }
}


BackendTrilinos::BackendTrilinos( const BackendTrilinos& tc )
    :
    super( tc ),
    M_options( tc.M_options ),
    M_prec_type( tc.M_prec_type ),
    M_local_solver( tc.M_local_solver ),
    M_local_solver_type( tc.M_local_solver_type ),
    M_Prec( tc.M_Prec )

{
}



// -- SETTING OF OPTIONS --
void
BackendTrilinos::set_maxiter ( int maxiter )
{
    M_options.set( "max_iter", maxiter );
}

void
BackendTrilinos::set_tol ( double tol )
{
    M_options.set( "tol", tol );
}

void
BackendTrilinos::set_drop ( double drop )
{
    M_options.set( "drop", drop );
    M_options.set( "fact: drop tolerance", drop );
}

void
BackendTrilinos::set_overlap ( int overlap )
{
    M_options.set( "overlap", overlap );
    M_options.set( "partitioner: overlap", overlap );
}

void
BackendTrilinos::set_restart ( int restart )
{
    M_options.set( "kspace", restart );
}

void
BackendTrilinos::set_fillin ( int fillin )
{
    M_options.set( "ilut_fill", fillin );
    M_options.set( "fact: level-of-fill", fillin );
}

void
BackendTrilinos::set_order ( int order )
{
    M_options.set( "poly_order", order );
}

void
BackendTrilinos::set_overlap_type( std::string str )
{
    M_options.set( "type_overlap", str );
}

void
BackendTrilinos::set_residual ( std::string str )
{
    M_options.set( "conv", str );
}

void
BackendTrilinos::set_localparts ( int localparts )
{
    M_options.set( "partitioner: local parts", localparts );
}

void
BackendTrilinos::set_solver( std::string str )
{
    if ( str == "gmres" )
        M_options.set( "solver", "AZ_gmres" );

    else if ( str == "gmres-cond" )
        M_options.set( "solver", "AZ_gmres_condnum" );

    else if ( str == "cg" )
        M_options.set( "solver", "AZ_cg" );

    else if ( str == "cg-cond" )
        M_options.set( "solver", "AZ_cg_condnum" );

    else if ( str == "bicgstab" )
        M_options.set( "solver", "AZ_bicgstab" );
}

void
BackendTrilinos::set_prec_local_solver( std::string str )
{
    M_local_solver = str;
}

void
BackendTrilinos::set_prec_local_solver_type( std::string str )
{
    M_local_solver_type = str;
}

void
BackendTrilinos::set_prec( std::string str )
{
    if ( str == "id" )
        M_options.set( "precond", "none" );

    else if ( str == "block-jacobi" )
        M_options.set( "precond", "Jacobi" );

    else if ( str == "domain-decomposition" )
        M_options.set( "precond", "dom_decomp" );

    else if ( str == "neumann" )
        M_options.set( "precond", "Neumann" );

    else if ( str == "ifpack" )
        M_prec_type = "Ifpack";

    else if ( str == "ml" )
        M_prec_type = "MultiLevel";
}

void
BackendTrilinos::set_verbose( const std::string verb )
{
    if ( verb == "none" )
        M_options.set( "output", "none" );

    else if ( verb == "summary" )
        M_options.set( "output", "summary" );

    else if ( verb == "last" )
        M_options.set( "output", "last" );

    else
        M_options.set( "output", "all" );
}

void
BackendTrilinos::set_options( list_type opts )
{
    M_options = opts;
}

BackendTrilinos::list_type
BackendTrilinos::get_options()
{
    return M_options;
}



BackendTrilinos::solve_return_type
BackendTrilinos::solve( base_sparse_matrix_ptrtype const& A,
                        base_sparse_matrix_ptrtype const& B,
                        base_vector_ptrtype& x,
                        base_vector_ptrtype const& b )
{
#if 0
    BackendTrilinos opt;
    opt.set_verbose( "none" );
    opt.set_tol( 1e-14 );
    opt.set_solver( "gmres" );
    opt.set_residual( "rhs" );
#endif
    operator_ptrtype prec = IfpackPrec( B, this->get_options() );

    op_mat_ptrtype matrixOp = op_mat_ptrtype( new op_mat_type( A, this->get_options(), "A", prec ) );

    int iter = matrixOp->ApplyInverse( b, x );
    return boost::make_tuple( true, iter, 1e-16 );
}

// -- LINEAR ALGEBRA INTERFACE --
void
BackendTrilinos::applyMatrix( sparse_matrix_type const& _A,
                              const Vector<value_type>& _x,
                              vector_type& _b )
{
    epetra_sparse_matrix_type const& A( dynamic_cast<epetra_sparse_matrix_type const&>( _A ) );

    epetra_vector_type& b = dynamic_cast<epetra_vector_type&>( _b );

    try
    {
        epetra_vector_type const& x = dynamic_cast<epetra_vector_type const&>( _x );
        A.mat().Apply( x.vec(), b.vec() );
    }

    catch ( std::bad_cast& )
    {

        Epetra_Map map( _x.size(), _x.localSize(), 0, Epetra_MpiComm( _x.comm() ) );
        epetra_vector_type v_ep( map );

        //v_ep = x;
        for ( size_type lid =0; lid < _x.localSize(); ++lid )
        {
            size_type gid = _x.firstLocalIndex()+lid;
            v_ep( gid  ) = _x( gid );
        }

        A.mat().Apply( v_ep.vec(), b.vec() );
    }
}

BackendTrilinos::real_type
BackendTrilinos::dot( const vector_type& _f,
                      const vector_type& _x ) const
{
    epetra_vector_type const& x = dynamic_cast<epetra_vector_type const&>( _x );
    epetra_vector_type const& f = dynamic_cast<epetra_vector_type const&>( _f );
    value_type result;
    f.vec().Dot( x.vec(), &result );
    return result;
}

#endif // FEELPP_HAS_TRILINOS_EPETRA
} // Feel

