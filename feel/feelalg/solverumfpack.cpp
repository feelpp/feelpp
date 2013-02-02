/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-04

  Copyright (C) 2007, 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2004 EPFL


  This library is free software; you can redlstribute it and/or
  modlfy it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is dlstributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file solverumfpack.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-04
 */
#include <boost/timer.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/solverumfpack.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#if defined(FEELPP_HAS_UMFPACK)

namespace Feel
{
class SolverUMFPACK::Pimpl
{
public:

    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::column_major> matrix_type;

    Pimpl()
        :
        nr( 0 ),
        nc( 0 ),
        Ap(),
        Ai(),
        Ax()
    {}

    Pimpl( Pimpl const& __p )
        :
        nr( __p.nr ),
        nc( __p.nc ),
        Ap( __p.Ap ),
        Ai( __p.Ai ),
        Ax( __p.Ax )
    {}

    int nr;
    int nc;
    std::vector<int> Ap;
    std::vector<int> Ai;
    std::vector<double> Ax;
};
SolverUMFPACK::SolverUMFPACK()
    :
    M_p( new Pimpl ),
    M_matrix_reset( true ),
    M_matrix_values_reset( true ),
    M_symbolic( 0 ),
    M_numeric( 0 ),
    M_Control( new double[UMFPACK_CONTROL] ),
    M_Info( new double[UMFPACK_INFO] )
{
    ::umfpack_di_defaults( M_Control );
    M_Control[UMFPACK_PRL] = 6;
    //M_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    //M_Control[UMFPACK_AMD_DENSE] = 1;
}
SolverUMFPACK::SolverUMFPACK( SolverUMFPACK const & umfpackSolver )
    :
    M_p( umfpackSolver.M_p ),
    M_matrix_reset( umfpackSolver.M_matrix_reset ),
    M_matrix_values_reset( umfpackSolver.M_matrix_values_reset ),
    M_symbolic( 0 ),
    M_numeric( 0 ),
    M_Control( umfpackSolver.M_Control ),
    M_Info( umfpackSolver.M_Info )
{

}
SolverUMFPACK::~SolverUMFPACK()
{
    if ( M_numeric )
        ::umfpack_di_free_numeric( &M_numeric );

    if ( M_symbolic )
        ::umfpack_di_free_symbolic( &M_symbolic );

}

void
SolverUMFPACK::reportInfo()
{
    umfpack_di_report_info( M_Control, M_Info );
}
void
SolverUMFPACK::reportStatus( int status )
{
    umfpack_di_report_status( M_Control, status );
}
void
SolverUMFPACK::setStrategy( int strategy )
{
    if  ( strategy == ( UMFPACK_STRATEGY_AUTO ||
                        UMFPACK_STRATEGY_UNSYMMETRIC ||
                        UMFPACK_STRATEGY_SYMMETRIC ||
                        UMFPACK_STRATEGY_2BY2 ) )
        M_Control [UMFPACK_STRATEGY] = strategy;

    else
        M_Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_AUTO;
}
void
SolverUMFPACK::setMatrix( const matrix_type& m )
{
    DVLOG(2) << "[SolverUMFPACK::setMatrix] set matrix (nr="
                  << m.nrows() << ", nc=" << m.ncols() << ", nz=" << m.nz() << ")\n";
    boost::timer ti;
    int *Map = 0;
    M_p->nr = m.nrows();
    M_p->nc = m.ncols();


    M_p->Ap.resize( m.ncols()+1 );
    M_p->Ai.resize( m.nz() );
    M_p->Ax.resize( m.nz() );

    int status = umfpack_di_triplet_to_col ( m.nrows(), m.ncols(), m.nz(),
                 m.Ti(), m.Tj(), m.Tx(),
                 &M_p->Ap[0], &M_p->Ai[0], &M_p->Ax[0],
                 Map ) ;

    if ( status != UMFPACK_OK )
    {
        reportInfo();
        reportStatus( status );
        Error() << "[SolverUMFPACK::setMatrix] set matrix failed\n";

    }

    DVLOG(2) << "[SolverUMFPACK::setMatrix] set matrix done in " << ti.elapsed() << "\n";
}

void
SolverUMFPACK::solve( array_type& __X, array_type const& __B )
{

    prepareSolve();

    boost::timer ti;

    DVLOG(2) << "[SolverUMFPACK::solve] solve A x = b using UMFPACK version " << UMFPACK_VERSION << "\n";
    int status = umfpack_di_solve( UMFPACK_A,
                                   &M_p->Ap[0],
                                   &M_p->Ai[0],
                                   &M_p->Ax[0],
                                   boost::addressof( __X[0] ),
                                   boost::addressof( __B[0] ),
                                   M_numeric,
                                   M_Control,
                                   M_Info );

    if ( status != UMFPACK_OK )
    {
        reportInfo();
        reportStatus( status );
        Error() << "[SolverUMFPACK::solve] solve failed\n";
    }

    DVLOG(2) << "[SolverUMFPACK::solve] solve in " << ti.elapsed() << "\n";
}
void SolverUMFPACK::prepareSolve()
{
    boost::timer ti;
    DVLOG(2) << "[SolverUMFPACK::solve] preparesolve starts\n";

    if ( M_matrix_reset )
    {
        if ( M_symbolic )
        {
            DVLOG(2) << "[SolverUMFPACK::prepareSolve] Destroying symbolic factorization\n";

            umfpack_di_free_symbolic( &M_symbolic );
            M_symbolic = 0;
        }

        DVLOG(2) << "[SolverUMFPACK::prepareSolve] computing symbolic factorization\n";
        int status = umfpack_di_symbolic( M_p->nr,
                                          M_p->nc,
                                          &M_p->Ap[0],
                                          &M_p->Ai[0],
                                          &M_p->Ax[0],
                                          &M_symbolic,
                                          M_Control,
                                          M_Info );

        if ( status != UMFPACK_OK )
        {
            reportInfo();
            reportStatus( status );
            Error() << "[SolverUMFPACK::prepareSolve] symbolic factorization failed\n";
        }
    }

    if ( M_matrix_reset || M_matrix_values_reset )
    {
        if ( M_numeric )
        {
            DVLOG(2) << "[SolverUMFPACK::prepareSolve] Destroying numeric factorization\n";
            umfpack_di_free_numeric( &M_numeric );
            M_numeric = 0;
        }

        DVLOG(2) << "[SolverUMFPACK::prepareSolve] computing numeric factorization\n";
        int status = umfpack_di_numeric( &M_p->Ap[0],
                                         &M_p->Ai[0],
                                         &M_p->Ax[0],
                                         M_symbolic, &M_numeric,
                                         M_Control,
                                         M_Info );

        if ( status != UMFPACK_OK )
        {
            reportInfo();
            reportStatus( status );
            Error() << "[SolverUMFPACK::prepareSolve] numeric factorization failed\n";
        }
    }

    M_matrix_reset = false;
    M_matrix_values_reset = false;
    DVLOG(2) << "[SolverUMFPACK::solve] preparesolve done in " << ti.elapsed() << "\n";
}

}
#endif // FEELPP_HAS_UMFPACK
