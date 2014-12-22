/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
              Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2007-01-18

   Copyright (C) 2007 EPFL
   Copyright (C) 2008-2011 Université Joseph Fourier (Grenoble 1)


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
   \file backendtrilinos.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-08-14
*/
#ifndef _BACKENDTRILINOS_HPP_
#define _BACKENDTRILINOS_HPP_


#include <boost/program_options/variables_map.hpp>
#include <feel/feelconfig.h>

#if defined ( FEELPP_HAS_TRILINOS_EPETRA )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#include <feel/feelalg/matrixepetra.hpp>
#include <feel/feelalg/vectorepetra.hpp>
#include <feel/feelalg/preconditionerifpack.hpp>
#include <feel/feelalg/preconditionerml.hpp>
#include <feel/feelalg/solverlineartrilinos.hpp>
#include <feel/feelalg/operatortrilinos.hpp>
#include <feel/feelalg/datamap.hpp>
#include <feel/feelalg/backend.hpp>


#include <Teuchos_ParameterList.hpp>

namespace Feel
{
namespace po = boost::program_options;

// forward declararion
template<class T> class Backend;
class OperatorMatrix;

/**
 * \class BackendTrilinos
 * \brief backend interface for trilinos framework
 *
 */
class BackendTrilinos : public Backend<double>
{
    typedef Backend<double> super;
public:

    // -- TYPEDEFS --
    typedef super::value_type value_type;

    typedef super::vector_ptrtype base_vector_ptrtype;
    typedef super::sparse_matrix_ptrtype base_sparse_matrix_ptrtype;
    typedef super::vector_type vector_type;
    typedef super::sparse_matrix_type sparse_matrix_type;

    typedef sparse_matrix_type::graph_type graph_type;
    typedef sparse_matrix_type::graph_ptrtype graph_ptrtype;

    typedef MatrixEpetra epetra_sparse_matrix_type;
    typedef VectorEpetra<value_type> epetra_vector_type;
    typedef Epetra_Operator operator_type;

    typedef OperatorMatrix op_mat_type;
    typedef boost::shared_ptr<OperatorMatrix> op_mat_ptrtype;

    typedef boost::shared_ptr<epetra_sparse_matrix_type> epetra_sparse_matrix_ptrtype;
    typedef boost::shared_ptr<epetra_vector_type> epetra_vector_ptrtype;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    typedef super::solve_return_type solve_return_type;
    typedef super::nl_solve_return_type nl_solve_return_type;

    typedef std::map< size_type, std::vector< size_type > > procdist_type;

    typedef Teuchos::ParameterList list_type;


    // -- CONSTRUCTOR --
    BackendTrilinos( WorldComm const& worldComm=Environment::worldComm())
        :
        super( worldComm ),
        M_options(),
        M_prec_type( "" ),
        M_Prec()
    {
        this->M_backend = BackendType::BACKEND_TRILINOS;
        set_maxiter( 1000 );
        set_tol( 1e-10 );
    }


    BackendTrilinos( po::variables_map const& vm, std::string const& prefix = "", WorldComm const& worldComm=Environment::worldComm() );

    BackendTrilinos( const BackendTrilinos& tc );


    // -- FACTORY METHODS --
    static Epetra_Map epetraMap( DataMap const& dmap )
    {
        std::vector<int> e( dmap.nMyElements() );
        std::copy( dmap.myGlobalElements().begin(),
                   dmap.myGlobalElements().end(),
                   e.begin() );
        return Epetra_Map( -1, dmap.nMyElements(), e.data(), 0, Epetra_MpiComm( dmap.comm() ) );
    }


    static Epetra_Map epetraMapStatic( DataMap const& dmap )
    {
        return Epetra_Map( dmap.nGlobalElements(), dmap.nMyElements(), 0, Epetra_MpiComm( dmap.comm() ) );
    }


    template<typename DomainSpace, typename DualImageSpace>
    sparse_matrix_ptrtype newMatrix( DomainSpace const& Xh,
                                     DualImageSpace const& Yh,
                                     size_type matrix_properties = NON_HERMITIAN,
                                     bool init = true )
    {
        return newMatrix( Xh->map(), Yh->map(), matrix_properties );
    }

    sparse_matrix_ptrtype
    newMatrix( const size_type m,
               const size_type n,
               const size_type m_l,
               const size_type n_l,
               const size_type nnz=30,
               const size_type noz=10,
               size_type matrix_properties = NON_HERMITIAN )
    {
        sparse_matrix_ptrtype  mat( new epetra_sparse_matrix_type( m,n,m_l,n_l,nnz,noz ) );
        mat->setMatrixProperties( matrix_properties );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( const size_type m,
               const size_type n,
               const size_type m_l,
               const size_type n_l,
               graph_ptrtype const & graph,
               size_type matrix_properties = NON_HERMITIAN )
    {
        sparse_matrix_ptrtype  mat( new epetra_sparse_matrix_type( m,n,m_l,n_l,30,10 ) );
        mat->setMatrixProperties( matrix_properties );
        return mat;
    }

    sparse_matrix_ptrtype newMatrix( DataMap const& domainmap,
                                     DataMap const& imagemap,
                                     size_type matrix_properties = NON_HERMITIAN,
                                     bool init = true )
    {
        Epetra_Map erowmap = BackendTrilinos::epetraMap( imagemap );
        Epetra_Map ecolmap = BackendTrilinos::epetraMapStatic( domainmap );
        Epetra_Map edomainmap = BackendTrilinos::epetraMap( imagemap );

        //std::cout << "Rowmap: " << erowmap << "\n";
        //std::cout << "Colmap: " << ecolmap << "\n";
        //std::cout << "Domainmap: " << edomainmap << "\n";

        //std::cout << "Is matrix rectangular? " << !erowmap.SameAs( ecolmap ) << "\n";

        if ( !erowmap.SameAs( ecolmap ) )
        {
            Epetra_Map eimagemap = BackendTrilinos::epetraMap( domainmap );
            //std::cout << "Imagemap: " << eimagemap << "\n";
            auto A= sparse_matrix_ptrtype( new epetra_sparse_matrix_type( erowmap, ecolmap, edomainmap, eimagemap ) );
            A->setMatrixProperties( matrix_properties );
            return A;
        }

        else
        {
            auto A= sparse_matrix_ptrtype( new epetra_sparse_matrix_type( erowmap, ecolmap ) );
            A->setMatrixProperties( matrix_properties );
            return A;
        }
    }

    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l )
    {
        sparse_matrix_ptrtype  mat( new epetra_sparse_matrix_type( m,n,m_l,n_l,0,0 ) );
        return mat;

    }
    sparse_matrix_ptrtype
    newZeroMatrix( DataMap const& domainmap, DataMap const& imagemap )
    {
        Epetra_Map erowmap = BackendTrilinos::epetraMap( imagemap );
        Epetra_Map ecolmap = BackendTrilinos::epetraMapStatic( domainmap );

        auto A= sparse_matrix_ptrtype( new epetra_sparse_matrix_type( erowmap, ecolmap ) );
        //A->setMatrixProperties( matrix_properties );
        return A;
    }


    template<typename SpaceT>
    vector_ptrtype newVector( SpaceT const& space )
    {
        return newVector( space->map() );
    }
    vector_ptrtype newVector( DataMap const& domainmap )
    {
        Epetra_Map emap = BackendTrilinos::epetraMap( domainmap );

        return vector_ptrtype( new epetra_vector_type( emap ) );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local )
    {
        return vector_ptrtype( new epetra_vector_type( /*n, n_local*/ ) );
#warning to fix!
    }


    static operator_ptrtype IfpackPrec( sparse_matrix_ptrtype const& M, list_type options, std::string precType = "Amesos" )
    {
        PreconditionerIfpack P( options, precType );

        P.buildPreconditioner( M );

        return P.getPrec();
    }


    static operator_ptrtype MLPrec( sparse_matrix_ptrtype& M, list_type options )
    {
        PreconditionerML P( options );

        P.buildPreconditioner( M );

        return P.getPrec();
    }

#if 0
    template< typename element_type >
    static void Epetra2Ublas( vector_ptrtype const& u, element_type& x )
    {
        epetra_vector_ptrtype const& _v( dynamic_cast<epetra_vector_ptrtype const&>( u ) );
        Epetra_Map v_map( _v->Map() );

        vector_type v = *u;

        //DVLOG(2) << "Initial EpetraVector " << v << "\n";

        const size_type L = v.localSize();

        for ( size_type i=0; i<L; i++ )
        {
            DVLOG(2) << "x(" << x.firstLocalIndex() + i  << ")="
                           << "v[" << v_map.GID( i ) << "] = "
                           << v( i ) << "\n";

            x( x.firstLocalIndex() + i ) = v( i );
        }

        DVLOG(2) << "Epetra2Ublas:" << x << "\n";
    }


    template< typename element_type >
    static void Ublas2Epetra( element_type const& x, vector_ptrtype& v )
    {
        epetra_vector_type& _v( dynamic_cast<epetra_vector_type&>( *v ) );
        Epetra_Map v_map( _v.Map() );

        DVLOG(2) << "Local size of ublas vector" << x.localSize() << "\n";
        DVLOG(2) << "Local size of epetra vector" << v->localSize() << "\n";

        const size_type L = v->localSize();

        for ( size_type i=0; i<L; i++ )
        {
            DVLOG(2) << "v[" << v_map.GID( i ) << "] = "
                           << "x(" << x.firstLocalIndex() + i  << ")="
                           << x( x.firstLocalIndex() + i ) << "\n";

            v->set( v_map.GID( i ), x( x.firstLocalIndex() + i ) );
        }
    }
#endif

    template< int index, typename spaceT >
    static Epetra_MultiVector getComponent( spaceT const& Xh, Epetra_MultiVector const& sol )
    {
        Epetra_Map componentMap ( epetraMap( Xh->template functionSpace<index>()->map() ) );
        Epetra_Map globalMap ( epetraMap( Xh->map() ) );

        //DVLOG(2) << "Component map: " << componentMap << "\n";

        Epetra_MultiVector component( componentMap, 1 );

        int Length = component.MyLength();

        int shift = Xh->nDofStart( index );

        for ( int i=0; i < Length; i++ )
        {
            int compGlobalID = componentMap.GID( i );

            if ( compGlobalID >= 0 )
            {
                int compLocalID = componentMap.LID( compGlobalID );

                int localID = globalMap.LID( compGlobalID+shift );
                //                         int globalID = globalMap.GID(localID);

                DVLOG(2) << "[MyBackend] Copy entry sol[" << localID << "]=" <<  sol[0][localID]
                               << " to component[" << compLocalID << "]\n";

                component[0][compLocalID] = sol[0][localID];

                DVLOG(2) << component[0][compLocalID] << "\n";
            }
        }

        return component;
    }


    template< int index, typename spaceT >
    static void UpdateComponent( spaceT const& Xh, Epetra_MultiVector& sol, Epetra_MultiVector& comp )
    {
        Epetra_Map componentMap ( epetraMap( Xh->template functionSpace<index>()->map() ) );
        Epetra_Map globalMap ( epetraMap( Xh->map() ) );

        int shift = Xh->nDofStart( index );

        int Length = comp.MyLength();

        for ( int i=0; i < Length; i++ )
        {
            int compGlobalID = componentMap.GID( i );

            if ( compGlobalID >= 0 )
            {
                int compLocalID = componentMap.LID( compGlobalID );

                int localID = globalMap.LID( compGlobalID+shift );
                //                         int globalID = globalMap.GID(localID);

                DVLOG(2) << "Copy entry component[" << compLocalID << "] to sol[" << localID << "]="
                               <<  sol[0][localID]
                               << "]\n";

                sol[0][localID] = comp[0][compLocalID] ;

                DVLOG(2) << comp[0][compLocalID] << "\n";
            }
        }
    }



    // -- SETTING OF OPTIONS --
    void set_maxiter ( int maxiter );

    void set_tol ( double tol );

    void set_drop ( double drop );

    void set_overlap ( int overlap );

    void set_restart ( int restart );

    void set_fillin ( int fillin );

    void set_order ( int order );

    void set_overlap_type( std::string str );

    void set_residual ( std::string str );

    void set_localparts ( int localparts );

    void set_solver( std::string str );

    void set_prec_local_solver( std::string str );

    void set_prec_local_solver_type( std::string str );

    void set_prec( std::string str );

    void set_verbose( const std::string verb="none" );

    void set_options( list_type opts );

    list_type get_options();

    void set_symmetric( bool /*is_sym*/ ) {}


    // -- LINEAR ALGEBRA INTERFACE --
    static void
    applyMatrix( sparse_matrix_type const& A,
                 const Vector<value_type>& x,
                 vector_type& b );

    real_type dot( const vector_type& f, const vector_type& x ) const;

    void prod( sparse_matrix_type const& A,
               vector_type const& x,
               vector_type& b, bool transpose ) const
    {
        FEELPP_ASSERT( !transpose ).warn( "not implemented yet" );
        epetra_sparse_matrix_type const& _A = dynamic_cast<epetra_sparse_matrix_type const&>( A );
        epetra_vector_type const& _x = dynamic_cast<epetra_vector_type const&>( x );
        epetra_vector_type& _b = dynamic_cast<epetra_vector_type&>( b );
        _A.mat().Apply( _x.vec(), _b.vec() );
    }

    solve_return_type solve( base_sparse_matrix_ptrtype const& A,
                             base_sparse_matrix_ptrtype const& B,
                             base_vector_ptrtype& x,
                             base_vector_ptrtype const& b );

    bool converged()
    {
        return true;
    }

private:

    mpi::communicator M_comm;

    list_type M_options;

    std::string M_prec_type, M_local_solver, M_local_solver_type;

    operator_ptrtype M_Prec;

}; // class BackendTrilinos

/**
 * Create a program_options options description for trilinos backend
 *
 * \param prefix prefix string to differentiate trilinos solvers
 */
po::options_description backendtrilinos_options( std::string const& prefix );


} // Feel

#endif // FEELPP_HAS_TRILINOS_EPETRA


#endif /* _BACKENDTRILINOS_HPP_ */

