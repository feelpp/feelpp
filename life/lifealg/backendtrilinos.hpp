/* -*- mode: c++ -*-

   This file is part of the Life library

   Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
              Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2007-01-18

   Copyright (C) 2007 EPFL
   Copyright (C) 2008 Université Joseph Fourier (Grenoble 1)


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
/**
   \file backendtrilinos.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-08-14
*/
#ifndef _BACKENDTRILINOS_HPP_
#define _BACKENDTRILINOS_HPP_


#include <boost/program_options/variables_map.hpp>
#include <lifeconfig.h>

#if defined ( HAVE_TRILINOS_EPETRA )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#include <life/lifealg/matrixepetra.hpp>
#include <life/lifealg/vectorepetra.hpp>
#include <life/lifealg/preconditionerifpack.hpp>
#include <life/lifealg/preconditionerml.hpp>
#include <life/lifealg/solverlineartrilinos.hpp>
#include <life/lifealg/operatortrilinos.hpp>
#include <life/lifealg/datamap.hpp>
#include <life/lifealg/backend.hpp>


#include <Teuchos_ParameterList.hpp>

namespace Life
{
namespace po = boost::program_options;

// forward declararion
template<class T> class Backend;

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
    BackendTrilinos()
        :
        super(),
        M_options(),
        M_prec_type( "" ),
        M_Prec()
    {
        set_maxiter( 1000 );
        set_tol( 1e-10 );
    }


    BackendTrilinos( po::variables_map const& vm, std::string const& prefix = "" );

    BackendTrilinos( const BackendTrilinos& tc );


    // -- FACTORY METHODS --
    static Epetra_Map epetraMap( DataMap const& dmap )
    {
        std::vector<int> e( dmap.nMyElements() );
        std::copy( dmap.myGlobalElements().begin(),
                   dmap.myGlobalElements().end(),
                   e.begin() );
        return Epetra_Map( -1, dmap.nMyElements(), e.data(), 0, Epetra_MpiComm( Application::comm() ) );
    }


    static Epetra_Map epetraMapStatic( DataMap const& dmap )
    {
        return Epetra_Map( dmap.nGlobalElements(), dmap.nMyElements(), 0, Epetra_MpiComm( Application::comm() ) );
    }


    template<typename DomainSpace, typename DualImageSpace>
    sparse_matrix_ptrtype newMatrix( DomainSpace const& Xh,
                                     DualImageSpace const& Yh )
    {
        return newMatrix( Xh->map(), Yh->map() );
    }
    sparse_matrix_ptrtype newMatrix( DataMap const& domainmap,
                                     DataMap const& imagemap )
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
        		return sparse_matrix_ptrtype( new epetra_sparse_matrix_type( erowmap, ecolmap, eimagemap, edomainmap ) );
		}
	else
		return sparse_matrix_ptrtype( new epetra_sparse_matrix_type( erowmap, ecolmap ) );
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


    static operator_ptrtype IfpackPrec( sparse_matrix_ptrtype const& M, list_type options )
    {
        PreconditionerIfpack P( options );

        P.buildPreconditioner( M );

        return P.getPrec();
    }


    static operator_ptrtype MLPrec( sparse_matrix_ptrtype& M, list_type options )
    {
        PreconditionerML P( options );

        P.buildPreconditioner( M );

        return P.getPrec();
    }

    template< typename element_type >
    static void Epetra2Ublas( vector_ptrtype const& u, element_type& x )
    {
        epetra_vector_ptrtype const& _v( dynamic_cast<epetra_vector_ptrtype const&>( u ) );
        Epetra_Map v_map( _v->Map() );

        vector_type v = *u;

        //Debug(10003) << "Initial EpetraVector " << v << "\n";

        const size_type L = v.localSize();

        for ( size_type i=0; i<L; i++)
            {
                Debug(10003) << "x(" << x.firstLocalIndex() + i  << ")="
                             << "v[" << v_map.GID(i) << "] = "
                             << v(i) << "\n";

                x( x.firstLocalIndex() + i ) = v(i);
            }

        Debug(10003) << "Epetra2Ublas:" << x << "\n";
    }


    template< typename element_type >
    static void Ublas2Epetra( element_type const& x, vector_ptrtype& v )
    {
        epetra_vector_type& _v( dynamic_cast<epetra_vector_type&>( *v ) );
        Epetra_Map v_map( _v.Map() );

        Debug(10002) << "Local size of ublas vector" << x.localSize() << "\n";
        Debug(10002) << "Local size of epetra vector" << v->localSize() << "\n";

        const size_type L = v->localSize();

        for ( size_type i=0; i<L; i++ )
            {
                Debug(10002) << "v[" << v_map.GID(i) << "] = "
                             << "x(" << x.firstLocalIndex() + i  << ")="
                             << x( x.firstLocalIndex() + i ) << "\n";

                v->set(v_map.GID(i), x( x.firstLocalIndex() + i ) );
            }
    }


    template< int index, typename spaceT >
    static Epetra_MultiVector getComponent( spaceT const& Xh, Epetra_MultiVector const& sol )
    {
        Epetra_Map componentMap ( epetraMap( Xh->template functionSpace<index>()->map() ) );
        Epetra_Map globalMap ( epetraMap( Xh->map() ) );

        //Debug(10006) << "Component map: " << componentMap << "\n";

        Epetra_MultiVector component( componentMap, 1 );

        int Length = component.MyLength();

        int shift = Xh->nDofStart(index);

        for ( int i=0; i < Length; i++)
            {
                int compGlobalID = componentMap.GID(i);

                if ( compGlobalID >= 0 )
                    {
                        int compLocalID = componentMap.LID(compGlobalID);

                        int localID = globalMap.LID(compGlobalID+shift);
//                         int globalID = globalMap.GID(localID);

                        Debug(10006) << "[MyBackend] Copy entry sol[" << localID << "]=" <<  sol[0][localID]
                                     << " to component[" << compLocalID << "]\n";

                        component[0][compLocalID] = sol[0][localID];

                        Debug(10006) << component[0][compLocalID] << "\n";
                    }
            }

        return component;
    }


    template< int index, typename spaceT >
    static void UpdateComponent( spaceT const& Xh, Epetra_MultiVector& sol, Epetra_MultiVector& comp )
    {
        Epetra_Map componentMap ( epetraMap( Xh->template functionSpace<index>()->map() ) );
        Epetra_Map globalMap ( epetraMap( Xh->map() ) );

        int shift = Xh->nDofStart(index);

        int Length = comp.MyLength();

        for ( int i=0; i < Length; i++)
            {
                int compGlobalID = componentMap.GID(i);

                if ( compGlobalID >= 0 )
                    {
                        int compLocalID = componentMap.LID(compGlobalID);

                        int localID = globalMap.LID(compGlobalID+shift);
//                         int globalID = globalMap.GID(localID);

                        Debug(10006) << "Copy entry component[" << compLocalID << "] to sol[" << localID << "]="
                                     <<  sol[0][localID]
                                     << "]\n";

                        sol[0][localID] = comp[0][compLocalID] ;

                        Debug(10006) << comp[0][compLocalID] << "\n";
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
               vector_type& b ) const
    {
        epetra_sparse_matrix_type const& _A = dynamic_cast<epetra_sparse_matrix_type const&>( A );
        epetra_vector_type const& _x = dynamic_cast<epetra_vector_type const&>( x );
        epetra_vector_type& _b = dynamic_cast<epetra_vector_type&>( b );
        _A.mat().Apply( _x.vec(), _b.vec() );
    }

    solve_return_type solve( base_sparse_matrix_ptrtype const& A,
                             base_sparse_matrix_ptrtype const& B,
                             base_vector_ptrtype& x,
                             base_vector_ptrtype const& b );

    bool converged() { return true; }

private:

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


} // Life

#endif // HAVE_TRILINOS_EPETRA


#endif /* _BACKENDTRILINOS_HPP_ */

