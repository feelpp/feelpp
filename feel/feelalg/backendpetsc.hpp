/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-25

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file backendpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-25
 */

#ifndef _BACKENDPETSC_HPP_
#define _BACKENDPETSC_HPP_

#include <boost/program_options/variables_map.hpp>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/solverlinearpetsc.hpp>
#include <feel/feelalg/solvernonlinearpetsc.hpp>
#include <feel/feelalg/backend.hpp>



namespace Feel
{
namespace po = boost::program_options;





#if defined( FEELPP_HAS_PETSC_H )
/**
 * \class BackendPetsc
 *
 * this class provides an interface to the PETSC linear algebra library
 */
template<typename T>
class BackendPetsc : public Backend<T>
{
    typedef Backend<T> super;
public:

    // -- TYPEDEFS --
    typedef typename super::value_type value_type;

    /* matrix */
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef MatrixPetsc<value_type> petsc_sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> petsc_sparse_matrix_ptrtype;
    typedef MatrixPetscMPI<value_type> petscMPI_sparse_matrix_type;
    //typedef boost::shared_ptr<sparse_matrix_type> petscMPI_sparse_matrix_ptrtype;

    typedef typename sparse_matrix_type::graph_type graph_type;
    typedef typename sparse_matrix_type::graph_ptrtype graph_ptrtype;

    /* vector */
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;
    typedef VectorPetsc<value_type> petsc_vector_type;
    typedef boost::shared_ptr<vector_type> petsc_vector_ptrtype;
    typedef VectorPetscMPI<value_type> petscMPI_vector_type;

    typedef typename super::solve_return_type solve_return_type;
    typedef typename super::nl_solve_return_type nl_solve_return_type;

    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;


    // -- CONSTRUCTOR --
    BackendPetsc( WorldComm const& worldComm=Environment::worldComm() )
        :
        super( worldComm ),
        M_solver_petsc( worldComm )
        //M_nl_solver_petsc( "",worldComm )
    {
        this->M_backend = BackendType::BACKEND_PETSC;
    }

    BackendPetsc( po::variables_map const& vm, std::string const& prefix = "",
                  WorldComm const& worldComm=Environment::worldComm() )
        :
        super( vm, prefix, worldComm ),
        M_solver_petsc( vm, worldComm )
        //M_nl_solver_petsc( prefix,worldComm )
    {
        this->M_backend = BackendType::BACKEND_PETSC;

        std::string _prefix = prefix;

        if ( !_prefix.empty() )
            _prefix += "-";
    }
    ~BackendPetsc();
    void clear();

    sparse_matrix_ptrtype 
    newMatrix() 
    {
        sparse_matrix_ptrtype mat;
        if ( this->comm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( this->comm() ) );
        else // seq
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( this->comm() ) );
        return mat;
    }

    // -- FACTORY METHODS --
    template<typename DomainSpace, typename DualImageSpace>
    static sparse_matrix_ptrtype newMatrix( DomainSpace const& Xh,
                                            DualImageSpace const& Yh,
                                            size_type matrix_properties = NON_HERMITIAN )
    {
        auto s = stencil( _test=Yh,_trial=Xh );

        sparse_matrix_ptrtype mat;
        if ( Yh->worldComm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( Yh->dof(),Xh->dof() ) );
        else // seq
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( Yh->dof(),Xh->dof() ) );

        mat->setMatrixProperties( matrix_properties );
        mat->init( Yh->nDof(), Xh->nDof(),
                   Yh->nLocalDofWithoutGhost(), Xh->nLocalDofWithoutGhost(),
                   s->graph() );
        //Yh->nLocalDof(), Xh->nLocalDof() );
#if 0
        auto nSpace = DomainSpace::nSpaces;

        std::vector < std::vector<int> > is( nSpace );
        uint cptSpaces=0;

        //boost::tuple< typename DomainSpace::functionspace_vector_type, uint, std::vector < std::vector<int> > > hola;
        //        auto result = boost::make_tuple(Xh->functionSpaces(),cptSpaces,is);
        auto result = boost::make_tuple( cptSpaces,is );
        boost::fusion::fold( Xh->functionSpaces(), result, computeNDofForEachSpace() );

        for ( uint i = 0; i<nSpace; i++ )
        {
            //is[i].resize()
        }

#endif

        return mat;
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
        sparse_matrix_ptrtype mat;

        if ( this->comm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type );
        else
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type );

        mat->setMatrixProperties( matrix_properties );
        mat->init( m,n,m_l,n_l,nnz,noz );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( datamap_ptrtype const& domainmap,
               datamap_ptrtype const& imagemap,
               size_type matrix_properties = NON_HERMITIAN,
               bool init = true )
    {
        sparse_matrix_ptrtype mat;

        if ( imagemap->worldComm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( imagemap,domainmap,imagemap->worldComm() ) );
        else
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( imagemap,domainmap,imagemap->worldComm() ) );

        mat->setMatrixProperties( matrix_properties );

        if ( init )
        {
            mat->init( imagemap->nDof(), domainmap->nDof(),
                       imagemap->nLocalDofWithoutGhost(), domainmap->nLocalDofWithoutGhost() );
        }

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
        sparse_matrix_ptrtype mat;
        auto const& mapGraphRow = graph->mapRowPtr();
        auto const& mapGraphCol = graph->mapColPtr();
        if ( this->comm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( mapGraphRow,mapGraphCol,mapGraphRow->worldComm() ) );
        else
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( mapGraphRow,mapGraphCol,mapGraphRow->worldComm() ) ) ;

        mat->setMatrixProperties( matrix_properties );
        //mat->init( m,n,m_l,n_l,graph );
        mat->init( mapGraphRow->nDof(), mapGraphCol->nDof(),
                   mapGraphRow->nLocalDofWithoutGhost(), mapGraphCol->nLocalDofWithoutGhost(),
                   graph );


        return mat;
    }

   sparse_matrix_ptrtype
   newZeroMatrix( datamap_ptrtype const& domainmap,
                  datamap_ptrtype const& imagemap )
    {
        graph_ptrtype sparsity_graph( new graph_type( imagemap, domainmap ) );
        sparsity_graph->zero();
        sparsity_graph->close();

        sparse_matrix_ptrtype mat;
        if ( imagemap->worldComm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( imagemap,domainmap,imagemap->worldComm() ) );
        else
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( imagemap,domainmap,imagemap->worldComm() ) );

        mat->init( imagemap->nDof(), domainmap->nDof(),
                   imagemap->nLocalDofWithoutGhost(), domainmap->nLocalDofWithoutGhost(),
                   sparsity_graph );

        return mat;
    }

    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l )
    {
        graph_ptrtype sparsity_graph( new graph_type( 0,0,m_l-1,0,n_l-1 ) );
        sparsity_graph->zero();
        sparsity_graph->close();

        sparse_matrix_ptrtype mat;

        if ( this->comm().globalSize()>1 )
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type );
        else
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type );

        //mat->setMatrixProperties( matrix_properties );
        mat->init( m, n, m_l, n_l, sparsity_graph );

        return mat;
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( SpaceT const& space )
    {
        if ( space->worldComm().globalSize()>1 )
            return vector_ptrtype( new petscMPI_vector_type( space->dof() ) );
        else
            return vector_ptrtype( new petsc_vector_type( space->dof() ) );
    }

    vector_ptrtype newVector( datamap_ptrtype const& dm )
    {
        if ( this->comm().globalSize()>1 ) return vector_ptrtype( new petscMPI_vector_type( dm ) );

        else return vector_ptrtype( new petsc_vector_type( dm ) );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local )
    {
        return vector_ptrtype( new petsc_vector_type( n, n_local, this->comm() ) );
    }


    void set_symmetric( bool /*is_sym*/ ) {}

    // -- LINEAR ALGEBRA INTERFACE --
    template <class Vector>
    static void applyMatrix( sparse_matrix_type const& A,
                             const Vector& x,
                             vector_type& b )
    {
        // petsc mat/vec here
    }

    void prod( sparse_matrix_type const& A,
               vector_type const& x,
               vector_type& b, bool transpose ) const
    {
        int ierr = 0;
        petsc_sparse_matrix_type const& _A = dynamic_cast<petsc_sparse_matrix_type const&>( A );
        petsc_vector_type const& _x = dynamic_cast<petsc_vector_type const&>( x );
        petsc_vector_type const& _b = dynamic_cast<petsc_vector_type const&>( b );
        if(!transpose) {
            if ( _A.mapCol().worldComm().globalSize() == x.map().worldComm().globalSize() )
            {
                //std::cout << "BackendPetsc::prod STANDART"<< std::endl;
                ierr = MatMult( _A.mat(), _x.vec(), _b.vec() );
                CHKERRABORT( _A.comm().globalComm(),ierr );
            }
            else
            {
                //std::cout << "BackendPetsc::prod with convert"<< std::endl;
                auto x_convert = petscMPI_vector_type(_A.mapColPtr());
                x_convert.duplicateFromOtherPartition(x);
                x_convert.close();
                ierr = MatMult( _A.mat(), x_convert.vec(), _b.vec() );
                CHKERRABORT( _A.comm().globalComm(),ierr );
            }
        }
        else {
            if ( _A.mapRow().worldComm().globalSize() == x.map().worldComm().globalSize() )
            {
                //std::cout << "BackendPetsc::prod STANDART"<< std::endl;
                ierr = MatMultTranspose( _A.mat(), _x.vec(), _b.vec() );
                CHKERRABORT( _A.comm().globalComm(),ierr );
            }
            else
            {
                //std::cout << "BackendPetsc::prod with convert"<< std::endl;
                auto x_convert = petscMPI_vector_type(_A.mapRowPtr());
                x_convert.duplicateFromOtherPartition(x);
                x_convert.close();
                ierr = MatMultTranspose( _A.mat(), x_convert.vec(), _b.vec() );
                CHKERRABORT( _A.comm().globalComm(),ierr );
            }
        }
        b.close();
    }

    /**
     * get the matrix \c M whose diagonal is \c -v
     */
    int diag( vector_type const& v, sparse_matrix_type& M ) const;

    /**
     * @return the vector \c v with diagonal of \c M
     */
    int diag( sparse_matrix_type const& M, vector_type& v ) const;

    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             vector_type const& b );

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& B,
                             vector_ptrtype& x,
                             vector_ptrtype const& b );

    /**
     * assemble \f$C=P^T A P\f$
     */
    int PtAP( sparse_matrix_ptrtype const& A,
              sparse_matrix_ptrtype const& P,
              sparse_matrix_ptrtype& C ) const;
    /**
     * assemble \f$C=P A P^T\f$
     */
    int PAPt( sparse_matrix_ptrtype const& A,
              sparse_matrix_ptrtype const& P,
              sparse_matrix_ptrtype& C ) const;
    
    template <class Vector>
    static value_type dot( const vector_type& f,
                           const Vector& x )
    {
        value_type result( 0 );

        // petsc dot here

        return result;
    }

    /**
     * @return the linear solver (const version)
     */
    SolverLinearPetsc<double> const& linearSolver() const { return M_solver_petsc; }
    /**
     * @return the linear solver 
     */
    SolverLinearPetsc<double> & linearSolver() { return M_solver_petsc; }
    
private:

    SolverLinearPetsc<double> M_solver_petsc;
    //SolverNonLinearPetsc<double> M_nl_solver_petsc;

}; // class BackendPetsc




#if !defined(FEELPP_BACKEND_PETSC_NOEXTERN)
extern template class BackendPetsc<double>;
extern template class BackendPetsc<std::complex<double>>;
#endif

#endif // FEELPP_HAS_PETSC_H
} // Feel

#endif /* _BACKENDPETSC_HPP_ */
