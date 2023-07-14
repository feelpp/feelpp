/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
template<typename T, typename SizeT = uint32_type>
class FEELPP_EXPORT BackendPetsc : public Backend<T,SizeT>
{
    typedef Backend<T> super;
public:

    // -- TYPEDEFS --
    typedef typename super::value_type value_type;
    using size_type = typename super::size_type;
    /* matrix */
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef MatrixPetsc<value_type> petsc_sparse_matrix_type;
    typedef std::shared_ptr<sparse_matrix_type> petsc_sparse_matrix_ptrtype;
    typedef MatrixPetscMPI<value_type> petscMPI_sparse_matrix_type;
    //typedef std::shared_ptr<sparse_matrix_type> petscMPI_sparse_matrix_ptrtype;

    typedef typename sparse_matrix_type::graph_type graph_type;
    typedef typename sparse_matrix_type::graph_ptrtype graph_ptrtype;

    /* vector */
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;
    typedef VectorPetsc<value_type> petsc_vector_type;
    typedef std::shared_ptr<vector_type> petsc_vector_ptrtype;
    typedef VectorPetscMPI<value_type> petscMPI_vector_type;
    typedef std::shared_ptr<petscMPI_vector_type> petscMPI_vector_ptrtype;

    typedef typename super::solve_return_type solve_return_type;
    typedef typename super::nl_solve_return_type nl_solve_return_type;

    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;


    // -- CONSTRUCTOR --
    explicit BackendPetsc( worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() )
        :
        super( worldComm ),
        M_solver_petsc( worldComm )
        //M_nl_solver_petsc( "",worldComm )
    {
        this->M_backend = BackendType::BACKEND_PETSC;
    }

    BackendPetsc( po::variables_map const& vm, std::string const& prefix = "",
                  worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() )
        :
        super( vm, prefix, worldComm ),
        M_solver_petsc( prefix, worldComm, vm )
        //M_nl_solver_petsc( prefix,worldComm )
    {
        this->M_backend = BackendType::BACKEND_PETSC;
    }
    ~BackendPetsc() override;
    void clear() override;


    /**
     * convert a vector into a backend pointer vector
     */
    vector_ptrtype toBackendVectorPtr( vector_type const& v  ) override;

    /**
     * convert a pointer vector into a backend pointer vector
     */
    vector_ptrtype toBackendVectorPtr( vector_ptrtype const& v  ) override;


    sparse_matrix_ptrtype
    newMatrix() override
    {
        sparse_matrix_ptrtype mat;
        if ( this->comm().globalSize()>1 )
            mat = std::make_shared<petscMPI_sparse_matrix_type>( this->worldCommPtr() );
        else // seq
            mat = std::make_shared<petsc_sparse_matrix_type>( this->worldCommPtr() );
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
            mat = std::make_shared<petscMPI_sparse_matrix_type>( Yh->dof(),Xh->dof() );
        else // seq
            mat = std::make_shared<petsc_sparse_matrix_type>( Yh->dof(),Xh->dof() );

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
               size_type matrix_properties = NON_HERMITIAN ) override
    {
        sparse_matrix_ptrtype mat;

        if ( this->comm().globalSize()>1 )
            mat = std::make_shared<petscMPI_sparse_matrix_type>();
        else
            mat = std::make_shared<petsc_sparse_matrix_type>();

        mat->setMatrixProperties( matrix_properties );
        mat->init( m,n,m_l,n_l,nnz,noz );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( datamap_ptrtype const& domainmap,
               datamap_ptrtype const& imagemap,
               size_type matrix_properties = NON_HERMITIAN,
               bool init = true ) override
    {
        sparse_matrix_ptrtype mat;

        if ( imagemap->worldComm().globalSize()>1 )
            mat = std::make_shared<petscMPI_sparse_matrix_type>( imagemap,domainmap,imagemap->worldCommPtr() );
        else
            mat = std::make_shared<petsc_sparse_matrix_type>( imagemap,domainmap,imagemap->worldCommPtr() );

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
               size_type matrix_properties = NON_HERMITIAN ) override
    {
        sparse_matrix_ptrtype mat;
        auto const& mapGraphRow = graph->mapRowPtr();
        auto const& mapGraphCol = graph->mapColPtr();
        if ( this->comm().globalSize()>1 )
            mat = std::make_shared<petscMPI_sparse_matrix_type>( mapGraphRow,mapGraphCol,mapGraphRow->worldCommPtr() );
        else
            mat = std::make_shared<petsc_sparse_matrix_type>( mapGraphRow,mapGraphCol,mapGraphRow->worldCommPtr() );

        mat->setMatrixProperties( matrix_properties );
        //mat->init( m,n,m_l,n_l,graph );
        mat->init( mapGraphRow->nDof(), mapGraphCol->nDof(),
                   mapGraphRow->nLocalDofWithoutGhost(), mapGraphCol->nLocalDofWithoutGhost(),
                   graph );


        return mat;
    }

   sparse_matrix_ptrtype
   newZeroMatrix( datamap_ptrtype const& domainmap,
                  datamap_ptrtype const& imagemap ) override
    {
        graph_ptrtype sparsity_graph( new graph_type( imagemap, domainmap ) );
        sparsity_graph->zero();
        sparsity_graph->close();

        sparse_matrix_ptrtype mat;
        if ( imagemap->worldComm().globalSize()>1 )
            mat = std::make_shared<petscMPI_sparse_matrix_type>( imagemap,domainmap,imagemap->worldCommPtr() );
        else
            mat = std::make_shared<petsc_sparse_matrix_type>( imagemap,domainmap,imagemap->worldCommPtr() );

        mat->init( imagemap->nDof(), domainmap->nDof(),
                   imagemap->nLocalDofWithoutGhost(), domainmap->nLocalDofWithoutGhost(),
                   sparsity_graph );

        return mat;
    }

    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l ) override
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

    sparse_matrix_ptrtype newIdentityMatrix( datamap_ptrtype const& domainmap, datamap_ptrtype const& imagemap ) override
    {
        graph_ptrtype sparsity_graph( new graph_type( imagemap,domainmap ) );
        sparsity_graph->addMissingZeroEntriesDiagonal();
        sparsity_graph->close();
        sparse_matrix_ptrtype mat = this->newMatrix(0,0,0,0,sparsity_graph);
        auto vecDiag = this->newVector( imagemap );
        vecDiag->setConstant( 1. );
        mat->setDiagonal( vecDiag );
        return mat;
    }

    template <typename SpaceT>
    static auto newVector( SpaceT const& space )
    {
        using underlying_t = typename std::remove_reference<decltype(*space)>::type;
        static_assert(std::is_member_function_pointer<decltype(&underlying_t::dof)>::value, "Underlying type of SpaceT must have a dof() member function");

        return newVector(space->dof());
    }

    /**
     * @brief create a new vector of \c N VectorPetsc
     *
     * @tparam SpaceT function space type
     * @param space function space
     * @param N number of vectors
     * @return std::vector<vector_ptrtype>
     */
    template <typename SpaceT, typename = std::enable_if_t<is_functionspace_v<decay_type<SpaceT>>>>
    static auto newVectors( SpaceT const& space, size_t vector_count ) 
    {
        using underlying_t = typename std::remove_reference<decltype(*space)>::type;
        static_assert(std::is_member_function_pointer<decltype(&underlying_t::dof)>::value, "Underlying type of SpaceT must have a dof() member function");

        return newVectors(space->dof());
    }

    /**
     * @brief create a new vector of \c N VectorPetsc
     * 
     * @param dm Datamap
     * @param vector_count size of the vector
     * @return vector_ptrtype 
     */
    std::vector<vector_ptrtype> newVectors( datamap_ptrtype const& dm, const size_type vector_count ) override
    {
        std::vector<vector_ptrtype> result;

        // Reserve memory in advance for efficiency
        result.reserve(vector_count);

        try
        {
            for ( size_t i = 0; i < vector_count; ++i )
            {
                result.push_back(newVector(dm));
            }     
        }
        catch(const std::bad_alloc& e) // catches if new allocation fails
        {
            LOG(WARNING) << fmt::format("[BackendPetsc::newVectors] Memory allocation failed: {} ",e.what()) << '\n';
            throw;
        }
        return result;
    }

    /**
     * @brief create a new VectorPetsc
     * 
     * @param dm Datamap
     * @return vector_ptrtype 
     */
    vector_ptrtype newVector( datamap_ptrtype const& dm ) override
    {
        if ( this->comm().globalSize()>1 ) 
            return std::make_shared<petscMPI_vector_type>( dm );
        else 
           return std::make_shared<petsc_vector_type>( dm );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local ) override
    {
        return std::make_shared<petsc_vector_type>( n, n_local, this->worldCommPtr() );
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

    /**
     * apply matrix vector product with the matrix A and vector x
     * and stores the result in b.
     */
    void prod( sparse_matrix_type const& A,
               vector_type const& x,
               vector_type& b, bool transpose ) const override;

    /**
     * get the matrix \c M whose diagonal is \c -v
     */
    int diag( vector_type const& v, sparse_matrix_type& M ) const override;

    /**
     * @return the vector \c v with diagonal of \c M
     */
    int diag( sparse_matrix_type const& M, vector_type& v ) const override;

    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             vector_type const& b );

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& B,
                             vector_ptrtype& x,
                             vector_ptrtype const& b ) override;

    /**
     * assemble \f$C=P^T A P\f$
     */
    int PtAP( sparse_matrix_ptrtype const& A,
              sparse_matrix_ptrtype const& P,
              sparse_matrix_ptrtype& C ) const override;
    /**
     * assemble \f$C=P A P^T\f$
     */
    int PAPt( sparse_matrix_ptrtype const& A,
              sparse_matrix_ptrtype const& P,
              sparse_matrix_ptrtype& C ) const override;

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
extern template class BackendPetsc<double, uint32_type>;
extern template class BackendPetsc<std::complex<double>,uint32_type>;
#endif

#endif // FEELPP_HAS_PETSC_H
} // Feel

#include <feel/feelalg/topetsc.hpp>

#endif /* _BACKENDPETSC_HPP_ */
