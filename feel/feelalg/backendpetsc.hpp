/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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

    // -- CONSTRUCTOR --
    BackendPetsc( WorldComm const& worldComm=Environment::worldComm() )
        :
        super( worldComm ),
        M_solver_petsc( worldComm ),
        M_nl_solver_petsc( worldComm )
    {}

    BackendPetsc( po::variables_map const& vm, std::string const& prefix = "",
                  WorldComm const& worldComm=Environment::worldComm() )
        :
        super( vm, prefix, worldComm ),
        M_solver_petsc( vm, worldComm ),
        M_nl_solver_petsc( worldComm )
    {
        std::string _prefix = prefix;

        if ( !_prefix.empty() )
            _prefix += "-";
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
            mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( Yh->map(),Xh->map() ) );
        else // seq
            mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( Yh->map(),Xh->map() ) );

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

        if ( this->comm().globalSize()>1 ) mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type );

        else mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type );

        mat->setMatrixProperties( matrix_properties );
        mat->init( m,n,m_l,n_l,nnz,noz );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( DataMap const& domainmap,
               DataMap const& imagemap,
               size_type matrix_properties = NON_HERMITIAN,
               bool init = true )
    {
        sparse_matrix_ptrtype mat;

        if ( imagemap.worldComm().globalSize()>1 ) mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( imagemap,domainmap,imagemap.worldComm() ) );

        else mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( imagemap,domainmap,imagemap.worldComm() ) );

        mat->setMatrixProperties( matrix_properties );

        if ( init )
        {
            mat->init( imagemap.nDof(), domainmap.nDof(),
                       imagemap.nLocalDofWithoutGhost(), domainmap.nLocalDofWithoutGhost() );
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

        if ( this->comm().globalSize()>1 ) mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type );

        else mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type );

        mat->setMatrixProperties( matrix_properties );
        mat->init( m,n,m_l,n_l,graph );
        return mat;
    }

    sparse_matrix_ptrtype
    newZeroMatrix( DataMap const& domainmap,
                   DataMap const& imagemap )
    {
        graph_ptrtype sparsity_graph( new graph_type( 0,
                                      0, imagemap.nLocalDofWithoutGhost()-1,
                                      0, domainmap.nLocalDofWithoutGhost()-1,
                                      imagemap.worldComm() ) );
        sparsity_graph->zero();
        sparsity_graph->close();

        sparse_matrix_ptrtype mat;

        if ( imagemap.worldComm().globalSize()>1 ) mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type( imagemap,domainmap,imagemap.worldComm() ) );

        else mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type( imagemap,domainmap,imagemap.worldComm() ) );

        //mat->setMatrixProperties( matrix_properties );
        mat->init( imagemap.nDof(), domainmap.nDof(),
                   imagemap.nLocalDofWithoutGhost(), domainmap.nLocalDofWithoutGhost(),
                   sparsity_graph );

        return mat;
        //return newZeroMatrix(imagemap.nDof(), domainmap.nDof(),
        //                     imagemap.nLocalDofWithoutGhost(), domainmap.nLocalDofWithoutGhost() );
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

        if ( this->comm().globalSize()>1 ) mat = sparse_matrix_ptrtype( new petscMPI_sparse_matrix_type );

        else mat = sparse_matrix_ptrtype( new petsc_sparse_matrix_type );

        //mat->setMatrixProperties( matrix_properties );
        mat->init( m, n, m_l, n_l, sparsity_graph );

        return mat;
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( SpaceT const& space )
    {
        if ( space->worldComm().globalSize()>1 ) return vector_ptrtype( new petscMPI_vector_type( space->map() ) );

        else return vector_ptrtype( new petsc_vector_type( space->map() ) );

        //return this->newVector(space->map());
        //return vector_ptrtype( new vector_type( space->nDof(), space->nLocalDof() ) );
    }

    vector_ptrtype newVector( DataMap const& dm )
    {
        if ( this->comm().globalSize()>1 ) return vector_ptrtype( new petscMPI_vector_type( dm ) );

        else return vector_ptrtype( new petsc_vector_type( dm ) );

        //return vector_ptrtype( new petsc_vector_type( dm.nGlobalElements(), dm.nMyElements() ) );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local )
    {
        return vector_ptrtype( new petsc_vector_type( n, n_local ) );
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
               vector_type& b ) const
    {
        int ierr = 0;
        petsc_sparse_matrix_type const& _A = dynamic_cast<petsc_sparse_matrix_type const&>( A );
        petsc_vector_type const& _x = dynamic_cast<petsc_vector_type const&>( x );
        petsc_vector_type const& _b = dynamic_cast<petsc_vector_type const&>( b );
        if ( _A.mapCol().worldComm().globalSize() == x.map().worldComm().globalSize() )
        {
            //std::cout << "BackendPetsc::prod STANDART"<< std::endl;
            ierr = MatMult( _A.mat(), _x.vec(), _b.vec() );
            CHKERRABORT( _A.comm().globalComm(),ierr );
        }
        else
        {
            //std::cout << "BackendPetsc::prod with convert"<< std::endl;
            auto x_convert = petscMPI_vector_type(_A.mapCol());
            x_convert.duplicateFromOtherPartition(x);
            x_convert.close();
            ierr = MatMult( _A.mat(), x_convert.vec(), _b.vec() );
            CHKERRABORT( _A.comm().globalComm(),ierr );
        }
        b.close();
    }


    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             vector_type const& b );

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& B,
                             vector_ptrtype& x,
                             vector_ptrtype const& b );

    template <class Vector>
    static value_type dot( const vector_type& f,
                           const Vector& x )
    {
        value_type result( 0 );

        // petsc dot here

        return result;
    }

private:

    SolverLinearPetsc<double> M_solver_petsc;
    SolverNonLinearPetsc<double> M_nl_solver_petsc;

}; // class BackendPetsc


template<typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_ptrtype const& A,
                        sparse_matrix_ptrtype const& B,
                        vector_ptrtype& x,
                        vector_ptrtype const& b )
{
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if (!M_solver_petsc.initialized())
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance=this->rTolerance(),
                                  _atolerance=this->aTolerance(),
                                  _dtolerance=this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    auto res = M_solver_petsc.solve( *A, *B, *x, *b, this->rTolerance(), this->maxIterations(), this->transpose() );
    Debug( 7005 ) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    Debug( 7005 ) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.get<0>() ) std::cerr<< "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve


template<typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_type const& A,
                        vector_type& x,
                        vector_type const& b )
{
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if (!M_solver_petsc.initialized())
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance=this->rTolerance(),
                                  _atolerance=this->aTolerance(),
                                  _dtolerance=this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    auto res = M_solver_petsc.solve( A, x, b, this->rTolerance(), this->maxIterations() );
    Debug( 7005 ) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    Debug( 7005 ) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.get<0>() ) std::cerr<< "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve

po::options_description backendpetsc_options( std::string const& prefix = "" );

#endif // FEELPP_HAS_PETSC_H
} // Feel

#endif /* _BACKENDPETSC_HPP_ */
