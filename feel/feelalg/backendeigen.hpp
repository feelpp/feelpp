/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-06-20

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file backendeigen.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-06-20
 */

#ifndef _BACKENDEIGEN_HPP_
#define _BACKENDEIGEN_HPP_

#include <boost/make_shared.hpp>
#include <boost/program_options/variables_map.hpp>
#include <feel/feelcore/feelpetsc.hpp>
#undef MatType
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feelalg/matrixeigensparse.hpp>
#include <feel/feelalg/vectoreigen.hpp>



#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <feel/feelalg/backend.hpp>
namespace Feel
{
namespace po = boost::program_options;


/**
 * \class BackendEigen
 *
 * this class provides an interface to the EIGEN linear algebra library
 */
template<typename T, int _Options = 0>
class BackendEigen : public Backend<T>
{
    typedef Backend<T> super;
public:

    // -- TYPEDEFS --
    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;

    static const bool IsDense = (_Options == 1);
    static const bool IsSparse = (_Options == 0);

    /* matrix */
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename mpl::if_<mpl::bool_<IsDense>,
                              mpl::identity<MatrixEigenDense<value_type> >,
                              mpl::identity<MatrixEigenSparse<value_type> > >::type::type eigen_sparse_matrix_type;
    typedef boost::shared_ptr<eigen_sparse_matrix_type> eigen_sparse_matrix_ptrtype;

    typedef typename sparse_matrix_type::graph_type graph_type;
    typedef typename sparse_matrix_type::graph_ptrtype graph_ptrtype;

    /* vector */
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;
    typedef VectorEigen<value_type> eigen_vector_type;
    typedef boost::shared_ptr<vector_type> eigen_vector_ptrtype;

    typedef typename super::solve_return_type solve_return_type;
    typedef typename super::nl_solve_return_type nl_solve_return_type;

    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;

    // -- CONSTRUCTOR --
    BackendEigen();
    BackendEigen( WorldComm const& worldComm=Environment::worldComm() );
    BackendEigen( po::variables_map const& vm, std::string const& prefix = "",
                  WorldComm const& worldComm=Environment::worldComm() );

    // -- FACTORY METHODS --
    sparse_matrix_ptrtype
    newMatrix()
    {
        auto A= boost::make_shared<eigen_sparse_matrix_type>(0,0,this->comm());
        return A;
    }
    template<typename DomainSpace, typename DualImageSpace>
    static sparse_matrix_ptrtype newMatrix( boost::shared_ptr<DomainSpace> const& space1,
                                            boost::shared_ptr<DualImageSpace> const& space2,
                                            size_type matrix_properties = NON_HERMITIAN )
    {
        Context ctx( matrix_properties );
        //if ( ctx.test( DENSE ) )
        {
            //auto A= sparse_matrix_ptrtype( new eigen_sparse_matrix_type( space1->nDof(), space2->nDof() ) );
            //auto A= sparse_matrix_ptrtype( new eigen_sparse_matrix_type( space1->nDof(), space2->nDof(), this->comm() ) );
            auto A= boost::make_shared<eigen_sparse_matrix_type>( space1->nDof(), space2->nDof(), space1->map()->worldComm() );
            A->setMatrixProperties( matrix_properties );
            return A;
        }
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
        //sparse_matrix_ptrtype mat( new eigen_sparse_matrix_type( m,n ) );
        sparse_matrix_ptrtype mat( new eigen_sparse_matrix_type( m,n,this->comm() ) );
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
        //sparse_matrix_ptrtype mat( new eigen_sparse_matrix_type( m,n ) );
        sparse_matrix_ptrtype mat( new eigen_sparse_matrix_type( m,n,this->comm() ) );
        mat->setMatrixProperties( matrix_properties );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( datamap_ptrtype const& d1, datamap_ptrtype const& d2, size_type matrix_properties = NON_HERMITIAN, bool init = true )
    {
        auto A = sparse_matrix_ptrtype( new eigen_sparse_matrix_type( d1->nGlobalElements(),
                                                                      d2->nGlobalElements(),
                                                                      this->comm() ) );
        A->setMatrixProperties( matrix_properties );
        return A;
    }

    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l )
    {
        auto A = sparse_matrix_ptrtype( new eigen_sparse_matrix_type( m, n ) );
        //A->setMatrixProperties( matrix_properties );
        return A;
    }


    sparse_matrix_ptrtype
    newZeroMatrix( datamap_ptrtype const& d1, datamap_ptrtype const& d2 )
    {
        auto A = sparse_matrix_ptrtype( new eigen_sparse_matrix_type( d1->nGlobalElements(), d2->nGlobalElements() ) );
        //A->setMatrixProperties( matrix_properties );
        return A;
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( boost::shared_ptr<SpaceT> const& space )
    {
        //return vector_ptrtype( new eigen_vector_type( space->nDof() ) );
        return vector_ptrtype( new eigen_vector_type( space->nDof(), space->map()->comm() ) );
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( SpaceT const& space )
    {
        //return vector_ptrtype( new eigen_vector_type( space.nDof() ) );
        return vector_ptrtype( new eigen_vector_type( space.nDof(), space.map()->comm() ) );
    }

    vector_ptrtype
    newVector( datamap_ptrtype const& d )
    {
        return vector_ptrtype( new eigen_vector_type( d->nGlobalElements(), this->comm() ) );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local )
    {
        return vector_ptrtype( new eigen_vector_type( n, this->comm() ) );
    }


    // -- LINEAR ALGEBRA INTERFACE --
    void prod( sparse_matrix_type const& A,
               vector_type const& x,
               vector_type& b, bool transpose ) const
    {
        eigen_sparse_matrix_type const& _A = dynamic_cast<eigen_sparse_matrix_type const&>( A );
        eigen_vector_type const& _x = dynamic_cast<eigen_vector_type const&>( x );
        eigen_vector_type& _b = dynamic_cast<eigen_vector_type&>( b );
        if(!transpose)
            _b.vec() = _A.mat()*_x.vec();
        else
            _b.vec() = _A.mat().adjoint()*_x.vec();
    }

    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             const vector_type& b )
        {
            return this->solve( A, x, b, mpl::bool_<IsDense>() );
        }

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             vector_ptrtype& x,
                             const vector_ptrtype& b )
    {
        return this->solve( *A, *x, *b );
    }

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& P,
                             vector_ptrtype& x,
                             const vector_ptrtype& b )
    {
        return this->solve( *A, *x, *b );
    }


    value_type dot( const vector_type& f,
                   const vector_type& x ) const
    {
        eigen_vector_type const& _f = dynamic_cast<eigen_vector_type const&>( f );
        eigen_vector_type const& _x = dynamic_cast<eigen_vector_type const&>( x );
        return _f.vec().dot( _x.vec() );
    }

private:
    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             const vector_type& b,
                             mpl::bool_<true> );
    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             const vector_type& b,
                             mpl::bool_<false> );

}; // class BackendEigen

// -- CONSTRUCTOR --
template<typename T, int _Options>
BackendEigen<T,_Options>::BackendEigen( WorldComm const& _worldComm )
    :
    super(_worldComm)
{
    if ( IsSparse )
        this->M_backend = BackendType::BACKEND_EIGEN;
    else
        this->M_backend = BackendType::BACKEND_EIGEN_DENSE;
}

template<typename T, int _Options>
BackendEigen<T,_Options>::BackendEigen( po::variables_map const& vm, std::string const& prefix, WorldComm const& _worldComm  )
    :
    super( vm, prefix, _worldComm )
{
    if ( IsSparse )
        this->M_backend = BackendType::BACKEND_EIGEN;
    else
        this->M_backend = BackendType::BACKEND_EIGEN_DENSE;

    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";
}



template<typename T, int _Options>
typename BackendEigen<T,_Options>::solve_return_type
BackendEigen<T,_Options>::solve( sparse_matrix_type const& _A,
                                 vector_type& _x,
                                 const vector_type& _b,
                                 mpl::bool_<true>)
{
    bool reusePC = ( this->precMatrixStructure() == SAME_PRECONDITIONER );

    eigen_sparse_matrix_type const& A( dynamic_cast<eigen_sparse_matrix_type const&>( _A ) );
    eigen_vector_type      & x( dynamic_cast<eigen_vector_type      &>( _x ) );
    eigen_vector_type const& b( dynamic_cast<eigen_vector_type const&>( _b ) );
    x.vec() = A.mat().lu().solve(b.vec());

    return solve_return_type( boost::make_tuple(true,1,1e-10) );
} // BackendEigen::solve

template<typename T, int _Options>
typename BackendEigen<T,_Options>::solve_return_type
BackendEigen<T,_Options>::solve( sparse_matrix_type const& _A,
                                 vector_type& _x,
                                 const vector_type& _b,
                                 mpl::bool_<false>)
{
    bool reusePC = ( this->precMatrixStructure() == SAME_PRECONDITIONER );

    eigen_sparse_matrix_type const& A( dynamic_cast<eigen_sparse_matrix_type const&>( _A ) );
    eigen_vector_type      & x( dynamic_cast<eigen_vector_type      &>( _x ) );
    eigen_vector_type const& b( dynamic_cast<eigen_vector_type const&>( _b ) );

    //x.vec()=A.mat().template fullPivLu().solve(b.vec());
    Eigen::SimplicialLDLT<typename eigen_sparse_matrix_type::matrix_type> solver;
    solver.compute(A.mat());
    x.vec() = solver.solve(b.vec());

    // if(solver.info()!=Eigen::Succeeded) {
    //     // solving failed
    //     return boost::make_tuple(false,1,1e-10);;
    // }
    return solve_return_type( boost::make_tuple(true,1,1e-10) );
} // BackendEigen::solve




po::options_description backendeigen_options( std::string const& prefix = "" );

} // Feel

#endif /* _BACKENDEIGEN_HPP_ */
