//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 23 Sep 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_BILINEARFORMBASE_HPP
#define FEELPP_BILINEARFORMBASE_HPP 1

#include <future>


#include <feel/feelconfig.h>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelvf/block.hpp>




namespace Feel
{

//!
//! Base class for Bilinear Forms
//! handle algebraic representation and multithreading
//!
template<typename T=double>
class BilinearFormBase : public CommObject
{
public:
    using super = CommObject;
    using list_block_type = Feel::vf::list_block_type;
    using value_type = T;
    using pre_solve_type = typename Backend<value_type>::pre_solve_type;
    using post_solve_type = typename Backend<value_type>::post_solve_type;

    //typedef ublas::compressed_matrix<value_type, ublas::row_major> csr_matrix_type;
    typedef MatrixSparse<value_type> matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;
    static const bool is_row_major = true;//matrix_type::is_row_major;

    using size_type =  typename matrix_type::size_type;
    
    typedef typename mpl::if_<mpl::equal_to<mpl::bool_<is_row_major>, mpl::bool_<true> >,
                              mpl::identity<ublas::row_major>,
                              mpl::identity<ublas::column_major> >::type::type layout_type;


    BilinearFormBase() = default;
    
    template<typename FE1,  typename FE2>
    BilinearFormBase( std::string name,
                      FE1 const& Xh,
                      FE2 const& Yh,
                      matrix_ptrtype& __M,
                      size_type rowstart = 0,
                      size_type colstart = 0,
                      bool build = true,
                      bool do_threshold = false,
                      value_type threshold = type_traits<value_type>::epsilon(),
                      size_type graph_hints = Pattern::COUPLED );

    template<typename FE1,  typename FE2>
    BilinearFormBase( std::string name,
                      FE1 const& Xh,
                      FE2 const& Yh,
                      matrix_ptrtype& __M,
                      list_block_type const& __lb = {},
                      size_type rowstart = 0,
                      size_type colstart = 0,
                      bool do_threshold = false,
                      value_type threshold = type_traits<value_type>::epsilon(),
                      size_type graph_hints = Pattern::COUPLED );

    BilinearFormBase( BilinearFormBase const& __vf );
    BilinearFormBase( BilinearFormBase && __vf ) = default;
    ~BilinearFormBase() override
        {
            //toc(M_name, FLAGS_v > 0 );
        }

    /**
     * copy operator
     */
    BilinearFormBase&
    operator=( BilinearFormBase const& form );

    BilinearFormBase& operator+=( BilinearFormBase& a )
        {
            if ( this == &a )
            {
                M_matrix->scale( 2.0 );
                return *this;
            }
            M_matrix->addMatrix( 1.0, a.M_matrix );
            return *this;
        }

    BilinearFormBase& add( double alpha, BilinearFormBase&  a )
        {
            M_matrix->addMatrix( alpha, a.M_matrix );
            return *this;
        }
    
    virtual void push_back( std::future<void>&& f ) { M_fut_assign.push_back( std::forward<std::future<void>>( f ) ); }
    virtual void get()
        {
            //std::cout << "-- before fut: " << M_fut_assign.size() << std::endl;
            for( auto& f : M_fut_assign )
                f.get();
            //M_fut_assign.clear();
            //std::cout << "-- after fut: " << M_fut_assign.size() << std::endl;
        }

    /** @name Accessors
     */
    //@{

    //!
    //! @return the name of the bilinear form
    //!
    std::string const& name() const { return M_name; }
    
    /**
     * return the pattern
     */
    size_type pattern() const
    {
        return M_pattern;
    }

    /**
     * \return true if the pattern is coupled with respect to the components,
     * false otherwise
     */
    bool isPatternCoupled() const
    {
        Feel::Context ctx( M_pattern );
        return ctx.test( Pattern::COUPLED );
    }

    /**
     * \return true if the pattern is the default one, false otherwise
     */
    bool isPatternDefault() const
    {
        Feel::Context ctx( M_pattern );
        return ctx.test( Pattern::DEFAULT );
    }

    /**
     * \return true if the pattern adds the neighboring elements, false otherwise
     */
    bool isPatternNeighbor() const
    {
        Feel::Context ctx( M_pattern );
        return ctx.test( Pattern::EXTENDED );
    }
    bool isPatternExtended() const
    {
        Feel::Context ctx( M_pattern );
        return ctx.test( Pattern::EXTENDED );
    }

    bool isPatternSymmetric() const
    {
        Feel::Context ctx( M_pattern );
        return ctx.test( Pattern::PATTERN_SYMMETRIC );
    }

    /**
     * \return the matrix associated to the bilinear form
     */
    matrix_type const& matrix() const
    {
        return *M_matrix;
    }

    matrix_type& matrix()
    {
        return *M_matrix;
    }

    matrix_ptrtype const& matrixPtr() const
    {
        return M_matrix;
    }

    matrix_ptrtype& matrixPtr()
    {
        return M_matrix;
    }

    list_block_type const& blockList() const
    {
        return M_lb;
    }

    size_type rowStartInMatrix() const
    {
        return M_row_startInMatrix;
    }

    size_type colStartInMatrix() const
    {
        return M_col_startInMatrix;
    }
    //!
    //! @return number of non-zero entries
    //!
    std::size_t nnz() const
    {
        return M_matrix->nnz();
    }
    /**
     * @brief set the bilinear form to zero
     * @details set the bilinear form and its
     * algebraic representation to zero
     */
    void zero()
    {
        M_matrix->zero();
    }
    /**
     * \return the threshold
     */
    value_type threshold() const
    {
        return M_threshold;
    }

    /**
     * \return \c true if threshold applies, false otherwise
     */
    bool doThreshold( value_type const& v ) const
    {
        return ( math::abs( v ) > M_threshold );
    }

    /**
     * return true if do threshold. false otherwise
     */
    bool doThreshold() const
    {
        return M_do_threshold;
    }

    /**
     * \return the mapping from test dof id to container id with global process numbering
     */
    std::vector<size_type> const& dofIdToContainerIdTest() const { return *M_dofIdToContainerIdTest; }
    /**
     * \return the mapping from trial dof id to container id with global process numbering
     */
    std::vector<size_type> const& dofIdToContainerIdTrial() const { return *M_dofIdToContainerIdTrial; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set a threshold value for the matrix entries associated with the
     * bilinear form
     */
    void setThreshold( value_type eps )
    {
        M_threshold = eps;
    }

    /**
     * set the threshold strategy, true threshold the matrix entries,
     * false do not threshold
     */
    void setDoThreshold( bool do_threshold )
    {
        M_do_threshold = do_threshold;
    }
    /**
     * set mapping from test dof id to container id with global process numbering
     */
    void setDofIdToContainerIdTest( std::vector<size_type> const& gpmap ) { M_dofIdToContainerIdTest = std::addressof( gpmap ); }
    /**
     * set mapping from trial dof id to container id with global process numbering
     */
    void setDofIdToContainerIdTrial( std::vector<size_type> const& gpmap ) { M_dofIdToContainerIdTrial = std::addressof( gpmap ); }

    //@}

    /** @name  Methods
     */
    //@{


    // close matrix
    void close()
        {
            this->get(); // futures
            M_matrix->close();
        }

    bool closed() const noexcept
        {
            return M_matrix->closed();
        }

    /**
     * Diagonalize representation(matrix) associated to the \p
     * BilinearForm at selected dofs \p dofs by putting 0 on
     * the extra diagonal terms and 1 on the diagonal.
     *
     * If \p ON_ELIMINATION_KEEP_DIAGONAL is set in \p on_context then
     * the diagonal value of the matrix is kept and the right habd
     * side \p rhs is modified accordingly.
     */
    void zeroRows( std::vector<int> const& __dofs,
                   Vector<value_type> const& __values,
                   Vector<value_type>& rhs,
                   Feel::Context const& on_context,
                   double value_on_diagonal );

    /**
     * add value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void add( size_type i,  size_type j,  value_type const& v )
    {
        if ( M_do_threshold )
        {
            if ( doThreshold( v ) )
                M_matrix->add( i+this->rowStartInMatrix(),
                                j+this->colStartInMatrix(),
                                v );
        }

        else
            M_matrix->add( i+this->rowStartInMatrix(),
                            j+this->colStartInMatrix(),
                            v );

    }
    /**
     * add value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void addMatrix( int* rows, int nrows,
                    int* cols, int ncols,
                    value_type* data,
                    size_type K  = 0,
                    size_type K2 = invalid_v<size_type> );


    /**
     * set value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void set( size_type i,  size_type j,  value_type const& v )
    {
        M_matrix->set( i, j, v );
    }

    void addToNOz( size_type i, size_type n ) 
    {
        M_n_oz[i] += n;
    }
    void addToNNz( size_type i, size_type n )
    {
        M_n_nz[i] += n;
    }
    size_type nOz( size_type i ) const
    {
        return M_n_oz[i];
    }
    size_type nNz( size_type i ) const
    {
        return M_n_nz[i];
    }

    template<typename X1, typename X2>
    void allocateMatrix( std::shared_ptr<X1> const& x1, std::shared_ptr<X2> const& x2 )
    {
        M_matrix = backend()->newMatrix( _test=x1, _trial=x2 );
    }
    bool isMatrixAllocated() const
    {
        return (bool)M_matrix;
    }
    BOOST_PARAMETER_MEMBER_FUNCTION( ( typename Backend<value_type>::solve_return_type ),
                                     solve,
                                     tag,
                                     ( required
                                       ( in_out( solution ),* )
                                       ( rhs, * ) )
                                     ( optional
                                       ( name,           ( std::string ), "" )
                                       ( kind,           ( std::string ), soption(_prefix=name,_name="backend") )
                                       ( rebuild,        ( bool ), boption(_prefix=name,_name="backend.rebuild") )
                                       ( pre, (pre_solve_type), pre_solve_type() )
                                       ( post, (post_solve_type), post_solve_type() )
                                         ) )
        {
            this->close();
            return Feel::backend( _name=name, _kind=kind, _rebuild=rebuild,
                                  _worldcomm=this->worldCommPtr() )->solve( _matrix=this->matrixPtr(),
                                                            _rhs=rhs.vectorPtr(),
                                                            _solution=solution,
                                                            _pre=pre,
                                                            _post=post
                                                            );
        }

    BOOST_PARAMETER_MEMBER_FUNCTION( ( typename Backend<value_type>::solve_return_type ),
                                     solveb,
                                     tag,
                                     ( required
                                       ( in_out( solution ),* )
                                       ( rhs, * )
                                       ( backend, *) )
                                     ( optional
                                       ( prec,           ( preconditioner_ptrtype ),
                                         preconditioner( _prefix=backend->prefix(),
                                                         _matrix=this->matrixPtr(),
                                                         _pc=backend->pcEnumType()/*LU_PRECOND*/,
                                                         _pcfactormatsolverpackage=backend->matSolverPackageEnumType(),
                                                         _backend=backend ) )
                                         ) )
        {
            this->close();
            return backend->solve( _matrix=this->matrixPtr(), _rhs=rhs.vectorPtr(),
                                   _solution=solution, _prec = prec );
        }

    //@}

protected:
    std::string M_name;
    size_type M_pattern;

    matrix_ptrtype M_matrix;

    list_block_type M_lb;
    size_type M_row_startInMatrix,M_col_startInMatrix;

    bool M_do_build;
    bool M_do_threshold;
    value_type M_threshold;

    std::vector<size_type> M_n_nz;
    std::vector<size_type> M_n_oz;

    std::vector<size_type> const* M_dofIdToContainerIdTest;
    std::vector<size_type> const* M_dofIdToContainerIdTrial;

    std::vector<std::future<void>> M_fut_assign;
    std::mutex b_mutex;
    
};

template<typename T>
template<typename FE1,  typename FE2>
BilinearFormBase<T>::BilinearFormBase( std::string name,
                                       FE1 const& Xh,
                                       FE2 const& Yh,
                                       matrix_ptrtype& __M,
                                       size_type rowstart,
                                       size_type colstart,
                                       bool build,
                                       bool do_threshold,
                                       value_type threshold,
                                       size_type graph_hints )
:
    super( Xh->worldCommPtr() ),
    M_pattern( graph_hints ),
    M_matrix( __M ),
    M_lb{},
    M_row_startInMatrix( rowstart ),
    M_col_startInMatrix( colstart ),
    M_do_build( build ),
    M_do_threshold( do_threshold ),
    M_threshold( threshold )
{
    boost::timer tim;
    DVLOG(2) << "begin constructor with default listblock\n";

    if ( !Xh->worldComm().isActive() ) return;

    if ( !M_matrix ) M_matrix = backend()->newMatrix( _test=Xh, _trial=Yh );

    M_lb.push_back( Feel::vf::Block ( 0, 0, 0, 0 ) );
    auto dmTest = M_matrix->mapRowPtr();
    auto dmTrial = M_matrix->mapColPtr();
    this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
    this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );

    DVLOG(2) << " - form init in " << tim.elapsed() << "\n";
    DVLOG(2) << "begin constructor with default listblock done\n";
}

template<typename T>
template<typename FE1,  typename FE2>
BilinearFormBase<T>::BilinearFormBase( std::string name,
                                       FE1 const& Xh,
                                       FE2 const& Yh,
                                       matrix_ptrtype& __M,
                                       list_block_type const& __lb,
                                       size_type rowstart,
                                       size_type colstart,
                                       bool do_threshold,
                                       value_type threshold,
                                       size_type graph_hints )
:
    super( Xh->worldCommPtr() ),
    M_name( name ),
    M_pattern( graph_hints ),
    M_matrix( __M ),
    M_lb( __lb ),
    M_row_startInMatrix( rowstart ),
    M_col_startInMatrix( colstart ),
    M_do_build( false ),
    M_do_threshold( do_threshold ),
    M_threshold( threshold )
{
    if ( !Xh->worldComm().isActive() ) return;

    if ( !M_matrix ) M_matrix = backend()->newMatrix( _test=Xh, _trial=Yh );
    auto dmTest = M_matrix->mapRowPtr();
    auto dmTrial = M_matrix->mapColPtr();
    this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
    this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );
}


#if !defined(FEELPP_BILINEARFORMBASE_NOEXTERN)
extern template class BilinearFormBase<double>;
//extern template class Backend<std::complex<double>>;
#endif


} // namespace Feel



#endif
