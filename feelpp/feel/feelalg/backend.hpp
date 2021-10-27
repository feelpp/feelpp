/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-12-23

  Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-present Feel++ Consortium

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
   \file backend.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-12-23
 */
#ifndef Backend_H
#define Backend_H 1

#include <functional>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/fusion/include/fold.hpp>
#include <boost/smart_ptr/enable_shared_from_this.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/matrixblock.hpp>
//#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/datamap.hpp>

#include <feel/feelalg/solverlinear.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feeldiscr/functionspacebase.hpp>

#include <feel/feelalg/matrixshell.hpp>
#include <feel/feelalg/matrixshellsparse.hpp>



//#include <feel/feelvf/vf.hpp>
//#include <boost/fusion/support/pair.hpp>
//#include <boost/fusion/container.hpp>
//#include <boost/fusion/sequence.hpp>
//#include <boost/fusion/algorithm.hpp>

//namespace fusion = boost::fusion;

//#include <feel/feelvf/bilinearform.hpp>
#include <feel/feelvf/pattern.hpp>
#include <feel/feelvf/block.hpp>

#include <feel/feelalg/nullspace.hpp>

namespace Feel
{
/*enum  Pattern
{
DEFAULT   = 1 << 0,
EXTENDED  = 1 << 1,
COUPLED   = 1 << 2,
SYMMETRIC = 1 << 3
};
*/
///! \cond detail
namespace detail
{
template<typename T>
FEELPP_EXPORT std::shared_ptr<DataMap<>> datamap( T const& t, mpl::true_ )
{
    return t->mapPtr();
}
template<typename T>
FEELPP_EXPORT std::shared_ptr<DataMap<>> datamap( T const& t, mpl::false_ )
{
    return t.mapPtr();
}
template<typename T>
FEELPP_EXPORT std::shared_ptr<DataMap<>> datamap( T const& t )
{
    return datamap( t, Feel::detail::is_shared_ptr<T>() );
}

template<typename T>
#if BOOST_VERSION >= 105300
FEELPP_EXPORT typename boost::detail::sp_dereference< typename T::element_type >::type
#else
FEELPP_EXPORT typename T::reference
#endif
ref( T t, mpl::true_ )
{
    return *t;
}
template<typename T>
FEELPP_EXPORT T& ref( T& t, mpl::false_ )
{
    return t;
}
template<typename T>
FEELPP_EXPORT auto ref( T& t ) -> decltype( ref( t, Feel::detail::is_shared_ptr<T>() ) )
{
    return ref( t, Feel::detail::is_shared_ptr<T>() );
}


}
///! \endcond detail

/**
 * default pre/post solve function, a no-op
 */
inline void default_prepost_solve( vector_ptrtype rhs, vector_ptrtype sol )
{

}

template<typename T> class MatrixBlockBase;
template<int NR, int NC, typename T> class MatrixBlock;
template<typename T, typename SizeT> class VectorBlockBase;
template<int NR, typename T, typename SizeT> class VectorBlock;

template<typename T> class BlocksBaseSparseMatrix;
template<typename T, typename SizeT> class BlocksBaseVector;

class BackendBase : public CommObject
{
public:
    using super = CommObject;
    BackendBase( worldcomm_ptr_t const& w ) : super( w ) {}
    BackendBase( BackendBase const& ) = default;
    BackendBase( BackendBase && ) = default;
    ~BackendBase() override = default;
};
/**
 * \class Backend
 * \brief base class for all linear algebra backends
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename T, typename SizeT = uint32_type>
class FEELPP_EXPORT Backend : public BackendBase, public std::enable_shared_from_this<Backend<T,SizeT>>
{
public:


    /** @name Typedefs
     */
    //@{
    using super = BackendBase;
    typedef T value_type;
    using size_type = SizeT;
    typedef typename type_traits<T>::real_type real_type;

    typedef Vector<value_type> vector_type;
    typedef std::shared_ptr<vector_type> vector_ptrtype;
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef std::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef MatrixShell<value_type> shell_matrix_type;
    typedef std::shared_ptr<shell_matrix_type> shell_matrix_ptrtype;

    typedef typename sparse_matrix_type::graph_type graph_type;
    typedef typename sparse_matrix_type::graph_ptrtype graph_ptrtype;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_ptrtype ptrtype;

    typedef SolverNonLinear<value_type,size_type> solvernonlinear_type;
    typedef std::shared_ptr<solvernonlinear_type> solvernonlinear_ptrtype;

    typedef typename SolverLinear<real_type,size_type>::solve_return_type solve_return_type;
    typedef typename solvernonlinear_type::solve_return_type nl_solve_return_type;

    typedef DataMap<size_type> datamap_type;
    typedef std::shared_ptr<datamap_type> datamap_ptrtype;

    typedef typename datamap_type::indexsplit_type indexsplit_type;
    typedef typename datamap_type::indexsplit_ptrtype indexsplit_ptrtype;

    using pre_solve_type = std::function<void(vector_ptrtype,vector_ptrtype)>;
    using post_solve_type = std::function<void(vector_ptrtype,vector_ptrtype)>;
    using update_nlsolve_type = typename solvernonlinear_type::update_iteration_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Backend( worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );
    Backend( po::variables_map const& vm, std::string const& prefix = "", worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );
    Backend( Backend const& ) = default;
    Backend( Backend && ) = default;
    ~Backend() override;


    /**
     * Builds a \p Backend, if Petsc is available, use Petsc by
     * default, otherwise use GMM which is distributed with feel
     */
    static FEELPP_DEPRECATED backend_ptrtype build(
#if defined( FEELPP_HAS_PETSC_H )
        BackendType = BACKEND_PETSC
#else
        BackendType = BACKEND_GMM
#endif
        , worldcomm_ptr_t const& worldComm=Environment::worldCommPtr()
    );

    /**
     * build a backend
     * \param kind the type of backend
     * \param prefix the name of the backend
     * \param worldcomm the communicator
     * \return the backend
     */
    static backend_ptrtype build( std::string const& kind, std::string const& prefix = "", worldcomm_ptr_t const& worldComm=Environment::worldCommPtr(), po::variables_map const& vm = Environment::vm() );


    static FEELPP_DEPRECATED backend_ptrtype build( po::variables_map const& vm, std::string const& prefix = "", worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );

    /**
     * build a backend
     * \param bt the type of backend
     * \param prefix the name of the backend
     * \param worldcomm the communicator
     * \return the backend
     */
    static FEELPP_DEPRECATED backend_ptrtype build( BackendType bt, std::string const& prefix = "", worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );

    /**
     * convert a vector into a backend pointer vector
     */
    virtual vector_ptrtype toBackendVectorPtr( vector_type const& v  )
    {
        vector_ptrtype _newvec;
        return _newvec;
    }

    /**
     * convert a pointer vector into a backend pointer vector
     */
    virtual vector_ptrtype toBackendVectorPtr( vector_ptrtype const& v  )
    {
        vector_ptrtype _newvec;
        return _newvec;
    }
    /**
     * convert a pointer vector into a backend pointer vector
     * apply a dynamic_pointer_cast is necessary for shared<VectorUblas>
     */
    template<typename VecType>
    vector_ptrtype toBackendVectorPtr( std::shared_ptr<VecType> const& v  )
    {
        vector_ptrtype vcast = std::dynamic_pointer_cast< vector_type >( v );
        return this->toBackendVectorPtr( vcast );
    }

    /**
     * instantiate a new sparse matrix
     */
    virtual sparse_matrix_ptrtype newMatrix() = 0;

    /**
     * instantiate a new sparse vector
     */
    virtual sparse_matrix_ptrtype newMatrix( const size_type m,
            const size_type n,
            const size_type m_l,
            const size_type n_l,
            const size_type nnz=30,
            const size_type noz=10,
            size_type prop = NON_HERMITIAN ) = 0;

    /**
     * instantiate a new sparse vector
     */
    virtual sparse_matrix_ptrtype newMatrix( const size_type m,
            const size_type n,
            const size_type m_l,
            const size_type n_l,
            graph_ptrtype const & graph,
            size_type matrix_properties = NON_HERMITIAN ) = 0;

    /**
     * instantiate a new sparse vector
     */
    sparse_matrix_ptrtype newMatrix( const size_type m,
                                     const size_type n,
                                     const size_type m_l,
                                     const size_type n_l,
                                     graph_ptrtype const & graph,
                                     indexsplit_ptrtype const& indexSplit,
                                     size_type matrix_properties = NON_HERMITIAN )
    {
        auto mat = this->newMatrix( m,n,m_l,n_l,graph,matrix_properties );
        mat->setIndexSplit( indexSplit );
        return mat;
    }


    /**
     * instantiate a new sparse vector
     */
    virtual sparse_matrix_ptrtype newMatrix( datamap_ptrtype const& dm1,
                                             datamap_ptrtype const& dm2,
                                             size_type prop = NON_HERMITIAN,
                                             bool init = true ) = 0;

    /**
     * instantiate a new sparse vector
     */
    sparse_matrix_ptrtype newMatrix( datamap_ptrtype const& domainmap,
                                     datamap_ptrtype const& imagemap,
                                     graph_ptrtype const & graph,
                                     size_type matrix_properties = NON_HERMITIAN,
                                     bool init = true )
    {
        auto mat = this->newMatrix( domainmap,imagemap, matrix_properties, false );

        if ( init ) mat->init( imagemap->nDof(), domainmap->nDof(),
                                   imagemap->nLocalDofWithoutGhost(), domainmap->nLocalDofWithoutGhost(),
                                   graph );

        mat->zero();
        // todo!
        //mat->setIndexSplit( trial->dofIndexSplit() );
        return mat;
    }


    /**
     * instantiate a new sparse vector
     */
    virtual
    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l ) =0;

    virtual sparse_matrix_ptrtype newZeroMatrix( datamap_ptrtype const& dm1, datamap_ptrtype const& dm2 ) = 0;

    virtual sparse_matrix_ptrtype newIdentityMatrix( datamap_ptrtype const& dm1, datamap_ptrtype const& dm2 )
    {
        CHECK( false ) << "Not Implemented in base class!";
        return sparse_matrix_ptrtype{};
    }

    /**
     * helper function
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( sparse_matrix_ptrtype ),
                                     newMatrix,
                                     tag,
                                     ( required
                                       ( trial,*( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> > ) )
                                       ( test,*( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> > ) ) )
                                     ( optional
                                       ( pattern,( size_type ),Pattern::COUPLED )
                                       ( properties,( size_type ),NON_HERMITIAN )
                                       ( buildGraphWithTranspose, ( bool ),false )
                                       ( pattern_block,    *, ( BlocksStencilPattern(1,1,size_type( Pattern::HAS_NO_BLOCK_PATTERN ) ) ) )
                                       ( diag_is_nonzero,  *( boost::is_integral<mpl::_> ), true )
                                       ( verbose,   ( bool ), M_verbose )
                                       ( collect_garbage, *( boost::is_integral<mpl::_> ), true )
                                     ) )

#endif
        template <typename TrialType,typename TestType>
        sparse_matrix_ptrtype newMatrix( NA::arguments<
                                         typename na::trial::template required_as_t<std::shared_ptr<TrialType>>,
                                         typename na::test::template required_as_t<std::shared_ptr<TestType>>,
                                         typename na::pattern::template required_as_t<size_type>,
                                         typename na::properties::template required_as_t<size_type>,
                                         typename na::buildGraphWithTranspose::template required_as_t<bool>,
                                         typename na::pattern_block::template required_as_t<BlocksStencilPattern const&>,
                                         typename na::diag_is_nonzero::template required_as_t<bool>,
                                         typename na::verbose::template required_as_t<bool>,
                                         typename na::collect_garbage::template required_as_t<bool>
                                         > && args )
    {
        auto && trial = args.get(_trial);
        auto && test = args.get(_test);
        size_type pattern = args.get(_pattern);
        size_type properties = args.get(_properties);
        bool buildGraphWithTranspose = args.get(_buildGraphWithTranspose);
        BlocksStencilPattern const& pattern_block = args.get(_pattern_block);
        bool diag_is_nonzero = args.get(_diag_is_nonzero);
        bool verbose = args.get(_verbose);
        bool collect_garbage = args.get(_collect_garbage);

        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::newMatrix begin" );
        }

        if ( !this->comm().isActive() ) return sparse_matrix_ptrtype();

        //auto mat = this->newMatrix( trial->map(), test->map(), properties, false );
        auto mat = this->newMatrix( trial->dofOnOff(), test->dofOn(), properties, false );

        if ( this->type() == BackendType::BACKEND_EIGEN_DENSE ||
             this->type() == BackendType::BACKEND_EIGEN )
        {
            mat->init( test->nDof(), trial->nDof(),
                       test->nLocalDofWithoutGhost(), trial->nLocalDofWithoutGhost() );
        }
        else
        {
            if ( !buildGraphWithTranspose )
            {
                tic();
                auto s = stencil( _test=test,
                                  _trial=trial,
                                  _pattern=pattern,
                                  _pattern_block=pattern_block,
                                  _diag_is_nonzero=diag_is_nonzero,
                                  _collect_garbage=collect_garbage);
                toc( "Backend::newMatrix:: build stencil", FLAGS_v > 0 );
                tic();
                mat->init( test->nDof(), trial->nDof(),
                           test->nLocalDofWithoutGhost(), trial->nLocalDofWithoutGhost(),
                           s->graph() );
                toc( "Backend::newMatrix:: initialize matrix", FLAGS_v > 0 );
            }
            else
            {
                auto s = stencil( _test=trial,
                                  _trial=test,
                                  _pattern=pattern,
                                  _pattern_block=pattern_block.transpose(),
                                  _diag_is_nonzero=false,// because transpose(do just after)
                                  _close=false,
                                  _collect_garbage=collect_garbage );
                // get the good graph
                auto graph = s->graph()->transpose(false);
                if ( diag_is_nonzero )
                    graph->addMissingZeroEntriesDiagonal();
                graph->close();

                //maybe do that
                //stencilManagerGarbage(boost::make_tuple( trial, test, pattern, pattern_block.transpose().getSetOfBlocks(), false/*diag_is_nonzero*/));
                //now save the good graph in the StencilManager(only if entry is empty)
                stencilManagerAdd(boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero), graph);

                mat->init( test->nDof(), trial->nDof(),
                           test->nLocalDofWithoutGhost(), trial->nLocalDofWithoutGhost(),
                           graph );
            }
            tic();
            mat->zero();
            mat->setIndexSplit( trial->dofIndexSplit() );
            toc("Backend::newMatrix:: zero out matrix + set split", FLAGS_v > 0 );
        }

        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::newMatrix end" );
        }
        return mat;
    }

    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
        sparse_matrix_ptrtype newMatrix( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... )
            .add_default_arguments( NA::make_default_argument( _pattern, Pattern::COUPLED ),
                                    NA::make_default_argument( _properties, NON_HERMITIAN ),
                                    NA::make_default_argument( _buildGraphWithTranspose, false ),
                                    NA::make_default_argument( _pattern_block, BlocksStencilPattern(1,1,size_type( Pattern::HAS_NO_BLOCK_PATTERN )) ),
                                    NA::make_default_argument( _diag_is_nonzero, true ),
                                    NA::make_default_argument( _verbose, M_verbose ),
                                    NA::make_default_argument( _collect_garbage, true )
                                    );


        using trial_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype( args.get(_trial) )>>>;
        using test_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype( args.get(_test) )>>>;

        return newMatrix<trial_type,test_type>( std::move( args ) );
    }


    template<typename DomainSpace, typename ImageSpace>
    sparse_matrix_ptrtype newMatrix( DomainSpace const& dm, ImageSpace const& im, sparse_matrix_ptrtype const& M, size_type prop = NON_HERMITIAN  )
    {
        sparse_matrix_ptrtype m = newMatrix( dm, im, prop  );
        m->init( im->nDof(), dm->nDof(), im->nLocalDof(), dm->nLocalDof(), M->graph() );
        return m;
    }

    //!
    //! create a new block Matrix from a set of blocks
    //!
    sparse_matrix_ptrtype newBlockMatrixImpl( BlocksBaseSparseMatrix<value_type> const & b,
                                              bool copy_values=true,
                                              bool diag_is_nonzero=true )
    {
        typedef MatrixBlockBase<value_type> matrix_block_type;
        typedef std::shared_ptr<matrix_block_type> matrix_block_ptrtype;

        matrix_block_ptrtype mb;
        if ( b.isClosed() )
        {
            mb.reset( new matrix_block_type( b, *this, copy_values, diag_is_nonzero ) );
        }
        else
        {
            BlocksBaseSparseMatrix<value_type> copyBlock( b );
            copyBlock.close();
            mb.reset( new matrix_block_type( copyBlock, *this, copy_values, diag_is_nonzero ) );
        }
        return mb->getSparseMatrix();
    }

    //!
    //! create a new block Matrix from a set of graph blocks
    //!
    sparse_matrix_ptrtype newBlockMatrixImpl( BlocksBaseGraphCSR const & b,
                                              bool copy_values=true,
                                              bool diag_is_nonzero=true )
    {
        typedef MatrixBlockBase<value_type> matrix_block_type;
        typedef std::shared_ptr<matrix_block_type> matrix_block_ptrtype;

        matrix_block_ptrtype mb;
        if ( b.isClosed() )
        {
            mb.reset( new matrix_block_type( b, *this, diag_is_nonzero ) );
        }
        else
        {
            BlocksBaseGraphCSR copyBlock( b );
            copyBlock.close();
            mb.reset( new matrix_block_type( copyBlock, *this, diag_is_nonzero ) );
        }
        return mb->getSparseMatrix();
    }

    /**
     * instantiate a new block matrix sparse
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( sparse_matrix_ptrtype ),
                                     newBlockMatrix,
                                     tag,
                                     ( required
                                       ( block,* )
                                     )
                                     ( optional
                                       ( copy_values,*( boost::is_integral<mpl::_> ),true )
                                       ( diag_is_nonzero,  *( boost::is_integral<mpl::_> ), true )
                                     )
                                   )
    {
        return newBlockMatrixImpl( block,copy_values,diag_is_nonzero );
    }
#endif
    template <typename ... Ts>
    sparse_matrix_ptrtype newBlockMatrix( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && block = args.get(_block );
        bool copy_values = args.get_else( _copy_values, true );
        bool diag_is_nonzero = args.get_else( _diag_is_nonzero, true );
        return newBlockMatrixImpl( block,copy_values,diag_is_nonzero );
    }

    /**
     * instantiate a new block matrix sparse
     */
    template<typename TB>
        vector_ptrtype newBlockVectorImpl( BlocksBaseVector<TB,size_type> const & b,
                                           bool copy_values=true )
    {
        using vector_block_type = VectorBlockBase<TB,size_type>;
        return std::make_shared<vector_block_type>( b, *this, copy_values )->getVector();
    }

    /**
     * instantiate a new block matrix sparse
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( vector_ptrtype ),
                                     newBlockVector,
                                     tag,
                                     ( required
                                       ( block,* )
                                     )
                                     ( optional
                                       ( copy_values,*( boost::is_integral<mpl::_> ),true )
                                     )
                                   )
    {
        return newBlockVectorImpl( block,copy_values );
    }
#endif
    template <typename ... Ts>
    vector_ptrtype newBlockVector( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && block = args.get(_block );
        bool copy_values = args.get_else( _copy_values, true );
        return newBlockVectorImpl( block,copy_values );
    }


    /**
     * instantiate a new zero matrix
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( sparse_matrix_ptrtype ),
                                     newZeroMatrix,
                                     tag,
                                     ( required
                                       ( test,*( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> >) )
                                       ( trial,*( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> >) )
                                     )
                                   )
    {
        return this->newZeroMatrix( trial->dofOnOff(), test->dofOn() );
    }
#endif
    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) == 2 && ( NA::is_named_argument_v<Ts> && ...) > >
    sparse_matrix_ptrtype newZeroMatrix( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && test = args.get(_test );
        auto && trial = args.get(_trial );
        return this->newZeroMatrix( trial->dofOnOff(), test->dofOn() );
    }

    /**
     * instantiate a new vector
     */
    virtual vector_ptrtype newVector( datamap_ptrtype const& dm ) = 0;

    /**
     * instantiate a new vector
     */
    virtual vector_ptrtype newVector( const size_type n, const size_type n_local ) = 0;

    /**
     * helper function
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( vector_ptrtype ),
                                     newVector,
                                     tag,
                                     ( required
                                       ( test,*( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> >) )
                                     )
                                   )
    {
        if ( !this->comm().isActive() ) return vector_ptrtype();

        return this->newVector( test->dof() );
    }
#endif
    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
    vector_ptrtype newVector( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && test = args.get(_test );
        if ( !this->comm().isActive() )
            return vector_ptrtype();
        return this->newVector( test->dof() );
    }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return prefix of backend
     */
    std::string prefix() const
    {
        return M_prefix;
    }

    /**
     * \return the type of linear solver
     */
    std::string kspType() const
    {
        return M_ksp;
    }

    /**
     * \return the type of non linear solver
     */
    std::string snesType() const
    {
        return snesTypeConvertEnumToStr( this->nlSolver()->getType() );
    }

    /**
     * \return the type of preconditioner
     */
    std::string pcType() const
    {
        return M_pc;
    }

    /**
     * return true if the null space is the constant values, false otherwise
     */
    bool hasConstantNullSpace() const
    {
        return M_constant_null_space;
    }

    /**
     * \return the type of fieldsplitType
     */
    std::string fieldsplitType() const
    {
        return M_fieldSplit;
    }

    /**
     * \return enum pc type from options
     **/
    PreconditionerType pcEnumType() const;

    /**
     * \return enum solver type from options
     **/
    SolverType kspEnumType() const;

    /**
     * \return enum snes solver type from string
     */
    SolverNonLinearType snesEnumType() const;

    /**
     * \return enum fieldsplit type from options
     **/
    FieldSplitType fieldSplitEnumType() const;

    /**
     * \return the type of pcFactorMatSolverPackageType
     */
    std::string pcFactorMatSolverPackageType() const
    {
        return M_pcFactorMatSolverPackage;
    }

    /**
     * \return enum MatSolverPackage type from options
     **/
    MatSolverPackageType matSolverPackageEnumType() const;

    /**
     * \return the type of preconditioner associated to the matrix
     */
    MatrixStructure precMatrixStructure() const
    {
        return M_prec_matrix_structure;
    }

    /**
     * \return the relative tolerance
     */
    real_type rTolerance() const
    {
        return M_rtolerance;
    }

    /**
     * \return the relative tolerance SNES
     */
    real_type rToleranceSNES() const
    {
        return M_rtoleranceSNES;
    }

    /**
     * \return the divergence tolerance
     */
    real_type dTolerance() const
    {
        return M_dtolerance;
    }

    /**
     * \return the SNES step length tolerance
     */
    real_type sToleranceSNES() const
    {
        return M_stoleranceSNES;
    }

    /**
     * \return the absolute tolerance
     */
    real_type aTolerance() const
    {
        return M_atolerance;
    }

    /**
     * \return the SNES absolute tolerance
     */
    real_type aToleranceSNES() const
    {
        return M_atoleranceSNES;
    }

    /**
     * \return the maximum number of iterations
     */
    size_type maxIterations() const
    {
        return M_maxitKSP;
    }
    size_type maxIterationsKSP() const
    {
        return M_maxitKSP;
    }
    size_type maxIterationsKSPinSNES() const
    {
        return M_maxitKSPinSNES;
    }
    size_type maxIterationsSNES() const
    {
        return M_maxitSNES;
    }
    size_type maxIterationsKSPReuse() const
    {
        return M_maxitKSPReuse;
    }
    size_type maxIterationsKSPinSNESReuse() const
    {
        return M_maxitKSPinSNESReuse;
    }
    size_type maxIterationsSNESReuse() const
    {
        return M_maxitSNESReuse;
    }


    /**
     * \return the KSP relative tolerance in SNES
     */
    real_type rtoleranceKSPinSNES() const
    {
        return M_rtoleranceKSPinSNES;
    }


    bool converged() const
    {
        return M_converged;
    }

    size_type nIterations() const
    {
        return M_iteration;
    }

    bool transpose() const
    {
        return M_transpose;
    }

    bool showKSPMonitor() const { return M_showKSPMonitor; }
    bool showKSPConvergedReason() const { return M_showKSPConvergedReason; }

    bool reusePrec() const { return M_reuse_prec; }
    bool reuseJac() const { return M_reuse_jac; }

    bool reusePrecRebuildAtFirstNewtonStep() const { return M_reusePrecRebuildAtFirstNewtonStep; }
    bool reuseJacRebuildAtFirstNewtonStep() const { return M_reuseJacRebuildAtFirstNewtonStep; }

    BackendType type() const { return M_backend; }

    static std::string enumToKind( BackendType bt )
        {
            if ( bt == BACKEND_EIGEN ) return "eigen";
            if ( bt == BACKEND_EIGEN_DENSE ) return "eigen_dense";
#if FEELPP_HAS_PETSC
            if ( bt == BACKEND_PETSC ) return "petsc";
#endif
            LOG(WARNING) << "Unknown backend, setting up for 'none'";
            return "none";
        }

    /**
     * @return the current datamap of the backend
     */
    datamap_ptrtype dataMap() const { return M_datamap; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set tolerances: relative tolerance \p rtol, divergence tolerance \p dtol
     * and absolute tolerance \p atol
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     setTolerances,
                                     tag,
                                     ( required
                                       ( rtolerance, ( real_type ) )
                                     )
                                     ( optional
                                       ( maxit,      ( size_type ), 1000 )
                                       ( atolerance, ( real_type ),    1e-50 )
                                       ( dtolerance, ( real_type ),    1e5 )
                                     ) )
    {
        M_rtolerance = rtolerance;
        M_dtolerance = dtolerance;
        M_atolerance = atolerance;
        M_maxitKSP = maxit;
    }
#endif
    template <typename ... Ts>
    void setTolerances( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        real_type rtolerance = args.get(_rtolerance );
        size_type maxit = args.get_else( _maxit, 1000 );
        real_type atolerance = args.get_else( _atolerance, 1e-50 );
        real_type dtolerance = args.get_else( _dtolerance, 1e5 );
        M_rtolerance = rtolerance;
        M_dtolerance = dtolerance;
        M_atolerance = atolerance;
        M_maxitKSP = maxit;
    }
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     setTolerancesSNES,
                                     tag,
                                     ( required
                                       ( rtolerance, ( double ) )
                                     )
                                     ( optional
                                       ( maxit,      ( size_type ), 50 )
                                       ( atolerance, ( double ),    1e-50 )
                                       ( stolerance, ( double ),    1e-8 )
                                     ) )
    {
        M_rtoleranceSNES = rtolerance;
        M_stoleranceSNES = stolerance;
        M_atoleranceSNES = atolerance;
        M_maxitSNES = maxit;
    }
#endif
    template <typename ... Ts>
    void setTolerancesSNES( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        real_type rtolerance = args.get(_rtolerance );
        size_type maxit = args.get_else( _maxit, 50 );
        real_type atolerance = args.get_else( _atolerance, 1e-50 );
        real_type stolerance = args.get_else( _stolerance, 1e-8 );
        M_rtoleranceSNES = rtolerance;
        M_stoleranceSNES = stolerance;
        M_atoleranceSNES = atolerance;
        M_maxitSNES = maxit;
    }

    /**
     * set solver: krylov subspace method and preconditioners
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     setSolverType,
                                     tag,
                                     ( required
                                       ( ksp, ( std::string ) )
                                     )
                                     ( optional
                                       ( pc,      ( std::string ), "lu" )
                                       ( constant_null_space,      ( bool ), false )
                                       ( pcfactormatsolverpackage,  ( std::string ), "petsc" )
                                     ) )
    {
        M_ksp = ksp;
        M_pc = pc;
        M_pcFactorMatSolverPackage = pcfactormatsolverpackage;
        M_constant_null_space = constant_null_space;
    }
#endif
    template <typename ... Ts>
    void setSolverType( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        std::string const& ksp = args.get(_ksp );
        std::string const& pc = args.get_else( _pc, "lu" );
        bool constant_null_space = args.get_else( _constant_null_space, false );
        std::string const& pcfactormatsolverpackage = args.get_else( _pcfactormatsolverpackage, "petsc" );
        M_ksp = ksp;
        M_pc = pc;
        M_pcFactorMatSolverPackage = pcfactormatsolverpackage;
        M_constant_null_space = constant_null_space;
    }

    /**
     * set the type of preconditioner associated to the matrix
     */
    void setPrecMatrixStructure( MatrixStructure mstruct )
    {
        M_prec_matrix_structure = mstruct;
    }

    /**
     * \return the non linear solver
     */
    solvernonlinear_ptrtype nlSolver()
    {
        return M_nlsolver;
    }
    solvernonlinear_ptrtype const& nlSolver() const
    {
        return M_nlsolver;
    }

    void setTranspose( bool transpose )
    {
        M_transpose = transpose;
    }

    void setShowKSPMonitor( bool b ) { M_showKSPMonitor=b; }
    void setShowKSPConvergedReason( bool b ) { M_showKSPConvergedReason=b; }

    void setReusePrec( bool b ) { M_reuse_prec=b; }
    void setReuseJac( bool b) { M_reuse_jac=b; }

    void setReusePrecRebuildAtFirstNewtonStep(bool b) { M_reusePrecRebuildAtFirstNewtonStep=b; }
    void setReuseJacRebuildAtFirstNewtonStep(bool b) { M_reuseJacRebuildAtFirstNewtonStep=b; }

    /**
     * @brief set the current datamap of the backend
     * this is used for example in the pre/post solve functions
     * to pass on the parallel data layout
     */
    void setDataMap( datamap_ptrtype dm ) { M_datamap = dm; }

    //@}

    /** @name  Methods
     */
    //@{

    //!
    //! build a new backend with the same properties
    //!
    backend_ptrtype clone()
    {
        return build( enumToKind( this->type() ), this->prefix(), this->worldCommPtr() );
    }

    /**
     * clean up
     */
    virtual void clear();

    /**
     * \return \f$ r = x^T * y \f$
     */
    virtual value_type dot( vector_type const& x, vector_type const& y ) const;


    /**
     * \return \f$ r = x^T * y \f$
     */
    value_type dot( vector_ptrtype const& x, vector_ptrtype const& y ) const
    {
        return this->dot( *x, *y );
    }
    /**
     * \return \f$ y = A * x \f$
     */
    virtual void prod( sparse_matrix_type const& A, vector_type const& x, vector_type& y, bool transpose = false ) const = 0;

    /**
     * \return \f$ y = A * x \f$
     */
    void prod( sparse_matrix_ptrtype const& A, vector_ptrtype const& x, vector_ptrtype& y, bool transpose = false ) const
    {
        this->prod( *A, *x, *y, transpose );
    }

    /**
     * get the matrix \c M whose diagonal is \c v
     */
    virtual int diag( vector_ptrtype const& v, sparse_matrix_ptrtype& M ) const
        {
            return diag( *v, *M );
        }

    /**
     * get the matrix \c M whose diagonal is \c v
     */
    virtual int diag( vector_type const& v, sparse_matrix_type& M ) const
        {
            CHECK(0) << "Invalid call to diag(v,M). Not implemented in Backend base class";
            return 0;
        }

    /**
     * @return the vector \c v with diagonal of \c M
     */
    virtual int diag( sparse_matrix_ptrtype const& M, vector_ptrtype& v ) const
        {
            return diag( *M, *v );
        }
    /**
     * @return the vector \c v with diagonal of \c M
     */
    virtual int diag( sparse_matrix_type const& M, vector_type& v ) const
        {
            CHECK(0) << "Invalid call to diag(M,v). Not implemented in Backend base class";
            return 0;
        }

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation
     * of the inverse of \f$A\f$. this interface uses the
     * boost.parameter library to ease the function usage
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param reuse_prec if true use adaptive preconditioning strategy
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( solve_return_type ),
                                     solve,
                                     tag,
                                     ( required
                                       //( matrix,*(mpl::or_<sparse_matrix_ptrtype,shell_matrix_ptrtype>) )
                                       ( matrix,(sparse_matrix_ptrtype) )
                                       //( in_out( solution ),*( mpl::or_<mpl::or_<boost::is_convertible<mpl::_,vector_type&>,boost::is_convertible<mpl::_,vector_type> >,boost::is_convertible<mpl::_,vector_ptrtype> > ) )
                                       ( in_out( solution ),* )
                                       ( rhs,( vector_ptrtype ) ) )
                                     ( optional
                                       //(prec,(sparse_matrix_ptrtype), matrix )
                                       ( prec,( preconditioner_ptrtype ), preconditioner_ptrtype() )
                                       ( null_space,( NullSpace<value_type> ), NullSpace<value_type>() )
                                       ( near_null_space,( NullSpace<value_type> ), NullSpace<value_type>() )
                                       ( maxit,( size_type ), M_maxitKSP/*1000*/ )
                                       ( rtolerance,( double ), M_rtolerance/*1e-13*/ )
                                       ( atolerance,( double ), M_atolerance/*1e-50*/ )
                                       ( dtolerance,( double ), M_dtolerance/*1e5*/ )
                                       ( reuse_prec,( bool ), M_reuse_prec )
                                       ( transpose,( bool ), false )
                                       ( close,( bool ), true )
                                       ( constant_null_space,( bool ), M_constant_null_space/*false*/ )
                                       ( pre, (pre_solve_type), pre_solve_type() )
                                       ( post, (post_solve_type), post_solve_type() )
                                       ( pc,( std::string ),M_pc/*"lu"*/ )
                                       ( ksp,( std::string ),M_ksp/*"gmres"*/ )
                                       ( pcfactormatsolverpackage,( std::string ), M_pcFactorMatSolverPackage )
                                       ( verbose,   ( bool ), M_verbose )
                                     )
                                   )
#endif
    template <typename ... Ts>
    solve_return_type solve( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        sparse_matrix_ptrtype matrix = args.get(_matrix);
        auto && solution = args.get(_solution);
        vector_ptrtype rhs = args.get(_rhs);
        preconditioner_ptrtype prec = args.get_else(_prec,preconditioner_ptrtype{});
        NullSpace<value_type> const& null_space = args.get_else(_null_space,NullSpace<value_type>{} );
        NullSpace<value_type> const& near_null_space = args.get_else(_near_null_space,NullSpace<value_type>{} );
        size_type maxit = args.get_else(_maxit,M_maxitKSP );
        double rtolerance = args.get_else(_rtolerance,M_rtolerance );
        double atolerance = args.get_else(_atolerance,M_atolerance );
        double dtolerance = args.get_else(_dtolerance,M_dtolerance );
        bool reuse_prec = args.get_else(_reuse_prec, M_reuse_prec);
        bool transpose = args.get_else(_transpose,false);
        bool close = args.get_else(_close,true);
        bool constant_null_space = args.get_else(_constant_null_space,M_constant_null_space);
        pre_solve_type pre = args.get_else(_pre,pre_solve_type{});
        post_solve_type post = args.get_else(_post,post_solve_type{});
        std::string const& pc = args.get_else(_pc,M_pc);
        std::string const& ksp = args.get_else(_ksp,M_ksp);
        std::string const& pcfactormatsolverpackage = args.get_else(_pcfactormatsolverpackage,M_pcFactorMatSolverPackage);
        bool verbose = args.get_else(_verbose, M_verbose );


        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::solve begin" );
        }
        this->setTolerances( _dtolerance=dtolerance,
                             _rtolerance=rtolerance,
                             _atolerance=atolerance,
                             _maxit=maxit );

        this->setSolverType( _pc=pc, _ksp=ksp,
                             _constant_null_space=constant_null_space,
                             _pcfactormatsolverpackage = pcfactormatsolverpackage );

        // preconditioner
        if ( !prec )
        {
            if ( !M_preconditioner )
                M_preconditioner = Feel::preconditioner( _prefix=this->prefix(),_matrix=matrix,_pc=this->pcEnumType()/*LU_PRECOND*/,
                                                         _pcfactormatsolverpackage=this->matSolverPackageEnumType(), _backend=this->shared_from_this(),
                                                         _worldcomm=this->comm() );
            else
            {
                M_preconditioner->setType( this->pcEnumType() );
                M_preconditioner->setMatSolverPackageType( this->matSolverPackageEnumType() );
                M_preconditioner->setMatrix( matrix );
            }
        }
        else
            this->attachPreconditioner( prec );

        // attach null space (or near null space for multigrid) in backend
        auto mynullspace = std::make_shared<NullSpace<value_type>>(this->shared_from_this(),null_space);
        auto myNearNullSpace = std::make_shared<NullSpace<value_type>>(this->shared_from_this(),near_null_space);
        if ( mynullspace->size() > 0 )
        {
            this->attachNullSpace( mynullspace );
            if ( myNearNullSpace->size() > 0 )
                this->attachNearNullSpace( myNearNullSpace );
            else
                this->attachNearNullSpace( mynullspace );
        }
        else if ( myNearNullSpace->size() > 0 )
        {
            this->attachNearNullSpace( myNearNullSpace );
        }

        // make sure matrix and rhs are closed
        if ( close )
        {
            matrix->close();
            rhs->close();
        }

        // print them in matlab format
        if ( !M_export.empty() )
        {
            matrix->printMatlab( M_export+"_A.m" );
            rhs->printMatlab( M_export+"_b.m" );
        }
        // set pre/post solve functions
        M_post_solve = post;
        M_pre_solve = pre;
        auto dm = Feel::detail::datamap( solution );
        this->setDataMap( dm );

        vector_ptrtype _sol( this->toBackendVectorPtr( solution ) );
        bool needToCopySolution = false;
        if( !_sol )
        {
            _sol = this->newVector( dm );
            *_sol = Feel::detail::ref( solution );
            needToCopySolution = true;
        }
        vector_ptrtype _rhs( this->toBackendVectorPtr( rhs ) );
        CHECK( _rhs ) << "converstion to backend vector of rhs fails";

        this->setTranspose( transpose );
        solve_return_type ret;

        if ( reuse_prec == false )
        {
            ret = solve( matrix, matrix, _sol, _rhs );
        }

        else
            ret = solve( matrix, matrix, _sol, _rhs, reuse_prec );

        _sol->close();

        if ( needToCopySolution )
        {
            Feel::detail::ref( solution ) = *_sol;
            Feel::detail::ref( solution ).close();
        }

        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::solve end" );
        }
        return ret;
    }

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation of
     * the inverse of \f$A\f$.
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
    virtual solve_return_type solve( sparse_matrix_ptrtype const& A,
                                     sparse_matrix_ptrtype const& P,
                                     vector_ptrtype& x,
                                     vector_ptrtype const& b
                                   ) = 0;

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation
     * of the inverse of \f$A\f$ using an adaptive preconditioning
     * strategy.
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& P,
                             vector_ptrtype& x,
                             vector_ptrtype const& b,
                             bool reuse_prec
                           );


    /**
     * solve for \f$P F(x)=0 b\f$
     */
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( nl_solve_return_type ),
                                     nlSolve,
                                     tag,
                                     ( required
                                       //( in_out( solution ),*( mpl::or_<boost::is_convertible<mpl::_,vector_type&>,boost::is_convertible<mpl::_,vector_ptrtype> > ) ) )
                                       ( in_out( solution ),*))
                                     ( optional
                                       ( jacobian,( sparse_matrix_ptrtype ), sparse_matrix_ptrtype() )
                                       ( residual,( vector_ptrtype ), vector_ptrtype() )
                                       //(prec,(sparse_matrix_ptrtype), jacobian )
                                       ( prec,( preconditioner_ptrtype ), preconditioner_ptrtype() )
                                       ( null_space,( NullSpace<value_type> ), NullSpace<value_type>() )
                                       ( near_null_space,( NullSpace<value_type> ), NullSpace<value_type>() )
                                       ( maxit,( size_type ), M_maxitSNES/*50*/ )
                                       ( rtolerance,( double ), M_rtoleranceSNES/*1e-8*/ )
                                       ( atolerance,( double ), M_atoleranceSNES/*1e-50*/ )
                                       ( stolerance,( double ), M_stoleranceSNES/*1e-8*/ )
                                       ( reuse_prec,( bool ), M_reuse_prec )
                                       ( reuse_jac,( bool ), M_reuse_jac )
                                       ( transpose,( bool ), false )
                                       ( pre, (pre_solve_type), pre_solve_type() )
                                       ( post, (post_solve_type), post_solve_type() )
                                       ( update, (update_nlsolve_type), update_nlsolve_type() )
                                       ( pc,( std::string ),M_pc/*"lu"*/ )
                                       ( ksp,( std::string ),M_ksp/*"gmres"*/ )
                                       ( pcfactormatsolverpackage,( std::string ), M_pcFactorMatSolverPackage )
                                       ( verbose,   ( bool ), M_verbose )
                                     )
                                   )
#endif
    template <typename ... Ts>
    nl_solve_return_type nlSolve( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && solution = args.get(_solution);
        sparse_matrix_ptrtype jacobian = args.get_else(_jacobian,sparse_matrix_ptrtype{});
        vector_ptrtype residual = args.get_else(_residual,vector_ptrtype{});
        preconditioner_ptrtype prec = args.get_else(_prec,preconditioner_ptrtype{});
        NullSpace<value_type> const& null_space = args.get_else(_null_space,NullSpace<value_type>{} );
        NullSpace<value_type> const& near_null_space = args.get_else(_near_null_space,NullSpace<value_type>{} );
        size_type maxit = args.get_else(_maxit,M_maxitSNES );
        double rtolerance = args.get_else(_rtolerance,M_rtoleranceSNES );
        double atolerance = args.get_else(_atolerance,M_atoleranceSNES );
        double stolerance = args.get_else(_stolerance,M_stoleranceSNES );
        bool reuse_prec = args.get_else(_reuse_prec, M_reuse_prec );
        bool reuse_jac = args.get_else(_reuse_jac, M_reuse_jac );
        bool transpose = args.get_else(_transpose, false );
        pre_solve_type pre = args.get_else(_pre,pre_solve_type{});
        post_solve_type post = args.get_else(_post,post_solve_type{});
        update_nlsolve_type update = args.get_else(_update,update_nlsolve_type{});
        std::string const& pc = args.get_else(_pc,M_pc);
        std::string const& ksp = args.get_else(_ksp,M_ksp);
        std::string const& pcfactormatsolverpackage = args.get_else(_pcfactormatsolverpackage,M_pcFactorMatSolverPackage);
        bool verbose = args.get_else(_verbose, M_verbose );

        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::nlSolve begin" );
        }
        this->setTolerancesSNES( _stolerance=stolerance,
                                 _rtolerance=rtolerance,
                                 _atolerance=atolerance,
                                 _maxit=maxit );
        this->setSolverType( _pc=pc, _ksp=ksp,
                             _pcfactormatsolverpackage = pcfactormatsolverpackage );

        auto dm = Feel::detail::datamap( solution );
        this->setDataMap( dm );

        // preconditioner
        if ( !prec )
        {
            if ( !M_preconditioner )
                M_preconditioner = Feel::preconditioner( _prefix=this->prefix(),_pc=this->pcEnumType()/*LU_PRECOND*/,
                                                         _pcfactormatsolverpackage=this->matSolverPackageEnumType(), _backend=this->shared_from_this(),
                                                         _worldcomm=this->comm() );
            else
            {
                M_preconditioner->setType( this->pcEnumType() );
                M_preconditioner->setMatSolverPackageType( this->matSolverPackageEnumType() );
            }
        }
        else
            this->attachPreconditioner( prec );

        vector_ptrtype _sol( this->toBackendVectorPtr( solution ) );
        bool needToCopySolution = false;
        if( !_sol )
        {
            _sol = this->newVector( dm );
            *_sol = Feel::detail::ref( solution );
            _sol->close();
            needToCopySolution = true;
        }

        this->setTranspose( transpose );
        solve_return_type ret;

        // residual vector
        vector_ptrtype _res;
        if ( !residual )
            _res = this->newVector( dm );
        else
            _res = this->toBackendVectorPtr( residual );

        //this->nlSolver()->setPrefix( this->prefix() );
        if ( !jacobian )
        {
            this->nlSolver()->jacobian( _sol, jacobian );
            jacobian->close();
        }

        if ( !this->nlSolver()->initialized() )
            this->nlSolver()->attachPreconditioner( M_preconditioner );

        // attach null space (or near null space for multigrid) in backend
        auto mynullspace = std::make_shared<NullSpace<value_type>>(this->shared_from_this(),null_space);
        auto myNearNullSpace = std::make_shared<NullSpace<value_type>>(this->shared_from_this(),near_null_space);
        if ( mynullspace->size() > 0 )
        {
            this->attachNullSpace( mynullspace );
            if ( myNearNullSpace->size() > 0 )
                this->attachNearNullSpace( myNearNullSpace );
            else
                this->attachNearNullSpace( mynullspace );
        }
        else if ( myNearNullSpace->size() > 0 )
        {
            this->attachNearNullSpace( myNearNullSpace );
        }
        this->nlSolver()->setPreSolve( pre );
        this->nlSolver()->setPostSolve( post );
        this->nlSolver()->setUpdateIteration( update );

        //if ( reuse_prec == false && reuse_jac == false )
        //    ret = nlSolve( jacobian, _sol, residual, rtolerance, maxit );
        //else
        ret = nlSolve( jacobian, _sol, _res, rtolerance, maxit, reuse_prec, reuse_jac );
        _sol->close();

        if ( needToCopySolution )
        {
            Feel::detail::ref( solution ) = *_sol;
            Feel::detail::ref( solution ).close();
        }
        if ( verbose )
        {
            Environment::logMemoryUsage( "backend::nlSolve end" );
        }
        return ret;
    }

    /**
     * solve for the nonlinear problem \f$F( u ) = 0\f$
     */
    virtual FEELPP_DEPRECATED nl_solve_return_type nlSolve( sparse_matrix_ptrtype& A,
                                                            vector_ptrtype& x,
                                                            vector_ptrtype& b,
                                                            const double, const int );

    /**
     * solve for the nonlinear problem \f$F( u ) = 0\f$ with an
     * adaptive strategy to reuse the preconditioner
     */
    virtual nl_solve_return_type nlSolve( sparse_matrix_ptrtype& A,
                                          vector_ptrtype& x,
                                          vector_ptrtype& b,
                                          const double, const int,
                                          bool reusePC, bool reuseJAC );

    /**
     * assemble \f$C=P^T A P\f$
     */
    virtual int PtAP( sparse_matrix_ptrtype const& A,
                       sparse_matrix_ptrtype const& P,
                       sparse_matrix_ptrtype & C
                       ) const;

    /**
     * assemble \f$C=P A P^T\f$
     */
    virtual int PAPt( sparse_matrix_ptrtype const& A,
                      sparse_matrix_ptrtype const& P,
                      sparse_matrix_ptrtype& C ) const;

    /**
     * attach the default preconditioner
     */
    void attachPreconditioner()
    {
    }

    /**
     * Attaches a Preconditioner object to be used by the solver
     */
    void attachPreconditioner( preconditioner_ptrtype preconditioner )
    {
        // if ( M_preconditioner && M_preconditioner != preconditioner )
        //     M_preconditioner->clear();
        M_preconditioner = preconditioner;
    }

    /**
     * @return the preconditioner attached to the backend
     */
    preconditioner_ptrtype preconditioner()
    {
        return M_preconditioner;
    }


    void attachNullSpace( std::shared_ptr<NullSpace<value_type> > nullSpace )
    {
        M_nullSpace = nullSpace;
    }
    void attachNearNullSpace( std::shared_ptr<NullSpace<value_type> > nearNullSpace )
    {
        M_nearNullSpace = nearNullSpace;
    }

    /**
     * register a backend observer for the delete signal of backend
     */
    template<typename Observer>
    void
    addDeleteObserver( Observer const& obs )
        {
            M_deleteObservers.connect( obs );
        }
    /**
     * register a backend observer for the delete signal of backend that is a
     * shared_ptr<>
     */
    template<typename Observer>
    void
    addDeleteObserver( std::shared_ptr<Observer> const& obs )
        {
            M_deleteObservers.connect(boost::bind(&Observer::operator(), obs));
        }
    /**
     * send the delete signal to all observers
     */
    void
    sendDeleteSignal()
        {
            M_deleteObservers();
        }

    /**
     * \return the pre solve function
     */
    pre_solve_type preSolve() { return M_pre_solve; }

    /**
     * call the pre solve function with \p x as the rhs and \p y as the solution
     */
    void preSolve(vector_ptrtype x, vector_ptrtype y) { return M_pre_solve(x,y); }

    /**
     * \return the post solve function
     */
    post_solve_type postSolve() { return M_post_solve; }

    /**
     * call the post solve function with \p x as the rhs and \p y as the solution
     */
    void postSolve(vector_ptrtype x, vector_ptrtype y) { return M_post_solve(x,y); }

    //@}



protected:
    preconditioner_ptrtype M_preconditioner;
    std::shared_ptr<NullSpace<value_type> > M_nullSpace, M_nearNullSpace;
private:

    void start();

    void stop();

    void reset();

private:
    po::variables_map M_vm;

protected:
    BackendType M_backend;

private:
    std::string M_prefix;

    solvernonlinear_ptrtype M_nlsolver;

    MatrixStructure M_prec_matrix_structure;

    double M_totalSolveIter;
    double M_lastSolveIter;
    double M_firstSolveTime;
    double M_residual;
    double M_rtolerance;
    double M_dtolerance;
    double M_atolerance;
    double M_rtoleranceSNES, M_stoleranceSNES, M_atoleranceSNES;
    double M_rtoleranceKSPinSNES;

    bool M_reuse_prec;
    bool M_reuse_jac;
    bool M_reusePrecIsBuild,M_reusePrecRebuildAtFirstNewtonStep;
    bool M_reuseJacIsBuild,M_reuseJacRebuildAtFirstNewtonStep;
    size_t M_nUsePC;
    bool   M_converged;
    bool   M_reusePC;
    bool   M_reusedPC;
    bool   M_reuseFailed;
    boost::timer M_timer;
    bool   M_transpose;
    size_type    M_maxitKSP, M_maxitKSPinSNES, M_maxitSNES;
    size_type    M_maxitKSPReuse, M_maxitKSPinSNESReuse, M_maxitSNESReuse;
    size_type    M_iteration;
    std::string M_export;
    std::string M_ksp;
    std::string M_pc;
    std::string M_fieldSplit;
    std::string M_pcFactorMatSolverPackage;
    bool M_constant_null_space;
    bool M_showKSPMonitor;
    bool M_showKSPConvergedReason;
    //std::map<std::string,boost::tuple<std::string,std::string> > M_sub;
    pre_solve_type M_pre_solve;
    post_solve_type M_post_solve;
    datamap_ptrtype M_datamap;
    boost::signals2::signal<void()> M_deleteObservers;

    bool M_verbose;
};


typedef Backend<double,uint32_type> backend_type;
typedef std::shared_ptr<backend_type> backend_ptrtype;

typedef Backend<std::complex<double>> c_backend_type;
typedef std::shared_ptr<c_backend_type> c_backend_ptrtype;

template<typename T = double, typename SizeT = uint32_type>
using backend_ptr_t = std::shared_ptr<Backend<T,SizeT>>;

namespace detail
{
template<typename T>
class BackendManagerImpl:
        public std::map<boost::tuple<std::string,std::string,int>, typename Backend<T>::ptrtype >,
        public boost::noncopyable
{
public:
    typedef typename Backend<T>::ptrtype value_type;
    typedef boost::tuple<std::string,std::string,int> key_type;
    typedef std::map<key_type, value_type> backend_manager_type;

};
template<typename T>
struct BackendManager : public  Feel::Singleton<BackendManagerImpl<T>> {};

template<typename T>
struct BackendManagerDeleterImpl
{
    void operator()() const
        {
            VLOG(2) << "[BackendManagerDeleter] clear BackendManager Singleton: " << Feel::detail::BackendManager<T>::instance().size() << "\n";
            Feel::detail::BackendManager<T>::instance().clear();
            VLOG(2) << "[BackendManagerDeleter] clear BackendManager done\n";
        }
};
template<typename T>
struct BackendManagerDeleter
    : public  Feel::Singleton<BackendManagerDeleterImpl<T>>
{};




template<typename T>
typename Backend<T>::ptrtype
backend_impl( std::string const& name, std::string const& kind, bool rebuild, worldcomm_ptr_t const& worldcomm )
{
    // register the BackendManager into Feel::Environment so that it gets the
    // BackendManager is cleared up when the Environment is deleted
    static bool observed=false;
    if ( !observed )
    {
        Environment::addDeleteObserver( Feel::detail::BackendManagerDeleter<T>::instance() );
        observed = true;
    }

    auto git = Feel::detail::BackendManager<T>::instance().find( boost::make_tuple( kind, name, worldcomm->globalSize() ) );

    if (  git != Feel::detail::BackendManager<T>::instance().end() && ( rebuild == false ) )
    {
        DVLOG(2) << "[backend] found backend name=" << name << " kind=" << kind << " rebuild=" << rebuild << " worldcomm->globalSize()=" << worldcomm->globalSize() << "\n";
        return git->second;
    }

    else
    {
        if (  git != Feel::detail::BackendManager<T>::instance().end() && ( rebuild == true ) )
            git->second->clear();

        DVLOG(2) << "[backend] building backend name=" << name << " kind=" << kind << " rebuild=" << rebuild << " worldcomm->globalSize()=" << worldcomm->globalSize() << "\n";

        typename Backend<T>::ptrtype b;
        b = Feel::Backend<T>::build( kind, name, worldcomm );
        DVLOG(2) << "[backend] storing backend in singleton" << "\n";
        Feel::detail::BackendManager<T>::instance().operator[]( boost::make_tuple( kind, name, worldcomm->globalSize() ) ) = b;
        return b;
    }

}

} // detail


#if 0
BOOST_PARAMETER_FUNCTION(
                         ( backend_ptrtype ), // return type
                         backend,           // 2. function name
                         tag,               // 3. namespace of tag types
                         ( optional
                           ( name,           ( std::string ), "" )
                           ( kind,           ( std::string ), soption(_prefix=name,_name="backend"))
                           ( rebuild,        ( bool ), boption(_prefix=name,_name="backend.rebuild") )
                           ( worldcomm,      (worldcomm_ptr_t), Environment::worldCommPtr() )
                           ) )
{
    return Feel::detail::backend_impl<double>( name, kind, rebuild, worldcomm);
}


/*
 * Complex backend
 */
BOOST_PARAMETER_FUNCTION(
                         ( c_backend_ptrtype ), // return type
                         cbackend,           // 2. function name
                         tag,               // 3. namespace of tag types
                         ( optional
                           ( name,           ( std::string ), "" )
                           ( kind,           ( std::string ), soption(_prefix=name,_name="backend"))
                           ( rebuild,        ( bool ), boption(_prefix=name,_name="backend.rebuild") )
                           ( worldcomm,      (worldcomm_ptr_t), Environment::worldCommPtr() )
                           ) )
{
    return Feel::detail::backend_impl<std::complex<double>>( name, kind, rebuild, worldcomm);
}
#endif

template <typename ... Ts>
backend_ptrtype backend( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    std::string const& name = args.get_else(_name,"");
    std::string const& kind = args.get_else(_kind,soption(_prefix=name,_name="backend"));
    bool rebuild = args.get_else(_rebuild,boption(_prefix=name,_name="backend.rebuild"));
    worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr());
    return Feel::detail::backend_impl<double>( name, kind, rebuild, worldcomm);
}
template <typename ... Ts>
c_backend_ptrtype cbackend( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    std::string const& name = args.get_else(_name,"");
    std::string const& kind = args.get_else(_kind,soption(_prefix=name,_name="backend"));
    bool rebuild = args.get_else(_rebuild,boption(_prefix=name,_name="backend.rebuild"));
    worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr());
    return Feel::detail::backend_impl<std::complex<double>>( name, kind, rebuild, worldcomm);
}

template<typename T>
bool isMatrixInverseSymmetric ( std::shared_ptr<MatrixSparse<T> >& A, std::shared_ptr<MatrixSparse<T> >& At, bool print=false  )
{
#if FEELPP_HAS_PETSC
    auto u = Backend<T>::build( BACKEND_PETSC, A->comm() )->newVector(A->size1(), A->size1());
    auto v = Backend<T>::build( BACKEND_PETSC, A->comm() )->newVector(A->size1(), A->size1());

    auto res_u = Backend<T>::build( BACKEND_PETSC, A->comm() )->newVector(A->size1(), A->size1());
    auto res_v = Backend<T>::build( BACKEND_PETSC, A->comm() )->newVector(A->size1(), A->size1());


    for (size_type i = 0; i < u->size(); i++)
    {
        u->set(i,(double(std::rand())/double(RAND_MAX)));
    }

    for (size_type i = 0; i < v->size(); i++)
    {
        v->set(i,(double(std::rand())/double(RAND_MAX)));
    }

    Backend<T>::build( BACKEND_PETSC, A->comm() )->solve(_matrix=A,
                                                         _solution=res_u,
                                                         _rhs=u,
                                                         _pcfactormatsolverpackage="mumps"
                                                         );



    Backend<T>::build( BACKEND_PETSC, A->comm() )->solve(_matrix=At,
                                                         _solution=res_v,
                                                         _rhs=v,
                                                         _pcfactormatsolverpackage="mumps"
                                                         );



    T val1 = inner_product(res_u,v);
    T val2 = inner_product(res_v,u);

    T res = math::abs(val1-val2);

    if ((res >= 1e-12) && print)
    {
        std::cout<<"-----------Subdomain "<< Environment::worldComm().rank() <<"-----------\n";
        std::cout<<"--Val1= "<< val1 <<"\n";
        std::cout<<"--Val2= "<< val2 <<"\n";
        std::cout<<"--|Val1-val2|= "<< res <<"\n";
        std::cout<<"---------------------------------\n";
    }

    return  res < 1e-12;
#else
    LOG(WARNING) << "isMatrixInverseSymmetric: Petsc is not available. This function will always return false.";
    return false;
#endif

}

#if !defined(FEELPP_BACKEND_NOEXTERN)
extern template class Backend<double>;
//extern template class Backend<std::complex<double>>;
#endif

}

#include <feel/feelalg/vectorblock.hpp>

#endif /* Backend_H */
