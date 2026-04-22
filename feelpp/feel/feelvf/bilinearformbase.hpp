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
#ifndef FEELPP_VF_BILINEARFORMBASE_H
#define FEELPP_VF_BILINEARFORMBASE_H

#include <future>


#include <feel/feelconfig.h>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelvf/block.hpp>
#include <feel/feelvf/dirichletconstraints.hpp>




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
    using vector_type = typename matrix_type::vector_type;
    using vector_ptrtype = typename matrix_type::vector_ptrtype;
    static inline const bool is_row_major = true;//matrix_type::is_row_major;
    using deferred_dirichlet_set_type = vf::DeferredDirichletSet<value_type>;

    struct DeferredDirichletState
    {
        vf::DeferredDirichletPolicy dirichletPolicy = vf::DeferredDirichletPolicy::automatic;
        deferred_dirichlet_set_type pendingConstraints;
        deferred_dirichlet_set_type appliedConstraints;
        matrix_ptrtype constrainedMatrix;
        vector_ptrtype constrainedVector;
        void const* constrainedVectorSource = nullptr;
        std::size_t constrainedVectorSourceRevision = 0;
    };
    using deferred_dirichlet_state_type = DeferredDirichletState;
    using deferred_dirichlet_state_ptrtype = std::shared_ptr<deferred_dirichlet_state_type>;

    using size_type =  typename matrix_type::size_type;
    
    using layout_type = mp11::mp_if_c<is_row_major, ublas::row_major,ublas::column_major>;

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
            //toc(M_name, Environment::logVerbosityLevel() > 0 );
        }

    /**
     * copy operator
     */
    BilinearFormBase&
    operator=( BilinearFormBase const& form );

    BilinearFormBase& operator+=( BilinearFormBase const& a )
        {
            this->invalidateMaterializedDeferredDirichlet();
            if ( this == &a )
            {
                M_matrix->scale( 2.0 );
                return *this;
            }
            M_matrix->addMatrix( 1.0, a.M_matrix );
            return *this;
        }

    BilinearFormBase& add( double alpha, BilinearFormBase const&  a )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->addMatrix( alpha, a.M_matrix );
            return *this;
        }
    
    BilinearFormBase& operator-=( BilinearFormBase const& a )
        {
            this->invalidateMaterializedDeferredDirichlet();
            if ( this == &a )
            {
                M_matrix->zero();
                return *this;
            }
            M_matrix->addMatrix( -1.0, a.M_matrix );
            return *this;
        }
    /**
     * @brief operator *=
     * 
     * @param __a bilinear form
     * @return Bilinear& 
     */
    BilinearFormBase& operator*=( value_type const& alpha )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->scale( alpha );
            return *this;
        }
    /**
     * @brief operator/= by a scalar
     * 
     * @param __a scalar to divide by
     * @return BilinearForm& 
     */
    BilinearFormBase& operator/=( value_type const& alpha )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->scale( 1.0/alpha );
            return *this;
        }
    //! scale the form
    BilinearFormBase& scale( double alpha )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->scale( alpha );
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
        this->invalidateMaterializedDeferredDirichlet();
        return *M_matrix;
    }

    matrix_ptrtype const& matrixPtr() const
    {
        return M_matrix;
    }

    matrix_ptrtype& matrixPtr()
    {
        this->invalidateMaterializedDeferredDirichlet();
        return M_matrix;
    }

    auto l1Norm() const
    {
        return M_matrix->l1Norm();
    }

    auto linftyNorm() const
    {
        return M_matrix->linftyNorm();
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
        this->clearDeferredDirichlet();
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


    // Close the assembled base operator only. Deferred Dirichlet constraints
    // remain logical form state until materialized explicitly.
    void close()
        {
            this->get(); // futures
            M_matrix->close();
        }

    void closeBaseOperator()
        {
            this->close();
        }

    bool closed() const noexcept
        {
            return M_matrix->closed();
        }

    bool baseOperatorClosed() const noexcept
        {
            return this->closed();
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

    bool useDeferredDirichlet() const noexcept
        {
            return vf::usesDeferredDirichlet( this->dirichletPolicy() );
        }

    vf::DeferredDirichletPolicy dirichletPolicy() const noexcept
        {
            return M_dirichletState ? M_dirichletState->dirichletPolicy :
                vf::DeferredDirichletPolicy::automatic;
        }

    void setDirichletPolicy( vf::DeferredDirichletPolicy policy ) noexcept
        {
            if ( !M_dirichletState && policy == vf::DeferredDirichletPolicy::automatic )
                return;
            if ( !M_dirichletState )
                M_dirichletState = std::make_shared<deferred_dirichlet_state_type>();
            M_dirichletState->dirichletPolicy = policy;
        }

    void setUseDeferredDirichlet( bool value ) noexcept
        {
            this->setDirichletPolicy( value ? vf::DeferredDirichletPolicy::deferred :
                                             vf::DeferredDirichletPolicy::immediate );
        }

    BilinearFormBase& deferDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::deferred );
            return *this;
        }

    BilinearFormBase& autoDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::automatic );
            return *this;
        }

    BilinearFormBase& immediateDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::immediate );
            return *this;
        }

    bool hasDeferredDirichlet() const noexcept
        {
            return M_dirichletState && !M_dirichletState->pendingConstraints.empty();
        }

    bool hasPendingDirichletConstraints() const noexcept
        {
            return this->hasDeferredDirichlet();
        }

    bool hasDirichletConstraints() const noexcept
        {
            return M_dirichletState &&
                   ( !M_dirichletState->pendingConstraints.empty() ||
                     !M_dirichletState->appliedConstraints.empty() );
        }

    bool supportsConstrainedOperatorView() const noexcept
        {
            return true;
        }

    bool hasMaterializedConstrainedOperator() const noexcept
        {
            return M_dirichletState && static_cast<bool>( M_dirichletState->constrainedMatrix );
        }

    bool shouldDeferDirichlet( Feel::Context const& on_context ) const noexcept
        {
            return vf::shouldDeferDirichlet( this->dirichletPolicy(), on_context );
        }

    void deferZeroRows( std::vector<int> const& dofs,
                        std::vector<value_type> const& values,
                        Feel::Context const& on_context,
                        double value_on_diagonal,
                        std::uint8_t entity_priority = vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::unspecified ) )
        {
            auto& state = this->ensureDeferredDirichletState();
            state.pendingConstraints.append( dofs, values, on_context, value_on_diagonal, entity_priority );
            this->invalidateMaterializedDeferredDirichlet();
        }

    void clearDeferredDirichlet() noexcept
        {
            if ( M_dirichletState )
            {
                M_dirichletState->pendingConstraints.clear();
                M_dirichletState->appliedConstraints.clear();
                M_dirichletState->constrainedMatrix.reset();
                M_dirichletState->constrainedVector.reset();
                M_dirichletState->constrainedVectorSource = nullptr;
                M_dirichletState->constrainedVectorSourceRevision = 0;
            }
        }

    deferred_dirichlet_state_ptrtype const& deferredDirichletStatePtr() const noexcept
        {
            return M_dirichletState;
        }

    void shareDeferredDirichletState( deferred_dirichlet_state_ptrtype state )
        {
            M_dirichletState = std::move( state );
        }

    matrix_ptrtype const& baseMatrixPtr() const noexcept
        {
            return M_matrix;
        }

    matrix_ptrtype& baseMatrixPtr() noexcept
        {
            this->invalidateMaterializedDeferredDirichlet();
            return M_matrix;
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto constrainedSystem( RhsType& rhs )
        {
            return this->materializeConstrainedSystem( rhs.vectorPtr() );
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto constrainedSystem( RhsType const& rhs )
        {
            return const_cast<BilinearFormBase*>( this )->materializeConstrainedSystem( rhs.vectorPtr() );
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto activeSystem( RhsType& rhs )
        {
            if ( this->hasDirichletConstraints() )
                return this->constrainedSystem( rhs );
            return std::pair{ this->baseMatrixPtr(), rhs.vectorPtr() };
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto activeSystem( RhsType const& rhs )
        {
            if ( this->hasDirichletConstraints() )
                return const_cast<BilinearFormBase*>( this )->materializeConstrainedSystem( rhs.vectorPtr() );
            return std::pair{ this->baseMatrixPtr(), rhs.vectorPtr() };
        }

    matrix_ptrtype constrainedMatrixPtr()
        {
            return this->materializeConstrainedMatrix();
        }

    matrix_ptrtype constrainedMatrixPtr() const
        {
            return const_cast<BilinearFormBase*>( this )->constrainedMatrixPtr();
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    matrix_ptrtype constrainedMatrixPtr( RhsType& rhs )
        {
            return this->constrainedSystem( rhs ).first;
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    matrix_ptrtype constrainedMatrixPtr( RhsType const& rhs )
        {
            return const_cast<BilinearFormBase*>( this )->constrainedMatrixPtr( rhs );
        }

    matrix_ptrtype activeMatrixPtr()
        {
            return this->hasDirichletConstraints() ? this->constrainedMatrixPtr() : this->baseMatrixPtr();
        }

    matrix_ptrtype activeMatrixPtr() const
        {
            return this->hasDirichletConstraints() ? this->constrainedMatrixPtr() : this->baseMatrixPtr();
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    matrix_ptrtype activeMatrixPtr( RhsType& rhs )
        {
            return this->activeSystem( rhs ).first;
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    matrix_ptrtype activeMatrixPtr( RhsType const& rhs )
        {
            return const_cast<BilinearFormBase*>( this )->activeMatrixPtr( rhs );
        }

    void materializeConstrainedOperator()
        {
            if ( !this->hasDirichletConstraints() )
                return;
            (void)this->constrainedMatrixPtr();
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    void materializeConstrainedOperator( RhsType& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;
            (void)this->constrainedSystem( rhs );
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    void materializeConstrainedOperator( RhsType const& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;
            (void)this->constrainedSystem( rhs );
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto constrainedVectorPtr( RhsType& rhs )
        {
            return this->constrainedSystem( rhs ).second;
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto constrainedVectorPtr( RhsType const& rhs )
        {
            return this->constrainedSystem( rhs ).second;
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto activeVectorPtr( RhsType& rhs )
        {
            return this->activeSystem( rhs ).second;
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto activeVectorPtr( RhsType const& rhs )
        {
            return this->activeSystem( rhs ).second;
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    void applyDeferredDirichlet( RhsType& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;

            auto constrainedVector = this->constrainedVectorPtr( rhs );
            vf::copyVectorValues( rhs.vectorPtr(), constrainedVector );
            auto& state = this->ensureDeferredDirichletState();
            state.constrainedVectorSource = static_cast<void const*>( rhs.vectorPtr().get() );
            state.constrainedVectorSourceRevision = rhs.vectorPtr()->revision();
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    void applyDeferredDirichlet( RhsType const& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;

            auto constrainedVector = this->constrainedVectorPtr( rhs );
            vf::copyVectorValues( rhs.vectorPtr(), constrainedVector );
            auto& state = this->ensureDeferredDirichletState();
            state.constrainedVectorSource = static_cast<void const*>( rhs.vectorPtr().get() );
            state.constrainedVectorSourceRevision = rhs.vectorPtr()->revision();
        }

protected:
    deferred_dirichlet_state_type& ensureDeferredDirichletState()
        {
            if ( !M_dirichletState )
                M_dirichletState = std::make_shared<deferred_dirichlet_state_type>();
            return *M_dirichletState;
        }

    void invalidateMaterializedDeferredDirichlet() noexcept
        {
            if ( M_dirichletState )
            {
                M_dirichletState->constrainedMatrix.reset();
                M_dirichletState->constrainedVector.reset();
                M_dirichletState->constrainedVectorSource = nullptr;
                M_dirichletState->constrainedVectorSourceRevision = 0;
            }
        }

    deferred_dirichlet_set_type allDeferredDirichletConstraints() const
        {
            deferred_dirichlet_set_type constraints;
            if ( !M_dirichletState )
                return constraints;

            constraints.append( M_dirichletState->appliedConstraints );
            constraints.append( M_dirichletState->pendingConstraints );
            return constraints;
        }

    void promotePendingDeferredDirichlet()
        {
            auto& state = this->ensureDeferredDirichletState();
            state.appliedConstraints.append( state.pendingConstraints );
            state.pendingConstraints.clear();
        }

    template<typename VectorPtrType>
    auto materializeConstrainedSystem( VectorPtrType const& rhsVector )
    {
        if ( !this->hasDirichletConstraints() )
            return std::pair{ M_matrix, rhsVector };

        auto& state = this->ensureDeferredDirichletState();
        auto const rhsSource = static_cast<void const*>( rhsVector.get() );
        auto const rhsRevision = rhsVector->revision();
        if ( state.constrainedMatrix &&
             state.constrainedVector &&
             state.pendingConstraints.empty() &&
             state.constrainedVectorSource == rhsSource &&
             state.constrainedVectorSourceRevision == rhsRevision )
        {
            return std::pair{ state.constrainedMatrix, state.constrainedVector };
        }

        this->get();
        M_matrix->closeIfNeeded();
        if ( !rhsVector->closed() )
            rhsVector->close();

        auto constrainedMatrix = M_matrix->clone();
        auto constrainedVector = vf::cloneVectorWithValues( rhsVector );
        auto const deferredEntries = this->allDeferredDirichletConstraints().entries();
        vf::applyDeferredDirichletEntries( deferredEntries, constrainedMatrix, constrainedVector );

        constrainedMatrix->close();
        if ( !constrainedVector->closed() )
            constrainedVector->close();

        state.constrainedMatrix = constrainedMatrix;
        state.constrainedVector = constrainedVector;
        state.constrainedVectorSource = rhsSource;
        state.constrainedVectorSourceRevision = rhsRevision;
        this->promotePendingDeferredDirichlet();
        return std::pair{ state.constrainedMatrix, state.constrainedVector };
    }

    matrix_ptrtype materializeConstrainedMatrix()
        {
            if ( !this->hasDirichletConstraints() )
                return M_matrix;

            auto& state = this->ensureDeferredDirichletState();
            if ( state.constrainedMatrix && state.pendingConstraints.empty() )
                return state.constrainedMatrix;

            this->get();
            M_matrix->closeIfNeeded();

            auto constrainedMatrix = M_matrix->clone();
            auto dummyRhs = Feel::backend( _worldcomm=this->worldCommPtr() )->newVector( M_matrix->mapRowPtr() );
            dummyRhs->zero();
            dummyRhs->close();

            auto const deferredEntries = this->allDeferredDirichletConstraints().entries();
            vf::applyDeferredDirichletEntries( deferredEntries, constrainedMatrix, dummyRhs );

            constrainedMatrix->close();
            state.constrainedMatrix = constrainedMatrix;
            state.constrainedVector.reset();
            state.constrainedVectorSource = nullptr;
            state.constrainedVectorSourceRevision = 0;
            this->promotePendingDeferredDirichlet();
            return state.constrainedMatrix;
        }

    template<typename VectorPtrType>
    auto materializeConstrainedVector( VectorPtrType const& rhsVector )
        {
            if ( !this->hasDirichletConstraints() )
                return rhsVector;
            auto [constrainedMatrix, constrainedVector] = this->materializeConstrainedSystem( rhsVector );
            return constrainedVector;
        }

public:

    /**
     * add value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void add( size_type i,  size_type j,  value_type const& v )
    {
        this->invalidateMaterializedDeferredDirichlet();
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
        this->invalidateMaterializedDeferredDirichlet();
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
        this->invalidateMaterializedDeferredDirichlet();
        M_matrix = backend()->newMatrix( _test=x1, _trial=x2 );
    }
    bool isMatrixAllocated() const
    {
        return (bool)M_matrix;
    }

    template <typename ... Ts>
    typename Backend<value_type>::solve_return_type solve( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            auto && solution = args.get(_solution);
            auto && rhs = args.get(_rhs);
            std::string const& name = args.get_else(_name,"");
            std::string const& kind = args.get_else_invocable(_kind,[&name]() { return soption(_prefix=name,_name="backend"); } );
            bool rebuild = args.get_else_invocable(_rebuild,[&name]() { return boption(_prefix=name,_name="backend.rebuild"); } );
            pre_solve_type pre = args.get_else(_pre,pre_solve_type());
            post_solve_type post = args.get_else(_post,post_solve_type());

            if ( !this->hasDirichletConstraints() )
                this->closeBaseOperator();
            auto [matrix, rhsVector] = this->activeSystem( rhs );
            auto solveBackend = Feel::backend( _name=name, _kind=kind, _rebuild=rebuild,
                                               _worldcomm=this->worldCommPtr() );
            return solveBackend->solve( _matrix=matrix,
                                        _auxiliary_matrix=this->baseMatrixPtr(),
                                        _rhs=rhsVector,
                                        _solution=solution,
                                        _pre=pre,
                                        _post=post
                                        );
        }

    template <typename ... Ts>
    typename Backend<value_type>::solve_return_type solveb( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            auto && solution = args.get(_solution);
            auto && rhs = args.get(_rhs);
            auto && backend = args.get(_backend);
            if ( !this->hasDirichletConstraints() )
                this->closeBaseOperator();
            auto [matrix, rhsVector] = this->activeSystem( rhs );
            preconditioner_ptrtype prec = args.get_else_invocable(_prec, [&backend,&matrix](){ return preconditioner( _prefix=backend->prefix(),
                                                                                                                      _matrix=matrix,
                                                                                                                      _pc=backend->pcEnumType(),
                                                                                                                      _pcfactormatsolverpackage=backend->matSolverPackageEnumType(),
                                                                                                                      _backend=backend ); } );
            return backend->solve( _matrix=matrix, _rhs=rhsVector,
                                   _solution=solution, _prec = prec,
                                   _auxiliary_matrix=this->baseMatrixPtr() );
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

    deferred_dirichlet_state_ptrtype M_dirichletState;

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
    M_threshold( threshold ),
    M_dirichletState( std::make_shared<deferred_dirichlet_state_type>() )
{
    DVLOG(2) << "begin constructor with default listblock\n";

    if ( !Xh->worldComm().isActive() ) return;

    if ( !M_matrix ) M_matrix = backend()->newMatrix( _test=Xh, _trial=Yh );

    M_lb.push_back( Feel::vf::Block ( 0, 0, 0, 0 ) );
    auto dmTest = M_matrix->mapRowPtr();
    auto dmTrial = M_matrix->mapColPtr();
    this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
    this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );

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
    M_threshold( threshold ),
    M_dirichletState( std::make_shared<deferred_dirichlet_state_type>() )
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
