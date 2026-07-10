/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#ifndef FEELPP_VF_BLOCKFORMS_H
#define FEELPP_VF_BLOCKFORMS_H

#include <feel/feelvf/form.hpp>
#include <feel/feelalg/productspaceconcepts.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/matrixcondensed.hpp>
#include <feel/feelalg/vectorcondensed.hpp>
#include <feel/feeldiscr/csrgraphblocks.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/dirichletconstraints.hpp>

#include <ranges>
#include <type_traits>


namespace Feel {

//!
//! forward declarations of @c BlockBilinearForm and @c blocform2()
//!
template<typename TestPS, typename TrialPS = TestPS>
class BlockBilinearForm;


template<typename PS>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps );

template<typename PS, typename BackendT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, BackendT&& b );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, solve::strategy s, BackendT&& b,
            size_type pattern = Pattern::COUPLED,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, solve::strategy s, BackendT&& b,
            std::vector<PatternSizeT> const& patterns,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b,
            size_type pattern = Pattern::COUPLED,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b,
            std::vector<PatternSizeT> const& patterns,
            RangeMapT r = stencilRangeMap() );

template<typename PS,typename T>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, condensed_matrix_ptr_t<T> & m );

namespace detail
{
template<typename T>
struct is_backend_shared_ptr : std::false_type {};

template<typename T>
struct is_backend_shared_ptr<std::shared_ptr<T>> : std::bool_constant<std::is_base_of_v<BackendBase, T>> {};

template<typename T>
inline constexpr bool is_backend_like_v =
    std::is_base_of_v<BackendBase, decay_type<T>> ||
    is_backend_shared_ptr<decay_type<T>>::value;

template<typename T>
struct is_std_vector : std::false_type {};

template<typename T, typename Allocator>
struct is_std_vector<std::vector<T, Allocator>> : std::true_type {};

template<typename T>
inline constexpr bool is_std_vector_v = is_std_vector<decay_type<T>>::value;

template<typename T>
concept BlockformStaticProductSpaces =
    FoldableProductSpacesType<T>;

template<BlockformStaticProductSpaces ProductSpacesT, typename N>
int
blockformFlattenedBlockIndex( ProductSpacesT const& ps, N n, int subIndex )
{
    int index = 0;
    int position = 0;
    int const target = int( n );
    auto const& spaces = ps.tupleSpaces();
    hana::for_each( spaces, [&]( auto const& space )
                    {
                        if ( position < target )
                            index += Feel::detail::csrGraphBlockCount( space );
                        ++position;
                    } );
    return index + subIndex;
}

template<typename... Ts>
concept Blockform2NamedArguments =
    ( sizeof...(Ts) > 0 ) &&
    ( NA::is_named_argument_v<std::decay_t<Ts>> && ... ) &&
    NA::arguments<std::decay_t<Ts>...>::template has_t<na::test>::value;

template<typename... Ts>
concept Blockform2NamedArgumentsWithTrial =
    Blockform2NamedArguments<Ts...> &&
    NA::arguments<std::decay_t<Ts>...>::template has_t<na::trial>::value;

template<typename ArgsT, typename TestT, typename TrialT>
auto blockform2Named( ArgsT& args, TestT&& test, TrialT&& trial )
{
    auto s = args.get_else( _strategy, solve::strategy::monolithic );
    auto&& b = args.get_else_invocable( _backend, []() { return Feel::backend(); } );
    auto&& r = args.get_else( _range, stencilRangeMap() );
    auto&& pattern = args.get_else( _pattern, Pattern::COUPLED );

    if constexpr ( is_std_vector_v<decltype( pattern )> )
    {
        return BlockBilinearForm<TestT, TrialT>( std::forward<TestT>( test ),
                                                 std::forward<TrialT>( trial ),
                                                 s, std::forward<decltype( b )>( b ),
                                                 pattern, r );
    }
    else
    {
        return BlockBilinearForm<TestT, TrialT>( std::forward<TestT>( test ),
                                                 std::forward<TrialT>( trial ),
                                                 s, std::forward<decltype( b )>( b ),
                                                 static_cast<size_type>( pattern ), r );
    }
}
}

//!
//! forward declarations of @c BlockLinearForm and @c blockform1()
//!
template<typename PS>
class BlockLinearForm;

template<typename PS>
BlockLinearForm<PS>
blockform1( PS&& ps );

template<typename PS,typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, BackendT&& b );

template<typename PS,typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, solve::strategy s, BackendT&& b );


template<typename PS,typename T>
BlockLinearForm<PS>
blockform1( PS&& ps, condensed_vector_ptr_t<T> v );


/**
 * Handles bilinear form over a product of spaces
 */
template<typename TestPS, typename TrialPS>
class BlockBilinearForm
{
public :
    using value_type = typename Feel::decay_type<TestPS>::value_type;
    using condensed_matrix_type = MatrixCondensed<value_type>;
    using size_type = typename condensed_matrix_type::size_type;
    using condensed_matrix_ptrtype = std::shared_ptr<condensed_matrix_type>;
    using test_product_space_t = decay_type<TestPS>;
    using trial_product_space_t = decay_type<TrialPS>;
    using product_space_t = test_product_space_t;
    using vector_type = Vector<value_type>;
    using vector_ptrtype = typename vector_type::clone_ptrtype;
    using deferred_dirichlet_set_type = vf::DeferredDirichletSet<value_type>;
    using sparse_matrix_ptrtype = typename condensed_matrix_type::sparse_matrix_ptrtype;

    template<typename SpacePtrType>
    class RowDirichletView
    {
    public:
        using parent_type = BlockBilinearForm<TestPS, TrialPS>;
        using space_ptrtype = std::decay_t<SpacePtrType>;

        RowDirichletView( parent_type& parent, space_ptrtype space, int rowstart )
            :
            M_parent( parent ),
            M_space( std::move( space ) ),
            M_rowDofIdToContainerId( parent.rowDirichletDofIdToContainerId( rowstart ) )
        {
            if ( M_parent.M_matrix->staticCondensation() )
            {
                auto const maxContainerId = std::ranges::max( M_rowDofIdToContainerId );
                M_rowContainerIdToDofId.assign( maxContainerId + 1, invalid_v<size_type> );
                for ( size_type dofId = 0; dofId < M_rowDofIdToContainerId.size(); ++dofId )
                {
                    auto const containerId = M_rowDofIdToContainerId[dofId];
                    DCHECK( containerId < M_rowContainerIdToDofId.size() ) << "invalid container id";
                    M_rowContainerIdToDofId[containerId] = dofId;
                }
            }
        }

        template<typename ExprT>
        RowDirichletView& operator+=( Expr<ExprT> const& expr )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView::operator+=" );
            expr.assemble( M_space, M_space, *this );
            return *this;
        }

        condensed_matrix_type& matrix()
        {
            return *M_parent.M_matrix;
        }

        std::vector<size_type> const& dofIdToContainerIdTrial() const
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView::dofIdToContainerIdTrial" );
            return M_rowDofIdToContainerId;
        }

        std::vector<size_type> const& dofIdToContainerIdTest() const
        {
            return M_rowDofIdToContainerId;
        }

        bool shouldDeferDirichlet( Feel::Context const& on_context ) const
        {
            return M_parent.shouldDeferDirichlet( on_context );
        }

        void deferZeroRows( std::vector<int> const& dofs,
                            std::vector<value_type> const& values,
                            Feel::Context const& on_context,
                            double value_on_diagonal,
                            std::uint8_t entity_priority = vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::unspecified ) )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView::deferZeroRows" );
            if ( M_parent.M_matrix->staticCondensation() )
            {
                std::vector<int> remappedDofs;
                remappedDofs.reserve( dofs.size() );
                for ( int dof : dofs )
                {
                    CHECK( dof >= 0 && static_cast<size_type>( dof ) < M_rowContainerIdToDofId.size() )
                        << "invalid deferred Dirichlet dof " << dof;
                    auto const mappedDof = M_rowContainerIdToDofId[dof];
                    CHECK( mappedDof != invalid_v<size_type> ) << "missing inverse dof mapping for deferred Dirichlet dof " << dof;
                    remappedDofs.push_back( static_cast<int>( mappedDof ) );
                }
                M_parent.deferZeroRows( remappedDofs, values, on_context, value_on_diagonal, entity_priority );
                return;
            }

            M_parent.deferZeroRows( dofs, values, on_context, value_on_diagonal, entity_priority );
        }

        void zeroRows( std::vector<int> const& dofs,
                       Vector<value_type> const& values,
                       Vector<value_type>& rhs,
                       Feel::Context const& on_context,
                       double value_on_diagonal )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView::zeroRows" );
            M_parent.invalidateMaterializedDeferredDirichlet();
            M_parent.close();
            M_parent.M_matrix->zeroRows( dofs, values, rhs, on_context, value_on_diagonal );
        }

        void set( size_type i, size_type j, value_type const& value )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView::set" );
            M_parent.invalidateMaterializedDeferredDirichlet();
            M_parent.M_matrix->set( i, j, value );
        }

    private:
        void checkSameTestTrialFunctionSpace( char const* where ) const
        {
            M_parent.checkSameTestTrialFunctionSpace( where );
        }

        parent_type& M_parent;
        space_ptrtype M_space;
        std::vector<size_type> const& M_rowDofIdToContainerId;
        std::vector<size_type> M_rowContainerIdToDofId;
    };

    BlockBilinearForm() = default;
    BlockBilinearForm( BlockBilinearForm const& ) = default;
    BlockBilinearForm( BlockBilinearForm && ) = default;
    
    template<typename T>
        requires StaticProductSpacesType<T>
    BlockBilinearForm( T&& ps )
        :
        M_test_ps(std::forward<T>(ps)),
        M_trial_ps(M_test_ps),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED), backend(), false ) )
        {}

    template<typename T>
        requires DynamicProductSpaceType<T>
    BlockBilinearForm( T&& ps )
        :
        M_test_ps(std::forward<T>(ps)),
        M_trial_ps(M_test_ps),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED), backend(), false ) )
        {}    

    template<typename TestT, typename TrialT>
        requires ( StaticProductSpacesType<TestT> || DynamicProductSpaceType<TestT> ) &&
                 ( StaticProductSpacesType<TrialT> || DynamicProductSpaceType<TrialT> )
    BlockBilinearForm( TestT&& testPs, TrialT&& trialPs )
        :
        M_test_ps(std::forward<TestT>(testPs)),
        M_trial_ps(std::forward<TrialT>(trialPs)),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED), backend(), false ) )
        {}
    
    template<typename T,typename BackendT, typename RangeMapT = StencilRangeMap0Type>
        requires StaticProductSpacesType<T> && Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( T&& ps, BackendT&& b, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<T>(ps)),
        M_trial_ps(M_test_ps),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED, r), std::forward<BackendT>(b), false ) )
        {}

    template<typename TestT, typename TrialT, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
        requires ( StaticProductSpacesType<TestT> || DynamicProductSpaceType<TestT> ) &&
                 ( StaticProductSpacesType<TrialT> || DynamicProductSpaceType<TrialT> ) &&
                 Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( TestT&& testPs, TrialT&& trialPs, BackendT&& b, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<TestT>(testPs)),
        M_trial_ps(std::forward<TrialT>(trialPs)),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED, r), std::forward<BackendT>(b), false ) )
        {}

    template<typename T, typename BackendT, typename RangeMapT>
        requires StaticProductSpacesType<T> && Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( T&& ps, solve::strategy s, BackendT&& b, size_type pattern = Pattern::COUPLED, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<T>(ps)),
        M_trial_ps(M_test_ps),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:pattern,r),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {}
    template<typename TestT, typename TrialT, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
        requires ( StaticProductSpacesType<TestT> || DynamicProductSpaceType<TestT> ) &&
                 ( StaticProductSpacesType<TrialT> || DynamicProductSpaceType<TrialT> ) &&
                 Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( TestT&& testPs, TrialT&& trialPs, solve::strategy s, BackendT&& b, size_type pattern = Pattern::COUPLED, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<TestT>(testPs)),
        M_trial_ps(std::forward<TrialT>(trialPs)),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:pattern,r),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {
            this->checkStrategySupportsTestTrialSpaces( s, "BlockBilinearForm constructor" );
        }
    template<typename T, typename BackendT, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
        requires StaticProductSpacesType<T> && Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( T&& ps, solve::strategy s, BackendT&& b, std::vector<PatternSizeT> const& patterns, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<T>(ps)),
        M_trial_ps(M_test_ps),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?std::vector<PatternSizeT>( patterns.size(), static_cast<PatternSizeT>( Pattern::ZERO ) ):patterns, r),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {}

    template<typename TestT, typename TrialT, typename BackendT, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
        requires ( StaticProductSpacesType<TestT> || DynamicProductSpaceType<TestT> ) &&
                 ( StaticProductSpacesType<TrialT> || DynamicProductSpaceType<TrialT> ) &&
                 Feel::detail::is_backend_like_v<BackendT>
    BlockBilinearForm( TestT&& testPs, TrialT&& trialPs, solve::strategy s, BackendT&& b, std::vector<PatternSizeT> const& patterns, RangeMapT r = stencilRangeMap() )
        :
        M_test_ps(std::forward<TestT>(testPs)),
        M_trial_ps(std::forward<TrialT>(trialPs)),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?std::vector<PatternSizeT>( patterns.size(), static_cast<PatternSizeT>( Pattern::ZERO ) ):patterns, r),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {
            this->checkStrategySupportsTestTrialSpaces( s, "BlockBilinearForm constructor" );
        }

    BlockBilinearForm(product_space_t&& ps, condensed_matrix_ptrtype & m)
        :
        M_test_ps(ps),
        M_trial_ps(M_test_ps),
        M_matrix(m)
        {}

    BlockBilinearForm(product_space_t const& ps, condensed_matrix_ptrtype & m)
        :
        M_test_ps(ps),
        M_trial_ps(M_test_ps),
        M_matrix(m)
        {}
    BlockBilinearForm& operator=( BlockBilinearForm const& a )
        {
            if ( this == &a )
                return *this;

            bool same_spaces = ( M_test_ps == a.M_test_ps ) && ( M_trial_ps == a.M_trial_ps );
            M_test_ps = a.M_test_ps;
            M_trial_ps = a.M_trial_ps;
            if ( !this->isMatrixAllocated() || !same_spaces )
            {
                this->allocateMatrix( a.M_matrix->solveStrategy(), a.M_matrix->backend() );
                M_matrix->setBackend( a.M_matrix->backend()->clone() );
            }
            M_matrix->zero();
            M_matrix->addMatrix( 1.,(MatrixSparse<value_type> const&)*a.M_matrix->getSparseMatrix() );
            M_pendingDirichlet = a.M_pendingDirichlet;
            M_appliedDirichlet = a.M_appliedDirichlet;
            M_constrainedMatrix = a.M_constrainedMatrix;
            M_dirichletPolicy = a.M_dirichletPolicy;
            
            return *this;
        }
    BlockBilinearForm& operator=( BlockBilinearForm && a ) = default;
    BlockBilinearForm& operator+=( BlockBilinearForm& a )
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

    //!
    //! allocate algebraic representation of the bilinear form
    //! @param s the solve strategy (monolithic, static condensation or local)
    //! @param b the algebraic backend to use (petsc or eigen)
    //!
    void allocateMatrix( solve::strategy s = solve::strategy::monolithic, backend_ptrtype const& b = backend() )
        {
            this->checkStrategySupportsTestTrialSpaces( s, "allocateMatrix" );
            M_matrix = std::make_shared<condensed_matrix_type>( s, csrGraphBlocks(M_test_ps, M_trial_ps, (s==solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), b, (s==solve::strategy::static_condensation)?false:true );
            this->clearDeferredDirichlet();
        }
    //!
    //! @return true if allocated, false otherwise
    //!
    bool isMatrixAllocated() const
        {
            return (bool)M_matrix;
        }
#if 0
    template<typename N1,typename N2>
    decltype(auto) operator()( N1 n1, N2 n2 )
        {
            cout << "filling out matrix block (" << n1 << "," << n2 << ")\n";

        }
#endif
    template<typename N1, typename N2>
    decltype(auto) operator()( N1 n1, N2 n2, int s1 = 0, int s2 = 0 )
        {
            this->invalidateMaterializedDeferredDirichlet();
            auto&& testSpaces = remove_shared_ptr_f( M_test_ps );
            auto&& trialSpaces = remove_shared_ptr_f( M_trial_ps );
            #if 0
            hana::if_( hana::bool_<Feel::is_shared_ptr_v<PS>>{},
                                     []( auto&& x ) { return *x; },
                                     []( auto&& x ) { return x; } )(M_test_ps);
            #endif
            int const rowIndex = Feel::detail::blockformFlattenedBlockIndex( testSpaces, n1, s1 );
            int const colIndex = Feel::detail::blockformFlattenedBlockIndex( trialSpaces, n2, s2 );
            auto test_space = hana::at( testSpaces.tupleSpaces(), n1 );
            auto trial_space = hana::at( trialSpaces.tupleSpaces(), n2 );

            return hana::eval_if(std::is_base_of<ProductSpaceBase,decay_type<decltype(test_space)>>{},
                                 [&]( auto _ ) { return hana::eval_if( std::is_base_of<ProductSpaceBase,decay_type<decltype(trial_space)>>{},
                                                                 [&] (auto _) {
                                                                     LOG(INFO) << "filling out dyn matrix block (" << rowIndex << "," << colIndex  << ")\n";
                                                                     return form2(_test=(*_(test_space))[s1],_trial=(*_(trial_space))[s2],
                                                                                  _name="bilinearform.a"s+"("+std::to_string(rowIndex)+","s+std::to_string(colIndex)+")"s,
                                                                                  _matrix=M_matrix->block(rowIndex,colIndex), _rowstart=rowIndex, _colstart=colIndex );
                                                                 },
                                                                 [&] (auto _){
                                                                     LOG(INFO) << "filling out dyn matrix block (" << rowIndex << "," << colIndex  << ")\n";
                                                                     return form2(_test=(*_(test_space))[s1],_trial=_(trial_space),
                                                                                  _name="bilinearform.a"s+"("+std::to_string(rowIndex)+","s+std::to_string(colIndex)+")"s,
                                                                                  _matrix=M_matrix->block(rowIndex,colIndex), _rowstart=rowIndex, _colstart=colIndex );
                                                                 }); },
                                 [&]( auto _ ) { return hana::eval_if( std::is_base_of<ProductSpaceBase,decay_type<decltype(trial_space)>>{},
                                                                 [&] (auto _) {
                                                                     LOG(INFO) << "filling out dyn matrix block (" << rowIndex << "," << colIndex  << ")\n";
                                                                     return form2(_test=_(test_space),_trial=(*_(trial_space))[s2],
                                                                                  _name="bilinearform.a"s+"("+std::to_string(rowIndex)+","s+std::to_string(colIndex)+")"s,
                                                                                  _matrix=M_matrix->block(rowIndex,colIndex), _rowstart=rowIndex, _colstart=colIndex );
                                                                 },
                                                                 [&] (auto _){
                                                                     LOG(INFO) << "filling out dyn matrix block (" << rowIndex << "," << colIndex  << ")\n";
                                                                     return form2(_test=_(test_space),_trial=_(trial_space),
                                                                                  _name="bilinearform.a"s+"("+std::to_string(rowIndex)+","s+std::to_string(colIndex)+")"s,
                                                                                  _matrix=M_matrix->block(rowIndex,colIndex), _rowstart=rowIndex, _colstart=colIndex );
                                                                 }); });


        }

    template<typename N1>
    decltype(auto) row( N1 n1, int s1 = 0 )
        requires Feel::detail::BlockformStaticProductSpaces<test_product_space_t>
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::row" );
            this->invalidateMaterializedDeferredDirichlet();
            auto&& spaces = remove_shared_ptr_f( M_test_ps );
            auto space = hana::at( spaces.tupleSpaces(), n1 );
            int const rowIndex = Feel::detail::blockformFlattenedBlockIndex( spaces, n1, s1 );

            if constexpr ( std::is_base_of_v<ProductSpaceBase, decay_type<decltype(space)>> )
            {
                auto rowSpace = (*space)[s1];
                return RowDirichletView<decltype( rowSpace )>( *this, rowSpace, rowIndex );
            }
            else
            {
                return RowDirichletView<decltype( space )>( *this, space, rowIndex );
            }
        }

    decltype(auto) operator()( int n1, int n2 )
        requires DynamicProductSpaceType<test_product_space_t> && DynamicProductSpaceType<trial_product_space_t>
        {
            this->invalidateMaterializedDeferredDirichlet();
            cout << "filling out matrix block (" << n1 << "," << n2 << ")\n";
            return form2(_test=M_test_ps[n1],_trial=M_trial_ps[n2], _matrix=M_matrix, _rowstart=int(n1), _colstart=int(n2) );
        }

    template<typename T>
    void setFunctionSpace( T&& ps )
        {
            M_test_ps = std::forward<T>(ps);
            M_trial_ps = M_test_ps;
        }
    template<typename TestT, typename TrialT>
    void setFunctionSpaces( TestT&& testPs, TrialT&& trialPs )
        {
            M_test_ps = std::forward<TestT>(testPs);
            M_trial_ps = std::forward<TrialT>(trialPs);
        }
    template<typename BackendT>
    void setStrategy( BackendT&& b )
        {
            
            M_matrix = std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_test_ps, M_trial_ps, Pattern::COUPLED), std::forward<BackendT>(b), false );
            this->clearDeferredDirichlet();
        }
    template<typename BackendT>
    void setStrategy( solve::strategy s, BackendT&& b, size_type pattern = Pattern::COUPLED )
        {
            this->checkStrategySupportsTestTrialSpaces( s, "setStrategy" );
            M_matrix =  std::make_shared<condensed_matrix_type>( s,
                                                                   csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:pattern),
                                                                   std::forward<BackendT>(b),
                                                                   (s>=solve::strategy::static_condensation)?false:true );
            this->clearDeferredDirichlet();
        }
    template<typename BackendT, typename PatternSizeT>
    void setStrategy( solve::strategy s, BackendT&& b, std::vector<PatternSizeT> const& patterns )
        {
            this->checkStrategySupportsTestTrialSpaces( s, "setStrategy" );
            M_matrix =  std::make_shared<condensed_matrix_type>( s,
                                                                   csrGraphBlocks(M_test_ps, M_trial_ps, (s>=solve::strategy::static_condensation)?std::vector<PatternSizeT>( patterns.size(), static_cast<PatternSizeT>( Pattern::ZERO ) ):patterns),
                                                                   std::forward<BackendT>(b),
                                                                   (s>=solve::strategy::static_condensation)?false:true );
            this->clearDeferredDirichlet();
        }
    // Close the assembled base block operator only. Deferred Dirichlet
    // constraints stay in form state until materialized explicitly.
    void close()
        {
            M_matrix->close();
        }
    void closeBaseOperator()
        {
            this->close();
        }
    bool baseOperatorClosed() const noexcept
        {
            return M_matrix->closed();
        }
    void zero()
        {
            M_matrix->zero();
            this->clearDeferredDirichlet();
        }
    void zero(int n1, int n2 )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->zero( n1, n2 );
        }
    void transpose(int n1, int n2 )
        {
            this->invalidateMaterializedDeferredDirichlet();
            M_matrix->transposeBlock( n1, n2 );
        }
    //!
    //! @return the number of non-zero entries in matrix representation
    //!
    std::size_t nnz() const { return M_matrix->nnz(); }
    
    void syncLocalMatrix()
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::syncLocalMatrix" );
            int s = M_test_ps.numberOfSpaces();
            int n = 0;
            auto pst = M_test_ps.tupleSpaces();
            auto cp = hana::cartesian_product( hana::make_tuple( pst, pst ) );
            int nstatic = hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(hana::back(pst))>>{},
                                    [s] (auto&& x ) { return s-hana::back(std::forward<decltype(x)>(x))->numberOfSpaces()+1; },
                                    [s] (auto&& x ) { return s; } )( pst );
            hana::for_each( cp, [&]( auto const& e )
                            {
                                int r = n/nstatic;
                                int c = n%nstatic;
                                auto test_space = e[0_c];
                                auto trial_space = e[1_c];
                                DVLOG( 1 ) << "syncLocalMatrix("<< r << ","<< c <<")dim "<< test_space->mesh()->dimension()<<" , " << trial_space->mesh()->dimension() <<"\n";
                                M_matrix->sc(r,c)->syncLocalMatrix( test_space,trial_space );
                                ++n;
                            });
        }
    void deferZeroRows( std::vector<int> const& dofs,
                        std::vector<value_type> const& values,
                        Feel::Context const& on_context,
                        double value_on_diagonal,
                        std::uint8_t entity_priority = vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::unspecified ) )
        {
            M_pendingDirichlet.append( dofs, values, on_context, value_on_diagonal, entity_priority );
            this->invalidateMaterializedDeferredDirichlet();
        }
    vf::DeferredDirichletPolicy dirichletPolicy() const noexcept
        {
            return M_dirichletPolicy;
        }
    void setDirichletPolicy( vf::DeferredDirichletPolicy policy ) noexcept
        {
            M_dirichletPolicy = policy;
        }
    bool useDeferredDirichlet() const noexcept
        {
            return vf::usesDeferredDirichlet( this->dirichletPolicy() );
        }
    void setUseDeferredDirichlet( bool value ) noexcept
        {
            this->setDirichletPolicy( value ? vf::DeferredDirichletPolicy::deferred :
                                             vf::DeferredDirichletPolicy::immediate );
        }
    BlockBilinearForm& deferDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::deferred );
            return *this;
        }
    BlockBilinearForm& autoDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::automatic );
            return *this;
        }
    BlockBilinearForm& immediateDirichlet() noexcept
        {
            this->setDirichletPolicy( vf::DeferredDirichletPolicy::immediate );
            return *this;
        }
    bool shouldDeferDirichlet( Feel::Context const& on_context ) const noexcept
        {
            return vf::shouldDeferDirichlet( this->dirichletPolicy(),
                                             on_context,
                                             static_cast<bool>( M_matrix ) && !M_matrix->localSolve() );
        }
    bool hasDeferredDirichlet() const noexcept
        {
            return !M_pendingDirichlet.empty();
        }
    bool hasPendingDirichletConstraints() const noexcept
        {
            return this->hasDeferredDirichlet();
        }
    bool hasDirichletConstraints() const noexcept
        {
            return !M_pendingDirichlet.empty() || !M_appliedDirichlet.empty();
        }
    bool supportsConstrainedOperatorView() const noexcept
        {
            return static_cast<bool>( M_matrix ) && M_matrix->monolithic();
        }
    bool hasMaterializedConstrainedOperator() const noexcept
        {
            return static_cast<bool>( M_constrainedMatrix );
        }
    void clearDeferredDirichlet() noexcept
        {
            M_pendingDirichlet.clear();
            M_appliedDirichlet.clear();
            M_constrainedMatrix.reset();
            M_constrainedVector.reset();
            M_constrainedVectorSource = nullptr;
            M_constrainedVectorSourceRevision = 0;
        }
    template<typename CondensedFormT, typename CondensedRhsT>
    void applyDeferredDirichlet( CondensedFormT& condensedForm, CondensedRhsT& condensedRhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;

            auto [matrix, rhsVector] = this->closeCondensedSystem( condensedForm, condensedRhs );
            this->applyDeferredDirichletToClosedCondensedSystem( condensedForm, matrix, rhsVector );
        }
    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    void applyDeferredDirichlet( RhsType& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;

            auto constrainedVector = this->constrainedVectorPtr( rhs );
            auto rhsVectorHandle = rhs.vectorPtr();
            vf::copyVectorValues( rhsVectorHandle->getVector(), constrainedVector );
            M_constrainedVectorSource = static_cast<void const*>( rhsVectorHandle->getVector().get() );
            M_constrainedVectorSourceRevision = rhsVectorHandle->getVector()->revision();
        }
    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    void applyDeferredDirichlet( RhsType const& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;

            auto rhsVectorHandle = rhs.vectorPtr();
            auto constrainedVector = this->constrainedVectorPtr( rhs );
            vf::copyVectorValues( rhsVectorHandle->getVector(), constrainedVector );
            M_constrainedVectorSource = static_cast<void const*>( rhsVectorHandle->getVector().get() );
            M_constrainedVectorSourceRevision = rhsVectorHandle->getVector()->revision();
        }
    sparse_matrix_ptrtype baseMatrixPtr() const
        {
            return M_matrix->getSparseMatrix();
        }
    sparse_matrix_ptrtype baseMatrixPtr()
        {
            this->invalidateMaterializedDeferredDirichlet();
            return M_matrix->getSparseMatrix();
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto constrainedSystem( RhsType& rhs )
        {
            return this->materializeConstrainedSystem( rhs.vectorPtr()->getVector() );
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto constrainedSystem( RhsType const& rhs )
        {
            return const_cast<BlockBilinearForm*>( this )->materializeConstrainedSystem( rhs.vectorPtr()->getVector() );
        }

    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    auto activeSystem( RhsType& rhs )
        {
            CHECK( this->supportsConstrainedOperatorView() || !this->hasDirichletConstraints() )
                << "activeSystem() is only available for monolithic block solves when Dirichlet constraints are present";
            if ( this->hasDirichletConstraints() )
                return this->constrainedSystem( rhs );
            return std::pair{ this->baseMatrixPtr(), rhs.vectorPtr()->getVector() };
        }

    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    auto activeSystem( RhsType const& rhs )
        {
            CHECK( this->supportsConstrainedOperatorView() || !this->hasDirichletConstraints() )
                << "activeSystem() is only available for monolithic block solves when Dirichlet constraints are present";
            if ( this->hasDirichletConstraints() )
                return const_cast<BlockBilinearForm*>( this )->materializeConstrainedSystem( rhs.vectorPtr()->getVector() );
            return std::pair{ this->baseMatrixPtr(), rhs.vectorPtr()->getVector() };
        }

    sparse_matrix_ptrtype constrainedMatrixPtr()
        {
            return this->materializeConstrainedMatrix();
        }
    sparse_matrix_ptrtype constrainedMatrixPtr() const
        {
            return const_cast<BlockBilinearForm*>( this )->constrainedMatrixPtr();
        }
    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    sparse_matrix_ptrtype constrainedMatrixPtr( RhsType& rhs )
        {
            return this->constrainedSystem( rhs ).first;
        }
    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    sparse_matrix_ptrtype constrainedMatrixPtr( RhsType const& rhs )
        {
            return const_cast<BlockBilinearForm*>( this )->constrainedMatrixPtr( rhs );
        }
    sparse_matrix_ptrtype activeMatrixPtr()
        {
            CHECK( this->supportsConstrainedOperatorView() || !this->hasDirichletConstraints() )
                << "activeMatrixPtr() is only available for monolithic block solves when Dirichlet constraints are present";
            return this->hasDirichletConstraints() ? this->constrainedMatrixPtr() : this->baseMatrixPtr();
        }
    sparse_matrix_ptrtype activeMatrixPtr() const
        {
            CHECK( this->supportsConstrainedOperatorView() || !this->hasDirichletConstraints() )
                << "activeMatrixPtr() is only available for monolithic block solves when Dirichlet constraints are present";
            return this->hasDirichletConstraints() ? this->constrainedMatrixPtr() : this->baseMatrixPtr();
        }
    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    sparse_matrix_ptrtype activeMatrixPtr( RhsType& rhs )
        {
            return this->activeSystem( rhs ).first;
        }
    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    sparse_matrix_ptrtype activeMatrixPtr( RhsType const& rhs )
        {
            return const_cast<BlockBilinearForm*>( this )->activeMatrixPtr( rhs );
        }
    void materializeConstrainedOperator()
        {
            if ( !this->hasDirichletConstraints() )
                return;
            CHECK( this->supportsConstrainedOperatorView() )
                << "materializeConstrainedOperator() is only available for monolithic block solves";
            (void)this->constrainedMatrixPtr();
        }
    template<typename RhsType>
        requires requires( RhsType& rhs ) { rhs.vectorPtr(); }
    void materializeConstrainedOperator( RhsType& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;
            CHECK( this->supportsConstrainedOperatorView() )
                << "materializeConstrainedOperator() is only available for monolithic block solves";
            (void)this->constrainedSystem( rhs );
        }
    template<typename RhsType>
        requires requires( RhsType const& rhs ) { rhs.vectorPtr(); }
    void materializeConstrainedOperator( RhsType const& rhs )
        {
            if ( !this->hasDirichletConstraints() )
                return;
            CHECK( this->supportsConstrainedOperatorView() )
                << "materializeConstrainedOperator() is only available for monolithic block solves";
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
    sparse_matrix_ptrtype matrixPtr() const { return M_matrix; }
    sparse_matrix_ptrtype matrixPtr() { this->invalidateMaterializedDeferredDirichlet(); return M_matrix; }
    condensed_matrix_type const& matrix() const { return *M_matrix; }
    condensed_matrix_type& matrix() { this->invalidateMaterializedDeferredDirichlet(); return *M_matrix; }
    test_product_space_t testFunctionSpace() const { return M_test_ps; }
    trial_product_space_t trialFunctionSpace() const { return M_trial_ps; }
    product_space_t functionSpace() const { return M_test_ps; }
    bool isRectangular() const { return !this->isSquareBlockForm(); }
    auto l1Norm() const { return M_matrix->l1Norm(); }
    auto linftyNorm() const { return M_matrix->linftyNorm(); }
    using pre_solve_type = typename Backend<value_type>::pre_solve_type;
    using post_solve_type = typename Backend<value_type>::post_solve_type;

    template <typename ... Ts>
    typename Backend<double>::solve_return_type solve( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            auto && solution = args.get(_solution);
            auto && rhs = args.get(_rhs);
            bool condense = args.get_else(_condense,false);
            auto && condenser = args.get_else(_condenser,condenser_poisson() );
            bool local = args.get_else(_local,false);
            std::string const& name = args.get_else(_name,"" );
            std::string const& kind = args.get_else_invocable(_kind,[&name](){ return soption(_prefix=name,_name="backend"); } );
            bool rebuild = args.get_else_invocable(_rebuild,[&name](){ return boption(_prefix=name,_name="backend.rebuild"); } );
            pre_solve_type pre = args.get_else(_pre, pre_solve_type() );
            post_solve_type post = args.get_else(_post,post_solve_type() );

            this->checkSquareBlockForm( "BlockBilinearForm::solve" );

            if constexpr ( StokesCondenserTag<decltype( condenser )> )
                return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
            else
            {
                if ( condense )
                {
                    this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solve(_condense=true)" );
                    return solveImplCondense( M_test_ps, solution, rhs, name, kind, rebuild, pre, post, condenser );
                }
                if ( local )
                {
                    this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solve(_local=true)" );
                    return solveImplLocal( M_test_ps, solution, rhs, name, kind, rebuild, pre, post );
                }
                return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
            }
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplLocal( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                    bool rebuild, pre_solve_type pre, post_solve_type post )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solveImplLocal" );
            auto sc = M_matrix->sc();
            tic();
            cout << " . starting local Solve" << std::endl;
            auto vsc = rhs.vectorPtr()->sc();
            sc->localSolve ( vsc, solution);
            cout << " . local Solve done" << std::endl;
            toc("blockform.local.localsolve",Environment::logVerbosityLevel()>0);
            typename Backend<double>::solve_return_type r;
            return r;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires DynamicProductSpaceType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplLocal( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                    bool rebuild, pre_solve_type pre, post_solve_type post )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solveImplLocal" );
            typename Backend<double>::solve_return_type r;
            return r;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t, typename CT>
        requires DynamicProductSpaceType<PS_t> && CondenserTag<CT>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post, CT ct )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solveImplCondense" );
            if constexpr ( Sb9CondenserTag<CT> )
            {
                return solveImplCondenseTwoFieldInternal( ps, solution, rhs, name, kind, rebuild, pre, post );
            }
            return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplCondenseTwoFieldInternal( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                                       bool rebuild, pre_solve_type pre, post_solve_type post )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solveImplCondenseTwoFieldInternal" );
            static_assert( Sb9CondensableProductElement<Solution_t>,
                           "SB9 static condensation expects a 2-field product element with plain local field blocks" );
            CHECK( std::decay_t<Solution_t>::nspaces == 2 ) << "SB9 static condensation expects a 2-field product space";

            auto& eu = solution( 0_c );
            auto sc = M_matrix->sc();

            tic();
            auto psS = product( eu.functionSpace() );
            auto S = blockform2( psS, solve::strategy::monolithic, backend() );
            auto V = blockform1( psS, solve::strategy::monolithic, backend() );
            toc( "blockform.sc.space", Environment::logVerbosityLevel() > 0 );

            tic();
            this->syncLocalMatrix();
            sc->condense( rhs.vectorPtr()->sc(), solution, S, V );
            toc( "blockform.sc.condense", Environment::logVerbosityLevel() > 0 );

            this->prepareCondensedSystemForSolve( S, V );

            tic();
            auto U = psS.element();
            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm( name, "sc" ), _kind=kind,
                              _rebuild=rebuild, _pre=pre, _post=post );
            toc( "blockform.sc.solve", Environment::logVerbosityLevel() > 0 );

            solution( 0_c ) = U( 0_c );

            tic();
            sc->localSolve( rhs.vectorPtr()->sc(), solution );
            toc( "blockform.sc.localsolve", Environment::logVerbosityLevel() > 0 );
            return r;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t, typename CT>
        requires FoldableProductSpacesType<PS_t> && CondenserTag<CT>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post, CT ct )
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::solveImplCondense" );
            if constexpr ( Sb9CondenserTag<CT> )
            {
                return solveImplCondenseTwoFieldInternal( ps, solution, rhs, name, kind, rebuild, pre, post );
            }
            return solveImplCondense( ps, solution, rhs, name, kind, rebuild, pre, post,hana::integral_constant<int,decltype(hana::size( M_test_ps.tupleSpaces() ))::value>() );
        }
    template <typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImpl( Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind = "petsc",
               bool rebuild = false, pre_solve_type pre = pre_solve_type(), post_solve_type post = post_solve_type() )
        {
            if ( !this->hasDirichletConstraints() )
            {
                this->closeBaseOperator();
                rhs.vectorPtr()->close();
            }
            auto [matrix, rhsVector] = this->activeSystem( rhs );
            auto U = backend()->newBlockVector(_block=solution, _copy_values=false);
            auto solveBackend = backend( _name=name, _kind=kind, _rebuild=rebuild,
                                         _worldcomm=Environment::worldCommPtr() );
            tic();
            auto r1 = solveBackend->solve( _matrix=matrix,
                                           _auxiliary_matrix=this->baseMatrixPtr(),
                                           _rhs=rhsVector,
                                           _solution=U,
                                           _pre=pre,
                                           _post=post
                                           );
            toc("blockform.monolithic",Environment::logVerbosityLevel()>0);
            if ( Environment::isSequential() && boption("exporter.matlab") )
            {
                M_matrix->getSparseMatrix()->printMatlab("A.m");
                rhs.vectorPtr()->getVector()->printMatlab("b.m");
            } 
            solution.localize(U);
            return r1;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,1> )
        {
            return typename Backend<double>::solve_return_type{};
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,2> )
        {
            return typename Backend<double>::solve_return_type{};
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,3> )
        {
#if 1
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();
            tic();
            auto psS = product( e3.functionSpace() );
            toc("blockform.sc.space",Environment::logVerbosityLevel()>0);
            tic();
            auto S = blockform2( psS, solve::strategy::monolithic, backend(), Pattern::HDG  );
            toc("blockform.sc.bilinearform",Environment::logVerbosityLevel()>0);
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto V = blockform1( psS, solve::strategy::monolithic, backend() );
            
            tic();
            this->syncLocalMatrix();
            toc("blockform.sc.sync", Environment::logVerbosityLevel()>0);
            tic();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", Environment::logVerbosityLevel()>0);
            this->prepareCondensedSystemForSolve( S, V );
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;
            auto U = psS.element();

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );
            //auto r = backend(_name=prefixvm(name,"sc"),_rebuild=rebuild)->solve( _matrix=S.matrixPtr(), _rhs=V.vectorPtr(), _solution=e3);
            solution(2_c)=U(0_c);
            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", Environment::logVerbosityLevel()>0);

#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->localSolve ( rhs.vectorPtr()->sc(), solution);
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",Environment::logVerbosityLevel()>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m");
#endif


            return r;
#else
            return {};
#endif
        }

    //!
    //! solve using static condensation in the case of 2 trace spaces
    //!
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,4> )
        {
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();

            auto Th = product2( M_test_ps[3_c], M_test_ps[2_c] );
            auto S = blockform2(Th, solve::strategy::monolithic, backend(), Pattern::HDG);
            auto V = blockform1(Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", Environment::logVerbosityLevel()>0);
            this->prepareCondensedSystemForSolve( S, V );
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );

            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", Environment::logVerbosityLevel()>0);

            solution(2_c)=U(0_c);
            for( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution(3_c,i)=U(1_c,i);
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_test_ps[3_c]->numberOfSpaces());
            sc->localSolve ( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",Environment::logVerbosityLevel()>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m"); 
#endif

            return r;
        }
 
    //!
    //! solve using static condensation in the case of 2 trace spaces
    //!
    template <typename PS_t, typename Solution_t, typename Rhs_t>
        requires StaticProductSpacesType<PS_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,5> )
        {
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();

            auto Th = product2( M_test_ps[3_c], M_test_ps[4_c], ps[2_c] );
            std::vector<size_type> patterns = {Pattern::HDG,Pattern::HDG,Pattern::ZERO,
                                               Pattern::HDG,Pattern::HDG,Pattern::COUPLED,
                                               Pattern::ZERO,Pattern::COUPLED,Pattern::COUPLED};
            auto S = blockform2(Th, solve::strategy::monolithic, backend(), patterns);
            auto V = blockform1(Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", Environment::logVerbosityLevel()>0);
            this->prepareCondensedSystemForSolve( S, V );
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );

            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", Environment::logVerbosityLevel()>0);

            solution(2_c)=U(0_c);
            for( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution(3_c,i)=U(1_c,i);
            for( int i = 0; i < Th[2_c]->numberOfSpaces(); ++i )
                solution(4_c,i) = U(2_c,i);
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_test_ps[3_c]->numberOfSpaces());
            sc->localSolve ( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",Environment::logVerbosityLevel()>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m"); 
#endif

            return r;
        }

        //!
        //! solve using static condensation in the case of 2 trace spaces
        //!
        template <typename PS_t, typename Solution_t, typename Rhs_t>
            requires StaticProductSpacesType<PS_t>
        typename Backend<double>::solve_return_type
        solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                           bool rebuild, pre_solve_type pre, post_solve_type post, hana::integral_constant<int, 3*5> )
        {
#if 0            
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution( 2_c );
            auto& e2 = solution( 1_c );
            auto& e1 = solution( 0_c );

            auto sc = M_matrix->sc();

            auto Th = product2( M_test_ps[3_c], M_test_ps[4_c], ps[2_c] );
            std::vector<size_type> patterns = { Pattern::HDG, Pattern::HDG, Pattern::ZERO,
                                                Pattern::HDG, Pattern::HDG, Pattern::COUPLED,
                                                Pattern::ZERO, Pattern::COUPLED, Pattern::COUPLED };
            auto S = blockform2( Th, solve::strategy::monolithic, backend(), patterns );
            auto V = blockform1( Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense( rhs.vectorPtr()->sc(), solution, S, V );
            toc( "blockform.sc.condense", Environment::logVerbosityLevel() > 0 );
            this->prepareCondensedSystemForSolve( S, V );
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution = U, _rhs = V, _name = prefixvm( name, "sc" ), _rebuild = rebuild ); //, _condense=true );

            cout << " . Solve done" << std::endl;
            toc( "blockform.sc.solve", Environment::logVerbosityLevel() > 0 );

            solution( 2_c ) = U( 0_c );
            for ( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution( 3_c, i ) = U( 1_c, i );
            for ( int i = 0; i < Th[2_c]->numberOfSpaces(); ++i )
                solution( 4_c, i ) = U( 2_c, i );
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_test_ps[3_c]->numberOfSpaces() );
            sc->localSolve( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc( "blockform.sc.localsolve", Environment::logVerbosityLevel() > 0 );
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m");
#endif

            return r;
#else
            return typename Backend<double>::solve_return_type{};
#endif
        }
private:
        template<typename CondensedFormT, typename CondensedRhsT>
        auto closeCondensedSystem( CondensedFormT& condensedForm, CondensedRhsT& condensedRhs )
        {
            condensedForm.close();
            condensedRhs.close();
            return std::pair{ condensedForm.matrixPtr(), condensedRhs.vectorPtr()->getVector() };
        }

        template<typename CondensedFormT>
        void applyDeferredDirichletToClosedCondensedSystem( CondensedFormT& condensedForm,
                                                            sparse_matrix_ptrtype const& matrix,
                                                            vector_ptrtype const& rhsVector )
        {
            auto const deferredEntries = this->allDeferredDirichletConstraints().entries();
            vf::applyDeferredDirichletEntries( deferredEntries, matrix, rhsVector );
            this->promotePendingDeferredDirichlet();
            condensedForm.close();
            if ( !rhsVector->closed() )
                rhsVector->close();
        }

        template<typename CondensedFormT, typename CondensedRhsT>
        auto prepareCondensedSystemForSolve( CondensedFormT& condensedForm, CondensedRhsT& condensedRhs )
        {
            auto [matrix, rhsVector] = this->closeCondensedSystem( condensedForm, condensedRhs );
            if ( this->hasDirichletConstraints() )
                this->applyDeferredDirichletToClosedCondensedSystem( condensedForm, matrix, rhsVector );
            return std::pair{ matrix, rhsVector };
        }

        void invalidateMaterializedDeferredDirichlet() noexcept
        {
            M_constrainedMatrix.reset();
            M_constrainedVector.reset();
            M_constrainedVectorSource = nullptr;
            M_constrainedVectorSourceRevision = 0;
        }

        deferred_dirichlet_set_type allDeferredDirichletConstraints() const
        {
            deferred_dirichlet_set_type constraints;
            constraints.append( M_appliedDirichlet );
            constraints.append( M_pendingDirichlet );
            return constraints;
        }

        void promotePendingDeferredDirichlet()
        {
            M_appliedDirichlet.append( M_pendingDirichlet );
            M_pendingDirichlet.clear();
        }

        template<typename VectorPtrType>
        auto materializeConstrainedSystem( VectorPtrType const& rhsVector )
        {
            CHECK( M_matrix->monolithic() ) << "constrainedVectorPtr() is only available for monolithic block solves";
            if ( !this->hasDirichletConstraints() )
                return std::pair{ this->baseMatrixPtr(), rhsVector };

            auto const rhsSource = static_cast<void const*>( rhsVector.get() );
            auto const rhsRevision = rhsVector->revision();
            if ( M_constrainedMatrix &&
                 M_constrainedVector &&
                 M_pendingDirichlet.empty() &&
                 M_constrainedVectorSource == rhsSource &&
                 M_constrainedVectorSourceRevision == rhsRevision )
            {
                return std::pair{ M_constrainedMatrix, M_constrainedVector };
            }

            this->close();
            if ( !rhsVector->closed() )
                rhsVector->close();

            auto constrainedMatrix = this->baseMatrixPtr()->clone();
            auto constrainedVector = vf::cloneVectorWithValues( rhsVector );
            auto const deferredEntries = this->allDeferredDirichletConstraints().entries();
            vf::applyDeferredDirichletEntries( deferredEntries, constrainedMatrix, constrainedVector );

            constrainedMatrix->close();
            if ( !constrainedVector->closed() )
                constrainedVector->close();

            M_constrainedMatrix = constrainedMatrix;
            M_constrainedVector = constrainedVector;
            M_constrainedVectorSource = rhsSource;
            M_constrainedVectorSourceRevision = rhsRevision;
            this->promotePendingDeferredDirichlet();
            return std::pair{ M_constrainedMatrix, M_constrainedVector };
        }

        sparse_matrix_ptrtype materializeConstrainedMatrix()
        {
            CHECK( M_matrix->monolithic() ) << "constrainedMatrixPtr() is only available for monolithic block solves";
            if ( !this->hasDirichletConstraints() )
                return this->baseMatrixPtr();
            if ( M_constrainedMatrix && M_pendingDirichlet.empty() )
                return M_constrainedMatrix;

            this->close();
            auto constrainedMatrix = this->baseMatrixPtr()->clone();
            auto dummyRhs = Feel::backend( _worldcomm=Environment::worldCommPtr() )->newVector( this->baseMatrixPtr()->mapRowPtr() );
            dummyRhs->zero();
            dummyRhs->close();
            auto const deferredEntries = this->allDeferredDirichletConstraints().entries();
            vf::applyDeferredDirichletEntries( deferredEntries, constrainedMatrix, dummyRhs );
            constrainedMatrix->close();

            M_constrainedMatrix = constrainedMatrix;
            M_constrainedVector.reset();
            M_constrainedVectorSource = nullptr;
            M_constrainedVectorSourceRevision = 0;
            this->promotePendingDeferredDirichlet();
            return M_constrainedMatrix;
        }

        template<typename VectorPtrType>
        auto materializeConstrainedVector( VectorPtrType const& rhsVector )
        {
            if ( !this->hasDirichletConstraints() )
                return rhsVector;
            auto [constrainedMatrix, constrainedVector] = this->materializeConstrainedSystem( rhsVector );
            return constrainedVector;
        }

        bool isSquareBlockForm() const
        {
            if ( M_test_ps.nDof() != M_trial_ps.nDof() )
                return false;
            if ( M_matrix )
                return M_matrix->mapRow().nDof() == M_matrix->mapCol().nDof();
            return true;
        }

        bool hasSameTestTrialFunctionSpace() const
        {
            if constexpr ( std::is_same_v<test_product_space_t, trial_product_space_t> )
                return M_test_ps == M_trial_ps;
            else
                return false;
        }

        void checkSquareBlockForm( char const* where ) const
        {
            CHECK( this->isSquareBlockForm() )
                << where << " requires equal row and column dof counts";
        }

        void checkSameTestTrialFunctionSpace( char const* where ) const
        {
            CHECK( this->hasSameTestTrialFunctionSpace() )
                << where << " requires identical test/trial product spaces for now; matching dof counts are not enough";
        }

        std::vector<size_type> const& rowDirichletDofIdToContainerId( int rowstart ) const
        {
            this->checkSameTestTrialFunctionSpace( "BlockBilinearForm::RowDirichletView" );
            return M_matrix->mapRowPtr()->dofIdToContainerId( rowstart );
        }

        void checkStrategySupportsTestTrialSpaces( solve::strategy s, char const* where ) const
        {
            CHECK( s == solve::strategy::monolithic || this->hasSameTestTrialFunctionSpace() )
                << where << " supports non-monolithic strategies only with identical test/trial product spaces for now";
        }

        deferred_dirichlet_set_type M_pendingDirichlet;
        deferred_dirichlet_set_type M_appliedDirichlet;
        sparse_matrix_ptrtype M_constrainedMatrix;
        vector_ptrtype M_constrainedVector;
        void const* M_constrainedVectorSource = nullptr;
        std::size_t M_constrainedVectorSourceRevision = 0;
        vf::DeferredDirichletPolicy M_dirichletPolicy = vf::DeferredDirichletPolicy::automatic;
        test_product_space_t M_test_ps;
        trial_product_space_t M_trial_ps;
        condensed_matrix_ptrtype M_matrix;
};

template<typename PS>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps) );
}

template<typename PS, typename BackendT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, BackendT&& b )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps), std::forward<BackendT>(b) );
}

template<typename PS, typename BackendT, typename RangeMapT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS && ps, solve::strategy s, BackendT&& b, size_type pattern, RangeMapT r )
{
    return BlockBilinearForm<PS>( std::forward<PS>( ps ), s, std::forward<BackendT>(b), pattern, r );
}

template<typename PS, typename BackendT, typename PatternSizeT, typename RangeMapT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS && ps, solve::strategy s, BackendT&& b, std::vector<PatternSizeT> const& patterns, RangeMapT r )
{
    return BlockBilinearForm<PS>( std::forward<PS>( ps ), s, std::forward<BackendT>(b), patterns, r );
}

template<typename PS, typename BackendT, typename RangeMapT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b, size_type pattern, RangeMapT r )
{
    return BlockBilinearForm<PS>( ps, s, std::forward<BackendT>(b), pattern, r );
}

template<typename PS, typename BackendT, typename PatternSizeT, typename RangeMapT>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b, std::vector<PatternSizeT> const& patterns, RangeMapT r )
{
    return BlockBilinearForm<PS>( ps, s, std::forward<BackendT>(b), patterns, r );
}

template<typename PS,typename T>
    requires ( !NA::is_named_argument_v<std::decay_t<PS>> )
BlockBilinearForm<PS>
blockform2( PS&& ps, condensed_matrix_ptr_t<T> & m )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps), m );
}

template<typename ... Ts>
    requires Feel::detail::Blockform2NamedArgumentsWithTrial<Ts...>
auto
blockform2( Ts&& ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& test = args.get( _test );
    auto&& trial = args.get( _trial );
    return Feel::detail::blockform2Named( args, test, trial );
}

template<typename ... Ts>
    requires ( Feel::detail::Blockform2NamedArguments<Ts...> &&
               !Feel::detail::Blockform2NamedArgumentsWithTrial<Ts...> )
auto
blockform2( Ts&& ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& test = args.get( _test );
    return Feel::detail::blockform2Named( args, test, test );
}
/**
 * Handles linear form over a product of spaces
 */
template<typename PS>
class BlockLinearForm
{
public :
    using value_type = typename decay_type<PS>::value_type;
    using product_space_t = decay_type<PS>;
    using condensed_vector_type = VectorCondensed<value_type>;
    using condensed_vector_ptrtype = std::shared_ptr<condensed_vector_type>;
    using vector_ptrtype = condensed_vector_ptrtype;
    
    BlockLinearForm() = default;
    BlockLinearForm( BlockLinearForm const& ) = default;

    template<typename T, typename BackendT>
        requires StaticProductSpacesType<T>
    BlockLinearForm( T&& ps, solve::strategy s, BackendT&& b )
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(s, blockVector(M_ps), std::forward<BackendT>(b), false))
        {}
    template<typename T>
        requires StaticProductSpacesType<T>
    BlockLinearForm(T&& ps)
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), backend(), false))
        {}
    template<typename T>
        requires DynamicProductSpaceType<T>
    BlockLinearForm(T&& ps)
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), backend(), false))
        {}    
    template<typename T, typename BackendT>
        requires StaticProductSpacesType<T>
    BlockLinearForm(T&& ps, BackendT&& b)
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), std::forward<BackendT>(b), false))
        {}
    template<typename T>
    BlockLinearForm(T&& ps, condensed_vector_ptrtype v )
        :
        M_ps(std::forward<T>(ps)),
        M_vector(v)
        {}
    BlockLinearForm( BlockLinearForm&& ) = default;
    BlockLinearForm& operator=( BlockLinearForm && lf ) = default;
    BlockLinearForm& operator=( BlockLinearForm const& lf )
        {
            if ( this == &lf )
                return *this;

            bool same_spaces = ( M_ps == lf.M_ps );
            M_ps = lf.M_ps;
            if ( !M_vector || !same_spaces )
            {
                M_vector = std::make_shared<condensed_vector_type>( lf.M_vector->solveStrategy(), blockVector(M_ps), lf.M_vector->backend(), false );
                M_vector->setBackend( lf.M_vector->backend()->clone() );
            }
            M_vector->zero();
            M_vector->add( 1., *lf.M_vector->getVector() );
            
            return *this;
        }

#if 0
    template<typename N1>
    decltype(auto) operator()( N1 n1 )
        {
            cout << "filling out vector block (" << n1 << ")\n";
            return form1(_test=M_ps[n1],_vector=M_vector, _rowstart=int(n1) );
        }
#endif
    template<typename N1>
    decltype(auto) operator()( N1 n1, int s = 0 )
        {
            int n = 0;
            auto&& spaces=hana::if_( hana::bool_<Feel::is_shared_ptr_v<PS>>{},
                                     []( auto&& x ) { return *x; },
                                     []( auto&& x ) { return x; } )(M_ps);
            auto space = hana::at( spaces.tupleSpaces(), n1 );

            return hana::eval_if(std::is_base_of<ProductSpaceBase,decay_type<decltype(space)>>{},
                                 [&] (auto _) {
                                     VLOG(2) << "filling out dyn vector block (" << int(n1) + s  << ") condense=" << M_vector->staticCondensation() << "\n";
                                     auto vec = M_vector->monolithic() ? M_vector->vectorPtr() : M_vector->block( int(n1)+s );
                                     return form1(_test=(*_(space))[s],_vector=vec, _rowstart=int(n1)+s );
                                 },
                                 [&] (auto _){
                                     VLOG(2) << "filling out vector block (" << n1  << ") condense=" << M_vector->staticCondensation() << "\n";
                                     auto vec = M_vector->monolithic() ? M_vector->vectorPtr() : M_vector->block( int(n1) );
                                     return form1(_test=_(space),_vector=vec, _rowstart=int(n1) );
                                 });
        }

    decltype(auto) operator()( int n1 )
        {
            VLOG(2) << "filling out vector block (" << n1 << ") condense=" << M_vector->staticCondensation() << "\n";
            auto vec = M_vector->monolithic() ? M_vector->vectorPtr() : M_vector->block( int(n1) );
            return form1(_test=M_ps[n1],_vector=vec, _rowstart=int(n1) );
        }
    template<typename T>
    void setFunctionSpace( T&& ps )
        {
            M_ps = std::forward<T>(ps);
        }
    template<typename BackendT>
    void setStrategy( BackendT&& b )
        {
            M_vector = std::make_shared<condensed_vector_type>(blockVector(M_ps), std::forward<BackendT>(b), false);
            
            
        }
    template<typename BackendT>
    void setStrategy( solve::strategy s, BackendT&& b )
        {
            M_vector = std::make_shared<condensed_vector_type>(s, blockVector(M_ps), std::forward<BackendT>(b), false);
        }
    void close()
        {
            M_vector->close();
        }
    product_space_t functionSpace() const { return M_ps; }

    condensed_vector_ptrtype const& vectorPtr() const { return M_vector; }
    condensed_vector_ptrtype vectorPtr() { return M_vector; }
    condensed_vector_type const& vector() const { return *M_vector; }
    condensed_vector_type& vector() { return *M_vector; }
    auto sum() const { return M_vector->sum(); }
    auto min() const { return M_vector->min(); }
    auto max() const { return M_vector->max(); }
    auto l1Norm() const { return M_vector->l1Norm(); }
    auto l2Norm() const { return M_vector->l2Norm(); }
    auto linftyNorm() const { return M_vector->linftyNorm(); }

    /**
     * set linear form to 0
     */
    void zero() { M_vector->zero(); }

    //!
    //! zero block @param n1
    //!
    void zeroBlock( int n1 ) { M_vector->zeroBlock( n1 ); }

    BlockLinearForm& operator+=( BlockLinearForm const& l )
        {
            if ( this == &l )
            {
                M_vector->scale( 2. );
                return *this;
            }

            *M_vector += *l.M_vector;

            return *this;
        }
    product_space_t M_ps;
    condensed_vector_ptrtype M_vector;
};

template<typename PS>
using blockform1_t = BlockLinearForm<PS>;
template<typename TestPS, typename TrialPS = TestPS>
using blockform2_t = BlockBilinearForm<TestPS, TrialPS>;


template<typename PS>
BlockLinearForm<PS>
blockform1( PS&& ps )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps) );
}

template<typename PS, typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, BackendT&& b )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), std::forward<BackendT>(b) );
}

template<typename PS, typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, solve::strategy s, BackendT&& b )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), s, std::forward<BackendT>(b) );
}

template<typename PS,typename T>
BlockLinearForm<PS>
blockform1( PS&& ps, condensed_vector_ptr_t<T> v )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), v );
}


}
#endif
