/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/** @file tabulation.hpp
 * @brief Requested reference tables, reusable storage and modal-to-nodal adapters.
 */
#ifndef FEELPP_FEELPOLY_TABULATION_HPP
#define FEELPP_FEELPOLY_TABULATION_HPP 1
#include <feel/feelpoly/polynomialcoefficients.hpp>
#include <type_traits>

namespace Feel
{

/** @brief Reference derivative tables requested independently of polynomial degree. */
enum class Evaluation : unsigned
{
    Values = 1,    /**< Function values. */
    Gradients = 2, /**< All reference first partials. */
    Hessians = 4,  /**< All ordered reference second partials. */
    All = 7        /**< Values, gradients and Hessians. */
};

/** @brief Combine reference-table requests. */
[[nodiscard]] constexpr Evaluation operator|( Evaluation a, Evaluation b ) noexcept
{
    return static_cast<Evaluation>( static_cast<unsigned>( a ) | static_cast<unsigned>( b ) );
}

/** @brief Whether a request contains the specified table group. */
[[nodiscard]] constexpr bool includes( Evaluation request, Evaluation group ) noexcept
{
    return ( static_cast<unsigned>( request ) & static_cast<unsigned>( group ) ) == static_cast<unsigned>( group );
}

namespace detail
{

/** @brief Reject empty requests and unknown bits before dispatch or allocation. */
[[nodiscard]] constexpr bool validEvaluation( Evaluation request ) noexcept
{
    return static_cast<unsigned>( request ) > 0 && static_cast<unsigned>( request ) <= 7;
}

/** @brief Highest derivative algebra needed by a requested table collection. */
[[nodiscard]] constexpr unsigned derivativeOrder( Evaluation request ) noexcept
{
    if ( includes( request, Evaluation::Hessians ) )
    {
        return 2;
    }

    if ( includes( request, Evaluation::Gradients ) )
    {
        return 1;
    }

    return 0;
}

/** @brief Number of scalar matrix tables in the requested collection. */
template <int Dim>
[[nodiscard]] constexpr std::size_t derivativeTableCount( Evaluation request ) noexcept
{
    return includes( request, Evaluation::Values ) + Dim * includes( request, Evaluation::Gradients ) + Dim * Dim * includes( request, Evaluation::Hessians );
}

/** @brief Visit value, coordinate-ordered gradient and row-ordered Hessian slots.
 * @param operation Callable templated on derivative order, consuming (slot,i,j).
 */
template <int Dim, Evaluation Request, typename Operation>
constexpr void forEachDerivative( Operation&& operation )
{
    static_assert( validEvaluation( Request ) );
    std::size_t slot = 0;

    // Keep values first, then coordinate gradients, then ordered Hessian pairs.
    if constexpr ( includes( Request, Evaluation::Values ) )
    {
        operation.template operator()<0>( slot++, 0, 0 );
    }

    if constexpr ( includes( Request, Evaluation::Gradients ) )
    {
        for ( int i = 0; i < Dim; ++i )
        {
            operation.template operator()<1>( slot++, i, 0 );
        }
    }

    if constexpr ( includes( Request, Evaluation::Hessians ) )
    {
        for ( int i = 0; i < Dim; ++i )
        {
            for ( int j = 0; j < Dim; ++j )
            {
                operation.template operator()<2>( slot++, i, j );
            }
        }
    }
}

/** @brief Exponent of one coordinate in a requested derivative multi-index.
 * @pre The slot and axis are valid indices in the requested tables and reference coordinates.
 */
template <int Dim, Evaluation Request>
[[nodiscard]] constexpr unsigned derivativeExponent( std::size_t slot, int axis ) noexcept
{
    static_assert( validEvaluation( Request ) );

    // Remove a leading value slot before indexing derivative groups.
    if constexpr ( includes( Request, Evaluation::Values ) )
    {
        if ( slot == 0 ) return 0;
        --slot;
    }

    // A first partial contributes one derivative along its selected coordinate.
    if constexpr ( includes( Request, Evaluation::Gradients ) )
    {
        if constexpr ( includes( Request, Evaluation::Hessians ) )
        {
            if ( slot < Dim ) return slot == static_cast<std::size_t>( axis );
            slot -= Dim;
        }
        else
            return slot == static_cast<std::size_t>( axis );
    }

    // Ordered second partials contribute once for each matching direction.
    if constexpr ( includes( Request, Evaluation::Hessians ) )
        return ( slot / Dim == static_cast<std::size_t>( axis ) ) +
               ( slot % Dim == static_cast<std::size_t>( axis ) );
    else
        return 0;
}

/** @brief Expansion degree, which can differ from an FE family's index. */
template <class Source>
[[nodiscard]] uint16_type polynomialOrder( Source const& source )
{
    // FE family indices may differ from the polynomial expansion degree.
    if constexpr ( requires { source.polynomialDegree(); } )
        return source.polynomialDegree();
    else if constexpr ( requires { source.order(); } )
        return source.order();
    else
        return staticOrderValue<Source>();
}

/** @brief Native modal kernels retain their own degree and normalization policy. */
template <typename Source, typename = void>
struct TabulationKernel
{
    using type = Source;
};

/** @brief FE/polynomial facades select the effective degree independently of their seed. */
template <typename Source>
struct TabulationKernel<Source, std::void_t<typename Source::basis_type>>
{
    using type = typename Source::basis_type::template TabulationKernelFor<orderSpecification<Source>()>;
};

/** @brief Selected native kernel for a modal or coefficient-transformed source. */
template <typename Source>
using TabulationKernelFor = typename TabulationKernel<Source>::type;

/** Whether functions are represented by rows of a modal coefficient matrix. */
template <typename Source>
inline constexpr bool transformedTabulation = requires( Source const& source ) { source.coefficients(); };

/** Stored field components are independent of degree and reference directions. */
template <typename Source>
inline constexpr int tabulationComponents = []
{
    if constexpr ( transformedTabulation<Source> )
    {
        return Source::nComponents;
    }
    else
    {
        return 1;
    }
}();

/** Function rows are independent of degree for a transformed polynomial set. */
template <typename Source>
inline constexpr int tabulationFunctionExtent = transformedTabulation<Source> ? Dynamic : TabulationKernelFor<Source>::basisExtent;
} // namespace detail

/** @brief Capability of a primitive kernel or facade backed by native reference tabulation. */
template <typename Source>
concept NativeTabulationSource = requires { Source::dimension; Source::basisExtent; typename Source::Values; typename Source::template Workspace<0>; } ||
                                 ( requires { Source::nDim; } && Source::nDim > 0 && requires { typename Source::basis_type::template TabulationKernelFor<orderSpecification<Source>()>; } );

/** Borrowed selection of tables whose owner outlives the view. */
template <typename Owner, Evaluation Request>
class ReferenceTabulationView;

/** @brief Own requested reference tables with independent function and sample axes.
 * @tparam T Scalar type.
 * @tparam Dim Reference dimension.
 * @tparam OrderSpec Polynomial degree specification, independent of derivative order.
 * @tparam Functions Number of interleaved function/component rows, or Dynamic.
 * @tparam Samples Number of sample columns, or Dynamic.
 * @tparam Request Compile-time selection of reference tables.
 *
 * Each row is function*nComponents+component; every matrix has sample columns.
 * Gradients add one reference direction; Hessians use ordered (i,j), i*Dim+j.
 * Views and references expire on owner mutation, movement or destruction.
 * There are no retained views of source coefficients or input coordinates.
 */
template <typename T, int Dim, int OrderSpec, int Functions, int Samples, Evaluation Request,
          int Components = 1, int Layout = Eigen::RowMajor>
class ReferenceTabulation : public OrderBase<OrderSpec>
{
    static_assert( Dim >= 0 && Dim <= 3 && Components > 0 && ::Feel::detail::validEvaluation( Request ) );

  public:
    /** Scalar type shared by all owned tables. */
    using value_type = T;

    /** Requested table groups. */
    static constexpr Evaluation request = Request;

    /** Reference dimension, independent of component count. */
    static constexpr int dimension = Dim;

    /** Unmodified polynomial-order specification. */
    static constexpr int order_spec = OrderSpec;

    /** Interleaved row extent, including all stored field components. */
    static constexpr int rowsAtCompileTime = Functions;

    /** Stored field components per function, independent of derivative directions. */
    static constexpr int components = Components;

    /** Independent sample-column extent. */
    static constexpr int samplesAtCompileTime = Samples;

    /** Number of requested matrix owners. */
    static constexpr std::size_t tableCount = ::Feel::detail::derivativeTableCount<Dim>( Request );
    static_assert( tableCount > 0, "a zero-dimensional reference requires a value table" );

    /** Checked owner; the complete table collection shares the inline byte budget. */
    using Buffer = ::Feel::detail::EigenMatrixBuffer<T, Functions, Samples, ::Feel::detail::defaultInlineMatrixBytes / tableCount, Layout>;

    /** @brief Allocate initialized tables at a checked polynomial degree and shape. */
    explicit ReferenceTabulation( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec },
                                  Eigen::Index functions = Functions == Dynamic ? 0 : Functions,
                                  Eigen::Index samples = Samples == Dynamic ? 0 : Samples )
        : OrderBase<OrderSpec>( degree ),
          M_tables( ::Feel::detail::makeMatrixBuffers<Buffer, tableCount>( checkedRows( functions ), samples ) ) {}

    /** @brief Actual number of interleaved function/component rows. */
    [[nodiscard]] Eigen::Index rows() const noexcept { return M_tables[0].rows(); }

    /** @brief Actual number of sample columns. */
    [[nodiscard]] Eigen::Index cols() const noexcept { return M_tables[0].cols(); }

    /** @brief Actual number of functions, excluding the interleaved component axis. */
    [[nodiscard]] Eigen::Index functionCount() const noexcept { return rows() / Components; }

    /** @brief Borrow one checked table in the documented derivative-slot order. */
    [[nodiscard]] Buffer& table( std::size_t slot ) & { return M_tables.at( slot ); }

    /** @brief Borrow immutable table storage while this owner remains unchanged. */
    [[nodiscard]] Buffer const& table( std::size_t slot ) const& { return M_tables.at( slot ); }

    /** Prevent references into temporary owners. */
    Buffer&& table( std::size_t ) && = delete;

    /** Prevent references into const temporary owners. */
    Buffer const&& table( std::size_t ) const&& = delete;

    /** @brief Borrow the value table; throw if values were not requested. */
    [[nodiscard]] Buffer& values() & { return table( valueSlot() ); }

    /** @brief Borrow immutable values. */
    [[nodiscard]] Buffer const& values() const& { return table( valueSlot() ); }

    /** Prevent a value view into a temporary owner. */
    Buffer&& values() && = delete;

    /** Prevent a const value view into a temporary owner. */
    Buffer const&& values() const&& = delete;

    /** @brief Borrow one checked reference gradient direction. */
    [[nodiscard]] Buffer& gradient( int direction ) & { return table( gradientSlot( direction ) ); }

    /** @brief Borrow one immutable reference gradient direction. */
    [[nodiscard]] Buffer const& gradient( int direction ) const& { return table( gradientSlot( direction ) ); }

    /** Prevent a gradient view into a temporary owner. */
    Buffer&& gradient( int ) && = delete;

    /** Prevent a const gradient view into a temporary owner. */
    Buffer const&& gradient( int ) const&& = delete;

    /** @brief Borrow one checked ordered reference Hessian direction pair. */
    [[nodiscard]] Buffer& hessian( int i, int j ) & { return table( hessianSlot( i, j ) ); }

    /** @brief Borrow one immutable ordered reference Hessian direction pair. */
    [[nodiscard]] Buffer const& hessian( int i, int j ) const& { return table( hessianSlot( i, j ) ); }

    /** Prevent a Hessian view into a temporary owner. */
    Buffer&& hessian( int, int ) && = delete;

    /** Prevent a const Hessian view into a temporary owner. */
    Buffer const&& hessian( int, int ) const&& = delete;

    /** @brief Borrow output pointers during one native kernel call. */
    [[nodiscard]] auto pointers() & { return ::Feel::detail::matrixBufferPointers( M_tables ); }

    /** Prevent pointers into a temporary owner. */
    auto pointers() && = delete;

    /** @brief Select existing owners without constructing or copying coefficient storage. */
    template <Evaluation Selected>
    [[nodiscard]] auto select() & { return ReferenceTabulationView<ReferenceTabulation, Selected>{ *this }; }

    /** Prevent a selection into a temporary owner. */
    template <Evaluation Selected>
    auto select() && = delete;

  private:
    /** @brief Reject incomplete function/component blocks before allocation. */
    [[nodiscard]] static Eigen::Index checkedRows( Eigen::Index rows )
    {
        if ( rows < 0 || rows % Components != 0 ) throw std::invalid_argument( "tabulation rows must contain complete component blocks" );
        return rows;
    }

    /** @brief Check that value storage exists. */
    [[nodiscard]] static std::size_t valueSlot()
    {
        if constexpr ( !includes( Request, Evaluation::Values ) ) throw std::logic_error( "values were not requested" );
        return 0;
    }

    /** @brief Check a direction independently of the chosen table group. */
    static void validateDirection( int direction )
    {
        if ( direction < 0 || direction >= Dim ) throw std::out_of_range( "invalid tabulation direction" );
    }

    /** @brief Check and locate a first partial. */
    [[nodiscard]] static std::size_t gradientSlot( int direction )
    {
        validateDirection( direction );
        if constexpr ( !includes( Request, Evaluation::Gradients ) ) throw std::logic_error( "gradients were not requested" );
        return includes( Request, Evaluation::Values ) + direction;
    }

    /** @brief Check and locate an ordered second partial. */
    [[nodiscard]] static std::size_t hessianSlot( int i, int j )
    {
        validateDirection( i );
        validateDirection( j );
        if constexpr ( !includes( Request, Evaluation::Hessians ) ) throw std::logic_error( "Hessians were not requested" );
        return includes( Request, Evaluation::Values ) + Dim * includes( Request, Evaluation::Gradients ) + i * Dim + j;
    }

    /** Sole owners of the requested reference tables. */
    std::array<Buffer, tableCount> M_tables;
};

/** @brief Borrow a compile-time subset of a longer-lived owning table collection.
 * @note The view owns only pointers. Owner movement, destruction or table resizing
 * invalidates borrowed coefficient views; the owner must outlive this selection.
 */
template <typename Owner, Evaluation Request>
class ReferenceTabulationView
{
    static_assert( includes( Owner::request, Request ) );

  public:
    /** Borrowed matrix owner type, retaining all declared extents and layout. */
    using Buffer = typename Owner::Buffer;

    /** Scalar type retained from the owner. */
    using value_type = typename Owner::value_type;

    /** Interleaved component count retained from the owner. */
    static constexpr int components = Owner::components;

    /** Requested groups selected from the owner. */
    static constexpr Evaluation request = Request;

    /** Reference dimension. */
    static constexpr int dimension = Owner::dimension;

    /** Number of selected matrix pointers. */
    static constexpr std::size_t tableCount = ::Feel::detail::derivativeTableCount<dimension>( Request );

    /** @brief Borrow the selected owners without allocating storage. */
    explicit ReferenceTabulationView( Owner& owner )
        : M_owner( owner )
    {
        // Borrow the selected canonical buffers without copying their coefficients.
        ::Feel::detail::forEachDerivative<dimension, Request>( [&]<unsigned Derivative>( std::size_t slot, int i, int j )
                                                               {
            if constexpr ( Derivative == 0 ) M_tables[slot] = &owner.values();
            if constexpr ( Derivative == 1 ) M_tables[slot] = &owner.gradient( i );
            if constexpr ( Derivative == 2 ) M_tables[slot] = &owner.hessian( i, j ); } );
    }

    /** @brief Polynomial degree retained by the owner. */
    [[nodiscard]] auto order() const noexcept { return M_owner.order(); }

    /** @brief Borrow one selected matrix owner. */
    [[nodiscard]] Buffer& table( std::size_t slot ) const { return *M_tables.at( slot ); }

    /** @brief Borrow output pointers during one native kernel call. */
    [[nodiscard]] auto pointers() const noexcept { return M_tables; }

  private:
    /** Owner whose lifetime and degree constrain this view. */
    Owner& M_owner;

    /** Pointers to the selected owners; no coefficient storage is duplicated. */
    std::array<Buffer*, tableCount> M_tables;
};

/** @brief Reusable modal recurrence and coefficient-transform scratch.
 * @tparam Source Modal kernel, polynomial or FE basis.
 * @tparam Request Compile-time table selection.
 * @tparam Samples Independent sample extent, or Dynamic.
 * @note Scratch is private and never aliases caller outputs or coefficients.
 */
template <typename Source, Evaluation Request, int Samples = Dynamic>
    requires NativeTabulationSource<Source>
class TabulationWorkspace : public OrderBase<orderSpecification<Source>()>
{

    /** Native numerical kernel selected independently of the facade seed degree. */
    using Kernel = ::Feel::detail::TabulationKernelFor<Source>;

    /** Scalar used by the native modal buffers. */
    using T = typename Kernel::Values::value_type;

    /** Modal output owned only for a coefficient-transformed source. */
    using ModalTables = ReferenceTabulation<T, Kernel::dimension, orderSpecification<Source>(), Kernel::basisExtent, Samples, Request>;

  public:
    /** Source policy whose degree and coefficient representation bind this workspace. */
    using source_type = Source;

    /** Requested table groups, retained into the numerical kernel. */
    static constexpr Evaluation request = Request;

    /** @brief Bind recurrence storage to one effective degree. */
    explicit TabulationWorkspace( RuntimeOrder degree = RuntimeOrder{ orderSpecification<Source>() == Dynamic ? 0 : orderSpecification<Source>() } )
        : OrderBase<orderSpecification<Source>()>( degree ),
          M_workspace( degree ),
          M_modal( makeModalStorage( degree ) ) {}

  private:
    /** The common adapter is the only consumer of private modal scratch. */
    template <typename S, typename P, typename O, typename W, Evaluation R, int N>
    friend void tabulateInto( S const&, P const&, O&, TabulationWorkspace<W, R, N>& );

    /** @brief Borrow native recurrence scratch during a tabulation call. */
    [[nodiscard]] auto& recurrence() & noexcept { return M_workspace; }

    /** @brief Borrow private modal output while coefficient multiplication consumes it. */
    [[nodiscard]] auto& modal() & noexcept
        requires ::Feel::detail::transformedTabulation<Source>
    {
        return M_modal[0];
    }

  private:
    /** Highest requested derivative order selects the compact numerical algebra. */
    typename Kernel::template Workspace<::Feel::detail::derivativeOrder( Request ), Samples> M_workspace;

    /** Modal staging exists only when C multiplies the modal tables. */
    using ModalStorage = std::array<ModalTables, ::Feel::detail::transformedTabulation<Source> ? 1 : 0>;

    /** @brief Construct modal owners directly, avoiding empty-owner replacement. */
    [[nodiscard]] static ModalStorage makeModalStorage( RuntimeOrder degree )
    {
        // Primitive kernels write directly into caller output and need no modal owners.
        if constexpr ( ::Feel::detail::transformedTabulation<Source> )
            return { ModalTables{ degree, Kernel{ degree }.basisCount(), Samples == Dynamic ? 0 : Samples } };
        else
            return {};
    }

    /** Private modal owners are constructed once at the selected degree. */
    ModalStorage M_modal;
};

/** @brief Allocate initialized tabulation storage without evaluating the source.
 * @tparam Request Requested reference derivative groups.
 * @tparam Samples Independent sample extent.
 * @tparam Functions Optional independently known function/component row extent.
 */
template <Evaluation Request, int Samples = Dynamic, int Functions = Dynamic, typename Source>
    requires NativeTabulationSource<Source>
[[nodiscard]] auto makeTabulation( Source const& source, Eigen::Index samples )
{
    // Keep degree, component rows and sample columns as independent type policies.
    using Kernel = ::Feel::detail::TabulationKernelFor<Source>;
    constexpr int rows = Functions != Dynamic ? Functions : ::Feel::detail::tabulationFunctionExtent<Source>;
    constexpr int components = ::Feel::detail::tabulationComponents<Source>;
    using Tables = ReferenceTabulation<typename Kernel::Values::value_type, Kernel::dimension, orderSpecification<Source>(), rows, Samples, Request, components>;

    // Transformed functions need not have the same cardinality as their modal basis.
    RuntimeOrder const degree{ ::Feel::detail::polynomialOrder( source ) };
    Eigen::Index count;
    if constexpr ( ::Feel::detail::transformedTabulation<Source> )
        count = source.coefficients().rows();
    else
        count = Kernel{ degree }.basisCount();

    return Tables{ degree, count, samples };
}

/** @brief Evaluate requested reference tables into reusable caller-owned output.
 * @param source Modal kernel or functions represented by C(function/component,mode).
 * @param points Column-oriented reference coordinates, including expressions or legacy matrices.
 * @param output Owning tables with the same degree and request as the workspace.
 * @param workspace Private modal scratch reused across point batches.
 * @note Coordinates may alias output. Source coefficients must remain unchanged
 * and their storage must be disjoint from destination owners during this call.
 */
template <typename Source, typename Points, typename Tables, typename WorkspaceSource, Evaluation Request, int Samples>
void tabulateInto( Source const& source, Points const& points, Tables& output,
                   TabulationWorkspace<WorkspaceSource, Request, Samples>& workspace )
{
    // Reject incompatible storage and kernel policies at the common API boundary.
    using Kernel = ::Feel::detail::TabulationKernelFor<Source>;
    using T = typename Kernel::Values::value_type;
    static_assert( std::is_same_v<Kernel, ::Feel::detail::TabulationKernelFor<WorkspaceSource>> );
    static_assert( std::is_same_v<typename Tables::value_type, T> );
    static_assert( Tables::components == ::Feel::detail::tabulationComponents<Source> );
    static_assert( ::Feel::detail::transformedTabulation<Source> == ::Feel::detail::transformedTabulation<WorkspaceSource> );
    static_assert( includes( Request, Tables::request ) && Tables::dimension == Kernel::dimension );
    static_assert( ::Feel::detail::derivativeOrder( Request ) == ::Feel::detail::derivativeOrder( Tables::request ) );

    // A retained workspace and its output must describe the source's effective degree.
    auto const degree = ::Feel::detail::polynomialOrder( source );
    if ( output.order() != degree || workspace.order() != degree )
        throw std::invalid_argument( "tabulation degree disagrees with source or workspace" );

    // Borrow compatible coordinates or materialize an input expression once.
    ::Feel::detail::withEigenMatrix<T>( points, [&]( auto const& coordinates )
                                        {
        Kernel kernel{ RuntimeOrder{ degree } };

        // Evaluate modal tables before resizing outputs that may own the coordinates.
        if constexpr ( ::Feel::detail::transformedTabulation<Source> )
        {
            auto coefficients = source.coefficients();
            if ( coefficients.cols() != kernel.basisCount() )
                throw std::invalid_argument( "tabulation coefficients disagree with expansion cardinality" );

            auto modal = workspace.modal().template select<Tables::request>();
            kernel.tabulateInto( modal, coordinates, workspace.recurrence() );

            // Apply the same N = C Phi transform to each requested derivative table.
            for ( std::size_t slot = 0; slot < Tables::tableCount; ++slot )
            {
                auto& target = output.table( slot );
                target.resizeForOverwrite( coefficients.rows(), coordinates.cols() );
                auto matrix = target.matrix();
                ::Feel::detail::multiplyPolynomialMatricesInto( matrix, coefficients, modal.table( slot ).matrix() );
            }
        }

        // Primitive modal sources already use the requested output representation.
        else kernel.tabulateInto( output, coordinates, workspace.recurrence() ); } );
}

/** @brief Allocate and evaluate one compile-time reference-table request.
 * @note The owning result retains polynomial order and independent declared point extents.
 */
template <Evaluation Request = Evaluation::All, typename Source, typename Points>
    requires NativeTabulationSource<Source>
[[nodiscard]] auto tabulate( Source const& source, Points const& points )
{
    using Kernel = ::Feel::detail::TabulationKernelFor<Source>;

    // Preserve known sample extents in both the owning result and temporary scratch.
    return ::Feel::detail::withEigenMatrix<typename Kernel::Values::value_type>( points, [&]( auto const& coordinates )
                                                                                 {
        constexpr int samples = std::remove_cvref_t<decltype( coordinates )>::ColsAtCompileTime;
        auto result = makeTabulation<Request, samples>( source, coordinates.cols() );
        TabulationWorkspace<Source, Request, samples> workspace{ RuntimeOrder{ ::Feel::detail::polynomialOrder( source ) } };

        tabulateInto( source, coordinates, result, workspace );
        return result; } );
}

/** @brief Dispatch a runtime request once to a compile-time request tag.
 * @param request Runtime selection of reference table groups.
 * @param operation Callable consuming an integral_constant<Evaluation,Request> tag.
 * @throws std::invalid_argument For an empty or unknown request mask.
 */
template <typename Operation>
void withEvaluationRequest( Evaluation request, Operation&& operation )
{
    // Resolve the runtime mask once before entering a specialized numerical kernel.
    auto consume = [&]<unsigned Mask>()
    { operation( std::integral_constant<Evaluation, static_cast<Evaluation>( Mask )>{} ); };

    switch ( static_cast<unsigned>( request ) )
    {
    case 1:
        consume.template operator()<1>();
        break;
    case 2:
        consume.template operator()<2>();
        break;
    case 3:
        consume.template operator()<3>();
        break;
    case 4:
        consume.template operator()<4>();
        break;
    case 5:
        consume.template operator()<5>();
        break;
    case 6:
        consume.template operator()<6>();
        break;
    case 7:
        consume.template operator()<7>();
        break;
    default:
        throw std::invalid_argument( "invalid reference tabulation request" );
    }
}

/** @brief Consume owning tables after one runtime request dispatch.
 * @note Borrowed table views must not escape the callback's owning result lifetime.
 */
template <typename Source, typename Points, typename Operation>
void withTabulation( Source const& source, Points const& points, Evaluation request, Operation&& operation )
{
    withEvaluationRequest( request, [&]( auto selection )
                           { auto tables = tabulate<decltype( selection )::value>( source, points ); operation( tables ); } );
}
} // namespace Feel
#endif // FEELPP_FEELPOLY_TABULATION_HPP
