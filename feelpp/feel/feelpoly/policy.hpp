/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2011-2016 Feel++ Consortium
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file policy.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-12-03
 */
#ifndef FEELPP_FEELPOLY_POLICY_HPP
#define FEELPP_FEELPOLY_POLICY_HPP 1


#include <Eigen/Core>
#include <vector>
#include <feel/feelpoly/eigenextent.hpp>
#include <feel/feelpoly/componentlayout.hpp>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/traits.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/concepts.hpp>
namespace Feel
{
namespace ublas = boost::numeric::ublas;

namespace fem
{
/** Mapping classification used by finite-element geometry. */
enum TransformationType
{
    LINEAR,    ///< Affine geometry transformation.
    BILINEAR,  ///< Tensor-product bilinear geometry transformation.
    NONLINEAR  ///< General higher-order geometry transformation.
};
/** Compatibility spelling for the established geometry API. */
using transformation_type = TransformationType;
}


/**
 * Policy for \c Scalar polynomials or polynomial set of dimension
 * \p Dim
 * \note \c Scalar can be seen as rank 0 tensor polynomials
 */
template<uint16_type Dim>
struct Scalar : public ScalarBase, public ::Feel::detail::IdentityCoefficientLayout
{
    static constexpr uint16_type rank = 0;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = true;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = 1;
    static constexpr uint16_type nComponents1 = 1;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = 1;

};


/**
 * Policy for \c Vectorial polynomials or polynomial sets of dimension
 * \p Dim
 * \note \c Vectorial can be seen as rank 1 Tensor polynomials
 */
template<uint16_type Dim>
struct Vectorial : public VectorialBase, public ::Feel::detail::ComponentCoefficientLayout<Dim>
{
    static constexpr uint16_type rank = 1;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = true;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents1;

};

/**
 * Policy for \c Scalar polynomials or polynomial set of dimension
 * \p Dim
 */
namespace detail
{
template<uint16_type N, uint16_type M>
struct Field : public ::Feel::detail::IdentityCoefficientLayout
{
    static constexpr uint16_type rank = ( M > 1 );
    static constexpr uint16_type nDim = N;
    static constexpr uint16_type nVariables = N;

    static constexpr bool is_scalar = ( M==1 );
    static constexpr bool is_vectorial = ( N==M );
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = M;
    static constexpr uint16_type nComponents1 = M;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = 1;

};
}

template<uint16_type M>
struct Field
{
    template <uint16_type Nvar>
    struct Apply
    {
        typedef detail::Field<Nvar,M> type;
    };
    /** Compatibility spelling for existing template-policy clients. */
    template<uint16_type Nvar> using apply = Apply<Nvar>;

};

/**
 * Policy for rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor2 : public Tensor2Base, public ::Feel::detail::ComponentCoefficientLayout<Dim*Dim>
{
    static constexpr uint16_type rank = 2;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = true;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents2;

};

struct Tensor2SymmBase : Tensor2Base {};

/**
 * Policy for symmetric rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor2Symm : public Tensor2SymmBase, public ::Feel::detail::ComponentCoefficientLayout<Dim*Dim>
{
    static constexpr uint16_type rank = 2;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = true;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents2;

};

/**
 * Identify policies derived from Tensor2SymmBase through a standard Boolean trait.
 */
template<typename T>
struct IsSymmetric
    : bool_c<std::is_base_of_v<Tensor2SymmBase,T>>
{
};

/** Compatibility spelling for the established symmetry trait. */
template<typename T> using is_symm = IsSymmetric<T>;

/**
 * helper constant that is true if the T derives from Tensor2SymmBase, false
 * orherwise
 */
template<typename T>
constexpr bool is_symm_v = std::is_base_of_v<Tensor2SymmBase,T>;

template<typename T>
using is_scalar_t = bool_c<std::is_base_of_v<ScalarBase,T>>;
template<typename T>
constexpr bool is_scalar_v = std::is_base_of_v<ScalarBase,T>;
template<typename T>
using is_vectorial_t = bool_c<std::is_base_of_v<VectorialBase,T>>;
template<typename T>
constexpr bool is_vectorial_v = std::is_base_of_v<VectorialBase,T>;
template<typename T>
using is_tensor2_t = bool_c<std::is_base_of_v<Tensor2Base,T>>;
template<typename T>
constexpr bool is_tensor2_v = std::is_base_of_v<Tensor2Base,T>;

/**
 * Policy for rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor3
{
    static constexpr uint16_type rank = 3;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = true;

    static constexpr uint16_type nComponents = nDim*nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = nDim;
    static constexpr uint16_type nComponentsLast = nComponents3;
};

template<int Dim>
struct ListReturnTypes
{
    using return_types = type_list<Scalar<Dim>, Vectorial<Dim>, Tensor2<Dim>, Tensor3<Dim> >;
};
template<typename T1, typename T2>
struct ReturnSelect
{
    using type = if_t<std::is_same_v<T1, T2>,
                      T1,
                      if_t<(T1::rank > T2::rank), T1, T2>>;
};

enum EnumIndex
{
    GLOBAL_COMPONENT,
    COMPONENT_IN_COMPONENT,
    GLOBAL_FUNCTION_INDEX,
    PER_COMPONENT_FUNCTION_INDEX,
    COMPONENT_IN_COMPONENT_FUNCTION_INDEX,
    FUNCTION_INDEX
};

inline constexpr int_c<GLOBAL_COMPONENT> INDEX_GLOBAL_COMPONENT{};
inline constexpr int_c<COMPONENT_IN_COMPONENT> INDEX_COMPONENT_IN_COMPONENT{};
inline constexpr int_c<GLOBAL_FUNCTION_INDEX> INDEX_GLOBAL_FUNCTION_INDEX{};
inline constexpr int_c<PER_COMPONENT_FUNCTION_INDEX> INDEX_PER_COMPONENT_FUNCTION_INDEX{};
inline constexpr int_c<COMPONENT_IN_COMPONENT_FUNCTION_INDEX> INDEX_COMPONENT_IN_COMPONENT_FUNCTION_INDEX{};
inline constexpr int_c<FUNCTION_INDEX> INDEX_FUNCTION_INDEX{};


/**
 * Get the component type out the available types
 * \code
 * typedef typename Component<Vectorial<3> >::type component_type;
 * // component_type should be of type \c Scalar<3>
 * \endcode
 */
template<typename T>
struct GetComponent
{
    static constexpr uint16_type nDim = T::nDim;
    using type = if_t<std::is_same_v<T, Scalar<nDim> >,
                      Scalar<nDim>,
                      if_t<std::is_same_v<T, Vectorial<nDim> >,
                           Vectorial<nDim>,
                           if_t<std::is_same_v<T, Tensor2<nDim> >,
                                Tensor2<nDim>,
                                if_t<std::is_same_v<T, Tensor2Symm<nDim> >,
                                     Tensor2Symm<nDim>,
                                     mp::mp_void<>>>>>;
};

/**
 * Get the next rank
 * \code
 * typedef typename RankUp<Scalar<3> >::type rank_type;
 * // rank_type should be of type \c Vectorial<3>
 * \endcode
 */
template<typename T>
struct RankUp
{
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index + 1>;
};

template<typename T>
struct RankUp2
{
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index + 2>;
};

template<typename T>
struct RankDown2
{
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index - 2>;
};

template<typename T>
struct RankSame
{
    static constexpr uint16_type nDim = T::nDim;
    typedef T type;
};
template<typename T>
struct Rank0
{
    static constexpr uint16_type nDim = T::nDim;
    typedef Scalar<nDim> type;
};
template<typename T>
struct Rank1
{
    static constexpr uint16_type nDim = T::nDim;
    typedef Vectorial<nDim> type;
};

template<typename T>
struct RankDown
{
    static constexpr uint16_type nDim = T::nDim;
    using types = if_t<std::is_base_of_v<Tensor2SymmBase,T>,
                       type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2Symm<nDim>, Tensor3<nDim>>,
                       type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim>>>;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_eval_if_c<(index == 0),
                                  Scalar<nDim>,
                                  mp::mp_at,
                                  types, mp::mp_size_t<index - 1>>;
};
template<typename T>
using rankdown_t = typename RankDown<T>::type;

template<typename T>
struct RankCurl
{
    static constexpr uint16_type nDim = T::nDim;
    typedef typename if_t<( nDim == 3 ), RankSame<T>, RankDown<T> >::type type;
    static constexpr uint16_type value = ( nDim == 3 ) ? 3 : 1;
};

/**
 * Policy for orthogonal polynomial set which can be \p normalized or
 * not.  This can be used as a policy for orthogonal polynomial sets
 * in order to select either the normalized or un-normalized version
 */
template<bool normalized>
struct Normalized
{
    static constexpr bool is_normalized = normalized;
};

/**
 * Storage Policy using ublas for numerical type \p T
 *
 *
 */
template<typename T>
struct StorageUBlas
{
    typedef T value_type;
    typedef ublas::vector<value_type> vector_type;
    typedef ublas::matrix<value_type> matrix_type;
    typedef ublas::vector<matrix_type> vector_matrix_type;
    typedef ublas::vector<vector_matrix_type> vector_vector_matrix_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename matrix_node<value_type>::type points_type;
    typedef typename node<value_type>::type node_type;

};

/**
 * Storage Policy using Eigen for numerical type \p T
 */
template<typename T>
struct StorageEigen
{
    using value_type = T;
    using vector_type = Eigen::Matrix<value_type, Eigen::Dynamic, 1>;
    using matrix_type = typename ::Feel::detail::EigenMatrixBuffer<value_type, Dynamic, Dynamic,
        ::Feel::detail::defaultInlineMatrixBytes, Eigen::ColMajor>::matrix_type;
    using vector_matrix_type = std::vector<matrix_type>;
    using vector_vector_matrix_type = std::vector<vector_matrix_type>;
    using matrix_node_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
    using points_type = matrix_node_type;
    using node_type = Eigen::Matrix<value_type, Eigen::Dynamic, 1>;
    /** Order-aware scalar-basis storage using the shared extent/allocation bridge.
     * ExpansionOrderSpec is the modal expansion degree, including any family shift.
     * Spatial dimension and cell identity are supplied by the layer owning them.
     */
    template<int Dim, int ExpansionOrderSpec, ::Feel::detail::NodalCellKind CellKind>
    using ForBasis = ::Feel::detail::PolynomialStorage<value_type, Dim, ExpansionOrderSpec, CellKind>;
};

} // Feel

#endif /* FEELPP_FEELPOLY_POLICY_HPP */
