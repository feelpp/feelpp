/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_TRAITS_HPP
#define FEELPP_TOOLBOXES_CORE_TRAITS_HPP 1

namespace Feel
{
namespace FeelModels
{

template<typename RangeType>
struct RangeTraits
{
    using element_iterator = typename RangeType::iterator_t;
    using range_elt_type = typename RangeType::element_t;
    using const_t = std::remove_reference_t<range_elt_type>;
    using element_type = std::remove_const_t<const_t>;
    //using the_face_element_type = std::remove_const_t<const_t>;
    //using element_type = typename the_face_element_type::super2::template Element<the_face_element_type>::type;
};
template<typename RangeType, typename ExprType>
struct ExprTraits
{
    typedef typename RangeTraits<RangeType>::element_type element_type;

    typedef typename element_type::gm_type gm_type;
    typedef typename gm_type::template Context</*ExprType::context|vm::JACOBIAN,*/ element_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename ExprType::template tensor<map_gmc_type> eval_expr_type;
    typedef typename eval_expr_type::shape shape;
};

template<typename MeshElementType, typename ExprType>
struct ExprTraitsFromMeshElement
{
    typedef MeshElementType element_type;

    typedef typename element_type::gm_type gm_type;
    typedef typename gm_type::template Context</*ExprType::context|vm::JACOBIAN,*/ element_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename ExprType::template tensor<map_gmc_type> eval_expr_type;
    typedef typename eval_expr_type::shape shape;
};

template<typename ContextType, typename ExprType>
struct ExprTraitsFromContext
{
    using gmc_type = typename ContextType::mapped_type::first_type::element_type::gmc_type;
    using gmc_ptrtype = std::shared_ptr<gmc_type>;
    using map_gmc_type = fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> >;
    using eval_expr_type = typename ExprType::template tensor<map_gmc_type>;
    using shape = typename eval_expr_type::shape;
};

} // namespace FeelModels
} // namespace Feel

#endif
