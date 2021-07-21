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
    typedef typename boost::tuples::template element<1, RangeType>::type element_iterator;
    typedef typename boost::unwrap_reference<typename element_iterator::value_type>::type range_elt_type;
    typedef typename boost::remove_reference<range_elt_type>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type element_type;
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
    using gmc_type = typename ContextType::mapped_type::element_type::gmc_type;
    using gmc_ptrtype = std::shared_ptr<gmc_type>;
    using map_gmc_type = fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> >;
    using eval_expr_type = typename ExprType::template tensor<map_gmc_type>;
    using shape = typename eval_expr_type::shape;
};

} // namespace FeelModels
} // namespace Feel

#endif
