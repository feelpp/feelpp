
#include <feel/feelmodels/modelmesh/ale.cpp>
#include <feel/feelmodels/modelmesh/ale_impl.cpp>
#include <feel/feelmodels/modelmesh/meshale.cpp>
#include <feel/feelmodels/modelmesh/metricmeshadaptation.cpp>

namespace Feel
{
namespace FeelModels
{

template class ALE< Simplex<2,1>, 1 >;
template class ALE_IMPL::ALE< Simplex<2,1>, 1 >;
template class MeshALE< Simplex<2,1> >;
template class MetricMeshAdaptation<Simplex<2,1> >;
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template class ALE< Simplex<2,1>, 2 >;
template class ALE_IMPL::ALE< Simplex<2,1>, 2 >;
template class MeshALE< Simplex<2,2> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template class ALE< Simplex<2,1>, 3 >;
template class ALE_IMPL::ALE< Simplex<2,1>, 3 >;
template class MeshALE< Simplex<2,3> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template class ALE< Simplex<2,1>, 4 >;
template class ALE_IMPL::ALE< Simplex<2,1>, 4 >;
template class MeshALE< Simplex<2,4> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 5 )
template class ALE< Simplex<2,1>, 5 >;
template class ALE_IMPL::ALE< Simplex<2,1>, 5 >;
template class MeshALE< Simplex<2,5> >;
#endif // FEELPP_MESH_MAX_ORDER

} // namespace FeelModels
} // namespace Feel

