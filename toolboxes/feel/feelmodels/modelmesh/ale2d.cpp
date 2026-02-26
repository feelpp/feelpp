
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

} // namespace FeelModels
} // namespace Feel
