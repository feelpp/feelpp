
#include <feel/feelmodels/modelmesh/winslow.cpp>

namespace Feel
{
namespace FeelModels
{

template class Winslow< Mesh<Simplex<2,1> >, 1 >;
template class Winslow< Mesh<Hypercube<2,1> >, 1 >;

} // namespace FeelModels
} // namespace Feel
