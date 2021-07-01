
#include <feel/feelcrb/toolboxmor.hpp>
#include <feel/feelcrb/crbplugin.hpp>

namespace Feel
{
using base_type = bases<Lagrange<1, Scalar, Continuous, PointSetFekete> >;
using space_type = FunctionSpace<Mesh<Simplex<2> >, base_type>;
// FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor, ToolboxMor<FunctionSpace<Mesh<Simplex<2> >, base_type>>, toolboxmor )

template class ToolboxMor<space_type>;
}
