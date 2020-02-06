
#include <toolboxmor_impl.cpp>
#include <feel/feelcrb/crbplugin.hpp>

namespace Feel
{
using mesh_type = Mesh<Simplex<2>>;
using base_type = bases<Lagrange<1, Scalar, Continuous, PointSetFekete> >;
using space_type = FunctionSpace<mesh_type, base_type>;
FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor, ToolboxMor<FunctionSpace<mesh_type, base_type>>, toolboxmor )

template class ToolboxMor<space_type>;
}
