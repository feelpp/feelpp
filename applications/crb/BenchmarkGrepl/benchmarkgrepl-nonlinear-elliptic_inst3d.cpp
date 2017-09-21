#include "benchmarkgrepl-nonlinear-elliptic_impl.cpp"

namespace Feel
{

FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic1_3D, BenchmarkGreplNonlinearElliptic<1 BOOST_PP_COMMA() 3>, benchmarkgreplnonlinearelliptic1_3d )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic2_3D, BenchmarkGreplNonlinearElliptic<2 BOOST_PP_COMMA() 3>, benchmarkgreplnonlinearelliptic2_3d )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic3_3D, BenchmarkGreplNonlinearElliptic<3 BOOST_PP_COMMA() 3>, benchmarkgreplnonlinearelliptic3_3d )

template class BenchmarkGreplNonlinearElliptic<1,3>;
template class BenchmarkGreplNonlinearElliptic<2,3>;
template class BenchmarkGreplNonlinearElliptic<3,3>;

}
