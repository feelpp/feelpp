#include "benchmarkgrepl-nonlinear-elliptic_impl.cpp"

namespace Feel
{

FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic1_2D, BenchmarkGreplNonlinearElliptic<1 BOOST_PP_COMMA() 2>, benchmarkgreplnonlinearelliptic1_2d )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic2_2D, BenchmarkGreplNonlinearElliptic<2 BOOST_PP_COMMA() 2>, benchmarkgreplnonlinearelliptic2_2d )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic3_2D, BenchmarkGreplNonlinearElliptic<3 BOOST_PP_COMMA() 2>, benchmarkgreplnonlinearelliptic3_2d )

template class BenchmarkGreplNonlinearElliptic<1,2>;
template class BenchmarkGreplNonlinearElliptic<2,2>;
template class BenchmarkGreplNonlinearElliptic<3,2>;

}
