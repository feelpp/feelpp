#include "benchmarkgrepl-nonlinear-elliptic_impl.cpp"

namespace Feel
{

template class BenchmarkGreplNonlinearElliptic<1,2>;
template class BenchmarkGreplNonlinearElliptic<2,2>;
template class BenchmarkGreplNonlinearElliptic<3,2>;

FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_2dP1, BenchmarkGreplNonlinearElliptic<1 BOOST_PP_COMMA() 2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,2dP1) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_2dP2, BenchmarkGreplNonlinearElliptic<2 BOOST_PP_COMMA() 2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,2dP2) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_2dP3, BenchmarkGreplNonlinearElliptic<3 BOOST_PP_COMMA() 2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,2dP3) )

}
