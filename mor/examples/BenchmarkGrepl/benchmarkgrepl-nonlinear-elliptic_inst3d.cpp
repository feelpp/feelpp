#include "benchmarkgrepl-nonlinear-elliptic_impl.cpp"

namespace Feel
{

FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_3dP1, BenchmarkGreplNonlinearElliptic<1 BOOST_PP_COMMA() 3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,3dP1) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_3dP2, BenchmarkGreplNonlinearElliptic<2 BOOST_PP_COMMA() 3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,3dP2) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplNonlinearElliptic_3dP3, BenchmarkGreplNonlinearElliptic<3 BOOST_PP_COMMA() 3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,3dP3) )

template class BenchmarkGreplNonlinearElliptic<1,3>;
template class BenchmarkGreplNonlinearElliptic<2,3>;
template class BenchmarkGreplNonlinearElliptic<3,3>;

}
