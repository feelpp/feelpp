
#define  FEELP_INSTANTIATE
#include <feel/feeldiscr/pch.hpp>

namespace Feel {
template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<0,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<1,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<2,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<3,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<0,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<1,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<2,Scalar>>>;
template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<3,Scalar>>>;
}