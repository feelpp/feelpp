
#include <feel/feelmodels/fsi/fsi.hpp>

namespace Feel
{
namespace FeelModels
{

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDEStrongDirichletBC_Fluid( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        this->log("FSI","updateLinearPDEStrongDirichletBC_Fluid", "start" );

        auto mesh = M_fluidModel->mesh();
        auto Xh = M_fluidModel->spaceVelocityPressure();
        auto rowStartInMatrix = M_fluidModel->rowStartInMatrix();
        auto colStartInMatrix = M_fluidModel->colStartInMatrix();
        auto rowStartInVector = M_fluidModel->rowStartInVector();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=rowStartInMatrix,
                                                  _colstart=colStartInMatrix );
        auto const& u = M_fluidModel->fieldVelocity();
        bilinearForm +=
            on( _range=markedfaces(mesh, M_fluidModel->markersNameMovingBoundary()),
                _element=u, _rhs=F,
                _expr=idv(M_fluidModel->meshVelocity2()) );

        this->log("FSI","updateLinearPDEStrongDirichletBC_Fluid", "finish" );
}
}

} // namespace FeelModels
} // namespace Feel
