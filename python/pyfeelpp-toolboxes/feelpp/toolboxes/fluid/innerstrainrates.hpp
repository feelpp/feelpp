
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


// /!\ Pourvoir chosir la bonne dimension 
using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;
using mesh_t = Mesh<Simplex<2, 1>>;







template<std::size_t residualType, typename FluidMechanics, typename Pchv, typename space>
auto
innerStrainRates(FluidMechanics &t, Pchv &u1, Pchv &u2, space &Xh)
{ 
    //auto Xh = Pdh<1>(u1.functionSpace()->mesh());;
    auto r = Xh->element();

    std::cout << "Computing strain rates" << std::endl;
    auto e1 = sym(gradv(u1));
    auto e2 = sym(gradv(u2));
  
    std::cout << "Computing inner product of strain rates" << std::endl;
    auto e1_inner_e2 = inner(e1,e2);

    auto nDof_u1 = u1.functionSpace()->nDof();
    auto nDof_u2 = u2.functionSpace()->nDof();
    auto size_u1 = u1.size();
    auto size_u2 = u2.size();
    std::cout << "nDof of u1: " << nDof_u1 << std::endl;
    std::cout << "nDof of u2: " << nDof_u2 << std::endl;
    std::cout << "size of u1: " << size_u1 << std::endl;
    std::cout << "size of u2: " << size_u2 << std::endl;

    auto mesh = Xh->mesh();//u1.functionSpace()->mesh();
    auto numGlobalPoints = mesh->numGlobalPoints();
    auto numGlobalElements = mesh->numGlobalElements(); 
    auto numGlobalFaces = mesh->numGlobalFaces();
    auto numGlobalEdges = mesh->numGlobalEdges();
    auto measure = mesh->measure();
    std::cout << "Mesh information:" << std::endl;
    std::cout << "Number of global points: " << numGlobalPoints << std::endl;
    std::cout << "Number of global elements: " << numGlobalElements << std::endl;
    std::cout << "Number of global faces: " << numGlobalFaces << std::endl;
    std::cout << "Number of global edges: " << numGlobalEdges << std::endl;
    std::cout << "Measure of the mesh: " << measure << std::endl;

    std::cout << "Evaluating inner product of strain rates" << std::endl;
    r.on( _range=elements(mesh), _expr=e1_inner_e2 );
    std::cout << "Inner product of strain rates evaluated" << std::endl;
    return r;

}


template<std::size_t residualType, typename FluidMechanics, typename Pdh>
void
saveinnerStrainRates(FluidMechanics &t, Pdh &e1_inner_e2, const std::string& path)
{
    e1_inner_e2.saveHDF5( path );
}

template<std::size_t residualType, typename FluidMechanics, typename Pchv>
void
saveVelocity(FluidMechanics &t, Pchv &u, const std::string& path)
{
    u.saveHDF5( path );
}
