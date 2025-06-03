
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


// /!\ Pourvoir chosir la bonne dimension 
using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;
using mesh_t = Mesh<Simplex<2, 1>>;







template<std::size_t residualType, typename FluidMechanics, typename Pchv>
auto
innerStrainRates(FluidMechanics &t, Pchv &u1, Pchv &u2)
{ 
    auto Xh = Pdh<1>(u1.functionSpace()->mesh());;
    auto r = Xh->element();

    std::cout << "Computing strain rates" << std::endl;
    auto e1 = sym(gradv(u1));
    auto e2 = sym(gradv(u2));
  
    std::cout << "Computing inner product of strain rates" << std::endl;
    auto e1_inner_e2 = inner(e1,e2);


    std::cout << "Evaluating inner product of strain rates" << std::endl;
    r.on( _range=elements(u1.functionSpace()->mesh()), _expr=e1_inner_e2 );
    std::cout << "Inner product of strain rates evaluated" << std::endl;
    return r;

}


template<std::size_t residualType, typename FluidMechanics, typename Pdh>
void
saveinnerStrainRates(FluidMechanics &t, Pdh &e1_inner_e2, const std::string& path)
{
    e1_inner_e2.saveHDF5( path );
}
