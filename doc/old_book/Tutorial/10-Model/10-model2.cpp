#include <feel/feel.hpp>
#include <feel/feelmodels/modelproperties.hpp>
int main(int argc, char**argv )
{
  using namespace Feel;
  po::options_description laplacianoptions( "Laplacian options" );
  laplacianoptions.add_options()
    ("myVerbose", po::value< bool >()-> default_value( true ), "Display information during execution")
    ;
  Environment env( _argc=argc, _argv=argv,
      _desc=laplacianoptions,
      _about=about(_name="aniso_laplacian",
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org"));
  ModelProperties model; // Will load --mod-file
  map_scalar_field<2> bc_u { model.boundaryConditions().getScalarFields<2>("heat","dirichlet") };
  ModelMaterials materials = model.materials();
  if(boption("myVerbose") && Environment::isMasterRank() )
    std::cout << "Model " << Environment::expand( soption("mod-file")) << " loaded." << std::endl;
  auto f = expr( soption(_name="functions.f"), "f" );
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<MODEL_DIM>>);
  auto Vh = Pch<2>( mesh );
  auto u = Vh->element();
  auto v = Vh->element();
  auto k11 = Vh->element();
  auto k12 = Vh->element();
  auto k22 = Vh->element();
#if MODEL_DIM == 3
  auto k13 = Vh->element();
  auto k23 = Vh->element();
  auto k33 = Vh->element();
#endif
  auto a = form2( _trial=Vh, _test=Vh);
  auto l = form1( _test=Vh );
  l = integrate(_range=elements(mesh),_expr=f*id(v));
  for(auto it : materials)
  {
    auto mat = material(it);
    if(boption("myVerbose") && Environment::isMasterRank() )
      std::cout << "[Materials] - Laoding data for " << it.second.name() << " that apply on marker " << it.first  << " with diffusion coef [" 
#if MODEL_DIM == 3
        << "[" << it.second.k11() << "," << it.second.k12() << "," << it.second.k13() << "],"
        << "[" << it.second.k12() << "," << it.second.k22() << "," << it.second.k23() << "],"
        << "[" << it.second.k13() << "," << it.second.k23() << "," << it.second.k33() << "]]" 
#else
        << "[" << it.second.k11() << "," << it.second.k12() << "],"
        << "[" << it.second.k12() << "," << it.second.k22() << "]]"
#endif
        << std::endl;
    k11.on(_range=markedelements(mesh,it.first),_expr=cst(it.second.k11()));
    k12.on(_range=markedelements(mesh,it.first),_expr=cst(it.second.k12()));
    k22.on(_range=markedelements(mesh,it.first),_expr=cst(it.second.k22()));
#if MODEL_DIM == 3
    k13 += vf::project(_space=Vh,_range=markedelements(mesh,marker(it)),_expr=mat.k13());
    k23 += vf::project(_space=Vh,_range=markedelements(mesh,marker(it)),_expr=mat.k23());
    k33 += vf::project(_space=Vh,_range=markedelements(mesh,marker(it)),_expr=mat.k33());
#endif
  }
#if MODEL_DIM == 2
  a += integrate(_range=elements(mesh),_expr=inner(mat<MODEL_DIM,MODEL_DIM>(idv(k11), idv(k12), idv(k12), idv(k22) )*trans(gradt(u)),trans(grad(v))) );
#else
  a += integrate(_range=elements(mesh),_expr=inner(mat<MODEL_DIM,MODEL_DIM>(idv(k11), idv(k12), idv(k13), idv(k12), idv(k22), idv(k23), idv(k31), idv(k32), idv(k33))*trans(gradt(u)),trans(grad(v))) );
#endif
  for(auto it : bc_u){
    if(boption("myVerbose") && Environment::isMasterRank() )
      std::cout << "[BC] - Applying " << it.second << " on " << it.first << std::endl;
    a+=on(_range=markedfaces(mesh,it.first), _rhs=l, _element=u, _expr=it.second );
  }
  a.solve(_rhs=l,_solution=u);
  auto e = exporter( _mesh=mesh );
  for(int i = 0; i < 3; i ++){
    for(auto const &it : model.postProcess()["Fields"] )
    {
      if(it == "diffused") 
        e->step(i)->add("diffused",u);
      else if(it == "k11")
        e->step(i)->add("k11",k11);
      else if(it == "k12")
        e->step(i)->add("k12",k12);
      else if(it == "k11")
        e->step(i)->add("k22",k22);
#if MODEL_DIM == 3
      else if(it == "k13")
        e->step(i)->add("k13",k13);
      else if(it == "k11")
        e->step(i)->add("k23",k23);
      else if(it == "k33")
        e->step(i)->add("k33",k33);
#endif
    }
    e->save();
  }
  return 0;
}
