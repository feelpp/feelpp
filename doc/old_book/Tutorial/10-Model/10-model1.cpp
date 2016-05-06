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
