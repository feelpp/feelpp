#include <feel/feel.hpp>

#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv

using namespace Feel;
int main( int argc, char** argv )
{
  // Initialize Feel++ Environment
  Environment env( _argc=argc, _argv=argv,
      _desc=feel_options(),
      _about=about( _name="test_ams" ,
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org" ) );

  /// Todo : add option regul & backend;

  auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

  auto Xh = Ned1h<0>( mesh );
  auto XhL = Pchv<1>( mesh );

  auto u = Xh->element();
  auto v = Xh->element();

  auto e=exporter(_mesh = mesh );

  auto a = form2(_Xh, Xh);
  auto l = form1( Xh );

  a = integrate(_range=elements(mesh),
      _expr = trans(curlt_op(u))*curl_op(v)
      + doption("regul")*inner(idt(u), id(u))
      );
  l = integrate(_range=elements(mesh),
      _expr = inner(expr<3,1>(soption("functions.j")),id(v)));
  a += on(_range=boundaryelements(mesh),
      _expr=expr<3,1>(soption("functions.u")));

  auto prec = preconditioner(_pc=pcTypeConvertStrToEnum(soption("myBackend.pc-type")),
      _backend=backend(_name="myBackend"),
      _prefix="myBackend",
      _matrix=a.matrixPtr()
      );
  if(soption("myBackend.pc-type") == "ams")
    auto Igrad   = Grad( _domainSpace=XhL, _imageSpace=Xh);

  auto ozz = Xh->element();
  auto zoz = Xh->element();
  auto zzo = Xh->element();
  ozz.on(_range=elements(Xh->mesh()),_expr=vec(cst(1),cst(0),cst(0)));
  zoz.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(1),cst(0)));
  zzo.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
  M_ozz = backend(_name="myBackend")->newVector(Xh); *M_ozz = ozz; M_ozz->close();
  M_zoz = backend(_name="myBackend")->newVector(Xh); *M_zoz = zoz; M_zoz->close();
  M_zzo = backend(_name="myBackend")->newVector(Xh); *M_zzo = zzo; M_zzo->close();


  prec->attachAuxiliarySparseMatrix("G",Igrad.matPtr());
  prec->attachAuxiliaryVector("Px",M_ozz);
  prec->attachAuxiliaryVector("Py",M_zoz);
  prec->attachAuxiliaryVector("Pz",M_zzo);
  //prec->auxiliarySparseMatrix("G")->printMatlab("G.m");
  //prec->auxiliaryVector("Px")->printMatlab("Px.m");
  //prec->auxiliaryVector("Py")->printMatlab("Py.m");
  //prec->auxiliaryVector("Pz")->printMatlab("Pz.m");
  if(boption(_name="setAlphaBeta",_prefix="myBackend"))
  {
    auto a_alpha = form2(_test=XhLv, _trial=XhLv);
    auto b_alpha = form1(_test=XhLv);
    for(auto it : model.materials() )
      a_alpha += integrate(_range=markedelements(XhLv->mesh(),marker(it)), _expr=cst(1.)/idv(M_mu_r)*inner(gradt(u), grad(u)));
    for(auto const & it : m_dirichlet)
      a_alpha += on(_range=markedfaces(XhLv->mesh(),it.first),_element=u, _expr=it.second, _rhs=b_alpha, _type="elimination_symmetric");
    prec->attachAuxiliarySparseMatrix("a_alpha",a_alpha.matrixPtr());

    auto uu = XhL->element();
    auto a_beta = form2(_test=XhL, _trial=XhL);
    auto b_beta = form1(_test=XhL);
    for(auto it : model.materials() )
    {
      std::string key = "Materials."+marker(it)+".mu_opt";
      if(doption("relax") > 0.)
        a_beta += integrate(_range=markedelements(XhL->mesh(),marker(it)), _expr=M_mu_0*doption("relax")*inner(grad(uu),gradt(uu)));
      else
      {
        std::string key = "Materials."+marker(it)+".mu_opt";
        a_beta += integrate(_range=markedelements(XhL->mesh(),marker(it)), _expr=(M_mu_0/expr(model.getEntry(key)))*inner(grad(uu),gradt(uu)));
      }
    }
    for(auto const & it : m_dirichlet)
      a_beta += on(_range=markedfaces(XhL->mesh(),it.first),_element=uu, _expr=cst(0.), _rhs=b_beta, _type="elimination_symmetric");
    prec->attachAuxiliarySparseMatrix("a_beta",a_beta.matrixPtr());
  }
  if(doption("relax") == 0)
    prec->attachAuxiliarySparseMatrix("a_beta",NULL);
}
a.solveb(_rhs=l, _solution=u, _backend=backend(_name="myBackend"), _prec=prec);

}
