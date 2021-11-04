#define BOOST_TEST_MODULE test_ams
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelvf/vf.hpp>

#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv

using namespace Feel;
inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_precAFP" );
    opts.add_options()
    ( "regul", po::value<double>()->default_value( 1. ), "rot(1/mu rot(u) + regul u = j" )
    ( "mu", po::value<double>()->default_value( 1. ), "rot(1/mu rot(u) + regul u = j" )
    ( "ms.setAlphaBeta", po::value<bool>()->default_value( false ), "[ams] set Alpha and Beta" )
    ;
    return opts.add( Feel::feel_options() )
        .add(Feel::backend_options("ms"));
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_ams"), makeOptions() )

BOOST_AUTO_TEST_SUITE( test_ams )
BOOST_AUTO_TEST_CASE( test_0 )
{

  auto mesh = loadMesh(_mesh = new Mesh<Simplex<FEELPP_DIM>> );

  auto Xh = Ned1h<0>( mesh );
  auto XhL = Pch<1>( mesh );
  auto XhLv= Pchv<1>( mesh );

  auto u = Xh->element();
  auto v = Xh->element();

  auto e=exporter(_mesh = mesh );

  auto a = form2(_test=Xh, _trial=Xh);
  auto l = form1( _test=Xh );
  
  ModelProperties model;
  map_vector_field<FEELPP_DIM,1,2> m_dirichlet_u { model.boundaryConditions().getVectorFields<FEELPP_DIM> ( "u", "Dirichlet" )};
  map_vector_field<FEELPP_DIM,1,2> m_weak_u { model.boundaryConditions().getVectorFields<FEELPP_DIM> ( "u", "Weakdir" ) };

  a = integrate(_range=elements(mesh),
      _expr = 
        doption("mu")*trans(curlt_op(u))*curl_op(v)
      + doption("regul")*inner(idt(u), id(u))
      );
  l = integrate(_range=elements(mesh),
      _expr = inner(expr<3,1>(soption("functions.j")),id(v)));
  for(auto it : m_dirichlet_u)
  {
    LOG(INFO) << it.second.first << " on " << it.first << "\n";
  a += on(_range=markedfaces(mesh,it.first),
      _expr=it.second.first,
      _rhs=l,
      _element=u,
      _type="elimination_symmetric"
      );
  }
  auto prec = preconditioner(_pc=pcTypeConvertStrToEnum(soption("ms.pc-type")),
      _backend=backend(_name="ms"),
      _prefix="ms",
      _matrix=a.matrixPtr()
      );
  if(soption("ms.pc-type") == "ams")
  {
    auto Igrad   = Grad( _domainSpace=XhL, _imageSpace=Xh);

    auto ozz = Xh->element();
    auto zoz = Xh->element();
    auto zzo = Xh->element();
    ozz.on(_range=elements(Xh->mesh()),_expr=vec(cst(1),cst(0),cst(0)));
    zoz.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(1),cst(0)));
    zzo.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
    auto M_ozz = backend(_name="ms")->newVector(Xh); *M_ozz = ozz; M_ozz->close();
    auto M_zoz = backend(_name="ms")->newVector(Xh); *M_zoz = zoz; M_zoz->close();
    auto M_zzo = backend(_name="ms")->newVector(Xh); *M_zzo = zzo; M_zzo->close();


    prec->attachAuxiliarySparseMatrix("G",Igrad.matPtr());
    prec->attachAuxiliaryVector("Px",M_ozz);
    prec->attachAuxiliaryVector("Py",M_zoz);
    prec->attachAuxiliaryVector("Pz",M_zzo);
    
    if(boption(_name="setAlphaBeta",_prefix="ms"))
    {
      auto a_alpha = form2(_test=XhLv, _trial=XhLv);
      auto b_alpha = form1(_test=XhLv);
      a_alpha = integrate(_range=elements(XhLv->mesh()), _expr=cst(1.)/doption("mu")*inner(gradt(u), grad(u)));
      a_alpha += on(_range=boundaryfaces(XhLv->mesh()),_element=u, _expr=expr<3,1>(soption("functions.u")), _rhs=b_alpha, _type="elimination_symmetric");
      a_alpha.matrixPtr()->close();
      prec->attachAuxiliarySparseMatrix("a_alpha",a_alpha.matrixPtr());

      auto uu = XhL->element();
      auto a_beta = form2(_test=XhL, _trial=XhL);
      auto b_beta = form1(_test=XhL);
      a_beta = integrate(_range=elements(XhL->mesh()), _expr=doption("regul")*inner(grad(uu),gradt(uu)));
      
      a_beta += on(_range=boundaryfaces(XhL->mesh()),_element=uu, _expr=cst(0.), _rhs=b_beta, _type="elimination_symmetric");
      a_beta.matrixPtr()->close();
      prec->attachAuxiliarySparseMatrix("a_beta",a_beta.matrixPtr());
    }
    if(doption("regul") == 0)
      prec->attachAuxiliarySparseMatrix("a_beta",NULL);
  }
  a.solveb(_rhs=l, _solution=u, _backend=backend(_name="ms"), _prec=prec);
  
  if(boption("exporter.export"))
  {
    auto ue = vf::project(_space=Xh,_range=elements(mesh), _expr=expr<3,1>(soption("functions.u")));
    auto e=exporter(_mesh = mesh );
    e->add("u",u);
    e->add("ue",ue);
    e->save();
  }

}
BOOST_AUTO_TEST_SUITE_END()
