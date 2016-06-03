#include <feel/feel.hpp>
#include <feel/feeldiscr/ned1h.hpp>

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

inline
AboutData
makeAbout()
{
#if FEELPP_DIM==2
    AboutData about( "precAFP2D" ,
                     "precAFP2D" ,
                     "0.1",
                     "test precAFP2D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#else
    AboutData about( "precAFP3D" ,
                     "precAFP3D" ,
                     "0.1",
                     "test precAFP3D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#endif

    about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );

    return about;
}

int main( int argc, char** argv )
{
  // Initialize Feel++ Environment
  Environment env( _argc=argc, _argv=argv,
      _desc=makeOptions(),
      _about=makeAbout()
      );

  /// Todo : add option regul & backend;

  auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

  auto Xh = Ned1h<0>( mesh );
  auto XhL = Pch<1>( mesh );
  auto XhLv= Pchv<1>( mesh );

  auto u = Xh->element();
  auto v = Xh->element();

  auto e=exporter(_mesh = mesh );

  auto a = form2(_test=Xh, _trial=Xh);
  auto l = form1( Xh );

  a = integrate(_range=elements(mesh),
      _expr = 
        doption("mu")*trans(curlt_op(u))*curl_op(v)
      + doption("regul")*inner(idt(u), id(u))
      );
  l = integrate(_range=elements(mesh),
      _expr = inner(expr<3,1>(soption("functions.j")),id(v)));
  a += on(_range=boundaryfaces(mesh),
      _expr=expr<3,1>(soption("functions.u")),
      _rhs=l,
      _element=u,
      _type="elimination_symmetric"
      );

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
      prec->attachAuxiliarySparseMatrix("a_alpha",a_alpha.matrixPtr());

      auto uu = XhL->element();
      auto a_beta = form2(_test=XhL, _trial=XhL);
      auto b_beta = form1(_test=XhL);
      a_beta = integrate(_range=elements(XhL->mesh()), _expr=doption("regul")*inner(grad(uu),gradt(uu)));
      
      a_beta += on(_range=boundaryfaces(XhL->mesh()),_element=uu, _expr=cst(0.), _rhs=b_beta, _type="elimination_symmetric");
      prec->attachAuxiliarySparseMatrix("a_beta",a_beta.matrixPtr());
    }
    if(doption("regul") == 0)
      prec->attachAuxiliarySparseMatrix("a_beta",NULL);
  }
  a.solveb(_rhs=l, _solution=u, _backend=backend(_name="ms"), _prec=prec);

}
