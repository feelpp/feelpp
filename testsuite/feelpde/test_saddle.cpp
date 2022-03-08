#define BOOST_TEST_MODULE test_saddle
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelpde/preconditionerblockms.hpp>
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
    ( "mu", po::value<double>()->default_value( 1. ), "rot(1/mu rot(u) + grad(p) = j; div(u) = 0" )
    ( "ms.11.setAlphaBeta", po::value<bool>()->default_value( false ), "[ams] set Alpha and Beta" )
    ;
    return opts.add( Feel::feel_options() )
        .add(Feel::blockms_options("ms"))
        .add(Feel::backend_options("ms"))
        ;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_saddle"), makeOptions() )

BOOST_AUTO_TEST_SUITE( test_saddle )
BOOST_AUTO_TEST_CASE( test_0 )
{
    typedef double value_type;
typedef Simplex<3> convex_type;
typedef Mesh<convex_type> mesh_type;
typedef std::shared_ptr<mesh_type> mesh_ptrtype;
//! Hcurl space
typedef Nedelec<0,NedelecKind::NED1 > curl_basis_type;
typedef FunctionSpace<mesh_type, Feel::detail::bases<curl_basis_type>, value_type,Feel::Periodicity<Feel::NoPeriodicity>, Feel::mortars<Feel::NoMortar>> curl_space_type;
typedef std::shared_ptr<curl_space_type> curl_space_ptrtype;
typedef typename curl_space_type::element_type curl_element_type;
//! Pch space
typedef Lagrange<1, Scalar> lag_basis_type; 
typedef FunctionSpace<mesh_type, Feel::detail::bases<lag_basis_type>, value_type,Feel::Periodicity<Feel::NoPeriodicity>, Feel::mortars<Feel::NoMortar>> lag_space_type;
typedef std::shared_ptr<lag_space_type> lag_space_ptrtype;
typedef typename lag_space_type::element_type lag_element_type;
//! Comp space 
typedef FunctionSpace<mesh_type, Feel::detail::bases<curl_basis_type,lag_basis_type>, value_type,Feel::Periodicity<Feel::NoPeriodicity>, Feel::mortars<Feel::NoMortar>> comp_space_type;
typedef std::shared_ptr<comp_space_type> comp_space_ptrtype;
typedef typename comp_space_type::element_type comp_element_type;

  auto mesh = loadMesh(_mesh = new Mesh<Simplex<FEELPP_DIM>> );

  auto Xh   = comp_space_type::New( mesh );

  auto U = Xh->element();
  auto V = Xh->element();
  auto u = U.element<0>(); 
  auto p = U.element<1>(); 
  auto v = U.element<0>(); 
  auto q = U.element<1>(); 

  auto a = form2(_test=Xh, _trial=Xh);
  auto l = form1( _test=Xh );
    
  ModelProperties model;

  a = integrate(_range=elements(mesh),
      _expr = 
        doption("mu")*trans(curlt_op(u))*curl_op(v)
      + inner(trans(id (v)), gradt(p)) // grad(p)
      + inner(trans(idt(u)), grad (q))  // div(u) = 0
      );
  l = integrate(_range=elements(mesh),_expr = inner(expr<3,1>(soption("functions.j")),id(v)));
  
  map_vector_field<FEELPP_DIM,1,2> m_dirichlet_u { model.boundaryConditions().getVectorFields<FEELPP_DIM> ( "u", "Dirichlet" )};
  map_vector_field<FEELPP_DIM,1,2> m_weak_u { model.boundaryConditions().getVectorFields<FEELPP_DIM> ( "u", "Weakdir" ) };
  
  map_scalar_field<2> m_dirichlet_p { model.boundaryConditions().getScalarFields ( "phi", "Dirichlet" )};
  map_scalar_field<2> m_weak_p { model.boundaryConditions().getScalarFields ( "phi", "Weakdir" )};
  /*
   * BC(u)
   */
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
  /*
   * BC(p)
   */
  for(auto it : m_dirichlet_p)
  {
    LOG(INFO) << it.second.first << " on " << it.first << "\n";
  a += on(_range=markedfaces(mesh,it.first),
      _expr=it.second.first,
      _rhs=l,
      _element=p,
      _type="elimination_symmetric"
      );
  }

  auto prec = preconditioner(_backend=backend(_name="ms"),_pc=pcTypeConvertStrToEnum(soption("ms.pc-type")), _prefix="ms",_matrix=a.matrixPtr());
  if(soption("ms.pc-type") == "blockms" )
  {
    auto precBMS = std::make_shared<PreconditionerBlockMS<comp_space_type>>(
        U.functionSpace(),
        model,
        "ms",
        a.matrixPtr(), 1.);
    prec->attachInHousePreconditioners("blockms",precBMS);
    //preconditioner(_backend=backend(_name="ms"), _pc=FEELPP_BLOCKMS_PRECOND, _prefix="ms")->attachInHousePreconditioners("blockms",prec);
  }
  a.solveb(_rhs=l, _solution=U, _backend=backend(_name="ms"),_prec=prec);
  if(boption("exporter.export"))
  {
    auto ue = vf::project(_space=Xh->functionSpace<0>(),_range=elements(mesh), _expr=expr<3,1>(soption("functions.u")));
    auto e=exporter(_mesh = mesh );
    e->add("u",u);
    e->add("ue",ue);
    e->add("p",p);
    e->save();
  }


}
BOOST_AUTO_TEST_SUITE_END()
