/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-28

  Copyright (C) 2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file magnetostatic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-28
 */
#include <feel/feel.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpde/boundaryconditions.hpp>

#if FEELPP_DIM==2
#define curl_op curlx
#define culrt_op curlxt
#else
#define curl_op curl
#define curlt_op curlt
#endif

int main(int argc, char**argv )
{
    using namespace Feel;
    po::options_description magnetooptions( "Magnetostatic options" );
    magnetooptions.add_options()
    ( "weakdirichlet", po::value<bool>()->default_value( false ), "false if strong Dirichlet, true if weak Dirichlet treatment" )
    ;

    /// [marker1]
    Environment env( _argc=argc, _argv=argv,
                     _desc=magnetooptions,
                     _about=about(_name="magnetostatic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef Mesh<Simplex<FEELPP_DIM>> mesh_type;
    typedef Nedelec<0, NedelecKind::NED1> curl_basis_type;
    typedef Lagrange<1,Scalar> lag_basis_type;
    typedef FunctionSpace<mesh_type, bases<curl_basis_type, lag_basis_type> > space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;

    BoundaryConditions bcs;
    map_vector_field<FEELPP_DIM,1,2> m_dirichlet { bcs.getVectorFields<FEELPP_DIM> ( "u", "Dirichlet")};
    map_vector_field<FEELPP_DIM,1,2> m_weak { bcs.getVectorFields<FEELPP_DIM> ( "u", "Dirichlet_w")};
    map_vector_field<FEELPP_DIM,1,2> m_phi { bcs.getVectorFields<FEELPP_DIM> ( "phi", "Dirichlet")};

    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = space_type::New(mesh);

    auto a = form2( _trial=Xh, _test=Xh );
    auto f = expr<FEELPP_DIM,1>(soption("functions.f"), "f"); // rhs
    //std::cout << "curl(curl(E)) + E = " << curl(curl(e)) << "\n";
    auto U = Xh->element();
    auto u = U.element<0>();
    auto phi =  U.element<1>();
    auto V = Xh->element(e);
    auto v = V.element<0>();
    auto psi =  V.element<1>();

    double penaldir=doption("parameters.d");

    std::cout << "Weakdirichlet = " << (boption("weakdirichlet") ? "true":"false") << std::endl;
    std::cout << "Penaldir      = " <<  penaldir << std::endl;

    a = integrate(_range=elements(mesh), _expr=c*trans(idt(u))*id(v)+curlt_op(u)*curl_op(v));

    if ( boption("weakdirichlet") ) //weak Dirichlet
        a += integrate(boundaryfaces(mesh), -curlt_op(u)*(cross(N(),id(u)) )
                       - curl_op(u)*(cross(N(),idt(u)) )
                       + penaldir*trans(cross(idt(u),N()))*cross(id(u),N())/hFace() );
    auto l = form1( _test=Nh );
    l = integrate( _range=elements(mesh), _expr=trans(f)*id(v));
    if ( boption("weakdirichlet") ) //weak Dirichlet
    {
      l += integrate(boundaryfaces(mesh), - curl_op(u)*(cross(N(),e) )
                   + penaldir*trans(cross(e,N()))*cross(id(u),N())/hFace() );
    }
    else //strong Dirichlet
    {
#if 0
        for(auto const & it:  m_dirichlet)
        {
            a += on( _range=markedfaces(mesh,it.first), _rhs=l, _element=u, _expr=it.second);
        }
#endif
        a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u,_expr=e );
    }

    a.solve(_rhs=l,_solution=u);
    auto erru = normL2( _range=elements(mesh), _expr=idv(u)-e );
    auto errinterp = normL2( _range=elements(mesh), _expr=idv(v)-e );
    if ( Environment::isMasterRank() )
    {
        std::cout << "L2 error (u-e) = " << erru << "\n";
        std::cout << "L2 error (interp_curl(e) -e) = " << errinterp << "\n";
    }

    auto b = form2( _test=Xh, _trial=Xh );
    b = integrate(_range=elements(mesh), _expr=trans(idt(w))*id(w));
    auto ll = form1( _test=Xh );
    ll = integrate( _range=elements(mesh), _expr=trans(idv(u))*id(w));
    b.solve( _solution=z, _rhs=ll, _rebuild=true );

    auto ex = exporter( _mesh=mesh );
    ex->add( "u", u );
    ex->add( "Ihe", v );
    ex->add( "e", w );
    ex->add( "ul2", z );
    ex->save();
}
