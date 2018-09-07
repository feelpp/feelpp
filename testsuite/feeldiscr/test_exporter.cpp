/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 janv. 2016
 
 Copyright (C) 2016 Feel++ Consortium
 
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

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/filters.hpp>
#include <feel/feeldiscr/projector.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace Feel
{
  inline
    Feel::po::options_description
    makeOptions()
    {
      return Feel::feel_options() ;
    }

  template<int Dim, int Order>
  class cv : public Simget
  {
      typedef Simget super;
  public:
      typedef double value_type;

      typedef Backend<value_type> backend_type;
      typedef std::shared_ptr<backend_type> backend_ptrtype;
      //** Mesh **
      typedef Mesh<Simplex<Dim>> mesh_type;
      typedef std::shared_ptr<mesh_type> mesh_ptrtype;

      //** Basis **
      //Scalar
      typedef bases<Lagrange<Order,Scalar>> basis_type;
      typedef std::shared_ptr<basis_type> basis_ptrtype;

      //Scalar P0
      typedef bases<Lagrange<0,Scalar,Discontinuous>> basis_pm1_type;
      typedef std::shared_ptr<basis_pm1_type> basis_pm1_ptrtype;

      //Vectorial P0
      typedef bases<Lagrange<0,Vectorial,Discontinuous>> basis_pm1vec_type;
      typedef std::shared_ptr<basis_pm1vec_type> basis_pm1vec_ptrtype;

      //Scalar space
      typedef FunctionSpace<mesh_type,basis_type> space_type;
      typedef std::shared_ptr<space_type> space_ptrtype;

      //Scalar space P0
      typedef FunctionSpace<mesh_type,basis_pm1_type> space_pm1_type;
      typedef std::shared_ptr<space_pm1_type> space_pm1_ptrtype;

      //Vectorial space P0
      typedef FunctionSpace<mesh_type,basis_pm1vec_type> space_pm1vec_type;
      typedef std::shared_ptr<space_pm1vec_type> space_pm1vec_ptrtype;

      /**
       * Constructor
       */
      cv() : super()
          {}

      void run();
  };

template<int Dim, int Order>
void
cv<Dim,Order>::run()
{
    auto mesh = loadMesh(_mesh = new mesh_type ,_update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    auto Xh = space_type::New( mesh );
    auto Xhm1 = space_pm1_type::New( mesh );
    auto Xhm1vec = space_pm1vec_type::New( mesh );

    auto phi = Xh->element();
    auto gphi = Xhm1vec->element();
    auto gphi0 = Xhm1->element();
    phi = vf::project(Xh,elements(mesh), Px()*Px() + Py()*Py() );
    // P0 vector
    gphi = vf::project(Xhm1vec, elements(mesh), trans( gradv(phi) ) );
    // P0 scalar
    auto gphi_compx = gphi.comp( ( ComponentType )0 );
    gphi0 = vf::project(Xhm1, elements(mesh), idv(gphi_compx) );

    auto   e = exporter(_mesh = mesh, _name = "myExporter");
    e->add( "phi", phi );
    e->add( "grad_phi", gphi );
    e->add( "grad_phi0", gphi0 );
    e->save();
}
} // Feel

int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="test_exporter",
                                  _author="Cemosis",
                                  _email="vincent.huber@cemosis.fr"));

    Application app2;
    app2.add(new cv<2,2>);
    app2.run();
}

