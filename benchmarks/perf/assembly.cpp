/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Pierre Jolivet <pierre.jolivet@imag.fr>
       Date: 2014-12-16

  Copyright (C) 2014 Université de Strasbourg
  Copyright (C) 2014 Université de Grenoble

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
   \file assembly.cpp
   \author Pierre Jolivet <pierre.jolivet@imag.fr>
   \date 2014-12-16
 */
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#define MAX_ORDER 30

namespace Feel
{
template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type, uint16_type OrderBis, template<uint16_type> class TypeBis, typename std::enable_if<!std::is_same<Type<Dim>, Vectorial<Dim>>::value && OrderBis == MAX_ORDER && std::is_same<Type<Dim>, Scalar<Dim>>::value>::type* = nullptr>
static inline void assemble(boost::shared_ptr<Mesh<Simplex<Dim>>>& mesh, double* vec) {
    boost::timer time;
    auto Vh = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>>>::New(_mesh = mesh);
    vec[2] = time.elapsed();
    vec[0] = mesh->numGlobalElements();
    vec[1] = Vh->nDof();
    time.restart();
    auto A = backend()->newMatrix(Vh, Vh);
    auto v = Vh->element();
    auto u = Vh->element();
    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(mesh), _expr = gradt(u) * trans(grad(v)));
    vec[3] = time.elapsed();
    time.restart();
    auto f = backend()->newVector(Vh);
    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(mesh), _expr = id(v));
    a += on(_range = markedfaces(mesh, "Dirichlet"), _rhs = l, _element = u, _expr = cst(0.0));
    vec[4] = time.elapsed();
}
template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type, uint16_type OrderBis, template<uint16_type> class TypeBis, typename std::enable_if<std::is_same<Type<Dim>, Vectorial<Dim>>::value && OrderBis == MAX_ORDER && std::is_same<TypeBis<Dim>, Scalar<Dim>>::value>::type* = nullptr>
static inline void assemble(boost::shared_ptr<Mesh<Simplex<Dim>>>& mesh, double* vec) {
    boost::timer time;
    auto Vh = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>>>::New(_mesh = mesh);
    vec[2] = time.elapsed();
    vec[0] = mesh->numGlobalElements();
    vec[1] = Vh->nDof();
    auto E = 1e+8;
    auto nu = 0.25;
    auto mu = E / (2 * (1 + nu));
    auto lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    time.restart();
    auto A = backend()->newMatrix(Vh, Vh);
    auto v = Vh->element();
    auto u = Vh->element();
    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(mesh),
                  _expr = lambda * divt(u) * div(v) +
                          2 * mu * trace(trans(sym(gradt(u))) * sym(grad(u))));
    vec[3] = time.elapsed();
    time.restart();
    auto f = backend()->newVector(Vh);
    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(mesh), _expr = -1e+3 * trans(oneY()) * id(v));
    a += on(_range = markedfaces(mesh, "Dirichlet"), _rhs = l, _element = u, _expr = zero<2, 1>());
    vec[4] = time.elapsed();
}
template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type, uint16_type OrderBis, template<uint16_type> class TypeBis, typename std::enable_if<std::is_same<Type<Dim>, Vectorial<Dim>>::value && OrderBis != MAX_ORDER && std::is_same<TypeBis<Dim>, Scalar<Dim>>::value>::type* = nullptr>
static inline void assemble(boost::shared_ptr<Mesh<Simplex<Dim>>>& mesh, double* vec) {
    boost::timer time;
    auto Vh = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>, Lagrange<OrderBis, TypeBis>>>::New(_mesh = mesh,
                                                                                   _worldscomm = std::vector<WorldComm>(2, mesh->worldComm()),
                                                                                   _extended_doftable = std::vector<bool>(2, false));
    vec[2] = time.elapsed();
    vec[0] = mesh->numGlobalElements();
    vec[1] = Vh->nDof();
    auto A = backend()->newMatrix(Vh, Vh);
    auto f = backend()->newVector(Vh);
    auto U = Vh->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(mesh),
                  _expr = trans(oneY()) * id(u));
    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(mesh),
                  _expr = trace(gradt(u) * trans(grad(u))) - div(u) * idt(p) - divt(u) * id(p));
    a += on(_range = markedfaces(mesh, "Dirichlet"), _rhs = l, _element = u, _expr = zero<2, 1>());
    
}

template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type, uint16_type OrderBis = MAX_ORDER, template<uint16_type> class TypeBis = Scalar>
class Assembly : public Simget
{
public:
    void run();
}; // Assembly

template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type, uint16_type OrderBis, template<uint16_type> class TypeBis>
void
Assembly<Dim, Order, Type, OrderBis, TypeBis>::run()
{
    double hSize = doption("gmsh.hsize");
    int level = std::max(doption("parameters.l"), 1.0);
    std::vector<double> timers(5 * level, 0);
    for(int i = 0; i < level; ++i) {
        boost::shared_ptr<Mesh<Simplex<Dim>>> mesh;
        mesh = createGMSHMesh(_mesh = new Mesh<Simplex<Dim>>,
                              _update = MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                              _desc = domain(_name = "hypercube_" + std::to_string(i + 1), _shape = "hypercube", _h = hSize / std::pow(2.0, i),
                                             _xmin = 0.0, _xmax = 10.0,
                                             _ymin = 0.0, _ymax = 1.0,
                                             _zmin = 0.0, _zmax = 1.0));
        mesh->addMarkerName("Dirichlet", Dim == 2 ? 1 : 19, Dim == 2 ? 1 : 2);
        assemble<Dim, Order, Type, OrderBis, TypeBis>(mesh, &(timers[5 * i]));
    }
    if(!std::is_same<Type<Dim>, Vectorial<Dim>>::value && OrderBis == MAX_ORDER && std::is_same<Type<Dim>, Scalar<Dim>>::value)
        std::cout << "hsize\t\tnelements\tnDof\t\tFunctionSpace\tform2\t\tform1" << std::endl;
    for(int i = 0; i < level; ++i) {
        std::cout.width(16);
        std::cout << std::left << hSize / std::pow(2.0, i);
        std::cout.width(16);
        std::cout << std::left << int(timers[5 * i + 0]);
        std::cout.width(16);
        std::cout << std::left << int(timers[5 * i + 1]);
        std::cout.width(16);
        std::cout << std::left << timers[5 * i + 2];
        std::cout.width(16);
        std::cout << std::left << timers[5 * i + 3];
        std::cout << std::left << timers[5 * i + 4] << std::endl;
    }
} // Assembly::run

} // Feel

int main(int argc, char** argv) {
    using namespace Feel;

    /**
     * Initialize Feel++ Environment
     */
    Environment env(_argc = argc, _argv = argv,
                    _about = about(_name = "assembly",
                                   _author = "Feel++ Consortium",
                                   _email = "feelpp-devel@feelpp.org"));
    Application app;
    app.add(new Assembly<2, 2, Scalar>());
    app.add(new Assembly<2, 2, Vectorial>());
    app.add(new Assembly<2, 3, Vectorial, 2, Scalar>());
    app.run();
}
