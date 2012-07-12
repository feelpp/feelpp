/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-07-11

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file curvature.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#ifndef __FEELPP_BENCH_CURVATURE_HPP
#define __FEELPP_BENCH_CURVATURE_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/simget.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
//#include <feel/feelpoly/crouzeixraviart.hpp>

#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

#define quad 10

namespace Feel
{
/**
 * \class Curvature class
 * \brief compute the curvature of a shape for by different methods
 *
 */
template<int Dim,
         typename BasisU,
         typename BasisU_Vec,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Curvature
    :
public Simget
{
    typedef Simget super;
public:

    enum shape_type {circle=1 /*,ellipse, star*/ };

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #
    typedef BasisU basis_u_type;
    typedef bases<basis_u_type> basis_type;


    typedef BasisU_Vec basis_u_Vec_type;
    typedef bases<basis_u_Vec_type> basis_Vec_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef FunctionSpace<mesh_type, basis_Vec_type> space_Vec_type;
    typedef boost::shared_ptr<space_Vec_type> space_Vec_ptrtype;
    typedef typename space_Vec_type::element_type element_Vec_type;


    /* export */
    typedef Exporter<mesh_type> export_type;

    Curvature( std::string const& basis_name,
               po::variables_map const& vm, AboutData const& ad )
        :
        super( vm, ad ),
        M_backend(),
        M_basis_name( basis_name ),
        exporter()
    {
        // mu = this->vm()["mu"].template as<value_type>();
        // penalbc = this->vm()["bccoeff"].template as<value_type>();
    }


    std::string name() const
    {
        return M_basis_name;
    }

    /**
     * run the convergence test
     */
    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N )
    {
        run();
    }


private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& Delta, element_type& k_l2, element_type& k_smooth, element_type& k_nod, element_type& k_hess);

private:

    backend_ptrtype M_backend;
    std::string M_basis_name;

    boost::shared_ptr<export_type> exporter;

    // bench parameters
    double x0, y0, Radius;
    element_type init_shape;

}; // Curvature


    template<int Dim, typename BasisU, typename BasisU_Vec, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Curvature<Dim, BasisU, BasisU_Vec, Entity>::run()
{
    using namespace Feel::vf;

    static double pi = M_PI;

    int nparts = Environment::worldComm().size();
    bool prepare = this->vm()["benchmark.prepare"].template as<bool>();
    if ( prepare )
        nparts = this->vm()["benchmark.partitions"].template as<int>();

    // give as parameter when other shapes are avaliable
    shape_type shape = circle;

    if ( this->vm().count( "nochdir" ) == false )
    {
        this->changeRepository( boost::format( "perf/%1%/%2%/%3%/h_%4%/l_%5%/parts_%6%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % M_basis_name
                                % meshSizeInit()
                                % level()
                                % nparts );
    }

    //! init backend
    M_backend = backend_type::build( this->vm() );

    exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );

    boost::mpi::timer t;

    double shear = this->vm()["shear"].template as<value_type>();
    bool recombine = this->vm()["recombine"].template as<bool>();

    /*
     * First we create the mesh, in the case of quads we wish to have
     * non-regular meshes to ensure that we don't have some super-convergence
     * associated to the mesh. Hence we use recombine=true to recombine
     * triangles generated by a Delaunay algorithm into quads
     */
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % Dim % 1 ).str() ,
                                              _shape="hypercube",
                                              _usenames=true,
                                              _dim=Dim,
                                              _h=meshSizeInit(),
                                              _shear=shear,
                                              _xmin=-1.,_xmax=1.,
                                              _ymin=-1.,_ymax=1. ),
                                _refine=level(),
                                _partitions=nparts );

    // space
    space_ptrtype Xh = space_type::New( mesh );
    space_Vec_ptrtype Xh_Vec = space_Vec_type::New( mesh );

    // backends
    auto backend_l2 = backend_type::build( this->vm() );
    auto backend_l2Vec = backend_type::build( this->vm() );
    auto backend_l2Smooth = backend_type::build( this->vm() );
    auto backend_l2SmoothVec = backend_type::build( this->vm() );

    // projectors
    double diffnum = meshSizeInit() * 0.01;
    auto l2p = projector(Xh , Xh, backend_l2, L2);
    auto l2pVec = projector(Xh_Vec, Xh_Vec, backend_l2Vec, L2);
    auto smooth = projector(Xh , Xh, backend_l2Smooth, DIFF, diffnum, 20);
    auto smoothVec = projector(Xh_Vec, Xh_Vec, backend_l2SmoothVec, DIFF, diffnum, 20);

    t.restart();
    if ( prepare ) return;

    switch (shape)
        {
        case circle :
            {

                Radius = 0.5;
                x0 = 0;
                y0 = 0;

                auto X = Px() - x0;
                auto Y = Py() - y0;

                // exact signed distance function to a circle
                // might be better with a l2 projection ???
                init_shape = vf::project(Xh, elements(mesh),
                                   sqrt( X*X + Y*Y ) - Radius );
                break;
            }//circle
        } //switch

    // "delta like" function
    auto eps = 1.5 * vf::h();
    auto Delta = l2p->project(vf::chi( idv(init_shape)<-eps )*vf::constant(0.0)
                              +
                              vf::chi( idv(init_shape)>=-eps )*vf::chi( idv(init_shape)<=eps )*
                              1/(2*eps) *( 1 + cos(pi*idv(init_shape)/eps) )
                              +
                              vf::chi(idv(init_shape)>eps)*vf::constant(0.0) );


    /* +++++++++++++++ compute quantities with different projections +++++++++++ */

    /* ------------------ nodal projection ---------------- */
    auto n_nod = vf::project(Xh_Vec, elements(mesh),
                             trans(gradv(init_shape)) /
                             sqrt( gradv(init_shape) * trans(gradv(init_shape))));
    auto k_nod = vf::project(Xh, elements(mesh),
                             divv(n_nod) );


    /* ------------------ L2 projection ---------------- */
    auto n_l2 = l2pVec->project( gradv(init_shape) /
                                 sqrt( gradv(init_shape) * trans(gradv(init_shape))), _quad=_Q<quad>() );
    auto k_l2 = l2p->project( divv(n_l2), _quad=_Q<quad>() );


    /* ------------------ smooth projection ---------------- */
    auto n_smooth = smoothVec->project( gradv(init_shape) /
                                        sqrt( gradv(init_shape) * trans(gradv(init_shape))), _quad=_Q<quad>() );
    auto k_smooth = smooth->project( divv(n_smooth),_quad=_Q<quad>() );


    /* ------------------ from hessian (ok if |grad(init_shape)| = 1 ) ---------------- */
    auto k_hess = vf::project(Xh, elements(mesh), trace( hessv(init_shape) ) );


    // +++++++++++++++++++ error computation ++++++++++++++++++++++
    double perimeter = integrate(elements(mesh), idv(Delta), _Q<quad>()).evaluate()(0,0);
    double error_perimeter = std::sqrt( (perimeter - 2 * pi * Radius)*(perimeter - 2 * pi * Radius) );
    M_stats.put( "e.l2.perim", error_perimeter);

    double int_modgradphi = integrate(elements(mesh), sqrt( gradv(init_shape) * trans(gradv(init_shape))), _quad = _Q<quad>() ).evaluate()(0,0);
    int_modgradphi /= integrate(elements(mesh), cst(1.)).evaluate()(0,0);
    double error_modgraphi = std::sqrt( (int_modgradphi - 1.)*(int_modgradphi - 1.) );
    M_stats.put( "e.l2.modgradphi", error_modgraphi);


    double error_nod = integrate(elements(mesh),
                               (idv(k_nod) -  1 / Radius) * (idv(k_nod) -  1 / Radius) * idv(Delta),
              _Q<quad>() ).evaluate()(0,0) / perimeter ;
    error_nod = std::sqrt(error_nod);
    M_stats.put( "e.nod.k", error_nod);

    double error_l2 = integrate(elements(mesh),
                                (idv(k_l2) -  1 / Radius) * (idv(k_l2) -  1 / Radius) * idv(Delta),
              _Q<quad>() ).evaluate()(0,0) / perimeter ;
    error_l2 = std::sqrt(error_l2);
    M_stats.put( "e.l2.k", error_l2);

    double error_smooth = integrate(elements(mesh),
                             (idv(k_smooth) - 1 / Radius) * (idv(k_smooth) - 1 / Radius) * idv(Delta),
              _Q<quad>() ).evaluate()(0,0) / perimeter ;
    error_smooth = std::sqrt(error_smooth);
    M_stats.put( "e.sm.k", error_smooth);


    double error_hess = integrate(elements(mesh),
                             (idv(k_hess) - 1 / Radius) * (idv(k_hess) - 1 / Radius) * idv(Delta),
              _Q<quad>() ).evaluate()(0,0) / perimeter ;
    error_hess = std::sqrt(error_hess);
    M_stats.put( "e.hs.k", error_hess);

    exportResults(Delta, k_l2, k_smooth, k_nod, k_hess);

} // Curvature::run



template<int Dim, typename BasisU, typename BasisU_Vec, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Curvature<Dim, BasisU, BasisU_Vec, Entity>::exportResults( element_type& Delta, element_type& k_l2, element_type& k_smooth, element_type& k_nod, element_type& k_hess)
{
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( Delta.functionSpace()->mesh() );
        exporter->step( 0 )->add( "Delta", Delta );
        exporter->step( 0 )->add("k_l2", k_l2);
        exporter->step( 0 )->add("k_smooth", k_smooth);
        exporter->step( 0 )->add("k_nod", k_nod);
        exporter->save();
    }
} // Curvature::export
} // Feel

#endif // __FEELPP_BENCH_CURVATURE_HPP

