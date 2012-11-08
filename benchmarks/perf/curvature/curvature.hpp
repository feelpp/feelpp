/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __FEELPP_BENCH_CURVATURE_HPP
#define __FEELPP_BENCH_CURVATURE_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/feel.hpp>

#define PROJ_INT_BY_PART 1

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

    // ---------------- P0 space -----------------
    typedef bases<Lagrange<0, Scalar, Discontinuous> > basisP0_type;
    typedef FunctionSpace<mesh_type, basisP0_type> spaceP0_type;
    typedef boost::shared_ptr<spaceP0_type> spaceP0_ptrtype;
    typedef typename spaceP0_type::element_type elementP0_type;

    // ---------------- P1 space -----------------
    typedef bases<Lagrange<1> > basisP1_type;
    typedef FunctionSpace<mesh_type, basisP1_type> spaceP1_type;
    typedef boost::shared_ptr<spaceP1_type> spaceP1_ptrtype;
    typedef typename spaceP1_type::element_type elementP1_type;


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


    /* spaces and interpolators for high order visu  */
    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<mesh_type>::element_const_iterator,
                         typename MeshTraits<mesh_type>::element_const_iterator> range_visu_ho_type;

    // ----------------------- operator interpolation -------------
    typedef OperatorInterpolation<space_type, //espace depart
                                  spaceP1_type, //espace arrivee
                                  range_visu_ho_type> op_inte_N_to_P1_type;

    typedef boost::shared_ptr<op_inte_N_to_P1_type> op_inte_N_to_P1_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Curvature( std::string const& basis_name )
        :
        super(),
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

private:

    backend_ptrtype M_backend;
    std::string M_basis_name;

    boost::shared_ptr<export_type> exporter;

    // bench parameters
    double x0, y0, Radius;
    element_type init_shape;
    space_ptrtype Xh;
    mesh_ptrtype mesh;

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
    mesh = createGMSHMesh( _mesh=new mesh_type,
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

    M_stats.put( "t.init.mesh",t.elapsed() );t.restart();

    size_type gnelts=0;
    mpi::all_reduce( this->comm(), mesh->numElements() , gnelts, [] ( size_type x, size_type y ) {return x + y;} );

    // space
    Xh = space_type::New( mesh );
    M_stats.put( "n.space.nelts", gnelts );
    M_stats.put( "n.space.nlocalelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nDof() );
    M_stats.put( "n.space.nlocaldof",Xh->nLocalDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";

    space_Vec_ptrtype Xh_Vec = space_Vec_type::New( mesh );

    M_stats.put( "n.spacev.nelts", gnelts );
    M_stats.put( "n.spacev.nlocalelts",Xh_Vec->mesh()->numElements() );
    M_stats.put( "n.spacev.ndof",Xh_Vec->nDof() );
    M_stats.put( "n.spacev.nlocaldof",Xh_Vec->nLocalDof() );
    M_stats.put( "t.init.spacev",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    spaceP0_ptrtype Xh_P0 = spaceP0_type::New( mesh );
    M_stats.put( "n.spacep0.nelts", gnelts );
    M_stats.put( "n.spacep0.nlocalelts",Xh_P0->mesh()->numElements() );
    M_stats.put( "n.spacep0.ndof",Xh_P0->nDof() );
    M_stats.put( "n.spacep0.nlocaldof",Xh_P0->nLocalDof() );
    M_stats.put( "t.init.spacep0",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    // backends
    auto backend_l2 = backend_type::build( this->vm(), "projections" );
    auto backend_l2Vec = backend_type::build( this->vm(),  "projections" );
    auto backend_l2Smooth = backend_type::build( this->vm(),  "projections" );
    auto backend_l2SmoothVec = backend_type::build( this->vm(),  "projections" );

    // projectors

    double diffnum = std::pow(meshSizeInit(), std::max(1,BasisU::nOrder-1)) * 0.001 ;

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
                auto shape_expr = (vf::sqrt( X*X + Y*Y ) - Radius);

                // exact signed distance function to a circle
                // might be better with a l2 projection ???
                init_shape = vf::project(Xh, elements(mesh), shape_expr );

                break;
            }//circle
        } //switch

    auto X = Px() - x0;
    auto Y = Py() - y0;

    auto shape_expr = (vf::sqrt( X*X + Y*Y ) - Radius);

    // "delta like" function
    double order = std::max(1,BasisU::nOrder-1);
    // auto eps = 1.5 * vf::pow(vf::h(),order);
    // auto Delta =  ( vf::chi( shape_expr<-eps )*vf::constant(0.0)
    //                 +
    //                 vf::chi( shape_expr>=-eps )*vf::chi( shape_expr<=eps )*
    //                 1/(2*eps) *( 1 + cos(pi*shape_expr/eps) )
    //                 +
    //                 vf::chi(shape_expr>eps)*vf::constant(0.0) );
    auto eps = 1.5 * h();
    auto Delta =  ( vf::chi( shape_expr<-eps )*vf::constant(0.0)
                    +
                    vf::chi( shape_expr>=-eps )*vf::chi( shape_expr<=eps )*
                    1/(2*eps) *( 1 + cos(pi*shape_expr/eps) )
                    +
                    vf::chi(shape_expr>eps)*vf::constant(0.0) );

    LOG(INFO) << "project Delta...\n";
    auto Delta_proj = vf::project(Xh, elements(mesh),Delta );//l2p->project(_expr=Delta);

    LOG(INFO) << "project Delta done...\n";



    auto marker_delta = vf::project(Xh_P0, elements(mesh),
                                    vf::chi( idv(Delta_proj)>0.0001) * vf::cst(1.) );
    mesh->updateMarker2(marker_delta);


    /* +++++++++++++++ compute quantities with different projections +++++++++++ */

    /* ------------------ nodal projection ---------------- */
    auto n_nod = vf::project(Xh_Vec, elements(mesh),
                             trans(gradv(init_shape)) /
                             sqrt( gradv(init_shape) * trans(gradv(init_shape))));
    auto k_nod = vf::project(Xh, elements(mesh),
                             divv(n_nod) );

    auto nk_nod = vf::project(Xh_Vec, elements(mesh),
                             trans(gradv(k_nod)) /
                             sqrt( gradv(k_nod) * trans(gradv(k_nod))));
    auto kk_nod = vf::project(Xh, elements(mesh),
                             divv(nk_nod) );


    auto phi_nod = init_shape;


    double max_modgradphi=0.01;

    /* ------------------ L2 projection ---------------- */
    auto n_l2 = l2pVec->project( gradv(init_shape) /
                                 vf::max(sqrt( gradv(init_shape) * trans(gradv(init_shape))), max_modgradphi) );
    auto k_l2 = l2p->project( divv(n_l2) );
    auto nk_l2 = l2pVec->project( gradv(k_l2) /
                                 vf::max(sqrt( gradv(k_l2) * trans(gradv(k_l2))), max_modgradphi) );
    auto kk_l2 = l2p->project( divv(nk_l2) );

    auto phi_l2 = l2p->project( shape_expr );


    /* L2 without using projector directly */
    auto backend_int = backend_type::build( this->vm() );
    auto n_l2_bis = Xh_Vec->element();

    auto D = backend_int->newMatrix(Xh_Vec, Xh_Vec);
    auto F = backend_int->newVector(Xh_Vec);

    form2(Xh_Vec, Xh_Vec, D) = integrate(elements(mesh), trans(idt(n_l2_bis)) * id(n_l2_bis) );

    form1(Xh_Vec, F) = integrate(elements(mesh), gradv(init_shape) * trans(grad(init_shape))
                                 / vf::max(sqrt( gradv(init_shape) * trans(gradv(init_shape))), max_modgradphi)  );

    D->close();
    F->close();

    backend_int->solve(D, n_l2_bis, F);

    backend_int = backend_type::build( this->vm() );
    auto k_l2_bis = Xh->element();

    D = backend_int->newMatrix(Xh, Xh);
    F = backend_int->newVector(Xh);

    form2(Xh, Xh, D) = integrate(elements(mesh), idt(k_l2_bis) * id(k_l2_bis) );
    form1(Xh, F) = integrate(elements(mesh), divv(n_l2_bis) * id(k_l2_bis) );

    D->close();
    F->close();

    backend_int->solve(D, k_l2_bis, F);


    /* ------------------ smooth projection ---------------- */
    auto n_smooth = smoothVec->project( gradv(init_shape) /
                                        vf::max(sqrt( gradv(init_shape) * trans(gradv(init_shape))), max_modgradphi) );
    auto k_smooth = smooth->project( divv(n_smooth) );
    auto phi_smooth = smooth->project( shape_expr );

    auto nk_smooth = smoothVec->project( gradv(k_smooth) /
                                        vf::max(sqrt( gradv(k_smooth) * trans(gradv(k_smooth))), max_modgradphi) );
    auto kk_smooth = smooth->project( divv(nk_smooth) );


    /* ------------------ from hessian (ok if |grad(init_shape)| = 1 ) ---------------- */
    auto k_hess = vf::project(Xh, elements(mesh), trace( hessv(init_shape) ) );
    auto kk_hess = vf::project(Xh, elements(mesh), trace( hessv(k_hess) ) );


    #if PROJ_INT_BY_PART
    /* instead of solving : int_omega k * v = int_omega div( grad(phi) / |grad phi|)
       solve : int k * v = - int_omega grad(phi) / |grad phi| . grad(v) + int_Gamma grad(phi) / |grad phi| . N * v
    */
    backend_int = backend_type::build( this->vm() );
    auto k_int = Xh->element();
    auto kk_int = Xh->element();

    D = backend_int->newMatrix(Xh, Xh);
    F = backend_int->newVector(Xh);

    form2(Xh, Xh, D) = integrate(elements(mesh), idt(k_int) * id(k_int) * vf::max(sqrt( gradv(init_shape) * trans(gradv(init_shape))), max_modgradphi) );
    form1(Xh, F) = integrate(elements(mesh),
                               - gradv(init_shape) * trans(grad(init_shape)) );

    form1(Xh, F) += integrate(boundaryfaces(mesh),
                                 gradv(init_shape) * N() * id(init_shape) );

    D->close();
    F->close();

    backend_int->solve(D, k_int, F);

    form2(Xh, Xh, D) = integrate(elements(mesh), idt(kk_int) * id(kk_int) * vf::max(sqrt( gradv(k_int) * trans(gradv(k_int))), max_modgradphi) );
    form1(Xh, F) = integrate(elements(mesh),
                               - gradv(k_int) * trans(grad(k_int)) );

    form1(Xh, F) += integrate(boundaryfaces(mesh),
                                 gradv(k_int) * N() * id(k_int) );

    D->close();
    F->close();

    backend_int->solve(D, kk_int, F);


    #endif

    // +++++++++++++++++++ error computation ++++++++++++++++++++++
    double perimeter = integrate( _range=marked2elements(mesh, 1.), _expr=Delta ).evaluate()(0,0);
    LOG(INFO) << "perimeter = " << perimeter << "\n";

    double error_perimeter = math::abs(perimeter - 2 * pi * Radius);
    M_stats.put( "e.l2.perim", error_perimeter);
    M_stats.put( "d.value.double.perimeter", perimeter);
    LOG(INFO) << "e.l2.perim = " << error_perimeter << "\n";

    double int_modgradphi = integrate(marked2elements(mesh, 1.), sqrt( gradv(init_shape) * trans(gradv(init_shape))) ).evaluate()(0,0);
    int_modgradphi /= integrate(marked2elements(mesh, 1.), cst(1.)).evaluate()(0,0);

    double error_modgraphi = math::abs(int_modgradphi - 1.);
    M_stats.put( "e.l2.modgradphi", error_modgraphi);
    LOG(INFO) << "e.l2.modgradphi = " << error_modgraphi << "\n";
    M_stats.put( "d.value.double.modgraphi", int_modgradphi);

    /*  nodal  */
    double error_nod = integrate(marked2elements(mesh, 1.),
                               (idv(k_nod) -  1 / Radius) * (idv(k_nod) -  1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_nod = std::sqrt(error_nod);
    M_stats.put( "e.nod.k", error_nod);
    LOG(INFO) << "e.nod.k = " << error_nod << "\n";

    double error_nod_kk = integrate(marked2elements(mesh, 1.),
                               (idv(kk_nod) -  1 / Radius) * (idv(kk_nod) -  1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_nod_kk = std::sqrt(error_nod_kk);
    M_stats.put( "e.nod.kk", error_nod_kk);
    LOG(INFO) << "e.nod.kk = " << error_nod_kk << "\n";

    double error_nod_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_nod) - shape_expr)*(idv(phi_nod) - shape_expr) ).evaluate()(0,0);
    error_nod_phi = std::sqrt(error_nod_phi);
    M_stats.put( "e.nod.phi", error_nod_phi);


    /* l2 */
    double error_l2 = integrate(marked2elements(mesh, 1.),
                                (idv(k_l2) -  1 / Radius) * (idv(k_l2) -  1 / Radius) * Delta).evaluate()(0,0) / perimeter ;
    error_l2 = std::sqrt(error_l2);
    M_stats.put( "e.l2.k", error_l2);
    LOG(INFO) << "e.l2.k = " << error_l2 << "\n";

    double error_l2_kk = integrate(marked2elements(mesh, 1.),
                                (idv(kk_l2) -  1 / Radius) * (idv(kk_l2) -  1 / Radius) * Delta).evaluate()(0,0) / perimeter ;
    error_l2_kk = std::sqrt(error_l2_kk);
    M_stats.put( "e.l2.kk", error_l2_kk);
    LOG(INFO) << "e.l2.kk = " << error_l2_kk << "\n";

    double error_l2_bis = integrate(marked2elements(mesh, 1.),
                                (idv(k_l2_bis) -  1 / Radius) * (idv(k_l2_bis) -  1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_l2_bis = std::sqrt(error_l2_bis);
    M_stats.put( "e.l2.k_bis", error_l2_bis);

    double error_l2_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_l2) - shape_expr)*(idv(phi_l2) - shape_expr) ).evaluate()(0,0);
    error_l2_phi = std::sqrt(error_l2_phi);
    M_stats.put( "e.l2.phi", error_l2_phi);


    #if PROJ_INT_BY_PART
    double error_l2_int = integrate(marked2elements(mesh, 1.),
                                    (idv(k_int) -  1 / Radius) * (idv(k_int) -  1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_l2_int = std::sqrt( error_l2_int );
    M_stats.put( "e.l2.kk_int", error_l2_int);

    double error_l2_int_kk = integrate(marked2elements(mesh, 1.),
                                    (idv(kk_int) -  1 / Radius) * (idv(kk_int) -  1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_l2_int_kk = std::sqrt( error_l2_int_kk );
    M_stats.put( "e.l2.kk_int", error_l2_int_kk);

    #endif

    /* smooth */
    double error_smooth = integrate(marked2elements(mesh, 1.),
                                    (idv(k_smooth) - 1 / Radius) * (idv(k_smooth) - 1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_smooth = std::sqrt(error_smooth);
    M_stats.put( "e.sm.k", error_smooth);
    LOG(INFO) << "e.sm.k = " << error_smooth << "\n";

    double error_smooth_kk = integrate(marked2elements(mesh, 1.),
                                    (idv(kk_smooth) - 1 / Radius) * (idv(kk_smooth) - 1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_smooth_kk = std::sqrt(error_smooth_kk);
    M_stats.put( "e.sm.kk", error_smooth_kk);
    LOG(INFO) << "e.sm.kk = " << error_smooth_kk << "\n";

    double error_smooth_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_smooth) - shape_expr)*(idv(phi_smooth) - shape_expr) ).evaluate()(0,0);
    error_smooth_phi = std::sqrt(error_smooth_phi);
    M_stats.put( "e.smooth.phi", error_smooth_phi);


    /* hess */
    double error_hess = integrate(marked2elements(mesh, 1.),
                             (idv(k_hess) - 1 / Radius) * (idv(k_hess) - 1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_hess = std::sqrt(error_hess);
    M_stats.put( "e.hs.kproj", error_hess);
    LOG(INFO) << "e.hs.kproj = " << error_hess << "\n";

    double error_hess_kk = integrate(marked2elements(mesh, 1.),
                             (idv(kk_hess) - 1 / Radius) * (idv(kk_hess) - 1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_hess_kk = std::sqrt(error_hess_kk);
    M_stats.put( "e.hs.kkproj", error_hess_kk);
    LOG(INFO) << "e.hs.kkproj = " << error_hess_kk << "\n";

    auto k_hess_proj = trace( hessv(init_shape) );
    double error_hessp = integrate(marked2elements(mesh, 1.),
                                   (k_hess_proj - 1 / Radius) * (k_hess_proj - 1 / Radius) * Delta ).evaluate()(0,0) / perimeter ;
    error_hessp = std::sqrt(error_hessp);
    M_stats.put( "e.hs.k", error_hessp);
    LOG(INFO) << "e.hs.k = " << error_hessp << "\n";

    std::cout<<"exporting ...\n";

    if ( exporter->doExport() )
    {

        if (BasisU::nOrder>1)
            {
                // high order visu
                op_inte_N_to_P1_ptrtype op_inte_N_to_P1;
                OperatorLagrangeP1<space_type> * opLagP1;

                auto backend_oplag = backend_type::build( this->vm() );
                opLagP1 = new OperatorLagrangeP1<space_type> (Xh, backend_oplag);
                auto Xh_P1 = spaceP1_type::New(_mesh=opLagP1->mesh() );
                op_inte_N_to_P1 = opInterpolation(_domainSpace = Xh, _imageSpace = Xh_P1);

                auto delta_p1 = Xh_P1->element();
                auto k_l2_p1 = Xh_P1->element();
                auto k_smooth_p1 = Xh_P1->element();
                auto k_nod_p1 = Xh_P1->element();
                auto k_hess_p1 = Xh_P1->element();
                auto k_int_p1 = Xh_P1->element();

                auto kk_nod_p1 = Xh_P1->element();
                auto kk_l2_p1  = Xh_P1->element();
                auto kk_smooth_p1 = Xh_P1->element();
                auto kk_int_p1 = Xh_P1->element();
                auto kk_hess_p1 = Xh_P1->element();

                op_inte_N_to_P1->apply(Delta_proj, delta_p1);
                op_inte_N_to_P1->apply(k_l2, k_l2_p1);
                op_inte_N_to_P1->apply(k_smooth, k_smooth_p1);
                op_inte_N_to_P1->apply(k_nod, k_nod_p1);
                op_inte_N_to_P1->apply(k_hess, k_hess_p1);
                op_inte_N_to_P1->apply(k_int, k_int_p1);

                op_inte_N_to_P1->apply(kk_l2, kk_l2_p1);
                op_inte_N_to_P1->apply(kk_smooth, kk_smooth_p1);
                op_inte_N_to_P1->apply(kk_nod, kk_nod_p1);
                op_inte_N_to_P1->apply(kk_hess, kk_hess_p1);
                op_inte_N_to_P1->apply(kk_int, kk_int_p1);

                exporter->step( 0 )->setMesh( opLagP1->mesh() );
                exporter->step( 0 )->add( "Delta", delta_p1 );
                exporter->step( 0 )->add("k_l2", k_l2_p1);
                exporter->step( 0 )->add("k_smooth",k_smooth_p1);
                exporter->step( 0 )->add("k_nod", k_nod_p1);
                exporter->step( 0 )->add("k_hess", k_hess_p1);
                exporter->step( 0 )->add("k_int", k_int_p1);


                exporter->step( 0 )->add("kk_nod", kk_nod_p1);
                exporter->step( 0 )->add("kk_l2", kk_l2_p1);
                exporter->step( 0 )->add("kk_smooth", kk_smooth_p1);
                exporter->step( 0 )->add("kk_int", kk_int_p1);
                exporter->step( 0 )->add("kk_hess", kk_hess_p1);


                exporter->step( 0 )->add("marker_delta", marker_delta);
                exporter->step( 0 )->add("n_l2", n_l2);

                exporter->save();
            }
        else
            {
                exporter->step( 0 )->setMesh( mesh );
                exporter->step( 0 )->add( "Delta", Delta_proj );
                exporter->step( 0 )->add("k_l2", k_l2);
                exporter->step( 0 )->add("k_smooth", k_smooth);
                exporter->step( 0 )->add("k_nod", k_nod);
                exporter->step( 0 )->add("k_hess", k_hess);
                exporter->step( 0 )->add("k_int", k_int);
                exporter->step( 0 )->add("marker_delta", marker_delta);

                exporter->step( 0 )->add("kk_l2", kk_l2);
                exporter->step( 0 )->add("kk_smooth", kk_smooth);
                exporter->step( 0 )->add("kk_int", kk_int);
                exporter->step( 0 )->add("kk_hess", kk_hess);

                exporter->step( 0 )->add("n_l2", n_l2);

                exporter->save();
            }

    }

} // Curvature::run


} // Feel

#endif // __FEELPP_BENCH_CURVATURE_HPP

