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


    // ---------------- P1 Vec space -----------------
    typedef bases<Lagrange<1, Vectorial> > basisP1_Vec_type;
    typedef FunctionSpace<mesh_type, basisP1_Vec_type> spaceP1_Vec_type;
    typedef boost::shared_ptr<spaceP1_Vec_type> spaceP1_Vec_ptrtype;
    typedef typename spaceP1_Vec_type::element_type elementP1_Vec_type;


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

    typedef OperatorInterpolation<space_Vec_type, //espace depart
                                  spaceP1_Vec_type, //espace arrivee
                                  range_visu_ho_type> op_inte_N_to_P1_Vec_type;


    typedef boost::shared_ptr<op_inte_N_to_P1_Vec_type> op_inte_N_to_P1_Vec_ptrtype;

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

    bool prepare = boption("benchmark.prepare");
    if ( prepare )
        nparts = ioption("benchmark.partitions");

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
    M_backend = backend_type::build( soption("backend") );

    exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );

    boost::mpi::timer t;

    double shear = this->vm()["shear"].template as<value_type>();
    bool recombine = boption("recombine");

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
                                              _h= meshSizeInit()/(std::pow(2, level())) ,
                                              _shear=shear,
                                              _xmin=-1.,_xmax=1.,
                                              _ymin=-1.,_ymax=1. ),
                           //                                _refine=level(),
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
    auto backend_l2          = backend_type::build( soption("backend"), "projections" );
    auto backend_l2Vec       = backend_type::build( soption("backend"), "projections" );
    auto backend_l2Smooth    = backend_type::build( soption("backend"), "projections" );
    auto backend_l2SmoothVec = backend_type::build( soption("backend"), "projections" );

    // projectors


    //    double diffnum = (meshSizeInit() / std::pow(2,level())) / (BasisU::nOrder*BasisU::nOrder*BasisU::nOrder);
    //    double diffnum = (meshSizeInit() / std::pow(2,level())) / (BasisU::nOrder*BasisU::nOrder); // with this one I can get order 1  for P1 and P2

    //    double diffnum = (meshSizeInit() / std::pow(2,level())) / (BasisU::nOrder+1); // 0 !
    double diffnum = (meshSizeInit() / std::pow(2,level()-1)) / (BasisU::nOrder*2); // order 1 for P1, P2 and P3

    //    double diffnum = (meshSizeInit() / std::pow(2,level())) / BasisU::nOrder; // order 0

        //std::pow((meshSizeInit() / std::pow(2,level())) , BasisU::nOrder) ;
    // std::pow(meshSizeInit(), std::max(1,BasisU::nOrder-1)) ;

    auto l2p = opProjection(Xh , Xh, _type=L2);
    auto l2pVec = opProjection(Xh_Vec, Xh_Vec, _type=L2);
    auto smooth = projector(Xh , Xh, backend_l2Smooth, DIFF, diffnum, 20);
    auto smoothVec = projector(Xh_Vec, Xh_Vec, backend_l2SmoothVec,DIFF, diffnum, 20);
    auto cip = opProjection(Xh , Xh, _type=CIP);
    auto cipVec = opProjection(Xh_Vec , Xh_Vec, _type=CIP);

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

    double max_modgradphi=0.9;

    auto modgradphi = sqrt( gradv(init_shape) * trans(gradv(init_shape)));

    /* +++++++++++++++ compute quantities with different projections +++++++++++ */

    // auto modgradphituned = sqrt( gradv(init_shape) * trans(gradv(init_shape))) * chi( abs(idv(init_shape)) <= 1.5*h())
    //     + chi( abs(idv(init_shape)) > 1.5*h() );

    /* ------------------ nodal projection ---------------- */
    auto n_nod = vf::project(Xh_Vec, elements(mesh),
                             trans(gradv(init_shape)) / modgradphi );
    auto k_nod = vf::project(Xh, elements(mesh),
                             divv(n_nod) );
    auto nk_nod = vf::project(Xh_Vec, elements(mesh),
                             trans(gradv(k_nod)) /
                             sqrt( gradv(k_nod) * trans(gradv(k_nod))));
    auto kk_nod = vf::project(Xh, elements(mesh),
                             divv(nk_nod) );

    auto phi_nod = init_shape;


    /* ------------------ smooth projection ---------------- */
    auto n_smooth = smoothVec->project( trans(gradv(init_shape)) / modgradphi );

    auto k_smooth = smooth->project( divv(n_smooth) );
    auto phi_smooth = smooth->project( shape_expr );

    auto nk_smooth = smoothVec->project( trans(gradv(k_smooth)) /
                                        vf::max(sqrt( gradv(k_smooth) * trans(gradv(k_smooth))), max_modgradphi) );
    auto kk_smooth = smooth->project( divv(nk_smooth) );


    /* ------------------ L2 projection ---------------- */
    auto n_l2 = l2pVec->project( trans(gradv(init_shape)) / modgradphi );
    auto k_l2 = l2p->project( divv(n_l2) );
    auto nk_l2 = l2pVec->project( trans(gradv(k_l2)) /
                                 vf::max(sqrt( gradv(k_l2) * trans(gradv(k_l2))), max_modgradphi) );
    auto kk_l2 = l2p->project( divv(nk_l2) );

    auto phi_l2 = l2p->project( shape_expr );



    /* ------------------ from hessian (ok if |grad(init_shape)| = 1 ) ---------------- */
    auto k_hess = vf::project(Xh, elements(mesh), trace( hessv(init_shape) ) );
    auto kk_hess = vf::project(Xh, elements(mesh), trace( hessv(k_hess) ) );


    /*  ------------------ integrate by part ----------------------------------   */
    auto n_int = l2pVec->derivate( idv(init_shape) / modgradphi );
    auto k_int = l2p->derivate( trans(idv(n_int)) );
    auto nk_int = l2pVec->derivate( idv(k_int) );
    auto kk_int = l2p->derivate( trans(idv((nk_int))) );


    /* ------------------ optimal projection ---------------- */
    /*           try to get the better combinaison            */
    auto n_opt = n_nod;
    auto k_opt = l2p->derivate( trans(idv(n_opt)) );
    auto nk_opt = vf::project( Xh_Vec, elements(mesh),
                               trans( gradv( k_opt ) )
                               / sqrt( gradv(k_opt) * trans(gradv(k_opt)) ) );

    auto kk_opt = l2p->derivate( trans( idv(nk_opt) ) );

    auto Radius_expr = sqrt( Px() * Px() + Py() * Py() );


    /* ------------------ projection L2 with CIP stabilization ---------------- */
    auto n_cip = cipVec->project( trans(gradv(init_shape)) / modgradphi );
    auto k_cip = cip->project( divv(n_cip) );
    auto nk_cip = cipVec->project( trans(gradv(k_cip) )/
                                 vf::max(sqrt( gradv(k_cip) * trans(gradv(k_cip))), max_modgradphi) );
    auto kk_cip = cip->project( divv(nk_cip) );


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
                               (idv(k_nod) -  1 / Radius_expr) * (idv(k_nod) -  1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_nod = std::sqrt(error_nod);
    M_stats.put( "e.nod.k", error_nod);
    LOG(INFO) << "e.nod.k = " << error_nod << "\n";

    double error_nod_kk = integrate(marked2elements(mesh, 1.),
                               (idv(kk_nod) -  1 / Radius_expr) * (idv(kk_nod) -  1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_nod_kk = std::sqrt(error_nod_kk);
    if (error_nod_kk != error_nod_kk) // for P1 error_nod_kk = Nan and put crashes (get<double>)
        M_stats.put( "e.nod.kk", 1.);
    else
        M_stats.put( "e.nod.kk", error_nod_kk);
    LOG(INFO) << "e.nod.kk = " << error_nod_kk << "\n";

    double error_nod_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_nod) - shape_expr)*(idv(phi_nod) - shape_expr) ).evaluate()(0,0);
    error_nod_phi = std::sqrt(error_nod_phi);
    M_stats.put( "e.nod.phi", error_nod_phi);


    /* l2 */
    double error_l2 = integrate(marked2elements(mesh, 1.),
                                (idv(k_l2) -  1 / Radius_expr) * (idv(k_l2) -  1 / Radius_expr) * Delta).evaluate()(0,0) / perimeter ;
    error_l2 = std::sqrt(error_l2);
    M_stats.put( "e.l2.k", error_l2);
    LOG(INFO) << "e.l2.k = " << error_l2 << "\n";

    double error_l2_kk = integrate(marked2elements(mesh, 1.),
                                (idv(kk_l2) -  1 / Radius_expr) * (idv(kk_l2) -  1 / Radius_expr) * Delta).evaluate()(0,0) / perimeter ;
    error_l2_kk = std::sqrt(error_l2_kk);
    M_stats.put( "e.l2.kk", error_l2_kk);
    LOG(INFO) << "e.l2.kk = " << error_l2_kk << "\n";

    double error_l2_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_l2) - shape_expr)*(idv(phi_l2) - shape_expr) ).evaluate()(0,0);
    error_l2_phi = std::sqrt(error_l2_phi);
    M_stats.put( "e.l2.phi", error_l2_phi);


    double error_l2_int = integrate(marked2elements(mesh, 1.),
                                    (idv(k_int) -  1 / Radius_expr) * (idv(k_int) -  1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_l2_int = std::sqrt( error_l2_int );
    M_stats.put( "e.int.k", error_l2_int);

    double error_l2_int_kk = integrate(marked2elements(mesh, 1.),
                                    (idv(kk_int) -  1 / Radius_expr) * (idv(kk_int) -  1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_l2_int_kk = std::sqrt( error_l2_int_kk );
    M_stats.put( "e.int.kk", error_l2_int_kk);


    /* smooth */
    double error_smooth = integrate(marked2elements(mesh, 1.),
                                    (idv(k_smooth) - 1 / Radius_expr) * (idv(k_smooth) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_smooth = std::sqrt(error_smooth);
    M_stats.put( "e.sm.k", error_smooth);
    LOG(INFO) << "e.sm.k = " << error_smooth << "\n";

    double error_smooth_kk = integrate(marked2elements(mesh, 1.),
                                    (idv(kk_smooth) - 1 / Radius_expr) * (idv(kk_smooth) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_smooth_kk = std::sqrt(error_smooth_kk);
    M_stats.put( "e.sm.kk", error_smooth_kk);
    LOG(INFO) << "e.sm.kk = " << error_smooth_kk << "\n";

    double error_smooth_phi = integrate(marked2elements(mesh, 1.),
                                     (idv(phi_smooth) - shape_expr)*(idv(phi_smooth) - shape_expr) ).evaluate()(0,0);
    error_smooth_phi = std::sqrt(error_smooth_phi);
    M_stats.put( "e.smooth.phi", error_smooth_phi);


    /* hess */
    double error_hess = integrate(marked2elements(mesh, 1.),
                             (idv(k_hess) - 1 / Radius_expr) * (idv(k_hess) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_hess = std::sqrt(error_hess);
    M_stats.put( "e.hs.kproj", error_hess);
    LOG(INFO) << "e.hs.kproj = " << error_hess << "\n";

    double error_hess_kk = integrate(marked2elements(mesh, 1.),
                             (idv(kk_hess) - 1 / Radius_expr) * (idv(kk_hess) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_hess_kk = std::sqrt(error_hess_kk);
    M_stats.put( "e.hs.kkproj", error_hess_kk);
    LOG(INFO) << "e.hs.kkproj = " << error_hess_kk << "\n";

    auto k_hess_proj = trace( hessv(init_shape) );
    double error_hessp = integrate(marked2elements(mesh, 1.),
                                   (k_hess_proj - 1 / Radius_expr) * (k_hess_proj - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_hessp = std::sqrt(error_hessp);
    M_stats.put( "e.hs.k", error_hessp);
    LOG(INFO) << "e.hs.k = " << error_hessp << "\n";

    /* opt */
    double error_opt = integrate(marked2elements(mesh, 1.),
                             (idv(k_opt) - 1 / Radius_expr) * (idv(k_opt) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_opt = std::sqrt(error_opt);
    M_stats.put( "e.opt.k", error_opt);
    LOG(INFO) << "e.opt.k = " << error_opt << "\n";

    double error_opt_kk = integrate(marked2elements(mesh, 1.),
                             (idv(kk_opt) - 1 / Radius_expr) * (idv(kk_opt) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_opt_kk = std::sqrt(error_opt_kk);
    M_stats.put( "e.opt.kk", error_opt_kk);
    LOG(INFO) << "e.opt.kk = " << error_opt_kk << "\n";


    /* CIP */
    double error_cip = integrate(marked2elements(mesh, 1.),
                                    (idv(k_cip) - 1 / Radius_expr) * (idv(k_cip) - 1 / Radius_expr) * Delta ).evaluate()(0,0) / perimeter ;
    error_cip = std::sqrt(error_cip);
    M_stats.put( "e.ci.k", error_cip);
    LOG(INFO) << "e.ci.k = " << error_cip << "\n";

    std::cout<<"exporting ...\n";

    if ( exporter->doExport() )
    {

        if (BasisU::nOrder>1)
            {
                // high order visu
                op_inte_N_to_P1_ptrtype op_inte_N_to_P1;
                OperatorLagrangeP1<space_type> * opLagP1;

                op_inte_N_to_P1_Vec_ptrtype op_inte_N_to_P1_Vec;
                OperatorLagrangeP1<space_Vec_type> * opLagP1Vec;

                auto backend_oplag = backend_type::build( soption("backend") );
                auto backend_oplagVec = backend_type::build( soption("backend") );
                opLagP1 = new OperatorLagrangeP1<space_type> (Xh, backend_oplag);
                opLagP1Vec = new OperatorLagrangeP1<space_Vec_type> (Xh_Vec, backend_oplagVec);

                auto Xh_P1 = spaceP1_type::New(_mesh=opLagP1->mesh() );
                op_inte_N_to_P1 = opInterpolation(_domainSpace = Xh, _imageSpace = Xh_P1);

                auto Xh_P1_Vec = spaceP1_Vec_type::New(_mesh=opLagP1Vec->mesh() );
                op_inte_N_to_P1_Vec = opInterpolation(_domainSpace = Xh_Vec, _imageSpace = Xh_P1_Vec);


                auto delta_p1 = Xh_P1->element();
                auto k_l2_p1 = Xh_P1->element();
                auto k_smooth_p1 = Xh_P1->element();
                auto k_nod_p1 = Xh_P1->element();
                auto k_hess_p1 = Xh_P1->element();
                auto k_int_p1 = Xh_P1->element();
                auto k_opt_p1 = Xh_P1->element();
                auto k_cip_p1 = Xh_P1->element();

                auto nk_l2_p1 = Xh_P1_Vec->element();
                auto nk_nod_p1 = Xh_P1_Vec->element();
                auto nk_smooth_p1 = Xh_P1_Vec->element();
                auto nk_int_p1 = Xh_P1_Vec->element();
                auto nk_opt_p1 = Xh_P1_Vec->element();
                auto nk_cip_p1 = Xh_P1_Vec->element();


                op_inte_N_to_P1->apply(Delta_proj, delta_p1);
                op_inte_N_to_P1->apply(k_l2, k_l2_p1);
                op_inte_N_to_P1->apply(k_smooth, k_smooth_p1);
                op_inte_N_to_P1->apply(k_nod, k_nod_p1);
                op_inte_N_to_P1->apply(k_hess, k_hess_p1);
                op_inte_N_to_P1->apply(k_int, k_int_p1);
                op_inte_N_to_P1->apply(k_opt, k_opt_p1);
                op_inte_N_to_P1->apply(k_cip, k_cip_p1);

                op_inte_N_to_P1_Vec->apply( nk_l2, nk_l2_p1);
                op_inte_N_to_P1_Vec->apply( nk_nod, nk_nod_p1);
                op_inte_N_to_P1_Vec->apply( nk_smooth, nk_smooth_p1);
                op_inte_N_to_P1_Vec->apply( nk_int, nk_int_p1);
                op_inte_N_to_P1_Vec->apply( nk_opt, nk_opt_p1);
                op_inte_N_to_P1_Vec->apply( nk_cip, nk_cip_p1);

                exporter->step( 0 )->setMesh( opLagP1->mesh() );
                exporter->step( 0 )->add( "Delta", delta_p1 );
                exporter->step( 0 )->add("k_l2", k_l2_p1);
                exporter->step( 0 )->add("k_smooth",k_smooth_p1);
                exporter->step( 0 )->add("k_nod", k_nod_p1);
                exporter->step( 0 )->add("k_hess", k_hess_p1);
                exporter->step( 0 )->add("k_int", k_int_p1);
                exporter->step( 0 )->add("k_opt", k_opt_p1);
                exporter->step( 0 )->add("k_cip", k_cip_p1);

                exporter->step( 0 )->add("nk_l2", nk_l2_p1);
                exporter->step( 0 )->add("nk_nod", nk_nod_p1);
                exporter->step( 0 )->add("nk_smooth", nk_smooth_p1);
                exporter->step( 0 )->add("nk_int", nk_int_p1);
                exporter->step( 0 )->add("nk_opt", nk_opt_p1);
                exporter->step( 0 )->add("nk_cip", nk_cip_p1);

                exporter->step( 0 )->add("marker_delta", marker_delta);
                exporter->step( 0 )->add("n_l2", n_l2);

                exporter->step(0)->add("modgradphi", l2p->project( modgradphi ));

                exporter->step( 0 )->add("init_shape", init_shape);

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
                exporter->step( 0 )->add("k_cip", k_cip);
                exporter->step( 0 )->add("marker_delta", marker_delta);

                exporter->step( 0 )->add("n_l2", n_l2);
                exporter->step( 0 )->add("n_smooth", n_smooth);

                exporter->step( 0 )->add("init_shape", init_shape);
                exporter->save();
            }
    }

} // Curvature::run


} // Feel

#endif // __FEELPP_BENCH_CURVATURE_HPP

