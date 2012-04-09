/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-11

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file levelset.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-11
 */

#ifndef _LEVELSET_HPP_
#define _LEVELSET_HPP_

#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/advreact.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/backend_adaptive_reuse_pc.hpp>

#include "reinit_fms.hpp"
#include "reinit_ilp.hpp"
#include "indicator.hpp"

namespace Feel
{
class LevelSet
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type Dim = 3; // was template parameter

    static const uint16_type pOrder = 1;

    typedef double value_type;

    /* entity */
#define ENTITY Simplex

    /* mesh */
    typedef Mesh<GeoEntity<ENTITY<Dim, 1, Dim> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /* bases */
    typedef fusion::vector<fem::Lagrange<Dim, pOrder,
            Scalar, Continuous,
            double, ENTITY> >
            basis_p_type;
    typedef fusion::vector<fem::Lagrange<Dim, 0,
            Scalar, Discontinuous,
            double, ENTITY> >
            basis_i_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_p_type, value_type> space_p_type;
    typedef FunctionSpace<mesh_type, basis_i_type, value_type> space_i_type;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef space_p_type::element_type element_p_type;
    typedef space_i_type::element_type element_i_type;

    /* backends */
    typedef Backend<value_type> backendS_type;
    typedef boost::shared_ptr<backendS_type> backendS_ptrtype;
    //typedef BackendAdaptiveReusePC<backendS_type> backend_type;
    typedef backendS_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /* quadrature for postprocessing */
    static const uint16_type imOrder = 2*pOrder;
    typedef IM<Dim, imOrder, value_type, ENTITY> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<mesh_type>::timeset_type timeset_type;

    LevelSet( int argc, char** argv, AboutData const& ad );

    LevelSet( int argc,
              char** argv,
              AboutData const& ad,
              po::options_description const& od );

    LevelSet( LevelSet const& tc );


    /**
     * alias for run()
     */
    void operator()();

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( int iter,
                        double time,
                        element_p_type& phi,
                        element_p_type& psi,
                        element_i_type& kappa,
                        element_p_type& vx,
                        element_p_type& vy );

    double mass( const element_p_type& psi );

    void advReactUpdateCN( AdvReact<space_p_type, imOrder, ENTITY>& advReact,
                           double dt,
                           double theta,
                           double sign,
                           const element_p_type& vx,
                           const element_p_type& vy,
                           const element_p_type& phi,
                           bool updateStabilization );

    void advReactUpdateBdf2( AdvReact<space_p_type, imOrder, ENTITY>& advReact,
                             double dt,
                             double sign,
                             const element_p_type& vx,
                             const element_p_type& vy,
                             const element_p_type& phi,
                             const element_p_type& phio,
                             bool updateStabilization );

    void statsAfterReinit( const element_p_type& psi,
                           const element_p_type& phi,
                           double massBefore,
                           double mass0 );

private:

    double M_meshSize;
    double M_domainSize;

    boost::shared_ptr<export_type> M_exporter;
    export_type::timeset_ptrtype M_timeSet;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;

    im_type M_im;
}; // LevelSet

} // Feel

#endif /* _LEVELSET_HPP_ */
