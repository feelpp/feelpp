/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-01-25

  Copyright (C) 2007 EPFL

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
   \file kovasznay.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-01-25
 */

#ifndef _KOVASZNAY_HPP_
#define _KOVASZNAY_HPP_

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/backend_adaptive_reuse_pc.hpp>
#include <feel/feeldiscr/oseen.hpp>



namespace Feel
{
class Kovasznay
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type uOrder = 2;
    static const uint16_type pOrder = 2;

    static const int Dim = 2; // was a template parameter

    typedef double value_type;

    /* entity */
#define ENTITY Simplex

    /* mesh */
    typedef Mesh<ENTITY<Dim > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /* bases */
    typedef bases<Lagrange<uOrder, Vectorial, Continuous>,
            Lagrange<pOrder, Scalar, Continuous> > basis_type;
    typedef bases<Lagrange<0, Scalar, Discontinuous> > basis_i_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    typedef element_type::sub_element<0>::type::functionspace_type space_U_type;
    typedef boost::shared_ptr<space_U_type> space_U_ptrtype;
    typedef space_U_type::element_type element_U_type;

    typedef element_type::sub_element<1>::type::functionspace_type space_p_type;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef space_p_type::element_type element_p_type;

    typedef FunctionSpace<mesh_type, basis_i_type, value_type> space_i_type;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef space_i_type::element_type element_i_type;

    /* quadrature for postprocessing */
    static const uint16_type imOrder = 3*uOrder-1; // 5 / 20
    typedef IM<Dim, imOrder, value_type, ENTITY> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    // linear algebra backends
    //typedef BackendAdaptiveReusePC<Backend<value_type> > backend_type;
    typedef Backend<value_type > backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    Kovasznay( int argc, char** argv, AboutData const& ad );

    Kovasznay( int argc,
               char** argv,
               AboutData const& ad,
               po::options_description const& od );

    Kovasznay( Kovasznay const& tc );

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( int iter,
                        element_U_type& U,
                        element_p_type& p );

    void oseenUpdateInit( Oseen<space_type, imOrder, ENTITY>& oseen,
                          element_U_type& U );

    void oseenUpdateIter( Oseen<space_type, imOrder, ENTITY>& oseen,
                          element_U_type& U,
                          bool updateStabilization );

private:

    double M_meshSize;
    double M_nu;

    boost::shared_ptr<export_type> M_exporter;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;

    im_type M_im;

    double M_uErrorL2;
    double M_uErrorH1;
    double M_pErrorL2;
    double M_divError;
}; // Kovasznay

} // Feel

#endif /* _KOVASZNAY_HPP_ */
