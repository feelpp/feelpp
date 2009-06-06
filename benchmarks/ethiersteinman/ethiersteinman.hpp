/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-04-03

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
   \file ethiersteinman.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-04-03
 */

#ifndef _ETHIER_STEINMAN_HPP_
#define _ETHIER_STEINMAN_HPP_

#include <boost/timer.hpp>

#include <life/lifecore/applicationserial.hpp>
#include <life/lifealg/backendgmm.hpp>
#include <life/lifealg/backend_adaptive_reuse_pc.hpp>
#include <life/lifepoly/im.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifevf/vf.hpp>
#include <life/lifefilters/exporterensight.hpp>

#include "oseen.hpp"

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description
        ethiersteinmanoptions("EthierSteinman options");
    ethiersteinmanoptions.add_options()
        ("dt", Life::po::value<double>()->default_value( 0.025 ),
         "time step value")
        ("ft", Life::po::value<double>()->default_value( 0.1 ),
         "Final time value")
        ("mu", Life::po::value<double>()->default_value( 1.0 ),
         "viscosity of fluid")
        ("hsize", Life::po::value<double>()->default_value( 0.1 ),
         "first h value to start convergence")
        ("bccoeffdiff", Life::po::value<double>()->default_value( 100.0 ),
         "coefficient for diffusive weak Dirichlet boundary terms")
        ("bccoeffconv", Life::po::value<double>()->default_value( 100.0 ),
         "coefficient for convective weak Dirichlet boundary terms")
        ("fixpointtol", Life::po::value<double>()->default_value( 0.025 ),
         "Convergence tolerance for fixed point sub-iterations")
        ("maxsubiter", Life::po::value<int>()->default_value( 1 ),
         "maximal number of fixed point sub-iterations")
        ("divdivcoeff", Life::po::value<double>()->default_value( 0.0 ),
         "coefficient for ( div u , div v ) term in matrix")
        ("stabcoeff-div", Life::po::value<double>()->default_value( 0.0 ),
         "interior penalty stabilization coefficient for divergence")
        ("stabcoeff-p", Life::po::value<double>()->default_value( 0.1 ),
         "interior penalty stabilization coefficient for pressure")
        ("stabtheta", Life::po::value<double>()->default_value( -1.0 ),
         "stabilization decoupling coefficient (no decoupling for values < 0)")
        ("epscompress", Life::po::value<double>()->default_value( 0.0 ),
         "coefficient of compressibility term")
        ("export", Life::po::value<int>()->default_value( 0 ),
         "stride for result export (0=no export)")
        ("bdf1", "use BDF1 time stepping")
        ;

    Life::po::options_description solveroptions("algebraic solver options");
    solveroptions.add_options()
        ("tolerance", Life::po::value<double>()->default_value( 2.e-10 ),
         "solver tolerance")
        ("verbose", Life::po::value<int>()->default_value( 0 ),
         "(=0,1,2) print solver iterations")
        ("maxiter", Life::po::value<int>()->default_value( 1000 ),
         "set maximum number of iterations")
        ("fillin", Life::po::value<int>()->default_value( 20 ),
         "fill-in for incomplete factorizations")
        ("threshold", Life::po::value<double>()->default_value( 1.e-4 ),
         "threshold for incomplete factorizations")
        ("solver", Life::po::value<std::string>()->default_value( "gmres" ),
         "solver type")
        ("precond", Life::po::value<std::string>()->default_value( "ilutp" ),
         "preconditioner type")
        ("restart", Life::po::value<int>()->default_value( 200 ),
         "gmres restart")
        ("direct", "use direct solver")
      ;
    return ethiersteinmanoptions.add( solveroptions );
}

inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "ethiersteinman" ,
                            "ethiersteinman" ,
                            "0.1",
                            "Ethier Steinman Benchmark",
                            Life::AboutData::License_GPL,
                            "Copyright (c) 2007 EPFL");

    about.addAuthor("Christoph Winkelmann", "developer",
                    "christoph.winkelmann@epfl.ch", "");
    return about;

}

namespace Life
{
class EthierSteinman
    :
        public ApplicationSerial
{
    typedef ApplicationSerial super;
public:

    // -- TYPEDEFS --
    static const uint16_type uOrder = 1;
    static const uint16_type pOrder = 1;

    static const int Dim = 2; // was a template parameter

    typedef double value_type;

    /* entity */
#define ENTITY Simplex

    /* mesh */
    typedef Mesh<GeoEntity<ENTITY<Dim, 1,Dim> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /* bases */
    typedef fusion::vector<fem::Lagrange<Dim, uOrder,
                                         Vectorial, Continuous,
                                         double, ENTITY> > basis_U_type;
    typedef fusion::vector<fem::Lagrange<Dim, uOrder,
                                         Scalar, Continuous,
                                         double, ENTITY> > basis_u_type;
    typedef fusion::vector<fem::Lagrange<Dim, pOrder,
                                         Scalar, Continuous,
                                         double, ENTITY> > basis_p_type;
    typedef fusion::vector<fem::Lagrange<Dim, 0,
                                         Scalar, Discontinuous,
                                         double, ENTITY> > basis_i_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_U_type, value_type> space_U_type;
    typedef FunctionSpace<mesh_type, basis_u_type, value_type> space_u_type;
    typedef FunctionSpace<mesh_type, basis_p_type, value_type> space_p_type;
    typedef FunctionSpace<mesh_type, basis_i_type, value_type> space_i_type;
    typedef boost::shared_ptr<space_U_type> space_U_ptrtype;
    typedef boost::shared_ptr<space_u_type> space_u_ptrtype;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef space_U_type::element_type element_U_type;
    typedef space_u_type::element_type element_u_type;
    typedef space_p_type::element_type element_p_type;
    typedef space_i_type::element_type element_i_type;

    /* quadrature for postprocessing */
    static const uint16_type imOrder = 3*uOrder+2; // 5 / 20
    typedef IM<Dim, imOrder, value_type, ENTITY> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    // linear algebra backends
    typedef BackendAdaptiveReusePC<BackendGmm<value_type> > backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    EthierSteinman( int argc, char** argv, AboutData const& ad );

    EthierSteinman( int argc,
                    char** argv,
                    AboutData const& ad,
                    po::options_description const& od );

    EthierSteinman( EthierSteinman const& tc );

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize, double R=0.0 );

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

    void oseenUpdateInit( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                          mesh_ptr_type mesh,
                          value_type dt );

    void oseenUpdateBdf1( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                          value_type dt,
                          element_u_type& uxn,
                          element_u_type& uyn,
                          element_u_type& uzn,
                          element_u_type& ux,
                          element_u_type& uy,
                          element_u_type& uz,
                          element_u_type& ux0,
                          element_u_type& uy0,
                          element_u_type& uz0,
                          bool updateStabilization,
                          element_p_type& pn );

    void oseenUpdateBdf2Trans( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                               value_type dt,
                               element_u_type& uxn,
                               element_u_type& uyn,
                               element_u_type& uzn,
                               element_u_type& ux,
                               element_u_type& uy,
                               element_u_type& uz,
                               element_u_type& uxo,
                               element_u_type& uyo,
                               element_u_type& uzo,
                               element_u_type& ux0,
                               element_u_type& uy0,
                               element_u_type& uz0,
                               bool updateStabilization,
                               element_p_type& pn );

    void oseenUpdateBdf2( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                          value_type dt,
                          element_u_type& uxn,
                          element_u_type& uyn,
                          element_u_type& uzn,
                          element_u_type& ux,
                          element_u_type& uy,
                          element_u_type& uz,
                          element_u_type& uxo,
                          element_u_type& uyo,
                          element_u_type& uzo,
                          element_u_type& ux0,
                          element_u_type& uy0,
                          element_u_type& uz0,
                          bool updateStabilization,
                          element_p_type& pn );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( int iter,
                        double time,
                        element_U_type& U,
                        element_u_type& ux,
                        element_u_type& uy,
                        element_u_type& uz,
                        element_p_type& p,
                        OperatorLinear<space_u_type, space_u_type, backend_type>& massU,
                        OperatorLinear<space_u_type, space_u_type, backend_type>& laplaceU,
                        OperatorLinear<space_p_type, space_p_type, backend_type>& massP,
                        element_u_type& ux0,
                        element_u_type& uy0,
                        element_u_type& uz0,
                        element_p_type& p0 );

private:

    double M_meshSize;
    double M_bcCoeffDiff;
    double M_bcCoeffConv;
    double M_mu;

    boost::shared_ptr<export_type> M_exporter;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;

    im_type M_im;

    double M_uErrorL2;
    double M_uErrorH1;
    double M_pErrorL2;
    double M_divError;
}; // EthierSteinman

} // Life

#endif /* _ETHIER_STEINMAN_HPP_ */
