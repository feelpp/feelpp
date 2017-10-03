/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-16-12

       Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
       Copyright (C) CNRS

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file estimators
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \date 2016-25-10
 */

#ifndef __ESTIMATORS_HPP
#define __ESTIMATORS_HPP 1

/** include predefined feel command line options */
#include <feel/options.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <feel/feeldiscr/region.hpp>

/** include integration methods */
#include <feel/feelpoly/im.hpp>

/** include gmsh mesh importer */
#include <feel/feelfilters/gmsh.hpp>

/** include exporter factory class */
#include <feel/feelfilters/exporter.hpp>

/** include  polynomialset header */
#include <feel/feelpoly/polynomialset.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <feel/feelvf/vf.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/elementdiv.hpp>

/** include modelproperties **/
#include <feel/feelmodels/modelproperties.hpp>

#include <fstream>

namespace Feel
{
    template<int Dim, int OrderV, int OrderT, int G_order>
    class Estimators
    {
    public:
        //! numerical type is double
        typedef double value_type;

        static const bool is_P1V = (OrderV==1);
        static const bool is_P1T = (OrderT==1);

        typedef Eigen::Matrix<double, Dim, 1> vectorN_type;
        typedef Eigen::Matrix<double, Dim, Dim> matrixN_type;

        //! geometry entities type composing the mesh, here Simplex in Dimension Dim (Order G_order)
        typedef Simplex<Dim,G_order> convex_type;
        //! mesh type
        typedef Mesh<convex_type> mesh_type;
        //! mesh shared_ptr<> type
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        //POTENTIAL
        using V_space_type = Pch_type<mesh_type,OrderV>;
        using V_space_ptrtype =  boost::shared_ptr<V_space_type>;
        typedef typename V_space_type::element_type V_element_type;

        //TEMPERATURE
        using  T_space_type = Pch_type<mesh_type,OrderT>;
        typedef boost::shared_ptr<T_space_type> T_space_ptrtype;
        typedef typename T_space_type::element_type T_element_type;

        //! Scalar P0 space
        using p0_space_type = Pdh_type<mesh_type,0>;
        typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
        typedef typename p0_space_type::element_type p0_element_type;

        //! Scalar P1 space
        using p1_space_type = Pch_type<mesh_type,1>;
        typedef boost::shared_ptr<p1_space_type> p1_space_ptrtype;
        typedef typename p1_space_type::element_type p1_element_type;

        //! Vectorial P1 space
        using p1vec_space_type = Pchv_type<mesh_type,1>;
        typedef boost::shared_ptr<p1vec_space_type> p1vec_space_ptrtype;
        typedef typename p1vec_space_type::element_type p1vec_element_type;

        //! Model properties type
        typedef ModelProperties model_prop_type;
        typedef std::shared_ptr<model_prop_type> model_prop_ptrtype;
        //! Material properties type
        typedef ModelMaterials material_map_type;
        //! Boundary conditions type
        typedef BoundaryConditions condition_map_type;

        template<typename element_type>
        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> zz_estimator(const element_type& U, const mesh_ptrtype& mesh);

        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_V(const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir, const value_type& penaldir);
        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_V(const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir, const value_type& penaldir,
                                                                           mpl::bool_<true>); //Adapted for P1 order
        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_V(const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir, const value_type& penaldir,
                                                                           mpl::bool_<false>);

        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_T(const T_element_type& temperature,
                                                                           const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const material_map_type& materials,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir);

        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_T(const T_element_type& temperature,
                                                                           const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const material_map_type& materials,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir, mpl::bool_<true>); //P1

        FEELPP_DONT_INLINE
        boost::tuple<double, double, p0_element_type> residual_estimator_T(const T_element_type& temperature, const V_element_type& potential,
                                                                           const mesh_ptrtype& mesh,
                                                                           const material_map_type& materials,
                                                                           const model_prop_ptrtype& M_modelProps,
                                                                           const bool& weakdir, mpl::bool_<false>);

        // Give H1 semi norm error bound
        template<typename element_type>
        FEELPP_DONT_INLINE
        std::vector<p0_element_type> anisotropic_coeffs(const element_type& U, const mesh_ptrtype& mesh,
                                                        std::vector<vectorN_type> measures, std::vector<matrixN_type> directions);

        FEELPP_DONT_INLINE
        boost::tuple<double, p0_element_type> anisotropic_estimator_V(const V_element_type& potential,
                                                                      const mesh_ptrtype& mesh,
                                                                      const material_map_type& M_matProps,
                                                                      const model_prop_ptrtype& M_modelProps,
                                                                      const bool& weakdir, const value_type& penaldir,
                                                                      std::vector<vectorN_type> measures, std::vector<matrixN_type> directions);

        FEELPP_DONT_INLINE
        boost::tuple<double, p0_element_type> anisotropic_estimator_T(const T_element_type& temperature,
                                                                      const V_element_type& potential,
                                                                      const mesh_ptrtype& mesh,
                                                                      const material_map_type& M_matProps,
                                                                      const model_prop_ptrtype& M_modelProps,
                                                                      const bool& weakdir, const value_type& penaldir,
                                                                      std::vector<vectorN_type> measures, std::vector<matrixN_type> directions);
    };
}

template<int Dim, int OrderV, int OrderT, int G_order>
template<typename element_type>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::zz_estimator(const element_type& U, const mesh_ptrtype& mesh)
{
    //using namespace Feel::vf;

    std::vector<double> res(2);

    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh );

    auto GhUh = div( vf::sum( P1hvec, trans(gradv(U))*vf::meas()), vf::sum( P1hvec, vf::meas()*vf::one()) );
    auto eta_k_U = integrate(elements(mesh), trans(idv(GhUh) - trans(gradv(U)))*(idv(GhUh) - trans(gradv(U))), _Q<10>() ).broken(P0h).sqrt();

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    //auto npen_V=vf::project(P0h, elements(mesh), vf::nPEN() );
    //auto H1estimator_U = eta_k_U;

    auto step1_H1 = eta_k_U.pow(2);
    auto step2_H1 = step1_H1.sum();
    auto estimatorH1_U = math::sqrt(step2_H1);

    auto step1_L2 = element_product(eta_k_U,h);
    auto step2_L2 = step1_L2.pow(2);
    auto step3_L2 = step2_L2.sum();
    auto estimatorL2_U = math::sqrt(step3_L2);

    res[0] = estimatorL2_U;
    res[1] = estimatorH1_U;

    return boost::make_tuple( res[0], res[1], eta_k_U);
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_V(const V_element_type& potential, const mesh_ptrtype& mesh,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir, const value_type& penaldir)
{
    return residual_estimator_V(potential, mesh, M_modelProps, weakdir, penaldir, mpl::bool_< is_P1V >() );
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type >
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_V(const V_element_type& potential, const mesh_ptrtype& mesh,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir, const value_type& penaldir,
                                                                     mpl::bool_<true>)
{
    using namespace Feel::vf;

    std::vector<double> res(2);

    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto V = potential;

    auto jump_V = jumpv(gradv(V)); // error on jump

    //// Residual estimator for Dim=2
    auto eta_k_V = integrate( internalfaces(mesh), 0.25*h()*jump_V*jump_V ).broken(P0h).sqrt();

    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Dirichlet = g - idv(V);
                eta_k_V += integrate( markedfaces(mesh, marker),
                                      h()*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(V)*N();
                eta_k_V += integrate( markedfaces(mesh, marker),
                                      h()*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(V) + g2 - gradv(V)*N();
                eta_k_V +=  integrate( markedfaces(mesh, marker),
                                       h()*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        }
    }

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    auto H1estimator_V = eta_k_V;

    auto eta_k_L2_V = element_product(H1estimator_V,h);

    // norm of errors
    double estimatorL2_V=math::sqrt((eta_k_L2_V.pow(2)).sum());
    double estimatorH1_V=math::sqrt((eta_k_V.pow(2)).sum());

    res[0] = estimatorL2_V; // L2 norm
    res[1] = estimatorH1_V; // H1 norm

    return boost::make_tuple(res[0], res[1], eta_k_V);

}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_V(const V_element_type& potential,
                                                                     const mesh_ptrtype& mesh,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir, const value_type& penaldir,
                                                                     mpl::bool_<false>)
{
    using namespace Feel::vf;

    std::vector<double> res(2);

    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto V = potential;

    auto laplacian_V = trace(hessv(V)); //error from equation
    auto jump_V = jumpv(gradv(V)); // error on jump

    auto eta_k_V = integrate( internalfaces(mesh), 0.25*h()*jump_V*jump_V ).broken(P0h).sqrt();
    eta_k_V += integrate(elements(mesh), h()*h()*laplacian_V*laplacian_V ).broken(P0h).sqrt();

    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Dirichlet = g - idv(V);
                eta_k_V += integrate( markedelements(mesh, marker),
                                      h()*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(V)*N();
                eta_k_V += integrate( markedelements(mesh, marker),
                                      h()*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(V) + g2 - gradv(V)*N();
                eta_k_V +=  integrate( markedelements(mesh, marker),
                                       h()*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        }
    }

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    auto H1estimator_V = eta_k_V;

    auto eta_k_L2_V = element_product(H1estimator_V,h);

    // norm of errors
    double estimatorL2_V=math::sqrt((eta_k_L2_V.pow(2)).sum());
    double estimatorH1_V=math::sqrt((eta_k_V.pow(2)).sum());

    res[0] = estimatorL2_V; // L2 norm
    res[1] = estimatorH1_V; // H1 norm

    return boost::make_tuple(res[0], res[1], eta_k_V);

}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_T(const T_element_type& temperature,
                                                                     const V_element_type& potential,
                                                                     const mesh_ptrtype& mesh,
                                                                     const material_map_type& M_matProps,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir)
{
    return residual_estimator_T(temperature, potential, mesh, M_matProps, M_modelProps, weakdir, mpl::bool_< is_P1T >() );
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_T(const T_element_type& temperature,
                                                                     const V_element_type& potential,
                                                                     const mesh_ptrtype& mesh,
                                                                     const material_map_type& M_matProps,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir, mpl::bool_<true>)
{
    using namespace Feel::vf;

    std::vector<double> res(2);

    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto T = temperature;
    auto V = potential;

    auto jump_T = jumpv(gradv(T));

    auto eta_k_T = integrate(internalfaces(mesh), 0.25*h()*jump_T*jump_T, _Q<10>()).broken( P0h ).sqrt();

    auto itField = M_modelProps->boundaryConditions().find( "temperature");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Dirichlet = g - idv(T);
                eta_k_T += integrate( markedelements(mesh, marker),
                                      h()*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(T)*N();
                eta_k_T += integrate( markedelements(mesh, marker),
                                      h()*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(T) + g2 - gradv(T)*N();
                eta_k_T +=  integrate( markedelements(mesh, marker),
                                       h()*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        }
    }

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    auto H1estimator_T = eta_k_T;

    auto eta_k_H1_T = element_product(H1estimator_T,h.pow(OrderT-1));
    auto estimatorH1_T=math::sqrt((eta_k_H1_T.pow(2)).sum());

    auto eta_k_L2_T = element_product(H1estimator_T,h.pow(OrderT));
    auto estimatorL2_T=math::sqrt((eta_k_L2_T.pow(2)).sum());

    res[0] = estimatorL2_T; // L2 norm
    res[1] = estimatorH1_T; // H1 norm

    return boost::make_tuple(res[0], res[1], eta_k_T);
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::residual_estimator_T(const T_element_type& temperature,
                                                                     const V_element_type& potential,
                                                                     const mesh_ptrtype& mesh,
                                                                     const material_map_type& M_matProps,
                                                                     const model_prop_ptrtype& M_modelProps,
                                                                     const bool& weakdir, mpl::bool_<false>)
{
    using namespace Feel::vf;

    std::vector<double> res(2);

    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto T = temperature;
    auto V = potential;
    auto Xh_V = potential.functionSpace();
    auto mesh_V = Xh_V->mesh();

    auto jump_T = jumpv(gradv(T));

    auto eta_k_T = integrate(internalfaces(mesh), 0.25*h()*jump_T*jump_T, _Q<10>()).broken( P0h ).sqrt();

    double M_sigmaMax = 0;
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto sigmaProj = vf::project( _range=markedelements(mesh_V, marker),
                                       _space=Xh_V,
                                       _expr=material.getScalar("sigma0") );
        auto norm = sigmaProj.linftyNorm();
        if ( norm > M_sigmaMax )
            M_sigmaMax = norm;
    }

    // TO BE CHANGED: only for marker in M_mesh_V
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        auto alpha = material.getDouble("alpha");
        auto T0 = material.getDouble("T0");
        auto k0 = material.getDouble("k0");
        auto k = material.getScalar("k", {"T"}, {idv(T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});
        auto sigma0 = material.getDouble("sigma0")/M_sigmaMax;
        auto sigma = material.getScalar("sigma", {"T"}, {idv(T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

        auto current_laplacian = (sigma/k)*gradv(V)*trans(gradv(V))+trace(hessv(T));
        eta_k_T += integrate( markedelements( mesh, marker ),
                              h()*h()*current_laplacian*current_laplacian ).broken(P0h).sqrt();
    };

    auto itField = M_modelProps->boundaryConditions().find( "Temperature");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Dirichlet = g - idv(T);
                eta_k_T += integrate( markedelements(mesh, marker),
                                      h()*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(T)*N();
                eta_k_T += integrate( markedelements(mesh, marker),
                                      h()*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(T) + g2 - gradv(T)*N();
                eta_k_T +=  integrate( markedelements(mesh, marker),
                                       h()*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        }
    }

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    auto H1estimator_T = eta_k_T;

    auto eta_k_H1_T = element_product(H1estimator_T,h.pow(OrderT-1));
    auto estimatorH1_T=math::sqrt((eta_k_H1_T.pow(2)).sum());

    auto eta_k_L2_T = element_product(H1estimator_T,h.pow(OrderT));
    auto estimatorL2_T=math::sqrt((eta_k_L2_T.pow(2)).sum());

    res[0] = estimatorL2_T; // L2 norm
    res[1] = estimatorH1_T; // H1 norm

    return boost::make_tuple(res[0], res[1], eta_k_T);
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::anisotropic_estimator_V(const V_element_type& potential,
                                                                        const mesh_ptrtype& mesh,
                                                                        const material_map_type& M_matProps,
                                                                        const model_prop_ptrtype& M_modelProps,
                                                                        const bool& weakdir, const value_type& penaldir,
                                                                        std::vector<vectorN_type> measures, std::vector<matrixN_type> directions)
{
    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p1_space_ptrtype P1h = p1_space_type::New( mesh );
    p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh );
    auto V = potential;

    std::vector<p0_element_type> aniso_coeffs = anisotropic_coeffs(V, mesh, measures, directions);
    p0_element_type jump_BC_coeff = aniso_coeffs[0];
    p0_element_type anisoCoeff = aniso_coeffs[1];

    // Error on equation
    p0_element_type eta_k_V = integrate( elements(mesh), cst(0.) ).broken(P0h);

    // TO BE CHANGED: only for marker in M_mesh_V
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        if(!is_P1V)
            {
                auto div_sigma_grad = trace(hessv(V));
                eta_k_V += integrate( markedelements(mesh, marker), idv(anisoCoeff)*div_sigma_grad*div_sigma_grad).broken(P0h).sqrt();
            }
    };

    // Jump
    auto jump_V = jumpv( gradv(V) );
    eta_k_V += integrate( internalfaces(mesh), idv(anisoCoeff)*0.25*idv(jump_BC_coeff)*jump_V*jump_V ).broken( P0h ).sqrt();

    // Error on boundary conditions
    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());

                auto current_Dirichlet = g - idv(V);

                eta_k_V += integrate( markedelements(mesh, marker),
                                      idv(anisoCoeff)*idv(jump_BC_coeff)*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
             }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(V)*N();
                eta_k_V += integrate( markedelements(mesh, marker),
                                      idv(anisoCoeff)*idv(jump_BC_coeff)*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(V) + g2 - gradv(V)*N();
                eta_k_V +=  integrate( markedelements(mesh, marker),
                                       idv(anisoCoeff)*idv(jump_BC_coeff)*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        };
    }

    double res = math::sqrt(eta_k_V.sum());
    return boost::make_tuple( res, eta_k_V );
}

template<int Dim, int OrderV, int OrderT, int G_order>
boost::tuple<double, typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::anisotropic_estimator_T(const T_element_type& temperature,
                                                                        const V_element_type& potential,
                                                                        const mesh_ptrtype& mesh,
                                                                        const material_map_type& M_matProps,
                                                                        const model_prop_ptrtype& M_modelProps,
                                                                        const bool& weakdir, const value_type& penaldir,
                                                                        std::vector<vectorN_type> measures, std::vector<matrixN_type> directions)
{
    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p1_space_ptrtype P1h = p1_space_type::New( mesh );
    p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh );
    auto V = potential;
    auto Xh_V = potential.functionSpace();
    auto mesh_V = Xh_V->mesh();

    auto T = temperature;

    std::vector<p0_element_type> aniso_coeffs = anisotropic_coeffs(T, mesh, measures, directions);
    p0_element_type jump_BC_coeff = aniso_coeffs[0];
    p0_element_type anisoCoeff = aniso_coeffs[1];

    // Error on equation
    p0_element_type eta_k_T = integrate( elements(mesh), cst(0.) ).broken(P0h);

    double M_sigmaMax = 0;
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto sigmaProj = vf::project( _range=markedelements(mesh_V, marker),
                                      _space=Xh_V,
                                       _expr=material.getScalar("sigma0") );
        auto norm = sigmaProj.linftyNorm();
        if ( norm > M_sigmaMax )
            M_sigmaMax = norm;
    }

    // TO BE CHANGED: only for marker in M_mesh_V
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        auto alpha = material.getDouble("alpha");
        auto T0 = material.getDouble("T0");
        auto k0 = material.getDouble("k0");
        auto k = material.getScalar("k", {"T"}, {idv(T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});
        auto sigma0 = material.getDouble("sigma0")/M_sigmaMax;
        auto sigma = material.getScalar("sigma", {"T"}, {idv(T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

        if(!is_P1T)
            {
                auto div_grad_T = trace(hessv(T));
                auto rhs = sigma/k*gradv(V)*trans(gradv(V));
                eta_k_T += integrate( markedelements(mesh, marker),
                                      idv(anisoCoeff)*(div_grad_T + rhs)*(div_grad_T + rhs) ).broken(P0h).sqrt();
            }
    };

    // Jump
    auto jump_T = jumpv( gradv(T) );
    eta_k_T += integrate( internalfaces(mesh), idv(anisoCoeff)*0.25*idv(jump_BC_coeff)*jump_T*jump_T ).broken( P0h ).sqrt();

    // Error on boundary conditions
     auto itField = M_modelProps->boundaryConditions().find( "Temperature");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Dirichlet = g - idv(T);
                eta_k_T += integrate( markedelements(mesh, marker),
                                      idv(anisoCoeff)*idv(jump_BC_coeff)*current_Dirichlet*current_Dirichlet ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());

                auto current_Neumann = g - gradv(T)*N();
                eta_k_T += integrate( markedelements(mesh, marker),
                                      idv(anisoCoeff)*idv(jump_BC_coeff)*current_Neumann*current_Neumann ).broken(P0h).sqrt();
            }
        }
        itType = mapField.find( "Robin");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());

                auto current_Robin = g1*idv(T) + g2 - gradv(T)*N();
                eta_k_T +=  integrate( markedelements(mesh, marker),
                                       idv(anisoCoeff)*idv(jump_BC_coeff)*current_Robin*current_Robin ).broken(P0h).sqrt();
            }
        }
    }

    double res = math::sqrt(eta_k_T.sum());
    return boost::make_tuple( res, eta_k_T );
}

template<int Dim, int OrderV, int OrderT, int G_order>
template<typename element_type>
std::vector<typename Feel::Estimators<Dim, OrderV, OrderT, G_order>::p0_element_type>
Feel::Estimators<Dim, OrderV, OrderT, G_order>::anisotropic_coeffs(const element_type& U, const mesh_ptrtype& mesh,
                                                                   std::vector<vectorN_type> measures, std::vector<matrixN_type> directions)
{
    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p1_space_ptrtype P1h = p1_space_type::New( mesh );
    p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh );

    auto dofptItP1 = P1h->dof()->dofPointBegin();
    auto dofptEnP1 = P1h->dof()->dofPointEnd();

    // Build G matrix
    // /!\ Store only superior part (index shift)
    std::vector< std::vector<p0_element_type> > ZZ_P0(Dim);
    std::vector< std::vector<p1_element_type> > ZZ_P1(Dim);
    for(int i=0; i<Dim; i++)
        {
            ZZ_P0[i].resize(Dim - i, P0h->element());
            ZZ_P1[i].resize(Dim - i, P1h->element());
        }

    // Calculation of eta_i^ZZ on P0 space
    auto Gh_dUdi_vec = div( vf::sum( P1hvec, trans(gradv(U))*vf::meas()), vf::sum( P1hvec, vf::meas()*vf::one()) );

    // G(0,0)
    ZZ_P0[0][0] = integrate(elements(mesh), (dxv(U)-idv(Gh_dUdi_vec[Component::X]))*(dxv(U)-idv(Gh_dUdi_vec[Component::X]))).broken(P0h);
    // G(1,1)
    ZZ_P0[1][0] = integrate(elements(mesh), (dyv(U)-idv(Gh_dUdi_vec[Component::Y]))*(dyv(U)-idv(Gh_dUdi_vec[Component::Y]))).broken(P0h);
    // G(0,1)
    ZZ_P0[0][1] = integrate(elements(mesh), (dxv(U)-idv(Gh_dUdi_vec[Component::X]))*(dyv(U)-idv(Gh_dUdi_vec[Component::Y]))).broken(P0h);

    if(Dim == 3)
        {
            // G(2,2)
            ZZ_P0[2][0] = integrate(elements(mesh), (dzv(U)-idv(Gh_dUdi_vec[Component::Z]))*(dzv(U)-idv(Gh_dUdi_vec[Component::Z]))).broken(P0h);
            // G(0,2)
            ZZ_P0[0][2] = integrate(elements(mesh), (dxv(U)-idv(Gh_dUdi_vec[Component::X]))*(dzv(U)-idv(Gh_dUdi_vec[Component::Z]))).broken(P0h);
            // G(1,2)
            ZZ_P0[1][1] = integrate(elements(mesh), (dyv(U)-idv(Gh_dUdi_vec[Component::Y]))*(dzv(U)-idv(Gh_dUdi_vec[Component::Z]))).broken(P0h);
        }

    // Projectioon on P1 space (needed because measures and directions are only available on P1 space)
    for(int i=0; i<Dim; i++)
        {
            for(int j=0; j<Dim-i; j++)
                {
                    ZZ_P1[i][j] = vf::project(P1h, elements(mesh), idv(ZZ_P0[i][j]) );
                    //std::cout << "ZZ_P1 = " << ZZ_P1[i][j] << std::endl;
                }
        }

    p1_element_type minMeasP1(P1h, "minMeas");
    p1_element_type anisoCoeffP1(P1h, "anisoCoeff");

    for ( ; dofptItP1 != dofptEnP1; dofptItP1++)
        {
            auto dofptCoordP1 = get<0>(dofptItP1->second);
            auto dofptIdP1 = get<1>(dofptItP1->second);

            // Minimum of measures (lambda)
            vectorN_type currentMeas = measures[dofptIdP1];
            auto minMeas = currentMeas.minCoeff();
            minMeasP1[dofptIdP1] = minMeas;

            // Matrix G (P1) : G(i,j) = \int_K eta_i_ZZ * \int_K eta_j_ZZ
            matrixN_type GkP1;
            for(int i=0; i<Dim; i++)
                {
                    for(int j=0; j<Dim-i; j++)
                        {
                            GkP1(i,j+i) = ZZ_P1[i][j][dofptIdP1];
                            if(j != 0)
                                GkP1(j+i,i) = ZZ_P1[i][j][dofptIdP1];
                        }
                }

            // Matrix with real eigenvectors (P1)
            matrixN_type R = matrixN_type::Zero();
            for (int j=0; j<Dim; j++)
                for (int i=0; i<Dim; i++)
                    R(i,j)=real(directions[dofptIdP1](i,j));

            // Coeff \sum_i \lambda_i^2 (R_i^T G R_i)
            for(int i=0; i<Dim; i++)
                {
                    double rGr = (R.row(i)*GkP1)*R.col(i);
                    anisoCoeffP1[dofptIdP1] += measures[dofptIdP1](i)*measures[dofptIdP1](i)*math::abs(rGr);
                }
        }

    // Projection on P0 space
    p0_element_type minMeasP0 = vf::project(P0h, elements(mesh), idv(minMeasP1) );
    p0_element_type oneP0 = vf::project(P0h, elements(mesh), cst(1.) );
    p0_element_type jump_BC_coeff = div( oneP0, minMeasP0 ); //jump_coeff = 1/min(measures)
    p0_element_type anisoCoeffP0 = vf::project(P0h, elements(mesh), idv(anisoCoeffP1) );

    std::vector<p0_element_type> res(2);
    res[0] = jump_BC_coeff;
    res[1] = anisoCoeffP0;
    return res;

}

#endif //__ESTIMATORS_HPP
