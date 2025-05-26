/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-28

  Copyright (C) 2011-2014 Feel++ Consortium

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
//#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE test_bdf2
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelts/bdfadaptive.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
//using Feel::project;



std::map<int, std::string> exact_expr_k = {
    {1, "x+1 + 2*t :t:x"},
    {2, "x+1 + 2*t + 3*t^2 :t:x"},
    {3, "x+1 + 2*t + 3*t^2 + 4*t^3 :t:x"},
    {4, "x+1 + 2*t + 3*t^2 + 4*t^3 + 5*t^4 :t:x"}
};

std::map<int, std::string> exact_d1_k = {
    {1, "2 :t:x"},
    {2, "2 + 6*t :t:x"},
    {3, "2 + 6*t + 12*t^2 :t:x"},
    {4, "2 + 6*t + 12*t^2 + 20*t^3 :t:x"}
};

std::map<int, std::string> exact_d2_k= {
    {1, "0 :t:x"},
    {2, "6 :t:x"},
    {3, "6 + 24*t :t:x"},
    {4, "6 + 24*t + 60*t^2 :t:x"}
};

std::map<int, std::string> exprs_km1 = {
    {1, "1 :t:x"},
    {2, "1 + 2*t :t:x"},
    {3, "1 + 2*t + 3*t^2 :t:x"},
    {4, "1 + 2*t + 3*t^2 + 4*t^3 :t:x"}
};
std::map<int, std::string> d1_exprs_km1 = {
    {1, "0 :t:x"},
    {2, "2 :t:x"},
    {3, "2 + 6*t :t:x"},
    {4, "2 + 6*t + 12*t^2 :t:x"}
};
std::map<int, std::string> d2_exprs_km1 = {
    {1, "0 :t:x"},
    {2, "0 :t:x"},
    {3, "6 :t:x"},
    {4, "6 + 24*t :t:x"}
};

void check_bdf( std::map<int,std::string> exprs, std::map<int,std::string> d1_exprs, std::map<int,std::string> d2_exprs, bool check_extrapolation = false )
{
    using namespace Feel;
    using namespace Feel::vf;

    for (int bdf_order = 1; bdf_order <= 4; ++bdf_order)
    {
        BOOST_TEST_MESSAGE(fmt::format("[order {}] Testing BDF order {}", bdf_order, bdf_order));

        auto mesh = createGMSHMesh(_mesh=new Mesh<Simplex<1>>(), _desc=domain(_name="unit", _dim=1));
        auto Xh = Pch<1>(mesh);

        double dt = 0.25;
        int n_history = bdf_order + 2;
        std::vector<double> times(n_history);
        for (int i = 0; i < n_history; ++i)
            times[i] = 1.0 - i * dt;
        BOOST_TEST_MESSAGE(fmt::format("[order {}] Time steps initialized", bdf_order));
        std::vector<decltype(Xh->element())> states(n_history);
        for (int i = 0; i < n_history; ++i)
        {
            auto ue_expr = expr(exprs[bdf_order]);
            ue_expr.setParameterValues({{"t", times[i]}});
            states[i] = Xh->element();
            states[i].on(_range=elements(mesh), _expr=ue_expr);
        }
        BOOST_TEST_MESSAGE(fmt::format("[order {}] States initialized", bdf_order));
        auto thebdf = bdf(_space=Xh, _name="bdf_test", _order=bdf_order, _time_step=dt);
        if (bdf_order == 1)
            thebdf->setHistory(states[0], states[1]);
        else if (bdf_order == 2)
            thebdf->setHistory(states[0], states[1], states[2], states[3]);
        else if (bdf_order == 3)
            thebdf->setHistory(states[0], states[1], states[2], states[3], states[4]);
        else if (bdf_order == 4)
            thebdf->setHistory(states[0], states[1], states[2], states[3], states[4], states[5]);
        BOOST_TEST_MESSAGE(fmt::format("[order {}] BDF initialized", bdf_order));

        // Check first derivative at t^{n}
        auto du_dt = thebdf->firstDerivative();
        auto d1_expr = expr(d1_exprs[bdf_order]);
        d1_expr.setParameterValues({{"t", times[0]}});
        auto d1_ref = project(_space=Xh, _expr=d1_expr);
        auto err1 = normL2(_range=elements(mesh), _expr=idv(du_dt) - idv(d1_ref));
        BOOST_CHECK_SMALL(err1, 1e-12);
        BOOST_TEST_MESSAGE(fmt::format("[order {}] First derivative passed", bdf_order));

        // check polyDeriv()
        auto  dtu = Xh->element();
        auto ue_expr = expr(exprs[bdf_order]);
        ue_expr.setParameterValues({{"t", times[0]+dt}});
        d1_expr.setParameterValues({{"t", times[0]+dt}});
        d1_ref = project(_space=Xh, _expr=ue_expr);
        dtu.on(_range=elements(mesh), _expr=thebdf->polyDerivCoefficient(0)*idv(d1_ref)-idv(thebdf->polyDeriv()));
        auto dtu_ref = project(_space=Xh, _expr=d1_expr);
        dtu.printMatlab("dtu");
        dtu_ref.printMatlab("dtu_ref");
        auto err_poly = normL2(_range=elements(mesh), _expr=idv(dtu) - idv(dtu_ref));
        BOOST_CHECK_SMALL(err_poly, 1e-12);
        BOOST_TEST_MESSAGE(fmt::format("[order {}] Poly derivative passed", bdf_order));
        // Check second derivative if available
        if (bdf_order >= 2)
        {
            auto d2u_dt2 = thebdf->secondDerivative();
            auto d2_expr = expr(d2_exprs[bdf_order]);
            d2_expr.setParameterValues({{"t", times[0]}});
            auto d2_ref = project(_space=Xh, _expr=d2_expr);
            auto err2 = normL2(_range=elements(mesh), _expr=idv(d2u_dt2) - idv(d2_ref));
            BOOST_CHECK_SMALL(err2, 1e-12);
            BOOST_TEST_MESSAGE(fmt::format("[order {}] Second derivative passed", bdf_order));
        }
        // Check extrapolation
        if (check_extrapolation)
        {
            double t = times[0]+dt;
            auto extrap = thebdf->extrapolation();
            auto extrap_expr = expr(exprs[bdf_order]);
            extrap_expr.setParameterValues({{"t", t}});
            auto extrap_ref = project(_space=Xh, _expr=extrap_expr);
            auto err_extrap = normL2(_range=elements(mesh), _expr=idv(extrap) - idv(extrap_ref));
            BOOST_CHECK_SMALL(err_extrap, 1e-12);
            BOOST_TEST_MESSAGE(fmt::format("[order {}] Extrapolation passed", bdf_order));
        }
    }
}

template<int Dim>
class TestPDEWithBDF
{
public:
    // pass maps to exact expressions
    using expr_map_t = std::map<int, std::string>;
    expr_map_t exprs_;
    expr_map_t d1_exprs_;
    expr_map_t d2_exprs_;
    TestPDEWithBDF(expr_map_t const& exprs,
                   expr_map_t const& d1_exprs,
                   expr_map_t const& d2_exprs)
        :
        exprs_(exprs),
        d1_exprs_(d1_exprs),
        d2_exprs_(d2_exprs)
    {
    }
    void run(int bdf_order)
    {
        using namespace Feel;
        using namespace Feel::vf;

        BOOST_TEST_MESSAGE(fmt::format("Running PDE test with BDF order {}", bdf_order));

        auto mesh = createGMSHMesh(_mesh=new Mesh<Simplex<Dim>>,
                                   _desc=domain(_name=(boost::format("unit-%1%d") % Dim).str(), _dim=Dim));
        auto Xh = Pch<2>(mesh);
        auto u = Xh->element();
        auto ue = Xh->element();
        auto v = Xh->element();
        auto solution = Xh->element();

        auto mybdf = bdf(_space=Xh, _name="mybdf", _order=bdf_order, _time_step=0.1, _final_time=0.1);

        std::string g = exprs_[bdf_order];
        auto ue_g = expr(g);
        auto fe   = -laplacian(ue_g) + diff(ue_g, "t");

        auto e = exporter(_mesh=mesh);

        std::vector<decltype(Xh->element())> states(mybdf->priorTimes().size());
        for (auto const& time : mybdf->priorTimes())
        {
            if (Environment::worldComm().isMasterRank())
                BOOST_TEST_MESSAGE( fmt::format("Initializing t={} (index {})", time.second, time.first));

            ue_g.setParameterValues({{"t", time.second}});
            states[time.first] = Xh->element();
            states[time.first].on(_range=elements(mesh), _expr=ue_g);
        }

        BOOST_TEST_MESSAGE(fmt::format("BDF initialized with order {} time: {}", bdf_order, mybdf->time()));

        auto d1_expr = expr(d1_exprs_[bdf_order]);
        d1_expr.setParameterValues({{"t", mybdf->time()}});
        fe.setParameterValues({{"t", mybdf->time()}});
        BOOST_CHECK_CLOSE(fe.evaluate()(0, 0), d1_expr.evaluate()(0, 0), 1e-12);


        ue_g.setParameterValues({{"t", mybdf->time()}});
        solution.on(_range=elements(mesh), _expr=ue_g);
        ue.on(_range=elements(mesh), _expr=ue_g);
        auto error = project(_space=Xh, _expr=idv(ue) - idv(solution));
        e->step(0)->add("exact", ue);
        e->step(0)->add("solution", solution);
        e->step(0)->add("error", error);
        e->save();

        double maxerror = error.linftyNorm();
        if (Environment::worldComm().isMasterRank())
            std::cout << "Initial max error: " << maxerror << "\n";

        for (mybdf->start(states); !mybdf->isFinished(); mybdf->next(solution))
        {
            BOOST_TEST_MESSAGE(fmt::format("BDF step at time {}", mybdf->time()));
            ue_g.setParameterValues({{"t", mybdf->time()}});
            fe.setParameterValues({{"t", mybdf->time()}});
            d1_expr.setParameterValues({{"t", mybdf->time()}});
            BOOST_CHECK_CLOSE(fe.evaluate()(0, 0), d1_expr.evaluate()(0, 0), 1e-12);

            auto ft = form1(_test=Xh);
            ft = integrate(_range=elements(mesh), _expr=(fe + idv(mybdf->polyDeriv())) * id(u));

            auto at = form2(_test=Xh, _trial=Xh);
            at = integrate(_range=elements(mesh),
                           _expr=gradt(u) * trans(grad(v)) +
                                 mybdf->polyDerivCoefficient(0) * idt(u) * id(u));

            at += on(_range=boundaryfaces(mesh), _element=solution, _rhs=ft, _expr=ue_g);
            at.solve(_solution=solution, _rhs=ft);

            ue = project(_space=Xh, _expr=ue_g);
            auto error = project(_space=Xh, _expr=idv(ue) - idv(solution));
            //auto error = project(_space=Xh, _expr=idv(mybdf->firstDerivative(ue)) - (mybdf->polyDerivCoefficient(0)*idv(ue) - idv(mybdf->polyDeriv())));
            double maxerror = error.linftyNorm();
            ue.printMatlab("ue");
            solution.printMatlab("solution");
            error.printMatlab("error");
            if (Environment::worldComm().isMasterRank())
                std::cout << "Max error at time " << mybdf->time() << ": " << maxerror << "\n";

            BOOST_CHECK_SMALL(maxerror, 1e-9);

            e->step(mybdf->time())->add("exact", ue);
            e->step(mybdf->time())->add("solution", solution);
            e->step(mybdf->time())->add("error", error);
            e->save();
        }
    }
};



FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( bdf2 )



BOOST_AUTO_TEST_CASE( test_bdf_order_km1 )
{
    check_bdf( exprs_km1, d1_exprs_km1, d2_exprs_km1, true );

}
BOOST_AUTO_TEST_CASE( test_bdf_order_k )
{
    check_bdf( exact_expr_k, exact_d1_k, exact_d2_k, false );
}

BOOST_AUTO_TEST_CASE( test_bdf_pde_km1 )
{
    TestPDEWithBDF<2> test1(exprs_km1, d1_exprs_km1, d2_exprs_km1);
    test1.run(1);
    test1.run(2);
    test1.run(3);
    test1.run(4);
}

BOOST_AUTO_TEST_CASE( test_bdf_pde_k )
{
    TestPDEWithBDF<2> test1(exact_expr_k, exact_d1_k, exact_d2_k);
    test1.run(1);
    test1.run(2);
    test1.run(3);
    test1.run(4);
}

BOOST_AUTO_TEST_CASE( test_bdfadaptive_constant_dt_coefficients )
{
    using namespace Feel;

    // simple 1D P1 space
    auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<1>>(), _desc=domain(_name="unit",_dim=1) );
    auto Xh = Pch<1>( mesh );

    constexpr double dt = 0.25;
    constexpr double t0 = 1.0;
    // classical BDF weights for k=1..4
    std::map<int,std::vector<double>> expected = {
        {1, {  1.0,   -1.0           } },
        {2, {  3./2., -2.0,  1./2.   } },
        {3, { 11./6., -3.0,  3./2., -1./3.} },
        {4, { 25./12.,-4.0,  3.0,  -4./3., 1./4. } }
    };

    for ( int k = 1; k <= 4; ++k )
    {
        // build time stamps [ t0, t0-dt, … ]
        std::vector<double> times( k+1 );
        for ( int i = 0; i <= k; ++i ) times[i] = t0 - i*dt;

        // dummy history (values not used for coeffs)
        std::vector<decltype(Xh->element())> H( k+1 );
        for ( int i = 0; i <= k; ++i ) H[i] = Xh->element();

        auto bdfA = std::make_shared< BdfAdaptive<decltype(Xh)::element_type> >( Xh, "adaptive_coeffs" );
        bdfA->setTimeStamps( times );
        bdfA->setHistory( H );
        bdfA->setCurrentTimeStep( dt );

        BOOST_TEST_MESSAGE( "[order " << k << "] checking classical BDF weights" );
        for ( int j = 0; j <= k; ++j )
        {
            // polyDerivCoefficient() gives α_j / dt, so multiply back by dt
            double alpha = bdfA->polyDerivCoefficient(j) * dt;
            double alpha_ref = expected[k][j];
            BOOST_CHECK_CLOSE( alpha, alpha_ref, 1e-12 );
        }
    }
}
#if 0
BOOST_AUTO_TEST_CASE(test_adaptive_error_estimation_constant_dt)
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = createGMSHMesh(_mesh=new Mesh<Simplex<1>>(), _desc=domain(_name="unit", _dim=1));
    auto Xh = Pch<1>(mesh);

    constexpr int order = 3;
    constexpr double dt = 0.25;

    // Time stamps: t_0, t_1, ..., t_k
    std::vector<double> times = {
        1.0, 0.75, 0.5, 0.25
    };

    auto u_expr = expr("1 + 2*t + 3*t^2 :t");

    std::vector<decltype(Xh->element())> states;
    for (auto t : times)
    {
        auto u = Xh->element();
        u_expr.setParameterValues({{"t", t}});
        u.on(_range=elements(mesh), _expr=u_expr);
        states.push_back(u);
    }

    auto bdf = std::make_shared<BdfAdaptive<decltype(Xh)::element_type>>(Xh, "adaptive_bdf");

    bdf->setTimeStamps(times);
    bdf->setHistory(states);

    double t_next = times[0] + dt;

    double err = bdf->estimateError(t_next);
    BOOST_CHECK_SMALL(err, 1e-12);

    auto adapt = bdf->checkAndAdapt(t_next, /*tol=*/1e-6);
    BOOST_CHECK(adapt.accepted);
    BOOST_TEST_MESSAGE(fmt::format("estimated error = {:.3e}, suggested dt = {:.3e}", adapt.err, adapt.new_dt));
}

BOOST_AUTO_TEST_CASE(test_adaptive_error_estimation_variable_dt)
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = createGMSHMesh(_mesh=new Mesh<Simplex<1>>(), _desc=domain(_name="unit", _dim=1));
    auto Xh = Pch<1>(mesh);

    constexpr int order = 3;
    std::vector<double> times = {
        1.0, 0.7, 0.4, 0.2  // non-uniform
    };

    auto u_expr = expr("1 + 2*t + 3*t^2 :t");

    std::vector<decltype(Xh->element())> states;
    for (auto t : times)
    {
        auto u = Xh->element();
        u_expr.setParameterValues({{"t", t}});
        u.on(_range=elements(mesh), _expr=u_expr);
        states.push_back(u);
    }

    auto bdf = std::make_shared<BdfAdaptive<decltype(Xh)::element_type>>(Xh, "adaptive_bdf");
    bdf->setTimeStamps(times);
    bdf->setHistory(states);

    double t_next = 1.2;  // forward step

    double err = bdf->estimateError(t_next);
    BOOST_CHECK_SMALL(err, 1e-12);

    auto adapt = bdf->checkAndAdapt(t_next, 1e-6);
    BOOST_CHECK(adapt.accepted);
    BOOST_TEST_MESSAGE(fmt::format("non-uniform dt error = {:.3e}, dt = {:.3e}", adapt.err, adapt.new_dt));
}
#endif
BOOST_AUTO_TEST_SUITE_END()
