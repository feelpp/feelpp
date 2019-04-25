#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/empiricalquadrature.hpp>

#include <algorithm>

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */

using namespace Feel;
using namespace Feel::vf;

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="biotsavart-nonlinear"),
                     _desc=feel_options().add(eq_options()) );

    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    int numGlobalElts = mesh->numGlobalElements();
    int numElts = mesh->numElements();
    Feel::cout << "#elts = " << numGlobalElts << " (" << numElts << ")" << std::endl;

    auto Dmu = ParameterSpace<>::New(0);
    Dmu->setDimension(1);
    auto mu_min = Dmu->element();
    mu_min << 4e3;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 6e3;
    Dmu->setMax( mu_max );

    auto mu = Dmu->element();
    mu << 4.3e3;

    if( boption("eq.test") )
    {
        auto r = markedelements(mesh,"Cu");

        // electric problem
        auto Xh = Pch<1>(mesh, r);
        auto u = Xh->element();
        auto lhs = form2(_test=Xh,_trial=Xh);
        auto rhs = form1(_test=Xh);
        lhs = integrate(_range=r, _expr=inner(gradt(u),grad(u)) );
        lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V0"), _expr=cst(0.));
        lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V1"), _expr=cst(9.));
        lhs.solve(_rhs=rhs, _solution=u);

        // biotsavart
        auto coeff = 1/(4*M_PI);
        auto mu0 = 4*M_PI; //SI unit : H.m-1 = m.kg.s-2.A-2
        // auto mu = 5.8e3;
        auto dist = inner( -P(), -P(),
                           mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
        auto ex = -mu0*coeff*cst_ref(mu(0))*cross(trans(gradv(u)),P())/(dist*dist*dist);

        auto eq = EmpiricalQuadrature(r, mu, Dmu, 2);
        eq.addExpression(ex);
        tic();
        int err = eq.offline();
        toc("offline");
        mu << 5.5e3;
        tic();
        auto a = integrate(_range=r, _expr=ex).evaluate()(2,0);
        toc("integrate");
        double d;
        if( !err )
        {
            tic();
            d = eq.evaluate();
            toc("online");
            Feel::cout << tc::green << "a = " << a << std::endl
                       << "d = " << d << tc::reset << std::endl;
        }
        else
            Feel::cout << "error " << err << std::endl;

    }
    else
    {
        // auto r = elements(mesh);
        // auto ex = expr(soption("functions.f"))*cst_ref(mu(0));
        // auto eq = EmpiricalQuadrature(r, ex, mu, Dmu);
        // tic();
        // int err = eq.offline();
        // toc("offline");
        // double d;
        // if( !err )
        // {
        //     tic();
        //     mu << 5.5e3;
        //     d = eq.evaluate();
        //     toc("online");
        // }

        // // auto q = _Q(ex.polynomialOrder());

        // tic();
        // auto a = integrate(_range=r, _expr=ex).evaluate()(0,0);
        // toc("integrate");
        // Feel::cout << tc::green << "a = " << a << std::endl
        //            << "d = " << d << tc::reset << std::endl;
    }
    // // quadrature
    // using iim_type = vf::detail::integrate_im_type<decltype(r),decltype(ex),decltype(q),decltype(q)>;
    // auto ims = iim_type::im(q,q,ex);
    // auto im = ims.first;
    // auto pts = im.points();
    // int npts = im.nPoints();
    // for ( uint32_type k = 0; k < npts; ++k )
    //     Feel::cout << "p = (" << im.point(k)(0) << "," << im.point(k)(1) << ")\t"
    //                << "w[" << k << "] = " << im.weights()[k] << std::endl;

    // auto eltForInit = boost::unwrap_ref(*boost::get<1>(r));

    // auto gm = mesh->gm();
    // auto geopc = gm->preCompute( pts );
    // auto ctx = gm->template context<vm::POINT|vm::JACOBIAN|decltype(ex)::context>( eltForInit, geopc );
    // auto expr_evaluator = ex.evaluator( mapgmc(ctx) );

    // int m = 1, j = 1;
    // int N = numElts*npts;
    // double delta = a*doption("parameters.d");
    // std::vector<std::vector<double> > eval(m*j, std::vector<double>(N, 0.0));
    // double res = 0.0;
    // int mj = 0;
    // int n = 0;
    // tic();
    // for( auto const& eltWrap : r )
    // {
    //     auto const& elt = unwrap_ref( eltWrap );
    //     if ( elt.processId() != Environment::rank() )
    //         continue;
    //     ctx->update( elt );
    //     expr_evaluator.update( vf::mapgmc( ctx ) );
    //     for ( uint16_type q=0;q<ctx->nPoints();++q )
    //     {
    //         eval[mj][n] = expr_evaluator.evalq( 2,0,q )*ctx->J(q);
    //         res += im.weights()[q]*eval[mj][n];
    //         n++;
    //     }
    // }
    // double b = 0.0;
    // mpi::all_reduce(Environment::worldComm().globalComm(), res, b, std::plus<double>());
    // toc("quadrature");
    // Feel::cout << tc::green << "b = " << b << tc::reset << std::endl;

    // tic();
    // glp_prob *lp;
    // int ia[1+m*j*N], ja[1+m*j*N];
    // double ar[1+m*j*N];
    // lp = glp_create_prob();
    // glp_set_prob_name(lp, "test");
    // glp_set_obj_dir(lp, GLP_MIN);
    // glp_add_rows(lp, m*j);
    // glp_set_row_name(lp, 1, "x_11");
    // glp_set_row_bnds(lp, 1, GLP_DB, res-delta, res+delta);
    // glp_add_cols(lp, N);
    // for( int i = 1; i <= N; ++i)
    // {
    //     glp_set_col_name(lp, i, (boost::format("p_%1%")%i).str().c_str() );
    //     glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
    //     glp_set_obj_coef(lp, i, 1.0);
    // }
    // int k = 1;
    // for( int i = 1; i <=m*j; ++i)
    // {
    //     for( int j = 1; j <= N; ++j)
    //     {
    //         ia[k] = i;
    //         ja[k] = j;
    //         ar[k] = eval[i-1][j-1];
    //         k++;
    //     }
    // }
    // glp_load_matrix(lp, m*j*N, ia, ja, ar);
    // glp_simplex(lp, NULL);
    // toc("init");
    // tic();
    // double z = glp_get_obj_val(lp);
    // toc("lp");
    // std::vector<std::pair<int,double> > weights;
    // for( int i = 1; i <= N; ++i)
    //     if( glp_get_col_prim(lp, i) > 1e-12 )
    //         weights.push_back(std::make_pair(i-1, glp_get_col_prim(lp, i)) );
    // glp_delete_prob(lp);

    // Feel::cout << "non zero = " << weights.size() << "/" << N << std::endl;

    // res = 0.0;
    // tic();
    // for( int n = 0; n < weights.size(); ++n )
    //     res += weights[n].second*eval[0][weights[n].first];
    // double c = 0.0;
    // mpi::all_reduce(Environment::worldComm().globalComm(), res, c, std::plus<double>());
    // toc("eq");
    // Feel::cout << tc::green << "c = " << c << tc::reset << std::endl;

    return 0;
}
