#include "../feel/feelmodels/hdg/mixedpoisson.hpp"

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-model" ,
                     "mixed-poisson-model" ,
                     "0.1",
                     "Mixed-Poisson-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;
}

void printCvg( std::ostream& out, std::vector<double> h, std::vector<std::vector<double>> err )
{
    auto errP = err[0];
    auto errRelP = err[1];
    auto errU = err[2];
    auto errRelU = err[3];
    auto errPP = err[4];
    auto errRelPP = err[5];

    boost::format fmter("%1% %|14t|");
    boost::format fmter_endl("%1%\n");

    out << fmter % "h";
    out << fmter % "errP_L2";
    out << fmter % "errP_H1";
    out << fmter % "errRelP_L2";
    out << fmter % "errRelP_H1";
    out << fmter % "errU_L2";
    out << fmter % "errU_H1";
    out << fmter % "errRelU_L2";
    out << fmter % "errRelU_H1";
    out << fmter % "errPP_L2";
    out << fmter % "errPP_H1";
    out << fmter % "errRelPP_L2";
    out << fmter_endl % "errRelPP_H1";

    int maxiter = h.size();
    for( int i = 0; i < maxiter; ++i )
    {
        out << fmter % h[i];
        out << fmter % errP[2*i];
        out << fmter % errP[2*i+1];
        out << fmter % errRelP[2*i];
        out << fmter % errRelP[2*i+1];
        out << fmter % errU[2*i];
        out << fmter % errU[2*i+1];
        out << fmter % errRelU[2*i];
        out << fmter % errRelU[2*i+1];
        out << fmter % errPP[2*i];
        out << fmter % errPP[2*i+1];
        out << fmter % errRelPP[2*i];
        out << fmter_endl % errRelPP[2*i+1];
    }
}

template<typename ExpF>
std::tuple<double,double,double,double>
computeErrors( ExpF F, Expr<GinacExVF<15>> T, bool isScalar)
{
    int quadError = ioption("benchmark.quad");
    double nL2 = normL2( _range=F.functionSpace()->template rangeElements<0>(),
                         _expr=idv(F),
                         _quad=quadError, _quad1=quadError );
    double errL2 = normL2( _range=F.functionSpace()->template rangeElements<0>(),
                           _expr=idv(F)-T,
                           _quad=quadError, _quad1=quadError );
    double errRelL2 = errL2/nL2;
    double nH1 = normH1( _range=F.functionSpace()->template rangeElements<0>(),
                         _expr=idv(F), _grad_expr=gradv(F),
                         _quad=quadError, _quad1=quadError );
    double errH1 = normH1( _range=F.functionSpace()->template rangeElements<0>(),
                           _expr=idv(F)-T,
                           _grad_expr=gradv(F)-grad<3>(T),
                           _quad=quadError, _quad1=quadError );
    double errRelH1 = errH1/nH1;

    return std::make_tuple(errL2,errRelL2,errH1,errRelH1);
}

template<typename ExpF, int M>
std::tuple<double,double,double,double>
computeErrors( ExpF F, Expr<GinacMatrix<M,1,15>> T, bool isScalar)
{
    int quadError = ioption("benchmark.quad");
    double nL2 = normL2( _range=F.functionSpace()->template rangeElements<0>(),
                         _expr=idv(F),
                         _quad=quadError, _quad1=quadError );
    double errL2 = normL2( _range=F.functionSpace()->template rangeElements<0>(),
                           _expr=idv(F)-T,
                           _quad=quadError, _quad1=quadError );
    double errRelL2 = errL2/nL2;
    double nH1 = normH1( _range=F.functionSpace()->template rangeElements<0>(),
                         _expr=idv(F), _grad_expr=gradv(F),
                         _quad=quadError, _quad1=quadError );
    double errH1 = normH1( _range=F.functionSpace()->template rangeElements<0>(),
                           _expr=idv(F)-T,
                           _grad_expr=gradv(F)-grad(T),
                           _quad=quadError, _quad1=quadError );
    double errRelH1 = errH1/nH1;

    return std::make_tuple(errL2,errRelL2,errH1,errRelH1);
}

template<int Dim, int OrderT, int GOrder = 1>
int
runApplicationMixedPoisson( std::string  const& prefix )
{
    using namespace Feel;

    typedef FeelModels::MixedPoisson<Dim,OrderT,GOrder> mp_type;
    using mesh_type = typename mp_type::mesh_type;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    std::string p = "hdg.poisson";
    auto MP = mp_type::New(p);
    auto model = MP->modelProperties().models().model();
    auto paramValues = MP->modelProperties().parameters().toParameterValues();

    auto pString = model.ptree().template get<std::string>("p_expr");
    auto pExpr = expr<15>(pString);
    pExpr.setParameterValues(paramValues);
    auto uString = model.ptree().get("u_expr","");
    auto uExpr = expr<Dim,1,15>(uString);
    uExpr.setParameterValues(paramValues);

    int maxiter = ioption("benchmark.nlevels");
    double factor = doption("benchmark.refine");
    double size = doption("benchmark.hsize");

    exporter_ptrtype e( exporter_type::New("convergence") );

    std::vector<double> h(maxiter);
    std::vector<double> normP(2*maxiter), errP(2*maxiter), errRelP(2*maxiter);
    std::vector<double> normU(2*maxiter), errU(2*maxiter,0), errRelU(2*maxiter,0);
    std::vector<double> normPP(2*maxiter), errPP(2*maxiter,0), errRelPP(2*maxiter,0);

    for( int i = 0; i < maxiter; ++i )
    {
        h[i] = size;
        auto mesh = loadMesh( _mesh=new mesh_type, _h=size );
        if( e->doExport() )
            e->step(i)->setMesh(mesh);
        MP->init(mesh);
        MP->assembleAll();
        MP->solve();
        MP->postProcess();

        auto P_h = MP->potentialField();
        auto PP_h = MP->postPotentialField();
        auto Ph = MP->potentialSpace();
        auto P_ex = Ph->element(pExpr);
        auto errTP = computeErrors(P_h,pExpr,true);
        auto errTPP = computeErrors(PP_h,pExpr,true);
        errP[2*i] = std::get<0>(errTP);
        errP[2*i+1] = std::get<1>(errTP);
        errRelP[2*i] = std::get<2>(errTP);
        errRelP[2*i+1] = std::get<3>(errTP);
        errPP[2*i] = std::get<0>(errTPP);
        errPP[2*i+1] = std::get<1>(errTPP);
        errRelPP[2*i] = std::get<2>(errTPP);
        errRelPP[2*i+1] = std::get<3>(errTPP);
        if( e->doExport() )
        {
            e->step(i)->add("P_h", P_h);
            e->step(i)->add("PP_h", PP_h);
            e->step(i)->add("P_ex", P_ex);
        }

        if( !uString.empty() )
        {
            auto U_h = MP->fluxField();
            auto Uh = MP->fluxSpace();
            auto U_ex = Uh->element(uExpr);
            auto errTU = computeErrors(U_h,uExpr,false);
            errU[2*i] = std::get<0>(errTU);
            errU[2*i+1] = std::get<1>(errTU);
            errRelU[2*i] = std::get<2>(errTU);
            errRelU[2*i+1] = std::get<3>(errTU);
            if( e->doExport() )
            {
                e->step(i)->add("U_h", U_h);
                e->step(i)->add("U_ex", U_ex);
            }
        }

        if( e->doExport() )
            e->save();

        size /= factor;
    }

    std::vector<std::vector<double>> err = {errP, errRelP, errU, errRelU, errPP, errRelPP};
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        printCvg( std::cout, h, err );
        std::ofstream file ( soption("benchmark.filename") );
        if( file )
        {
            printCvg( file, h, err );
            file.close();
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description mpoptions( "hdg.poisson options" );
    mpoptions.add( FeelModels::makeMixedPoissonOptions("","hdg.poisson") );
    mpoptions.add_options()
        ("benchmark.filename", Feel::po::value<std::string>()->default_value( "cvg.csv" ), "file to save the results" )
        ("benchmark.quad", Feel::po::value<int>()->default_value( 5 ), "quadrature order for the error")
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=mpoptions,
                           _desc_lib=FeelModels::makeMixedPoissonLibOptions("","hdg.poisson").add(feel_options())
                           );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2>, hana::int_c<1> ),
                                             hana::make_tuple("P3", hana::int_c<3>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ),
                                             hana::make_tuple("P2G2", hana::int_c<2>, hana::int_c<2> ),
                                             hana::make_tuple("P3G2", hana::int_c<3>, hana::int_c<2> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ),
                                             hana::make_tuple("P2G2", hana::int_c<2>, hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ) );
#endif

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                    [&discretization,&dimension]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                                runApplicationMixedPoisson<_dim,_torder,_gorder>( "" );
                        } );

    return 0;
}
