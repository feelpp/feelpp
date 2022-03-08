#include "thermoelectric.hpp"

using namespace Feel;

void printCvg( std::ostream& out, std::vector<double> h, std::vector<std::vector<double>> err )
{
    auto errT = err[0];
    auto errRelT = err[1];
    auto errV = err[2];
    auto errRelV = err[3];
    auto errC = err[4];
    auto errRelC = err[5];

    boost::format fmter("%1% %|14t|");
    boost::format fmter_endl("%1%\n");

    out << fmter % "h";
    out << fmter % "errT_L2";
    out << fmter % "errT_H1";
    out << fmter % "errRelT_L2";
    out << fmter % "errRelT_H1";
    out << fmter % "errV_L2";
    out << fmter % "errV_H1";
    out << fmter % "errRelV_L2";
    out << fmter % "errRelV_H1";
    out << fmter % "errC_L2";
    out << fmter % "errC_H1";
    out << fmter % "errRelC_L2";
    out << fmter_endl % "errRelC_H1";

    int maxiter = h.size();
    for( int i = 0; i < maxiter; ++i )
    {
        out << fmter % h[i];
        out << fmter % errT[2*i];
        out << fmter % errT[2*i+1];
        out << fmter % errRelT[2*i];
        out << fmter % errRelT[2*i+1];
        out << fmter % errV[2*i];
        out << fmter % errV[2*i+1];
        out << fmter % errRelV[2*i];
        out << fmter % errRelV[2*i+1];
        out << fmter % errC[2*i];
        out << fmter % errC[2*i+1];
        out << fmter % errRelC[2*i];
        out << fmter_endl % errRelC[2*i+1];
    }
}

template <int Dim, int OrderT, int OrderG>
int
runApplicationThermoElectric()
{
    using model_type = ThermoElectricHDG<Dim, OrderT, OrderT, OrderG>;
    using mesh_type = typename model_type::mesh_type;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    std::shared_ptr<model_type> thermoElectric( new model_type() );
    auto heat = thermoElectric->thermo();
    auto heatModel = heat->modelProperties().models().model();
    auto heatParamValues = heat->modelProperties().parameters().toParameterValues();
    auto electric = thermoElectric->electro();
    auto electricModel = electric->modelProperties().models().model();
    auto electricParamValues = electric->modelProperties().parameters().toParameterValues();

    auto tString = heatModel.ptree().template get<std::string>("p_expr");
    auto tExpr = expr<15>(tString);
    tExpr.setParameterValues(heatParamValues);
    auto vString = electricModel.ptree().template get<std::string>("p_expr");
    auto vExpr = expr<15>(vString);
    vExpr.setParameterValues(electricParamValues);
    auto cString = electricModel.ptree().get("u_expr", "");
    auto cExpr = expr<Dim,1,15>(cString);
    cExpr.setParameterValues(electricParamValues);

    int maxiter = ioption("benchmark.nlevels");
    double factor = doption("benchmark.refine");
    double size = doption("benchmark.hsize");
    int quadError = ioption("benchmark.quad");

    exporter_ptrtype e( exporter_type::New("convergence") );

    std::vector<double> h(maxiter);
    std::vector<double> normT(2*maxiter), errT(2*maxiter), errRelT(2*maxiter);
    std::vector<double> normV(2*maxiter), errV(2*maxiter), errRelV(2*maxiter);
    std::vector<double> normC(2*maxiter, 0), errC(2*maxiter, 0), errRelC(2*maxiter, 0);

    for( int i = 0; i < maxiter; ++i )
    {
        h[i] = size;
        auto mesh = loadMesh( _mesh=new mesh_type, _h=size );
        if( e->doExport() )
            e->step(i)->setMesh(mesh);

        thermoElectric->init(mesh);
        thermoElectric->solve();

        auto T_h = heat->potentialField();
        auto Th = heat->potentialSpace();
        auto T_ex = Th->element(tExpr);
        normT[2*i] = normL2( _range=Th->template rangeElements<0>(),
                             _expr=idv(T_h),
                             _quad=quadError, _quad1=quadError );
        errT[2*i] = normL2( _range=Th->template rangeElements<0>(),
                            _expr=idv(T_h)-tExpr,
                            _quad=quadError, _quad1=quadError );
        errRelT[2*i] = errT[2*i]/normT[2*i];
        normT[2*i+1] = normH1( _range=Th->template rangeElements<0>(),
                               _expr=idv(T_h), _grad_expr=gradv(T_h),
                               _quad=quadError, _quad1=quadError );
        errT[2*i+1] = normH1( _range=Th->template rangeElements<0>(),
                              _expr=idv(T_h)-tExpr,
                              _grad_expr=gradv(T_h)-grad<mesh_type::nRealDim>(tExpr),
                              _quad=quadError, _quad1=quadError );
        errRelT[2*i+1] = errT[2*i+1]/normT[2*i+1];
        if( e->doExport() )
        {
            e->step(i)->add("T_h", T_h);
            e->step(i)->add("T_ex", T_ex);
        }

        auto V_h = electric->potentialField();
        auto Vh = electric->potentialSpace();
        auto V_ex = Vh->element(vExpr);
        normV[2*i] = normL2( _range=Vh->template rangeElements<0>(),
                             _expr=idv(V_h),
                             _quad=quadError, _quad1=quadError );
        errV[2*i] = normL2( _range=Vh->template rangeElements<0>(),
                            _expr=idv(V_h)-vExpr,
                            _quad=quadError, _quad1=quadError );
        errRelV[2*i] = errV[2*i]/normV[2*i];
        normV[2*i+1] = normH1( _range=Vh->template rangeElements<0>(),
                               _expr=idv(V_h),
                               _grad_expr=gradv(V_h),
                               _quad=quadError, _quad1=quadError );
        errV[2*i+1] = normH1( _range=Vh->template rangeElements<0>(),
                              _expr=idv(V_h)-idv(V_ex),
                              _grad_expr=gradv(V_h)-grad<mesh_type::nRealDim>(vExpr),
                              _quad=quadError, _quad1=quadError );
        errRelV[2*i+1] = errV[2*i+1]/normV[2*i+1];
        if( e->doExport() )
        {
            e->step(i)->add("V_h", V_h);
            e->step(i)->add("V_ex", V_ex);
        }

        auto C_h = electric->fluxField();
        if( !cString.empty() )
        {
            auto Ch = electric->fluxSpace();
            auto C_ex = Ch->element(cExpr);
            normC[2*i] = normL2( _range=Ch->template rangeElements<0>(),
                                 _expr=idv(C_h),
                                 _quad=quadError, _quad1=quadError );
            errC[2*i] = normL2( _range=Ch->template rangeElements<0>(),
                                _expr=idv(C_h)-cExpr,
                                _quad=quadError, _quad1=quadError );
            errRelC[2*i] = errC[2*i]/normC[2*i];
            normC[2*i+1] = normH1( _range=Ch->template rangeElements<0>(),
                                   _expr=idv(C_h), _grad_expr=gradv(C_h),
                                   _quad=quadError, _quad1=quadError );
            errC[2*i+1] = normH1( _range=Ch->template rangeElements<0>(),
                                  _expr=idv(C_h)-idv(C_ex),
                                  _grad_expr=gradv(C_h)-grad<mesh_type::nRealDim>(cExpr),
                                  _quad=quadError, _quad1=quadError );
            errRelC[2*i+1] = errC[2*i+1]/normC[2*i+1];
            if( e->doExport() )
            {
                e->step(i)->add("C_h", C_h);
                e->step(i)->add("C_ex", C_ex);
            }
        }

        if( e->doExport() )
            e->save();
        size /= factor;
    }

    std::vector<std::vector<double>> err = {errT, errRelT, errV, errRelV, errC, errRelC};
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

int main(int argc, char** argv) {
    po::options_description thermoelectricoptions( "application thermo-electric options" );
    thermoelectricoptions.add( makeThermoElectricHDGOptions() );
    thermoelectricoptions.add_options()
        ("benchmark.filename", Feel::po::value<std::string>()->default_value( "cvg.csv" ), "file to save the results" )
        ("benchmark.quad", Feel::po::value<int>()->default_value( 5 ), "quadrature order for the error")
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=thermoelectricoptions,
                     _about=about(_name="cvg_thermoelectric",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");
    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1G1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P2G1", hana::int_c<2>, hana::int_c<1> ),
                                             hana::make_tuple("P3G1", hana::int_c<3>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ),
                                             hana::make_tuple("P2G2", hana::int_c<2>, hana::int_c<2> ),
                                             hana::make_tuple("P3G2", hana::int_c<3>, hana::int_c<2> ) );

    int status = -1;
    std::vector<std::string> combinations;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status,&combinations]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                        constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>( hana::at_c<1>(d) ))>::value;
                        combinations.push_back(std::to_string(_dim)+","+_discretization);
                        if ( dimension == _dim && discretization == _discretization )
                            status = runApplicationThermoElectric<_dim,_torder,_gorder>();
                    } );
    if( status == -1 )
    {
        Feel::cout << "Wrong combination of dimension( " << dimension << ") and discretization (" << discretization << "). Possible combination:" << std::endl;
        for( auto const& s : combinations )
            Feel::cout << "\t("<< s << ")" << std::endl;
    }
    return status;
}
