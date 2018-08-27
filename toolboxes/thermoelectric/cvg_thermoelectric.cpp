#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>
#include <feel/feelfilters/loadmesh.hpp>

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

template <int OrderT>
void
runApplicationThermoElectric()
{
    constexpr int geoOrder = OrderT < 2 ? OrderT : 2;
    using model_heat_type = FeelModels::Heat< Simplex<FEELPP_DIM,geoOrder>,
                                              Lagrange<OrderT, Scalar,Continuous,PointSetFekete> >;
    using model_electric_type = FeelModels::Electric< Simplex<FEELPP_DIM,geoOrder>,
                                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> >;
    using model_type = FeelModels::ThermoElectric< model_heat_type,model_electric_type>;
    using mesh_type = typename model_type::mesh_type;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    std::shared_ptr<model_type> thermoElectric( new model_type("thermo-electric", false) );

    int maxiter = ioption("cvg.nb-iter");
    double factor = doption("cvg.factor");
    double size = doption("gmsh.hsize");

    exporter_ptrtype e( exporter_type::New("convergence") );

    std::vector<double> h(maxiter);
    std::vector<double> normT(2*maxiter), errT(2*maxiter), errRelT(2*maxiter);
    std::vector<double> normV(2*maxiter), errV(2*maxiter), errRelV(2*maxiter);
    std::vector<double> normC(2*maxiter, 0), errC(2*maxiter, 0), errRelC(2*maxiter, 0);

    for( int i = 0; i < maxiter; ++i )
    {
        h[i] = size;
        auto mesh = loadMesh( _mesh=new mesh_type, _h=size );
        thermoElectric->setMesh(mesh);
        thermoElectric->init();
        thermoElectric->printAndSaveInfo();
        thermoElectric->solve();

        auto model = thermoElectric->modelProperties().models().model();
        e->step(i)->setMesh(mesh);

        auto heat = thermoElectric->heatModel();
        auto T_h = heat->fieldTemperature();
        auto tExpr = model.ptree().template get<std::string>("t_expr");
        auto Th = heat->spaceTemperature();
        auto T_ex = Th->element(expr<15>(tExpr));
        normT[2*i] = normL2( T_h.functionSpace()->template rangeElements<0>(),
                             idv(T_h), _quad=_Q<15>(), _quad1=_Q<15>() );
        errT[2*i] = normL2( T_h.functionSpace()->template rangeElements<0>(),
                            idv(T_h)-idv(T_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
        errRelT[2*i] = errT[2*i]/normT[2*i];
        normT[2*i+1] = normH1( T_h.functionSpace()->template rangeElements<0>(),
                               _expr=idv(T_h), _grad_expr=gradv(T_h), _quad=_Q<15>(), _quad1=_Q<15>() );
        errT[2*i+1] = normH1( T_h.functionSpace()->template rangeElements<0>(),
                              _expr=idv(T_h)-idv(T_ex), _grad_expr=gradv(T_h)-gradv(T_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
        errRelT[2*i+1] = errT[2*i+1]/normT[2*i+1];
        e->step(i)->add("T_h", T_h);
        e->step(i)->add("T_ex", T_ex);

        auto electric = thermoElectric->electricModel();
        auto V_h = electric->fieldElectricPotential();
        auto vExpr = model.ptree().template get<std::string>("v_expr");
        auto Vh = electric->spaceElectricPotential();
        auto V_ex = Vh->element(expr<15>(vExpr));
        normV[2*i] = normL2( V_h.functionSpace()->template rangeElements<0>(),
                           idv(V_h), _quad=_Q<15>(), _quad1=_Q<15>() );
        errV[2*i] = normL2( V_h.functionSpace()->template rangeElements<0>(),
                          idv(V_h)-idv(V_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
        errRelV[2*i] = errV[2*i]/normV[2*i];
        normV[2*i+1] = normH1( V_h.functionSpace()->template rangeElements<0>(),
                               _expr=idv(V_h), _grad_expr=gradv(V_h), _quad=_Q<15>(), _quad1=_Q<15>() );
        errV[2*i+1] = normH1( V_h.functionSpace()->template rangeElements<0>(),
                              _expr=idv(V_h)-idv(V_ex), _grad_expr=gradv(V_h)-gradv(V_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
        errRelV[2*i+1] = errV[2*i+1]/normV[2*i+1];
        e->step(i)->add("V_h", V_h);
        e->step(i)->add("V_ex", V_ex);

        auto C_h = electric->fieldElectricField();
        auto cExpr = model.ptree().get("c_expr", "");
        if( !cExpr.empty() )
        {
            auto Ch = electric->spaceElectricField();
            auto C_ex = Ch->element(expr<3,1,15>(cExpr));
            normC[2*i] = normL2( C_h.functionSpace()->template rangeElements<0>(),
                               idv(C_h), _quad=_Q<15>(), _quad1=_Q<15>() );
            errC[2*i] = normL2( C_h.functionSpace()->template rangeElements<0>(),
                              idv(C_h)-idv(C_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
            errRelC[2*i] = errC[2*i]/normC[2*i];
            normC[2*i+1] = normH1( C_h.functionSpace()->template rangeElements<0>(),
                                   _expr=idv(C_h), _grad_expr=gradv(C_h), _quad=_Q<15>(), _quad1=_Q<15>() );
            errC[2*i+1] = normH1( C_h.functionSpace()->template rangeElements<0>(),
                                  _expr=idv(C_h)-idv(C_ex), _grad_expr=gradv(C_h)-gradv(C_ex), _quad=_Q<15>(), _quad1=_Q<15>() );
            errRelC[2*i+1] = errC[2*i+1]/normC[2*i+1];
            e->step(i)->add("C_h", C_h);
            e->step(i)->add("C_ex", C_ex);
        }

        e->save();
        size /= factor;
    }

    std::vector<std::vector<double>> err = {errT, errRelT, errV, errRelV, errC, errRelC};
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        printCvg( std::cout, h, err );
        std::ofstream file ( soption("cvg.filename") );
        if( file )
        {
            printCvg( file, h, err );
            file.close();
        }
    }
}

int main(int argc, char** argv) {
    po::options_description thermoelectricoptions( "application thermo-electric options" );
    thermoelectricoptions.add( toolboxes_options("thermo-electric") );
    thermoelectricoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2,P3 ")
        ("cvg.nb-iter", po::value<int>()->default_value(1), "number of convergence iteration")
        ("cvg.factor", po::value<double>()->default_value(2), "factor of mesh size by iteration")
        ("cvg.filename", po::value<std::string>()->default_value("cvg.dat"), "filename to store convergence results")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=thermoelectricoptions,
                     _about=about(_name="cvg_thermoelectric",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationThermoElectric<1>();
#if FEELPP_INSTANTIATION_ORDER_MAX >= 2
    else if ( feapprox == "P2" )
        runApplicationThermoElectric<2>();
#endif
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    else if ( feapprox == "P3" )
        runApplicationThermoElectric<3>();
#endif
    else
        CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
