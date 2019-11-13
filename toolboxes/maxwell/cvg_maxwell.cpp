#include <feel/feelmodels/maxwell/maxwell.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

void printCvg( std::ostream& out, std::vector<double> h, std::vector<std::map<std::string,std::map<std::string,double> > > errors )
{
    boost::format fmter("%1% %|14t|");

    out << fmter % "h";
    for( auto const& fp : errors[0] )
        for( auto const& np : errors[0][fp.first] )
            out << fmter % ("err"+fp.first+"_"+np.first);
    out << "\n";

    int maxiter = h.size();
    for( int i = 0; i < maxiter; ++i )
    {
        out << fmter % h[i];
        for( auto const& fp : errors[i] )
            for( auto const& np : errors[i][fp.first] )
                out << fmter % np.second;
        out << "\n";
    }
}

auto bExprDim(std::string ex, hana::int_<2> )
{
    return expr<15>(ex);
}

auto bExprDim(std::string ex, hana::int_<3> )
{
    return expr<3,1,15>(ex);
}

template<int Dim>
int
runApplicationMaxwell()
{
    using m_type = Maxwell<Simplex<Dim, 1> >;
    using mesh_type = typename m_type::mesh_type;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    std::shared_ptr<m_type> M( new m_type("maxwell") );
    auto model = M->modelProperties().models().model();
    auto paramValues = M->modelProperties().parameters().toParameterValues();

    auto aString = model.ptree().template get<std::string>("a_expr");
    auto aExpr = expr<Dim,1,15>(aString);
    aExpr.setParameterValues(paramValues);
    auto bString = model.ptree().get("b_expr","");
    constexpr int exprDim = Dim == 3 ? 3 : 1;
    auto bExpr = bExprDim(bString, hana::int_<Dim>() );
    bExpr.setParameterValues(paramValues);

    int maxiter = ioption("benchmark.nlevels");
    double factor = doption("benchmark.refine");
    double size = doption("benchmark.hsize");
    int quadError = ioption("benchmark.quad");

    exporter_ptrtype e( exporter_type::New("convergence") );

    std::vector<double> h(maxiter);
    std::vector<std::map<std::string,std::map<std::string,double> > > errors(maxiter, std::map<std::string,std::map<std::string,double> >());

    for( int i = 0; i < maxiter; ++i )
    {
        h[i] = size;
        auto mesh = loadMesh( _mesh=new mesh_type, _h=size );
        if( e->doExport() )
            e->step(i)->setMesh(mesh);

        M->setMesh(mesh);
        M->init();
        M->printAndSaveInfo();
        M->solve();

        auto rangeElements = M->rangeMeshElements();

        auto A_h = M->fieldMagneticPotential();
        auto Ah = M->spaceMagneticPotential();
        auto A_ex = Ah->element(aExpr);
        errors[i]["A"]["L2"] = normL2(_range=rangeElements,
                                      _expr=idv(A_h)-aExpr,
                                      _quad=quadError, _quad1=quadError );
        errors[i]["A"]["semiHcurl"] = normL2(_range=rangeElements,
                                             _expr=curlv(A_h)-bExpr,
                                             _quad=quadError, _quad1=quadError);
        errors[i]["A"]["Hcurl"] = errors[i]["A"]["L2"] + errors[i]["A"]["semiHcurl"];
        if( e->doExport() )
        {
            e->step(i)->add("A_h", A_h);
            e->step(i)->add("A_ex", A_ex);
        }

        auto B_h = M->fieldMagneticField();
        auto Bh = M->spaceMagneticField();
        auto B_ex = Bh->element(bExpr);
        errors[i]["B"]["L2"] = normL2(_range=rangeElements,
                                      _expr=idv(B_h)-bExpr,
                                      _quad=quadError, _quad1=quadError);
#if FEELPP_DIM==3
        errors[i]["B"]["semiHdiv"] = normL2(_range=rangeElements,
                                            _expr=divv(B_h),
                                            _quad=quadError, _quad1=quadError);
        errors[i]["B"]["Hdiv"] = errors[i]["B"]["semiHdiv"] + errors[i]["B"]["L2"];
#endif
        if( e->doExport() )
        {
            e->step(i)->add("B_h", B_h);
            e->step(i)->add("B_ex", B_ex);
        }

        if( e->doExport() )
            e->save();

        size /= factor;

        if( i == maxiter-1 )
        {
            double nAL2 = normL2( _range=rangeElements,
                                  _expr=aExpr,
                                  _quad=quadError, _quad1=quadError );
            double nBL2 = normL2( _range=rangeElements,
                                  _expr=bExpr,
                                  _quad=quadError, _quad1=quadError );
            for( int j = 0; j < maxiter; ++j)
            {
                errors[j]["A"]["L2_rel"] = errors[j]["A"]["L2"]/nAL2;
                errors[j]["A"]["Hcurl_rel"] = errors[j]["A"]["Hcurl"]/(nAL2+nBL2);
                errors[j]["B"]["L2_rel"] = errors[j]["B"]["L2"]/nBL2;
#if FEELPP_DIM==3
                errors[j]["B"]["Hdiv_rel"] = errors[j]["B"]["Hdiv"]/nBL2;
#endif
            }
        }
    }

    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        printCvg( std::cout, h, errors );
        std::ofstream file ( soption("benchmark.filename") );
        if( file )
        {
            printCvg( file, h, errors );
            file.close();
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description moptions( "maxwell options" );
    moptions.add( toolboxes_options("maxwell") );
    moptions.add_options()
        ("benchmark.filename", Feel::po::value<std::string>()->default_value( "cvg.csv" ), "file to save the results" )
        ("benchmark.quad", Feel::po::value<int>()->default_value( 5 ), "quadrature order for the error")
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=about(_name="feelpp_toolbox_maxwell",
                                        _author="Feel++ Consortium",
                                        _email="feelpp-devel@feelpp.org"),
                           _desc=moptions
                           );

    int dimension = ioption(_name="case.dimension");
    if( dimension == 2 )
        runApplicationMaxwell<2>();
    else
        runApplicationMaxwell<3>();

    return 0;
}
