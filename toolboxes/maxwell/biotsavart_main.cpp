#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelmodels/maxwell/biotsavart.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "biotsavart",
                     "biotsavart",
                     "0.1",
                     "Coupled-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    return about;
}

int main(int argc, char** argv)
{
    auto opt = feel_options();
    opt.add_options()
        ("cond-marker", po::value<std::vector<std::string>>()->default_value({{"omega"}}),"")
        ("mag-marker", po::value<std::vector<std::string>>()->default_value({{"Box"}}),"")
        ("j-path", po::value<std::string>()->default_value("current-density"), "" )
        ("mag-dim", po::value<int>()->default_value(3),"")
        ("valid-x", po::value<double>()->default_value(0),"")
        ("valid-y", po::value<double>()->default_value(0),"")
        ("valid-z", po::value<double>()->default_value(0),"")
        ("compute-A", po::value<bool>()->default_value(true),"")
        ("compute-B", po::value<bool>()->default_value(true),"")
        ("compute-error", po::value<bool>()->default_value(false),"")
        ("box-radius", po::value<double>()->default_value(0.5), "")
        ("box-center", po::value<std::vector<double>>()->default_value({{0,0,0}}), "")
        ("export-bs", po::value<bool>()->default_value(true), "")
        ;
    Environment env( _argc=argc,
                     _argv=argv,
                     _about=makeAbout(),
                     _desc=opt.add(biotsavart_options())
                     );

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> >);

    // auto bs = BiotSavart<2>(mesh, GeoTool::Node(0,0,0), 2);
    // auto bs = BiotSavart<1>(mesh, GeoTool::Node(0,0,0), GeoTool::Node(0,0,1));
    // auto bs = BiotSavart<1>(mesh);
    auto center = vdoption("box-center");
    GeoTool::Node c(center[0],center[1],center[2]);
    auto radius = doption("box-radius");
    auto bs = BiotSavart<3>(mesh, c, radius);
    // // auto bs = BiotSavart<3>(mesh, GeoTool::Node(0,0,0), GeoTool::Node(0,0,1),2,4);

    bs.init();

    auto vecMarkersCond = vsoption("cond-marker");
    auto markersCond = std::set<std::string>(vecMarkersCond.begin(),vecMarkersCond.end());
    auto jEx = expr<3,1>(soption("functions.j"));
    // auto Xh = Pdhv<1>(mesh);//, markedelements(mesh,markersCond));
    // auto j = Xh->element();
    // j.load(soption("j-path"));
    // auto jEx = idv(j);

    bs.compute(jEx, boption("compute-B"), boption("compute-A"), markersCond);
    auto b = bs.magneticFieldValue(doption("valid-x"),doption("valid-y"),doption("valid-z"));
    auto a = bs.magneticPotentialValue(doption("valid-x"),doption("valid-y"),doption("valid-z"));
    Feel::cout << "b=\n" << b << "\na=\n" << a << std::endl;

    // auto bb = biotsavartComputeB(jEx,elements(mesh), doption("valid-x"),doption("valid-y"),doption("valid-z"),"mm");
    // Feel::cout << "b on valid:\n" << bb << std::endl;

    auto bEx = expr<3,1>(soption("functions.b"));
    bEx.setParameterValues(std::make_pair("z",doption("valid-z")));
    auto bbb = bEx.evaluate();
    Feel::cout << "b exact:\n" << bbb << std::endl;

    if( boption("export-bs") )
    {
        auto meshM = bs.mesh();
        auto B = bs.magneticField();
        auto A = bs.magneticPotential();
        auto e = exporter(_mesh=meshM);
        e->add("A",A);
        e->add("B",B);
        e->save();
    }
    return 0;
}
