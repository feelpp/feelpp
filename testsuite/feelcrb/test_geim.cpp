#include <feel/feelcrb/geim.hpp>
//#include <geim.hpp>
#include <feel/feel.hpp>
#include <iostream>
#include <feel/feeldiscr/pch.hpp>



int main(int argc, char**argv)
{
    using namespace Feel;
    po::options_description geimoptions( "GEIM options" );

    Environment env( _argc=argc, _argv=argv,
                 _desc=geimoptions,
                 _about=about(_name="qs_test_eim",
                              _author="Feel++ Consortium",
                              _email="feelpp-devel@feelpp.org"));
    // *** GeimFunctionalModel ***
    //--- Space ---
    using MeshType = Mesh< Simplex<2> >;
    tic();
    auto mesh = unitSquare();
    //typedef decltype( mesh ) MeshType;
    toc("loadMesh");

    tic();
    auto Xh = Pch<2>( mesh );
    toc("Vh");
    GeimFunctionalModel< Pch_type<MeshType,2> > sigma( Xh );
    //--- radius ---
    std::cout << "rayon:" << sigma.radius()<< '\n';
    sigma.setRadius(2.5);
    std::cout << "rayon:" << sigma.radius()<< '\n';
    //--- center ---
    std::cout << "center:" << sigma.center()<< '\n';
    std::vector<double> x_m(2,1.0);
    sigma.setCenter(x_m);
    std::cout << "center:" << sigma.center()<< '\n';
    // --- Linear Form ---
    sigma.defineForm();
    std::cout << sigma.container() << '\n';
    // *** DictionnaryGeim ***
    DictionnaryGeim< Pch_type<MeshType,2> > dico( Xh );
    dico.addElement(Xh,sigma.center(),sigma.radius());
    std::vector<double> v1(2,0);
    dico.addElement(Xh,v1,0.5);
    std::cout << "size:"<<dico.size() << '\n';

    for ( int i = 0; i < dico.size(); i++)
    {
        std::cout << "form:"<< i <<"\n"<< dico[i].container() << '\n';
    }

    return 0;
}
