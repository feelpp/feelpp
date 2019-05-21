#include <pybind11/embed.h> // everything needed for embedding
#include <sstream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;
using namespace py::literals;

int main(int argc, char** argv ) {

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     //_desc=makeOptions(),
                     _about=about(_name="pyexpr",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "f=[x**2,y**2,z**2];\n"
         << "k=4;\n"
         << "kgrad_f=k*grad(f,[x,y,z]);\n"
         << "print('k*grad(f)=', kgrad_f)\n"
         << "div_kgrad_f=div(k*grad(f,[x,y,z]),[x,y,z]);\n"
         << "print('div(k*grad(f))=', div_kgrad_f)\n";

    auto _locals = py::dict("f"_a="[x**2,y**2,z**2]");
    std::map<std::string,std::string> locals;
    auto m = Feel::pyexpr( ostr.str(), {"kgrad_f","div_kgrad_f"}, locals );
}
