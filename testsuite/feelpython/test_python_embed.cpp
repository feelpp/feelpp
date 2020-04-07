
#include <pybind11/embed.h> // everything needed for embedding
#include <sstream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;
using namespace py::literals;

int main(int argc, char** argv ) {
    using namespace Feel;
    Environment env(  _argc=argc, _argv=argv);

    auto locals = py::dict("f"_a="[x**2,y**2,z**2]");
    py::module::import("sys").attr("path").cast<py::list>().append(Feel::Environment::expand("$top_srcdir/feelpp/feel/feelpython/"));
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "f=[x**2,y**2,z**2];\n"
         << "laplacian_f=toginac(laplacian(f,[x,y,z]),[x,y,z]);\n"
         << "print('lap(f)=', laplacian_f)\n";
    
    py::exec(ostr.str().c_str(), py::globals(), locals);

}
