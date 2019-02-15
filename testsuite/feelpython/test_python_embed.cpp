
#include <pybind11/embed.h> // everything needed for embedding
#include <sstream>

#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;
using namespace py::literals;

int main() {
    py::scoped_interpreter guard{};

    auto locals = py::dict("f"_a="[x**2,y**2,z**2]");
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "f=[x**2,y**2,z**2];\n"
         << "laplacian_f=toginac(laplacian(f,[x,y,z]),[x,y,z]);\n"
         << "print('lap(f)=', laplacian_f)\n";
    
    py::exec(ostr.str().c_str(), py::globals(), locals);

}
