#define BOOST_TEST_MODULE python testsuite
#include <feel/feelcore/testsuite.hpp>
#include <fmt/core.h>
#include <pybind11/embed.h> // everything needed for embedding

#include <feel/feelcore/environment.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;
using namespace py::literals;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( python )

BOOST_AUTO_TEST_CASE( embed )
{
    using namespace Feel;
    auto locals = py::dict("f"_a="[x**2,y**2,z**2]");
    py::module::import("sys").attr("path").cast<py::list>().append(Feel::Environment::expand("$top_srcdir/feelpp/feel/feelpython/"));
    std::string code = fmt::format(
        "try:\n"
        "   from sympy2ginac import *\n"
        "except ImportError:\n"
        "   from feelpp.sympy2ginac import *\n"
        "f=[x**2,y**2,z**2];\n"
        "laplacian_f=toginac(laplacian(f,[x,y,z]),[x,y,z]);\n"
        "print('lap(f)=', laplacian_f)\n");
    
    py::exec(code.c_str(), py::globals(), locals);
    //py::print(locals["laplacian_f"]);
    //BOOST_CHECK( py::str(locals["laplacian_f"]) == std::string{"{2,2,2}:x:y:z"} );
}
BOOST_AUTO_TEST_SUITE_END()