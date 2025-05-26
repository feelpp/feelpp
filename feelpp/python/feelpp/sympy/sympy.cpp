#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// C++ convenience: call the Python get_coefficients API
std::map<std::string,std::string>
python_get_coeffs(std::string mode,
                  std::map<std::string,std::string> const& kwargs)
{
    py::gil_scoped_acquire G;
    auto sympy = py::module_::import("feelpp.sympy.api");
    auto get_coeffs = sympy.attr("get_coefficients");
    // Note: all values in kwargs are strings; Python will sympify them
    py::dict py_kwargs;
    for(auto& [k,v] : kwargs)
        py_kwargs[k.c_str()] = v;
    // call: get_coefficients(mode, **py_kwargs)
    py::object out = get_coeffs(py::str(mode), py_kwargs);
    // convert returned dict[str,str] into C++ map
    return out.cast<std::map<std::string,std::string>>();
}

// then expose a module if you like, or just call python_get_coeffs from your app
PYBIND11_EMBEDDED_MODULE(feelpp_sympy, m) {
    m.def("get_coeffs", &python_get_coeffs,
          "Call the feelpp.sympy API");
}