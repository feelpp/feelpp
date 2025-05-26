#include <boost/ut.hpp>
#include <feelpp/sympy/api.hpp>

using namespace boost::ut;
namespace sym = feelpp::sympy::api;

int main() {
  suite sympy_api = [] {
    "Evaluate simple addition"_test = [] {
      auto result = sym::evaluate("1 + 2");
      expect(result == "3"_s);
    };

    "Evaluate symbolic derivative"_test = [] {
      auto deriv = sym::derivative("sin(x)", "x");
      // Expect the string "cos(x)" or equivalent
      expect(deriv == "cos(x)"_s);
    };

    "Matrix creation and determinant"_test = [] {
      // Suppose sym::matrix takes rows, cols and a list of entries
      sym::Matrix M = sym::matrix(2, 2, { "1", "2", "3", "4" });
      auto det = sym::determinant(M);
      expect(det == "-2"_s);
    };

    "Invalid expression throws"_test = [] {
      try {
        sym::evaluate("foobar@@");
        expect(false);  // Should never reach here
      } catch (const sym::SympyError& e) {
        expect(std::string{e.what()}.find("SyntaxError") != std::string::npos);
      }
    };
  };

  return ut::cfg<>.run({ sympy_api });
}