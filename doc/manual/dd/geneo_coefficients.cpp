#include "geneo_coefficients.hpp"

namespace Feel
{
    double kappa::operator()(uint16_type, uint16_type, boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const&) const {
        return x[0] > 0.5 && x[1] > 0.5 ? val : 1.0;
    }

    double stripes::operator()(uint16_type, uint16_type, boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const&) const {
        return ((x[1] > 0.05 && x[1] < 0.2) || (x[1] > 0.3 && x[1] < 0.45) || (x[1] > 0.55 && x[1] < 0.7) || (x[1] > 0.8 && x[1] < 0.95)) ? first : second;
    }
}
