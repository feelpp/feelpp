#include <feel/feelpoly/context.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Feel
{
struct kappa
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    value_type val = 1;
    double operator()(uint16_type, uint16_type, boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const&) const;
};

struct stripes
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    value_type first  = 1;
    value_type second = 1;
    double operator()(uint16_type, uint16_type, boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const&) const;
};
}
