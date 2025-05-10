#if !defined(FEELPP_FEEL_FEELALG_FMT_HPP)
#define FEELPP_FEEL_FEELALG_FMT_HPP 1

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <iostream>


FMT_BEGIN_NAMESPACE

template <typename Scalar, int Rank> struct formatter<Eigen::Tensor<Scalar, Rank>> : ostream_formatter {};

FMT_END_NAMESPACE



#if FMT_VERSION >= 100200

// Use the new nested formatter with fmt >= 10.2.0.
// This support nested Eigen types as well as padding/format specifiers.
#include <type_traits>

template <typename T>
struct fmt::formatter<T, std::enable_if_t<std::is_base_of<Eigen::DenseBase<T>, T>::value, char>>
    : fmt::nested_formatter<typename T::Scalar>
{
    auto format(T const& a, format_context& ctx) const
    {
        return this->write_padded(ctx, [&](auto out) {
            for (Eigen::Index ir = 0; ir < a.rows(); ir++) {
                for (Eigen::Index ic = 0; ic < a.cols(); ic++) {
                    out = fmt::format_to(out, "{} ", this->nested(a(ir, ic)));
                }
                if (ir + 1 < a.rows()) {
                    out = fmt::format_to(out, "\n");
                }
            }
            return out;
        });
    }
};

template <typename Derived>
struct fmt::is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};

#elif (FMT_VERSION >= 100000) || (FMT_VERSION >= 90000 && !defined(FMT_DEPRECATED_OSTREAM))

// fmt >= 10.x or fmt 9.x without deprecated ostream support.
#include <type_traits>

template <typename Derived>
struct fmt::formatter<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return m_underlying.parse(ctx);
    }

    template <typename FormatContext>
    auto format(const Derived& mat, FormatContext& ctx) const
    {
        auto out = ctx.out();

        for (Eigen::Index row = 0; row < mat.rows(); ++row) {
            for (Eigen::Index col = 0; col < mat.cols(); ++col) {
                out = fmt::format_to(out, "  ");
                out = m_underlying.format(mat.coeff(row, col), ctx);
            }

            if (row < mat.rows() - 1) {
                out = fmt::format_to(out, "\n");
            }
        }

        return out;
    }

private:
    fmt::formatter<typename Derived::Scalar, char> m_underlying;
};

template <typename Derived>
struct fmt::is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};

#else

// Include legacy ostr support

// clang-format off
#include <feel/feelcore/warnoff.h>
#include <fmt/ostream.h>
#include <feel/feelipre/warnon.h>
// clang-format on

template <typename Derived>
struct fmt::is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};

#endif


 #endif /* FEELPP_FEEL_FEELALG_FMT_HPP */