/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 Oct 2024

 Copyright (C) 2024 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#if !defined(FEELPP_FEEL_FEELALG_FMT_HPP)
#define FEELPP_FEEL_FEELALG_FMT_HPP 1

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#if defined(FEELPP_HAS_EIGEN3)
#include <Eigen/Core>
#include <Eigen/Dense>
#if __has_include(<Eigen/CXX11/Tensor>)
#include <Eigen/CXX11/Tensor>
#endif
#endif

FMT_BEGIN_NAMESPACE

#if defined(FEELPP_HAS_EIGEN3)

#if __has_include(<Eigen/CXX11/Tensor>)
// Formatter for Eigen::Tensor
template <typename Scalar, int Rank> 
struct formatter<Eigen::Tensor<Scalar, Rank>> : ostream_formatter {};
#endif

#if FMT_VERSION >= 100200
// Use the new nested formatter with fmt >= 10.2.0.
// This supports nested Eigen types as well as padding/format specifiers.
#include <type_traits>

template <typename T>
struct formatter<T, std::enable_if_t<std::is_base_of<Eigen::DenseBase<T>, T>::value, char>>
    : nested_formatter<typename T::Scalar>
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
struct is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};

#elif (FMT_VERSION >= 100000) || (FMT_VERSION >= 90000 && !defined(FMT_DEPRECATED_OSTREAM))
// fmt >= 10.x or fmt 9.x without deprecated ostream support.
#include <type_traits>

template <typename Derived>
struct formatter<
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
    formatter<typename Derived::Scalar, char> m_underlying;
};

template <typename Derived>
struct is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};

#else
// Include legacy ostream support
template <typename Derived>
struct is_range<
    Derived,
    std::enable_if_t<std::is_base_of<Eigen::DenseBase<Derived>, Derived>::value, char>>
    : std::false_type
{
};
#endif

#endif // FEELPP_HAS_EIGEN3

FMT_END_NAMESPACE

#endif /* FEELPP_FEEL_FEELALG_FMT_HPP */
