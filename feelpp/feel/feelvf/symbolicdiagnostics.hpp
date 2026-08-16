/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
       Date: 2026-06-07

  Copyright (C) 2026 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file symbolicdiagnostics.hpp
   \author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
   \date 2026-06-07
 */
#ifndef FEELPP_VF_SYMBOLICDIAGNOSTICS_HPP
#define FEELPP_VF_SYMBOLICDIAGNOSTICS_HPP 1

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelvf/concepts.hpp>

namespace Feel
{
namespace vf
{

/**
 * @brief Symbolic differentiation diagnostic category.
 */
enum class SymbolicDifferentiationDiagnosticKind
{
    None,
    MissingDependencyQuery,
    MissingDiffOperation,
    UnsupportedNativeFunctor,
    UnsupportedPiecewiseFunction,
    MissingDomainAssumption,
    Exception
};

/**
 * @brief Convert a symbolic differentiation diagnostic category to text.
 */
inline std::string
toString( SymbolicDifferentiationDiagnosticKind kind )
{
    switch ( kind )
    {
    case SymbolicDifferentiationDiagnosticKind::None:
        return "none";
    case SymbolicDifferentiationDiagnosticKind::MissingDependencyQuery:
        return "missing dependency query";
    case SymbolicDifferentiationDiagnosticKind::MissingDiffOperation:
        return "missing symbolic diff operation";
    case SymbolicDifferentiationDiagnosticKind::UnsupportedNativeFunctor:
        return "unsupported native functor";
    case SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction:
        return "unsupported piecewise function";
    case SymbolicDifferentiationDiagnosticKind::MissingDomainAssumption:
        return "missing domain assumption";
    case SymbolicDifferentiationDiagnosticKind::Exception:
        return "exception";
    }
    return "unknown";
}

/**
 * @brief Stored source-location information for symbolic diagnostics.
 */
struct SymbolicDiagnosticSourceLocation
{
    std::string file;
    std::uint_least32_t line = 0;
    std::uint_least32_t column = 0;
    std::string function;

    /**
     * @brief Build a lightweight copy from std::source_location.
     */
    static SymbolicDiagnosticSourceLocation
    from( std::source_location const& location )
    {
        return { location.file_name(), location.line(), location.column(), location.function_name() };
    }

    /**
     * @brief Return a compact file:line description.
     */
    std::string
    summary() const
    {
        std::ostringstream os;
        os << file << ":" << line;
        if ( column != 0 )
            os << ":" << column;
        if ( !function.empty() )
            os << " in " << function;
        return os.str();
    }
};

/**
 * @brief Exception used when symbolic differentiation is known to be unsupported.
 */
class SymbolicDifferentiationError : public std::logic_error
{
public:
    SymbolicDifferentiationError( SymbolicDifferentiationDiagnosticKind kind,
                                  std::string message,
                                  std::source_location const& location = std::source_location::current() )
        :
        std::logic_error( message + " [source: " + SymbolicDiagnosticSourceLocation::from( location ).summary() + "]" ),
        M_kind( kind ),
        M_message( std::move( message ) ),
        M_location( SymbolicDiagnosticSourceLocation::from( location ) )
    {}

    /**
     * @brief Return the diagnostic category.
     */
    SymbolicDifferentiationDiagnosticKind kind() const noexcept { return M_kind; }

    /**
     * @brief Return the message without appended source information.
     */
    std::string const& message() const noexcept { return M_message; }

    /**
     * @brief Return the source location attached to the diagnostic.
     */
    SymbolicDiagnosticSourceLocation const& location() const noexcept { return M_location; }

private:
    SymbolicDifferentiationDiagnosticKind M_kind = SymbolicDifferentiationDiagnosticKind::Exception;
    std::string M_message;
    SymbolicDiagnosticSourceLocation M_location;
};

/**
 * @brief Report symbolic dependencies against an explicit candidate set.
 */
struct SymbolicDependencyReport
{
    bool supportsDependencyQuery = false;
    std::vector<std::string> checkedSymbols;
    std::vector<std::string> dependentSymbols;

    /**
     * @brief Return true if at least one checked symbol is a dependency.
     */
    bool hasDependency() const { return !dependentSymbols.empty(); }

    /**
     * @brief Return true if the expression depends on a checked symbol.
     */
    bool dependsOn( std::string const& symbol ) const
    {
        return std::find( dependentSymbols.begin(), dependentSymbols.end(), symbol ) != dependentSymbols.end();
    }

    /**
     * @brief Return a compact text summary.
     */
    std::string summary() const
    {
        std::ostringstream os;
        os << "supportsDependencyQuery=" << ( supportsDependencyQuery ? "true" : "false" )
           << ", dependentSymbols=[";
        for ( size_t i = 0; i < dependentSymbols.size(); ++i )
        {
            if ( i )
                os << ",";
            os << dependentSymbols[i];
        }
        os << "]";
        return os.str();
    }
};

/**
 * @brief Report whether an expression can be differentiated with respect to a symbol.
 */
struct SymbolicDifferentiabilityReport
{
    std::string symbol;
    bool supportsDependencyQuery = false;
    bool hasDependency = false;
    bool supportsSymbolicDiff = false;
    bool differentiable = false;
    SymbolicDifferentiationDiagnosticKind kind = SymbolicDifferentiationDiagnosticKind::None;
    std::string message;

    /**
     * @brief Return true if symbolic differentiation succeeded.
     */
    bool ok() const { return differentiable && kind == SymbolicDifferentiationDiagnosticKind::None; }

    /**
     * @brief Return a compact text summary.
     */
    std::string summary() const
    {
        std::ostringstream os;
        os << "symbol=" << symbol
           << ", supportsDependencyQuery=" << ( supportsDependencyQuery ? "true" : "false" )
           << ", hasDependency=" << ( hasDependency ? "true" : "false" )
           << ", supportsSymbolicDiff=" << ( supportsSymbolicDiff ? "true" : "false" )
           << ", differentiable=" << ( differentiable ? "true" : "false" )
           << ", kind=" << toString( kind );
        if ( !message.empty() )
            os << ", message=" << message;
        return os.str();
    }
};

/**
 * @brief Describe symbolic dependencies against an explicit candidate list.
 */
template <Feel::VfExpr ExprT>
SymbolicDependencyReport
describeSymbolicDependencies( ExprT const& expr, std::vector<std::string> const& symbols )
{
    SymbolicDependencyReport report;
    report.checkedSymbols = symbols;
    if constexpr ( Feel::HasSymbolDependency<ExprT> )
    {
        report.supportsDependencyQuery = true;
        for ( auto const& symbol : symbols )
        {
            if ( expr.hasSymbolDependency( symbol ) )
                report.dependentSymbols.push_back( symbol );
        }
    }
    return report;
}

/**
 * @brief Describe symbolic dependencies against an explicit candidate list.
 */
template <Feel::VfExpr ExprT>
SymbolicDependencyReport
describeSymbolicDependencies( ExprT const& expr, std::initializer_list<std::string> symbols )
{
    return describeSymbolicDependencies( expr, std::vector<std::string>( symbols ) );
}

/**
 * @brief Describe whether an expression differentiates with respect to a symbol.
 *
 * @details The function intentionally attempts the derivative when the public
 * `diff<1>(symbol)` API is available. This makes runtime unsupported paths,
 * such as native math functors without chain rules, visible to users and tests.
 */
template <Feel::VfExpr ExprT>
SymbolicDifferentiabilityReport
describeDifferentiability( ExprT const& expr, std::string const& symbol )
{
    SymbolicDifferentiabilityReport report;
    report.symbol = symbol;

    if constexpr ( Feel::HasSymbolDependency<ExprT> )
    {
        report.supportsDependencyQuery = true;
        report.hasDependency = expr.hasSymbolDependency( symbol );
    }
    else
    {
        report.kind = SymbolicDifferentiationDiagnosticKind::MissingDependencyQuery;
        report.message = "expression does not provide hasSymbolDependency(symbol)";
    }

    if constexpr ( Feel::SymbolicallyDifferentiableExpr<ExprT> )
    {
        report.supportsSymbolicDiff = true;
        try
        {
            auto derivative = expr.template diff<1>( symbol );
            (void)derivative;
            report.differentiable = true;
            report.kind = SymbolicDifferentiationDiagnosticKind::None;
            report.message.clear();
        }
        catch ( SymbolicDifferentiationError const& e )
        {
            report.differentiable = false;
            report.kind = e.kind();
            report.message = e.what();
        }
        catch ( std::exception const& e )
        {
            report.differentiable = false;
            report.kind = SymbolicDifferentiationDiagnosticKind::Exception;
            report.message = e.what();
        }
    }
    else
    {
        report.supportsSymbolicDiff = false;
        report.differentiable = false;
        report.kind = SymbolicDifferentiationDiagnosticKind::MissingDiffOperation;
        report.message = "expression does not provide diff<1>(symbol)";
    }

    return report;
}

namespace details
{
/**
 * @brief Build an unsupported native functor symbolic differentiation error.
 */
inline SymbolicDifferentiationError
unsupportedSymbolicDifferentiation( std::string const& functorKind, std::string const& functorName,
                                    std::source_location const& location = std::source_location::current() )
{
    return SymbolicDifferentiationError(
        SymbolicDifferentiationDiagnosticKind::UnsupportedNativeFunctor,
        std::string( "Feel::vf " ) + functorKind + " math functor '" + functorName +
            "' does not implement symbolic differentiation yet; use a GiNaC expr(...) law or add the native chain rule",
        location );
}

/**
 * @brief Build an unsupported piecewise or non-smooth function differentiation error.
 */
inline SymbolicDifferentiationError
unsupportedPiecewiseDifferentiation( std::string const& functionKind, std::string const& functionName,
                                     std::source_location const& location = std::source_location::current() )
{
    return SymbolicDifferentiationError(
        SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction,
        std::string( "Feel::vf " ) + functionKind + " function '" + functionName +
            "' is piecewise, non-smooth, or discontinuous; symbolic differentiation is disabled by default until an explicit piecewise AD policy is selected",
        location );
}

/**
 * @brief Build a missing domain-assumption differentiation error.
 */
inline SymbolicDifferentiationError
missingDomainAssumptionForPower( std::string const& message,
                                 std::source_location const& location = std::source_location::current() )
{
    return SymbolicDifferentiationError(
        SymbolicDifferentiationDiagnosticKind::MissingDomainAssumption,
        std::string( "Feel::vf pow differentiation requires an explicit domain assumption: " ) + message,
        location );
}
} // namespace details

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_SYMBOLICDIAGNOSTICS_HPP */
