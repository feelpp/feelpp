#pragma once

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <sstream>

namespace std
{
template<typename CharT>
struct formatter<GiNaC::ex, CharT>
{
    // We don’t support any format specifiers, so just skip parsing.
    constexpr auto parse(fmt::basic_format_parse_context<CharT>& ctx)
    {
        return ctx.begin();
    }
    template<typename FormatContext>
    auto format(GiNaC::ex const& e, FormatContext& ctx)
    {
        // Fallback to operator<< for GiNaC::ex
        std::basic_ostringstream<CharT> oss;
        oss << e;
        // Now format the resulting string
        return fmt::formatter<std::basic_string<CharT>, CharT>
                 ::format(oss.str(), ctx);
    }
};
} // namespace std
namespace fmt
{
template<typename CharT>
struct formatter<GiNaC::ex, CharT>
{
    // We don’t support any format specifiers, so just skip parsing.
    constexpr auto parse(fmt::basic_format_parse_context<CharT>& ctx)
    {
        return ctx.begin();
    }
    template<typename FormatContext>
    auto format(GiNaC::ex const& e, FormatContext& ctx)
    {
        // Fallback to operator<< for GiNaC::ex
        std::basic_ostringstream<CharT> oss;
        oss << e;
        // Now format the resulting string
        return fmt::formatter<std::basic_string<CharT>, CharT>
                 ::format(oss.str(), ctx);
    }
};
} // namespace std