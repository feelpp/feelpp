/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#pragma once

#include <feel/feelcore/feel.hpp>
#include <tabulate/table.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel
{

class TabulateInformationProperties
{
public:
    enum class VerboseLevel { NONE,SMALL,MEDIUM,FULL };

    TabulateInformationProperties( VerboseLevel vl = VerboseLevel::FULL );

    static void updateTerminalSize( bool updateEvenIfAlreadyDefined = false );
    static int hasTerminalSize() { return S_terminalSize? true:false; }
    static int terminalWidth() { return S_terminalSize->first; }
    static int terminalHeight() { return S_terminalSize->second; }
private:
    uint16_type M_information_level;
    VerboseLevel M_verbose_level;
    inline static std::optional<std::pair<int,int>> S_terminalSize;
};

namespace TabulateInformationTools
{
namespace FromJSON
{

void
addKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::vector<std::string> const& keys );

inline
void
addKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::initializer_list<std::string> const& keys )
{
    addKeyToValues( table, jsonInfo, tabInfoProp, std::vector<std::string>(keys) );
}

void
addAllKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

tabulate::Table
tabulateFunctionSpace( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

} // namespace FromJSON

void
mergeSections( tabulate::Table & tabInfoMerge, std::vector<std::pair<std::string,tabulate::Table>> tabInfos, bool hasTitle = false );

tabulate::Table
createSections( std::vector<std::pair<std::string,tabulate::Table>> tabInfos, std::string const& title = "" );

} // namespace TabulateInformationTools

} // namespace Feel
