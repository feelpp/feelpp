/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#pragma once

#include <feel/feelcore/feel.hpp>
#include <tabulate/table.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelmodels/modelcore/terminalsize.hpp>

namespace Feel
{

class TabulateInformationProperties
{
public:
    enum class VerboseLevel { NONE,SMALL,MEDIUM,FULL };

    TabulateInformationProperties( VerboseLevel vl = VerboseLevel::FULL );

    static TerminalProperties const& terminalProperties() { return S_terminalProperties; }
private:
    uint16_type M_information_level;
    VerboseLevel M_verbose_level;
    inline static TerminalProperties S_terminalProperties;
};

class TabulateInformations
{
public :
    using self_ptrtype = std::shared_ptr<TabulateInformations>;

    TabulateInformations() = default;
    TabulateInformations( TabulateInformations const& ) = default;
    TabulateInformations( TabulateInformations && ) = default;

    static self_ptrtype New( tabulate::Table const& table );

    virtual void updateForUse() = 0;

    friend std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti );

    std::vector<std::string> const& outputStringByLines() const { return M_outputStringByLines; }

    size_t outputStringWidth() const;
protected:
    //std::ostringstream M_ostr;
    std::vector<std::string> M_outputStringByLines;
};

std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti );

using tabulate_informations_ptr_t = typename TabulateInformations::self_ptrtype;

class TabulateInformationsTable : public TabulateInformations
{
public :
    explicit TabulateInformationsTable( tabulate::Table const& table ) : M_table( table ) {}
    TabulateInformationsTable( TabulateInformationsTable const& ) = default;
    TabulateInformationsTable( TabulateInformationsTable && ) = default;

    void updateForUse() override;

private:
    tabulate::Table M_table;
};

class TabulateInformationsSections : public TabulateInformations
{
public:
    TabulateInformationsSections() = default;
    TabulateInformationsSections( TabulateInformationsSections const& ) = default;
    TabulateInformationsSections( TabulateInformationsSections && ) = default;

    static std::shared_ptr<TabulateInformationsSections> New() { return std::make_shared<TabulateInformationsSections>(); }

    static std::shared_ptr<TabulateInformationsSections> cast( tabulate_informations_ptr_t t ) { return std::dynamic_pointer_cast<TabulateInformationsSections>(t); }

    void updateForUse() override;

    void add( std::string const& name, tabulate_informations_ptr_t t ) { M_subTab.push_back( std::make_pair( name, t ) ); }

    void erase( std::string const& name )
        {
            auto itFind = std::find_if( M_subTab.begin(), M_subTab.end(), [&name]( auto const& e ) { return name == e.first; } );
            if ( itFind != M_subTab.end() )
                M_subTab.erase( itFind );
        }
private:
    std::vector<std::pair<std::string,tabulate_informations_ptr_t>> M_subTab;
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

tabulate_informations_ptr_t
tabulateInformationsFunctionSpace( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

} // namespace FromJSON

void
mergeSections( tabulate::Table & tabInfoMerge, std::vector<std::pair<std::string,tabulate::Table>> tabInfos, bool hasTitle = false );

tabulate::Table
createSections( std::vector<std::pair<std::string,tabulate::Table>> tabInfos, std::string const& title = "" );

} // namespace TabulateInformationTools

} // namespace Feel
