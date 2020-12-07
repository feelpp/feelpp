/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/tabulateinformation.hpp>
#include <feel/feelmodels/modelcore/terminalsize.hpp>

namespace Feel
{

TabulateInformationProperties::TabulateInformationProperties( VerboseLevel vl )
    :
    M_information_level( 0 ),
    M_verbose_level( vl )
{
    this->updateTerminalSize();
}

void
TabulateInformationProperties::updateTerminalSize( bool updateEvenIfAlreadyDefined )
{
    if ( S_terminalSize && !updateEvenIfAlreadyDefined )
        return;
    S_terminalSize = get_terminal_size();
}

namespace TabulateInformationTools
{
namespace FromJSON
{
#if 0
bool is_array_of_primitive( nl::json const& jsonInfo )
{
    if ( !jsonInfo.is_array() )
        return false;
    for ( auto const& el : jsonInfo.items() )
        if ( !el.value().is_primitive() )
            return false;
    return true;
}

tabulate::Table tabulateArrayOfPrimitive( nl::json const& jsonInfo )
{
    //std::vector<std::string> arrayStringValues;
    std::vector<std::variant<std::string, tabulate::Table>> arrayStringValues; // will not compile on last tabulate
    for ( auto const& el : jsonInfo.items() )
    {
        auto const& arrayVal = el.value();
        if ( arrayVal.is_string() )
            arrayStringValues.push_back( arrayVal.get<std::string>() );
        else if ( arrayVal.is_number_integer() )
            arrayStringValues.push_back( std::to_string( arrayVal.get<int>() ) );
        else if ( arrayVal.is_number_float() )
            arrayStringValues.push_back( std::to_string( arrayVal.get<double>() ) );
        else
            CHECK( false ) << "not a primitive entry";
    }

    tabulate::Table table;
    table.add_row( arrayStringValues );
    table.format().hide_border_top().hide_border_bottom();
    return table;
}
#endif

void addKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::vector<std::string> const& keys )
{
    size_t nrow = 0; //std::distance(table.begin(),table.end());
    size_t max_size_col0 = 0, max_size_col1 = 0;
    for ( auto /*const*/& t : table )
    {
        CHECK( t.size() == 2 ) << "invalid table";
        ++nrow;
        if ( TabulateInformationProperties::hasTerminalSize() )
        {
            max_size_col0 = std::max( t[0].size(), max_size_col0 );
            max_size_col1 = std::max( t[1].size(), max_size_col1 );
        }
    }

    //for ( auto const& t : table ) ++nrow;
    for ( std::string const& key : keys )
    {
        if ( !jsonInfo.contains(key) )
            continue;
        auto const& valAtKey = jsonInfo.at(key);

        if ( !valAtKey.is_primitive() )
            continue;
        if ( valAtKey.is_string() )
            table.add_row( {key, valAtKey.get<std::string>()} );
        else if ( valAtKey.is_number_integer() )
            table.add_row( {key, std::to_string(valAtKey.get<int>())} );
        else if ( valAtKey.is_number_float() )
            table.add_row( {key, std::to_string(valAtKey.get<double>())} );
        else
            CHECK( false ) << "TODO missing primitive or not primitive";
        ++nrow;

        //if ( nrow > 1 )
        table[nrow-1].format().hide_border_top().hide_border_bottom();
        //table[nrow-1][0].format().column_separator(":");
        //table[nrow-1][1].format().column_separator(":");
        //table[nrow-1].format().hide_border_bottom();
        table[nrow-1][0].format().hide_border_left();
        table[nrow-1][1].format().hide_border_right();
        table[nrow-1][1].format().border_left(":");
        if ( TabulateInformationProperties::hasTerminalSize() )
        {
            max_size_col0 = std::max( table[nrow-1][0].size(), max_size_col0 );
            max_size_col1 = std::max( table[nrow-1][1].size(), max_size_col1 );
        }
    }
    //table.column(1).format().border_left(":");
    if ( TabulateInformationProperties::hasTerminalSize() )
    {
        int table_width = max_size_col0+max_size_col1+1;// table.shape().first;
        int table_width_max = (int)std::floor( (3./4)*TabulateInformationProperties::terminalWidth());
        if ( table_width > table_width_max )
        {
            table.column(0).format().width( max_size_col0 );
            table.column(1).format().width( table_width_max - (max_size_col0+1) );
        }
    }
}

void addAllKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    std::vector<std::string> keys;
    for ( auto const& el : jsonInfo.items() )
        keys.push_back( el.key() );
    addKeyToValues(table,jsonInfo,tabInfoProp,keys );
}


tabulate::Table tabulateFunctionSpace(  nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    tabulate::Table tabInfo;

    tabulate::Table tabInfoOthers;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );
    tabInfo.add_row({tabInfoOthers});

    tabulate::Table tabInfoBasis;
    tabInfoBasis.add_row({"Basis"});
    tabulate::Table tabInfoBasisEntries;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoBasisEntries, jsonInfo.at("basis"), tabInfoProp );
    tabInfoBasis.add_row({tabInfoBasisEntries});
    tabInfo.add_row({tabInfoBasis});

    tabulate::Table tabInfoDofTable;
    tabInfoDofTable.add_row({"Dof Table"});
    tabulate::Table tabInfoDofTableEntries;
    auto const& jsonInfoDofTable = jsonInfo.at("doftable");
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoDofTableEntries, jsonInfoDofTable, tabInfoProp );
    tabInfoDofTable.add_row({tabInfoDofTableEntries});

    if ( jsonInfoDofTable.contains( "nLocalDofWithoutGhost" ) )
    {
        tabulate::Table tabInfoDofTableEntriesDatByPartition;
        tabInfoDofTableEntriesDatByPartition.add_row({"partition id","nLocalDofWithGhost","nLocalDofWithoutGhost","nLocalGhost"});
        auto jarray_nLocalDofWithGhost = jsonInfoDofTable.at("nLocalDofWithGhost").items();
        auto jarray_nLocalDofWithoutGhost = jsonInfoDofTable.at("nLocalDofWithoutGhost").items();
        auto jarray_nLocalGhost = jsonInfoDofTable.at("nLocalGhost").items();
        int procId = 0;
        for (auto it1=jarray_nLocalDofWithGhost.begin(), it2 = jarray_nLocalDofWithoutGhost.begin(), it3=jarray_nLocalGhost.begin();
             it1 != jarray_nLocalDofWithGhost.end(); ++it1,++it2,++it3,++procId )
            tabInfoDofTableEntriesDatByPartition.add_row({ std::to_string(procId), it1.value().get<std::string>(), it2.value().get<std::string>(), it3.value().get<std::string>() });
        tabInfoDofTable.add_row({tabInfoDofTableEntriesDatByPartition});
    }
    tabInfo.add_row({tabInfoDofTable});

    tabInfo.format().hide_border();

    return tabInfo;
}

} // namespace FromJSON


void
mergeSections( tabulate::Table & tabInfoMerge, std::vector<std::pair<std::string,tabulate::Table>> tabInfos, bool hasTitle )
{
    size_t rowStartSections = hasTitle? 1 : 0;
    size_t nrow = 0; //std::distance(table.begin(),table.end());
    for ( auto const& t : tabInfoMerge ) ++nrow;
    for ( auto const& [tableName,tableInfo] : tabInfos )
    {
        tabulate::Table tabInfoSection;
        tabInfoSection.add_row({tableName});
        tabInfoSection.add_row({tableInfo});
        tabInfoMerge.add_row({tabInfoSection});

        if ( nrow > rowStartSections )
            tabInfoMerge[nrow].format().hide_border_top();
        ++nrow;
    }
}
tabulate::Table
createSections( std::vector<std::pair<std::string,tabulate::Table>> tabInfos, std::string const& title )
{
    tabulate::Table tabInfo;
    bool hasTitle = !title.empty();
    if ( hasTitle )
        tabInfo.add_row( {title} );
    mergeSections( tabInfo, tabInfos, hasTitle );
    return tabInfo;
}


} // namespace TabulateInformationTools
} //namespace Feel

