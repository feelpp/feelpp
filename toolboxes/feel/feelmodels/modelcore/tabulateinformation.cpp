/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/tabulateinformation.hpp>
#include <feel/feelmodels/modelcore/terminalsize.hpp>

#include <tabulate/asciidoc_exporter.hpp>
//#include <tabulate/markdown_exporter.hpp>


namespace Feel
{

TabulateInformationProperties::TabulateInformationProperties( VerboseLevel vl )
    :
    M_information_level( 0 ),
    M_verbose_level( vl )
{}

typename TabulateInformations::self_ptrtype TabulateInformations::New( tabulate::Table const& table )
{
    return std::make_shared<TabulateInformationsTable>( table );
}

std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti )
{
    return o << ti.exporterAscii();
}

std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAscii const& ea )
{
    for ( std::string const& s : ea.M_outputStringByLines )
        o << s << "\n";
    return o;
}

std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAsciiDoc const& ead )
{
    ead.M_tabulateInformations->exportAsciiDoc( o, ead.M_startLevelSection );
    return o;
}


std::vector<std::string>
TabulateInformationsTable::exportAscii() const
{
    std::ostringstream ostr;
    ostr << M_table;
    std::istringstream istr(ostr.str());
    std::string to;

    std::vector<std::string> outputStringByLines;
    while ( std::getline(istr,to,'\n') )
    {
        outputStringByLines.push_back( to );
    }
    return outputStringByLines;
}

void
TabulateInformationsTable::exportAsciiDoc( std::ostream &o, int levelSection ) const
{
    tabulate::AsciiDocExporter exporter;
    //tabulate::MarkdownExporter exporter;
    auto asciidoc = exporter.dump(const_cast<tabulate::Table&>(M_table));
    o << asciidoc;
}

std::vector<std::string>
TabulateInformationsSections::exportAscii() const
{
    auto gethorizontalLine = [] ( size_t maxLineSize, bool addLeftBorder, bool addRightBorder ) {
        std::string horizontalLine;
        if ( addLeftBorder )
            horizontalLine += "+-";
        horizontalLine += std::string(maxLineSize,'-');
        if ( addRightBorder )
            horizontalLine += "-+";
        return horizontalLine;
    };

    auto getLineInContext = []( std::string const& s, size_t maxLineSize, bool addLeftBorder, bool addRightBorder ) {

        std::string sModif;
        if ( s.size() < maxLineSize )
        {
            std::ostringstream ostr;
            if ( false )//if ( s.size() < (maxLineSize-1) )
            {
                ostr << std::setw(maxLineSize-1);
                sModif += " ";
            }
            else
                ostr << std::setw(maxLineSize);
            ostr  << std::left << s;
            sModif += ostr.str();
        }
        std::string const& sUsed = ( s.size() < maxLineSize )? sModif: s;

        std::string res;
        if ( addLeftBorder )
            res += "| ";
        res += sUsed;
        if ( addRightBorder )
            res += " |";
        return res;
    };

    auto outputStringWidth = []( std::vector<std::string> const& osbl ) {
        size_t res = 0;
        for (auto const& l : osbl )
            res = std::max( res, l.size() );
        return res;
    };

    std::vector<std::string> outputStringByLines;

    //bool hasDoFirstSection = false;
    for ( auto const& [name,st] : M_subTab )
    {
        //auto outputStringByLinesCurrent = st->exportAscii();
        auto exporterAsciiCurrent = st->exporterAscii();
        auto const& outputStringByLinesCurrent = exporterAsciiCurrent.outputStringByLines();
        if ( outputStringByLinesCurrent.empty() )
            continue;

        size_t maxLineSize = std::max( outputStringWidth( outputStringByLinesCurrent ), name.size() );
        bool addLeftBorder = !name.empty(), addRightBorder = addLeftBorder;
        std::string horizontalLine = gethorizontalLine( maxLineSize, addLeftBorder, addRightBorder );
        if ( !name.empty() )
        {
            outputStringByLines.push_back( horizontalLine );
            outputStringByLines.push_back( getLineInContext(name,maxLineSize,addLeftBorder, addRightBorder) );
            outputStringByLines.push_back( horizontalLine );
        }
        for ( std::string const& s : outputStringByLinesCurrent )
            outputStringByLines.push_back( getLineInContext(s,maxLineSize,addLeftBorder, addRightBorder) );
        if ( !name.empty() )
            outputStringByLines.push_back( horizontalLine );

        //hasDoFirstSection = true;
    }
    return outputStringByLines;
}


void
TabulateInformationsSections::exportAsciiDoc( std::ostream &o, int levelSection ) const
{
    // TODO
    std::string levelSectionStr = std::string(levelSection+1,'=');
    for ( auto const& [name,st] : M_subTab )
    {
        int currentLevelSection = levelSection;
        if ( !name.empty() )
        {
            o << levelSectionStr << " " << name << "\n";
            ++currentLevelSection;
        }
        o <<  st->exporterAsciiDoc( currentLevelSection ) << "\n";
    }
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
        if ( TabulateInformationProperties::terminalProperties().hasWidth() )
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
        if ( TabulateInformationProperties::terminalProperties().hasWidth() )
        {
            max_size_col0 = std::max( table[nrow-1][0].size(), max_size_col0 );
            max_size_col1 = std::max( table[nrow-1][1].size(), max_size_col1 );
        }
    }
    //table.column(1).format().border_left(":");
    if ( TabulateInformationProperties::terminalProperties().hasWidth() )
    {
        int table_width = max_size_col0+max_size_col1+1;// table.shape().first;
        int table_width_max = (int)std::floor( (3./4)*TabulateInformationProperties::terminalProperties().width());
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

tabulate_informations_ptr_t
tabulateInformationsFunctionSpace( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New();

    tabulate::Table tabInfoOthers;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );
    tabInfo->add( "", TabulateInformations::New( tabInfoOthers ) );

    tabulate::Table tabInfoBasisEntries;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoBasisEntries, jsonInfo.at("basis"), tabInfoProp );
    tabInfo->add( "Basis", TabulateInformations::New( tabInfoBasisEntries ) );

    auto tabInfoDofTable = TabulateInformationsSections::New();
    auto const& jsonInfoDofTable = jsonInfo.at("doftable");
    tabulate::Table tabInfoDofTableEntries;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoDofTableEntries, jsonInfoDofTable, tabInfoProp );
    tabInfoDofTable->add( "", TabulateInformations::New( tabInfoDofTableEntries ) );
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
        tabInfoDofTable->add( "", TabulateInformations::New( tabInfoDofTableEntriesDatByPartition ) );
    }
    tabInfo->add( "Dof Table", tabInfoDofTable );

    return tabInfo;
}

} // namespace FromJSON

} // namespace TabulateInformationTools
} //namespace Feel

