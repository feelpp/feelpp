/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelcore/tabulateinformations.hpp>

#include <tabulate/asciidoc_exporter.hpp>
//#include <tabulate/markdown_exporter.hpp>


namespace Feel
{

TabulateInformationProperties::TabulateInformationProperties( uint16_type vl )
    :
    M_verboseLevel( vl )
{}

typename TabulateInformations::self_ptrtype TabulateInformations::New( Feel::Table const& table, TabulateInformationProperties const& tip )
{
    return std::make_shared<TabulateInformationsTable>( table, tip.verboseLevel() );
}
typename TabulateInformations::self_ptrtype TabulateInformations::New( tabulate::Table const& table, TabulateInformationProperties const& tip )
{
    return std::make_shared<TabulateInformationsTableAlternative>( table, tip.verboseLevel() );
}

std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti )
{
    return o << ti.exporterAscii();
}

std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAscii const& ea )
{
    for ( auto const& s : ea.outputText() )
        o << s << "\n";
    return o;
}

std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAsciiDoc const& ead )
{
    ead.M_tabulateInformations->exportAsciiDoc( o, ead.M_startLevelSection );
    return o;
}


std::vector<Printer::OutputText>
TabulateInformationsTable::exportAscii() const
{
    return M_table.toOutputText();
}

void
TabulateInformationsTable::exportAsciiDoc( std::ostream &o, int levelSection ) const
{
    M_table.exportAsciiDoc( o );
}

std::vector<Printer::OutputText>
TabulateInformationsTableAlternative::exportAscii() const
{
    std::ostringstream ostr;
    ostr << M_table;
    std::istringstream istr(ostr.str());
    std::string to;

    std::vector<Printer::OutputText> otext;
    while ( std::getline(istr,to,'\n') )
    {
        otext.push_back( Printer::OutputText(to) );
    }
    return otext;
}

void
TabulateInformationsTableAlternative::exportAsciiDoc( std::ostream &o, int levelSection ) const
{
    tabulate::AsciiDocExporter exporter;
    //tabulate::MarkdownExporter exporter;
    auto asciidoc = exporter.dump(const_cast<tabulate::Table&>(M_table));
    o << asciidoc;
}

std::vector<Printer::OutputText>
TabulateInformationsSections::exportAscii() const
{
    uint16_type maxVerboseLevel = 1;
    Feel::Table t;
    t.format().setShowAllBorders( false ).setHasRowSeparator( false ).setAllPadding( 0 );
    for ( auto const& [name,st] : M_subTab )
    {
        if ( maxVerboseLevel != invalid_v<uint16_type> && st->verboseLevel() > maxVerboseLevel )
            continue;
        auto exporterAsciiCurrent = st->exporterAscii();
        auto const& otCurrent = exporterAsciiCurrent.outputText();
        if ( otCurrent.empty() )
            continue;

        Feel::Table t2;
         if ( !name.empty() )
             t2.add_row( {name} );
         else
             t2.format().setShowAllBorders( false ).setHasRowSeparator( false ).setAllPadding( 0 );

        t2.add_row( { otCurrent } );

        t.add_row( { t2.toOutputText() } );
    }
    return t.toOutputText();
}


void
TabulateInformationsSections::exportAsciiDoc( std::ostream &o, int levelSection ) const
{
    std::string levelSectionStr = std::string(levelSection+1,'=');
    for ( auto const& [name,st] : M_subTab )
    {
        int currentLevelSection = levelSection;
        if ( !name.empty() )
        {
            o << "\n" << levelSectionStr << " " << name << "\n";
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

void addKeyToValues( Feel::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::vector<std::string> const& keys )
{
    for ( std::string const& key : keys )
    {
        if ( !jsonInfo.contains(key) )
            continue;
        auto const& valAtKey = jsonInfo.at(key);

        if ( !valAtKey.is_primitive() || valAtKey.is_null() )
            continue;
        if ( valAtKey.is_string() )
            table.add_row( {key, valAtKey.get<std::string>()} );
        else if ( valAtKey.is_boolean() )
            table.add_row( {key, valAtKey.get<bool>()} );
        else if ( valAtKey.is_number_integer() )
            table.add_row( {key, valAtKey.get<int>()} );
        else if ( valAtKey.is_number_float() )
            table.add_row( {key, valAtKey.get<double>()} );
        else
            CHECK( false ) << "TODO missing primitive or not primitive at key " << key ;
    }

    // TODO : to improve (not very nice and perfect, should be done in printer)
    if ( TabulateInformationProperties::terminalProperties().hasWidth() )
    {
        int table_width_max = (int)std::floor( (3./4)*TabulateInformationProperties::terminalProperties().width());
        table.format().setWidthMax( table_width_max );
    }
}


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
        else if ( valAtKey.is_boolean() )
            table.add_row( {key, std::to_string(valAtKey.get<bool>())} );
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

void addAllKeyToValues( Feel::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    std::vector<std::string> keys;
    for ( auto const& el : jsonInfo.items() )
        keys.push_back( el.key() );
    addKeyToValues(table,jsonInfo,tabInfoProp,keys );
}
void addAllKeyToValues( tabulate::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    std::vector<std::string> keys;
    for ( auto const& el : jsonInfo.items() )
        keys.push_back( el.key() );
    addKeyToValues(table,jsonInfo,tabInfoProp,keys );
}

Feel::Table createTableFromArray( nl::json const& jsonInfo, bool applyDefaultFormat )
{
    Feel::Table tabInfo;
    if ( !jsonInfo.is_array() )
        return tabInfo;

    for ( auto const& el : jsonInfo.items() )
    {
        auto const& arrayVal = el.value();
        if ( arrayVal.is_string() )
            tabInfo.add_row( jsonInfo.get<std::vector<std::string>>() );
        break;
    }

    if ( applyDefaultFormat )
        tabInfo.format()
            .setShowAllBorders( false )
            .setColumnSeparator(",")
            .setHasRowSeparator( false );
    return tabInfo;
}

tabulate_informations_ptr_t
tabulateInformationsFunctionSpace( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    Feel::Table tabInfoOthers;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );
    tabInfoOthers.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    tabInfo->add( "", TabulateInformations::New( tabInfoOthers,tabInfoProp ) );

    Feel::Table tabInfoBasisEntries;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoBasisEntries, jsonInfo.at("basis"), tabInfoProp );
    tabInfoBasisEntries.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    tabInfo->add( "Basis", TabulateInformations::New( tabInfoBasisEntries,tabInfoProp ) );

    auto tabInfoDofTable = TabulateInformationsSections::New( tabInfoProp );
    auto const& jsonInfoDofTable = jsonInfo.at("doftable");
    Feel::Table tabInfoDofTableEntries;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoDofTableEntries, jsonInfoDofTable, tabInfoProp );
    tabInfoDofTableEntries.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    tabInfoDofTable->add( "", TabulateInformations::New( tabInfoDofTableEntries,tabInfoProp ) );
    if ( jsonInfoDofTable.contains( "nLocalDofWithoutGhost" ) )
    {
        Feel::Table tabInfoDofTableEntriesDatByPartition;
        tabInfoDofTableEntriesDatByPartition.format().setFirstRowIsHeader( true );
        tabInfoDofTableEntriesDatByPartition.add_row({"partition id","nLocalDofWithGhost","nLocalDofWithoutGhost","nLocalGhost"});
        auto jarray_nLocalDofWithGhost = jsonInfoDofTable.at("nLocalDofWithGhost").items();
        auto jarray_nLocalDofWithoutGhost = jsonInfoDofTable.at("nLocalDofWithoutGhost").items();
        auto jarray_nLocalGhost = jsonInfoDofTable.at("nLocalGhost").items();
        int procId = 0;
        for (auto it1=jarray_nLocalDofWithGhost.begin(), it2 = jarray_nLocalDofWithoutGhost.begin(), it3=jarray_nLocalGhost.begin();
             it1 != jarray_nLocalDofWithGhost.end(); ++it1,++it2,++it3,++procId )
            tabInfoDofTableEntriesDatByPartition.add_row({ procId, it1.value().get<int>(), it2.value().get<int>(), it3.value().get<int>() });
        tabInfoDofTable->add( "", TabulateInformations::New( tabInfoDofTableEntriesDatByPartition, tabInfoProp.newByIncreasingVerboseLevel() ) );
    }
    tabInfo->add( "Dof Table", tabInfoDofTable );

    return tabInfo;
}

tabulate_informations_ptr_t
tabulateInformationsSymbolsExpr( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, bool addNameExpr )
{
    Feel::Table tabInfo;
    if ( !jsonInfo.is_array() )
        return TabulateInformations::New( tabInfo, tabInfoProp );

    for ( auto const& el : jsonInfo.items() )
    {
        auto const& jse = el.value();
        if ( !jse.contains( "symbol" ) )
            continue;

        std::string name = jse.value<std::string>( "name", "" );
        std::string expr = jse.value( "expr", "" );
        if ( tabInfo.nRow() == 0 )
        {
            if ( addNameExpr )
                tabInfo.add_row({"Name","Expression","Symbol","Shape","Components"});
            else
                tabInfo.add_row({"Symbol","Shape","Components"});
            tabInfo.format().setFirstRowIsHeader( true );
        }

        if ( jse.contains( "components" ) )
        {
            auto const& jcomps = jse.at( "components" );
            Feel::Table tabInfoComp;
            tabInfoComp.add_row({"Symbol","Indices"});
            tabInfoComp.format().setFirstRowIsHeader( true );
            for ( auto const& el2 : jcomps.items() )
            {
                auto const& jsec = el2.value();
                auto const& idx = jsec.at("indices");
                std::string idxStr = std::to_string( idx[0].get<int>() ) + "," + std::to_string( idx[1].get<int>() );
                tabInfoComp.add_row({ jsec.at("symbol").get<std::string>(), idxStr });
            }
            if ( addNameExpr )
                tabInfo.add_row({ name, expr, jse.at("symbol").get<std::string>(), jse.at("shape").get<std::string>(), tabInfoComp });
            else
                tabInfo.add_row({ jse.at("symbol").get<std::string>(), jse.at("shape").get<std::string>(), tabInfoComp });
        }
        else
        {
            if ( addNameExpr )
                tabInfo.add_row({ name, expr, jse.at("symbol").get<std::string>(), jse.at("shape").get<std::string>(), "" });
            else
                tabInfo.add_row({ jse.at("symbol").get<std::string>(), jse.at("shape").get<std::string>(), "" });
        }
    }
    return TabulateInformations::New( tabInfo, tabInfoProp );
}


} // namespace FromJSON

} // namespace TabulateInformationTools
} //namespace Feel

