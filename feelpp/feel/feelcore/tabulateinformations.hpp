/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#pragma once

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/table.hpp>
#include <tabulate/table.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/terminalproperties.hpp>

namespace Feel
{

class TabulateInformationProperties
{
public:
    TabulateInformationProperties( uint16_type vl = 1 );
    TabulateInformationProperties( TabulateInformationProperties const& ) = default;
    TabulateInformationProperties( TabulateInformationProperties && ) = default;

    uint16_type verboseLevel() const { return M_verboseLevel; }

    TabulateInformationProperties newByIncreasingVerboseLevel( uint16_type vl = 1 ) const
        {
            TabulateInformationProperties tip( *this );
            tip.M_verboseLevel += vl;
            return tip;
        }

    static TerminalProperties const& terminalProperties() { return S_terminalProperties; }
private:
    uint16_type M_verboseLevel;
    inline static TerminalProperties S_terminalProperties;
};

/**
 * @brief base class that describe informations by tabulate/section design
 */
class TabulateInformations
{
public :
    using self_ptrtype = std::shared_ptr<TabulateInformations>;

    struct ExporterAscii
    {
        ExporterAscii( std::vector<Printer::OutputText> && ot, int maxWidth = -1 ) : M_outputText( std::move( ot ) )
            {
                if ( maxWidth > 0 )
                    for ( auto & ot : M_outputText )
                        ot.applyMaxWidth( maxWidth );
            }
        ExporterAscii( ExporterAscii const& ) = default;
        ExporterAscii( ExporterAscii && ) = default;
        friend std::ostream& operator<<( std::ostream& o, ExporterAscii const& ea );

        std::vector<Printer::OutputText> const& outputText() const { return M_outputText; }
    private :
        std::vector<Printer::OutputText> M_outputText;
    };

     struct ExporterAsciiDoc
     {
         ExporterAsciiDoc( std::shared_ptr<const TabulateInformations> ti, int startLevelSection )
             :
             M_tabulateInformations( std::const_pointer_cast<TabulateInformations>( ti) ),
             M_startLevelSection( startLevelSection )
             {}
         ExporterAsciiDoc( ExporterAsciiDoc const& ) = default;
         ExporterAsciiDoc( ExporterAsciiDoc && ) = default;
         friend std::ostream& operator<<( std::ostream& o, ExporterAsciiDoc const& ea );
     private :
         self_ptrtype M_tabulateInformations;
         int M_startLevelSection;
     };

    TabulateInformations( uint16_type vl = 1 ) : M_verboseLevel( vl ) {}
    TabulateInformations( TabulateInformations const& ) = default;
    TabulateInformations( TabulateInformations && ) = default;

    virtual ~TabulateInformations() {}

    uint16_type verboseLevel() const { return M_verboseLevel; }

    //! get an exporter in ascii format
    ExporterAscii exporterAscii( bool useTerminalWidth = false ) const
        {
            int maxWidth = -1;
            if ( useTerminalWidth && TabulateInformationProperties::terminalProperties().hasWidth() && TabulateInformationProperties::terminalProperties().width() > 0 )
                maxWidth = TabulateInformationProperties::terminalProperties().width();
            return ExporterAscii( this->exportAscii(), maxWidth );
        }
    //! get an exporter in asciidoc format
    ExporterAsciiDoc exporterAsciiDoc( int startLevelSection = 0 ) const { return ExporterAsciiDoc( this->shared_from_this_tabulate_informations(), startLevelSection ); }

    //! create an new tabulate informations from a table
    static self_ptrtype New( Feel::Table const& table, TabulateInformationProperties const& tip = TabulateInformationProperties{} );
    //! create an new tabulate informations from a table
    static self_ptrtype New( tabulate::Table const& table, TabulateInformationProperties const& tip = TabulateInformationProperties{} );

    friend std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti );
    friend std::ostream& operator<<( std::ostream& o, ExporterAsciiDoc const& ea );
protected :
    virtual std::shared_ptr<const TabulateInformations> shared_from_this_tabulate_informations() const = 0;
    virtual std::vector<Printer::OutputText> exportAscii() const = 0;
    virtual void exportAsciiDoc( std::ostream &o, int levelSection ) const = 0;
    uint16_type M_verboseLevel;
};

//! ostream operator
std::ostream& operator<<( std::ostream& o, TabulateInformations const& ti );
std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAscii const& ti );
std::ostream& operator<<( std::ostream& o, TabulateInformations::ExporterAsciiDoc const& ti );

//! types
using tabulate_informations_t = TabulateInformations;
using tabulate_informations_ptr_t = typename TabulateInformations::self_ptrtype;

/**
 * @brief describe informations in table
 */
class TabulateInformationsTable : public TabulateInformations,
                                  public std::enable_shared_from_this<TabulateInformationsTable>
{
public :
    explicit TabulateInformationsTable( Feel::Table const& table, uint16_type vl = 1 ) : TabulateInformations( vl ), M_table( table ) {}
    TabulateInformationsTable( TabulateInformationsTable const& ) = default;
    TabulateInformationsTable( TabulateInformationsTable && ) = default;

private:
    std::shared_ptr<const TabulateInformations> shared_from_this_tabulate_informations() const override { return std::dynamic_pointer_cast<const TabulateInformations>( this->shared_from_this() );  }
    std::vector<Printer::OutputText> exportAscii() const override;
    void exportAsciiDoc( std::ostream &o, int levelSection ) const override;
private:
    Feel::Table M_table;
};

/**
 * @brief describe informations in table
 */
class TabulateInformationsTableAlternative : public TabulateInformations,
                                             public std::enable_shared_from_this<TabulateInformationsTableAlternative>
{
public :
    explicit TabulateInformationsTableAlternative( tabulate::Table const& table, uint16_type vl = 1 ) : TabulateInformations( vl ), M_table( table ) {}
    TabulateInformationsTableAlternative( TabulateInformationsTableAlternative const& ) = default;
    TabulateInformationsTableAlternative( TabulateInformationsTableAlternative && ) = default;

private:
    std::shared_ptr<const TabulateInformations> shared_from_this_tabulate_informations() const override { return std::dynamic_pointer_cast<const TabulateInformations>( this->shared_from_this() );  }
    std::vector<Printer::OutputText> exportAscii() const override;
    void exportAsciiDoc( std::ostream &o, int levelSection ) const override;
private:
    tabulate::Table M_table;
};

/**
 * @brief describe informations by section
 */
class TabulateInformationsSections : public TabulateInformations,
                                     public std::enable_shared_from_this<TabulateInformationsSections>
{
public:
    TabulateInformationsSections( uint16_type vl = 1 ) : TabulateInformations( vl ) {}
    TabulateInformationsSections( TabulateInformationsSections const& ) = default;
    TabulateInformationsSections( TabulateInformationsSections && ) = default;

    //! create a new tabulate informations of section
    static std::shared_ptr<TabulateInformationsSections> New( TabulateInformationProperties const& tip = TabulateInformationProperties{} ) { return std::make_shared<TabulateInformationsSections>( tip.verboseLevel() ); }

    //! cast from the base class
    static std::shared_ptr<TabulateInformationsSections> cast( tabulate_informations_ptr_t t ) { return std::dynamic_pointer_cast<TabulateInformationsSections>(t); }

    //! add a section
    void add( std::string const& name, tabulate_informations_ptr_t t ) { M_subTab.push_back( std::make_pair( name, t ) ); }

    //! erase a section from the name
    void erase( std::string const& name )
        {
            auto itFind = std::find_if( M_subTab.begin(), M_subTab.end(), [&name]( auto const& e ) { return name == e.first; } );
            if ( itFind != M_subTab.end() )
                M_subTab.erase( itFind );
        }
private:
    std::shared_ptr<const TabulateInformations> shared_from_this_tabulate_informations() const override { return std::dynamic_pointer_cast<const TabulateInformations>( this->shared_from_this() );  }
    std::vector<Printer::OutputText> exportAscii() const override;
    void exportAsciiDoc( std::ostream &o, int levelSection ) const override;

private:
    std::vector<std::pair<std::string,tabulate_informations_ptr_t>> M_subTab;
};

//! types
using tabulate_informations_table_t = TabulateInformationsTable;
using tabulate_informations_table_ptr_t = std::shared_ptr<tabulate_informations_table_t>;
using tabulate_informations_sections_t = TabulateInformationsSections;
using tabulate_informations_sections_ptr_t = std::shared_ptr<tabulate_informations_sections_t>;


namespace TabulateInformationTools
{
namespace FromJSON
{

void
addKeyToValues( Feel::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::vector<std::string> const& keys );

inline
void
addKeyToValues( Feel::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::initializer_list<std::string> const& keys )
{
    addKeyToValues( table, jsonInfo, tabInfoProp, std::vector<std::string>(keys) );
}

void
addAllKeyToValues( Feel::Table &table, nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

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


Feel::Table createTableFromArray( nl::json const& jsonInfo, bool applyDefaultFormat = false );

tabulate_informations_ptr_t
tabulateInformationsFunctionSpace( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

tabulate_informations_ptr_t
tabulateInformationsSymbolsExpr( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, bool addNameExpr = false );


} // namespace FromJSON

} // namespace TabulateInformationTools

} // namespace Feel
