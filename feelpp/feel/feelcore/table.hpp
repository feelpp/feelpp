/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 18 Dec 2020

 Copyright (C) 2020 Feel++ Consortium

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
#pragma once

#include <iostream>
#include <variant>
#include <string>
#include <vector>
#include <optional>
#include <tuple>

#include <numeric>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{

namespace Font
{
enum class Align { left, right, center };

enum class FloatingPoint { defaultfloat, fixed, scientific, hexfloat };

enum class Color { none, grey, red, green, yellow, blue, magenta, cyan, white };

}

namespace Printer
{
class OutputText
{
public :
    OutputText() = default;
    explicit OutputText( std::string const& s ) { this->push_back( s ); }
    OutputText( OutputText && ) = default;
    OutputText( OutputText const& ) = default;

    OutputText& operator=( OutputText && ) = default;
    OutputText& operator=( OutputText const& ) = default;

    void push_back( std::string const& s ) { M_data.push_back( std::make_pair( s, std::make_tuple( Font::Color::none ) ) ); }
    void push_back( OutputText const& ot )
        {
            for ( auto const& d : ot.data() )
                M_data.push_back( d );
        }
    void push_front( std::string const& s ) { M_data.insert( M_data.begin(), std::make_pair( s, std::make_tuple( Font::Color::none ) ) ); }
    size_t size() const { return std::accumulate( M_data.begin(), M_data.end(), 0, []( size_t res, auto const& e ) { return res + e.first.size(); } ); }

    std::vector<std::pair<std::string,std::tuple<Font::Color>> > const& data() const { return M_data; }
    std::vector<std::pair<std::string,std::tuple<Font::Color>> > & data() { return M_data; }

    OutputText& operator<<( std::string const& s ) { this->push_back( s ); return *this; }
    OutputText& operator<<( OutputText const& ot ) { this->push_back( ot ); return *this; }

    void setColor( Font::Color c ) { for ( auto & [s,opt] : M_data ) std::get<0>( opt ) = c; }
private:
    std::vector<std::pair<std::string,std::tuple<Font::Color>> > M_data;
};

}


namespace TableImpl
{
class Cell;
struct Format;
}

class Table
{
public :
    using variant_value_type = std::variant<std::string,int,size_t,double,Feel::Table,Printer::OutputText,std::vector<Printer::OutputText>>;

    Table();
    //Table( Table const& ) = default;
    Table( Table const& t );
    Table( Table && ) = default;
    Table& operator=( Table const& ) = default;
    Table& operator=( Table && ) = default;

    //! number of row in table
    int nRow() const { return M_nRow; }

    //! number of column in table
    int nCol() const { return M_nCol; }

    //! resize the table
    void resize( int nRow, int nCol )
        {
            CHECK( nRow*nCol > M_nRow*M_nCol) << "TODO case reduce size";
            M_cells.resize( nRow*nCol );
            M_nRow = nRow;
            M_nCol = nCol;
        }

    //! get cell at indices i,j in table
    TableImpl::Cell const& operator()(int i,int j) const { return M_cells[i*M_nCol+j]; }

    //! get cell at indices i,j in table
    TableImpl::Cell & operator()(int i,int j) { return M_cells[i*M_nCol+j]; }

    //! add row in table
    void add_row( std::initializer_list<variant_value_type> const& rowValues ) { this->add_row<std::initializer_list<variant_value_type>>( rowValues ); }

    template <typename T>
    void
    add_row( T const& rowValues, std::enable_if_t< is_iterable_v<T> >* = nullptr );

    //! get table format
    TableImpl::Format& format() { return *M_format; }
    //! get table format
    TableImpl::Format const& format() const { return *M_format; }

    //! return the output text of the table representation (each element of the vector start at a new line)
    std::vector<Printer::OutputText> toOutputText() const { return this->toOutputText( this->format() ); }

    //! return the output text of the table representation (each element of the vector start at a new line)
    //! the \format arg override the format of the table
    std::vector<Printer::OutputText> toOutputText( TableImpl::Format const& format ) const;

    //! export into stream \o the asciidoc representation of the table
    void exportAsciiDoc( std::ostream &o ) const { this->exportAsciiDocImpl( o ); }

    //friend void updateOutputStreamOfCellUsingAsciiDoc( std::ostream &o, TableImpl::Cell const& c, TableImpl::Cell::Format const& format, std::string const& tableSeparator, int nestedTableLevel = 0 );
    //private :
    void exportAsciiDocImpl( std::ostream &o, std::string const& tableSeparator = "|", int nestedTableLevel = 0 ) const;

private :
    int M_nRow, M_nCol;
    std::vector<TableImpl::Cell> M_cells;
    std::shared_ptr/*unique_ptr*/<TableImpl::Format> M_format;
};

std::ostream&
operator<<(std::ostream& o, Printer::OutputText const& cb );

std::ostream&
operator<<(std::ostream& o, Table const& cb );

std::ostream&
operator<<( std::ostream& o, TableImpl::Cell const& c );


namespace TableImpl
{
class Cell
{
public:

    /**
     * Format of a Cell in a Table
     */
    struct Format
    {
        Format() = default;
        Format( Format const& ) = default;
        Format( Format && ) = default;
        Format& operator=( Format const& ) = default;
        Format& operator=( Format && ) = default;

        Font::Align fontAlign() const { return M_fontAlign? *M_fontAlign : Font::Align::left; }
        Format& setFontAlign( Font::Align a ) { M_fontAlign = a; return *this; }

        Font::Color fontColor() const { return M_fontColor? *M_fontColor : Font::Color::none; }
        Format& setFontColor( Font::Color a ) { M_fontColor = a; return *this; }

        int paddingLeft() const { return M_paddingLeft? *M_paddingLeft : 1; }
        int paddingRight() const { return M_paddingRight? *M_paddingRight : 1; }
        Format& setPaddingLeft( int v ) { if ( v < 0 ) return *this; M_paddingLeft = v; return *this; }
        Format& setPaddingRight( int v ) { if ( v < 0 ) return *this; M_paddingRight = v; return *this; }
        Format& setAllPadding( int v ) { if ( v < 0 ) return *this; M_paddingLeft = v; M_paddingRight = v; return *this; }

        Font::FloatingPoint floatingPoint() const { return M_fontFloatingPoint? *M_fontFloatingPoint : Font::FloatingPoint::scientific; }
        Format& setFloatingPoint( Font::FloatingPoint a ) { M_fontFloatingPoint = a; return *this; }

        int widthMax() const { return M_widthMax? *M_widthMax : -1; }
        Format& setWidthMax( int i ) { M_widthMax = i; return *this; }

        Format newFromParent( Format const& parentFormat ) const;

    private :
        std::optional<Font::Align> M_fontAlign;
        std::optional<Font::Color> M_fontColor;
        std::optional<Font::FloatingPoint> M_fontFloatingPoint;
        std::optional<int> M_paddingLeft, M_paddingRight;

        std::optional<int> M_widthMax;
    };
    using variant_value_type = Table::variant_value_type;

    Cell() = default;
    explicit Cell( variant_value_type const& v ) : M_value( v ) {}
    Cell( Cell const& ) = default;
    Cell( Cell && ) = default;

    Cell& operator=( Cell const& ) = default;
    Cell& operator=( Cell && ) = default;

    Format& format() { return M_format; }
    Format const& format() const { return M_format; }

    template <typename T>
    bool is_a() const { return std::holds_alternative<T>( this->M_value ); }

    template <typename T>
    T const& value() const
        {
            if ( !this->is_a<T>() )
                CHECK(false) << "value type of cell is not compatible";
            return std::get<T>( this->M_value );
        }

    std::vector<Printer::OutputText> toOutputText( Format const& format, bool enableWidthMax = true ) const;

    static void updateWidthAndHorizontalAlign( Printer::OutputText & ot, size_t width, Format const& format );

private:
    friend std::ostream& operator<<(std::ostream& o, Cell const& c );

private:
    variant_value_type M_value;
    Format M_format;
};


/**
 * Format of a Table
 * inherits of Cell::Format which define the defaut format of cell (can be override in each cell
 */
struct Format : public Cell::Format
{
    Format() = default;
    Format( Format const& ) = default;
    Format( Format && ) = default;
    Format& operator=( Format const& ) = default;
    Format& operator=( Format && ) = default;

    Format( Cell::Format const& f )
        : Cell::Format( f ) {}

    bool hasRowSeparator() const { return M_hasRowSeparator; }
    Format& setHasRowSeparator( bool b ) { M_hasRowSeparator = b; return *this; }
    bool hasColumnSeparator() const { return M_hasColumnSeparator; }
    Format& setHasColumnSeparator( bool b ) { M_hasColumnSeparator = b; return *this; }

    std::string const& rowSeparator() const { return M_rowSeparator; }
    Format& setRowSeparator( std::string const& s ) { M_rowSeparator = s; return *this; }
    std::string const& columnSeparator() const { return M_columnSeparator; }
    Format& setColumnSeparator( std::string const& s ) { M_columnSeparator = s; return *this; }

    bool showLeftBorder() const { return M_showLeftBorder; }
    bool showRightBorder() const { return M_showRightBorder; }
    bool showTopBorder() const { return M_showTopBorder; }
    bool showBottomBorder() const { return M_showBottomBorder; }
    Format& setShowLeftBorder( bool b ) { M_showLeftBorder = b; return *this; }
    Format& setShowRightBorder( bool b ) { M_showRightBorder = b; return *this; }
    Format& setShowTopBorder( bool b ) { M_showTopBorder = b; return *this; }
    Format& setShowBottomBorder( bool b ) { M_showBottomBorder = b; return *this; }
    Format& setShowAllBorders( bool b ) { this->setShowLeftBorder( b ); this->setShowRightBorder( b ); this->setShowTopBorder( b ); this->setShowBottomBorder( b ); return *this; }


    bool firstRowIsHeader() const { return M_firstRowIsHeader; }
    Format& setFirstRowIsHeader( bool b ) { M_firstRowIsHeader = b; return *this; }
    bool firstColumnIsHeader() const { return M_firstColumnIsHeader; }
    Format& setFirstColumnIsHeader( bool b ) { M_firstColumnIsHeader = b; return *this; }

    // override Cell::Format (return the good type)
    Format& setFontAlign( Font::Align a ) { Cell::Format::setFontAlign( a ); return *this; }
    Format& setFontColor( Font::Color a ) { Cell::Format::setFontColor( a ); return *this; }
    Format& setPaddingLeft( int v ) { Cell::Format::setPaddingLeft( v ); return *this; }
    Format& setPaddingRight( int v ) { Cell::Format::setPaddingRight( v ); return *this; }
    Format& setAllPadding( int v ) { Cell::Format::setAllPadding( v ); return *this; }
    Format& setFloatingPoint( Font::FloatingPoint a ) { Cell::Format::setFloatingPoint( a ); return *this; }
    Format& setWidthMax( int i ) { Cell::Format::setWidthMax( i ); return *this; }
private :

    bool M_hasRowSeparator = true;
    bool M_hasColumnSeparator = true;
    std::string M_rowSeparator = "-";
    std::string M_columnSeparator = "|";
    bool M_showLeftBorder = true;
    bool M_showRightBorder = true;
    bool M_showTopBorder = true;
    bool M_showBottomBorder = true;
    bool M_firstRowIsHeader = false;
    bool M_firstColumnIsHeader = false;
};

void updateOutputStreamOfCellUsingAsciiDoc( std::ostream &o, Cell const& c, Cell::Format const& format, std::string const& tableSeparator, int nestedTableLevel = 0 );


} // namespace TableImpl


template <typename T>
void
Table::add_row( T const& rowValues, std::enable_if_t< is_iterable_v<T> >* )
{
    size_t containerSize = std::distance( rowValues.begin(), rowValues.end() ); // rowValues.size()
    //size_t containerSize = std::distance( std::begin(rowValues), std::end( rowValues ) ); // rowValues.size()
    if ( M_nRow > 0 && containerSize != M_nCol )
        return;

    this->resize( M_nRow+1, containerSize );
    int j = 0, i= M_nRow-1;
    for ( auto const& rowValue : rowValues )
        this->operator()(i,j++) = TableImpl::Cell(rowValue);
}



} // namespace Feel
