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

namespace Feel
{
struct TableBase
{
    virtual ~TableBase() {}
};

namespace Font
{
enum class Align { left, right, center };

enum class FloatingPoint { defaultfloat, fixed, scientific, hexfloat };
}

class Table : public TableBase
{
public :
    class Cell
    {
    public:

        struct Format
        {
            Format() = default;
            Format( Format const& ) = default;
            Format( Format && ) = default;
            Format& operator=( Format const& ) = default;
            Format& operator=( Format && ) = default;

            Font::Align fontAlign() const { return M_fontAlign? *M_fontAlign : Font::Align::left; }
            Format& setFontAlign( Font::Align a ) { M_fontAlign = a; return *this; }

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
            std::optional<Font::FloatingPoint> M_fontFloatingPoint;
            std::optional<int> M_paddingLeft, M_paddingRight;

            std::optional<int> M_widthMax;
        };
        using variant_value_type = std::variant<std::string,int,size_t,double,Feel::TableBase>;

        Cell() = default;
        explicit Cell( variant_value_type const& v ) : M_value( v ) {}
        Cell( Cell const& ) = default;
        Cell( Cell && ) = default;

        Cell& operator=( Cell const& ) = default;
        Cell& operator=( Cell && ) = default;

        Format& format() { return M_format; }
        Format const& format() const { return M_format; }

        std::vector<std::string> toRawString( Format const& format ) const;

        static void updateStream( std::ostream &o, std::string const& input, size_t width, Format const& format );

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

    Table() : M_nRow(0), M_nCol(0) {}
    Table( Table const& ) = default;
    Table( Table && ) = default;

    using variant_value_type = typename Cell::variant_value_type;

    //! number of row in table
    int nRow() const { return M_nRow; }

    //! number of column in table
    int nCol() const { return M_nCol; }

    //! resize the table
    void resize( int nRow, int nCol )
        {
            // TODO case reduce size
            M_cells.resize( nRow*nCol );
            M_nRow = nRow;
            M_nCol = nCol;
        }

    //! get cell at indices i,j in table
    Cell const& operator()(int i,int j) const { return M_cells[i*M_nCol+j]; }

    //! get cell at indices i,j in table
    Cell & operator()(int i,int j) { return M_cells[i*M_nCol+j]; }

    //! add row in table
    void add_row( std::initializer_list<variant_value_type> const& rowValues );

    //! get table format
    Format& format() { return M_format; }
    //! get table format
    Format const& format() const { return M_format; }

    void exportAsciiDoc( std::ostream &o ) const;

    std::vector<std::string> toRawString( Cell::Format const& format ) const;

private :
    int M_nRow, M_nCol;
    std::vector<Cell> M_cells;
    Format M_format;
};


std::ostream&
operator<<(std::ostream& o, TableBase const& cb );

std::ostream&
operator<<( std::ostream& o, Table::Cell const& c );

} // namespace Feel
