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

#include <feel/feelcore/table.hpp>

#if 0
#include <iostream>
#include <variant>
#include <string>
#include <vector>
#include <sstream>
#include <optional>
#endif

#include <iomanip>
#include <algorithm>






// http://www.cplusplus.com/forum/general/256212/
namespace utility
{
// output stream manipulator to centre an item in a field.
// uses the field width and fill character specified by the stream.
// for the most part, obeys the format flags specified by the stream
template < typename T > struct centre
{
    explicit centre( const T& v ) : what(v) {}
    const T& what ;

    template < typename C, typename CT >
    friend std::basic_ostream<C,CT>& operator << ( std::basic_ostream<C,CT>& stm, centre c )
        {
            const auto width = stm.width() ;
            if( width < 2 ) return stm << c.what ;

            std::basic_ostringstream<C,CT> str_stm ;
            str_stm.flags( stm.flags() ) ;
            // TO DO: handle internal justification
            str_stm << std::left << c.what ;
            const auto string_rep = str_stm.str() ;

            if( string_rep.size() >= width ) return stm << string_rep ;

            const auto left = ( width - string_rep.size() ) / 2 ;
            const auto right = width - string_rep.size() - left ;
            const auto fill = stm.fill() ;

            str_stm.str({}) ;
            for( int i = 0 ; i < left ; ++i ) str_stm << fill ;
            str_stm << string_rep ;
            for( int i = 0 ; i < right ; ++i ) str_stm << fill ;

            stm.width(0) ;
            return stm << str_stm.str() ;
        }

    // TO DO: allow chaining std::put_money etc.
};
}




namespace Feel
{

#if 0
std::vector<std::string>
splitByLines( std::string const& input )
{
    std::istringstream istr(input);
    std::string to;

    std::vector<std::string> outputStringByLines;
    while ( std::getline(istr,to,'\n') )
    {
        outputStringByLines.push_back( to );
    }
    return outputStringByLines;
}
#endif
  

typename Table::Cell::Format
Table::Cell::Format::newFromParent( Format const& parentFormat ) const
{
    Format newFormat;
    if ( M_fontAlign || parentFormat.M_fontAlign )
        newFormat.setFontAlign( M_fontAlign? *M_fontAlign : *(parentFormat.M_fontAlign) );
    if ( M_fontFloatingPoint || parentFormat.M_fontFloatingPoint )
        newFormat.setFloatingPoint( M_fontFloatingPoint? *M_fontFloatingPoint : *(parentFormat.M_fontFloatingPoint) );
    if ( M_paddingLeft || parentFormat.M_paddingLeft )
        newFormat.setPaddingLeft( M_paddingLeft? *M_paddingLeft : *(parentFormat.M_paddingLeft) );
    if ( M_paddingRight || parentFormat.M_paddingRight )
        newFormat.setPaddingRight( M_paddingRight? *M_paddingRight : *(parentFormat.M_paddingRight) );

    if ( M_widthMax || parentFormat.M_widthMax )
        newFormat.setWidthMax( M_widthMax? *M_widthMax : *(parentFormat.M_widthMax) );

    return newFormat;
}

void
Table::Cell::updateStream( std::ostream &o, std::string const& input, size_t width, Format const& format )
{
    //std::cout << "updateStream : "<< input << std::endl;
    int widthInput = width - format.paddingLeft() - format.paddingRight();
    if ( widthInput <= 0 )
        return;

    std::string paddingLeftStr(format.paddingLeft(),' ');
    std::string paddingRightStr(format.paddingRight(),' ');

    o << paddingLeftStr
      << std::setw(widthInput);
    switch ( format.fontAlign() )
    {
    case Font::Align::left:
        o << std::left << input;
        break;
    case Font::Align::right:
        o << std::right << input;
        break;
    case Font::Align::center:
        o << utility::centre(input);
        break;
    }

    o << paddingRightStr;
}


void
Table::add_row( std::initializer_list<variant_value_type> const& rowValues )
{
    if ( M_nRow > 0 && rowValues.size() != M_nCol )
        return;

    this->resize( M_nRow+1, rowValues.size() );
    int j = 0, i= M_nRow-1;
    for ( auto const& rowValue : rowValues )
        this->operator()(i,j++) = Cell(rowValue);
}

std::vector<std::string>
Table::toRawString( Cell::Format const& format ) const
{
    auto const& c = *this;

    int nRow = this->nRow();
    int nCol = this->nCol();
    std::vector<std::vector<std::string>> allOutputStringByLine(nRow*nCol);
    std::vector<Table::Cell::Format> allNewFormat(nRow*nCol);
    for (int i=0;i<nRow;++i)
    {
        for (int j=0;j<nCol;++j)
        {
            allNewFormat[i*nCol+j] = c(i,j).format().newFromParent( this->format() );
            allOutputStringByLine[i*nCol+j] = c(i,j).toRawString( allNewFormat[i*nCol+j] );
        }
    }

    std::vector<size_t> columnsWidth(nCol,0);
    std::vector<size_t> rowsHeight(nRow,0);
    for (int j=0;j<nCol;++j)
        for (int i=0;i<nRow;++i)
        {
            auto const& outputStringByLine = allOutputStringByLine[i*nCol+j];
            if ( outputStringByLine.empty() )
                continue;

            auto const& theFormat = allNewFormat[i*nCol+j];
            int paddingWidth = theFormat.paddingLeft() + theFormat.paddingRight();

            columnsWidth[j] = std::max( columnsWidth[j],
                                        std::max_element(outputStringByLine.begin(),outputStringByLine.end(),
                                                         [](std::string const& a, std::string const& b) { return a.size() < b.size(); } )->size() + paddingWidth );

            rowsHeight[i] = std::max( rowsHeight[i], outputStringByLine.size() );
        }

    // resize row heigh
    for (int i=0;i<nRow;++i)
    {
        int rowHeight = rowsHeight[i];
        for (int j=0;j<nCol;++j)
        {
            allOutputStringByLine[i*nCol+j].resize( rowHeight, std::string{} );
        }
    }


    std::string const& rowSep = this->format().rowSeparator();
    std::string const& colSep = this->format().columnSeparator();
    bool addLeftBorder = this->format().showLeftBorder();
    bool addRightBorder = this->format().showRightBorder();

    std::string horizontalLine, horizontalLineHeaderSep;
    std::string horizontalLineColSep( colSep.size(), '+');
    std::string headerColumnSep = "||";
    std::string horizontalLineFistColSep( this->format().firstColumnIsHeader()? headerColumnSep.size() : colSep.size(), '+');
    if ( addLeftBorder )
    {
        horizontalLine += "+";
        horizontalLineHeaderSep += "+";
    }
    for (int k=0;k<columnsWidth.size();++k )
    {
        if ( k == 1 && nCol > 1 )
        {
            horizontalLine += horizontalLineFistColSep;
            horizontalLineHeaderSep += horizontalLineFistColSep;

        }
        else if ( k>0 )
        {
            horizontalLine += horizontalLineColSep;
            horizontalLineHeaderSep += horizontalLineColSep;
        }
        int colWidth = columnsWidth[k];
        horizontalLine += std::string(colWidth,'-');
        horizontalLineHeaderSep += std::string(colWidth,'=');
    }
    if ( addRightBorder )
    {
        horizontalLine += "+";
        horizontalLineHeaderSep += "+";
    }

    std::vector<std::string> output;
    if ( this->format().showTopBorder() )
        output.push_back( horizontalLine );
    for (int i=0;i<nRow;++i)
    {
        if ( nCol > 1 && i == 1 && this->format().firstRowIsHeader() )
        {
            output.push_back( horizontalLineHeaderSep );
        }
        else if ( i>0 && this->format().hasRowSeparator() )
            output.push_back( horizontalLine );

        for ( int i2 = 0; i2 <rowsHeight[i] ;++i2 )
        {
            std::ostringstream o;
            if ( addLeftBorder )
                o << "|";
            for (int j=0;j<nCol;++j)
            {
                std::string const& cijStr = allOutputStringByLine[i*nCol+j][i2];
                if ( nCol > 1 && j==1 && this->format().firstColumnIsHeader() )
                    o << headerColumnSep;
                else if ( j>0 && this->format().hasColumnSeparator() )
                    o << colSep;

                Table::Cell::updateStream( o, cijStr, columnsWidth[j], allNewFormat[i*nCol+j] );
            }
            if ( addRightBorder )
                o << "|";
            output.push_back( o.str() );
        }
    }
    if ( this->format().showBottomBorder() )
        output.push_back( horizontalLine );
    return output;
}


void
Table::exportAsciiDoc( std::ostream &o ) const
{
    o<< "\n";
    o << "[cols=\"";
    if ( this->format().firstColumnIsHeader() )
    {
        o << "h";
        if ( this->nCol() > 1 )
            o << "," << this->nCol()-1 << "*";
    }
    else
        o << this->nCol();
    o << "\"";

    if ( this->format().firstRowIsHeader() )
        o << ",options=\"header\"";
    o << "]" << "\n";
    o << "|===" << "\n";

    for (int i=0;i<this->nRow();++i)
    {
        for (int j=0;j<this->nCol();++j)
        {
            auto const& cell = this->operator()(i,j);
            auto format = cell.format().newFromParent( this->format() );

            Table::updateOutputStreamOfCellUsingAsciiDoc(o,cell,format);
            o << "\n";
        }
        o << "\n";
    }
    o << "|===" << "\n";
}

void
Table::updateOutputStreamOfCellUsingAsciiDoc( std::ostream &o, Cell const& c, Cell::Format const& format )
{
    switch ( format.fontAlign() )
    {
    case Font::Align::left:
        o << "<";
        break;
    case Font::Align::right:
        o << ">";
        break;
    case Font::Align::center:
        o << "^";
        break;
    }

    o << "|";
    auto strByLine = c.toRawString( format );
    for ( std::string const& s : strByLine )
        o << s;
}

std::vector<std::string>
Table::Cell::toRawString( Format const& format ) const
{
    std::ostringstream ostr;
    if ( std::holds_alternative<Feel::TableBase>( this->M_value ) )
    {
        std::cout << "TODO" << std::endl;
    }
    else if ( std::holds_alternative<double>( this->M_value ) )
    {
        double val  = std::get<double>(this->M_value);
        switch ( this->format().floatingPoint() )
        {
        case Font::FloatingPoint::fixed:
            ostr << std::fixed << val;
            break;
        case Font::FloatingPoint::scientific:
            ostr << std::scientific << val;
            break;
        case Font::FloatingPoint::hexfloat:
            ostr << std::hexfloat << val;
            break;
        case Font::FloatingPoint::defaultfloat:
            ostr << std::defaultfloat << val;
            break;
        };
    }
    else
        std::visit([&ostr](auto && val) { ostr << val; }, this->M_value );

    // auto res = splitByLines( ostr.str() );
    // return res;
    //return ostr.str();

    int witdhMax = format.widthMax();
#if 1
    std::istringstream istr(ostr.str());
    std::string to;

    std::vector<std::string> outputStringByLines;
    while ( std::getline(istr,to,'\n') )
    {
        if ( witdhMax < 0 ) // not enable
            outputStringByLines.push_back( to );
        else
        {
            size_t currentWidth = to.size();
            size_t pos = 0;
            while ( pos < currentWidth )
            {
                int widthUsed = std::min( currentWidth-pos, (size_t)witdhMax);
                outputStringByLines.push_back( to.substr( pos,widthUsed ) );
                pos += widthUsed;
            }
        }
    }
    return outputStringByLines;
#endif
}

std::ostream&
operator<<(std::ostream& o, TableBase const& cb )
{
    auto const& c = static_cast<Table const&>( cb );

    for ( std::string const& s : c.toRawString( c.format() ) )
        o << s << "\n";
    return o;
}

std::ostream&
operator<<( std::ostream& o, Table::Cell const& c )
{
    for ( std::string const& s : c.toRawString( c.format() ) )
        o << s << "\n";
    return o;
}


} // namespace Feel
