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

#include <iomanip>
#include <algorithm>

#include <feel/feelcore/termcolor.hpp>
#include <feel/feelcore/traits.hpp>




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
  


Table::Table( int nRow, int nCol, Table::MemoryLayout ml )
    :
    M_memoryLayout( ml ),
    M_nRow(nRow), M_nCol(nCol),
    M_format( new TableImpl::Format{} )
{
    M_cells.resize( nRow*nCol );
}

Table::Table( Table const& t )
    :
    M_memoryLayout( t.M_memoryLayout ),
    M_nRow( t.M_nRow ),
    M_nCol( t.M_nCol ),
    M_cells( t.M_cells ),
    M_format( new TableImpl::Format{t.format()} )
{}

void
Table::resize( int nRow, int nCol )
{
    if ( (nRow == M_nRow) && (nCol == M_nCol) )
        return;

    std::vector<TableImpl::Cell> cellsSave = std::move( M_cells );
    int nRowSave = M_nRow, nColSave = M_nCol;
    M_cells.resize( nRow*nCol );
    M_nRow = nRow;
    M_nCol = nCol;
    for (int i=0;i<std::min(M_nRow,nRowSave);++i)
    {
        for (int j=0;j<std::min(M_nCol,nColSave);++j)
        {
            int id = M_memoryLayout == MemoryLayout::RowMajor? i*nColSave+j : i+j*nRowSave;
            this->operator()(i,j) = std::move( cellsSave[id] );
        }
    }
}


std::vector<Printer::OutputText>
Table::toOutputText( TableImpl::Format const& format ) const
{
    auto const& c = *this;

    int nRow = this->nRow();
    int nCol = this->nCol();
    if ( nRow*nCol == 0 )
        return std::vector<Printer::OutputText>{};

    std::vector<std::vector<Printer::OutputText>> allOutputStringByLine(nRow*nCol);
    std::vector<TableImpl::Cell::Format> allNewFormat(nRow*nCol);
    for (int i=0;i<nRow;++i)
    {
        for (int j=0;j<nCol;++j)
        {
            allNewFormat[i*nCol+j] = c(i,j).format().newFromParent( format );
            allOutputStringByLine[i*nCol+j] = c(i,j).toOutputText( allNewFormat[i*nCol+j] );
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
                                                         [](auto const& a, auto const& b) { return a.size() < b.size(); } )->size() + paddingWidth );

            rowsHeight[i] = std::max( rowsHeight[i], outputStringByLine.size() );
        }

    // resize row heigh
    for (int i=0;i<nRow;++i)
    {
        int rowHeight = rowsHeight[i];
        for (int j=0;j<nCol;++j)
        {
            allOutputStringByLine[i*nCol+j].resize( rowHeight, Printer::OutputText{} );
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

    std::vector<Printer::OutputText> output;
    if ( this->format().showTopBorder() )
        output.push_back( Printer::OutputText{horizontalLine} );
    for (int i=0;i<nRow;++i)
    {
        if ( nCol > 1 && i == 1 && this->format().firstRowIsHeader() )
        {
            output.push_back( Printer::OutputText{horizontalLineHeaderSep} );
        }
        else if ( i>0 && this->format().hasRowSeparator() )
            output.push_back( Printer::OutputText{horizontalLine} );

        for ( int i2 = 0; i2 <rowsHeight[i] ;++i2 )
        {
            Printer::OutputText o;
            if ( addLeftBorder )
                o << "|";
            for (int j=0;j<nCol;++j)
            {
                if ( nCol > 1 && j==1 && this->format().firstColumnIsHeader() )
                    o << headerColumnSep;
                else if ( j>0 && this->format().hasColumnSeparator() )
                    o << colSep;

                TableImpl::Cell::updateWidthAndHorizontalAlign( allOutputStringByLine[i*nCol+j][i2], columnsWidth[j], allNewFormat[i*nCol+j] );
                o << allOutputStringByLine[i*nCol+j][i2];
            }
            if ( addRightBorder )
                o << "|";
            output.push_back( o );
        }
    }
    if ( this->format().showBottomBorder() )
        output.push_back( Printer::OutputText{horizontalLine} );
    return output;
}


void
Table::exportAsciiDocImpl( std::ostream &o, std::string const& tableSeparator, int nestedTableLevel ) const
{
    if ( this->empty() )
        return;
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
    o << tableSeparator << "===" << "\n";

    for (int i=0;i<this->nRow();++i)
    {
        for (int j=0;j<this->nCol();++j)
        {
            auto const& cell = this->operator()(i,j);
            auto format = cell.format().newFromParent( this->format() );
            TableImpl::updateOutputStreamOfCellUsingAsciiDoc(o,cell,format,tableSeparator,nestedTableLevel);
            o << "\n";
        }
        if ( i < (this->nRow()-1) )
            o << "\n";
    }
    o << tableSeparator << "===" << "\n";
}

void
Table::exportCSV( std::ostream &o ) const
{
    if ( this->empty() )
        return;

    for (int i=0;i<this->nRow();++i)
    {
        for (int j=0;j<this->nCol();++j)
        {
            auto const& c = this->operator()(i,j);
            auto cellFormatUsed = c.format().newFromParent( this->format() );
            auto ots = c.toOutputText( cellFormatUsed );
            if ( j > 0 )
                o << " , ";
            for ( auto const& ot : ots )
                o << ot;

        }
        o << "\n";
    }
}


std::ostream&
operator<<(std::ostream& o, Printer::OutputText const& cb )
{
    //namespace tc = termcolor;
    for ( auto const& d : cb.data() )
    {
        std::string const& txt = d.first;

        Font::Color color = std::get<0>(d.second);
        bool hasColor = color != Font::Color::none;
        if ( hasColor )
        {
            switch ( color )
            {
            case Font::Color::grey: o << termcolor::grey; break;
            case Font::Color::red: o << termcolor::red; break;
            case Font::Color::green: o << termcolor::green; break;
            case Font::Color::yellow: o << termcolor::yellow; break;
            case Font::Color::blue: o << termcolor::blue; break;
            case Font::Color::magenta: o << termcolor::magenta; break;
            case Font::Color::cyan: o << termcolor::cyan; break;
            case Font::Color::white: o << termcolor::white; break;
            default:
                break;
            }
        }
        o << txt;

        if ( hasColor )
            o << termcolor::reset;
    }
    return o;
}

std::ostream&
operator<<(std::ostream& o, Table const& c )
{
    for ( auto const& s : c.toOutputText( /*c.format()*/ ) )
        o << s << "\n";
    return o;
}

std::ostream&
operator<<( std::ostream& o, TableImpl::Cell const& c )
{
    for ( auto const& ot : c.toOutputText( c.format() ) )
        o << ot << "\n";
    return o;
}

namespace TableImpl
{

typename Cell::Format
Cell::Format::newFromParent( Format const& parentFormat ) const
{
    Format newFormat;
    if ( M_fontAlign || parentFormat.M_fontAlign )
        newFormat.setFontAlign( M_fontAlign? *M_fontAlign : *(parentFormat.M_fontAlign) );
    if ( M_fontColor || parentFormat.M_fontColor )
        newFormat.setFontColor( M_fontColor? *M_fontColor : *(parentFormat.M_fontColor) );
    if ( M_fontFloatingPoint || parentFormat.M_fontFloatingPoint )
        newFormat.setFloatingPoint( M_fontFloatingPoint? *M_fontFloatingPoint : *(parentFormat.M_fontFloatingPoint) );
    if ( M_fontFloatingPointPrecision || parentFormat.M_fontFloatingPointPrecision )
        newFormat.setFloatingPointPrecision( M_fontFloatingPointPrecision? *M_fontFloatingPointPrecision : *(parentFormat.M_fontFloatingPointPrecision) );
    if ( M_paddingLeft || parentFormat.M_paddingLeft )
        newFormat.setPaddingLeft( M_paddingLeft? *M_paddingLeft : *(parentFormat.M_paddingLeft) );
    if ( M_paddingRight || parentFormat.M_paddingRight )
        newFormat.setPaddingRight( M_paddingRight? *M_paddingRight : *(parentFormat.M_paddingRight) );

    if ( M_widthMax || parentFormat.M_widthMax )
        newFormat.setWidthMax( M_widthMax? *M_widthMax : *(parentFormat.M_widthMax) );

    return newFormat;
}

void
Cell::updateWidthAndHorizontalAlign( Printer::OutputText & ot, size_t width, Format const& format )
{
    int widthInput = width - format.paddingLeft() - format.paddingRight();
    if ( widthInput <= 0 )
        return;

    std::string paddingLeftStr(format.paddingLeft(),' ');
    std::string paddingRightStr(format.paddingRight(),' ');

    if ( ot.data().size() == 0 )
    {
        ot.push_back( std::string( width,' ' ) );
    }
    else if ( ot.data().size() == 1 )
    {
        std::string const& input = ot.data().front().first;
        std::ostringstream o;
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

        ot.data().front().first = o.str();
    }
    else
    {
        int widthMissing = widthInput - ot.size();

        switch ( format.fontAlign() )
        {
        case Font::Align::left:
        {
            std::string missingSpace( widthMissing, ' ' );
            ot.push_front( paddingLeftStr );
            ot.push_back( paddingRightStr + missingSpace );
            break;
        }
        case Font::Align::right:
        {
            std::string missingSpace( widthMissing, ' ' );
            ot.push_front( paddingLeftStr + missingSpace );
            ot.push_back( paddingRightStr  );
            break;
        }
        case Font::Align::center:
        {
            int widthMissingLeft = widthMissing/2;
            int widthMissingRight = widthMissing - widthMissingLeft;
            std::string missingSpaceLeft( widthMissingLeft, ' ' );
            std::string missingSpaceRight( widthMissingRight, ' ' );
            ot.push_front( paddingLeftStr + missingSpaceLeft );
            ot.push_back( paddingRightStr + missingSpaceRight );
            break;
        }
        }

    }
}

std::vector<Printer::OutputText>
Cell::toOutputText( Format const& format, bool enableWidthMax ) const
{
    std::ostringstream ostr;
    if ( std::holds_alternative<Feel::Table>( this->M_value ) )
    {
        auto const& c = std::get<Feel::Table>( this->M_value );
        return c.toOutputText( TableImpl::Format( c.format().newFromParent( format ) ) );
    }
    else if ( std::holds_alternative<Printer::OutputText>( this->M_value ) )
    {
        return std::vector<Printer::OutputText>(1,std::get<Printer::OutputText>( this->M_value ) );
    }
    else if ( std::holds_alternative<std::vector<Printer::OutputText>>( this->M_value ) )
    {
        return std::get<std::vector<Printer::OutputText>>( this->M_value );
    }
    else if ( std::holds_alternative<double>( this->M_value ) )
    {
        double val = std::get<double>(this->M_value);
        ostr << std::setprecision( format.floatingPointPrecision() );
        switch ( format.floatingPoint() )
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
        std::visit([&ostr](auto && val)
                   {
                       if constexpr ( !is_std_vector_v<std::decay_t<decltype(val)>> )
                                        ostr << val;
                   }, this->M_value );

    int witdhMax = enableWidthMax? format.widthMax() : -1;
#if 1
    std::istringstream istr(ostr.str());
    std::string to;

    std::vector<Printer::OutputText> outputStringByLines;

    while ( std::getline(istr,to,'\n') )
    {
        if ( witdhMax < 0 ) // not enable
            outputStringByLines.push_back( Printer::OutputText( to ) );
        else
        {
            size_t currentWidth = to.size();
            size_t pos = 0;
            while ( pos < currentWidth )
            {
                int widthUsed = std::min( currentWidth-pos, (size_t)witdhMax);
                outputStringByLines.push_back( Printer::OutputText( to.substr( pos,widthUsed ) ) );
                pos += widthUsed;
            }
        }
    }
    for ( auto & rs : outputStringByLines )
        rs.setColor( format.fontColor() );
    return outputStringByLines;
#endif
}

void updateOutputStreamOfCellUsingAsciiDoc( std::ostream &o, Cell const& c, Cell::Format const& format, std::string const& tableSeparator, int nestedTableLevel )
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

    if ( c.template is_a<Feel::Table>() && nestedTableLevel == 0 )
    {
        o << "a" << tableSeparator;
        c.template value<Feel::Table>().exportAsciiDocImpl( o, "!", nestedTableLevel+1 );
    }
    else
    {
        o << tableSeparator;
        for ( auto const& ot : c.toOutputText( format, false ) )
            for ( auto const& [s,p] : ot.data() )
                o << s;
    }
}


} // namespace TableImpl

} // namespace Feel
