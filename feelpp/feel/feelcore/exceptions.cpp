/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2022-02-20

 Copyright (C) 2022 Feel++ Consortium

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
#include <boost/core/demangle.hpp>
#include <fmt/core.h>
#include <fmt/color.h>
#include <feel/feelcore/environment.hpp>

namespace Feel
{
void
printGitReport()
{
    if ( GitMetadata::populated() )
    {
        if ( GitMetadata::anyUncommittedChanges() )
        {
            std::cerr << fmt::format( "WARN: there were uncommitted changes at build-time." ) << std::endl;
        }
        const char* str = R"({:*^30}
 - commit {} (HEAD)
 - describe {}
 - Author: {} <{}>
 - Date: {}
 - Subject: {}
 - Body: {}
{:*^30})";
        std::cout << fmt::format( str, " Git Report ",
                                  GitMetadata::commitSHA1(),
                                  GitMetadata::describe(),
                                  GitMetadata::authorName(), GitMetadata::authorEmail(),
                                  GitMetadata::commitDate(),
                                  GitMetadata::commitSubject(),
                                  GitMetadata::commitBody(),
                                  " End Git Report " )
                  << std::endl;
    }
    else
    {
        LOG(INFO) << "WARN: failed to get the current git state. Is this a git repo?" << std::endl;
    }
}
template <typename E>
void
print_and_trace( std::string const& s, E const& e )
{
    
    const boost::stacktrace::stacktrace* st = boost::get_error_info<traced>( e );
    if ( st )
    {
        fmt::print( "{:*^30}\n", " Stack Trace " );
        std::cerr << *st << '\n';
        fmt::print( "{:*^30}\n", " Stack Trace " );
    }
    printGitReport();
    fmt::print( fmt::emphasis::bold | fg( fmt::color::red ), s );
}
void handleExceptions()
{

    try
    {
        throw; // re-throw exception already in flight } 
    }
    catch( boost::bad_lexical_cast const& e )
    {
        print_and_trace( fmt::format( "[feel++.boost.bad_lexical_cast] {}, source type: {}, target type: {}", 
                                        e.what(), 
                                        boost::core::demangle(e.source_type().name()), 
                                        boost::core::demangle(e.target_type().name())), e );
    }
    catch(const boost::filesystem::filesystem_error& e)
    {
        if ( e.code() == boost::system::errc::permission_denied )
            print_and_trace( fmt::format( "[feel++.boost.filesystem.filesystem_error.permission_denied] {}, path1: {}, path2: {}\n", e.what(), e.path1().string(), e.path2().string() ), e );
        else
            print_and_trace( fmt::format( "[feel++.boost.filesystem.filesystem_error] {}, path1: {}, path2: {}\n", e.what(), e.path1().string(), e.path2().string() ), e );
    }
    catch (json::exception& e)
    {
        print_and_trace(fmt::format("[feel++.json.exception] {}\n",e.what()),e);
    }
    catch ( std::invalid_argument const& e )
    {
        print_and_trace( fmt::format( "[feelpp.std.invalid_argument] {}\n",e.what()), e );
    }
    catch (const std::runtime_error & e) 
    {
        print_and_trace( fmt::format( "[feel++.std.runtime_error] {}", e.what() ), e );
    }
    catch ( const boost::mpi::exception& e )
    {
        print_and_trace( fmt::format( "[feel++.boost.mpi.exception] {}, routine: {}, result_code: {}\n", e.what(), e.routine(), e.result_code() ), e );
    }
    catch ( const std::exception& e )
    {
        print_and_trace( fmt::format( "[feel++.std.exception] {}", e.what() ), e );
    }
    catch(...)
    {
        fmt::print( "[feel++.unknown.exception] caught unknown exception" );
    }
}    

} // namespace Feel
