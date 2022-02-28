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
#include <feel/feelcore/environment.hpp>

namespace Feel
{

void handleExceptions()
{
    try
    {
        throw; // re-throw exception already in flight } 
    }
    catch( boost::bad_lexical_cast const& e )
    {
        if ( Environment::isMasterRank() )
            fmt::print( "[feel++.boost.bad_lexical_cast] {}, source type: {}, target type: {}", 
                        e.what(), 
                        boost::core::demangle(e.source_type().name()), 
                        boost::core::demangle(e.target_type().name()) );
    }
    catch(const boost::filesystem::filesystem_error& e)
    {
        if ( e.code() == boost::system::errc::permission_denied )
            fmt::print( "[feel++.boost.filesystem.filesystem_error.permission_denied] {}\n", e.what() );
        else
            fmt::print( "[feel++.boost.filesystem.filesystem_error] {}\n", e.what() );
    }
    catch (json::exception& e)
    {
        fmt::print("[feel++.json.exception] {}\n",e.what());
    }
    catch ( std::invalid_argument const& e )
    {
        fmt::print( "[feelpp.std.invalid_argument] {}\n",e.what());
    }
    catch (const std::runtime_error & e) 
    {
        fmt::print( "[feel++.std.runtime_error] {}", e.what() );
    }
    catch ( const boost::mpi::exception& e )
    {
        fmt::print( "[feel++.boost.mpi.exception] {}, routine: {}, result_code: {}\n", e.what(), e.routine(), e.result_code() );
    }
    catch ( const std::exception& e )
    {
        fmt::print( "[feel++.std.exception] {}", e.what() );
    }
    catch(...)
    {
        fmt::print( "[feel++.unknown.exception] {}" );
    }
    LOG(WARNING) << fmt::format("[feel++] report status\n");
    if ( GitMetadata::populated() )
    {
        if ( GitMetadata::anyUncommittedChanges() )
        {
            std::cerr << fmt::format("WARN: there were uncommitted changes at build-time.") << std::endl;
        }
        const char* str = R"(commit {} (HEAD)
                             describe {}
                             Author: {} <{}>
                             Date: {}
                             Subject: {}
                             Body: {}
                            )";
        std::cout << fmt::format( str,
                        GitMetadata::commitSHA1(),
                        GitMetadata::describe(),
                        GitMetadata::authorName(),GitMetadata::authorEmail(),
                        GitMetadata::commitDate(),
                        GitMetadata::commitSubject(),
                        GitMetadata::commitBody() ) << std::endl;
    }
    else
    {
        std::cerr << "WARN: failed to get the current git state. Is this a git repo?" << std::endl;
    }
    LOG(WARNING) << fmt::format("[feel++] abort all mpi processes...\n");
    Environment::abort(EXIT_FAILURE);
}    

} // namespace Feel
