/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
      Date: 2022-02-19

 Copyright (C) 2022 Feel++ Consortium


 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_GIT_HPP
#define FEELPP_GIT_HPP 1 

#include <string>

namespace Feel
{

/**
 * @brief GitMetaData class holder
 * 
 * This class holds the git metadata if available
 */
class GitMetadata
{
  public:
    /**
     * @brief check if git data is available
     * if there wasn't a .git directory (e.g. downloaded source
     *  code without revision history).
     */
    static bool populated();

    /**
     * @brief Were there any uncommitted changes that won't be reflected in the CommitID?
     *
     */
    static bool anyUncommittedChanges();

    //! The commit author's name.
    static std::string authorName();
    //! The commit author's email.
    static std::string authorEmail();
    //! The commit SHA1.
    static std::string commitSHA1();
    //! The ISO8601 commit date.
    static std::string commitDate();
    //! The commit subject.
    static std::string commitSubject();
    //! The commit body.
    static std::string commitBody();
    //! The commit describe.
    static std::string describe();
};

} // namespace Feel
#endif