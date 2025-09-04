/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Thomas Saigre <thomas.saigre@cemosis.fr>
       Date: 2025-09-04

  Copyright (C) 2025 Université de Strasbourg

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
/**
   \file test_journal.cpp
   \author Thomas Saigre <thomas.saigre@cemosis.fr>
   \date 2025-09-04
 */


#define BOOST_TEST_MODULE test_journal
#include <feel/feel.hpp>
#include <feel/feelcore/testsuite.hpp>

using namespace Feel;

nl::json ptree_1 = {
    { "testsuite", {
        { "cases", { "case1", "case2", "case3" } }
    } }
};
nl::json ptree_2 = {
    { "testsuite", {
        { "config", true  }
    } }
};

nl::json merged_ptree = {
    {"testsuite", {
        { "cases", {"case1", "case2", "case3" } },
        { "config", true }
    } }
};

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_journal"), Feel::feel_options() )
BOOST_AUTO_TEST_SUITE( journal )

BOOST_AUTO_TEST_CASE( test_add_ptree )
{
    JournalManager::journalAdd( ptree_1 );

    const std::string filename = "journal_test_add_ptree.json";
    JournalManager::journalCheckpoint( true, filename );

    nl::json loaded_ptree;
    std::ifstream infile(filename);
    BOOST_REQUIRE(infile.is_open());
    infile >> loaded_ptree;
    infile.close();

    BOOST_CHECK_EQUAL(loaded_ptree["testsuite"].dump(), ptree_1["testsuite"].dump());
}

BOOST_AUTO_TEST_CASE( test_add_ptrees )
{
    JournalManager::journalAdd( ptree_1 );
    JournalManager::journalAdd( ptree_2 );

    const std::string filename = "journal_test_add_ptrees.json";
    JournalManager::journalCheckpoint( true, filename );

    nl::json loaded_ptree;
    std::ifstream infile(filename);
    BOOST_REQUIRE(infile.is_open());
    infile >> loaded_ptree;
    infile.close();

    BOOST_CHECK_EQUAL(loaded_ptree["testsuite"].dump(), merged_ptree["testsuite"].dump());
}

BOOST_AUTO_TEST_CASE( test_add_ptrees_not_merged )
{
    JournalManager::journalAdd( ptree_1 );
    JournalManager::journalAdd( ptree_2, false );

    const std::string filename = "journal_test_add_ptrees_not_merged.json";
    JournalManager::journalCheckpoint( true, filename );

    nl::json loaded_ptree;
    std::ifstream infile(filename);
    BOOST_REQUIRE(infile.is_open());
    infile >> loaded_ptree;
    infile.close();

    BOOST_CHECK_EQUAL(loaded_ptree["testsuite"].dump(), ptree_2["testsuite"].dump());
}


BOOST_AUTO_TEST_SUITE_END()