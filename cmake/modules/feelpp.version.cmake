###  feelpp.directive.cmake; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-21
#
#  Copyright (C) 2013 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

################################################################################
# set the version of feel++
################################################################################

OPTION(FEELPP_ENABLE_GIT "enable Feel++ looking up for git information" OFF)
SET(FEELPP_SCM "git")
FIND_PACKAGE(Git)
if(GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git )
  message(STATUS "Git found" )
  if ( NOT EXISTS ${PROJECT_SOURCE_DIR}/contrib/ginac/ginac )
    message(STATUS "initializing submodules...")
    execute_process(COMMAND git submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} )
    message(STATUS "updating submodules...")
    execute_process(COMMAND git submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} )
  endif()

  execute_process(
    COMMAND git rev-parse --verify -q --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE Project_WC_REVISION
    ERROR_VARIABLE GIT_COMMIT_error
    RESULT_VARIABLE GIT_COMMIT_result)
  if ( NOT ${GIT_COMMIT_result} EQUAL 0 )
    MESSAGE(SEND_ERROR "Command \"git rev-parse --verify -q --short git-svn\" failed with output:\n${GIT_COMMIT_error}")
  else()
    STRING(STRIP "${Project_WC_REVISION}" Project_WC_REVISION )
    MESSAGE(STATUS "Git commit: ${Project_WC_REVISION}")
  endif()
else()
  if(EXISTS ${PROJECT_SOURCE_DIR}/GITREVISION)
    file(READ ${PROJECT_SOURCE_DIR}/GITREVISION Project_WC_REVISION)
  endif()
endif()


set(FEELPP_VERSION_MAJOR "0")
set(FEELPP_VERSION_MINOR "100")
set(FEELPP_VERSION_MICRO "0")
set(FEELPP_REVISION "0" )
set(FEELPP_BUILDID "0" )
set(FEELPP_VERSION_PRERELEASE "-beta.7" )
if (FEELPP_ENABLE_GIT AND GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git )
  set(FEELPP_VERSION_METADATA "+${FEELPP_SCM}${Project_WC_REVISION}")
  set(FEELPP_REVISION "${Project_WC_REVISION}")
  set(FEELPP_BUILDID "${Project_WC_REVISION}")
  file(WRITE "${CMAKE_SOURCE_DIR}/GITREVISION" "${Project_WC_REVISION}")
else()
  set(FEELPP_VERSION_METADATA "")
  if ( EXISTS "${CMAKE_SOURCE_DIR}/GITREVISION" )
    file(READ "${CMAKE_SOURCE_DIR}/GITREVISION" FEELPP_REVISION)
    file(READ "${CMAKE_SOURCE_DIR}/GITREVISION" FEELPP_BUILDID)
  endif()
endif()
set(FEELPP_VERSION "(((${FEELPP_VERSION_MAJOR}) << 16) | ((${FEELPP_VERSION_MINOR}) << 9) | (${FEELPP_VERSION_MICRO}))" )
set(FEELPP_VERSION_STRING "${FEELPP_VERSION_MAJOR}.${FEELPP_VERSION_MINOR}.${FEELPP_VERSION_MICRO}${FEELPP_VERSION_PRERELEASE}${FEELPP_VERSION_METADATA}")
set(FEELPP_SHARED_VERSION "1.0.0" )
set(FEELPP_SHARED_SOVERSION "1" )
