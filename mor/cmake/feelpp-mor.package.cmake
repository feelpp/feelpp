###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2020-01-20
#
#  Copyright (C) 2013-2020 Feel++ Consortium
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
##
## Archive generation using cpack
##
if (UNIX)
  execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE FEELPP_SYSTEM_MACHINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE FEELPP_SYSTEM_MACHINE_error
    RESULT_VARIABLE FEELPP_SYSTEM_MACHINE_result)
endif()
SET(CPACK_PACKAGE_NAME "feelpp-mor")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++ MOR")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.adoc")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.adoc")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FEELPP_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FEELPP_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FEELPP_VERSION_MICRO}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "feelpp-mor")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
SET(CPACK_SYSTEM_NAME "${FEELPP_OS}-${FEELPP_SYSTEM_MACHINE}")

SET(CPACK_PACKAGE_NAME "feelpp-mor")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++ MOR")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${FEELPP_VERSION_MAJOR}.${FEELPP_VERSION_MINOR}.${FEELPP_VERSION_MICRO}${FEELPP_VERSION_PRERELEASE}${FEELPP_VERSION_METADATA}")

if (NOT GIT_FOUND)
  find_package(Git QUIET)
endif()

function(feelpp_source_regex_literal INPUT OUTPUT)
  set(_value "${INPUT}")
  string(REPLACE "[" "[[]" _value "${_value}")
  string(REPLACE "]" "[]]" _value "${_value}")
  string(REPLACE "." "[.]" _value "${_value}")
  string(REPLACE "+" "[+]" _value "${_value}")
  string(REPLACE "*" "[*]" _value "${_value}")
  string(REPLACE "?" "[?]" _value "${_value}")
  string(REPLACE "$" "[$]" _value "${_value}")
  string(REPLACE "(" "[(]" _value "${_value}")
  string(REPLACE ")" "[)]" _value "${_value}")
  string(REPLACE "{" "[{]" _value "${_value}")
  string(REPLACE "}" "[}]" _value "${_value}")
  string(REPLACE "|" "[|]" _value "${_value}")
  set("${OUTPUT}" "${_value}" PARENT_SCOPE)
endfunction()

function(feelpp_append_git_untracked_source_ignores VARIABLE)
  set(_patterns ${${VARIABLE}})

  if (NOT GIT_FOUND)
    set("${VARIABLE}" "${_patterns}" PARENT_SCOPE)
    return()
  endif()

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}" ls-files --others --exclude-standard --directory --no-empty-directory
    OUTPUT_VARIABLE _git_untracked
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _git_untracked_result
  )

  if (NOT _git_untracked_result EQUAL 0 OR "${_git_untracked}" STREQUAL "")
    set("${VARIABLE}" "${_patterns}" PARENT_SCOPE)
    return()
  endif()

  string(REPLACE "\n" ";" _git_untracked_list "${_git_untracked}")
  foreach(_git_untracked_path IN LISTS _git_untracked_list)
    if (_git_untracked_path STREQUAL "")
      continue()
    endif()

    feelpp_source_regex_literal("${_git_untracked_path}" _git_untracked_regex)
    if (_git_untracked_path MATCHES "/$")
      list(APPEND _patterns "/${_git_untracked_regex}")
    else()
      list(APPEND _patterns "/${_git_untracked_regex}$")
    endif()
  endforeach()

  list(LENGTH _git_untracked_list _git_untracked_count)
  message(STATUS "[cpack] excluding ${_git_untracked_count} git-untracked path(s) from the source archive")
  set("${VARIABLE}" "${_patterns}" PARENT_SCOPE)
endfunction()

SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR};/")
set(CPACK_SOURCE_IGNORE_FILES
  "/[.]git/"
  "/[.]svn/"
  "/[.]venv[^/]*/"
  "/[.]pytest_cache/"
  "/[.]mypy_cache/"
  "/[.]cache/"
  "/[.]idea/"
  "/[.]vscode/"
  "/__pycache__/"
  "/build/"
  "/build[-_.][^/]*/"
  "/install/"
  "/install[-_.][^/]*/"
  "/_dist/"
  "/packaging/"
  "/doc/analysis/"
  "/feelpp_pkg[.]egg-info/"
  "CMakeLists[.]txt[.]user$"
  "gmsh-(config[.]err|info[.]log)$"
  "[.](tar[.](gz|bz2|xz)|deb|dsc|changes|build|buildinfo)$"
)
feelpp_append_git_untracked_source_ignores(CPACK_SOURCE_IGNORE_FILES)

include( CPack )
