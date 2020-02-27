###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 28 Feb 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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
include (GNUInstallDirs)
Find_Package(Asciidoctor)

if ( ASCIIDOCTOR_FOUND )
  set (FEELPP_DOCDIR ${CMAKE_CURRENT_SOURCE_DIR}/doc)
  set (FEELPP_STYLESHEET ${FEELPP_DOCDIR}/stylesheet.css)
  set (FEELPP_TITLE ${PROJECT_NAME} ${FEELPP_PACKAGE_VERSION})
  set (FEELPP_A2P ${ASCIIDOCTOR_PDF_EXECUTABLE}  -amanmanual='${FEELPP_TITLE}')
  set (FEELPP_HAS_ASCIIDOCTOR_PDF 1)
  set (FEELPP_A2M ${ASCIIDOCTOR_EXECUTABLE} -b manpage -amanmanual='${FEELPP_TITLE}')
  set (FEELPP_A2H ${ASCIIDOCTOR_EXECUTABLE} -d manpage -b html5 -a stylesheeet=${FEELPP_STYLESHEET} -aversion-label=${PROJECT_NAME} -arevnumber=${FEELPP_PACKAGE_VERSION})
  set (FEELPP_A2M_STR "${ASCIIDOCTOR_EXECUTABLE} -b manpage -amanmanual='${FEELPP_TITLE}'")
  set (FEELPP_A2H_STR "${ASCIIDOCTOR_EXECUTABLE} -d manpage -b html5 -a stylesheeet=${FEELPP_STYLESHEET} -aversion-label=${PROJECT_NAME} -arevnumber=${FEELPP_PACKAGE_VERSION}")

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check.asciidoc.man.adoc "= toto(1)\n\n== NAME\n\nabc - def")
  execute_process(
    COMMAND ${FEELPP_A2M} check.asciidoc.man.adoc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}
    OUTPUT_FILE check.asciidoc.manpage.log
    ERROR_FILE check.asciidoc.manpage.err
    RESULT_VARIABLE ERROR_CODE
    )
  if(ERROR_CODE EQUAL "0")
    set(FEELPP_HAS_ASCIIDOCTOR_MANPAGE 1)
  else()
    unset(FEELPP_HAS_ASCIIDOCTOR_MANPAGE)
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check.asciidoc.html.adoc "")
  execute_process(
    COMMAND ${FEELPP_A2H} check.asciidoc.html.adoc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}
    OUTPUT_FILE check.asciidoc.html5.log
    ERROR_FILE check.asciidoc.html5.err
    RESULT_VARIABLE ERROR_CODE
    )
  if(ERROR_CODE EQUAL "0")
    set(FEELPP_HAS_ASCIIDOCTOR_HTML5 1)
  else()
    unset(FEELPP_HAS_ASCIIDOCTOR_HTML5)
  endif()
  if ( FEELPP_HAS_ASCIIDOCTOR_MANPAGE OR FEELPP_HAS_ASCIIDOCTOR_HTML5 OR FEELPP_HAS_ASCIIDOCTOR_PDF)
    set (FEELPP_HAS_ASCIIDOCTOR 1)
    set(FEELPP_ASCIIDOCTOR_ENABLED_BACKEND)
    if (FEELPP_HAS_ASCIIDOCTOR_MANPAGE)
      list(APPEND FEELPP_ASCIIDOCTOR_ENABLED_BACKEND "manpage")
    endif()
    if (FEELPP_HAS_ASCIIDOCTOR_HTML5)
      list(APPEND FEELPP_ASCIIDOCTOR_ENABLED_BACKEND "html5")
    endif()
    if (FEELPP_HAS_ASCIIDOCTOR_PDF)
      list(APPEND FEELPP_ASCIIDOCTOR_ENABLED_BACKEND "pdf")
    endif()
    SET( FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Asciidoctor(${FEELPP_ASCIIDOCTOR_ENABLED_BACKEND})" )
  endif()
endif( ASCIIDOCTOR_FOUND )
