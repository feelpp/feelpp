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
  set (FEELPP_A2M ${ASCIIDOCTOR_EXECUTABLE} -b manpage -amanmanual='${FEELPP_TITLE}')
  set (FEELPP_A2H ${ASCIIDOCTOR_EXECUTABLE} -d manpage -b html5 -a stylesheeet=${FEELPP_STYLESHEET} -aversion-label=${PROJECT_NAME} -arevnumber=${FEELPP_PACKAGE_VERSION})

  add_custom_target (man)
  add_custom_target (html)
  macro (feelpp_add_man NAME MAN SECT)
    message(STATUS "building manual page ${NAME}.${SECT}")
    add_custom_target(${NAME}.${SECT})
    add_custom_target(${NAME}.${SECT}.html)
    add_custom_command (
      TARGET ${NAME}.${SECT}
      #OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}
      COMMAND ${FEELPP_A2M} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT} ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
      )
    #add_custom_target(${NAME}.${SECT} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT})
    add_custom_command (
      TARGET ${NAME}.${SECT}.html
      #OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html
      COMMAND ${FEELPP_A2H} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
      DEPENDS ${FEELPP_STYLESHEET}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
      )
    #add_custom_target(${NAME}.${SECT}.html DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html)
    
    #set(FEELPP_MANS ${FEELPP_MANS} ${NAME}.${SECT})
    #set(FEELPP_HTMLS ${FEELPP_HTMLS} ${NAME}.${SECT}.html)
    add_dependencies(man ${NAME}.${SECT})
    add_dependencies(html ${NAME}.${SECT}.html)

    install (
      FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html
      DESTINATION ${CMAKE_INSTALL_DOCDIR}
      COMPONENT Bin
      )
    install (
      FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}
      DESTINATION ${CMAKE_INSTALL_MANDIR}/man${SECT}
      COMPONENT Bin
      )

  endmacro (feelpp_add_man)
else( ASCIIDOCTOR_FOUND )
  macro (feelpp_add_man NAME SECT)
    # no-op
  endmacro(feelpp_add_man)
endif( ASCIIDOCTOR_FOUND )

