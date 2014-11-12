###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-18
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

# Set to "ON" if you want the LaTeX documentation to be built.
set(DOXYFILE_LATEX 1)

find_program(DOXYGEN doxygen)
if ( DOXYGEN AND FEELPP_ENABLE_DOXYGEN)

  ## We extract the cmake option for the documentation
  execute_process(
    COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/extract_command.sh ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE FEELPP_ENABLE_LISTE
    )
   # doxygen is optional
   include(UseDoxygen OPTIONAL)

   # make sure that doxygen is not run byt default
   set_target_properties(doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)

   # add subpages to projects/modules
   configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc.in
     ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc
     )

   set(FEELPP_PROJECTS_DIR,"")
	 set(FEELPP_DOXYGEN_PROJECT "OFF")
   # if projects have been detected
   foreach(project ${projects})
     #get shortname
     GET_FILENAME_COMPONENT(dirname ${project} NAME)

     #define option to control doc generation
     string(TOUPPER ${dirname} PROJECT)
     if (FEELPP_ENABLE_${PROJECT}_DOCUMENTATION )
			 set(FEELPP_DOXYGEN_PROJECT "ON")
		 endif()
	 endforeach()
	 if(${FEELPP_DOXYGEN_PROJECT})
       FILE(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc "/** \\page Project Project\n")
		 endif()
   
	 foreach(project ${projects})
     #get shortname
     GET_FILENAME_COMPONENT(dirname ${project} NAME)

     #define option to control doc generation
     string(TOUPPER ${dirname} PROJECT)
     if (FEELPP_ENABLE_${PROJECT}_DOCUMENTATION )
       set(FEELPP_PROJECTS_DIR ${FEELPP_PROJECTS_DIR} "\"${CMAKE_SOURCE_DIR}/research/${project}/\" \\ \n")
       FILE(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc "   - \\subpage ${dirname}\n")
       MESSAGE(STATUS "Doxygen Documention for ${project} will be generated")
     endif()
   endforeach()
   FILE(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc "*/\n")

   #set(FEELPP_PROJECTS_DIR ${FEELPP_PROJECTS_DIR} "\"\" \n") #force empty directory
   STRING(REPLACE ";" " " _TMP_STR "${FEELPP_PROJECTS_DIR}")
   SET(FEELPP_PROJECTS_DIR "${_TMP_STR}")

   FILE(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/pages.doc "\n\n}\n")


   #
   configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
     ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
     )

   configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/latexheader.tex.in
     ${CMAKE_CURRENT_BINARY_DIR}/latexheader.tex
     )
   file(GLOB LOGOS "${CMAKE_CURRENT_SOURCE_DIR}/pngs/logos/*")
   file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html)
   foreach(logo ${LOGOS})
     GET_FILENAME_COMPONENT(logoname ${logo} NAME)
     configure_file( ${logo} "${CMAKE_CURRENT_BINARY_DIR}/html/${logoname}" COPYONLY)
   endforeach()
endif()
