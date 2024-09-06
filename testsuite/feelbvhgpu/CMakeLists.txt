###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2023-11-20
#
#  Copyright (C) 2023 Universite de Strasbourg
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


#set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")



project(test_bvhgpu)

set_directory_properties(PROPERTIES LABEL testbvhgpu )


if(FALSE)
   find_package(rocthrust REQUIRED)
   if (rocthrust_FOUND)
      message("{rocthrust_LIBRARIES}") 
   else()
      message(FATAL_ERROR "Could not find rocThrust")     
   endif()
endif()


if(TRUE)
   find_package(rocthrust REQUIRED CONFIG PATHS "/nvme0/lemoinep/soft/rocThrust")
   if (rocthrust_FOUND)
      message("{rocthrust_LIBRARIES}") 
      set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/nvme0/lemoinep/soft/rocThrust/include")
      list(APPEND include_dirs "/nvme0/lemoinep/soft/rocThrust/include")
   else()
      message(FATAL_ERROR "Could not find rocThrust")     
   endif()
endif()



include(CheckLanguage)
check_language(HIP)
if(CMAKE_HIP_COMPILER)
   enable_language(HIP)
endif()

find_package(hip)
feelpp_add_test(bvhgpu SRCS test_bvhgpu.cpp CFG test_bvhgpu.cfg GEO cubic.geo 
  LINK_LIBRARIES ${GPU_HIP_LIBRARIES} Feelpp::feelpp  hip::device hip::host hipblas roc::rocthrust)
  set_source_files_properties( test_bvhgpu.cpp 
  PROPERTIES LANGUAGE ${GPU_RUNTIME})


