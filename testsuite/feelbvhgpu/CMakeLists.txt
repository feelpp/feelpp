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

#feelpp_add_test(bvhgpu SRCS test_bvhgpu.cpp LINK_LIBRARIES Feelpp::feelpp_gpu )
feelpp_add_test(bvhgpu SRCS test_old.cpp LINK_LIBRARIES Feelpp::feelpp_gpu )
#feelpp_add_test(bvhgpu SRCS test_old_empty.cpp LINK_LIBRARIES Feelpp::feelpp_gpu )
set_target_properties(feelpp_test_bvhgpu PROPERTIES HIP_LANGUAGE CXX)


