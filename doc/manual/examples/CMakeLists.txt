# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2025-01-09
#
#  Copyright (C) 2025 Feel++ Consortium
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Custom repository examples
feelpp_add_application( example_custom_repository 
    SRCS example_custom_repository.cpp 
    LABELS doc 
    NO_TEST )

feelpp_add_application( example_custom_repository_env 
    SRCS example_custom_repository_env.cpp 
    LABELS doc 
    NO_TEST )

feelpp_add_application( example_custom_repository_cluster 
    SRCS example_custom_repository_cluster.cpp 
    LABELS doc 
    NO_TEST )

feelpp_add_application( test_feelio 
    SRCS test_feelio.cpp 
    LABELS doc 
    NO_TEST )
