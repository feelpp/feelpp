###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-12-08
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

if ( APPLE )
  if (EXISTS /usr/local/lib/petscdir/3.4.3/darwin-cxx-debug/ )
    set(ENV{PETSC_DIR} /usr/local/lib/petscdir/3.4.3/darwin-cxx-debug/)
    if (EXISTS /usr/local/lib/slepcdir/3.4.3/darwin-cxx-debug/ )
      set(ENV{SLEPC_DIR} /usr/local/lib/slepcdir/3.4.3/darwin-cxx-debug/)
    endif()
  endif()
  if (EXISTS /usr/local/opt/scotch5/)
    set(ENV{PTSCOTCH_DIR} /usr/local/opt/scotch5/)
  endif()
endif()

message(STATUS "Feel++ uses a Homebrew system")
