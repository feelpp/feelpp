###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-07-15
#
#  Copyright (C) 2014-2016 Feel++ Consortium
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

set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testopt")

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_PETSC_TAO)
if ( FEELPP_HAS_PETSC_TAO )
  feelpp_add_test( tao )
  feelpp_add_test( tao_functional )

  # All rank variants share the same generated-expression repository.  Keep
  # CTest from executing them concurrently, which can otherwise race while
  # compiling and loading the generated target expressions.
  set_property(
    TEST feelpp_test_tao_functional-np-1
    PROPERTY RESOURCE_LOCK feelpp_test_tao_functional_repository )
  if ( NProcs2 GREATER 1 )
    set_property(
      TEST feelpp_test_tao_functional-np-${NProcs2}
      PROPERTY RESOURCE_LOCK feelpp_test_tao_functional_repository )
  endif()

  # The standard helper registers np-1 and the configured half-core count.
  # Keep an explicit two-rank regression so TAO field and block layouts are
  # exercised on 1, 2, and 5 ranks in the macOS Spack development preset.
  if ( NProcs2 GREATER 1 AND NOT NProcs2 EQUAL 2 )
    add_test(
      NAME feelpp_test_tao_functional-np-2
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:feelpp_test_tao_functional>
              --log_level=all -- ${MPIEXEC_POSTFLAGS}
              --directory=testsuite/test_tao_functional --rm )
    set_property(
      TEST feelpp_test_tao_functional-np-2
      PROPERTY LABELS "testopt;testsuite;testsuite-feelpp" )
    set_property(
      TEST feelpp_test_tao_functional-np-2 PROPERTY PROCESSORS 2 )
    set_property(
      TEST feelpp_test_tao_functional-np-2
      PROPERTY TIMEOUT ${FEELPP_DEFAULT_TEST_TIMEOUT} )
    set_property(
      TEST feelpp_test_tao_functional-np-2
      PROPERTY RESOURCE_LOCK feelpp_test_tao_functional_repository )
    feelpp_append_solver_runtime_env( feelpp_test_tao_functional-np-2 )
  endif()
endif()

# for applications using NLOpt
if ( FEELPP_HAS_NLOPT )

  feelpp_add_test( nlopt )

endif()

if ( FEELPP_HAS_IPOPT )

  feelpp_add_test( ipopt )

endif()
