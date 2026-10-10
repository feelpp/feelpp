# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-2.1-or-later

set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;teststats")
feelpp_add_test( iterative_statistics NO_MPI_TEST )
feelpp_add_test( kll NO_MPI_TEST )
feelpp_add_test( otherstatistics NO_MPI_TEST )
