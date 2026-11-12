/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2005,2006 EPFL
 * SPDX-FileCopyrightText: 2006-2011 Université Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2011-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file vf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#ifndef FEELPP_VF_HPP
#define FEELPP_VF_HPP 1

/**
 * \brief allow automatic type naming of complex expression
 */
#define AUTO( a, b ) auto a = (b);

#include <feel/feelvf/vf_dsl.hpp>
#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/basis.hpp>
#include <feel/feelvf/field.hpp>
#include <feel/feelvf/facet.hpp>
#include <feel/feelvf/shellgeometric.hpp>
#include <feel/feelvf/sb9.hpp>
#include <feel/feelvf/tensorbasis.hpp>
#include <feel/feelvf/voigt.hpp>
#include <feel/feelvf/contractions.hpp>

#if defined( FEELPP_HAS_FFTW)
#include <feel/feelvf/msi.hpp>
#endif


#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/symbolicdiagnostics.hpp>

#include <boost/preprocessor/comparison/equal.hpp>

namespace Feel
{
using namespace vf;
}
/// \endcond

#endif /* __VF_H */
