//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feelcrb/biotsavartrb.hpp>
#include <feel/feelcrb/crbsaddlepoint.hpp>
#include "biotsavart.hpp"
#include "thermoelectric-nonlinear.hpp"

using namespace Feel;

int main( int argc, char** argv )
{
    using te_type = ThermoElectric;
    using BSModel_type = BiotSavart<te_type>;
    using BSRB_type = BiotSavartRB<BSModel_type>;
    using BSRB_ptrtype = boost::shared_ptr<BSRB_type>;

    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="biotsavart-nonlinear"),
                     _desc=BSModel_type::makeOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(crbSaddlePointOptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB")) );

    auto bs = BSRB_type::New( soption("biotsavart.basename"), crb::stage::offline );
    bs->offline();
    auto mu = bs->paramFromProperties();
    auto vtN = bs->onlineVT( mu );
    auto bn = bs->online( mu, vtN );
    auto B = bs->expansion(bn);
    auto VT = bs->expansionVT(vtN );
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();

    auto e = exporter( _mesh=bs->mesh(), _name="biotsavart" );
    e->add("B", B );
    e->add("V", V );
    e->add("T", T );
    e->save();

    bs->exportBasis();

    return 0;
}
