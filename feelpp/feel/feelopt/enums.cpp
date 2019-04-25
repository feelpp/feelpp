/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <feel/feelopt/enums.hpp>

namespace Feel {

std::map<std::string, nlopt::algorithm> nloptAlgoMap = {
    {"GN_DIRECT", nlopt::GN_DIRECT },
    {"GN_DIRECT_L", nlopt::GN_DIRECT_L },
    {"GN_DIRECT_L_RAND", nlopt::GN_DIRECT_L_RAND },
    {"GN_DIRECT_NOSCAL", nlopt::GN_DIRECT_NOSCAL },
    {"GN_DIRECT_L_NOSCAL", nlopt::GN_DIRECT_L_NOSCAL },
    {"GN_DIRECT_L_RAND_NOSCAL", nlopt::GN_DIRECT_L_RAND_NOSCAL },

    {"GN_ORIG_DIRECT", nlopt::GN_ORIG_DIRECT },
    {"GN_ORIG_DIRECT_L", nlopt::GN_ORIG_DIRECT_L },

    {"GD_STOGO", nlopt::GD_STOGO },
    {"GD_STOGO_RAND", nlopt::GD_STOGO_RAND },

    {"LD_LBFGS_NOCEDAL", nlopt::LD_LBFGS_NOCEDAL },

    {"LD_LBFGS", nlopt::LD_LBFGS },

    {"LN_PRAXIS", nlopt::LN_PRAXIS },

    {"LD_VAR1", nlopt::LD_VAR1 },
    {"LD_VAR2", nlopt::LD_VAR2 },

    {"LD_TNEWTON", nlopt::LD_TNEWTON },
    {"LD_TNEWTON_RESTART", nlopt::LD_TNEWTON_RESTART },
    {"LD_TNEWTON_PRECOND", nlopt::LD_TNEWTON_PRECOND },
    {"LD_TNEWTON_PRECOND_RESTART", nlopt::LD_TNEWTON_PRECOND_RESTART },

    {"GN_CRS2_LM", nlopt::GN_CRS2_LM },

    {"GN_MLSL", nlopt::GN_MLSL },
    {"GD_MLSL", nlopt::GD_MLSL },
    {"GN_MLSL_LDS", nlopt::GN_MLSL_LDS },
    {"GD_MLSL_LDS", nlopt::GD_MLSL_LDS },

    {"LD_MMA", nlopt::LD_MMA },

    {"LN_COBYLA", nlopt::LN_COBYLA },

    {"LN_NEWUOA", nlopt::LN_NEWUOA },
    {"LN_NEWUOA_BOUND", nlopt::LN_NEWUOA_BOUND },

    {"LN_NELDERMEAD", nlopt::LN_NELDERMEAD },
    {"LN_SBPLX", nlopt::LN_SBPLX },

    {"LN_AUGLAG", nlopt::LN_AUGLAG },
    {"LD_AUGLAG", nlopt::LD_AUGLAG },
    {"LN_AUGLAG_EQ", nlopt::LN_AUGLAG_EQ },
    {"LD_AUGLAG_EQ", nlopt::LD_AUGLAG_EQ },

    {"LN_BOBYQA", nlopt::LN_BOBYQA },

    {"GN_ISRES", nlopt::GN_ISRES },

    /* new variants that require local_optimizer to be set,
       not with older constants for backwards compatibility */
    {"AUGLAG", nlopt::AUGLAG },
    {"AUGLAG_EQ", nlopt::AUGLAG_EQ },
    {"G_MLSL", nlopt::G_MLSL },
    {"G_MLSL_LDS", nlopt::G_MLSL_LDS },

    {"LD_SLSQP", nlopt::LD_SLSQP },

    {"LD_CCSAQ", nlopt::LD_CCSAQ },

    {"GN_ESCH", nlopt::GN_ESCH },

    {"GN_AGS", nlopt::GN_AGS }
};

std::map<nlopt::result, std::string> nloptResultMap = {
    {nlopt::FAILURE, "NLOPT Generic Failure!"},
    {nlopt::INVALID_ARGS, "NLOPT Invalid arguments!"},
    {nlopt::OUT_OF_MEMORY, "NLOPT Out of memory!"},
    {nlopt::ROUNDOFF_LIMITED, "NLOPT Roundoff limited!"},
    {nlopt::FORCED_STOP, "NLOPT Forced stop!"},
    {nlopt::SUCCESS, "NLOPT Success"},
    {nlopt::STOPVAL_REACHED, "NLOPT Stop value reached!"},
    {nlopt::FTOL_REACHED, "NLOPT ftol reached!"},
    {nlopt::XTOL_REACHED, "NLOPT xtol reached!"},
    {nlopt::MAXEVAL_REACHED, "NLOPT Maximum number of evaluation reached!"},
    {nlopt::MAXTIME_REACHED, "NLOPT Maximum time reached"}
};

} // namespace Feel
