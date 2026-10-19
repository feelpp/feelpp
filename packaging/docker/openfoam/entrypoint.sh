#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

set -o pipefail
# OpenFOAM's setup probes optional libraries with nonzero return statuses.
# Supply an empty argument so it does not interpret the command as settings.
source "${FOAM_BASHRC}" '' || exit $?
set -e
export LD_LIBRARY_PATH="${PRECICE_OPENFOAM_TARGET_DIR}:${LD_LIBRARY_PATH:-}"
exec "$@"
