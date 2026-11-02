// -*- mode: hcl -*-
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

variable "REGISTRY" {
  default = "ghcr.io/feelpp"
}

variable "BASE_IMAGE" {
  default = "ubuntu:24.04"
}

variable "OPENFOAM_PACKAGE_VERSION" {
  default = "2606.0~rc2-1"
}

variable "BUILD_JOBS" {
  default = "4"
}

group "default" {
  targets = ["openfoam2606-precice"]
}

target "openfoam2606" {
  context = "."
  dockerfile = "Dockerfile"
  target = "openfoam"
  platforms = ["linux/amd64"]
  args = {
    BASE_IMAGE = BASE_IMAGE
    OPENFOAM_PACKAGE_VERSION = OPENFOAM_PACKAGE_VERSION
  }
  tags = ["${REGISTRY}/openfoam:2606-ubuntu24.04"]
}

target "openfoam2606-precice" {
  inherits = ["openfoam2606"]
  target = "openfoam-precice"
  args = {
    BUILD_JOBS = BUILD_JOBS
  }
  tags = ["${REGISTRY}/openfoam:2606-precice3.4.1-ubuntu24.04"]
}
