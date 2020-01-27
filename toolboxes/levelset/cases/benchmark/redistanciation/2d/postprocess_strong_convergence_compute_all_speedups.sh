#!/usr/bin/env bash
ROOT_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d/results
FILENAME_PATTERN=redistanciation_2d_time_per_core
SPEEDUP_SCRIPT=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d/postprocess_strong_convergence_speedup.py

find ${ROOT_DIR} -name "${FILENAME_PATTERN}*" -exec ${SPEEDUP_SCRIPT} {} \; 
