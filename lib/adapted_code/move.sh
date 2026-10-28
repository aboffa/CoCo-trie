#!/usr/bin/env bash
# Compatibility entry point for older build instructions. CMake now applies
# these adaptations to private copies under the build directory.
set -eu
printf '%s\n' 'No copying needed: configure with cmake -S . -B build to stage adapted headers.'
