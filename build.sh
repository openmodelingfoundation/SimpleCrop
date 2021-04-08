#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

mkdir -p src/_build
cd src/_build
cmake ../
make all
