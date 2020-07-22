#!/usr/bin/env bash

gfortran-8 -Wall -Wextra -g -pedantic -std=f2003 -Wimplicit-interface src/*Component.f03 src/*CLI.f03 src/Main.f03 -o test/sc
