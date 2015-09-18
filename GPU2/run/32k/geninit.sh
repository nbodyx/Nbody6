#!/bin/sh
rm -f fort.2
../nbody6.gpu < input > INIT
mv -f fort.2 fort.1
