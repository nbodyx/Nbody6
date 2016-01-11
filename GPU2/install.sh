#!/bin/sh
make clean
make -C ../Ncode -f ../GPU2/Makefile_gpu clean
rm run/nbody6.gpu run/nbody6.sse
rm -f params.h common6.h
ln -s ../Ncode/params.h ./
ln -s ../Ncode/common6.h ./
