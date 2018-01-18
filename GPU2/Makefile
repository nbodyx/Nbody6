#avx       = enable
chainmode = archain
#ksblock   = on

ifeq ($(avx), enable)
MKOPTIONS += avx=enable
endif
ifeq ($(chainmode), archain)
MKOPTIONS += chainmode=archain
endif
ifeq ($(ksblock), on)
MKOPTIONS += ksblock=on
endif

ifeq ($(chainmode), archain)
 ifeq ($(ksblock), on)
  BINNAME = nbody7b
 else
  BINNAME = nbody7
 endif
else
 ifeq ($(ksblock), on)
  BINNAME = nbody6b
 else
  BINNAME = nbody6
 endif
endif

FC       = gfortran
FFLAGS   = -O3 -Wall -march=core2
CXX      = g++ 
CXXFLAGS = -O3 -Wall -march=core2
ifeq ($(avx), enable)
FFLAGS   += -mavx
CXXFLAGS += -mavx
endif

CUDA_PATH = /usr/local/cuda
#SDK_PATH  = /usr/local/cuda_sdk
SDK_PATH  = /usr/local/cuda/samples

NVCC = $(CUDA_PATH)/bin/nvcc
NVCC += -arch=sm_20 -Xptxas -v
NVCC += -I $(SDK_PATH)/common/inc
ifeq ($(avx), enable)
NVCC += -Xcompiler "-O3 -Wall -fopenmp -mavx"
else
NVCC += -Xcompiler "-O3 -Wall -fopenmp"
endif
NVCC += -DWITH_CUDA5

RUNDIR = ./run

ifeq ($(avx), enable)
all: gpu avx
else
all: gpu sse
endif

gpu: gpunb.gpu.o gpupot.gpu.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o
	cp -f gpunb.gpu.o gpunb.o
	cp -f gpupot.gpu.o gpupot.o
	cp -f *.o ./Build
	make gpu -C ./Build -f $(shell pwd)/Makefile.build -j4 $(MKOPTIONS)
	mv -f ./Build/nbody6 $(RUNDIR)/$(BINNAME).gpu

avx: gpunb.avx.o gpupot.avx.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o
	cp -f gpunb.avx.o gpunb.o
	cp -f gpupot.avx.o gpupot.o
	cp -f *.o ./Build
	make sse -C ./Build -f $(shell pwd)/Makefile.build -j4 $(MKOPTIONS)
	mv -f ./Build/nbody6 $(RUNDIR)/$(BINNAME).avx

sse: gpunb.sse.o gpupot.sse.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o
	cp -f gpunb.sse.o gpunb.o
	cp -f gpupot.sse.o gpupot.o
	cp -f *.o ./Build
	make sse -C ./Build -f $(shell pwd)/Makefile.build -j4 $(MKOPTIONS)
	mv -f ./Build/nbody6 $(RUNDIR)/$(BINNAME).sse

cnbint.o: ./lib/cnbint.cpp
	$(CXX) $(CXXFLAGS) -c $^

wtime.o: ./lib/wtime.cpp
	$(CXX) $(CXXFLAGS) -c $^

cxvpred.o: ./lib/cxvpred.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $^

ifeq ($(avx), enable)
gpuirr.o: gpuirr.s
	$(CXX) -c $^ -o $@

gpuirr.s: gpuirr.avx.s
	grep -E -v 'v?movaps\s(%xmm[0-9]+),\s?\1$$' < $^ > $@ #Thanks, Sakuraba

gpuirr.avx.s: irrlib/gpuirr.avx.cpp
	$(CXX) $(CXXFLAGS) -mavx -fopenmp -S $^ -o $@

else
gpuirr.o: ./irrlib/gpuirr.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $^ -o $@

endif
intgrt.o: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -c -o $@ 

intgrt.s: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -S -o $@

gpunb.gpu.o: lib/gpunb.velocity.cu
	$(NVCC) $^ -c -o $@

gpunb.velocity.su: gpunb.velocity.cubin
	cuobjdump -sass $< | c++filt > $@

gpunb.velocity.cubin: lib/gpunb.velocity.cu
	$(NVCC) -cubin --ptxas-options=-v $<

gpupot.gpu.o: lib/gpupot.gpu.cu
	$(NVCC) $^ -c -o $@

gpunb.avx.o: lib/gpunb.avx.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpupot.avx.o: lib/gpupot.avx.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpunb.sse.o: lib/gpunb.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpupot.sse.o: lib/gpupot.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@


clean:
	rm -f *.o ./Build/*.o

