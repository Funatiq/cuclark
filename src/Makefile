CXXFLAGS = -O3
NVCC = nvcc #/usr/local/cuda/bin/nvcc
NVCCFLAGS = -arch=sm_30 -O3 #-D_FORCE_INLINES
# Checking for openmp
NB = $(shell echo |cpp -fopenmp -dM |grep -i open | wc -l)
ifeq ($(NB),1)
OPENMP = -Xcompiler -fopenmp
endif

CUCLARKCC = CuClarkDB.cu main.cc analyser.cc file.cc kmersConversion.cc
CUCLARK = $(CUCLARKCC) CuClarkDB.cuh CuCLARK_hh.hh dataType.hh hashTable_hh.hh HashTableStorage_hh.hh analyser.hh dataType.hh file.hh kmersConversion.hh

TPROGS = getTargetsDef getAccssnTaxID getfilesToTaxNodes #getGInTaxID
PROGS = cuCLARK cuCLARK-l $(TPROGS)

.PHONY: all clean target_definition debug

all: $(PROGS)

clean:
	rm -f $(PROGS) *.o
	
debug: NVCCFLAGS += -DTIME_DBLOADING -DDEBUG_DMEM #-DDEBUG_HMEM -DDEBUG_DB -DDEBUGQUERY -DDEBUG_KERNEL -DDEBUG_MERGE -DDEBUG_RESULT -DDEBUG_BATCH
debug: cuCLARK cuCLARK-l

cuCLARK: $(CUCLARK) parameters.hh
	$(NVCC) $(NVCCFLAGS) $(OPENMP) -o cuCLARK $(CUCLARKCC)

cuCLARK-l: $(CUCLARK) parameters_light_hh
	@mv parameters.hh parameters_full_hh
	@cp parameters_light_hh parameters.hh
	 $(NVCC) $(NVCCFLAGS) $(OPENMP) -o cuCLARK-l $(CUCLARKCC)
	@mv parameters_full_hh parameters.hh

target_definition: $(TPROGS)

getTargetsDef: getTargetsDef.cc file.cc file.hh
	$(CXX) $(CXXFLAGS) -o getTargetsDef getTargetsDef.cc file.cc

#getGInTaxID: getGInTaxID.cc file.cc file.hh
#	$(CXX) $(CXXFLAGS) -o getGInTaxID getGInTaxID.cc file.cc

getAccssnTaxID: getAccssnTaxID.cc file.cc file.hh
	$(CXX) $(CXXFLAGS) -o getAccssnTaxID getAccssnTaxID.cc file.cc

getfilesToTaxNodes: getfilesToTaxNodes.cc file.cc file.hh
	$(CXX) $(CXXFLAGS) -o getfilesToTaxNodes getfilesToTaxNodes.cc file.cc
