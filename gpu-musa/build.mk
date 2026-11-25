CXX ?= g++
ARCH ?= mp_31
SDK ?= /usr/local/musa
CXX_MAJOR := $(shell $(CXX) -dumpversion | cut -d. -f1)
TRIPLE := $(shell $(CXX) -dumpmachine)
HOST_FLAGS := -std=c++17 -O3 -DNDEBUG -DSMALL_BOX -DOPENCL_3_0 -DBOOST_TIMER_ENABLE_DEPRECATED -I$(BOOST)/include -Isrc/lib -Isrc/OpenCL/inc -I$(SDK)/include -Iport
# Jammy can have GCC 12 runtime files but only GCC 11 C++ headers.
DEVICE_FLAGS := -std=c++17 -O3 --offload-arch=$(ARCH) --musa-path=$(SDK) -isystem /usr/include/c++/$(CXX_MAJOR) -isystem /usr/include/$(TRIPLE)/c++/$(CXX_MAJOR) -funsafe-math-optimizations -ffp-contract=fast -Xclang -cl-single-precision-constant -DSMALL_BOX -Isrc/OpenCL/inc -Isrc/OpenCL/src/kernels
LIBS := -L$(BOOST)/lib -Wl,-rpath,'$$ORIGIN/lib' -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L$(SDK)/lib -Wl,-rpath,$(SDK)/lib -lmusart -lstdc++fs -lpthread -lm
SOURCES := $(wildcard src/lib/*.cpp) src/main/main.cpp runtime/musa_backend.cpp
HOST_OBJS := $(SOURCES:%.cpp=obj/%.o)
OBJS := $(HOST_OBJS) obj/kernel1.o obj/kernel2.o
ifneq ($(wildcard reference/grid_fallback.mu),)
OBJS += obj/kernel1_reference.o
endif

AutoDock-Vina-GPU-2-1: $(OBJS)
	$(CXX) -no-pie $^ $(LIBS) -o $@

obj/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(HOST_FLAGS) -MMD -MP -c $< -o $@

obj/kernel%.o: port/kernel%.mu port/opencl_device_compat.h port/vina_musa_launch.h $(wildcard src/OpenCL/src/kernels/*.cl) src/OpenCL/inc/kernel2.h
	@mkdir -p $(@D)
	$(SDK)/bin/mcc $(DEVICE_FLAGS) -c $< -o $@

obj/kernel1_reference.o: reference/grid_fallback.mu reference/kernel1.cl port/opencl_device_compat.h port/vina_musa_launch.h src/OpenCL/inc/kernel2.h
	@mkdir -p $(@D)
	$(SDK)/bin/mcc $(DEVICE_FLAGS) -Iport -Dvina_musa_launch_kernel1=vina_musa_launch_kernel1_reference -c $< -o $@

-include $(HOST_OBJS:.o=.d)
