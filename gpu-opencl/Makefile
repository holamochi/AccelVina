# BASE_BUILD is the directory produced by ../build_baseline.sh.
BASE_BUILD ?=
ifeq ($(strip $(BASE_BUILD)),)
$(error Set BASE_BUILD to a verified baseline build directory)
endif
BUILD_DIR ?= build
SRC_ROOT ?= .
CXX := g++
CPPFLAGS := -I$(BASE_BUILD)/deps/boost-headers/include \
            -I$(BASE_BUILD)/deps/OpenCL-Headers-2024.10.24 \
            -I$(SRC_ROOT)/lib -I$(SRC_ROOT)/OpenCL/inc \
            -DOPENCL_3_0 -DNVIDIA_PLATFORM -DSMALL_BOX -DBOOST_TIMER_ENABLE_DEPRECATED -DNDEBUG
CXXFLAGS := -O3 -std=gnu++17
LDLIBS := -L$(BASE_BUILD)/runtime/lib -l:libboost_program_options.so.1.84.0 \
          -l:libboost_system.so.1.84.0 -l:libboost_filesystem.so.1.84.0 \
          -l:libboost_thread.so.1.84.0 -l:libOpenCL.so.1 -lstdc++fs -lm -lpthread

sources := $(patsubst $(SRC_ROOT)/%,%,$(wildcard $(SRC_ROOT)/lib/*.cpp)) main/main.cpp OpenCL/src/wrapcl.cpp
objects := $(addprefix $(BUILD_DIR)/obj/,$(sources:.cpp=.o))

.PHONY: all
all: $(BUILD_DIR)/AutoDock-Vina-GPU-2-1

$(BUILD_DIR)/AutoDock-Vina-GPU-2-1: $(objects)
	$(CXX) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/obj/%.o: $(SRC_ROOT)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(objects:.o=.d)
