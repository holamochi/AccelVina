SHELL := /bin/bash
BOOST_ROOT ?= /usr
BOOST_LIBDIR ?= $(BOOST_ROOT)/lib/x86_64-linux-gnu
CXX := g++
CC := gcc
BUILD ?= build
CPPFLAGS := -I$(BOOST_ROOT)/include -Isrc/lib -DNDEBUG -DVERSION='"v1.2.7-profiled"'
PROFILE_DIR ?= $(abspath profiles)
PROFILE_FLAGS ?= -fprofile-use=$(PROFILE_DIR) -Wmissing-profile
STRICT := -O3 -flto=4 $(PROFILE_FLAGS) -fprofile-prefix-path=$(abspath $(BUILD)) -march=icelake-server -mprefer-vector-width=256 -mno-fma -fno-fast-math -ffp-contract=off -pthread
CXXFLAGS := $(STRICT) -std=c++11 -Wno-deprecated-declarations
CFLAGS := $(STRICT) -std=c11
LIBS := $(addprefix $(BOOST_LIBDIR)/libboost_,$(addsuffix .a,program_options filesystem thread serialization system chrono atomic))
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(wildcard src/lib/*.cpp) src/main/main.cpp) $(BUILD)/lib/packed_bfgs.o $(BUILD)/lib/c_kernel.o
.PHONY: all
all: $(BUILD)/vina
$(BUILD)/vina: $(OBJECTS)
	$(CXX) $(STRICT) -o $@ $(OBJECTS) $(LIBS)
$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
$(BUILD)/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@
-include $(OBJECTS:.o=.d)
$(BUILD)/lib/c_kernel.o: CFLAGS += -mavx512f -mavx512dq -mavx512vl -mno-fma -ffp-contract=off

# GCC LTO may omit the profile of a folded static initializer. Require the
# measured file of every translation unit while keeping those warnings in the log.
ifneq ($(findstring -fprofile-use,$(PROFILE_FLAGS)),)
.PHONY: check-profiles
$(OBJECTS): | check-profiles
check-profiles:
	@set -eu; for file in $(patsubst $(BUILD)/%.o,%.gcda,$(OBJECTS)); do \
	    path="$(PROFILE_DIR)/$${file//\//#}"; test -s "$$path" || { echo "Missing profile: $$path" >&2; exit 1; }; \
	done
endif
