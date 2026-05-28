# ========================================
#  Project Configuration
# ========================================
TARGET 		:= gencore
SRCS 		:= $(wildcard *.c)
OBJS 		:= $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

# ========================================
#  Compiler and Flags
# ========================================

CC         := gcc
CFLAGS     := -O3 -Wall -Wextra -Wpedantic -D_GNU_SOURCE
LDFLAGS    := -lm -pthread -lz

# ========================================
#  NUMA
# ========================================

NUMA_AVAILABLE 	:= 0
NUMA_INC 		:=
NUMA_LIB 		:=

ifneq ($(shell pkg-config --exists libnuma && echo yes),)
		NUMA_AVAILABLE := 1
		NUMA_INC := $(shell pkg-config --cflags libnuma)
		NUMA_LIB := $(shell pkg-config --libs libnuma)
	else
	ifneq ($(shell printf '\#include <numa.h>\n' | $(CC) -E - >/dev/null 2>&1 && echo yes),)
		NUMA_AVAILABLE := 1
		NUMA_LIB := -lnuma
	endif
endif

# Add macro and libs accordingly
CFLAGS   += -DNUMA_AVAILABLE=$(NUMA_AVAILABLE) $(NUMA_INC)
LDFLAGS  += $(NUMA_LIB)

# ========================================
#  External Libraries
# ========================================

# lcptools
LCPTOOLS_INC := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LIB := -L$(CURRENT_DIR)/lcptools/lib -llcptools -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib

# htslib
HTSLIB_INC := -I$(CURRENT_DIR)/htslib/include
HTSLIB_LIB := -L$(CURRENT_DIR)/htslib/lib -lhts -Wl,-rpath,$(CURRENT_DIR)/htslib/lib

# ========================================
#  Combined Flags
# ========================================

INCLUDES := $(LCPTOOLS_INC) $(HTSLIB_INC)
CXXLIBS  := $(LCPTOOLS_LIB) $(HTSLIB_LIB)

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(CXXLIBS) $(LDFLAGS)
	rm -f $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ========================================
#  Utility Targets
# ========================================

clean: 
	@echo "Cleaning"
	rm -f $(OBJS)
	rm -f $(TARGET)

install: clean install-htslib install-lcptools $(TARGET)

install-htslib:
	@echo "Installing htslib"
	cd htslib && \
	autoreconf -i && \
	./configure && \
	make && \
	make prefix=$(CURRENT_DIR)/htslib install

reinstall-htslib:
	@echo "Re-installing htslib"
	git submodule deinit -f -- htslib
	rm -rf htslib
	git submodule update --init --recursive

install-lcptools:
	@echo "Installing lcptool"
	cd lcptools && \
	make install PREFIX=$(CURRENT_DIR)/lcptools

recompile-lcptools:
	cd lcptools && \
	make uninstall PREFIX=$(CURRENT_DIR)/lcptools && \
	make install PREFIX=$(CURRENT_DIR)/lcptools
