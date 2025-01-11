# programs
TARGET := gencore
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

# compiler
GXX := gcc
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic

# object files that need lcptools
LCPTOOLS_CXXFLAGS := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LDFLAGS := -L$(CURRENT_DIR)/lcptools/lib -llcptools -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib -lz
HTSLIB_CXXFLAGS := -I$(CURRENT_DIR)/htslib/include
HTSLIB_LDFLAGS := -L$(CURRENT_DIR)/htslib/lib -lhts -Wl,-rpath,$(CURRENT_DIR)/htslib/lib -pthread

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) $(HTSLIB_LDFLAGS) -lm
	rm *.o

gencore.o: gencore.c
	$(GXX) $(CXXFLAGS) $(HTSLIB_CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

rfasta.o: rfasta.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

rfastq.o: rfastq.c
	$(GXX) $(CXXFLAGS) $(HTSLIB_CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

rload.o: rload.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

utils.o: utils.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@ -lz

init.o: init.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

%.o: %.c
	$(GXX) $(CXXFLAGS) -c $< -o $@

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

reinstall-lcptools:
	@echo "Re-installing lcptools"
	git submodule deinit -f -- lcptools
	rm -rf lcptools
	git submodule update --init --recursive

recompile-lcptools:
	cd lcptools && \
	make uninstall PREFIX=$(CURRENT_DIR)/lcptools && \
	make install PREFIX=$(CURRENT_DIR)/lcptools
