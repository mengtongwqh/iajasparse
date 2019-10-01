# COMPILER
export CC = gcc
export CXX = g++

# COMPILER FLAGS
FLAGS_DEBUG = -g -pedantic -Wall -Wextra -DDEBUG -fno-elide-constructors
FLAGS_RELEASE = -DNDEBUG -O2 -ffast-math -finline-functions
debug: CXXFLAGS   := $(CXXFLAGS) $(FLAGS_DEBUG)
debug: CFLAGS     := $(CFLAGS) $(FLAGS_DEBUG)
release: CXXFLAGS := $(CXXFLAGS) $(FLAGS_RELEASE)
release: CFLAGS   := $(CFLAGS) $(FLAGS_RELEASE)
export CXXFLAGS
export CFLAGS

# SUBDIRS TO BE BUILT
LIBDIR = ./lib
TESTDIR = ./test
TESTDIRS = $(wildcard $(TESTDIR)/*)
DIRS = $(LIBDIR) + $(TESTDIRS)

#TARGETS
all: clean debug

debug release:
	@echo "Compiling libiajasparse and test cases with flags: " $(CXXFLAGS)
	cd $(LIBDIR); $(MAKE);
	@for dir in $(DIRS); do cd $(PWD)/$$dir; $(MAKE); done
	@echo "Compilation completed."

clean:
	@for dir in $(DIRS); do cd $(PWD)/$$dir; $(MAKE) clean; done
