CXXFLAGS += -I../../include -g -Wall
SRCEXT = .cc

# NAME OF OUTPUT PROGRAM
TESTS = $(basename $(wildcard *$(SRCEXT)))

# IAJA LIBRARY
LIBIAJA = ../../lib/libiajasparse.a

.PHONY: all clean

all:
	@for prog in $(TESTS); \
		do \
		$(CXX) $$prog$(SRCEXT) $(LIBIAJA) $(CXXFLAGS) -o $$prog.x; \
		done
	@echo "done compiling all move tests"

clean:
	-$(RM) *.x *.o




