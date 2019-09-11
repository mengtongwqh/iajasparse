# NAME OF OUTPUT PROGRAM
PROGRAM = vector_test

# COMPILER FLAGS
SRCEXTS = .cc

# IAJA LIBRARY
LIBIAJA = ../../lib/libiajasparse.a
CXXFLAGS += -I../../include

SRC = $(wildcard *$(SRCEXTS))
OBJS = $(foreach x, $(SRCEXTS), \
	   $(patsubst %$(x), %.o, $(filter %$(x),$(SRC))))

.PHONY: run memcheck clean

all: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROGRAM).x $(OBJS) $(LIBIAJA)
objs: $(OBJS)

run: 
	./$(PROGRAM).x
memcheck:
	valgrind ./$(PROGRAM).x

clean:
	-$(RM) *.x *.o vgcore.*
