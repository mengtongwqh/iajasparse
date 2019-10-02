# NAME OF OUTPUT PROGRAM
PROGRAM = andre_std_tests

# COMPILER FLAGS
SRCEXTS = .cc .cpp

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
	./$(PROGRAM).x ichol pcg 20 1
memcheck:
	valgrind ./$(PROGRAM).x ichol pcg 20 1

clean:
	-$(RM) *.x *.o vgcore.* *.txt *.timer
