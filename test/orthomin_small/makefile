# NAME OF OUTPUT PROGRAM
PROGRAM = orthomin_small

# COMPILER FLAGS
SRCEXTS = .cc .cpp

# IAJA LIBRARY
LIBIAJA = ../../build/libiajasparse.a
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
