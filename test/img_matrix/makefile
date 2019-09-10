# IAJA LIBRARY
LIBIAJA = ../../lib/libiajasparse.a
NONLIN = img_matrix_nonlinear
LIN = img_matrix_linear
CXXFLAGS += -I../../include

SRC = img_matrix.cc
OBJS = img_matrix.o

.PHONY: run memcheck clean all

all: linear nonlinear

linear: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(LIN).x $< $(LIBIAJA) $(LIN).cc

nonlinear: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(NONLIN).x $< $(LIBIAJA) $(NONLIN).cc

objs: $(OBJS)

clean:
	-$(RM) *.x *.o img.m *.txt vgcore.* octave-workspace*
