# IAJA LIBRARY
LIBIAJA = ../../lib/libiajasparse.a
CXXFLAGS += -I../../include/

NONLIN = img_denoise_nonlinear
LIN = img_denoise_linear

.PHONY: all clean linear nonlinear

all: linear nonlinear

linear:
	$(CXX) $(CXXFLAGS)  $(LIN).cc $(LIBIAJA) -o $(LIN).x 

nonlinear:
	$(CXX) $(CXXFLAGS) -o $(NONLIN).x  $(NONLIN).cc $(LIBIAJA) 

clean:
	-$(RM) *.x *.o img.m *.txt vgcore.* octave-workspace*
