# =====================================
# BASIC DEFS
# =====================================

# OUPUT PROGRAM AND SOURCE DIRECTORIES
PROGRAM = libiajasparse.a
SRCDIRS = ../src
SRCEXTS = .c .cc

# COMPILER AND FLAGS
ARFLAGS = rsv
CXXFLAGS += -I../include


# =====================================
# BUILD PROCESS
# =====================================
SOURCES = $(foreach d, $(SRCDIRS), $(wildcard $(addprefix $(d)/*, $(SRCEXTS))))

OBJS = $(foreach x, $(SRCEXTS), \
	   $(patsubst %$(x), %.o, $(notdir $(filter %$(x),$(SOURCES)))))

# this is fine
DEPS := $(patsubst %.o,%.d,$(OBJS))

# commands to be executed for each build process
.PHONY = cleanobj all clean objs rebuild

all: $(PROGRAM)

# GENERATE DEPENDENCIES
%.d : $(SRCDIRS)/%.c
	$(CC) -MM -MD $(CFLAGS) $<
%.d : $(SRCDIRS)/%.cc
	$(CXX) -MM -MD $(CXXFLAGS) $<

# COMPILE OBJECTS
objs: $(OBJS)
%.o : $(SRCDIRS)/%.c
	$(CC)  -c -fPIC $(CFLAGS) $< -o $@
%.o : $(SRCDIRS)/%.cc
	$(CXX) -c -fPIC $(CXXFLAGS) $< -o $@


# PRODUCE PROGRAM
$(PROGRAM) : $(OBJS)
	$(AR) $(ARFLAGS) $(PROGRAM) $(OBJS)

-include $(DEPS)


# ifeq ($(strip $(SRCEXTS)), .c)  # C file
	# $(CC) -o $(PROGRAM) $(OBJS)
# else # CPP file
	# $(CXX) -o $(PROGRAM) $(OBJS)
# endif
# =====================================
# CLEANUP
# =====================================
rebuild: clean all
clean:
	-$(RM) *.o *.d *.a
	
