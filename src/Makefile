ROOT_DIR := .
IO_DIR := $(ROOT_DIR)/io
UTILS_DIR := $(ROOT_DIR)/utils

EXEC   := interactions
TARGETSRC   := main.c
LIBSRC := tree_conversion.c interactions.c $(UTILS_DIR)/utils.c  $(IO_DIR)/lhalotree/read_lhalotree.c $(UTILS_DIR)/progressbar.c 
TARGETOBJS := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)
INCL   := tree_conversion.h interactions.h $(IO_DIR)/io_trees.h $(IO_DIR)/ctree/ctree.h $(IO_DIR)/lhalotree/read_lhalotree.h  $(IO_DIR)/lhalotree/lhalotree.h Makefile $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/utils.h 

INCLUDE := -I$(UTILS_DIR) -I$(IO_DIR)
#USE-MPI := yes  # set this if you want to run in embarrassingly parallel
#OPT += -DSINGLE_TREE_ATATIME
OPT += -DFIX_FLYBYS

ifdef USE-MPI
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC := mpicc  # sets the C-compiler
else
    CC := gcc  # sets the C-compiler
endif

# GSL automatic detection
GSL_FOUND := $(shell gsl-config --version 2>/dev/null)
ifndef GSL_FOUND
  $(warning GSL not found in path - please install GSL before installing SAGE (or, update the PATH environment variable such that "gsl-config" is found))
  # if the automatic detection fails, set GSL_DIR appropriately
  GSL_DIR := /opt/local
  GSL_INCL := -I$(GSL_DIR)/include  
  GSL_LIBDIR := $(GSL_DIR)/lib
  # since GSL is not in PATH, the runtime environment might not be setup correctly either
  # therefore, adding the compiletime library path is even more important (the -Xlinker bit)
  GSL_LIBS := -L$(GSL_LIBDIR) -lgsl -lgslcblas -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 
else
  # GSL is probably configured correctly, pick up the locations automatically
  GSL_INCL := $(shell gsl-config --cflags)
  GSL_LIBDIR := $(shell gsl-config --prefix)/lib
  GSL_LIBS   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
endif

OPTIMIZE := -g -std=gnu11 -O3 -march=native
CFLAGS := -Wextra -Wshadow -Wall -g -Wsign-compare -Wall -Wextra -Wshadow -Wunused -fPIC -D_POSIX_SOURCE=200809L -D_GNU_SOURCE -D_DARWIN_C_SOURCE 
CFLAGS += $(OPT) $(OPTIMIZE) $(GSL_INCL)

CLINK  :=   -g -lrt -lm  $(GSL_LIBS) 


OBJS:=$(TARGETOBJS) $(LIBOBJS)

all: $(EXEC) $(OBJS) Makefile $(LIBSRC) $(TARGETSRC)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(CLINK) $(OBJS) $(LIBS)   -o  $(EXEC)


%.o:%.c $(INCL)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean celna celan

clean:
	rm -f $(OBJS) $(EXEC)
