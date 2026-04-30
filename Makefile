# -----------------------------------------------------------------
# opticx Makefile
# -----------------------------------------------------------------
# Author    : J. J. Esteve-Paredes (JJEP)
# Modified  : D. Hernangómez-Pérez (DH)
# Version   : 1.1
# Date      : 20.06.2025
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Compiler and flags
# -----------------------------------------------------------------
FC     = gfortran
#FFLAGS = -g -fcheck=all -Wall -O0
FFLAGS = -O2

# -----------------------------------------------------------------
# Optional: use MKL  
# Usage: make USE_MKL=1 
# Note: only tested with Ubuntu, install with sudo apt-get install intel-mkl
# -----------------------------------------------------------------
ifeq ($(USE_MKL),1)
LIBS = -lmkl_rt -fopenmp -lpthread -lm -ldl	
# MKLROOT ?= /opt/intel/oneapi/mkl/latest
# LIBS = -L$(MKLROOT)/lib/intel64 \
#       -Wl,--start-group \
#       -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread \
#       -Wl,--end-group -fopenmp -lpthread -lm -ldl
else
LIBS   = -lopenblas -fopenmp -lgfortran 
endif

# -----------------------------------------------------------------
# Directories
# -----------------------------------------------------------------
MAINDIR  = main
SRCDIR   = src
BINDIR   = bin
BUILDDIR = build

# -----------------------------------------------------------------
# Files
# -----------------------------------------------------------------
SRC_MODULES = $(wildcard $(SRCDIR)/*.f90)
OBJ_MODULES = $(SRC_MODULES:$(SRCDIR)/%.f90=$(BUILDDIR)/%.o)
SRC_MAIN    = $(MAINDIR)/opticx.f90
OBJ_MAIN    = $(BINDIR)/opticx.o
TARGET      = $(BINDIR)/opticx

# -----------------------------------------------------------------
# Build target
# -----------------------------------------------------------------
all: $(TARGET)
	rm -f $(OBJ_MAIN)

# -----------------------------------------------------------------
# Directory creation rules
# -----------------------------------------------------------------
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# -----------------------------------------------------------------
# Linking rules
# -----------------------------------------------------------------
# Use the -J flag to specify the directory for .mod files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# -----------------------------------------------------------------------------
# Compilation rules
# -----------------------------------------------------------------------------
# Rule for creating the executable
$(TARGET): $(OBJ_MODULES) $(OBJ_MAIN)
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(OBJ_MAIN) -o $(TARGET) $(LIBS)

$(BINDIR)/opticx.o: $(SRC_MAIN) | $(BINDIR) $(BUILDDIR)
	$(FC) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@


# -----------------------------------------------------------------
# Module dependencies
# -----------------------------------------------------------------
# DH: Add more dependencies as needed for other modules
#     $(BUILDDIR)/some_other_module.o: \
# 	  $(BUILDDIR)/dependency_module.o

$(BUILDDIR)/parser_wannier90_tb.o: \
	$(BUILDDIR)/parser_input_file.o 

$(BUILDDIR)/parser_optics_xatu_dim.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_input_file.o 

$(BUILDDIR)/exciton_envelopes.o: \
    $(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o

$(BUILDDIR)/bands.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o

$(BUILDDIR)/ome_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o  

$(BUILDDIR)/ome_ex.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/exciton_envelopes.o 

$(BUILDDIR)/ome.o: \
    $(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/ome_sp.o \
	$(BUILDDIR)/ome_ex.o

$(BUILDDIR)/optical_response.o: \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/sigma_first_sp.o \
	$(BUILDDIR)/sigma_first_ex.o \
	$(BUILDDIR)/sigma_second_sp.o \
	$(BUILDDIR)/sigma_second_ex.o

$(BUILDDIR)/sigma_first_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/ome.o 

$(BUILDDIR)/sigma_first_ex.o: \
	$(BUILDDIR)/constants_math.o 

$(BUILDDIR)/sigma_second_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/ome.o 

$(BUILDDIR)/sigma_second_ex.o: \
	$(BUILDDIR)/constants_math.o \
    $(BUILDDIR)/ome_ex.o \
	$(BUILDDIR)/sigma_second_sp.o 
# -----------------------------------------------------------------
# Clean
# -----------------------------------------------------------------
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(BINDIR)/opticx.o $(TARGET)
