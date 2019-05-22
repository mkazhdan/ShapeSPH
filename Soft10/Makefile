TARGET = Soft10
SOURCE_FILES = cospmls.cpp OURmods.cpp soft_fftw.cpp csecond.cpp OURperms.cpp soft_fftw_pc.cpp FFTcode.cpp permroots.cpp soft_fftw_wo.cpp fft_grids.cpp primitive.cpp soft_sym.cpp fft_grids_so3.cpp primitive_FST.cpp utils_so3.cpp FST_semi_memo.cpp rotate_so3.cpp utils_vec_cx.cpp FST_semi_memo_fftw.cpp rotate_so3_mem.cpp weights.cpp indextables.cpp seminaive.cpp wignerTransforms.cpp makeWigner.cpp seminaive_fftw.cpp wignerTransforms_fftw.cpp naive_synthesis.cpp so3_correlate_fftw.cpp wignerTransforms_sym.cpp newFCT.cpp so3_correlate_sym.cpp oddweights.cpp soft.cpp
SRC_DIR = ./
BIN_DIR = ../Bin/
OBJ_DIR = ../Objects/
PLATFORM = Linux

TARGET_LIB = lib$(TARGET).a

CFLAGS += -I. -I../Include -I/usr/include/ -std=c++14 -Wunused-result -O3 -DRELEASE -funroll-loops -ffast-math -DNDEBUG

CC  = gcc
CXX = g++
MD  = mkdir
AR  = ar

OBJECTS=$(addprefix $(OBJ_DIR)$(PLATFORM)/$(TARGET)/, $(addsuffix .o, $(basename $(SOURCE_FILES))))

all: $(BIN_DIR)$(PLATFORM)/
all: $(OBJ_DIR)$(PLATFORM)/$(TARGET)/
all: $(BIN_DIR)$(PLATFORM)/$(TARGET_LIB)

clean:
	rm -f $(BIN_DIR)$(PLATFORM)/$(TARGET_LIB)
	rm -f $(OBJECTS)

$(BIN_DIR)$(PLATFORM)/:
	$(MD) -p $(BIN_DIR)
	$(MD) -p $(BIN_DIR)$(PLATFORM)

$(OBJ_DIR)$(PLATFORM)/$(TARGET)/:
	$(MD) -p $(OBJ_DIR)
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/$(TARGET)/

$(BIN_DIR)$(PLATFORM)/$(TARGET_LIB): $(OBJECTS)
	$(AR) rcs $@ $(OBJECTS)

$(OBJ_DIR)$(PLATFORM)/$(TARGET)/%.o: $(SRC_DIR)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) $<
