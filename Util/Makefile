TARGET = Util
SOURCE_FILES = PlyFile.cpp Util.cpp
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
