TARGET1 = ShapeAlign
TARGET2 = ShapeDescriptor
TARGET3 = ShapeSymmetry
SOURCE_FILES1 = $(TARGET1)/ShapeAlign.cpp
SOURCE_FILES2 = $(TARGET2)/ShapeDescriptor.cpp
SOURCE_FILES3 = $(TARGET3)/ShapeSymmetry.cpp

SRC_DIR = ./
BIN_DIR = Bin/
OBJ_DIR = Objects/
PLATFORM = Linux

CFLAGS += -fopenmp -I. -IInclude -std=c++14 -Wunused-result -O3 -DRELEASE -funroll-loops -ffast-math -DNDEBUG -g
LFLAGS += -lgomp -L. -L$(BIN_DIR)$(PLATFORM) -lUtil -lSoft10 -ljpeg -O3 -lfftw3 -lfftw3f

CC  = gcc
CXX = g++
MD  = mkdir
AR  = ar

OBJECTS1=$(addprefix $(OBJ_DIR)$(PLATFORM)/, $(addsuffix .o, $(basename $(SOURCE_FILES1))))
OBJECTS2=$(addprefix $(OBJ_DIR)$(PLATFORM)/, $(addsuffix .o, $(basename $(SOURCE_FILES2))))
OBJECTS3=$(addprefix $(OBJ_DIR)$(PLATFORM)/, $(addsuffix .o, $(basename $(SOURCE_FILES3))))

all: $(BIN_DIR)$(PLATFORM)/
all: $(OBJ_DIR)$(PLATFORM)/$(TARGET1)/
all: $(OBJ_DIR)$(PLATFORM)/$(TARGET2)/
all: $(OBJ_DIR)$(PLATFORM)/$(TARGET3)/
all: $(BIN_DIR)$(PLATFORM)/$(TARGET1)
all: $(BIN_DIR)$(PLATFORM)/$(TARGET3)
all: $(BIN_DIR)$(PLATFORM)/$(TARGET2)

clean:
	rm -f $(BIN_DIR)$(PLATFORM)/$(TARGET1)
	rm -f $(BIN_DIR)$(PLATFORM)/$(TARGET2)
	rm -f $(BIN_DIR)$(PLATFORM)/$(TARGET3)
	rm -f $(OBJECTS1)
	rm -f $(OBJECTS2)
	rm -f $(OBJECTS3)
	make clean -C Util
	make clean -C Soft10

$(BIN_DIR)$(PLATFORM)/:
	$(MD) -p $(BIN_DIR)
	$(MD) -p $(BIN_DIR)$(PLATFORM)

$(OBJ_DIR)$(PLATFORM)/$(TARGET1)/:
	$(MD) -p $(OBJ_DIR)
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/$(TARGET1)/

$(OBJ_DIR)$(PLATFORM)/$(TARGET2)/:
	$(MD) -p $(OBJ_DIR)
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/$(TARGET2)/

$(OBJ_DIR)$(PLATFORM)/$(TARGET3)/:
	$(MD) -p $(OBJ_DIR)
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/
	$(MD) -p $(OBJ_DIR)$(PLATFORM)/$(TARGET3)/

$(BIN_DIR)$(PLATFORM)/$(TARGET1): $(OBJECTS1)
	make -C Util
	make -C Soft10
	$(CXX) -o $@ $(OBJECTS1) $(LFLAGS)

$(BIN_DIR)$(PLATFORM)/$(TARGET2): $(OBJECTS2)
	make -C Util
	make -C Soft10
	$(CXX) -o $@ $(OBJECTS2) $(LFLAGS)

$(BIN_DIR)$(PLATFORM)/$(TARGET3): $(OBJECTS3)
	make -C Util
	make -C Soft10
	$(CXX) -o $@ $(OBJECTS3) $(LFLAGS)
	
$(OBJ_DIR)$(PLATFORM)/$(TARGET1)/%.o: $(SRC_DIR)$(TARGET1)/%.cpp
	$(CXX) -c -o $@ $(CFLAGS) $<

$(OBJ_DIR)$(PLATFORM)/$(TARGET2)/%.o: $(SRC_DIR)$(TARGET2)/%.cpp
	$(CXX) -c -o $@ $(CFLAGS) $<

$(OBJ_DIR)$(PLATFORM)/$(TARGET3)/%.o: $(SRC_DIR)$(TARGET3)/%.cpp
	$(CXX) -c -o $@ $(CFLAGS) $<
