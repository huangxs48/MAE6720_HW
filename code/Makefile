CXX = g++
CXXFLAGS = -O -std=c++11 -I /usr/local/include
LDFLAGS =
LDLIBS = 

EXE_DIR := bin/
EXECUTABLE := $(EXE_DIR)main
SRC_FILES := $(wildcard *.cpp) \
             $(wildcard MatrixSolver/*.cpp) \
             $(wildcard PDE/*.cpp) \
             $(wildcard Mesh/*.cpp) \
             $(wildcard PDESolver/*.cpp) \
             $(wildcard Output/*.cpp) \
             $(wildcard BoundaryCondition/*.cpp) \

OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_PREFIX := src/
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

.PHONY : all dirs clean

all : dirs $(EXECUTABLE)

objs : dirs $(OBJ_FILES)

dirs : $(EXE_DIR) $(OBJ_DIR)


# Link objects into executable

$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $(OBJ_FILES) $(LDFLAGS) $(LDLIBS)

# Create objects from source files

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# Cleanup

clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
