# Makefile to compile the read_header program

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17

# Target executable
TARGET = reader

# Source files
SRCS = read_header.cpp
HEADERS = read_header.hpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule to compile the program
all: $(TARGET)

# Rule to compile the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)	

# Rule to compile the .cpp files into .o (object) files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) $(TARGET)
