# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2
LDFLAGS = -lhdf5_cpp -lhdf5

# Source files
SRCS = VoroLite.cpp Mesh.cpp Rays.cpp
OBJS = $(SRCS:.cpp=.o)

# Executable
TARGET = vorolite

# Default target
all: $(TARGET)

# Link object files into the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Optional: Show files being built
.PHONY: all clean
