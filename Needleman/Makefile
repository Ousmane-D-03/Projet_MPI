CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -fopenmp
LDFLAGS = -fopenmp

SOURCES = Needleman.cpp main_needleman.cpp
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = needleman_test

# <<< fichier utilisé automatiquement >>>
INPUT = data2000.fa

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(TARGET)
	./$(TARGET) $(INPUT)

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all test clean
