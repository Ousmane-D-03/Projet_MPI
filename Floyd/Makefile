CXX = mpic++
CFLAGS = -O2 -std=c++11 -DUSE_MPI -DWITH_GRAPHVIZ -fopenmp
INCLUDES = -I.
LDFLAGS = -lcgraph -fopenmp   

SRCS = FoydPar.cpp ForGraph.cpp Utils.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = mpi_floyd

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
