CXX = g++
CFLAGS = -O2 -std=c++11

INCLUDES = -I../Floyd
LDFLAGS =

SRCS = main_pam.cpp PAM.cpp
OBJS = $(SRCS:.cpp=.o)

all: pam

pam: $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

# build with MPI (if mpic++ available)
mpipam:
	mpic++ -DUSE_MPI -DWITH_GRAPHVIZ $(CFLAGS) $(INCLUDES) -o pam_mpi $(SRCS) -lcgraph

clean:
	rm -f $(OBJS) pam pam_mpi

.PHONY: all clean mpipam
