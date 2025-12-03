CXX = g++
CFLAGS = -O2 -std=c++11 -DWITH_GRAPHVIZ

INCLUDES = -I../Floyd
LDFLAGS = -lcgraph -lgvc  

SRCS = main_pam.cpp PAM.cpp ../Floyd/ForGraph.cpp 
OBJS = $(SRCS:.cpp=.o)

all: pam mpipam

pam: $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

# build with MPI (if mpic++ available)
mpipam:
	mpic++ -DUSE_MPI -DWITH_GRAPHVIZ $(CFLAGS) $(INCLUDES) -o pam_mpi $(SRCS) -lcgraph


test: pam mpipam
	@echo "Test PAM séquentiel"
	./pam ../Exemple2.dot 3
	@echo "---------------------------------"
	@echo "Test PAM MPI"
	mpirun -np 4 ./pam_mpi ../Exemple2.dot 3

clean:
	rm -f $(OBJS) pam pam_mpi

.PHONY: all clean mpipam


