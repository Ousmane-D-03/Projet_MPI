# Makefile pour ARN - VERSION CORRIGÉE avec Needleman
CXX = mpic++

# ⚠️ IMPORTANT : Ajouter -fopenmp dans CXXFLAGS ET LDFLAGS
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -fopenmp -DWITH_GRAPHVIZ
LDFLAGS = -fopenmp -lcgraph

# Pour activer Needleman, décommenter cette ligne :
# CXXFLAGS += -DUSE_NEEDLEMAN

# Fichiers source de base
SOURCES = ARNSequence.cpp main_arn.cpp \
          ../PAM/PAM.cpp \
          ../Floyd/FoydPar.cpp ../Floyd/ForGraph.cpp ../Floyd/Utils.cpp

# Si USE_NEEDLEMAN est défini, ajouter Needleman.cpp
ifeq ($(findstring USE_NEEDLEMAN,$(CXXFLAGS)),USE_NEEDLEMAN)
    SOURCES += ../Needleman/Needleman.cpp
    $(info ✅ Needleman activé)
endif

OBJECTS = $(SOURCES:.cpp=.o)
TARGET = arn_main

# Cible principale
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✅ Compilé: $(TARGET)"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Test simple
test: $(TARGET)
	@echo "\n=== Test ARN (4 procs) ==="
	@echo ">seq1\nACGTACGTACGTACGT\n>seq2\nACGTACGTACGTACGT\n>seq3\nTGCATGCATGCATGCA\n>seq4\nACGTACGTACGTACGT" > test.fasta
	mpirun -np 4 ./$(TARGET) test.fasta 10 2 test.dot

# Cible pour compiler avec Needleman
needleman: CXXFLAGS += -DUSE_NEEDLEMAN
needleman: clean all

clean:
	rm -f $(OBJECTS) $(TARGET) test.fasta test.dot *.dot ../Needleman/Needleman.o

.PHONY: all clean test needleman