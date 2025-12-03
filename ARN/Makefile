# Makefile pour le traitement des séquences d'ARN
CXX = mpic++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11
LDFLAGS =

# Fichiers source
SOURCES = ARNSequence.cpp main_arn.cpp ../PAM/PAM.cpp ../Floyd/FoydPar.cpp
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = arn_main

# Cibles principales
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Compilation réussie : $(TARGET)"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Tests
test: $(TARGET)
	@echo "Création de fichier FASTA de test..."
	@echo ">seq1" > test_sequences.fasta
	@echo "ACGTACGTACGTACGT" >> test_sequences.fasta
	@echo ">seq2" >> test_sequences.fasta
	@echo "ACGTACGTACGTACGT" >> test_sequences.fasta
	@echo ">seq3" >> test_sequences.fasta
	@echo "TGCATGCATGCATGCA" >> test_sequences.fasta
	@echo "Test avec epsilon = 10"
	./$(TARGET) test_sequences.fasta 10 3 test_graph.dot

# Nettoyage
clean:
	rm -f $(OBJECTS) $(TARGET) test_sequences.fasta test_graph.dot

.PHONY: all clean test
