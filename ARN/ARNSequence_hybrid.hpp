/**
 * @file ARNSequence_hybrid.hpp
 * @brief Version HYBRIDE avec fonction computeDistanceMatrix_Hybrid
 * À placer dans ARN/ARNSequence_hybrid.hpp (garde ARNSequence.hpp original)
 */

#ifndef ARNSEQUENCE_HYBRID_HPP
#define ARNSEQUENCE_HYBRID_HPP

#include <string>
#include <vector>
#include <iostream>

using namespace std;

struct ARNSeq {
    int id;
    string sequence;
    string label;
};

int levenshteinDistance(const string& seq1, const string& seq2);
int hammingDistance(const string& seq1, const string& seq2);
int readFASTAFile(const string& filename, vector<ARNSeq>& sequences);

int* computeDistanceMatrix(const vector<ARNSeq>& sequences, 
                           int (*distanceFunc)(const string&, const string&));

#ifdef USE_MPI
/**
 * @brief VERSION HYBRIDE : MPI + OpenMP
 */
int* computeDistanceMatrix_Hybrid(const vector<ARNSeq>& sequences, 
                                  int (*distanceFunc)(const string&, const string&),
                                  int rank, int nprocs);
#endif

int writeGraphDOT(const vector<ARNSeq>& sequences, int* distanceMatrix, 
                  int epsilon, const string& outputFile);

void printARNSeq(const ARNSeq& seq);
void printDistanceMatrix(int* distanceMatrix, int n);

#endif // ARNSEQUENCE_HYBRID_HPP
