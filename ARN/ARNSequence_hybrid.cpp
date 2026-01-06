/**
 * @file ARNSequence_hybrid.cpp
 * @brief Version HYBRIDE MPI + OpenMP
 * À placer dans ARN/ARNSequence_hybrid.cpp (garde ARNSequence.cpp original)
 */

#include "ARNSequence_hybrid.hpp"
#include <fstream>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <cmath>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

int levenshteinDistance(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();
    
    if (m == 0) return n;
    if (n == 0) return m;
    
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    
    for (int i = 0; i <= m; i++) dp[i][0] = i;
    for (int j = 0; j <= n; j++) dp[0][j] = j;
    
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int cost = (seq1[i - 1] == seq2[j - 1]) ? 0 : 1;
            dp[i][j] = min({
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost
            });
        }
    }
    
    return dp[m][n];
}

int hammingDistance(const string& seq1, const string& seq2) {
    if (seq1.length() != seq2.length()) return -1;
    
    int distance = 0;
    for (size_t i = 0; i < seq1.length(); i++) {
        if (seq1[i] != seq2[i]) distance++;
    }
    return distance;
}

int readFASTAFile(const string& filename, vector<ARNSeq>& sequences) {
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Erreur : impossible d'ouvrir " << filename << endl;
        return -1;
    }
    
    string line, currentLabel, currentSequence;
    int seqId = 0;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!currentSequence.empty()) {
                ARNSeq seq;
                seq.id = seqId++;
                seq.sequence = currentSequence;
                seq.label = currentLabel;
                sequences.push_back(seq);
                currentSequence.clear();
            }
            currentLabel = line.substr(1);
        } else {
            currentSequence += line;
        }
    }
    
    if (!currentSequence.empty()) {
        ARNSeq seq;
        seq.id = seqId++;
        seq.sequence = currentSequence;
        seq.label = currentLabel;
        sequences.push_back(seq);
    }
    
    file.close();
    cout << "Fichier FASTA lu : " << sequences.size() << " séquence(s)" << endl;
    return sequences.size();
}

int* computeDistanceMatrix(const vector<ARNSeq>& sequences, 
                           int (*distanceFunc)(const string&, const string&)) {
    int n = sequences.size();
    int* distanceMatrix = new int[n * n];
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                distanceMatrix[i * n + j] = 0;
            } else if (i < j) {
                int dist = distanceFunc(sequences[i].sequence, sequences[j].sequence);
                distanceMatrix[i * n + j] = dist;
                distanceMatrix[j * n + i] = dist;
            }
        }
    }
    
    return distanceMatrix;
}

#ifdef USE_MPI
/**
 * VERSION HYBRIDE : MPI distribue paires + OpenMP parallélise calculs
 */
int* computeDistanceMatrix_Hybrid(const vector<ARNSeq>& sequences, 
                                  int (*distanceFunc)(const string&, const string&),
                                  int rank, int nprocs) {
    int n = sequences.size();
    int total_pairs = (n * (n - 1)) / 2;
    
    int pairs_per_proc = total_pairs / nprocs;
    int remainder = total_pairs % nprocs;
    
    int my_start = rank * pairs_per_proc + min(rank, remainder);
    int my_end = my_start + pairs_per_proc + (rank < remainder ? 1 : 0);
    int my_count = my_end - my_start;
    
    auto pair_to_ij = [n](int pair_id) -> pair<int,int> {
        int i = (int)((2*n - 1 - sqrt((2*n-1)*(2*n-1) - 8*pair_id)) / 2);
        int j = pair_id - (i*n - i*(i+1)/2) + i + 1;
        return {i, j};
    };
    
    vector<int> local_distances(my_count);
    vector<pair<int,int>> local_indices(my_count);
    
    // ====== OPENMP : Calcul des distances locales ======
    #pragma omp parallel for schedule(dynamic, 32)
    for(int local_id = 0; local_id < my_count; ++local_id) {
        int pair_id = my_start + local_id;
        auto [i, j] = pair_to_ij(pair_id);
        
        int dist = distanceFunc(sequences[i].sequence, sequences[j].sequence);
        
        local_distances[local_id] = dist;
        local_indices[local_id] = {i, j};
    }
    
    // ====== MPI : Rassemblement ======
    int* distanceMatrix = nullptr;
    
    if(rank == 0) {
        distanceMatrix = new int[n * n];
        
        #pragma omp parallel for
        for(int i = 0; i < n; ++i) {
            distanceMatrix[i * n + i] = 0;
        }
    }
    
    vector<int> recvcounts(nprocs), displs(nprocs);
    MPI_Gather(&my_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        displs[0] = 0;
        for(int p = 1; p < nprocs; ++p) {
            displs[p] = displs[p-1] + recvcounts[p-1];
        }
    }
    
    vector<int> all_distances;
    if(rank == 0) all_distances.resize(total_pairs);
    
    MPI_Gatherv(local_distances.data(), my_count, MPI_INT,
                all_distances.data(), recvcounts.data(), displs.data(),
                MPI_INT, 0, MPI_COMM_WORLD);
    
    struct PairInt { int i; int j; };
    vector<PairInt> local_pairs(my_count);
    for(int k = 0; k < my_count; ++k) {
        local_pairs[k] = {local_indices[k].first, local_indices[k].second};
    }
    
    vector<PairInt> all_pairs;
    if(rank == 0) all_pairs.resize(total_pairs);
    
    MPI_Datatype MPI_PAIR_INT;
    MPI_Type_contiguous(2, MPI_INT, &MPI_PAIR_INT);
    MPI_Type_commit(&MPI_PAIR_INT);
    
    MPI_Gatherv(local_pairs.data(), my_count, MPI_PAIR_INT,
                all_pairs.data(), recvcounts.data(), displs.data(),
                MPI_PAIR_INT, 0, MPI_COMM_WORLD);
    
    MPI_Type_free(&MPI_PAIR_INT);
    
    if(rank == 0) {
        #pragma omp parallel for
        for(int k = 0; k < total_pairs; ++k) {
            int i = all_pairs[k].i;
            int j = all_pairs[k].j;
            int dist = all_distances[k];
            
            distanceMatrix[i * n + j] = dist;
            distanceMatrix[j * n + i] = dist;
        }
    }
    
    return distanceMatrix;
}
#endif

int writeGraphDOT(const vector<ARNSeq>& sequences, int* distanceMatrix, 
                  int epsilon, const string& outputFile) {
    ofstream file(outputFile);
    
    if (!file.is_open()) {
        cerr << "Erreur : impossible d'ouvrir " << outputFile << endl;
        return -1;
    }
    
    int n = sequences.size();
    
    file << "graph ARN {" << endl;
    file << "  rankdir=LR;" << endl;
    
    for (int i = 0; i < n; i++) {
        file << "  seq" << i << " [label=\"" << sequences[i].label << "\"];" << endl;
    }
    
    file << endl;
    
    int edgeCount = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int dist = distanceMatrix[i * n + j];
            if (dist < epsilon) {
                file << "  seq" << i << " -- seq" << j 
                     << " [weight=" << dist << ", label=\"" << dist << "\"];" << endl;
                edgeCount++;
            }
        }
    }
    
    file << "}" << endl;
    file.close();
    
    cout << "Graphe écrit : " << outputFile << endl;
    cout << "  Nœuds: " << n << ", Arêtes: " << edgeCount << endl;
    
    return 0;
}

void printARNSeq(const ARNSeq& seq) {
    cout << "ID: " << seq.id << " | Label: " << seq.label 
         << " | Séquence: " << seq.sequence.substr(0, 50);
    if (seq.sequence.length() > 50) cout << "...";
    cout << " (taille: " << seq.sequence.length() << ")" << endl;
}

void printDistanceMatrix(int* distanceMatrix, int n) {
    cout << "Matrice de distances (" << n << "x" << n << "):" << endl;
    
    cout << "    ";
    for (int j = 0; j < n; j++) cout << setw(6) << j;
    cout << endl;
    
    for (int i = 0; i < n; i++) {
        cout << setw(3) << i << " ";
        for (int j = 0; j < n; j++) {
            cout << setw(6) << distanceMatrix[i * n + j];
        }
        cout << endl;
    }
}
