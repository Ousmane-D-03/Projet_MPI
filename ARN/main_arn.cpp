// main_arn.cpp - Version CORRIGÉE
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "ARNSequence.hpp"

#ifdef USE_NEEDLEMAN
#include "../Needleman/Needleman.hpp"
#endif

#include "../PAM/PAM.hpp"
#include "../Floyd/FoydPar.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(argc < 4) {
        if(pid == 0) {
            cerr << "Usage: " << argv[0] 
                 << " <fichier_fasta> <epsilon> <k_clusters> [output_dot]" << endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    string fastaFile = argv[1];
    int epsilon = stoi(argv[2]);
    int k_clusters = stoi(argv[3]);
    string outputFile = (argc >= 5) ? argv[4] : "arn_graph.dot";

    double t_total_start = MPI_Wtime();

    vector<ARNSeq> sequences;
    int nbSeq = 0;
    int* distanceMatrix = nullptr;

    // ===== ÉTAPE 1 : Lecture et calcul distances (rank 0) =====
    if(pid == 0) {
        nbSeq = readFASTAFile(fastaFile, sequences);
        if(nbSeq <= 0) {
            cerr << "Erreur : impossible de lire le fichier FASTA" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cout << "Séquences lues : " << nbSeq << endl;

        double t_dist_start = MPI_Wtime();
        
        // Calcul de la matrice de distances
#ifdef USE_NEEDLEMAN
        cout << "Utilisation de Needleman-Wunsch..." << endl;
        
        ScoringParams params(1, -1, -3, -1);
        distanceMatrix = new int[nbSeq * nbSeq];
        
        // Trouver score max pour normalisation
        int max_score = nbSeq * params.match; // Score si séquences identiques
        
        for(int i = 0; i < nbSeq; ++i) {
            distanceMatrix[i*nbSeq + i] = 0; // Diagonale
            
            for(int j = i+1; j < nbSeq; ++j) {
                // Calculer score
                int score = needleman_wunsch_sequential(
                    sequences[i].sequence,
                    sequences[j].sequence,
                    params
                );
                
                // CONVERSION : Score → Distance
                // Option 1 : Distance = max_score - score
                int distance = max_score - score;
                
                // Option 2 : Si score négatif, prendre valeur absolue
                // int distance = abs(score);
                
                distanceMatrix[i*nbSeq + j] = distance;
                distanceMatrix[j*nbSeq + i] = distance; // Symétrie
            }
        }
#else
        cout << "Utilisation de Levenshtein..." << endl;
        distanceMatrix = computeDistanceMatrix(sequences, levenshteinDistance);
#endif
        
        double t_dist_end = MPI_Wtime();
        cout << "Temps calcul distances : " << (t_dist_end - t_dist_start) << " sec" << endl;
    }

    // ===== ÉTAPE 2 : Vérifications pour Floyd =====
    int p_sqrt = static_cast<int>(sqrt(nprocs));
    if(p_sqrt * p_sqrt != nprocs) {
        if(pid == 0) {
            cerr << "Erreur : nprocs doit être un carré parfait (1,4,9,16,...)" << endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    // CORRECTION 1 : Broadcaster nbSeq D'ABORD
    MPI_Bcast(&nbSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // CORRECTION 2 : Calculer et broadcaster block_size CORRECTEMENT
    int block_size = 0;
    if(pid == 0) {
        if(nbSeq % p_sqrt != 0) {
            cerr << "Erreur : nbSeq (" << nbSeq << ") doit être divisible par sqrt(nprocs) (" 
                 << p_sqrt << ")" << endl;
            MPI_Finalize();
            return EXIT_FAILURE;
        }
        block_size = nbSeq / p_sqrt;
    }
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Maintenant on peut allouer avec block_size correct
    int* D_local = new int[block_size * block_size];

    // ===== ÉTAPE 3 : Floyd-Warshall =====
    if(pid == 0) cout << "\nCalcul Floyd-Warshall..." << endl;
    
    double t_floyd_start = MPI_Wtime();
    
    decouperMatrice(distanceMatrix, D_local, nbSeq, block_size, p_sqrt, 0, pid);
    int* D_global = floydBlocsHybrid(D_local, nbSeq, p_sqrt, pid, 0);
    
    double t_floyd_end = MPI_Wtime();

    // ===== ÉTAPE 4 : PAM et écriture graphe (rank 0) =====
    if(pid == 0) {
        cout << "Temps Floyd : " << (t_floyd_end - t_floyd_start) << " sec" << endl;
        cout << "\nClustering PAM (k=" << k_clusters << ")..." << endl;

        double t_pam_start = MPI_Wtime();
        
        pam::Result res = pam::pam_sequential(
            nbSeq, 
            vector<int>(D_global, D_global + nbSeq*nbSeq), 
            k_clusters, 
            42
        );
        
        double t_pam_end = MPI_Wtime();
        
        cout << "Temps PAM : " << (t_pam_end - t_pam_start) << " sec" << endl;
        
        cout << "\nRésultats PAM :" << endl;
        cout << "  Coût : " << res.cost << endl;
        cout << "  Médoïdes : ";
        for(auto m : res.medoids) cout << m << " ";
        cout << endl;
        
        // Compter points par cluster
        vector<int> cluster_counts(k_clusters, 0);
        for(int i = 0; i < nbSeq; ++i) {
            if(res.membership[i] >= 0 && res.membership[i] < k_clusters) {
                cluster_counts[res.membership[i]]++;
            }
        }
        cout << "  Points/cluster : ";
        for(auto c : cluster_counts) cout << c << " ";
        cout << endl;

        // Écriture du graphe DOT
        cout << "\nÉcriture graphe..." << endl;
        if(writeGraphDOT(sequences, D_global, epsilon, outputFile) == 0) {
            cout << "Graphe écrit : " << outputFile << endl;
        }
        
        double t_total_end = MPI_Wtime();
        cout << "\n=== TEMPS TOTAL : " << (t_total_end - t_total_start) << " sec ===" << endl;

        delete[] distanceMatrix;
        delete[] D_global;
    }

    delete[] D_local;
    MPI_Finalize();
    return 0;
}