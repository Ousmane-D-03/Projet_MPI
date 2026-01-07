// main_arn_hybrid.cpp - Pipeline HYBRIDE MPI + OpenMP
// À placer dans ARN/main_arn_hybrid.cpp (garde main_arn.cpp original)

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ARNSequence_hybrid.hpp"
#include "../PAM/PAM.hpp"
#include "../Floyd/FoydPar.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    if(pid == 0) {
        cout << "=====================================" << endl;
        cout << "  PIPELINE HYBRIDE MPI + OpenMP" << endl;
        cout << "=====================================" << endl;
        cout << "MPI processes: " << nprocs << endl;
        cout << "MPI thread support: ";
        switch(provided) {
            case MPI_THREAD_SINGLE: cout << "SINGLE"; break;
            case MPI_THREAD_FUNNELED: cout << "FUNNELED"; break;
            case MPI_THREAD_SERIALIZED: cout << "SERIALIZED"; break;
            case MPI_THREAD_MULTIPLE: cout << "MULTIPLE"; break;
        }
        cout << endl;
#ifdef _OPENMP
        cout << "OpenMP threads/process: " << omp_get_max_threads() << endl;
#endif
        cout << "=====================================" << endl;
    }

    if(argc < 4) {
        if(pid == 0) {
            cerr << "Usage: " << argv[0] 
                 << " <fasta> <epsilon> <k> [output.dot] [omp_threads]" << endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    string fastaFile = argv[1];
    int epsilon = stoi(argv[2]);
    int k_clusters = stoi(argv[3]);
    string outputFile = (argc >= 5) ? argv[4] : "arn_graph.dot";
    int num_omp_threads = (argc >= 6) ? stoi(argv[5]) : 0;
    
    if(num_omp_threads > 0) {
#ifdef _OPENMP
        omp_set_num_threads(num_omp_threads);
#endif
    }

    double t_total_start = MPI_Wtime();

    // ===== ÉTAPE 1 : Lecture FASTA =====
    vector<ARNSeq> sequences;
    int nbSeq = 0;
    
    if(pid == 0) {
        cout << "\n[1/5] Lecture FASTA..." << endl;
        nbSeq = readFASTAFile(fastaFile, sequences);
        if(nbSeq <= 0) {
            cerr << "Erreur lecture FASTA" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cout << "      Séquences: " << nbSeq << endl;
    }
    
    MPI_Bcast(&nbSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Broadcast séquences
    if(pid != 0) sequences.resize(nbSeq);
    
    for(int i = 0; i < nbSeq; ++i) {
        int len = (pid == 0) ? sequences[i].sequence.length() : 0;
        int label_len = (pid == 0) ? sequences[i].label.length() : 0;
        
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&label_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if(pid != 0) {
            sequences[i].sequence.resize(len);
            sequences[i].label.resize(label_len);
            sequences[i].id = i;
        }
        
        MPI_Bcast(const_cast<char*>(sequences[i].sequence.data()), len, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(const_cast<char*>(sequences[i].label.data()), label_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    
    // ===== ÉTAPE 2 : Calcul distances HYBRIDE =====
    if(pid == 0) cout << "\n[2/5] Calcul distances (HYBRIDE)..." << endl;
    
    double t_dist_start = MPI_Wtime();
    int* distanceMatrix = computeDistanceMatrix_Hybrid(sequences, levenshteinDistance, pid, nprocs);
    double t_dist_end = MPI_Wtime();
    
    if(pid == 0) {
        cout << "      Temps: " << (t_dist_end - t_dist_start) << " sec" << endl;
    }
    
    // ===== ÉTAPE 3 : Vérifications Floyd =====
    int p_sqrt = static_cast<int>(sqrt(nprocs));
    
    if(p_sqrt * p_sqrt != nprocs) {
        if(pid == 0) {
            cerr << "ERREUR: nprocs doit être carré parfait (1,4,9,16...)" << endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    int block_size = 0;
    if(pid == 0) {
        if(nbSeq % p_sqrt != 0) {
            cerr << "ERREUR: nbSeq (" << nbSeq << ") divisible par sqrt(nprocs) (" 
                 << p_sqrt << ")" << endl;
            MPI_Finalize();
            return EXIT_FAILURE;
        }
        block_size = nbSeq / p_sqrt;
    }
    
    // CORRECTION CRITIQUE
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int* D_local = new int[block_size * block_size];
    
    // ===== ÉTAPE 4 : Floyd-Warshall =====
    if(pid == 0) cout << "\n[3/5] Floyd-Warshall..." << endl;
    
    double t_floyd_start = MPI_Wtime();
    decouperMatrice(distanceMatrix, D_local, nbSeq, block_size, p_sqrt, 0, pid);
    int* D_global = floydBlocsHybrid(D_local, nbSeq, p_sqrt, pid, 0);
    double t_floyd_end = MPI_Wtime();
    
    if(pid == 0) {
        cout << "      Temps: " << (t_floyd_end - t_floyd_start) << " sec" << endl;
    }
    
    // ===== ÉTAPE 5 : PAM =====
    if(pid == 0) {
        cout << "\n[4/5] Clustering PAM..." << endl;
        
        double t_pam_start = MPI_Wtime();
        
        pam::Result res = pam::pam_sequential(nbSeq, 
                                               vector<int>(D_global, D_global + nbSeq*nbSeq), 
                                               k_clusters, 42);
        
        double t_pam_end = MPI_Wtime();
        
        cout << "      Temps: " << (t_pam_end - t_pam_start) << " sec" << endl;
        cout << "      Coût: " << res.cost << endl;
        cout << "      Médoïdes: ";
        for(auto m : res.medoids) cout << m << " ";
        cout << endl;
        
        vector<int> cluster_counts(k_clusters, 0);
        for(int i = 0; i < nbSeq; ++i) {
            if(res.membership[i] >= 0 && res.membership[i] < k_clusters) {
                cluster_counts[res.membership[i]]++;
            }
        }
        cout << "      Points/cluster: ";
        for(auto c : cluster_counts) cout << c << " ";
        cout << endl;
        
        // ===== ÉTAPE 6 : Graphe DOT =====
        cout << "\n[5/5] Génération graphe..." << endl;
        if(writeGraphDOT(sequences, D_global, epsilon, outputFile) == 0) {
            cout << "      ✅ Graphe: " << outputFile << endl;
        }
        
        double t_total_end = MPI_Wtime();
        
        cout << "\n=====================================" << endl;
        cout << "  TEMPS TOTAL: " << (t_total_end - t_total_start) << " sec" << endl;
        cout << "=====================================" << endl;
        
        delete[] distanceMatrix;
        delete[] D_global;
    }
    
    delete[] D_local;
    MPI_Finalize();
    return 0;
}
