// main_pam_hybrid.cpp - Version HYBRIDE
// À placer dans PAM/main_pam_hybrid.cpp (garde main_pam.cpp original)

#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "PAM.hpp"

#ifdef WITH_GRAPHVIZ
#include "../Floyd/ForGraph.hpp"
#endif

using namespace std;

static void usage(const char* prog) {
    if (prog) cerr << "Usage: " << prog << " <graph.dot> <k> [seed] [num_threads]\n";
}

int main(int argc, char* argv[]) {
    int rank = 0, size = 1;
    
#ifdef USE_MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(rank == 0) {
        cout << "=== PAM HYBRIDE MPI + OpenMP ===" << endl;
        cout << "MPI processes: " << size << endl;
        cout << "MPI Thread support: ";
        switch(provided) {
            case MPI_THREAD_SINGLE: cout << "SINGLE" << endl; break;
            case MPI_THREAD_FUNNELED: cout << "FUNNELED" << endl; break;
            case MPI_THREAD_SERIALIZED: cout << "SERIALIZED" << endl; break;
            case MPI_THREAD_MULTIPLE: cout << "MULTIPLE" << endl; break;
        }
    }
#endif

    if (argc < 3) {
        if (rank == 0) usage(argv[0]);
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    char* dotfile = argv[1];
    int k = stoi(argv[2]);
    int seed = (argc >= 4) ? stoi(argv[3]) : 12345;
    int num_threads = (argc >= 5) ? stoi(argv[4]) : 0;
    
    if(num_threads > 0) {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
        if(rank == 0) cout << "OpenMP threads: " << num_threads << endl;
#endif
    } else {
#ifdef _OPENMP
        if(rank == 0) cout << "OpenMP threads: " << omp_get_max_threads() << endl;
#endif
    }

    int n = 0;
    vector<int> D;

    if (rank == 0) {
        #ifdef WITH_GRAPHVIZ
            map<string,int> nodes;
            int* mat_adj = lectureGraphe(dotfile, &n, &nodes);
            if (!mat_adj) {
                cerr << "Failed to read graph" << endl;
        #ifdef USE_MPI
                MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                return 1;
        #endif
            }
            int* Dk = MatDistance(n, mat_adj);
            D.resize(n*n);
            for (int i = 0; i < n*n; ++i) D[i] = Dk[i];
            delete[] mat_adj;
            delete[] Dk;
        #else
            FILE* f = fopen(dotfile, "r");
            if (!f) {
                cerr << "Failed to open distance file: " << dotfile << endl;
        #ifdef USE_MPI
                MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                return 1;
        #endif
            }
            if (fscanf(f, "%d", &n) != 1) {
                cerr << "Failed to read n" << endl;
                fclose(f);
        #ifdef USE_MPI
                MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                return 1;
        #endif
            }
            D.resize(n * n);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    int v;
                    if (fscanf(f, "%d", &v) != 1) {
                        cerr << "Failed to read distance at " << i << "," << j << endl;
                        fclose(f);
        #ifdef USE_MPI
                        MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                        return 1;
        #endif
                    }
                    D[i*n + j] = v;
                }
            }
            fclose(f);
        #endif
    }

#ifdef USE_MPI
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    
    if (n <= 0) {
        if (rank == 0) cerr << "n <= 0" << endl;
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if (rank != 0) D.resize(n*n);

#ifdef USE_MPI
    MPI_Bcast(D.data(), n*n, MPI_INT, 0, MPI_COMM_WORLD);
    
    double t_start = MPI_Wtime();
    pam::Result r = pam::pam_distributed(n, D, k, seed, rank, size);
    double t_end = MPI_Wtime();
#else
    double t_start = omp_get_wtime();
    pam::Result r = pam::pam_sequential(n, D, k, seed);
    double t_end = omp_get_wtime();
#endif

    if (rank == 0) {
        cout << "\n=== RÉSULTATS ===" << endl;
        cout << "Temps: " << (t_end - t_start) << " sec" << endl;
        cout << "Cost: " << r.cost << "\nMedoids:";
        for (int m: r.medoids) cout << " " << m;
        cout << "\n";
        
        vector<int> counts(k,0);
        for (int i = 0; i < (int)r.membership.size(); ++i) 
            if (r.membership[i] >= 0 && r.membership[i] < k) 
                counts[r.membership[i]]++;
        
        cout << "Points per cluster:";
        for (int c: counts) cout << " " << c;
        cout << "\n";
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
