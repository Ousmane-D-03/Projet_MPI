// main_pam.cpp
// Simple driver to run PAM. It can construct distance matrix from a .dot via existing Floyd code (requires Graphviz)

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <unistd.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "PAM.hpp"

#ifdef WITH_GRAPHVIZ
#include "../Floyd/ForGraph.hpp"
#endif

using namespace std;

static void usage(const char* prog) {
    if (prog) cerr << "Usage: " << prog << " <graph.dot> <k> [seed]\n";
}

int main(int argc, char* argv[]) {
    int rank = 0, size = 1;
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
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

    int n = 0;
    vector<int> D;

    // Two input modes:
    // - WITH_GRAPHVIZ defined at compile time: read .dot file and use Floyd code (Graphviz required)
    // - otherwise: input is a distance matrix text file: first line n, then n lines with n ints
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
                // read distance matrix from text file `dotfile`
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
                    cerr << "Failed to read n from distance file" << endl;
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

    // broadcast n then data
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
#endif

   
    

    // run PAM distributed
    pam::Result r = pam::pam_distributed(n, D, k, seed, rank, size);

    if (rank == 0) {
        cout << "Cost: " << r.cost << "\nMedoids:";
        for (int m: r.medoids) cout << " " << m;
        cout << "\n";
        // Print membership counts
        vector<int> counts(k,0);
        for (int i = 0; i < (int)r.membership.size(); ++i) if (r.membership[i] >= 0 && r.membership[i] < k) counts[r.membership[i]]++;
        cout << "Counts per medoid:";
        for (int c: counts) cout << " " << c;
        cout << "\n";

        // Detailed membership
        cout << "Membership:\n";
        for (int i = 0; i < (int)r.membership.size(); ++i) {
            cout << "Point " << i << " -> Medoid " << r.membership[i] << "\n";

        }
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
