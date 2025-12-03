// main_mpi.cpp
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include "ARNSequence.hpp"
#include "../PAM/PAM.hpp"   // adapter selon ton chemin
#include "../Floyd/FoydPar.hpp"       // ton fichier avec floydBlocsMPI
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

    vector<ARNSeq> sequences;
    int nbSeq = 0;
    int* distanceMatrix = nullptr;

    if(pid == 0) {
        nbSeq = readFASTAFile(fastaFile, sequences);
        if(nbSeq <= 0) {
            cerr << "Erreur : impossible de lire le fichier FASTA" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cout << "Nombre de séquences : " << nbSeq << endl;

        // calcul de la matrice de distances
        distanceMatrix = computeDistanceMatrix(sequences, levenshteinDistance);
    }

    // Vérifier que nbSeq est divisible par sqrt(nprocs)
    int p_sqrt = static_cast<int>(sqrt(nprocs));
    if(p_sqrt*p_sqrt != nprocs) {
        if(pid == 0) cerr << "Erreur : le nombre de processus doit être un carré parfait" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int block_size = (pid==0 ? nbSeq : 0)/p_sqrt; // taille d'un bloc par ligne/col
    MPI_Bcast(&nbSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* D_local = new int[block_size*block_size];

    // Découper et distribuer la matrice
    decouperMatrice(distanceMatrix, D_local, nbSeq, block_size, p_sqrt, 0, pid);

    // Appliquer Floyd–Warshall MPI
    int* D_global = floydBlocsMPI(D_local, nbSeq, p_sqrt, pid, 0);

    // Sur le root : appliquer PAM
    if(pid == 0) {
        cout << "Matrice des plus courts chemins calculée." << endl;

        // PAM séquentiel sur la matrice complète
        pam::Result res = pam::pam_sequential(nbSeq, vector<int>(D_global, D_global + nbSeq*nbSeq), k_clusters, 42);
        cout << "Résultat PAM :" << endl;
        cout << "Médoïdes : ";
        for(auto m : res.medoids) cout << m << " ";
        cout << endl;

        // Écriture du graphe DOT
        if(writeGraphDOT(sequences, D_global, epsilon, outputFile) == 0) {
            cout << "Graphe écrit dans " << outputFile << endl;
        }

        delete[] distanceMatrix;
        delete[] D_global;
    }

    delete[] D_local;
    MPI_Finalize();
    return 0;
}

