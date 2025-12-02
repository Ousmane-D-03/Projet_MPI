#include <mpi.h>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <iomanip>
#include "ForGraph.hpp"
#include "FoydPar.hpp"
#include "Utils.hpp"
#include <vector>
using namespace std;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc != 2) {
        if (pid == 0)
            cout << "Usage : ./main fichier.dot (graphe au format dot)" << endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char* file_name = argv[1];
    int nb_nodes = 0;
    map<string,int> my_nodes;
    int* D = nullptr;

    // Lecture du graphe sur le root
    if (pid == 0) {
        int* mat_adjacence = lectureGraphe(file_name, &nb_nodes, &my_nodes);
        if (!mat_adjacence) {
            cerr << "Erreur : impossible de lire le graphe !" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cout << "=== Matrice d'adjacence ===" << endl;
        affichage(mat_adjacence, nb_nodes, nb_nodes, 2);
        cout << endl;

        D = InitDk(nb_nodes, mat_adjacence);
        delete[] mat_adjacence;

        cout << "=== Matrice de distances initiale ===" << endl;
        affichage(D, nb_nodes, nb_nodes, 3);
        cout << endl;
    }

    // Diffusion du nombre de noeuds à tous les processus
    MPI_Bcast(&nb_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Vérification du nombre de processus carré parfait
    int p_sqrt = (int)sqrt(nprocs);
    if (p_sqrt * p_sqrt != nprocs || nb_nodes % p_sqrt != 0) {
        if (pid == 0)
            cerr << "Erreur : nprocs pas carré parfait ou nb_nodes non divisible" << endl;
        MPI_Finalize();
        return -1;
    }

    int block_size = nb_nodes / p_sqrt;
    int* D_local = new int[block_size * block_size];

    // ----- CALCUL SEQUENTIEL -----
    double t_seq_start = 0, t_seq_end = 0;
    if (pid == 0) {
        t_seq_start = MPI_Wtime();
        int* D_seq = MatDistance(nb_nodes, D);
        t_seq_end = MPI_Wtime();

        cout << "=== Matrice de distances (séquentiel) ===" << endl;
        affichage(D_seq, nb_nodes, nb_nodes, 3);
        delete[] D_seq;

        cout << "\nTemps séquentiel : " << (t_seq_end - t_seq_start) << " sec" << endl;
        cout << "\n" << endl;
    }

    // ----- CALCUL PARALLELE -----
    decouperMatrice(D, D_local, nb_nodes, block_size, p_sqrt, 0, pid);

    double t0 = MPI_Wtime();
    int* D_final = floydBlocsMPI(D_local, nb_nodes, p_sqrt, pid, 0);
    double t1 = MPI_Wtime();

    double local_time = t1 - t0, max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (pid == 0) {
        cout << "=== Matrice globale après Floyd par blocs MPI ===" << endl;
        affichage(D_final, nb_nodes, nb_nodes, 3);

        cout << "\nTemps parallèle : " << max_time << " sec" << endl;

        delete[] D_final;
        delete[] D;
    }

    delete[] D_local;
    MPI_Finalize();
    return 0;
}
