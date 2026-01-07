#include <mpi.h>
#include <omp.h>
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
    int pid, nprocs;
    int provided;
    
    // Initialisation MPI avec support multi-threading
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc < 2 || argc > 3) {
        if (pid == 0) {
            cout << "Usage : mpirun -np <P> ./main fichier.dot [num_threads]" << endl;
            cout << endl;
            cout << "Paramètres:" << endl;
            cout << "  <P>           : Nombre de processus MPI (carré parfait: 4, 9, 16...)" << endl;
            cout << "  fichier.dot   : Graphe au format DOT" << endl;
            cout << "  [num_threads] : Threads OpenMP par processus (défaut: auto)" << endl;
            cout << endl;
            cout << "Exemples:" << endl;
            cout << "  mpirun -np 4 ./main Exemple2.dot" << endl;
            cout << "  mpirun -np 4 ./main Exemple2.dot 2" << endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char* file_name = argv[1];
    int num_threads = 4; // Default value
    
    // Configuration OpenMP
    if(num_threads > 0){
        omp_set_num_threads(num_threads);
    }

    int nb_nodes = 0;
    map<string,int> my_nodes;
    int* D = nullptr;

    if (pid == 0) {
        cout << "╔═══════════════════════════════════════════════════════════╗" << endl;
        cout << "║     FLOYD-WARSHALL HYBRIDE MPI+OPENMP                     ║" << endl;
        cout << "╚═══════════════════════════════════════════════════════════╝" << endl;
        cout << endl;
        cout << "Configuration:" << endl;
        cout << "  Fichier       : " << file_name << endl;
        cout << "  Processus MPI : " << nprocs << endl;
        cout << "  Threads/proc  : " << (num_threads > 0 ? to_string(num_threads) : "auto") 
             << " (max: " << omp_get_max_threads() << ")" << endl;
        cout << "  Total workers : " << nprocs * omp_get_max_threads() << endl;
        cout << "  Niveau thread : " << provided;
        if(provided >= MPI_THREAD_FUNNELED) cout << " ✓";
        cout << endl << endl;
    }

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
        if (pid == 0) {
            cerr << "Erreur : nprocs=" << nprocs << " pas carré parfait" << endl;
            cerr << "   ou nb_nodes=" << nb_nodes << " non divisible par √P=" << p_sqrt << endl;
        }
        MPI_Finalize();
        return -1;
    }

    int block_size = nb_nodes / p_sqrt;
    int* D_local = new int[block_size * block_size];

    // ----- CALCUL SEQUENTIEL -----
    double t_seq_start = 0, t_seq_end = 0;
    if (pid == 0) {
        cout << "┌─────────────────────────────────────────────────────────┐" << endl;
        cout << "│ CALCUL SÉQUENTIEL (référence)                          │" << endl;
        cout << "└─────────────────────────────────────────────────────────┘" << endl;
        
        t_seq_start = MPI_Wtime();
        int* D_seq = MatDistance(nb_nodes, D);
        t_seq_end = MPI_Wtime();

        cout << "=== Matrice de distances (séquentiel) ===" << endl;
        affichage(D_seq, nb_nodes, nb_nodes, 3);
        delete[] D_seq;

        cout << "\n✓ Temps séquentiel : " << (t_seq_end - t_seq_start) << " sec" << endl;
        cout << endl;
    }

    // ----- CALCUL PARALLELE HYBRIDE -----
    if (pid == 0) {
        cout << "┌─────────────────────────────────────────────────────────┐" << endl;
        cout << "│ CALCUL HYBRIDE MPI+OPENMP                               │" << endl;
        cout << "└─────────────────────────────────────────────────────────┘" << endl;
        cout << "  Grille        : " << p_sqrt << " × " << p_sqrt << " processus" << endl;
        cout << "  Taille bloc   : " << block_size << " × " << block_size << endl;
        cout << "  Mode          : Hybride (MPI + OpenMP)" << endl;
        cout << endl;
    }
    
    decouperMatrice(D, D_local, nb_nodes, block_size, p_sqrt, 0, pid);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    
    int* D_final = floydBlocsHybrid(D_local, nb_nodes, p_sqrt, pid, 0, num_threads);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    double local_time = t1 - t0, max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (pid == 0) {
        cout << "=== Matrice globale après Floyd par blocs Hybride (MPI+OpenMP) ===" << endl;
        affichage(D_final, nb_nodes, nb_nodes, 3);

        cout << "\n✓ Temps parallèle : " << max_time << " sec" << endl;
        cout << endl;
        
        // Statistiques
        if(t_seq_end > t_seq_start && max_time > 0){
            double speedup = (t_seq_end - t_seq_start) / max_time;
            double efficiency = speedup / (nprocs * omp_get_max_threads()) * 100;
            
            cout << "╔═══════════════════════════════════════════════════════════╗" << endl;
            cout << "║                  RÉSUMÉ DES PERFORMANCES                  ║" << endl;
            cout << "╚═══════════════════════════════════════════════════════════╝" << endl;
            cout << "  Taille graphe      : " << nb_nodes << " noeuds" << endl;
            cout << "  Temps séquentiel   : " << (t_seq_end - t_seq_start) << " sec" << endl;
            cout << "  Temps parallèle    : " << max_time << " sec" << endl;
            cout << "  Speedup            : " << speedup << "x" << endl;
            cout << "  Efficacité         : " << efficiency << "%" << endl;
            cout << "  Configuration      : " << nprocs << " proc × " 
                 << omp_get_max_threads() << " threads = " 
                 << (nprocs * omp_get_max_threads()) << " workers" << endl;
            cout << endl;
        }

        delete[] D_final;
        delete[] D;
    }

    delete[] D_local;
    MPI_Finalize();
    return 0;
}
