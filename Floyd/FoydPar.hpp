#ifndef FOYDPAR_HPP
#define FOYDPAR_HPP

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "ForGraph.hpp"

using namespace std;

inline int& A(int* B, int b, int i, int j);

/**
 * @brief Découpe la matrice globale en blocs et les distribue
 * 
 * @param D Matrice globale (sur processus root)
 * @param D_local Bloc local (sortie)
 * @param n Taille de la matrice globale
 * @param block_size Taille d'un bloc
 * @param p_sqrt Racine carrée du nombre de processus
 * @param root Processus racine
 * @param pid Identifiant du processus courant
 */
void decouperMatrice(int* D, int* D_local,
                     int n, int block_size, int p_sqrt,
                     int root, int pid);

/**
 * @brief Rassemble les blocs distribués en matrice globale
 * 
 * @param D_local Bloc local
 * @param n Taille de la matrice globale
 * @param block_size Taille d'un bloc
 * @param p_sqrt Racine carrée du nombre de processus
 * @param root Processus racine
 * @param pid Identifiant du processus courant
 * @return int* Matrice globale (NULL si pid != root)
 */
int* rassemblerMatrice(int* D_local,
                       int n, int block_size, int p_sqrt,
                       int root, int pid);

void afficherBloc(int* D_local, int block_size,
                  int pid, int nprocs, const string &titre);

/**
 * @brief Algorithme de Floyd-Warshall par blocs (VERSION HYBRIDE MPI+OpenMP)
 * 
 * Utilise MPI pour la distribution des blocs et OpenMP pour paralléliser
 * les calculs au sein de chaque bloc.
 * 
 * @param D_local Bloc local de la matrice
 * @param nb_nodes Nombre de nœuds du graphe
 * @param p_sqrt Racine carrée du nombre de processus
 * @param pid Identifiant du processus courant
 * @param root Processus racine
 * @param num_threads Nombre de threads OpenMP (0 = automatique)
 * @return int* Matrice globale (NULL si pid != root)
 * 
 * @note Chaque processus MPI utilise num_threads threads OpenMP
 * @note Configuration optimale : p × t ≈ nombre de cœurs physiques
 */
int* floydBlocsHybrid(int* D_local,
                   int nb_nodes, int p_sqrt, int pid, int root, int num_threads);

#endif
