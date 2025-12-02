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

// Accès à un élément dans un bloc
inline int& A(int* B, int b, int i, int j);

// Découpe la matrice globale en blocs puis distribution MPI
void decouperMatrice(int* D, int* D_local,
                     int n, int block_size, int p_sqrt,
                     int root, int pid);

// Rassemble les blocs distribués en une matrice globale
int* rassemblerMatrice(int* D_local,
                       int n, int block_size, int p_sqrt,
                       int root, int pid);

// Affiche un bloc local pour debug
void afficherBloc(int* D_local, int block_size,
                  int pid, int nprocs, const string &titre);

// Algorithme de Floyd–Warshall par blocs en MPI
int* floydBlocsMPI(int* D_local,
                   int nb_nodes, int p_sqrt, int pid, int root);

#endif
