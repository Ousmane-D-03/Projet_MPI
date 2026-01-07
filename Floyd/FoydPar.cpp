#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "ForGraph.hpp"
#include "FoydPar.hpp"
#include "Utils.hpp"
using namespace std;

inline int& A(int* B, int b, int i, int j) { return B[i*b + j]; }

void decouperMatrice(int* D, int* D_local, int n, int block_size, int p_sqrt, int root, int pid) {
    if(pid==root){
        int* temp = new int[block_size*block_size];
        for(int bi=0; bi<p_sqrt; bi++)
            for(int bj=0; bj<p_sqrt; bj++){
                int dest = bi*p_sqrt + bj;
                for(int i=0;i<block_size;i++)
                    for(int j=0;j<block_size;j++)
                        temp[i*block_size+j] = D[(bi*block_size+i)*n + (bj*block_size+j)];
                if(dest==root)
                    copy(temp,temp+block_size*block_size,D_local);
                else
                    MPI_Send(temp, block_size*block_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
            }
        delete[] temp;
    } else {
        MPI_Recv(D_local, block_size*block_size, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int* rassemblerMatrice(int* D_local, int n, int block_size, int p_sqrt, int root, int pid)
{
    int bloc_elem = block_size * block_size;
    int nb_procs = p_sqrt * p_sqrt;

    int* gathered = nullptr;
    if (pid == root)
        gathered = new int[bloc_elem * nb_procs];

    MPI_Gather(D_local, bloc_elem, MPI_INT,
               gathered, bloc_elem, MPI_INT,
               root, MPI_COMM_WORLD);

    if (pid == root) {
        int* D = new int[n * n];
        for (int p = 0; p < nb_procs; p++) {
            int bi = p / p_sqrt;
            int bj = p % p_sqrt;
            int* src = gathered + p * bloc_elem;
            for (int i = 0; i < block_size; i++)
                for (int j = 0; j < block_size; j++)
                    D[(bi * block_size + i) * n + (bj * block_size + j)] =
                        src[i * block_size + j];
        }
        delete[] gathered;
        return D;
    }

    return nullptr;
}

void afficherBloc(int* D_local, int block_size, int pid, int nprocs, const string &titre){
    MPI_Barrier(MPI_COMM_WORLD);
    for(int p=0;p<nprocs;p++){
        if(p==pid){
            cout<<"==== "<<titre<<" (PID "<<pid<<") ====\n";
            for(int i=0;i<block_size;i++){
                for(int j=0;j<block_size;j++)
                    cout<<setw(3) <<(D_local[i*block_size+j]==INF?" ∞ ":to_string(D_local[i*block_size+j])+" ");
                cout<<"\n";
            }
            cout<<endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/**
 * @brief Floyd-Warshall par blocs avec MPI et OpenMP
 */
int* floydBlocsHybrid(int* D_local, int nb_nodes, int p_sqrt, int pid, int root, int num_threads){
    int block_size = nb_nodes/p_sqrt;
    int px = pid/p_sqrt;   // Position ligne du processus dans la grille
    int py = pid%p_sqrt;   // Position colonne du processus dans la grille
    omp_set_num_threads(num_threads);
    int* pivot = new int[block_size*block_size];
    int* row_block = new int[block_size*block_size];
    int* col_block = new int[block_size*block_size];

    // Pour chaque bloc diagonal (pivot)
    for(int k=0; k<p_sqrt; k++){
        int pivot_rank = k*p_sqrt + k;  // Processus qui possède le bloc diagonal

        // ======== PHASE 1 : Calcul du bloc pivot [k,k] ========
        if(pid == pivot_rank){
            // Floyd-Warshall standard sur le bloc diagonal
            #pragma omp parallel for collapse(2) 
            for(int kk=0; kk<block_size; kk++){
                for(int i=0; i<block_size; i++){
                    for(int j=0; j<block_size; j++){
                        if(A(D_local,block_size,i,kk) < INF && 
                           A(D_local,block_size,kk,j) < INF){
                            int new_dist = A(D_local,block_size,i,kk) + 
                                         A(D_local,block_size,kk,j);
                            if(new_dist < A(D_local,block_size,i,j)){
                                A(D_local,block_size,i,j) = new_dist;
                            }
                        }
                    }
                }
            }
            copy(D_local, D_local+block_size*block_size, pivot);
        }

        // Broadcast du bloc pivot à tous les processus
        MPI_Bcast(pivot, block_size*block_size, MPI_INT, pivot_rank, MPI_COMM_WORLD);

        // ======== PHASE 2 : Mise à jour blocs LIGNE k ========
        if(px == k && py != k){
            // Je suis dans la ligne k mais pas sur la diagonale
            #pragma omp parallel for collapse(2) 
            for(int kk=0; kk<block_size; kk++){
                for(int i=0; i<block_size; i++){
                    for(int j=0; j<block_size; j++){
                        if(pivot[i*block_size+kk] < INF && 
                           A(D_local,block_size,kk,j) < INF){
                            int new_dist = pivot[i*block_size+kk] + 
                                         A(D_local,block_size,kk,j);
                            if(new_dist < A(D_local,block_size,i,j)){
                                A(D_local,block_size,i,j) = new_dist;
                            }
                        }
                    }
                }
            }
        }
        
        // ======== PHASE 3 : Mise à jour blocs COLONNE k ========
        if(py == k && px != k){
            // Je suis dans la colonne k mais pas sur la diagonale
            #pragma omp parallel for collapse(2) 
            for(int kk=0; kk<block_size; kk++){
                for(int i=0; i<block_size; i++){
                    for(int j=0; j<block_size; j++){
                        if(A(D_local,block_size,i,kk) < INF && 
                           pivot[kk*block_size+j] < INF){
                            int new_dist = A(D_local,block_size,i,kk) + 
                                         pivot[kk*block_size+j];
                            if(new_dist < A(D_local,block_size,i,j)){
                                A(D_local,block_size,i,j) = new_dist;
                            }
                        }
                    }
                }
            }
        }

        // Synchronisation avant de broadcaster ligne et colonne
        MPI_Barrier(MPI_COMM_WORLD);

        // ======== PHASE 4 : Broadcast ligne k et colonne k ========
        // Chaque processus a besoin de deux blocs :
        // - Le bloc [k, py] (dans la ligne k, colonne py)
        // - Le bloc [px, k] (dans la ligne px, colonne k)
        
        // Broadcast de tous les blocs de la ligne k
        for(int col=0; col<p_sqrt; col++){
            int source = k*p_sqrt + col;  // Processus [k, col]
            
            int* temp_buf = new int[block_size*block_size];
            if(pid == source){
                copy(D_local, D_local+block_size*block_size, temp_buf);
            }
            MPI_Bcast(temp_buf, block_size*block_size, MPI_INT, source, MPI_COMM_WORLD);
            
            // Si c'est le bloc dont j'ai besoin, je le garde
            if(col == py){
                copy(temp_buf, temp_buf+block_size*block_size, row_block);
            }
            delete[] temp_buf;
        }
        
        // Broadcast de tous les blocs de la colonne k
        for(int row=0; row<p_sqrt; row++){
            int source = row*p_sqrt + k;  // Processus [row, k]
            
            int* temp_buf = new int[block_size*block_size];
            if(pid == source){
                copy(D_local, D_local+block_size*block_size, temp_buf);
            }
            MPI_Bcast(temp_buf, block_size*block_size, MPI_INT, source, MPI_COMM_WORLD);
            
            // Si c'est le bloc dont j'ai besoin, je le garde
            if(row == px){
                copy(temp_buf, temp_buf+block_size*block_size, col_block);
            }
            delete[] temp_buf;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // ======== PHASE 5 : Mise à jour AUTRES blocs ========
        if(px != k && py != k){
            // Maintenant row_block = bloc[k, py] et col_block = bloc[px, k]
            #pragma omp parallel for collapse(2) 
            for(int kk=0; kk<block_size; kk++){
                for(int i=0; i<block_size; i++){
                    for(int j=0; j<block_size; j++){
                        if(col_block[i*block_size+kk] < INF && 
                           row_block[kk*block_size+j] < INF){
                            int new_dist = col_block[i*block_size+kk] + 
                                         row_block[kk*block_size+j];
                            if(new_dist < A(D_local,block_size,i,j)){
                                A(D_local,block_size,i,j) = new_dist;
                            }
                        }
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    delete[] pivot;
    delete[] row_block;
    delete[] col_block;

    return rassemblerMatrice(D_local, nb_nodes, block_size, p_sqrt, root, pid);
}
