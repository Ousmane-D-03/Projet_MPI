#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "ForGraph.hpp"
#include "FoydPar.hpp"
#include "Utils.hpp"
#include <omp.h>

using namespace std;

inline int& A(int* B, int b, int i, int j) { return B[i*b + j]; }

/**
 * @brief Découpe la matrice globale en blocs et les distribue aux processus
 * @param D matrice globale
 * @param D_local bloc local
 * @param n taille de la matrice globale
 * @param block_size taille d’un bloc
 * @param p_sqrt racine carrée du nombre de processus
 * @param root processus racine
 * @param pid identifiant du processus courant
 */
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

/**
 * @brief Rassemble les blocs distribués en matrice globale
 * @param D_local bloc local
 * @param n taille de la matrice globale
 * @param block_size taille d’un bloc
 * @param p_sqrt racine carrée du nombre de processus
 * @param root processus racine
 * @param pid identifiant du processus courant
 * @return pointeur vers la matrice globale (NULL si non root)
 */
int* rassemblerMatrice(int* D_local,
                       int n, int block_size, int p_sqrt, int root, int pid)
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

/**
 * @brief Affiche un bloc local de matrice
 * @param D_local bloc local
 * @param block_size taille du bloc
 * @param pid id du processus
 * @param nprocs nombre total de processus
 * @param titre titre à afficher
 */
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
 * @brief Algorithme de Floyd–Warshall par blocs (MPI)
 * @param D_local bloc local de matrice
 * @param nb_nodes nombre de nœuds du graphe
 * @param p_sqrt racine carrée du nombre de processus
 * @param pid id du processus courant
 * @param root id du processus racine
 * @return pointeur vers la matrice globale (uniquement sur le root, NULL sinon)
 */
int* floydBlocsHybrid(int* D_local, int nb_nodes, int p_sqrt, int pid, int root){
    omp_set_num_threads(omp_get_max_threads());
    int block_size = nb_nodes/p_sqrt;
    int px = pid/p_sqrt;
    int py = pid%p_sqrt;

    int* pivot = new int[block_size*block_size];
    int* row_block = new int[block_size*block_size];
    int* col_block = new int[block_size*block_size];

    for(int k=0;k<p_sqrt;k++){
        int pivot_rank = k*p_sqrt + k;

        if(pid==pivot_rank){
            #pragma omp parallel for collapse(2) schedule(static)
            for(int i=0;i<block_size;i++)
                for(int j=0;j<block_size;j++)
                    for(int x=0;x<block_size;x++)
                        if(A(D_local,block_size,i,x)<INF && A(D_local,block_size,x,j)<INF)
                            A(D_local,block_size,i,j) = min(A(D_local,block_size,i,j),
                                                            A(D_local,block_size,i,x)+A(D_local,block_size,x,j));
            copy(D_local,D_local+block_size*block_size,pivot);
        }

        MPI_Comm row_comm, col_comm;
        MPI_Comm_split(MPI_COMM_WORLD, px, pid, &row_comm);
        MPI_Comm_split(MPI_COMM_WORLD, py, pid, &col_comm);

        MPI_Bcast(pivot, block_size*block_size, MPI_INT, k, row_comm);
        MPI_Bcast(pivot, block_size*block_size, MPI_INT, k, col_comm);

        if(px==k && py!=k){
            #pragma omp parallel for collapse(2) schedule(static)
            for(int i=0;i<block_size;i++)
                for(int j=0;j<block_size;j++)
                    for(int x=0;x<block_size;x++)
                        if(pivot[i*block_size+x]<INF && A(D_local,block_size,x,j)<INF)
                            A(D_local,block_size,i,j) = min(A(D_local,block_size,i,j),
                                                            pivot[i*block_size+x]+A(D_local,block_size,x,j));
        }
        if(py==k && px!=k){
            #pragma omp parallel for collapse(2) schedule(static)
            for(int i=0;i<block_size;i++)
                for(int j=0;j<block_size;j++)
                    for(int x=0;x<block_size;x++)
                        if(A(D_local,block_size,i,x)<INF && pivot[x*block_size+j]<INF)
                            A(D_local,block_size,i,j) = min(A(D_local,block_size,i,j),
                                                            A(D_local,block_size,i,x)+pivot[x*block_size+j]);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        int row_root=-1, col_root=-1;
        for(int p=0;p<p_sqrt*p_sqrt;p++){
            if(p/p_sqrt==k) row_root=p;
            if(p%p_sqrt==k) col_root=p;
        }
        if(px==k) copy(D_local,D_local+block_size*block_size,row_block);
        if(py==k) copy(D_local,D_local+block_size*block_size,col_block);

        MPI_Bcast(row_block, block_size*block_size, MPI_INT, row_root, MPI_COMM_WORLD);
        MPI_Bcast(col_block, block_size*block_size, MPI_INT, col_root, MPI_COMM_WORLD);

        if(px!=k && py!=k){
            #pragma omp parallel for collapse(2) schedule(static)
            for(int i=0;i<block_size;i++)
                for(int j=0;j<block_size;j++)
                    for(int x=0;x<block_size;x++)
                        if(col_block[i*block_size+x]<INF && row_block[x*block_size+j]<INF)
                            A(D_local,block_size,i,j) = min(A(D_local,block_size,i,j),
                                                            col_block[i*block_size+x]+row_block[x*block_size+j]);
        }

        MPI_Comm_free(&row_comm);
        MPI_Comm_free(&col_comm);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    delete[] pivot;
    delete[] row_block;
    delete[] col_block;

    // Rassemble les blocs sur le root pour obtenir la matrice globale
    return rassemblerMatrice(D_local, nb_nodes, block_size, p_sqrt, root, pid);
}


