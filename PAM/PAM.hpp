// PAM.hpp
// Implementation minimale de PAM (k-medoids) séquentiel + parallèle (MPI)
// Doxygen-style comments
#ifndef PAM_HPP
#define PAM_HPP

#include <vector>

namespace pam {

struct Result {
    /** Indices des médoines (indices de sommets choisis comme centres) */
    std::vector<int> medoids;
    /** Pour chaque point, indice du médoine (valeur entre 0 et k-1) */
    std::vector<int> membership;
    /** Coût total (somme des distances au médoine le plus proche) */
    long long cost;
};


/**
 * @brief Exécute PAM en mode distribué (chaque processus peut recevoir
 * la matrice D partiellement). Cette API est générale : si MPI n'est pas
 * activé la fonction peut appeler la version séquentielle.
 *
 * @param n Nombre de points
 * @param D Matrice des distances (taille n*n), stockage row-major
 * @param k Nombre de médoines
 * @param seed Graine aléatoire
 * @param rank Rang MPI (0 si non MPI)
 * @param size Nombre de processus MPI (1 si non MPI)
 * @return Result Résultat (sur le rang 0 si distribué)
 */
Result pam_distributed(int n, const std::vector<int>& D, int k, int seed, int rank, int size);

/**
 * @brief Version séquentielle de PAM (mono-processus).
 *
 * @param n Nombre de points
 * @param D Matrice des distances (taille n*n), stockage row-major
 * @param k Nombre de médoines
 * @param seed Graine aléatoire
 * @return Result Résultat complet (médoines, affectation, coût)
 */
Result pam_sequential(int n, const std::vector<int>& D, int k, int seed);

}

#endif
