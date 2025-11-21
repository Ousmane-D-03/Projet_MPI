// PAM.hpp
// Implementation minimale de PAM (k-medoids) séquentiel + parallèle (MPI)
// Doxygen-style comments
#ifndef PAM_HPP
#define PAM_HPP

#include <vector>

namespace pam {

struct Result {
    std::vector<int> medoids; // indices des medoids
    std::vector<int> membership; // pour chaque point, index du medoid
    long long cost;
};


// - n: nombre de points
// - D: vecteur length n*n (row-major)
// - k: nombre de medoids
// - seed: seed aléatoire
// - rank,size: fournis par MPI
// Retourne Result (valide sur tous les processus)
Result pam_distributed(int n, const std::vector<int>& D, int k, int seed, int rank, int size);

// Version séquentielle pratique (wrapper)
Result pam_sequential(int n, const std::vector<int>& D, int k, int seed);

}

#endif
