// PAM.cpp
#include "PAM.hpp"

#include <algorithm>
#include <random>
#include <limits>
#include <set>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif


using namespace std;

namespace pam {

/**
 * @brief calculer_affectation
 * Calcule pour chaque point de l'intervalle [start,end) :
 * - le médoïde le plus proche (index dans medoids)
 * - la distance au meilleur et au second meilleur médoïde
 *
 * @param n nombre de points total
 * @param D matrice des distances (row-major)
 * @param medoids indices des medoids
 * @param start indice de départ (inclus)
 * @param end indice de fin (exclus)
 * @param membership sortie : index du medoid (taille end-start)
 * @param bestDist sortie : distance au medoid le plus proche
 * @param secondBestDist sortie : distance au 2e meilleur
 */
static void calculer_affectation(int n, const vector<int>& D, const vector<int>& medoids,
                                 int start, int end,
                                 vector<int>& membership, vector<int>& bestDist, vector<int>& secondBestDist) {
    int k = (int)medoids.size();
    for (int ii = start; ii < end; ++ii) {
        int best = numeric_limits<int>::max();
        int second = numeric_limits<int>::max();
        int bestMed = -1;
        for (int m = 0; m < k; ++m) {
            int dist = D[ii*n + medoids[m]];
            if (dist < best) {
                second = best;
                best = dist;
                bestMed = m;
            } else if (dist < second) {
                second = dist;
            }
        }
        membership[ii - start] = bestMed;
        bestDist[ii - start] = best;
        secondBestDist[ii - start] = second;
    }
}

/**
 * @brief pam_sequential
 * Version simple et lisible de l'algorithme PAM (k-médoïdes).
 * Initialise k médoines aléatoirement puis tente des échanges
 * (swap) tant qu'ils réduisent le coût total.
 *
 * @param n nombre de points
 * @param D matrice des distances (row-major)
 * @param k nombre de médoines
 * @param seed graine aléatoire
 * @return Result résultat complet (médoines, affectation, coût)
 */
Result pam_sequential(int n, const vector<int>& D, int k, int seed) {
    Result res;
    if (k <= 0 || k > n) {
        cerr << "Invalid k" << endl;
        return res;
    }

    vector<int> medoids;
    medoids.reserve(k);
    std::mt19937 rng(seed);
    // simple random init without replacement
    vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < k; ++i) medoids.push_back(idx[i]);

    vector<int> membership(n);
    vector<int> bestDist(n), secondBestDist(n);

    // initial membership
    calculer_affectation(n, D, medoids, 0, n, membership, bestDist, secondBestDist);

    long long cost = 0;
    for (int i = 0; i < n; ++i) cost += bestDist[i];

    bool improved = true;
    while (improved) {
        improved = false;
        long long bestDelta = 0;
        int bestSwapMed = -1;
        int bestSwapCand = -1;

        // brute force all swaps
        set<int> medset(medoids.begin(), medoids.end());
        for (int mi = 0; mi < k; ++mi) {
            int med = medoids[mi];
            for (int cand = 0; cand < n; ++cand) {
                if (medset.count(cand)) continue;
                long long delta = 0;
                for (int i = 0; i < n; ++i) {
                    int distToCand = D[i*n + cand];
                    if (membership[i] == mi) {
                        // currently assigned to the med being swapped out
                        int newdist = min(distToCand, secondBestDist[i]);
                        delta += (long long)newdist - bestDist[i];
                    } else {
                        if (distToCand < bestDist[i]) delta += (long long)distToCand - bestDist[i];
                    }
                }
                if (delta < bestDelta) {
                    bestDelta = delta;
                    bestSwapMed = mi;
                    bestSwapCand = cand;
                }
            }
        }

        if (bestDelta < 0) {
            // apply swap
            medoids[bestSwapMed] = bestSwapCand;
            // recompute membership and dists
            calculer_affectation(n, D, medoids, 0, n, membership, bestDist, secondBestDist);
            cost += bestDelta;
            improved = true;
        }
    }

    res.medoids = medoids;
    res.membership = membership;
    res.cost = cost;
    return res;
}

/**
 * @brief pam_distributed
 * Version distribuée de PAM : chaque processus reçoit une partie des lignes
 * de la matrice D (scatterv) ; les calculs de delta sont réduits globalement
 * pour décider des échanges.
 *
 * @param n nombre de points
 * @param D matrice des distances complète (sur le rang 0) ou vide sur les autres
 * @param k nombre de médoines
 * @param seed graine aléatoire
 * @param rank rang MPI
 * @param size nombre de processus
 * @return Result résultat (valide sur le rang 0)
 */
#ifndef USE_MPI
Result pam_distributed(int n, const vector<int>& D, int k, int seed, int /*rank*/, int /*size*/) {
    return pam_sequential(n, D, k, seed);
}
#else
Result pam_distributed(int n, const vector<int>& D, int k, int seed, int rank, int size) {
    Result res;

    if (k <= 0 || k > n) {
        if (rank == 0) cerr << "Erreur : k invalide" << endl;
        return res;
    }

    // Distribution de la matrice D par lignes
    vector<int> localD;
    int rowsPerProc = n / size;
    int remainder = n % size;
    vector<int> sendCounts(size);
    vector<int> displs(size);
    int offset = 0;

    // remplir sendCounts et displs
    for (int i = 0; i < size; ++i) {
         sendCounts[i] = (i < remainder) ? (rowsPerProc + 1) * n : rowsPerProc * n;
        displs[i] = offset;
        offset += sendCounts[i];
    }

    int localRows = sendCounts[rank] / n;
    localD.resize(localRows * n);

    MPI_Scatterv(D.data(), sendCounts.data(), displs.data(), MPI_INT,
                 localD.data(), localRows * n, MPI_INT, 0, MPI_COMM_WORLD);

    // Initialisation aléatoire des medoids sur rank 0
    vector<int> medoids;
    if (rank == 0) {
        medoids.reserve(k);
        std::mt19937 rng(seed);
        vector<int> perm(n);
        for (int i = 0; i < n; ++i) perm[i] = i;
        shuffle(perm.begin(), perm.end(), rng);
        for (int i = 0; i < k; ++i) medoids.push_back(perm[i]);
    }

    // Broadcast des medoids à tous les processus
    medoids.resize(k);
    MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);

    // Tableaux d'affectation locaux
    vector<int> membership(localRows);
    vector<int> bestDist(localRows), secondBestDist(localRows);

    // Calcul initial des affectations locales
    calculer_affectation(n, localD, medoids, 0, localRows, membership, bestDist, secondBestDist);

    // Coût total local
    long long localCost = 0;
    for (int i = 0; i < localRows; ++i) localCost += bestDist[i];

    bool changed = true;
    while (changed) {
        changed = false;
        long long bestDelta = 0;
        int bestMedIndex = -1;
        int bestCand = -1;
        set<int> medset(medoids.begin(), medoids.end());

        for (int mi = 0; mi < k; ++mi) {
            for (int cand = 0; cand < n; ++cand) {
                if (medset.count(cand)) continue;

                long long delta = 0;
                for (int i = 0; i < localRows; ++i) {
                    int distToCand = localD[i*n + cand];
                    if (membership[i] == mi) {
                        int nd = min(distToCand, secondBestDist[i]);
                        delta += (long long)nd - bestDist[i];
                    } else {
                        if (distToCand < bestDist[i]) delta += (long long)distToCand - bestDist[i];
                    }
                }

                long long globalDelta = 0;
                MPI_Allreduce(&delta, &globalDelta, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

                if (globalDelta < bestDelta) {
                    bestDelta = globalDelta;
                    bestMedIndex = mi;
                    bestCand = cand;
                }
            }
        }

        if (bestDelta < 0) {
            // Appliquer l'échange
            medoids[bestMedIndex] = bestCand;
            MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);
            calculer_affectation(n, localD, medoids, 0, localRows, membership, bestDist, secondBestDist);
            localCost += bestDelta; // chaque proc ajoute la même valeur
            changed = true;
        }
    }

    // Rassembler le membership complet sur rank 0
    vector<int> fullMembership;
    vector<int> recvCounts(size);
    for (int i = 0; i < size; ++i) {
        recvCounts[i] = (i < remainder ? rowsPerProc + 1 : rowsPerProc);
        displs[i] = (i == 0 ? 0 : displs[i-1] + recvCounts[i-1]);
    }

    if (rank == 0) fullMembership.resize(n);
    MPI_Gatherv(membership.data(), localRows, MPI_INT,
                fullMembership.data(), recvCounts.data(), displs.data(),
                MPI_INT, 0, MPI_COMM_WORLD);

    // Calcul du coût total final
    long long totalCost = 0;
    MPI_Reduce(&localCost, &totalCost, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);  
    // Rassembler les résultats sur rank 0 
    if (rank == 0) {
        res.medoids = medoids;
        res.membership = fullMembership;
        res.cost = totalCost;
    }
    return res;
}
#endif

} // namespace pam
