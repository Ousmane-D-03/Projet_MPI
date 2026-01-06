// PAM_hybrid.cpp - Version HYBRIDE MPI + OpenMP
// À placer dans PAM/PAM_hybrid.cpp (garde PAM.cpp original)

#include "PAM.hpp"

#include <algorithm>
#include <random>
#include <limits>
#include <set>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace pam {

/**
 * @brief calculer_affectation - Version HYBRIDE avec OpenMP
 */
static void calculer_affectation(int n, const vector<int>& D, const vector<int>& medoids,
                                 int start, int end,
                                 vector<int>& membership, vector<int>& bestDist, vector<int>& secondBestDist) {
    int k = (int)medoids.size();
    
    #pragma omp parallel for schedule(static)
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
 * @brief pam_sequential - Version avec OpenMP
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
    vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < k; ++i) medoids.push_back(idx[i]);

    vector<int> membership(n);
    vector<int> bestDist(n), secondBestDist(n);

    calculer_affectation(n, D, medoids, 0, n, membership, bestDist, secondBestDist);

    long long cost = 0;
    #pragma omp parallel for reduction(+:cost)
    for (int i = 0; i < n; ++i) cost += bestDist[i];

    bool improved = true;
    while (improved) {
        improved = false;
        long long bestDelta = 0;
        int bestSwapMed = -1;
        int bestSwapCand = -1;

        set<int> medset(medoids.begin(), medoids.end());
        
        vector<int> candidates;
        for(int c = 0; c < n; ++c) {
            if(!medset.count(c)) candidates.push_back(c);
        }
        
        #pragma omp parallel
        {
            long long thread_best_delta = 0;
            int thread_best_mi = -1;
            int thread_best_cand = -1;
            
            #pragma omp for collapse(2) schedule(dynamic, 16)
            for (int mi = 0; mi < k; ++mi) {
                for (int c_idx = 0; c_idx < (int)candidates.size(); ++c_idx) {
                    int cand = candidates[c_idx];
                    long long delta = 0;
                    
                    for (int i = 0; i < n; ++i) {
                        int distToCand = D[i*n + cand];
                        if (membership[i] == mi) {
                            int newdist = min(distToCand, secondBestDist[i]);
                            delta += (long long)newdist - bestDist[i];
                        } else {
                            if (distToCand < bestDist[i]) 
                                delta += (long long)distToCand - bestDist[i];
                        }
                    }
                    
                    if (delta < thread_best_delta) {
                        thread_best_delta = delta;
                        thread_best_mi = mi;
                        thread_best_cand = cand;
                    }
                }
            }
            
            #pragma omp critical
            {
                if (thread_best_delta < bestDelta) {
                    bestDelta = thread_best_delta;
                    bestSwapMed = thread_best_mi;
                    bestSwapCand = thread_best_cand;
                }
            }
        }

        if (bestDelta < 0) {
            medoids[bestSwapMed] = bestSwapCand;
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

#ifndef USE_MPI
Result pam_distributed(int n, const vector<int>& D, int k, int seed, int , int ) {
    return pam_sequential(n, D, k, seed);
}
#else
Result pam_distributed(int n, const vector<int>& D, int k, int seed, int rank, int size) {
    Result res;

    if (k <= 0 || k > n) {
        if (rank == 0) cerr << "Erreur : k invalide" << endl;
        return res;
    }

    vector<int> localD;
    int rowsPerProc = n / size;
    int remainder = n % size;
    vector<int> sendCounts(size);
    vector<int> displs(size);
    int offset = 0;

    for (int i = 0; i < size; ++i) {
        sendCounts[i] = (i < remainder) ? (rowsPerProc + 1) * n : rowsPerProc * n;
        displs[i] = offset;
        offset += sendCounts[i];
    }

    int localRows = sendCounts[rank] / n;
    localD.resize(localRows * n);

    MPI_Scatterv(D.data(), sendCounts.data(), displs.data(), MPI_INT,
                 localD.data(), localRows * n, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> medoids;
    if (rank == 0) {
        medoids.reserve(k);
        std::mt19937 rng(seed);
        vector<int> perm(n);
        for (int i = 0; i < n; ++i) perm[i] = i;
        shuffle(perm.begin(), perm.end(), rng);
        for (int i = 0; i < k; ++i) medoids.push_back(perm[i]);
    }

    medoids.resize(k);
    MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> membership(localRows);
    vector<int> bestDist(localRows), secondBestDist(localRows);

    calculer_affectation(n, localD, medoids, 0, localRows, membership, bestDist, secondBestDist);

    long long localCost = 0;
    #pragma omp parallel for reduction(+:localCost)
    for (int i = 0; i < localRows; ++i) localCost += bestDist[i];

    bool changed = true;
    while (changed) {
        changed = false;
        
        set<int> medset(medoids.begin(), medoids.end());
        vector<int> candidates;
        for(int c = 0; c < n; ++c) {
            if(!medset.count(c)) candidates.push_back(c);
        }
        int num_candidates = candidates.size();
        
        vector<long long> all_deltas(k * num_candidates, 0);
        
        #pragma omp parallel for collapse(2) schedule(dynamic, 16)
        for(int mi = 0; mi < k; ++mi) {
            for(int c_idx = 0; c_idx < num_candidates; ++c_idx) {
                int cand = candidates[c_idx];
                long long delta = 0;
                
                for(int i = 0; i < localRows; ++i) {
                    int distToCand = localD[i * n + cand];
                    
                    if(membership[i] == mi) {
                        int nd = min(distToCand, secondBestDist[i]);
                        delta += (long long)nd - bestDist[i];
                    } else {
                        if(distToCand < bestDist[i]) {
                            delta += (long long)distToCand - bestDist[i];
                        }
                    }
                }
                
                all_deltas[mi * num_candidates + c_idx] = delta;
            }
        }
        
        vector<long long> global_deltas(k * num_candidates);
        MPI_Allreduce(all_deltas.data(), global_deltas.data(), 
                      k * num_candidates, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        
        long long bestDelta = 0;
        int bestMedIndex = -1;
        int bestCand = -1;
        
        #pragma omp parallel
        {
            long long thread_best_delta = 0;
            int thread_best_mi = -1;
            int thread_best_cand = -1;
            
            #pragma omp for collapse(2)
            for(int mi = 0; mi < k; ++mi) {
                for(int c_idx = 0; c_idx < num_candidates; ++c_idx) {
                    long long delta = global_deltas[mi * num_candidates + c_idx];
                    
                    if(delta < thread_best_delta) {
                        thread_best_delta = delta;
                        thread_best_mi = mi;
                        thread_best_cand = candidates[c_idx];
                    }
                }
            }
            
            #pragma omp critical
            {
                if(thread_best_delta < bestDelta) {
                    bestDelta = thread_best_delta;
                    bestMedIndex = thread_best_mi;
                    bestCand = thread_best_cand;
                }
            }
        }

        if (bestDelta < 0) {
            medoids[bestMedIndex] = bestCand;
            MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);
            calculer_affectation(n, localD, medoids, 0, localRows, membership, bestDist, secondBestDist);
            localCost += bestDelta;
            changed = true;
        }
    }

    vector<int> fullMembership;
    vector<int> recvCounts(size);
    vector<int> displs_membership(size);
    
    for (int i = 0; i < size; ++i) {
        recvCounts[i] = (i < remainder ? rowsPerProc + 1 : rowsPerProc);
    }
    
    displs_membership[0] = 0;
    for(int i = 1; i < size; ++i) {
        displs_membership[i] = displs_membership[i-1] + recvCounts[i-1];
    }

    if (rank == 0) fullMembership.resize(n);
    
    MPI_Gatherv(membership.data(), localRows, MPI_INT,
                fullMembership.data(), recvCounts.data(), displs_membership.data(),
                MPI_INT, 0, MPI_COMM_WORLD);

    long long totalCost = 0;
    MPI_Reduce(&localCost, &totalCost, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        res.medoids = medoids;
        res.membership = fullMembership;
        res.cost = totalCost;
    }
    return res;
}
#endif

} // namespace pam
