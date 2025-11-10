// PAM.cpp
#include "PAM.hpp"

#include <algorithm>
#include <random>
#include <limits>
#include <set>
#include <iostream>

using namespace std;

namespace pam {

static void compute_membership_and_dists(int n, const vector<int>& D, const vector<int>& medoids,
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
    compute_membership_and_dists(n, D, medoids, 0, n, membership, bestDist, secondBestDist);

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
            compute_membership_and_dists(n, D, medoids, 0, n, membership, bestDist, secondBestDist);
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
// Fallback: if MPI not enabled, just run sequential PAM
Result pam_replicated(int n, const vector<int>& D, int k, int seed, int /*rank*/, int /*size*/) {
    return pam_sequential(n, D, k, seed);
}
#else
// MPI-enabled implementation
#include <mpi.h>

Result pam_replicated(int n, const vector<int>& D, int k, int seed, int rank, int size) {
    // We distribute work by rows for computing deltas but D is fully available on each process.
    // We'll perform the same swap search but each process computes contributions for a subset of rows.
    Result res;
    if (k <= 0 || k > n) {
        if (rank == 0) cerr << "Invalid k" << endl;
        return res;
    }

    // initial medoids on rank 0
    vector<int> medoids;
    if (rank == 0) {
        vector<int> idx(n);
        for (int i = 0; i < n; ++i) idx[i] = i;
        mt19937 rng(seed);
        shuffle(idx.begin(), idx.end(), rng);
        for (int i = 0; i < k; ++i) medoids.push_back(idx[i]);
    }
    // broadcast medoids
    medoids.resize(k);
    MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);

    // prepare local buffers: partition rows among processes
    int rows_per = n / size;
    int rem = n % size;
    int start = rank * rows_per + min(rank, rem);
    int end = start + rows_per + (rank < rem ? 1 : 0);
    int local_n = end - start;

    vector<int> membership_local(local_n);
    vector<int> bestDist_local(local_n), secondBestDist_local(local_n);

    // initial membership and dists (local)
    compute_membership_and_dists(n, D, medoids, start, end, membership_local, bestDist_local, secondBestDist_local);

    // compute initial cost (reduce)
    long long local_cost = 0;
    for (int i = 0; i < local_n; ++i) local_cost += bestDist_local[i];
    long long cost = 0;
    MPI_Allreduce(&local_cost, &cost, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    bool improved = true;
    while (improved) {
        improved = false;
        long long bestDeltaGlobal = 0;
        int bestMedGlobal = -1;
        int bestCandGlobal = -1;

        set<int> medset(medoids.begin(), medoids.end());

        // For each medoid and candidate, compute local delta and reduce
        for (int mi = 0; mi < k; ++mi) {
            for (int cand = 0; cand < n; ++cand) {
                if (medset.count(cand)) continue;
                long long local_delta = 0;
                for (int ii = start; ii < end; ++ii) {
                    int i_local = ii - start;
                    int distToCand = D[ii*n + cand];
                    if (membership_local[i_local] == mi) {
                        int newdist = min(distToCand, secondBestDist_local[i_local]);
                        local_delta += (long long)newdist - bestDist_local[i_local];
                    } else {
                        if (distToCand < bestDist_local[i_local]) local_delta += (long long)distToCand - bestDist_local[i_local];
                    }
                }
                long long total_delta = 0;
                MPI_Allreduce(&local_delta, &total_delta, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                if (total_delta < bestDeltaGlobal) {
                    bestDeltaGlobal = total_delta;
                    bestMedGlobal = mi;
                    bestCandGlobal = cand;
                }
            }
        }

        if (bestDeltaGlobal < 0) {
            improved = true;
            // apply swap on all processes
            medoids[bestMedGlobal] = bestCandGlobal;
            MPI_Bcast(medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);
            // recompute local membership/dists
            compute_membership_and_dists(n, D, medoids, start, end, membership_local, bestDist_local, secondBestDist_local);
            // update cost
            local_cost = 0;
            for (int i = 0; i < local_n; ++i) local_cost += bestDist_local[i];
            MPI_Allreduce(&local_cost, &cost, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        }
    }

    // gather memberships to root (optional) - we return local membership only aggregated on rank 0
    if (rank == 0) {
        res.membership.assign(n, -1);
        // copy root part
        for (int i = start; i < end; ++i) res.membership[i] = membership_local[i - start];
        // receive from others
        for (int p = 1; p < size; ++p) {
            int pstart = p * rows_per + min(p, rem);
            int pend = pstart + rows_per + (p < rem ? 1 : 0);
            int plen = pend - pstart;
            if (plen > 0) {
                vector<int> buf(plen);
                MPI_Recv(buf.data(), plen, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < plen; ++i) res.membership[pstart + i] = buf[i];
            }
        }
    } else {
        // send local membership to root
        if (local_n > 0) MPI_Send(membership_local.data(), local_n, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // gather medoids to all (already consistent)
    res.medoids = medoids;
    res.cost = cost;
    return res;
}
#endif

} // namespace pam
