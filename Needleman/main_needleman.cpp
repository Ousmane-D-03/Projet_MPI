/**
 * @file main_needleman.cpp
 * @brief Programme de test pour Needleman-Wunsch
 */

#include "Needleman.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>

using namespace std;

int main() {
    cout << "Test Needleman-Wunsch\n" << endl;
    
    string seq1 = "ACGTACGT";
    string seq2 = "ACGTTAGC";
    
    ScoringParams params;
    
    cout << "Seq1: " << seq1 << endl;
    cout << "Seq2: " << seq2 << "\n" << endl;
    
    // Séquentiel
    auto start = chrono::high_resolution_clock::now();
    int score_seq = needleman_wunsch_sequential(seq1, seq2, params);
    auto end = chrono::high_resolution_clock::now();
    double time_seq = chrono::duration<double>(end - start).count();
    
    cout << "Sequential: score=" << score_seq 
         << ", time=" << fixed << setprecision(6) << time_seq << "s" << endl;
    
    // Parallèle
    start = chrono::high_resolution_clock::now();
    int score_par = needleman_wunsch_parallel(seq1, seq2, params, 4);
    end = chrono::high_resolution_clock::now();
    double time_par = chrono::duration<double>(end - start).count();
    
    cout << "Parallel:   score=" << score_par 
         << ", time=" << time_par << "s" << endl;
    
    if (score_seq == score_par) {
        cout << "\nOK: same result" << endl;
    } else {
        cout << "\nERROR: different results!" << endl;
    }
    
    return 0;
}
