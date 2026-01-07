/**
 * @file main_needleman.cpp
 * @brief Programme de test pour Needleman-Wunsch
 */

#include "Needleman.hpp"
#include "../Floyd/ForGraph.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <vector>
#include <omp.h>

using namespace std;

vector<string> Lecture(const string& file_name)
{
    vector<string> V;
    ifstream f(file_name);
    if (!f) {
        cerr << "Impossible d’ouvrir le fichier : " << file_name << endl;
        return V;
    }

    string line, seq;

    while (getline(f, line)) {

        if (line.empty()) continue;

        if (line[0] == '>') {          // nouvelle séquence
            if (!seq.empty()) {
                V.push_back(seq);
                seq.clear();
            }
        } 
        else {
            seq += line;
        }
    }

    if (!seq.empty())
        V.push_back(seq);

    return V;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " fichier.fasta\n";
        return 1;
    }

    vector<string> seqs = Lecture(argv[1]);

    if (seqs.size() < 2) {
        cerr << "Erreur : le fichier doit contenir au moins 2 séquences\n";
        return 1;
    }

    string seq1 = seqs[0];
    string seq2 = seqs[1];

    cout << "Test Needleman-Wunsch\n\n";
    cout << "Seq1: " << seq1 << endl;
    cout << "Seq2: " << seq2 << "\n" << endl;

    ScoringParams params;

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

    if (score_seq == score_par)
        cout << "\nOK: same result" << endl;
    else
        cout << "\nERROR: different results!" << endl;

    return 0;
}
