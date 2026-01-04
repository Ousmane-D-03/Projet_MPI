/**
 * @file Needleman.hpp
 * @brief Algorithme de Needleman-Wunsch pour l'alignement global de séquences
 * @author Projet MPI - Partie 2
 * @date 2025
 */

#ifndef NEEDLEMAN_HPP
#define NEEDLEMAN_HPP

#include <string>
#include <vector>

using namespace std;

/**
 * @brief Structure des paramètres de scoring
 */
struct ScoringParams {
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    
    ScoringParams() : match(1), mismatch(-1), gap_open(-3), gap_extend(-1) {}
    ScoringParams(int m, int mm, int go, int ge) 
        : match(m), mismatch(mm), gap_open(go), gap_extend(ge) {}
};

/**
 * @brief Calcule le score d'alignement (version séquentielle)
 * @param seq1 Première séquence
 * @param seq2 Deuxième séquence
 * @param params Paramètres de scoring
 * @return Score d'alignement optimal
 */
int needleman_wunsch_sequential(const string& seq1, const string& seq2, 
                                const ScoringParams& params = ScoringParams());

/**
 * @brief Calcule le score d'alignement (version parallèle OpenMP)
 * @param seq1 Première séquence
 * @param seq2 Deuxième séquence
 * @param params Paramètres de scoring
 * @param num_threads Nombre de threads (0 = défaut)
 * @return Score d'alignement optimal
 */
int needleman_wunsch_parallel(const string& seq1, const string& seq2,
                              const ScoringParams& params = ScoringParams(),
                              int num_threads = 0);

#endif // NEEDLEMAN_HPP
