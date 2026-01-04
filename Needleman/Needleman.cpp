/**
 * @file Needleman.cpp
 * @brief Implémentation de l'algorithme de Needleman-Wunsch
 */

#include "Needleman.hpp"
#include <algorithm>
#include <omp.h>

using namespace std;

inline int substitution_score(char a, char b, const ScoringParams& params) {
    return (a == b) ? params.match : params.mismatch;
}

int needleman_wunsch_sequential(const string& seq1, const string& seq2, 
                                const ScoringParams& params) {
    int m = seq1.length();
    int n = seq2.length();
    
    if (m == 0) return n * params.gap_open;
    if (n == 0) return m * params.gap_open;
    
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    vector<vector<int>> gap_state(m + 1, vector<int>(n + 1, 0));
    
    dp[0][0] = 0;
    for (int i = 1; i <= m; i++) {
        dp[i][0] = params.gap_open + (i - 1) * params.gap_extend;
        gap_state[i][0] = 1;
    }
    for (int j = 1; j <= n; j++) {
        dp[0][j] = params.gap_open + (j - 1) * params.gap_extend;
        gap_state[0][j] = 2;
    }
    
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match_score = dp[i-1][j-1] + substitution_score(seq1[i-1], seq2[j-1], params);
            
            int gap_penalty_vertical = (gap_state[i-1][j] == 1) 
                ? params.gap_extend : params.gap_open;
            int gap_seq2 = dp[i-1][j] + gap_penalty_vertical;
            
            int gap_penalty_horizontal = (gap_state[i][j-1] == 2)
                ? params.gap_extend : params.gap_open;
            int gap_seq1 = dp[i][j-1] + gap_penalty_horizontal;
            
            dp[i][j] = max({match_score, gap_seq2, gap_seq1});
            
            if (dp[i][j] == gap_seq2) {
                gap_state[i][j] = 1;
            } else if (dp[i][j] == gap_seq1) {
                gap_state[i][j] = 2;
            } else {
                gap_state[i][j] = 0;
            }
        }
    }
    
    return dp[m][n];
}

int needleman_wunsch_parallel(const string& seq1, const string& seq2,
                              const ScoringParams& params,
                              int num_threads) {
    int m = seq1.length();
    int n = seq2.length();
    
    if (m == 0) return n * params.gap_open;
    if (n == 0) return m * params.gap_open;
    
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    vector<vector<int>> gap_state(m + 1, vector<int>(n + 1, 0));
    
    dp[0][0] = 0;
    for (int i = 1; i <= m; i++) {
        dp[i][0] = params.gap_open + (i - 1) * params.gap_extend;
        gap_state[i][0] = 1;
    }
    for (int j = 1; j <= n; j++) {
        dp[0][j] = params.gap_open + (j - 1) * params.gap_extend;
        gap_state[0][j] = 2;
    }
    
    int total_diagonals = m + n - 1;
    
    for (int diag = 0; diag < total_diagonals; diag++) {
        int start_i = max(0, diag - n + 1);
        int end_i = min(m, diag + 1);
        
        #pragma omp parallel for schedule(static)
        for (int i = start_i; i < end_i; i++) {
            int j = diag - i + 1;
            
            if (j > 0 && j <= n) {
                int match_score = dp[i-1][j-1] + 
                    substitution_score(seq1[i-1], seq2[j-1], params);
                
                int gap_penalty_vertical = (gap_state[i-1][j] == 1)
                    ? params.gap_extend : params.gap_open;
                int gap_seq2 = dp[i-1][j] + gap_penalty_vertical;
                
                int gap_penalty_horizontal = (gap_state[i][j-1] == 2)
                    ? params.gap_extend : params.gap_open;
                int gap_seq1 = dp[i][j-1] + gap_penalty_horizontal;
                
                dp[i][j] = max({match_score, gap_seq2, gap_seq1});
                
                if (dp[i][j] == gap_seq2) {
                    gap_state[i][j] = 1;
                } else if (dp[i][j] == gap_seq1) {
                    gap_state[i][j] = 2;
                } else {
                    gap_state[i][j] = 0;
                }
            }
        }
    }
    
    return dp[m][n];
}
