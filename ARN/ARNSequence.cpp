/**
 * @file ARNSequence.cpp
 * @brief Implémentation des fonctions pour manipuler des séquences d'ARN
 * @author Projet MPI
 * @date 2025
 */

#include "ARNSequence.hpp"
#include <fstream>
#include <algorithm>
#include <cstring>
#include <iomanip>

/**
 * @brief Calcule la distance de Levenshtein entre deux séquences
 * 
 * Utilise la programmation dynamique pour calculer le nombre minimum
 * d'éditions (insertions, suppressions, substitutions) nécessaires.
 */
int levenshteinDistance(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Cas particuliers
    if (m == 0) return n;
    if (n == 0) return m;
    
    // Allocation dynamique de la matrice DP
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    
    // Initialisation
    for (int i = 0; i <= m; i++) {
        dp[i][0] = i;
    }
    for (int j = 0; j <= n; j++) {
        dp[0][j] = j;
    }
    
    // Remplissage de la matrice
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int cost = (seq1[i - 1] == seq2[j - 1]) ? 0 : 1;
            dp[i][j] = min({
                dp[i - 1][j] + 1,      // suppression
                dp[i][j - 1] + 1,      // insertion
                dp[i - 1][j - 1] + cost // substitution
            });
        }
    }
    
    return dp[m][n];
}

/**
 * @brief Calcule la distance de Hamming entre deux séquences
 * 
 * Compte le nombre de positions où les deux séquences diffèrent.
 * N'est défini que pour des séquences de même longueur.
 */
int hammingDistance(const string& seq1, const string& seq2) {
    if (seq1.length() != seq2.length()) {
        return -1; // Erreur : longueurs différentes
    }
    
    int distance = 0;
    for (size_t i = 0; i < seq1.length(); i++) {
        if (seq1[i] != seq2[i]) {
            distance++;
        }
    }
    
    return distance;
}

/**
 * @brief Lit un fichier FASTA contenant des séquences d'ARN
 * 
 * Lit le fichier ligne par ligne et construit les séquences.
 * Les lignes commençant par '>' sont des labels.
 */
int readFASTAFile(const string& filename, vector<ARNSeq>& sequences) {
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Erreur : impossible d'ouvrir le fichier " << filename << endl;
        return -1;
    }
    
    string line;
    string currentLabel;
    string currentSequence;
    int seqId = 0;
    
    while (getline(file, line)) {
        // Ignorer les lignes vides
        if (line.empty()) {
            continue;
        }
        
        // Si c'est un label
        if (line[0] == '>') {
            // Sauvegarder la séquence précédente si elle existe
            if (!currentSequence.empty()) {
                ARNSeq seq;
                seq.id = seqId++;
                seq.sequence = currentSequence;
                seq.label = currentLabel;
                sequences.push_back(seq);
                currentSequence.clear();
            }
            
            // Nouveau label
            currentLabel = line.substr(1); // Enlever le '>'
        } else {
            // Ajouter la ligne à la séquence courante
            currentSequence += line;
        }
    }
    
    // Ne pas oublier la dernière séquence
    if (!currentSequence.empty()) {
        ARNSeq seq;
        seq.id = seqId++;
        seq.sequence = currentSequence;
        seq.label = currentLabel;
        sequences.push_back(seq);
    }
    
    file.close();
    
    cout << "Fichier FASTA lu avec succès : " << sequences.size() << " séquence(s)" << endl;
    return sequences.size();
}

/**
 * @brief Calcule la matrice de distances entre toutes les paires de séquences
 * 
 * Utilise la fonction de distance fournie pour calculer toutes les distances.
 * La matrice est symétrique et les diagonales sont à zéro.
 */
int* computeDistanceMatrix(const vector<ARNSeq>& sequences, 
                           int (*distanceFunc)(const string&, const string&)) {
    int n = sequences.size();
    int* distanceMatrix = new int[n * n];
    
    // Initialiser et calculer les distances
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                distanceMatrix[i * n + j] = 0;
            } else if (i < j) {
                // Calculer la distance
                int dist = distanceFunc(sequences[i].sequence, sequences[j].sequence);
                distanceMatrix[i * n + j] = dist;
                distanceMatrix[j * n + i] = dist; // Symétrie
            }
        }
    }
    
    return distanceMatrix;
}

/**
 * @brief Écrit un graphe au format DOT (Graphviz)
 * 
 * Crée une arête entre deux séquences si leur distance < epsilon.
 * Le poids de l'arête est la distance elle-même.
 */
int writeGraphDOT(const vector<ARNSeq>& sequences, int* distanceMatrix, 
                  int epsilon, const string& outputFile) {
    ofstream file(outputFile);
    
    if (!file.is_open()) {
        cerr << "Erreur : impossible d'ouvrir le fichier " << outputFile << endl;
        return -1;
    }
    
    int n = sequences.size();
    
    file << "graph ARN {" << endl;
    file << "  rankdir=LR;" << endl;
    
    // Ajouter les nœuds
    for (int i = 0; i < n; i++) {
        file << "  seq" << i << " [label=\"" << sequences[i].label << "\"];" << endl;
    }
    
    file << endl;
    
    // Ajouter les arêtes (uniquement si distance < epsilon)
    int edgeCount = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int dist = distanceMatrix[i * n + j];
            if (dist < epsilon) {
                file << "  seq" << i << " -- seq" << j << " [weight=" << dist << ", label=\"" << dist << "\"];" << endl;
                edgeCount++;
            }
        }
    }
    
    file << "}" << endl;
    file.close();
    
    cout << "Graphe écrit dans " << outputFile << endl;
    cout << "  Nœuds : " << n << ", Arêtes : " << edgeCount << endl;
    
    return 0;
}

/**
 * @brief Affiche une séquence d'ARN
 */
void printARNSeq(const ARNSeq& seq) {
    cout << "ID: " << seq.id << " | Label: " << seq.label 
         << " | Séquence: " << seq.sequence.substr(0, 50);
    if (seq.sequence.length() > 50) {
        cout << "...";
    }
    cout << " (taille: " << seq.sequence.length() << ")" << endl;
}

/**
 * @brief Affiche la matrice de distances
 */
void printDistanceMatrix(int* distanceMatrix, int n) {
    cout << "Matrice de distances (" << n << "x" << n << "):" << endl;
    
    // Afficher l'en-tête
    cout << "    ";
    for (int j = 0; j < n; j++) {
        cout << setw(6) << j;
    }
    cout << endl;
    
    // Afficher les lignes
    for (int i = 0; i < n; i++) {
        cout << setw(3) << i << " ";
        for (int j = 0; j < n; j++) {
            cout << setw(6) << distanceMatrix[i * n + j];
        }
        cout << endl;
    }
}
