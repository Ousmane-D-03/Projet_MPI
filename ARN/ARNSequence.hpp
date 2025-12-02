/**
 * @file ARNSequence.hpp
 * @brief Structures et fonctions pour manipuler des séquences d'ARN
 * @author Projet MPI
 * @date 2025
 */

#ifndef ARNSEQUENCE_HPP
#define ARNSEQUENCE_HPP

#include <string>
#include <vector>
#include <iostream>

using namespace std;

/**
 * @brief Représentation d'une séquence d'ARN
 * 
 * Une séquence d'ARN est une suite de nucléotides (A, C, G, T/U)
 */
struct ARNSeq {
    int id;           /**< Identifiant unique de la séquence */
    string sequence;  /**< La séquence elle-même (A, C, G, T) */
    string label;     /**< Étiquette/nom de la séquence */
};

/**
 * @brief Calcule la distance de Levenshtein (edit distance) entre deux séquences
 * 
 * La distance de Levenshtein est le nombre minimum d'éditions (insertions, 
 * suppressions, substitutions) nécessaires pour transformer une séquence en une autre.
 * 
 * @param seq1 Première séquence d'ARN
 * @param seq2 Deuxième séquence d'ARN
 * @return int La distance de Levenshtein entre seq1 et seq2
 * 
 * @note Complexité : O(|seq1| * |seq2|)
 */
int levenshteinDistance(const string& seq1, const string& seq2);

/**
 * @brief Calcule la distance de Hamming entre deux séquences de même longueur
 * 
 * La distance de Hamming est le nombre de positions où les nucléotides diffèrent.
 * Cette distance n'est définie que pour des séquences de même longueur.
 * 
 * @param seq1 Première séquence d'ARN
 * @param seq2 Deuxième séquence d'ARN
 * @return int La distance de Hamming, ou -1 si les séquences ont des longueurs différentes
 * 
 * @note Complexité : O(max(|seq1|, |seq2|))
 */
int hammingDistance(const string& seq1, const string& seq2);

/**
 * @brief Lit un fichier contenant des séquences d'ARN (format FASTA)
 * 
 * Format FASTA attendu :
 * >label1
 * ACGTACGT...
 * >label2
 * ACGTACGT...
 * 
 * @param filename Chemin du fichier FASTA
 * @param sequences Vecteur de séquences (sortie)
 * @return int Nombre de séquences lues, ou -1 en cas d'erreur
 */
int readFASTAFile(const string& filename, vector<ARNSeq>& sequences);

/**
 * @brief Crée une matrice de distances entre toutes les paires de séquences
 * 
 * @param sequences Vecteur de séquences d'ARN
 * @param distanceFunc Fonction de distance à utiliser (e.g., levenshteinDistance)
 * @return int* Matrice n x n contenant les distances, où n est le nombre de séquences
 * 
 * @note La matrice est stockée de manière linéaire en mémoire (row-major)
 * @note La matrice est symétrique : mat[i*n + j] == mat[j*n + i]
 */
int* computeDistanceMatrix(const vector<ARNSeq>& sequences, 
                           int (*distanceFunc)(const string&, const string&));

/**
 * @brief Construit un graphe à partir des séquences en utilisant un seuil epsilon
 * 
 * Une arête existe entre deux séquences si leur distance est strictement inférieure à epsilon.
 * Les poids des arêtes sont les distances elles-mêmes.
 * 
 * @param sequences Vecteur de séquences d'ARN
 * @param distanceMatrix Matrice de distances entre les séquences
 * @param epsilon Seuil de distance pour créer une arête
 * @param outputFile Fichier de sortie au format .dot (Graphviz)
 * @return int 0 si succès, -1 en cas d'erreur
 * 
 * @note Le fichier généré peut être utilisé directement avec Floyd-Warshall
 */
int writeGraphDOT(const vector<ARNSeq>& sequences, int* distanceMatrix, 
                  int epsilon, const string& outputFile);

/**
 * @brief Affiche une séquence d'ARN
 * 
 * @param seq La séquence à afficher
 */
void printARNSeq(const ARNSeq& seq);

/**
 * @brief Affiche la matrice de distances
 * 
 * @param distanceMatrix Pointeur vers la matrice
 * @param n Nombre de séquences (taille de la matrice n x n)
 */
void printDistanceMatrix(int* distanceMatrix, int n);

#endif // ARNSEQUENCE_HPP
