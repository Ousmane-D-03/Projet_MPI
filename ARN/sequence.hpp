#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <string>
#include <vector>
#include <map>

/**
 * @brief Structure représentant une séquence d'ARN
 */
struct Sequence {
    std::string id;      // Identifiant de la séquence (ex: >seq1)
    std::string data;    // Données de la séquence (ACGT)
};

/**
 * @brief Lecture d'un fichier FASTA contenant des séquences d'ARN
 * @param filename Chemin du fichier FASTA
 * @return Vector de séquences
 */
std::vector<Sequence> read_fasta(const char* filename);

/**
 * @brief Distance de Hamming entre deux séquences (même longueur requise)
 * @param s1 Première séquence
 * @param s2 Deuxième séquence
 * @return Distance de Hamming (INF si longueurs différentes)
 */
int distance_hamming(const Sequence& s1, const Sequence& s2);

/**
 * @brief Distance d'édition (Levenshtein) entre deux séquences
 * @param s1 Première séquence
 * @param s2 Deuxième séquence
 * @return Distance d'édition
 */
int distance_edit(const Sequence& s1, const Sequence& s2);

/**
 * @brief Distance basée sur les k-mers
 * @param s1 Première séquence
 * @param s2 Deuxième séquence
 * @param k Taille des k-mers
 * @return Distance (0-100, multiplié par 100 pour avoir des entiers)
 */
int distance_kmer(const Sequence& s1, const Sequence& s2, int k);

/**
 * @brief Construction de la matrice de distances entre toutes les séquences
 * @param seqs Vector de séquences
 * @param dist_type Type de distance: "hamming", "edit", "kmer"
 * @param k Paramètre k pour distance kmer (ignoré pour autres distances)
 * @return Matrice de distances (n*n, allocation dynamique)
 */
int* build_distance_matrix(const std::vector<Sequence>& seqs, 
                          const std::string& dist_type = "edit",
                          int k = 3);

/**
 * @brief Filtrage de la matrice de distances pour créer une matrice d'adjacence
 * @param D Matrice de distances complète
 * @param n Nombre de séquences
 * @param epsilon Seuil de distance pour créer une arête
 * @return Matrice d'adjacence (arête existe si distance < epsilon)
 */
int* filter_graph(int* D, int n, int epsilon);

/**
 * @brief Génération de séquences aléatoires pour test
 * @param n Nombre de séquences
 * @param length Longueur de chaque séquence
 * @param filename Fichier de sortie FASTA
 * @param num_families Nombre de familles (pour clustering)
 */
void generate_test_sequences(int n, int length, const char* filename, int num_families = 3);

/**
 * @brief Affichage des statistiques sur la matrice de distances
 * @param D Matrice de distances
 * @param n Taille de la matrice
 */
void print_distance_stats(int* D, int n);

/**
 * @brief Export des résultats de clustering
 * @param sequences Vector de séquences
 * @param medoids Indices des médoïdes
 * @param membership Affectation de chaque séquence
 * @param cost Coût total du clustering
 * @param filename Fichier de sortie
 */
void export_clustering_results(const std::vector<Sequence>& sequences,
                               const std::vector<int>& medoids,
                               const std::vector<int>& membership,
                               long long cost,
                               const char* filename);

#endif