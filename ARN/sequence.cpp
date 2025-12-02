#include "sequence.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <set>
#include <cmath>
#include <iomanip>

#define INF 1000

using namespace std;

// ============================================================================
// Lecture FASTA
// ============================================================================

vector<Sequence> read_fasta(const char* filename) {
    vector<Sequence> sequences;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Erreur: impossible d'ouvrir " << filename << endl;
        return sequences;
    }
    
    Sequence current;
    string line;
    
    while (getline(file, line)) {
        // Ignorer les lignes vides
        if (line.empty()) continue;
        
        // Nouvelle séquence
        if (line[0] == '>') {
            // Sauvegarder la séquence précédente si elle existe
            if (!current.id.empty()) {
                sequences.push_back(current);
            }
            // Commencer une nouvelle séquence
            current.id = line.substr(1);  // Enlever le '>'
            current.data.clear();
        } else {
            // Ajouter à la séquence courante
            current.data += line;
        }
    }
    
    // Ajouter la dernière séquence
    if (!current.id.empty()) {
        sequences.push_back(current);
    }
    
    file.close();
    
    cout << "Lecture de " << sequences.size() << " séquences depuis " << filename << endl;
    return sequences;
}

// ============================================================================
// Distance de Hamming
// ============================================================================

int distance_hamming(const Sequence& s1, const Sequence& s2) {
    if (s1.data.length() != s2.data.length()) {
        return INF;  // Séquences de longueurs différentes
    }
    
    int dist = 0;
    for (size_t i = 0; i < s1.data.length(); ++i) {
        if (s1.data[i] != s2.data[i]) {
            dist++;
        }
    }
    return dist;
}

// ============================================================================
// Distance d'édition (Levenshtein)
// ============================================================================

int distance_edit(const Sequence& s1, const Sequence& s2) {
    int m = s1.data.length();
    int n = s2.data.length();
    
    // Tableau de programmation dynamique
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    
    // Initialisation
    for (int i = 0; i <= m; ++i) {
        dp[i][0] = i;  // Suppression de i caractères
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = j;  // Insertion de j caractères
    }
    
    // Remplissage du tableau
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int cost = (s1.data[i-1] == s2.data[j-1]) ? 0 : 1;
            
            dp[i][j] = min({
                dp[i-1][j] + 1,        // Suppression
                dp[i][j-1] + 1,        // Insertion
                dp[i-1][j-1] + cost    // Substitution
            });
        }
    }
    
    return dp[m][n];
}

// ============================================================================
// Distance k-mer
// ============================================================================

int distance_kmer(const Sequence& s1, const Sequence& s2, int k) {
    if ((int)s1.data.length() < k || (int)s2.data.length() < k) {
        return 100;  // Distance maximale
    }
    
    // Extraction des k-mers
    set<string> kmers1, kmers2;
    
    for (size_t i = 0; i <= s1.data.length() - k; ++i) {
        kmers1.insert(s1.data.substr(i, k));
    }
    
    for (size_t i = 0; i <= s2.data.length() - k; ++i) {
        kmers2.insert(s2.data.substr(i, k));
    }
    
    // Calcul de l'intersection
    set<string> intersection;
    set_intersection(kmers1.begin(), kmers1.end(),
                    kmers2.begin(), kmers2.end(),
                    inserter(intersection, intersection.begin()));
    
    // Calcul de l'union
    set<string> union_set;
    set_union(kmers1.begin(), kmers1.end(),
             kmers2.begin(), kmers2.end(),
             inserter(union_set, union_set.begin()));
    
    // Distance de Jaccard (1 - similarité) * 100
    if (union_set.empty()) return 100;
    
    double similarity = (double)intersection.size() / union_set.size();
    return (int)((1.0 - similarity) * 100);
}

// ============================================================================
// Construction de la matrice de distances
// ============================================================================

int* build_distance_matrix(const vector<Sequence>& seqs, 
                          const string& dist_type,
                          int k) {
    int n = seqs.size();
    int* D = new int[n * n];
    
    cout << "Construction de la matrice de distances (" << dist_type << ")..." << endl;
    
    // Fonction de distance à utiliser
    auto dist_func = distance_edit;
    if (dist_type == "hamming") {
        dist_func = distance_hamming;
    }
    
    int total_pairs = n * (n - 1) / 2;
    int computed = 0;
    
    for (int i = 0; i < n; ++i) {
        D[i*n + i] = 0;  // Distance à soi-même = 0
        
        for (int j = i + 1; j < n; ++j) {
            int dist;
            
            if (dist_type == "kmer") {
                dist = distance_kmer(seqs[i], seqs[j], k);
            } else {
                dist = dist_func(seqs[i], seqs[j]);
            }
            
            D[i*n + j] = dist;
            D[j*n + i] = dist;  // Symétrie
            
            computed++;
            if (computed % 1000 == 0 || computed == total_pairs) {
                cout << "  Progression: " << computed << "/" << total_pairs 
                     << " (" << (100*computed/total_pairs) << "%)\r" << flush;
            }
        }
    }
    cout << endl;
    
    return D;
}

// ============================================================================
// Filtrage avec epsilon
// ============================================================================

int* filter_graph(int* D, int n, int epsilon) {
    int* adjacence = new int[n * n]();  // Initialisation à 0
    
    int num_edges = 0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (D[i*n + j] < epsilon) {
                adjacence[i*n + j] = D[i*n + j];
                adjacence[j*n + i] = D[i*n + j];
                num_edges++;
            }
        }
    }
    
    cout << "Graphe filtré avec epsilon=" << epsilon << ": " 
         << num_edges << " arêtes créées" << endl;
    
    return adjacence;
}

// ============================================================================
// Génération de séquences de test
// ============================================================================

void generate_test_sequences(int n, int length, const char* filename, int num_families) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Erreur: impossible de créer " << filename << endl;
        return;
    }
    
    const char bases[] = "ACGT";
    mt19937 rng(42);
    uniform_int_distribution<int> base_dist(0, 3);
    uniform_int_distribution<int> mutation_dist(0, 99);
    
    // Générer des séquences "template" pour chaque famille
    vector<string> templates;
    for (int f = 0; f < num_families; ++f) {
        string tmpl;
        for (int j = 0; j < length; ++j) {
            tmpl += bases[base_dist(rng)];
        }
        templates.push_back(tmpl);
    }
    
    // Générer n séquences basées sur les templates
    for (int i = 0; i < n; ++i) {
        int family = i % num_families;
        string seq = templates[family];
        
        // Introduire des mutations aléatoires (10% de chance par base)
        for (size_t j = 0; j < seq.length(); ++j) {
            if (mutation_dist(rng) < 10) {
                seq[j] = bases[base_dist(rng)];
            }
        }
        
        // Écrire au format FASTA
        file << ">seq" << i << "_family" << family << "\n";
        
        // Découper en lignes de 80 caractères (standard FASTA)
        for (size_t j = 0; j < seq.length(); j += 80) {
            file << seq.substr(j, 80) << "\n";
        }
    }
    
    file.close();
    cout << "Généré " << n << " séquences dans " << filename 
         << " (" << num_families << " familles)" << endl;
}

// ============================================================================
// Statistiques sur la matrice de distances
// ============================================================================

void print_distance_stats(int* D, int n) {
    vector<int> distances;
    
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (D[i*n + j] < INF) {
                distances.push_back(D[i*n + j]);
            }
        }
    }
    
    if (distances.empty()) {
        cout << "Aucune distance valide." << endl;
        return;
    }
    
    sort(distances.begin(), distances.end());
    
    int min_dist = distances.front();
    int max_dist = distances.back();
    int median = distances[distances.size() / 2];
    
    double mean = 0;
    for (int d : distances) mean += d;
    mean /= distances.size();
    
    cout << "\n=== Statistiques des distances ===" << endl;
    cout << "  Nombre de paires: " << distances.size() << endl;
    cout << "  Min: " << min_dist << endl;
    cout << "  Max: " << max_dist << endl;
    cout << "  Médiane: " << median << endl;
    cout << "  Moyenne: " << fixed << setprecision(2) << mean << endl;
    
    // Quartiles
    int q1 = distances[distances.size() / 4];
    int q3 = distances[3 * distances.size() / 4];
    cout << "  Q1 (25%): " << q1 << endl;
    cout << "  Q3 (75%): " << q3 << endl;
    cout << endl;
}

// ============================================================================
// Export des résultats
// ============================================================================

void export_clustering_results(const vector<Sequence>& sequences,
                               const vector<int>& medoids,
                               const vector<int>& membership,
                               long long cost,
                               const char* filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Erreur: impossible de créer " << filename << endl;
        return;
    }
    
    file << "=== Résultats du Clustering PAM ===" << endl;
    file << "Coût total: " << cost << endl;
    file << "Nombre de clusters: " << medoids.size() << endl;
    file << endl;
    
    // Compter les membres de chaque cluster
    vector<int> counts(medoids.size(), 0);
    for (int m : membership) {
        if (m >= 0 && m < (int)medoids.size()) {
            counts[m]++;
        }
    }
    
    // Détails par cluster
    for (size_t m = 0; m < medoids.size(); ++m) {
        file << "--- Cluster " << m << " ---" << endl;
        file << "Médoïde: " << sequences[medoids[m]].id << endl;
        file << "Taille: " << counts[m] << " séquences" << endl;
        file << "Membres:" << endl;
        
        for (size_t i = 0; i < membership.size(); ++i) {
            if (membership[i] == (int)m) {
                file << "  - " << sequences[i].id << endl;
            }
        }
        file << endl;
    }
    
    file.close();
    cout << "Résultats exportés dans " << filename << endl;
}