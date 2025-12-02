#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "sequence.hpp"
#include "../Floyd/ForGraph.hpp"
#include "../PAM/PAM.hpp"

using namespace std;

void print_usage(const char* prog) {
    cout << "Usage: " << prog << " [options]\n\n";
    cout << "Options:\n";
    cout << "  --generate <n> <len> <output.fasta> [families]  Générer des séquences test\n";
    cout << "  --input <file.fasta>                             Fichier FASTA d'entrée\n";
    cout << "  --distance <type>                                Type de distance: hamming, edit, kmer (défaut: edit)\n";
    cout << "  --kmer <k>                                       Taille des k-mers (défaut: 3)\n";
    cout << "  --epsilon <value>                                Seuil pour filtrage du graphe (défaut: INF = pas de filtrage)\n";
    cout << "  --clusters <k>                                   Nombre de clusters pour PAM\n";
    cout << "  --seed <value>                                   Graine aléatoire (défaut: 12345)\n";
    cout << "  --output <file>                                  Fichier de sortie des résultats\n";
    cout << "  --no-floyd                                       Ne pas appliquer Floyd-Warshall\n";
    cout << "  --help                                           Afficher cette aide\n";
    cout << "\nExemples:\n";
    cout << "  # Générer 30 séquences de 100 bases (3 familles)\n";
    cout << "  " << prog << " --generate 30 100 test.fasta 3\n\n";
    cout << "  # Clustering avec distance d'édition\n";
    cout << "  " << prog << " --input test.fasta --distance edit --clusters 3\n\n";
    cout << "  # Avec filtrage du graphe (epsilon=20) et Floyd-Warshall\n";
    cout << "  " << prog << " --input test.fasta --epsilon 20 --clusters 3\n";
}

int main(int argc, char* argv[]) {
    int rank = 0, size = 1;
    
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    // Paramètres par défaut
    string input_file = "";
    string output_file = "clustering_results.txt";
    string distance_type = "edit";
    int k_kmer = 3;
    int epsilon = 1000;  // INF par défaut
    int k_clusters = 3;
    int seed = 12345;
    bool use_floyd = true;
    bool generate_mode = false;
    int gen_n = 0, gen_len = 0, gen_families = 3;
    string gen_output = "";

    // Parsing des arguments
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            if (rank == 0) print_usage(argv[0]);
#ifdef USE_MPI
            MPI_Finalize();
#endif
            return 0;
        }
        else if (arg == "--generate") {
            if (i + 3 >= argc) {
                if (rank == 0) cerr << "Erreur: --generate nécessite <n> <len> <output>\n";
#ifdef USE_MPI
                MPI_Finalize();
#endif
                return 1;
            }
            generate_mode = true;
            gen_n = stoi(argv[++i]);
            gen_len = stoi(argv[++i]);
            gen_output = argv[++i];
            if (i + 1 < argc && argv[i+1][0] != '-') {
                gen_families = stoi(argv[++i]);
            }
        }
        else if (arg == "--input" && i + 1 < argc) {
            input_file = argv[++i];
        }
        else if (arg == "--distance" && i + 1 < argc) {
            distance_type = argv[++i];
        }
        else if (arg == "--kmer" && i + 1 < argc) {
            k_kmer = stoi(argv[++i]);
        }
        else if (arg == "--epsilon" && i + 1 < argc) {
            epsilon = stoi(argv[++i]);
        }
        else if (arg == "--clusters" && i + 1 < argc) {
            k_clusters = stoi(argv[++i]);
        }
        else if (arg == "--seed" && i + 1 < argc) {
            seed = stoi(argv[++i]);
        }
        else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
        else if (arg == "--no-floyd") {
            use_floyd = false;
        }
    }

    // Mode génération
    if (generate_mode) {
        if (rank == 0) {
            generate_test_sequences(gen_n, gen_len, gen_output.c_str(), gen_families);
        }
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 0;
    }

    // Mode clustering
    if (input_file.empty()) {
        if (rank == 0) {
            cerr << "Erreur: fichier d'entrée requis (--input)\n";
            print_usage(argv[0]);
        }
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if (rank == 0) {
        cout << "\n=== Clustering de Séquences ARN ===" << endl;
        cout << "Fichier d'entrée: " << input_file << endl;
        cout << "Distance: " << distance_type;
        if (distance_type == "kmer") cout << " (k=" << k_kmer << ")";
        cout << endl;
        cout << "Epsilon: " << (epsilon < 1000 ? to_string(epsilon) : "INF (pas de filtrage)") << endl;
        cout << "Nombre de clusters: " << k_clusters << endl;
        cout << "Floyd-Warshall: " << (use_floyd ? "Oui" : "Non") << endl;
        cout << "Seed: " << seed << endl;
#ifdef USE_MPI
        cout << "Processus MPI: " << size << endl;
#endif
        cout << endl;
    }

    // ========================================================================
    // ÉTAPE 1: Lecture des séquences
    // ========================================================================
    
    vector<Sequence> sequences;
    if (rank == 0) {
        sequences = read_fasta(input_file.c_str());
        if (sequences.empty()) {
            cerr << "Erreur: aucune séquence lue\n";
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, 1);
#else
            return 1;
#endif
        }
    }

    int n = 0;
    if (rank == 0) {
        n = sequences.size();
    }

#ifdef USE_MPI
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    // ========================================================================
    // ÉTAPE 2: Calcul de la matrice de distances
    // ========================================================================
    
    int* D = nullptr;
    
    if (rank == 0) {
        D = build_distance_matrix(sequences, distance_type, k_kmer);
        print_distance_stats(D, n);
    }

    // ========================================================================
    // ÉTAPE 3: Filtrage du graphe (optionnel)
    // ========================================================================
    
    int* adjacence = nullptr;
    
    if (epsilon < 1000) {
        if (rank == 0) {
            adjacence = filter_graph(D, n, epsilon);
        }
    } else {
        // Pas de filtrage, on utilise D directement
        if (rank == 0) {
            adjacence = new int[n * n];
            for (int i = 0; i < n * n; ++i) {
                adjacence[i] = D[i];
            }
            cout << "Pas de filtrage (epsilon=INF), graphe complet utilisé" << endl;
        }
    }

    // ========================================================================
    // ÉTAPE 4: Floyd-Warshall 
    // ========================================================================
    
    int* D_all = nullptr;
    
    if (use_floyd && rank == 0) {
        cout << "\nCalcul de Floyd-Warshall..." << endl;
        D_all = MatDistance(n, adjacence);
        cout << "Matrice des plus courts chemins calculée" << endl;
    } else if (rank == 0) {
        // Sans Floyd, on utilise la matrice de distances directe
        D_all = new int[n * n];
        for (int i = 0; i < n * n; ++i) {
            D_all[i] = adjacence[i];
        }
    }

    // ========================================================================
    // ÉTAPE 5: Clustering PAM
    // ========================================================================
    
    pam::Result result;
    
#ifdef USE_MPI
    // Version parallèle MPI
    vector<int> D_vec;
    if (rank == 0) {
        D_vec.assign(D_all, D_all + n*n);
    } else {
        D_vec.resize(n*n);
    }
    
    // Broadcast de la matrice à tous les processus
    MPI_Bcast(D_vec.data(), n*n, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) cout << "\nExécution de PAM (MPI, " << size << " processus)..." << endl;
    
    result = pam::pam_distributed(n, D_vec, k_clusters, seed, rank, size);
    
#else
    // Version séquentielle
    if (rank == 0) {
        cout << "\nExécution de PAM (séquentiel)..." << endl;
        vector<int> D_vec(D_all, D_all + n*n);
        result = pam::pam_sequential(n, D_vec, k_clusters, seed);
    }
#endif

    // ========================================================================
    // ÉTAPE 6: Affichage et export des résultats
    // ========================================================================
    
    if (rank == 0) {
        cout << "\n=== RÉSULTATS ===" << endl;
        cout << "Coût total: " << result.cost << endl;
        
        // Compter les membres par cluster
        vector<int> counts(k_clusters, 0);
        for (int m : result.membership) {
            if (m >= 0 && m < k_clusters) {
                counts[m]++;
            }
        }
        
        cout << "\nMédoïdes:" << endl;
        for (int m = 0; m < k_clusters; ++m) {
            cout << "  Cluster " << m << ": " << sequences[result.medoids[m]].id 
                 << " (" << counts[m] << " membres)" << endl;
        }
        
        cout << "\nDétails par cluster:" << endl;
        for (int m = 0; m < k_clusters; ++m) {
            cout << "\n--- Cluster " << m << " ---" << endl;
            cout << "Médoïde: " << sequences[result.medoids[m]].id << endl;
            cout << "Membres (" << counts[m] << "):" << endl;
            
            int shown = 0;
            for (size_t i = 0; i < result.membership.size(); ++i) {
                if (result.membership[i] == m) {
                    cout << "  - " << sequences[i].id;
                    if (i == (size_t)result.medoids[m]) {
                        cout << " [MÉDOÏDE]";
                    }
                    cout << endl;
                    
                    shown++;
                    if (shown >= 10 && counts[m] > 10) {
                        cout << "  ... (" << (counts[m] - 10) << " autres)" << endl;
                        break;
                    }
                }
            }
        }
        
        // Export complet
        export_clustering_results(sequences, result.medoids, 
                                 result.membership, result.cost, 
                                 output_file.c_str());
        
        cout << "\n✓ Terminé avec succès!" << endl;
        
        // Nettoyage
        delete[] D;
        delete[] adjacence;
        delete[] D_all;
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
    
    return 0;
}