/**
/**
 * @file main_arn.cpp
 * @brief Programme principal pour traiter des séquences d'ARN et construire un graphe
 * @author Projet MPI
 * @date 2025
 * 
 * Ce programme lit des séquences d'ARN au format FASTA, calcule les distances
 * entre les séquences, et construit un graphe pondéré qui peut être utilisé
 * avec l'algorithme de Floyd-Warshall.
 * 
 * Usage: ./arn_main <fichier_fasta> <epsilon> [output_dot]
 */

#include <iostream>
#include <string>
#include <vector>
#include "ARNSequence.hpp"

using namespace std;

void printUsage(const char* programName) {
    cout << "Usage: " << programName << " <fichier_fasta> <epsilon> [output_dot]" << endl;
    cout << endl;
    cout << "Arguments:" << endl;
    cout << "  fichier_fasta  Fichier contenant les séquences d'ARN (format FASTA)" << endl;
    cout << "  epsilon        Seuil de distance pour créer une arête" << endl;
    cout << "  output_dot     Fichier de sortie au format Graphviz .dot (optionnel)" << endl;
    cout << endl;
    cout << "Exemple:" << endl;
    cout << "  " << programName << " sequences.fasta 15 sequences_graph.dot" << endl;
}

int main(int argc, char* argv[]) {
    
    if (argc < 3 || argc > 4) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    string fastaFile = argv[1];
    int epsilon;
    string outputFile = "arn_graph.dot";
    
    // Parser epsilon
    try {
        epsilon = stoi(argv[2]);
    } catch (...) {
        cerr << "Erreur : epsilon doit être un entier positif" << endl;
        return EXIT_FAILURE;
    }
    
    if (epsilon <= 0) {
        cerr << "Erreur : epsilon doit être strictement positif" << endl;
        return EXIT_FAILURE;
    }
    
    // Fichier de sortie optionnel
    if (argc == 4) {
        outputFile = argv[3];
    }
    
    cout << "=== Traitement des séquences d'ARN ===" << endl;
    cout << "Fichier FASTA : " << fastaFile << endl;
    cout << "Epsilon : " << epsilon << endl;
    cout << "Fichier de sortie : " << outputFile << endl << endl;
    
    // Lire les séquences
    vector<ARNSeq> sequences;
    int nbSeq = readFASTAFile(fastaFile, sequences);
    
    if (nbSeq <= 0) {
        cerr << "Erreur : impossible de lire le fichier FASTA" << endl;
        return EXIT_FAILURE;
    }
    
    // Afficher les séquences
    cout << endl << "Séquences lues :" << endl;
    for (const auto& seq : sequences) {
        printARNSeq(seq);
    }
    
    // Calculer la matrice de distances
    cout << endl << "Calcul de la matrice de distances (Levenshtein)..." << endl;
    int* distanceMatrix = computeDistanceMatrix(sequences, levenshteinDistance);
    
    // Afficher la matrice de distances
    cout << endl;
    printDistanceMatrix(distanceMatrix, sequences.size());
    
    // Écrire le graphe au format DOT
    cout << endl;
    if (writeGraphDOT(sequences, distanceMatrix, epsilon, outputFile) != 0) {
        cerr << "Erreur : impossible d'écrire le fichier DOT" << endl;
        delete[] distanceMatrix;
        return EXIT_FAILURE;
    }
    
    // Nettoyage
    delete[] distanceMatrix;
    
    cout << endl << "Processus terminé avec succès !" << endl;
    cout << "Utilisez le fichier " << outputFile << " avec Floyd-Warshall ou PAM" << endl;
    
    return EXIT_SUCCESS;
}
*/
// main_mpi.cpp
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include "ARNSequence.hpp"
#include "../PAM/PAM.hpp"   // adapter selon ton chemin
#include "../Floyd/FoydPar.hpp"       // ton fichier avec floydBlocsMPI
using namespace std;

void printUsage(const char* programName) {
    cout << "Usage: " << programName << " <fichier_fasta> <epsilon> [output_dot]" << endl;
    cout << endl;
    cout << "Arguments:" << endl;
    cout << "  fichier_fasta  Fichier contenant les séquences d'ARN (format FASTA)" << endl;
    cout << "  epsilon        Seuil de distance pour créer une arête" << endl;
    cout << "  output_dot     Fichier de sortie au format Graphviz .dot (optionnel)" << endl;
    cout << endl;
    cout << "Exemple:" << endl;
    cout << "  " << programName << " sequences.fasta 15 sequences_graph.dot" << endl;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc < 3 || argc > 4) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    string fastaFile = argv[1];
    int epsilon = stoi(argv[2]);
    int k_clusters = stoi(argv[3]);
    string outputFile = (argc >= 5) ? argv[4] : "arn_graph.dot";

    vector<ARNSeq> sequences;
    int nbSeq = 0;
    int* distanceMatrix = nullptr;

    if(pid == 0) {
        nbSeq = readFASTAFile(fastaFile, sequences);
        if(nbSeq <= 0) {
            cerr << "Erreur : impossible de lire le fichier FASTA" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cout << "Nombre de séquences : " << nbSeq << endl;

        // calcul de la matrice de distances
        distanceMatrix = computeDistanceMatrix(sequences, levenshteinDistance);
    }

    // Vérifier que nbSeq est divisible par sqrt(nprocs)
    int p_sqrt = static_cast<int>(sqrt(nprocs));
    if(p_sqrt*p_sqrt != nprocs) {
        if(pid == 0) cerr << "Erreur : le nombre de processus doit être un carré parfait" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int block_size = (pid==0 ? nbSeq : 0)/p_sqrt; // taille d'un bloc par ligne/col
    MPI_Bcast(&nbSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* D_local = new int[block_size*block_size];

    // Découper et distribuer la matrice
    decouperMatrice(distanceMatrix, D_local, nbSeq, block_size, p_sqrt, 0, pid);

    // Appliquer Floyd–Warshall MPI
    int* D_global = floydBlocsMPI(D_local, nbSeq, p_sqrt, pid, 0);

    // Sur le root : appliquer PAM
    if(pid == 0) {
        cout << "Matrice des plus courts chemins calculée." << endl;

        // PAM
        pam::Result res = pam::pam_sequential(nbSeq, vector<int>(D_global, D_global + nbSeq*nbSeq), k_clusters, 42);
        cout << "Résultat PAM :" << endl;
        cout << "Médoïdes : ";
        for(auto m : res.medoids) cout << m << " ";
        cout << endl;

        // Écriture du graphe DOT
        if(writeGraphDOT(sequences, D_global, epsilon, outputFile) == 0) {
            cout << "Graphe écrit dans " << outputFile << endl;
        }

        delete[] distanceMatrix;
        delete[] D_global;
    }

    delete[] D_local;
    MPI_Finalize();
    return 0;
}
