# Projet MPI : Floyd-Warshall + PAM + Application ARN

## Vue d'ensemble

Ce projet implémente une pipeline complète pour classifier des séquences d'ARN en utilisant des algorithmes parallèles (MPI) :

1. **Traitement des séquences d'ARN** (`ARN/`) : Lecture de fichiers FASTA et calcul des distances entre séquences
2. **Floyd-Warshall parallèle** (`Floyd/`) : Calcul des plus courtes distances dans un graphe pondéré
3. **PAM parallèle** (`PAM/`) : Clustering des nœuds en k partitions

## Architecture générale

```
Séquences d'ARN (FASTA)
        ↓
    [ARN Module]
        ↓
Graphe pondéré (.dot)
        ↓
    [Floyd-Warshall MPI]
        ↓
Matrice de distances
        ↓
    [PAM Clustering MPI]
        ↓
Partitions / Clusters
```

## Démarrage rapide

### Pipeline complet (recommandé)

```bash
cd Projet_MPI
./run_complete_pipeline.sh ../../../Data/petit.txt 15 3 4
```

**Arguments:**
- `../../../Data/petit.txt` : Fichier FASTA avec les séquences
- `15` : Epsilon (seuil pour créer une arête)
- `3` : Nombre de clusters pour PAM
- `4` : Nombre de processus MPI (doit être un carré parfait : 1, 4, 9, 16...)

### Étapes individuelles

#### 1. Traitement ARN
```bash
cd ARN
make
./arn_main ../../../Data/petit.txt 15 arn_graph.dot
```

#### 2. Floyd-Warshall
```bash
cd ../Floyd
make
mpirun -np 4 ./mpi_floyd ../ARN/arn_graph.dot
```

#### 3. PAM Clustering
```bash
cd ../PAM
make mpipam
mpirun -np 4 ./pam_mpi ../ARN/arn_graph.dot 3
```

## Structure du projet

```
Projet_MPI/
├── ARN/
│   ├── ARNSequence.hpp        # Structures et déclarations
│   ├── ARNSequence.cpp        # Implémentation (Levenshtein, Hamming, I/O FASTA)
│   ├── main_arn.cpp           # Programme principal
│   ├── Makefile               # Compilation
│   ├── README.md              # Documentation ARN
│   └── Doxyfile               # Configuration Doxygen
│
├── Floyd/
│   ├── ForGraph.hpp/cpp       # Lecture du graphe .dot (Graphviz)
│   ├── FoydPar.hpp/cpp        # Floyd parallèle par blocs (MPI)
│   ├── Utils.hpp/cpp          # Utilitaires d'affichage
│   ├── main.cpp               # Programme principal
│   ├── Makefile               # Compilation
│   ├── README.md              # Documentation Floyd
│   ├── Exemple3.dot           # Graphe de test
│   └── Doxyfile               # Configuration Doxygen
│
├── PAM/
│   ├── PAM.hpp/cpp            # Implémentation PAM (séquentiel + distribué)
│   ├── main_pam.cpp           # Programme principal
│   ├── Makefile               # Compilation
│   ├── README.md              # Documentation PAM
│   └── Doxyfile               # Configuration Doxygen
│
├── run_complete_pipeline.sh   # Script d'automatisation complète
├── README.md                  # Ce fichier
└── README_FULL.md             # Documentation technique détaillée
```

## Dépendances

### Obligatoires
- **Compilateur** : g++/mpic++ avec support C++11 ou supérieur
- **MPI** : OpenMPI ou MPICH
- **Make** : Système de compilation

### Optionnelles (mais recommandées)
- **Graphviz** : `libcgraph-dev` pour la lecture des fichiers .dot
- **Doxygen** : Pour générer la documentation API

### Installation (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install build-essential openmpi-bin libopenmpi-dev libgraphviz-dev doxygen graphviz
```

### Installation (Fedora/RHEL)
```bash
sudo dnf install gcc-c++ openmpi-devel graphviz-devel doxygen graphviz
```

### Installation (macOS)
```bash
brew install open-mpi graphviz doxygen
```

## Compilation

### Compiler tous les modules
```bash
make -C ARN
make -C Floyd
make -C PAM mpipam
```

### Nettoyer tous les modules
```bash
make -C ARN clean
make -C Floyd clean
make -C PAM clean
```

## Fonctionnalités principales

### ARN Module
- **Distance de Levenshtein** : Edit distance entre deux séquences (insertions, suppressions, substitutions)
- **Distance de Hamming** : Distance pour séquences de même longueur
- **Lecture FASTA** : Support du format standard bioinformatique
- **Construction de graphe** : Graphe pondéré basé sur seuil epsilon

### Floyd-Warshall
- **Découpage par blocs** : Parallélisation efficace (2D blocking)
- **MPI** : Distribution sur architecture mémoire distribuée
- **Support .dot** : Lecture/écriture au format Graphviz
- **Synchronisation** : Barrières MPI et communicateurs par ligne/colonne

### PAM (Partitioning Around Medoids)
- **k-médoïdes** : Clustering basé sur les distances
- **Version séquentielle** : Pour petits ensembles de données
- **Version MPI** : Parallélisation avec distribution par lignes
- **Optimisation** : Évaluation efficace des échanges (swaps)

## Exemple complet d'exécution

### 1. Préparer les données
```bash
cat > test_sequences.fasta << 'EOF'
>seq1
ACGTACGTACGTACGT
>seq2
ACGTACGTACGTACGT
>seq3
TGCATGCATGCATGCA
>seq4
TGCATGCATGCATGCA
EOF
```

### 2. Exécuter le pipeline
```bash
./run_complete_pipeline.sh test_sequences.fasta 10 2 4
```

### 3. Résultats dans `results_YYYYMMDD_HHMMSS/`
- `arn_graph_temp.dot` : Graphe généré
- `floyd_output.txt` : Matrice de distances
- `pam_output.txt` : Partitions finales

## Visualisation des résultats

### Graphe au format Graphviz
```bash
# PNG
dot -Tpng arn_graph.dot -o arn_graph.png
eog arn_graph.png  # affichage

# PDF
dot -Tpdf arn_graph.dot -o arn_graph.pdf
```

## Documentation technique

Pour la documentation détaillée :
- Consulter [README_FULL.md](README_FULL.md)
- Chaque module (`ARN/`, `Floyd/`, `PAM/`) contient son propre README

## Génération de la documentation API

```bash
cd ARN && doxygen Doxyfile && open html/index.html
cd ../Floyd && doxygen Doxyfile && open html/index.html
cd ../PAM && doxygen Doxyfile && open html/index.html
```

## Performances

### Complexité
- **Floyd-Warshall** : O(n³) séquentiel, O(n³/p) parallèle avec p processus
- **PAM** : O(n² × iterations) avec n = nombre de points
- **Levenshtein** : O(m₁ × m₂) avec m₁, m₂ = longueurs des séquences

### Notes
- Le nombre de processus MPI doit être un **carré parfait** pour Floyd (1, 4, 9, 16, 25...)
- Le nombre de nœuds du graphe doit être divisible par √(nombre de processus)
- Pour optimal performance, respecter : `nb_sequences = k × √nprocs` où k est entier

## Limitations connues

- Floyd par blocs requiert une distribution équitable des blocs
- PAM utilise une stratégie d'échange complète (brute-force) ; peut être lent pour n > 5000
- La version répliquée de PAM (chaque processus a D entière) limite à des matrices modérées

## Troubleshooting

### Erreur : "Erreur : nprocs pas carré parfait"
**Solution** : Utiliser un nombre de processus carré parfait (1, 4, 9, 16, 25...)

### Erreur : "Failed to open distance file"
**Solution** : Vérifier le chemin du fichier FASTA, utiliser des chemins absolus

### Graphviz non trouvé
**Solution** : `sudo apt-get install libgraphviz-dev`

### MPI non installé
**Solution** : `sudo apt-get install openmpi-bin libopenmpi-dev`

## Auteurs et remerciements

- Projet MPI - Université d'Orléans - M1 Informatique
- Année 2024-2025

## Licence

Ce projet est fourni à titre éducatif.

## Références

- Floyd-Warshall Algorithm : https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
- PAM Clustering : https://en.wikipedia.org/wiki/K-medoids
- Levenshtein Distance : https://en.wikipedia.org/wiki/Levenshtein_distance
- OpenMPI Documentation : https://www.open-mpi.org/doc/
- Graphviz : https://graphviz.org/