# Projet MPI : Floyd-Warshall + PAM + Application ARN

## Vue d'ensemble

Ce projet implémente une pipeline complète pour classifier des séquences d'ARN :

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
    [Floyd-Warshall]
        ↓
Matrice de distances
        ↓
    [PAM Clustering]
        ↓
Partitions / Clusters
```

## Démarrage rapide

### Pipeline complet (recommandé)

```bash
./run_complete_pipeline.sh sequences.fasta 15 3 4
```

Arguments:
- `sequences.fasta` : Fichier FASTA avec les séquences
- `15` : Epsilon (seuil pour créer une arête)
- `3` : Nombre de clusters
- `4` : Nombre de processus MPI

### Étapes individuelles

#### 1. Traitement ARN
```bash
cd ARN
make
./arn_main sequences.fasta 15 graph.dot
```

#### 2. Floyd-Warshall
```bash
cd Floyd
make
mpirun -np 4 ./floyd_par ../graph.dot
```

#### 3. PAM
```bash
cd PAM
make
mpirun -np 4 ./pam ../graph.dot 3
```

## Structure du projet

```
Projet_MPI/
├── ARN/                        # Module de traitement des séquences d'ARN
│   ├── ARNSequence.hpp/cpp     # Distances et I/O FASTA
│   ├── main_arn.cpp            # Programme principal
│   ├── Makefile
│   ├── README.md
│   └── Doxyfile
│
├── Floyd/                      # Algorithme Floyd-Warshall parallèle
│   ├── ForGraph.hpp/cpp        # Lecture du graphe .dot
│   ├── FoydPar.hpp/cpp         # Floyd parallèle (MPI)
│   ├── Utils.hpp/cpp           # Utilitaires
│   ├── main.cpp                # Programme principal
│   ├── Makefile
│   ├── README.md
│   └── Exemples de graphes (.dot)
│
├── PAM/                        # Algorithme PAM parallèle
│   ├── PAM.hpp/cpp             # Implémentation PAM
│   ├── main_pam.cpp            # Programme principal
│   ├── Makefile
│   ├── README.md
│   └── Doxyfile
│
├── run_complete_pipeline.sh    # Script d'automatisation complète
└── README.md                   # Ce fichier
```

## Fonctionnalités principales

### ARN Module
- **Distance de Levenshtein** : Edit distance entre deux séquences
- **Distance de Hamming** : Distance pour séquences de même longueur
- **Lecture FASTA** : Support du format standard bioinformatique
- **Construction de graphe** : Graphe pondéré basé sur seuil epsilon

### Floyd-Warshall
- **Découpage par blocs** : Parallélisation efficace
- **MPI** : Distribution sur architecture mémoire distribuée
- **Support .dot** : Lecture/écriture au format Graphviz

### PAM
- **k-médoïdes** : Clustering basé sur les distances
- **Distribution répliquée** : Chaque processus détient la matrice complète
- **Version MPI** : Parallélisation des calculs d'échanges

## Documentation

Chaque module contient :
- Un **README.md** spécifique avec instructions de compilation et utilisation
- Un **Doxyfile** pour générer la documentation API

Générer la documentation :
```bash
cd ARN && doxygen Doxyfile
cd ../Floyd && doxygen Doxyfile
cd ../PAM && doxygen Doxyfile
```

## Dépendances

- **Compilateur** : g++ avec support C++11
- **MPI** : OpenMPI ou MPICH
- **Graphviz** : libcgraph pour la lecture des fichiers .dot
- **Doxygen** (optionnel) : Pour générer la documentation

## Installation des dépendances

### Ubuntu/Debian
```bash
sudo apt-get install build-essential openmpi-bin libopenmpi-dev libgraphviz-dev doxygen
```

### Fedora/RHEL
```bash
sudo dnf install gcc-c++ openmpi-devel graphviz-devel doxygen
```

### macOS
```bash
brew install open-mpi graphviz doxygen
```

## Exemple complet

### Données de test

1. **Créer un fichier FASTA** :
```bash
cat > sequences.fasta << 'EOF'
>seq1
ACGTACGTACGTACGT
>seq2
ACGTACGTACGTACGT
>seq3
TGCATGCATGCATGCA
EOF
```

2. **Exécuter le pipeline** :
```bash
./run_complete_pipeline.sh sequences.fasta 10 2 4
```

3. **Résultats** dans `results_YYYYMMDD_HHMMSS/` :
   - `arn_graph_temp.dot` : Graphe généré
   - `floyd_output.txt` : Matrice de distances
   - `pam_output.txt` : Partitions finales

## Nettoyage

```bash
cd ARN && make clean
cd ../Floyd && make clean
cd ../PAM && make clean
```

## Notes

- La complexité de Floyd-Warshall est O(n³)
- PAM est heuristique et peut ne pas trouver le clustering optimal
- Pour de très grands graphes, considérez des optimisations ou algorithmes approchés

## Auteurs

- Projet MPI - Université d'Orléans - M1 Informatique
- Année 2025-2026
