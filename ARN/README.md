# Application Floyd + PAM aux Séquences d'ARN

## Description

Cette application combine les algorithmes Floyd-Warshall et PAM (k-médoïdes) pour effectuer un clustering de séquences d'ARN basé sur leurs distances.

## Pipeline Complet

```
Séquences FASTA → Matrice de Distances → Graphe Filtré → Floyd-Warshall → PAM → Clusters
```

## Prérequis

- **Compilateur C++** avec support C++11 (g++)
- **MPI** (OpenMPI ou MPICH) pour la version parallèle
- **GraphViz** (libcgraph) pour la lecture de graphes .dot
- Bibliothèques Floyd et PAM (dans les dossiers `../Floyd` et `../PAM`)


## Compilation

### Version Séquentielle
```bash
cd RNA
make rna_seq
```

### Version MPI Parallèle
```bash
make rna_mpi
```

### Tout compiler
```bash
make all
```

## Utilisation

### 1. Génération de Séquences Test

```bash
# Générer 30 séquences de 100 bases réparties en 3 familles
./rna_seq --generate 30 100 test_sequences.fasta 3

# Générer 100 séquences de 50 bases en 5 familles
./rna_seq --generate 100 50 large_test.fasta 5
```

### 2. Clustering de Séquences

#### Exemple Simple (distance d'édition)
```bash
./rna_seq --input test_sequences.fasta --clusters 3
```

#### Avec Tous les Paramètres
```bash
./rna_seq \
  --input test_sequences.fasta \
  --distance edit \
  --epsilon 20 \
  --clusters 3 \
  --seed 42 \
  --output results.txt
```

#### Sans Floyd-Warshall (distances directes)
```bash
./rna_seq \
  --input test_sequences.fasta \
  --distance hamming \
  --clusters 4 \
  --no-floyd
```

#### Avec Distance k-mer
```bash
./rna_seq \
  --input test_sequences.fasta \
  --distance kmer \
  --kmer 4 \
  --clusters 3
```

### 3. Version Parallèle MPI

```bash
# Avec 4 processus
mpirun -np 4 ./rna_mpi \
  --input test_sequences.fasta \
  --distance edit \
  --epsilon 20 \
  --clusters 3
```

## Options Détaillées

| Option | Description | Défaut |
|--------|-------------|--------|
| `--generate <n> <len> <output> [fam]` | Génère n séquences de longueur len | - |
| `--input <file>` | Fichier FASTA d'entrée | Requis |
| `--distance <type>` | Type: hamming, edit, kmer | edit |
| `--kmer <k>` | Taille des k-mers | 3 |
| `--epsilon <val>` | Seuil de filtrage du graphe | 1000 (INF) |
| `--clusters <k>` | Nombre de clusters PAM | 3 |
| `--seed <val>` | Graine aléatoire | 12345 |
| `--output <file>` | Fichier de résultats | clustering_results.txt |
| `--no-floyd` | Désactive Floyd-Warshall | Activé |
| `--help` | Affiche l'aide | - |

## Types de Distances

### 1. Distance de Hamming (`--distance hamming`)
- **Principe**: Nombre de positions différentes
- **Contrainte**: Séquences de même longueur
- **Usage**: Séquences alignées de même taille

**Exemple**:
```
s1 = ACGT
s2 = ACCT
d_H(s1, s2) = 1  (position 2: G ≠ C)
```

### 2. Distance d'Édition (`--distance edit`)
- **Principe**: Nombre d'insertions/suppressions/substitutions
- **Contrainte**: Aucune
- **Usage**: Séquences de longueurs variées

**Exemple**:
```
s1 = ACGT
s2 = ACCGT
d_edit(s1, s2) = 1  (insertion de C)
```

### 3. Distance k-mer (`--distance kmer`)
- **Principe**: Similarité des sous-chaînes de longueur k
- **Contrainte**: Longueur ≥ k
- **Usage**: Approximation rapide

**Exemple** (k=2):
```
s1 = ACGT → {AC, CG, GT}
s2 = ACCT → {AC, CC, CT}
Jaccard = |intersection|/|union| = 1/5 = 0.2
Distance = (1 - 0.2) * 100 = 80
```

## Paramètre Epsilon (ε)

Le paramètre epsilon contrôle la densité du graphe :

- **ε = 0**: Graphe vide (aucune arête)
- **ε petit** (ex: 5-10): Graphe sparse, seules les séquences très similaires sont connectées
- **ε moyen** (ex: 20-30): Bon compromis
- **ε = 1000 (INF)**: Graphe complet, toutes les séquences sont connectées

### Recommandations

| Type de Distance | Epsilon Suggéré |
|-----------------|----------------|
| Hamming | 5-15 |
| Edit | 10-30 |
| k-mer | 30-60 |

**Conseil**: Exécutez d'abord sans epsilon pour voir les statistiques de distances, puis choisissez un epsilon approprié.

## Interprétation des Résultats

### Sortie Console

```
=== RÉSULTATS ===
Coût total: 1245
Médoïdes:
  Cluster 0: seq0_family0 (10 membres)
  Cluster 1: seq1_family1 (9 membres)
  Cluster 2: seq2_family2 (11 membres)
```

### Fichier de Résultats (`clustering_results.txt`)

```
=== Résultats du Clustering PAM ===
Coût total: 1245
Nombre de clusters: 3

--- Cluster 0 ---
Médoïde: seq0_family0
Taille: 10 séquences
Membres:
  - seq0_family0
  - seq3_family0
  - seq6_family0
  ...
```

## Exemples Complets

### Exemple 1: Clustering Basique

```bash
# Générer des données
./rna_seq --generate 30 100 test.fasta 3

# Clustering avec distance d'édition
./rna_seq --input test.fasta --distance edit --clusters 3

# Résultats attendus: 3 clusters correspondant aux 3 familles
```

### Exemple 2: Graphe Filtré + Floyd

```bash
# Créer un graphe sparse avec epsilon=15
./rna_seq --input test.fasta --distance edit --epsilon 15 --clusters 3

# Floyd-Warshall calcule les plus courts chemins dans ce graphe
```

### Exemple 3: Comparaison de Distances

```bash
# Hamming
./rna_seq --input test.fasta --distance hamming --clusters 3 --output results_hamming.txt

# Édition
./rna_seq --input test.fasta --distance edit --clusters 3 --output results_edit.txt

# k-mer
./rna_seq --input test.fasta --distance kmer --kmer 4 --clusters 3 --output results_kmer.txt

# Comparer les coûts et la qualité
```

### Exemple 4: Grande Échelle avec MPI

```bash
# Générer un grand jeu de données
./rna_seq --generate 200 150 large.fasta 5

# Clustering parallèle avec 8 processus
mpirun -np 8 ./rna_mpi \
  --input large.fasta \
  --distance edit \
  --epsilon 30 \
  --clusters 5 \
  --output large_results.txt
```

## Format FASTA

Le programme accepte les fichiers FASTA standard :

```
>seq1_description
ACGTACGTACGT
ACGTACGTACGT
>seq2_description
GGTTAACCGGTT
>seq3_description
ACGTACGTACGT
```

**Règles**:
- Chaque séquence commence par `>` suivi d'un identifiant
- Les lignes suivantes contiennent la séquence (peuvent être multi-lignes)
- Alphabet: A, C, G, T (ou U pour ARN)

## Performances

### Complexités

| Étape | Complexité | Commentaire |
|-------|-----------|-------------|
| Calcul de D | O(n² × m²) | n=séquences, m=longueur (pour edit) |
| Floyd-Warshall | O(n³) | Peut être désactivé |
| PAM | O(k(n-k)² × iter) | k=clusters, iter=itérations |

### Temps d'Exécution Estimés (1 cœur)

| n (séquences) | m (longueur) | Distance | Temps |
|--------------|--------------|----------|-------|
| 30 | 100 | edit | ~5s |
| 100 | 100 | edit | ~1min |
| 200 | 150 | edit | ~10min |
| 30 | 100 | hamming | ~1s |
| 100 | 100 | kmer | ~10s |

**Note**: Floyd-Warshall ajoute ~O(n³) au temps total.

## Dépannage

### Erreur: "Failed to read graph"
- Vérifiez que le fichier FASTA existe
- Format FASTA valide requis

### Erreur: "Invalid k"
- k (nombre de clusters) doit être entre 1 et n (nombre de séquences)

### Distance INF avec Hamming
- Les séquences doivent avoir la même longueur
- Utilisez `--distance edit` pour des longueurs variables

### Graphe vide après filtrage
- Epsilon trop petit, aucune arête créée
- Augmentez epsilon ou vérifiez les statistiques de distances

### Segmentation fault
- Vérifiez que les bibliothèques Floyd et PAM sont compilées
- Assurez-vous que GraphViz est installé

## Tests Automatiques

```bash
# Test complet séquentiel
make test_seq

# Test complet MPI
make test_mpi
```

## Validation Biologique

Pour valider la qualité du clustering :

1. **Annotations**: Si disponibles, comparer avec des familles connues
2. **Silhouette Score**: Mesure de cohésion intra-cluster vs séparation inter-cluster
3. **Visualisation**: Export pour analyse visuelle (à implémenter)

## Références

- **Floyd-Warshall**: Algorithme des plus courts chemins
- **PAM**: Kaufman & Rousseeuw (1990) "Finding Groups in Data"
- **Distance d'édition**: Levenshtein (1966)
- **k-mer**: Approche standard en bioinformatique

## Auteurs et Support

Projet M1 Info - Université d'Orléans 2025-2026

Pour toute question ou bug, consultez la documentation du projet.

## Licence

Code fourni dans le cadre du projet académique.