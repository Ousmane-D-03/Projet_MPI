# Application ARN - Classification de Séquences d'ARN

## Vue d'ensemble

Ce module traite des séquences d'ARN et construit un graphe pondéré basé sur les distances entre séquences. Le graphe généré peut ensuite être utilisé avec :
- L'algorithme **Floyd-Warshall** pour calculer les plus courtes distances
- L'algorithme **PAM** pour partitionner les séquences en clusters

## Fonctionnalités

### Distances entre séquences

- **Distance de Levenshtein** (édition) : nombre minimum d'insertions, suppressions et substitutions
- **Distance de Hamming** : compte les positions différentes (séquences de même longueur)

### Construction du graphe

À partir d'un seuil epsilon (ε), deux séquences sont connectées si leur distance < ε.
Le poids de l'arête est la distance elle-même.

Ceci formalise l'énoncé :
```
∀vi, vj ∈ V  (vi, vj) ∈ E  ssi  d(vi, vj) < ε
∀e = (vi, vj) ∈ E  w(e) = d(vi, vj)
```

## Compilation

```bash
make
```

Fichier exécutable généré : `arn_main`

## Utilisation

### Syntaxe générale
```bash
./arn_main <fichier_fasta> <epsilon> [fichier_sortie.dot]
```

### Arguments
- `fichier_fasta` : Fichier contenant les séquences au format FASTA
- `epsilon` : Seuil de distance pour créer une arête (entier positif)
- `fichier_sortie.dot` : Fichier de sortie Graphviz (optionnel, défaut: `arn_graph.dot`)

### Exemple
```bash
./arn_main sequences.fasta 15 my_graph.dot
```

## Format d'entrée (FASTA)

Le fichier FASTA doit respecter le format suivant :
```
>label1
ACGTACGTACGT...
>label2
ACGTACGTACGT...
>label3
TGCATGCATGCA...
```

- Les lignes commençant par `>` contiennent les labels
- Les lignes suivantes contiennent la séquence (A, C, G, T)
- Les séquences peuvent s'étendre sur plusieurs lignes

## Exemple de test

```bash
make test
```

Cela crée un fichier `test_sequences.fasta` et exécute le programme.

## Fichier de sortie

Le fichier `.dot` généré peut être visualisé avec :
```bash
dot -Tpng arn_graph.dot -o arn_graph.png
```

Ou converti en PDF :
```bash
dot -Tpdf arn_graph.dot -o arn_graph.pdf
```

## Intégration avec Floyd et PAM

Après génération du fichier `.dot`, vous pouvez :

1. **Utiliser Floyd-Warshall**:
   ```bash
   cd ../Floyd
   mpirun -np 4 ./floyd_par ../ARN/arn_graph.dot
   ```

2. **Utiliser PAM pour le clustering**:
   ```bash
   cd ../PAM
   mpirun -np 4 ./pam ../ARN/arn_graph.dot 3
   ```

## Implémentation

### ARNSequence.hpp / ARNSequence.cpp
Fonctions principales :
- `levenshteinDistance()` : Distance d'édition
- `hammingDistance()` : Distance de Hamming
- `readFASTAFile()` : Lecture du fichier FASTA
- `computeDistanceMatrix()` : Calcul de toutes les distances
- `writeGraphDOT()` : Écriture du graphe au format DOT

### main_arn.cpp
Programme principal orchestrant le traitement complet.

## Notes

- Les distances de Levenshtein sont utilisées par défaut (efficaces pour des séquences de longueurs variables)
- La matrice de distances est stockée en mémoire avec une organisation row-major
- La complexité temporelle est O(n² × m²) où n est le nombre de séquences et m leur longueur moyenne
- Pour de très nombreuses séquences, considérez des optimisations ou une parallélisation

## Nettoyage

```bash
make clean
```

Supprime les fichiers compilés et les fichiers de test.
