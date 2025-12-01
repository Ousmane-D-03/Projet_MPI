Floyd–Warshall Parallèle en MPI (Bloc par Bloc)

Ce projet implémente l’algorithme de Floyd–Warshall en version parallèle distribuée à l’aide de MPI. La matrice de distances est découpée en blocs répartis sur une grille de processus, permettant l’exécution parallèle de l’algorithme.
Le programme inclut également une version séquentielle afin de comparer les performances.

Fonctionnalités

Lecture d’un graphe au format DOT.

Construction de la matrice des distances initiale.

Exécution de l’algorithme de Floyd–Warshall :

version séquentielle,

version parallèle par blocs avec MPI.

Diffusion des blocs pivot, ligne et colonne.

Reconstruction de la matrice finale avec MPI_Gather.

Génération d’un fichier CSV contenant les temps d’exécution moyens.

Pré-requis

OpenMPI, MPICH ou toute implémentation compatible MPI.

Compilateur C++ (g++ ou autre compatible C++11).

Un fichier ForGraph.hpp fournissant les fonctions :

lectureGraphe()

InitDk()

MatDistance()

La constante INF

Structure du projet
main.cpp
ForGraph.hpp
graph.dot
Makefile
resultats_floyd.csv (généré automatiquement)
README.md

Compilation

Grâce au Makefile :

make


Cela génère un exécutable nommé floyd (ou selon ton Makefile).

Exécution
Syntaxe générale
mpirun -np <nb_processus> ./floyd fichier.dot [nb_iterations]


<nb_processus> doit être un carré parfait (4, 9, 16, 25…)

fichier.dot est un graphe orienté pondéré

nb_iterations est optionnel (50 par défaut) et sert au benchmark

Exemples

Exécution avec 4 processus :

mpirun -np 4 ./floyd graph.dot


Exécution avec 16 processus et 100 itérations :

mpirun -np 16 ./floyd graph.dot 100

Contraintes importantes

Le nombre total de processus doit être un carré parfait :
1, 4, 9, 16, 25, 36, …

Le nombre de nœuds du graphe doit être divisible par √P.

Exemple :
Pour 16 processus → √16 = 4 → le nombre de nœuds doit être divisible par 4.

Sorties du programme
Affichage console

Matrice de distances initiale.

Matrice finale après exécution du Floyd–Warshall distribué.

Temps séquentiel moyen sur N itérations.

Temps parallèle moyen sur N itérations.

Fichier CSV généré

resultats_floyd.csv, format :

nb_processus,temps_seq,temps_par


Ce fichier peut être utilisé pour tracer des graphes de performance.