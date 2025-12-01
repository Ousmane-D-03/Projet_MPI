# Floyd–Warshall MPI (Bloc par Bloc)

Implémentation séquentielle et parallèle (MPI) de l’algorithme de Floyd–Warshall avec découpage en blocs.

## Compilation
make

## Exécution
mpirun -np <P> ./floyd graph.dot [iters]

- P doit être un carré parfait (4, 9, 16…)
- iters = nombre d'itérations (50 par défaut)
- Le nombre de nœuds du graphe doit être divisible par √P.

## Sorties
- Matrice finale
- Temps séquentiel / parallèle
- resultats_floyd.csv
