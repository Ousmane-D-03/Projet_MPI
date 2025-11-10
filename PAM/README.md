# PAM (k-medoids) — Séquentiel + MPI (répliqué)

Ce dossier contient une implémentation de l'algorithme PAM (Partitioning Around Medoids).

But : fournir une version simple répliquée (chaque processus détient la matrice D entière) utile
pour n modéré (ex : 500 / 2000) et qui s'intègre avec le code `Floyd/` pour générer la matrice
de distances à partir d'un fichier `.dot`.

Build
```
cd PAM
make
```

Run
```
# séquentiel (mpirun 1 proc) :
mpirun -np 1 ./pam path/to/graph.dot k [seed]

# parallèle (ex : 4 processus)
mpirun -np 4 ./pam path/to/graph.dot k [seed]
```

Notes
- Le code utilise la lecture `.dot` et la fonction `MatDistance` de `Floyd/` (Graphviz requis).
- Mode de distribution : répliqué pour D, partition des lignes pour le calcul des deltas.
- Cette version est brute-force (évalue tous les échanges) et convient pour des n modestes.
