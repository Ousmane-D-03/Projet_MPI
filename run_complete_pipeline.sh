#!/bin/bash

##
# @file run_complete_pipeline.sh
# @brief Script d'intégration complète : ARN → Floyd-Warshall → PAM
# @author Projet MPI
# @date 2025
#
# Ce script automatise l'ensemble du pipeline :
# 1. Traitement des séquences d'ARN
# 2. Construction du graphe
# 3. Calcul des plus courtes distances (Floyd-Warshall)
# 4. Clustering avec PAM
#

set -e  # Arrêter en cas d'erreur

# Couleurs pour l'affichage
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # Pas de couleur

# Fonctions utilitaires
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERREUR]${NC} $1"
}

print_step() {
    echo -e "\n${YELLOW}===== $1 =====${NC}\n"
}

# Vérifier les arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <fichier_fasta> <epsilon> [nb_clusters] [nb_processus]"
    echo ""
    echo "Arguments:"
    echo "  fichier_fasta   Fichier FASTA contenant les séquences d'ARN"
    echo "  epsilon         Seuil de distance pour le graphe"
    echo "  nb_clusters     Nombre de clusters pour PAM (optionnel, défaut: 3)"
    echo "  nb_processus    Nombre de processus MPI (optionnel, défaut: 1)"
    echo ""
    echo "Exemple:"
    echo "  $0 sequences.fasta 15 3 4"
    exit 1
fi

FASTA_FILE="$1"
EPSILON="$2"
NB_CLUSTERS="${3:-3}"
NB_PROCS="${4:-1}"

# Vérifier que le fichier FASTA existe
if [ ! -f "$FASTA_FILE" ]; then
    print_error "Le fichier FASTA '$FASTA_FILE' n'existe pas"
    exit 1
fi

# Noms des fichiers intermédiaires
GRAPH_DOT="arn_graph_temp.dot"
DISTANCES_FILE="distances_matrix.txt"
RESULT_DIR="results_$(date +%Y%m%d_%H%M%S)"

# Créer le répertoire de résultats
mkdir -p "$RESULT_DIR"

print_step "Étape 1: Traitement des séquences d'ARN"

# Aller dans le répertoire ARN
cd ARN || exit 1

# Compiler si nécessaire
if [ ! -f "arn_main" ]; then
    print_info "Compilation du module ARN..."
    make clean > /dev/null 2>&1
    make > /dev/null 2>&1
fi

# Vérifier que le fichier FASTA est accessible
if [ ! -f "../$FASTA_FILE" ]; then
    print_error "Fichier FASTA introuvable: ../$FASTA_FILE"
    exit 1
fi

print_info "Traitement: $FASTA_FILE avec epsilon=$EPSILON"
./arn_main "../$FASTA_FILE" "$EPSILON" "$GRAPH_DOT"

# Copier le graphe vers le répertoire de résultats
cp "$GRAPH_DOT" "../$RESULT_DIR/"

print_info "Graphe généré: $GRAPH_DOT"

# Revenir au répertoire parent
cd .. || exit 1

print_step "Étape 2: Calcul des plus courtes distances (Floyd-Warshall)"

# Aller dans le répertoire Floyd
cd Floyd || exit 1

# Compiler si nécessaire
if [ ! -f "floyd_par" ] && [ ! -f "main" ]; then
    print_info "Compilation du module Floyd..."
    make clean > /dev/null 2>&1
    make > /dev/null 2>&1
fi

print_info "Exécution Floyd-Warshall avec $NB_PROCS processus"

# Chercher l'exécutable
FLOYD_EXEC=""
for exec in "mpi_floyd" "floyd_par" "main"; do
    if [ -f "$exec" ]; then
        FLOYD_EXEC="$exec"
        break
    fi
done

if [ -z "$FLOYD_EXEC" ]; then
    print_error "Aucun exécutable Floyd trouvé (cherché: mpi_floyd, floyd_par, main)"
    exit 1
fi

print_info "Utilisation de l'exécutable: $FLOYD_EXEC"

if [ "$NB_PROCS" -eq 1 ]; then
    ./$FLOYD_EXEC "../$RESULT_DIR/$GRAPH_DOT" > "../$RESULT_DIR/floyd_output.txt" 2>&1
else
    mpirun -np "$NB_PROCS" ./$FLOYD_EXEC "../$RESULT_DIR/$GRAPH_DOT" > "../$RESULT_DIR/floyd_output.txt" 2>&1
fi

print_info "Résultats Floyd sauvegardés"

# Revenir au répertoire parent
cd .. || exit 1

print_step "Étape 3: Clustering avec PAM"

# Aller dans le répertoire PAM
cd PAM || exit 1

# Compiler si nécessaire
if [ ! -f "pam" ] && [ ! -f "pam_mpi" ]; then
    print_info "Compilation du module PAM..."
    make clean > /dev/null 2>&1
    make > /dev/null 2>&1
fi

print_info "Exécution PAM avec $NB_CLUSTERS clusters et $NB_PROCS processus"

# Chercher l'exécutable
PAM_EXEC=""
for exec in "pam_mpi" "pam"; do
    if [ -f "$exec" ]; then
        PAM_EXEC="$exec"
        break
    fi
done

if [ -z "$PAM_EXEC" ]; then
    print_error "Aucun exécutable PAM trouvé (cherché: pam_mpi, pam)"
    exit 1
fi

print_info "Utilisation de l'exécutable: $PAM_EXEC"

mpirun -np "$NB_PROCS" ./$PAM_EXEC "../$RESULT_DIR/$GRAPH_DOT" "$NB_CLUSTERS" > "../$RESULT_DIR/pam_output.txt" 2>&1

print_info "Résultats PAM sauvegardés"

# Revenir au répertoire parent
cd .. || exit 1

print_step "Pipeline terminé avec succès!"

echo -e "${GREEN}Résultats dans:${NC} $RESULT_DIR/"
echo ""
echo "Fichiers disponibles:"
ls -lah "$RESULT_DIR/"
echo ""
echo -e "${GREEN}Visualisation du graphe:${NC}"
echo "  dot -Tpng $RESULT_DIR/$GRAPH_DOT -o $RESULT_DIR/graph.png"
