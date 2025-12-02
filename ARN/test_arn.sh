#!/bin/bash

# Script de test complet pour l'application ARN

echo "=========================================="
echo "Test de l'Application ARN Clustering"
echo "=========================================="
echo ""

# Couleurs pour la sortie
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Fonction de test
run_test() {
    local test_name=$1
    local command=$2
    
    echo -e "${YELLOW}[TEST]${NC} $test_name"
    echo "Commande: $command"
    
    if eval $command; then
        echo -e "${GREEN}✓ SUCCÈS${NC}"
        echo ""
        return 0
    else
        echo -e "${RED}✗ ÉCHEC${NC}"
        echo ""
        return 1
    fi
}

# Compteurs
total_tests=0
passed_tests=0

# Vérifier que le programme existe
if [ ! -f "./arn_seq" ]; then
    echo -e "${RED}Erreur: arn_seq n'existe pas. Exécutez 'make' d'abord.${NC}"
    exit 1
fi

# ============================================================================
# Test 1: Génération de séquences
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Génération de 30 séquences (3 familles)" \
    "./arn_seq --generate 30 100 test_30_3fam.fasta 3"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 2: Clustering basique avec distance d'édition
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering basique (edit, k=3)" \
    "./arn_seq --input test_30_3fam.fasta --distance edit --clusters 3 --output results_test2.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 3: Clustering avec filtrage du graphe
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering avec epsilon=20" \
    "./arn_seq --input test_30_3fam.fasta --distance edit --epsilon 20 --clusters 3 --output results_test3.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 4: Distance de Hamming (séquences même longueur)
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering avec distance Hamming" \
    "./arn_seq --input test_30_3fam.fasta --distance hamming --clusters 3 --output results_test4.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 5: Distance k-mer
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering avec distance k-mer (k=4)" \
    "./arn_seq --input test_30_3fam.fasta --distance kmer --kmer 4 --clusters 3 --output results_test5.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 6: Sans Floyd-Warshall
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering sans Floyd-Warshall" \
    "./arn_seq --input test_30_3fam.fasta --distance edit --clusters 3 --no-floyd --output results_test6.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 7: Plus grand jeu de données
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Génération de 60 séquences (4 familles)" \
    "./arn_seq --generate 60 80 test_60_4fam.fasta 4"; then
    passed_tests=$((passed_tests + 1))
fi

total_tests=$((total_tests + 1))
if run_test "Clustering sur 60 séquences (k=4)" \
    "./arn_seq --input test_60_4fam.fasta --distance edit --epsilon 25 --clusters 4 --output results_test7.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test 8: Différentes graines aléatoires
# ============================================================================
total_tests=$((total_tests + 1))
if run_test "Clustering avec seed=42" \
    "./arn_seq --input test_30_3fam.fasta --distance edit --clusters 3 --seed 42 --output results_seed42.txt"; then
    passed_tests=$((passed_tests + 1))
fi

total_tests=$((total_tests + 1))
if run_test "Clustering avec seed=99" \
    "./arn_seq --input test_30_3fam.fasta --distance edit --clusters 3 --seed 99 --output results_seed99.txt"; then
    passed_tests=$((passed_tests + 1))
fi

# ============================================================================
# Test MPI (si disponible)
# ============================================================================
if command -v mpirun &> /dev/null; then
    echo ""
    echo "=========================================="
    echo "Tests MPI"
    echo "=========================================="
    echo ""
    
    if [ -f "./arn_mpi" ]; then
        total_tests=$((total_tests + 1))
        if run_test "Clustering MPI avec 2 processus" \
            "mpirun -np 2 ./arn_mpi --input test_30_3fam.fasta --distance edit --clusters 3 --output results_mpi2.txt"; then
            passed_tests=$((passed_tests + 1))
        fi
        
        total_tests=$((total_tests + 1))
        if run_test "Clustering MPI avec 4 processus" \
            "mpirun -np 4 ./arn_mpi --input test_30_3fam.fasta --distance edit --epsilon 20 --clusters 3 --output results_mpi4.txt"; then
            passed_tests=$((passed_tests + 1))
        fi
    else
        echo -e "${YELLOW}[INFO]${NC} arn_mpi non compilé, tests MPI ignorés"
    fi
else
    echo -e "${YELLOW}[INFO]${NC} MPI non disponible, tests MPI ignorés"
fi

# ============================================================================
# Résumé des tests
# ============================================================================
echo ""
echo "=========================================="
echo "RÉSUMÉ DES TESTS"
echo "=========================================="
echo -e "Tests réussis: ${GREEN}$passed_tests${NC}/$total_tests"

if [ $passed_tests -eq $total_tests ]; then
    echo -e "${GREEN}✓ Tous les tests ont réussi!${NC}"
    exit_code=0
else
    echo -e "${RED}✗ Certains tests ont échoué${NC}"
    exit_code=1
fi

# ============================================================================
# Vérification des fichiers de sortie
# ============================================================================
echo ""
echo "Fichiers générés:"
ls -lh test_*.fasta results_*.txt 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "Pour voir les résultats détaillés:"
echo "  cat results_test2.txt"
echo "  cat clustering_results.txt"

exit $exit_code