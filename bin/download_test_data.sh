#!/usr/bin/env bash

# ==============================================================================
# TELECHARGEMENT DU JEU DE DONNEES DE DEMONSTRATION : PBMC 3k
# ==============================================================================
#
# ROLE
#   Recupere le jeu de donnees PBMC 3k distribue publiquement par 10X Genomics
#   et le place a l'emplacement attendu par la configuration d'exemple.
#
#   Ce jeu de donnees sert de reference au pipeline : il permet a tout
#   utilisateur de verifier que son installation fonctionne avant d'engager ses
#   propres donnees, et de reproduire les metriques de validation annoncees dans
#   le README.
#
# USAGE
#   bash bin/download_test_data.sh
#
#   Puis, depuis la racine du projet :
#     nextflow run main.nf \
#         --config params/example_pbmc.yaml \
#         --input_path data/raw/pbmc3k/hg19 \
#         --output_path results
#
# VOLUME
#   Archive compressee : environ 7 Mo. Donnees decompressees : environ 25 Mo.
#
# DEPENDANCES
#   curl ou wget, tar
# ==============================================================================

# Arret immediat en cas d'erreur, de variable non definie, ou d'echec dans un
# enchainement de commandes par tube. Sans ces options, un telechargement
# echoue passerait inapercu et le script poursuivrait sur une archive absente.
set -euo pipefail

URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
DEST_DIR="data/raw/pbmc3k"
ARCHIVE="pbmc3k_filtered_gene_bc_matrices.tar.gz"

echo "=== Telechargement du jeu de donnees de demonstration PBMC 3k ==="

# Verification du repertoire d'execution : le script utilise des chemins
# relatifs et doit donc etre lance depuis la racine du projet.
if [ ! -f "main.nf" ]; then
    echo "ERREUR : ce script doit etre lance depuis la racine du projet." >&2
    echo "         Exemple : bash bin/download_test_data.sh" >&2
    exit 1
fi

mkdir -p "${DEST_DIR}"

# Court-circuit si les donnees sont deja presentes, afin d'eviter un
# telechargement inutile lors d'une reexecution du script.
if [ -f "${DEST_DIR}/hg19/matrix.mtx" ]; then
    echo "Les donnees sont deja presentes dans ${DEST_DIR}/hg19"
    echo "Supprimer ce repertoire pour forcer un nouveau telechargement."
    exit 0
fi

# Selection de l'outil de telechargement disponible sur le systeme.
echo "Telechargement depuis 10X Genomics..."
if command -v curl > /dev/null 2>&1; then
    curl -L -o "${DEST_DIR}/${ARCHIVE}" "${URL}"
elif command -v wget > /dev/null 2>&1; then
    wget -O "${DEST_DIR}/${ARCHIVE}" "${URL}"
else
    echo "ERREUR : ni curl ni wget n'est disponible sur ce systeme." >&2
    exit 1
fi

echo "Decompression..."
# L'option --strip-components retire le niveau filtered_gene_bc_matrices de
# l'arborescence de l'archive, de sorte que le repertoire hg19 se retrouve
# directement sous data/raw/pbmc3k, chemin attendu par la configuration.
tar -xzf "${DEST_DIR}/${ARCHIVE}" -C "${DEST_DIR}" --strip-components=1

rm -f "${DEST_DIR}/${ARCHIVE}"

# Verification que les trois fichiers attendus sont bien presents.
echo ""
echo "=== Verification ==="
for f in matrix.mtx barcodes.tsv genes.tsv; do
    if [ -f "${DEST_DIR}/hg19/${f}" ]; then
        echo "  OK   ${f}"
    else
        echo "  MANQUANT   ${f}" >&2
    fi
done

echo ""
echo "Donnees disponibles dans : ${DEST_DIR}/hg19"
echo ""
echo "Lancer l'analyse de demonstration :"
echo "  nextflow run main.nf \\"
echo "      --config params/example_pbmc.yaml \\"
echo "      --input_path ${DEST_DIR}/hg19 \\"
echo "      --output_path results"
