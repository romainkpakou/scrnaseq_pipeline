# Exemples de sorties

Ce repertoire contient des livrables produits par le pipeline, a titre
d'illustration. Ils permettent d'apprecier le rendu sans avoir a executer
l'analyse.

## report_pbmc3k.pdf

Rapport complet genere sur le jeu de donnees de demonstration PBMC 3k de
10X Genomics, avec la configuration `params/example_pbmc.yaml`.

Le document couvre l'integralite du flux de travail : controle qualite,
reduction de dimensionnalite, clustering, annotation cellulaire et analyse de
trajectoire. Il correspond aux metriques de validation annoncees dans le
README : 2638 cellules retenues, 9 types cellulaires annotes et un lignage
lymphocytaire T unique.

Une version HTML du meme rapport, avec table des matieres flottante, est
produite en parallele par le pipeline dans `results/report.html`.

## Reproduire ces sorties

    bash bin/download_test_data.sh

    nextflow run main.nf -profile conda \
        --config params/example_pbmc.yaml \
        --input_path data/raw/pbmc3k/hg19 \
        --output_path results
