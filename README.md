# scrnaseq_pipeline

Pipeline **Nextflow DSL2** reproductible et modulaire pour l'analyse complète de données single-cell RNA-seq (10X Genomics), depuis les matrices brutes jusqu'aux trajectoires cellulaires.

## Aperçu

Le pipeline enchaîne six modules paramétrables via un fichier YAML :

| Module | Étape | Sorties |
|--------|-------|---------|
| 1 | Contrôle qualité + normalisation | Objet filtré, figures QC, HVG |
| 2 | ACP + sélection des composantes | ElbowPlot, JackStraw (optionnel) |
| 3 | Clustering (Louvain) + UMAP | Clusters, projection UMAP, marqueurs canoniques |
| 4 | Marqueurs différentiels | FindAllMarkers, DotPlot, Heatmap |
| 5 | Annotation cellulaire automatisée | Types cellulaires (AddModuleScore), UMAP annoté |
| 6 | Trajectoires (Slingshot) | Pseudotemps, gènes dynamiques |

## Validation

Le pipeline a été validé sur le dataset **PBMC 3k** (10X Genomics) : chaque étape reproduit à l'identique une analyse Seurat/Slingshot menée manuellement (9 types cellulaires, proportions et signatures moléculaires exactes, lignage T CD4 naïf vers CD8 effecteur).

## Utilisation

```bash
nextflow run main.nf \
    --config params/example_pbmc.yaml \
    --input_path data/raw/pbmc3k/hg19 \
    --output_path results
```

## Prérequis

- Nextflow >= 22.04
- R >= 4.3 avec Seurat 5.x, Slingshot 2.x
- Java 17

## Architecture

Chaque module est composé d'un process Nextflow (`modules/`) appelant un script R paramétrable (`bin/`). La configuration d'analyse est centralisée dans un fichier YAML (`params/`), la base de marqueurs canoniques est externalisée (`assets/`).

## Auteur

Romain KPAKOU : M2 Bioinformatique, Biostatistique et Biologie Computationnelle.

## Licence

MIT
