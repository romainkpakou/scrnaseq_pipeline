# Changelog

Toutes les evolutions notables de ce projet sont documentees dans ce fichier.
Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.0.0/) et le
versionnage respecte [SemVer](https://semver.org/lang/fr/).

## [0.2.0] : 2026-07-23

### Ajoute
- Module 5 : annotation cellulaire automatisee par score de module (AddModuleScore)
  confronte a une base de marqueurs canoniques externalisee
- Module 6 : reconstruction de trajectoire par Slingshot et identification des
  genes dynamiques par correlation de Spearman au pseudotemps
- Module 7 : generation d'un rapport HTML parametrable, dont l'ensemble des
  chiffres et sections sont derives dynamiquement des sorties du pipeline
- Base de marqueurs canoniques PBMC humain (9 types cellulaires)
- Fichier `nextflow.config` : manifeste, profils standard, conda et slurm,
  rapports de trace, execution et chronologie
- Environnement Conda reproductible (`environment/environment.yml`)
- Documentation exhaustive de l'ensemble des scripts R et modules Nextflow

### Modifie
- Les scripts R sont desormais declares comme entrees des process. Toute
  modification de code invalide le cache du module concerne et de ses
  dependants, ce qui rend `-resume` fiable durant le developpement.

### Valide
Reproduction a l'identique d'une analyse Seurat et Slingshot menee manuellement
sur le jeu de donnees PBMC 3k : 2638 cellules retenues, 9 clusters aux effectifs
identiques, 3446 marqueurs differentiels, 9 types cellulaires annotes avec les
memes proportions, lignage T unique et pseudotemps median identique.

## [0.1.0] : 2026-07-22

### Ajoute
- Structure initiale du projet et configuration Git
- Module 1 : controle qualite et normalisation
- Module 2 : analyse en composantes principales et selection des dimensions
- Module 3 : clustering non supervise et projection UMAP
- Module 4 : identification des marqueurs differentiels
- Configuration d'analyse centralisee au format YAML
- Licence MIT
