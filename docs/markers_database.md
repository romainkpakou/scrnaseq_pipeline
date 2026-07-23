# Base de marqueurs canoniques

Le module d'annotation attribue une identite biologique a chaque cluster en
confrontant son profil d'expression a une base de marqueurs. Cette base est un
fichier YAML externe au pipeline, ce qui permet de l'adapter a n'importe quel
tissu ou espece sans modifier le code.

## Format

Chaque entree de premier niveau est un nom de type cellulaire. Ce nom sera
utilise tel quel pour etiqueter les clusters dans toutes les sorties.

    Nom du type cellulaire:
      markers: [GENE1, GENE2, GENE3, GENE4, GENE5]
      description: "Description libre, non exploitee par le code"

Exemple issu de la base PBMC humain fournie :

    CD4 T naive:
      markers: [CCR7, LEF1, TCF7, IL7R, MAL]
      description: "Lymphocytes T CD4+ naifs, quiescents, circulants"

    NK:
      markers: [GNLY, GZMB, NKG7, SPON2, AKR1C3]
      description: "Cellules Natural Killer"

Le champ description est purement documentaire : il aide le lecteur du fichier
sans intervenir dans le calcul.

## Regles de redaction

**Nomenclature des genes.** Les symboles doivent correspondre exactement a ceux
de la matrice d'expression. La convention differe selon l'espece : majuscules
chez l'humain (CD14, GNLY), initiale majuscule seule chez la souris (Cd14,
Gnly). Un marqueur mal orthographie est silencieusement ignore.

**Nombre de marqueurs.** Trois a huit marqueurs par type constituent un bon
compromis. En dessous de trois, le score devient instable et sensible au bruit.
Au dela de huit, l'ajout de marqueurs peu specifiques dilue le signal et
rapproche les scores de types distincts.

**Specificite.** Privilegier des marqueurs exprimes preferentiellement par le
type vise. Un marqueur partage entre plusieurs types de la base n'est pas
disqualifiant, mais reduit le pouvoir discriminant entre ces types. C'est
notamment le cas d'IL7R, present a la fois chez les lymphocytes T naifs et
memoires, deux populations formant un continuum biologique.

**Ordre des marqueurs.** Le premier marqueur de chaque liste est utilise comme
marqueur representatif dans la figure de synthese produite en fin d'annotation.
Placer en tete le marqueur le plus emblematique du type.

## Verification apres annotation

Les scores obtenus sont consignes dans `results/logs/annotation_log.txt`. Leur
lecture permet d'evaluer la qualite de la base :

- un score eleve traduit une signature nette et specifique
- un score faible signale une signature peu discriminante, frequente pour les
  types apparentes formant un continuum
- une etiquette `Unknown` indique qu'aucun type n'a atteint le seuil de
  confiance, soit que le type manque a la base, soit que le cluster correspond
  a un etat inattendu

La figure `18_dotplot_annotated.png` constitue le controle visuel : une
annotation coherente s'y lit comme une diagonale nette, chaque marqueur
representatif s'exprimant preferentiellement sur le type qu'il definit.

## Adapter la base a un nouveau tissu

Copier la base fournie et l'adapter :

    cp assets/canonical_markers/pbmc_markers_human.yaml \
       assets/canonical_markers/mon_tissu.yaml

Puis pointer vers elle dans la configuration d'analyse :

    annotation:
      perform: true
      markers_file: "assets/canonical_markers/mon_tissu.yaml"
      score_threshold: 0.0

Le seuil de confiance est exprime dans l'echelle des scores de module, centres
autour de zero. Un seuil a 0.0 accepte toute signature exprimee au dela du fond
attendu. L'augmenter rend l'annotation plus exigeante et produit davantage
d'etiquettes `Unknown`.
