01\_data-import
================

  - [Downloading data](#downloading-data)
  - [Decompress the data](#decompress-the-data)
  - [charger la base de donnée de
    silva](#charger-la-base-de-donnée-de-silva)
      - [Lier l’annotation d’espèce à l’annotation
        taxonomique](#lier-lannotation-despèce-à-lannotation-taxonomique)

# Downloading data

Télécharger les données sur la Home page de M. Maignien. Les
échatillions sont prélevés de la rade de Brest, de différentes
profondeurs, et différents date.

``` bash
wget https://pagesperso.univ-brest.fr/~maignien/teaching/M1-MFA/UE-Ecogenomique2/EcoG2_data_cc2.tar.gz
```

# Decompress the data

la commande suivante sert à compresser les data.

``` bash
tar xvzf EcoG2_data_cc2.tar.gz
```

# charger la base de donnée de silva

SILVA fournit des ensembles de données complets,contrôlés par la qualité
et régulièrement mis à jour de séquences d’ARN ribosomal (ARNr) alignées
de petite (16S / 18S, SSU) et de grande sous-unité (23S / 28S, LSU) pour
les trois domaines de la vie (bactéries, archées et eucarya).

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

## Lier l’annotation d’espèce à l’annotation taxonomique

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```
