02\_data-analysis
================

  - [appeler la librairie dada2](#appeler-la-librairie-dada2)
  - [Inspect read quality profiles](#inspect-read-quality-profiles)
  - [Filter and trim](#filter-and-trim)
  - [appretissage des erreurs](#appretissage-des-erreurs)
  - [sample inference](#sample-inference)
  - [mairged paired reads](#mairged-paired-reads)
  - [construire la table
    d’observation](#construire-la-table-dobservation)
  - [Remove chimeras](#remove-chimeras)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
  - [assignation de la taxonomie](#assignation-de-la-taxonomie)

# appeler la librairie dada2

``` r
library("dada2")
```

    ## Warning: multiple methods tables found for 'which'

Les données utilisées sont dans le fichier donnees que j’ai crée et j’ai
mis à l’intérieur toutes les données importées dedans. La commande
suivant permet de définir une variable path pour qu’elle pointe vers le
répertoire donnees où j’ai mis l’ensemble des données extrait sur la
machine. Pour l’instant, considérer simplement les fichiers fastq
appariés à traiter.

``` r
path <- "~/ecog2_cc2/donnees" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

cela permet de lire les noms des fichiers fastq et effectuer quelques
manipulations de chaînes pour obtenir des listes correspondantes des
fichiers fastq forward et reverse. On va créer une variable fnFs et on
lui assigne la valeur de résultat de la fonction sort() qui va classer
les résultats de la fonction list.files(), qui va lister le fichier
R1\_001.fastq. ensuite t va faire la même chose pour R2\_001.fastq. En
suite on va extrairer les sample names avec la fonction strsplit(), en
supposant que les noms de fichiers ont le format: SAMPLENAME\_XXX.fastq

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_1.fastq and SAMPLENAME_R2_1.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "R"), `[`, 1)
```

# Inspect read quality profiles

Nous commençons par visualiser les profils de qualité des Forward reads
en utilisant la fonction plotQualityProfile qui va permettre de tracer
un résumé visuel de la distribution des scores de qualité en fonction de
la position de la séquence pour le fichier: fnFs fastq d’entrée.

``` r
library(dada2)
plotQualityProfile(fnFs[1:2])
```

![](02_data-analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Les lectures forward sont de bonne qualité, avant 240. on va rogner
ensuite les derniers nucléotides pour éviter des erreurs moins bien
contrôlées qui peuvent s’y produire. on va tronquer donc les lectures
avant à la position 240 (en coupant les 10 derniers nucléotides).

on visualise le profil de qualité des reverse reads en utilisant la
fonction plotQualityProfile qui va permettre de tracer un résumé visuel
de la distribution des scores de qualité en fonction de la position de
la séquence pour le fichier: fnRs fastq d’entrée.

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_data-analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
Comentaire de graphe au dessus: sur ce graphe on remarque que la qualité
est moins bonne que celle du graphe des forwards, on tranque donc les
lectures inversées à la position 200 où la distribution de qualité se
bloque.

# Filter and trim

Attribuer les noms de fichiers aux fichiers fastq.gz filtrés. on crée
une variable filtFs et on met dedans \_F\_filt.fastq.gz et puis une
filtRs et on met dedans \_R\_filt.fastq.gz après on nome l’objet filtFs
en utilisant la fonction names on fait pareil pour filtRs et on leur
donne le nom sample.names qui est la valeur de la fonction.

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Nous utiliserons les paramètres de filtrage standard: maxN=0(DADA2 ne
nécessite aucun Ns) truncQ=2, rm.phix=TRUEet maxEE=2. Le paramètre
maxEEp définit le nombre maximum d ’«erreurs attendues» autorisées dans
une lecture. Ici, on crée une variale out et on lui assigne les valeurs
des résultats de la fonction filterAndTrim(), qui va filtrer et ajuster
les fichiers fnFs, filtFs, fnRs, filtRs de fastq d’entrée (peut être
compressé) en fonction de plusieurs critères: maxN=0, maxEE=c(2,2),
truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE , et génère des
fichiers fastq (compressés par défaut) contenant les lectures coupées
qui ont passé les filtrés. Des fichiers fastq revers et forward
correspondants peuvent être fournis en entrée, auquel cas le filtrage
est effectué sur les lectures avant et arrière indépendamment, et les
deux lectures doivent passer pour que la paire de lecture soit sortie.
En suite on va utiliser la fonction head pour avoir un apperçu de
l’objet out

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),trimLeft=21,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145448
    ## Station5_Fond1_11mars15_R1.fastq     175993    160423
    ## Station5_Fond2_10sept14_R1.fastq     197039    177018
    ## Station5_Fond2_11mars15_R1.fastq      87585     79989
    ## Station5_Fond3_10sept14_R1.fastq     117140    106150
    ## Station5_Median1_10sept14_R1.fastq   116519    106745

# appretissage des erreurs

dada2 calcul un model d’erreurs apartir des données de séquençage, cette
méthode sur les reads F Reverse. L’algorithme DADA2 utilise un modèle
d’erreur paramétrique ( err) et chaque jeu de données d’amplicon a un
ensemble différent de taux d’erreur. La learnErrors() méthode apprend ce
modèle d’erreur à partir des données, en alternant l’estimation des taux
d’erreur et l’inférence de la composition de l’échantillon jusqu’à ce
qu’ils convergent vers une solution cohérente conjointement.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 105752691 total bases in 482889 reads from 3 samples will be used for learning the error rates.

dada2 calcul un model d’erreurs apartir des données de séquençage, on
applique cette méthode sur les reads Reverse de la même manière que pour
errF.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100755162 total bases in 562878 reads from 4 samples will be used for learning the error rates.

La fonction plotErrors() va tracer la fréquence observée de chaque
transition (par exemple A-\> C) en fonction du score de qualité associé.
Il trace également les taux d’erreur estimés finaux (s’ils existent).
l’argument nominalQ=TRUE va permettre de tracer les taux d’erreur
attendus (ligne rouge) si les scores de qualité correspondent exactement
à leur définition nominale: Q = -10 log10 (p\_err).

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_data-analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Intérprétation de graphe précédent: graphe 1: la probabilité que A est A
doit être max ainsi de suite… les taux d’erreur pour chaque transition
sont indiquées sur les graphes (graphe 2: A–\>C etc) c2: probabilité
d’erreus de séq pour qu’un A soit un C. Ici, les taux d’erreur estimés
(ligne noire) correspondent bien aux taux observés (points), et les taux
d’erreur diminuent avec une qualité accrue comme prévu. Tout semble
raisonnable et nous procédons en toute confiance.

# sample inference

on crée une nouvelle variable dadaFs pour corriger les jeux de données
dada appliquée au donné Forward. La fonction dada() supprime toutes les
erreurs de séquençage pour révéler les membres de la communauté
séquencée. l’argument multithread=TRUE: le multithreading est activé
et le nombre de threads disponibles est automatiquement déterminé. Si un
entier est fourni, le nombre de threads à utiliser est défini en
transmettant l’argument à setThreadOptions.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 37907 unique sequences.
    ## Sample 2 - 160423 reads in 35863 unique sequences.
    ## Sample 3 - 177018 reads in 47212 unique sequences.
    ## Sample 4 - 79989 reads in 20356 unique sequences.
    ## Sample 5 - 106150 reads in 30255 unique sequences.
    ## Sample 6 - 106745 reads in 28836 unique sequences.
    ## Sample 7 - 98823 reads in 25824 unique sequences.
    ## Sample 8 - 107427 reads in 26733 unique sequences.
    ## Sample 9 - 71082 reads in 17976 unique sequences.
    ## Sample 10 - 78645 reads in 20422 unique sequences.
    ## Sample 11 - 91534 reads in 24487 unique sequences.

dada appliqué au donné reverse. on crée une nouvelle variable dadaRs et
on lui assigne la valeur de résultat de la fonction dada comme on a fait
pour les reverse à fin de corriger les jeux de données appliquées pour
les Reverse.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 45486 unique sequences.
    ## Sample 2 - 160423 reads in 41638 unique sequences.
    ## Sample 3 - 177018 reads in 55554 unique sequences.
    ## Sample 4 - 79989 reads in 23239 unique sequences.
    ## Sample 5 - 106150 reads in 34625 unique sequences.
    ## Sample 6 - 106745 reads in 31673 unique sequences.
    ## Sample 7 - 98823 reads in 29093 unique sequences.
    ## Sample 8 - 107427 reads in 28947 unique sequences.
    ## Sample 9 - 71082 reads in 21426 unique sequences.
    ## Sample 10 - 78645 reads in 22051 unique sequences.
    ## Sample 11 - 91534 reads in 28266 unique sequences.

dada a crée des objet de class dada: dadaFs et dadaRs on regarde ce qui
est dans le premier étagers du placard de dadaFs, on peut changer le 1
pour regarder à n’importe quel étage de dadaFS.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1010 sequence variants were inferred from 37907 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# mairged paired reads

on va merger Read 1 et read 2 pour obtenir les séquences entièrement
débruitées en utilisant la fonction mergePairs(). La fusion est
effectuée en alignant les lectures avant débruitées avec le complément
inverse des lectures inverses débruitées correspondantes, puis en
construisant les séquences «contig» fusionnées. verbose=TRUE:monter les
étape avec le texte head(mergers\[\[1\]\]):inspecter les résultats en
regardant les premières lignes

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 117318 paired-reads (in 5196 unique pairings) successfully merged out of 141000 (in 21451 pairings) input.

    ## 138940 paired-reads (in 4296 unique pairings) successfully merged out of 156462 (in 15709 pairings) input.

    ## 142188 paired-reads (in 6989 unique pairings) successfully merged out of 171439 (in 27056 pairings) input.

    ## 67622 paired-reads (in 2721 unique pairings) successfully merged out of 77764 (in 9556 pairings) input.

    ## 83613 paired-reads (in 3458 unique pairings) successfully merged out of 102224 (in 16304 pairings) input.

    ## 86212 paired-reads (in 3348 unique pairings) successfully merged out of 103447 (in 14293 pairings) input.

    ## 80661 paired-reads (in 2727 unique pairings) successfully merged out of 95866 (in 12350 pairings) input.

    ## 89385 paired-reads (in 3073 unique pairings) successfully merged out of 104354 (in 12135 pairings) input.

    ## 59716 paired-reads (in 1939 unique pairings) successfully merged out of 68711 (in 7974 pairings) input.

    ## 66157 paired-reads (in 1763 unique pairings) successfully merged out of 76701 (in 8283 pairings) input.

    ## 75048 paired-reads (in 3149 unique pairings) successfully merged out of 88514 (in 12054 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                sequence
    ## 1     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 2     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 3     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 4     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 5     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 6 TACGAGGGGTCCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTACGTAGGCGTTTTAATAAGTTGTATGTTAAATATCTTAGCTTAACTAAGAAAGTGCATACAAAACTGTTAAGATAGAGTTTGAGAGAGGAACGCAGAATTCATGGTGGAGCGGTGACATGCGTAGATATCATGAGGAAAGTCAAATGCGAAGGCAGCCTTCTGGCTCAAAACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTATTTGGTGCTGGGGGATTCGACCCTTTCAGTGCCGTAGCTAACGCGATAAATACTCCGCCTGGGGACTACGATCGCAAGATT
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      5170       1       2     29         0      0      2   TRUE
    ## 2      4129       2       1     29         0      0      2   TRUE
    ## 3      3757       3       1     29         0      0      2   TRUE
    ## 4      2481       1       1     29         0      0      2   TRUE
    ## 5      2182       2       2     29         0      0      2   TRUE
    ## 6      2132       5       9     25         0      0      1   TRUE

Donc là on a les reads fusionnés avec succès.

# construire la table d’observation

A partir des merged, on crée un nouvelle objet seqtab La fonction
makeSequenceTable(), va permettre de construire une table de séquence
(analogue à une table OTU) à partir de la liste d’échantillon mergers.
la fonction dim va permettre de récupérer l’objet seqtab.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]    11 19426

Ici on veut regarder la distribution de la longeur des séquences. la
fonction getSequences() va extraire les séquences de l’objet seqtab.
nchar() va prenddre le résulat de la fonction getSequences comme
argument et renvoie un vecteur dont les éléments contiennent les tailles
des éléments correspondants. table() va permettre ensuite de créer un
tableau de tout cela.

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  183   27  165  184 5608 3594 2312 2613 2738  126 1770 
    ##  376  377  378  382  386 
    ##   90    4    1    1    2

dans mon objet seqtab j’ai une séq qui fait 352 paires de bases, une
autre qui fait 353 etc

# Remove chimeras

Une chimère: ça se passe pendant l’amplification par PCR, donc par ex un
ADN 16s amplifier par un fragment reverse et forrward. Si on prend que
le forward,y aura élongation mais imaginant qu’elle s’arréte avant la
fin de la séq 16S. Donc on va avoir le fragment 16S qui n’a pas bouger
et un fragment non complet. donc après le 2ème cycle, le fragment non
complet va pouvoir s’hybrider avec un 16s d’une autre bactérie, et
l’élongation va continuer. on va avoir comme résultat au final un
fragment hybride qui provient du premier ARN 16 et du deuxième. cela
s’appelle chimère. On va éliminer ces chimères en utlisant la fonction
removeBimeraDenovo(), et on va donner la valeur de résultat de la
fonction à une nouvelle variable appelée seqtab.mochim le système va
regarder tout les séq rares dans le début contig correspont au premier
ARN et la fin au deuxième.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

    ## [1]   11 1557

Donc le résultat montre qui il a détecter 17869 (=19426-1557) chimère
sur 19426.

calcul de ratio chimère qui est égale à la somme des sequences se
trouvant dans l’objet seqtab.mochim (c’est l’objet après remove des
chimères) / somme des séquences se traouvent dans l’objet seqtab(avant
d’enlever les chimères)

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.7769154

y a (1-0.96)\*100 = 22.31% de séq chimériques dans notre jeu de données.

# Track reads through the pipeline

résumé des fichiers qualité. construire une table on crée un nouvel
objet qui est getN c’est une variable qui prend le rôle d’une fonction
apliquer la fonnction get N qui est la somme des get uniq de dadaFS la
table track va être la concaténation de tout ce qui est entre parentèse.
head(track) : pour visualiser le tableau

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##                             input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_   159971   145448    142931    143292 117318   87962
    ## Station5_Fond1_11mars15_   175993   160423    158128    158473 138940  111552
    ## Station5_Fond2_10sept14_   197039   177018    173601    174591 142188  103668
    ## Station5_Fond2_11mars15_    87585    79989     78618     78926  67622   54711
    ## Station5_Fond3_10sept14_   117140   106150    103806    104338  83613   64259
    ## Station5_Median1_10sept14_ 116519   106745    104811    105173  86212   65559

# assignation de la taxonomie

Il va regarder dans le base de données et à partir des séq qui sont
proches, et va assigner une taxonomie.

c’est une façon d’attribuer une taxonomie aux séquences. La fonction
assignTaxonomy() prend en entrée un ensemble de séquences à classer et
un ensemble d’apprentissage de séquences de référence avec une taxonomie
connue (silva en ce cas), et produit des affectations taxonomiques avec
au moins une minBootconfiance bootstrap.

``` r
library(dada2)
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

ajout d’espèces dans le répertoire taxa contenant les fichiers fastq.

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

Inspecter les affectations taxonomiques:

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                     Species
    ## [1,] "Clade I"          "Clade Ia"                NA     
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"    NA     
    ## [3,] "Clade I"          "Clade Ia"                NA     
    ## [4,] "Clade I"          "Clade Ia"                NA     
    ## [5,] "Clade II"         NA                        NA     
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina" NA

les Proteobacteria sont bien représentés parmi les taxons les plus
abondants dans ces échantillons de la rade de Brest. Pas d’attributions
d’espèces, parce qu’il est souvent impossible de faire des assignations
d’espèces sans ambiguïté à partir de sous-segments du gène 16S.

``` r
save.image(file="02_data-analysis")
```
