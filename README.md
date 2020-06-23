# phylogenomics_pipeline
Phylogeny inference based on genome-wide data is a routine step for newly sequenced genomes.After gene family clustering implemented in Orthofinder or OrthoMCL, we can reconstruct the phlogeny using both single copy nuclear genes or low-copy nuclear genes. We found that the available species tree derived from built-in program STAG of Orthofinder was not accurate among closely related species, ie,taxon within a family. three independent perl scripts 
a. phylogenomics_4dtv.pl. phylogeny inference using four-fold degenerated sites of concatenate loci.
b. phylogenomics_astral.pl phylogeny inference using raxml gene tree,and summary methods colaescent methods implented in ASTRAL-II based single copy nuclear gene.
c. phylogenomics_para.pl. phylogeny inference using low copy nuclear genes implemented in ASTRAL-Pro (Zhang et al 2019) 
