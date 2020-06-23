# phylogenomics_pipeline
Phylogeny inference based on genome-wide data is a routine step for newly sequenced genomes.After gene family clustering implemented in Orthofinder or OrthoMCL, we can reconstruct the phlogeny using both single copy nuclear genes or low-copy nuclear genes. We found that the available species tree derived from built-in program STAG of Orthofinder was not accurate among closely related species, ie,taxon within a family. three independent perl scripts 
a. phylogenomics_4dtv.pl. phylogeny inference using four-fold degenerated sites of concatenate loci.
b. phylogenomics_astral.pl phylogeny inference using raxml gene tree,and summary methods colaescent methods implented in ASTRAL-II based single copy nuclear gene.
c. phylogenomics_para.pl. phylogeny inference using low copy nuclear genes implemented in ASTRAL-Pro (Zhang et al 2019), scripts or software used in this pipleline cited as:

Mikita Suyama, David Torrents, and Peer Bork (2006) PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments.Nucleic Acids Res. 34, W609-W612.

Alexandros Stamatakis, RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies, Bioinformatics, Volume 30, Issue 9, 1 May 2014, Pages 1312–1313, https://doi.org/10.1093/bioinformatics/btu033

Siavash Mirarab, Tandy Warnow, ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes, Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages i44–i52, https://doi.org/10.1093/bioinformatics/btv234

Mikita Suyama, David Torrents, Peer Bork, PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments, Nucleic Acids Research, Volume 34, Issue suppl_2, 1 July 2006, Pages W609–W612, https://doi.org/10.1093/nar/gkl315

Salvador Capella-Gutiérrez, José M. Silla-Martínez, Toni Gabaldón, trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses, Bioinformatics, Volume 25, Issue 15, 1 August 2009, Pages 1972–1973, https://doi.org/10.1093/bioinformatics/btp348

Kazutaka Katoh, Daron M. Standley, MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability, Molecular Biology and Evolution, Volume 30, Issue 4, April 2013, Pages 772–780, https://doi.org/10.1093/molbev/mst010

Stajich JE, Block D, Boulez K, et al. The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 2002;12(10):1611-1618. doi:10.1101/gr.361602

Installation:
1. the pipeline phylogenomics_4dtvs.pl need BioPerl
  1a. >perl -MCPAN -e shell
      cpan>install Bundle::CPAN
      cpan>q
      
1b. upgrade cpan
    >cpan
    cpan>install Module::Build
    cpan>o conf prefer_installer MB
    cpan>o conf commit
    cpan>q
    
1c. install two modules
    cpan>install Bio::SeqIO
    cpan>install Bio::AlignIO
    cpan>install Bio::AlignI
other binary scripts  were ready for use.

the orhogroups_list can be obtained used the perl script summary.pl

