##############Genome-wide phylogenomics analysis##############
=pod
   The scripts implemented the phylogenomics analysis using ASTRAL-II based on gene trees using 
one-to one orthologs between related species
=cut
###################alignment block################################
 #!/usr/bin/perl -w
  use strict;
  use warnings;
  die "USAGE:\n $0 <pro_sequences> <cds_sequences> <orthologs_list>" unless(@ARGV==3);
  open(PRO, "<$ARGV[0]")||die "$!";
  open(CDS, "<$ARGV[1]")||die "$!";
  my %pro;my %cds;
  my $name1; my $seq1;
  my $name2; my $seq2;
  while(<PRO>){
       chomp $_;
       if(/>(.*)/){
          $name1=(split/\|/,$1)[1];
          $seq1=();
         }
       else{
         $seq1.=$_;
         $pro{$name1}=$seq1;
         }
  }
  close PRO;
   while(<CDS>){
        chomp $_;
        if(/>(.*)/){
          $name2=(split/\|/,$1)[1];
          $seq2=();
        }
        else{
          $seq2.=$_;
          $cds{$name2}=$seq2;
       }
  }
  close CDS;
###OrthG1200       9       9       ath|AT5G57040.1 cme|MELO3C015494P1 egr|Eucgr.C03409.1.p
### grm|Gorai.010G084800.1 mec|Mecan_P01885.1 mtr|Medtr4g125860.1 pca|Potri.018G056200
  my $orthid;
  my $proid;
  my $species;
  my @orth_files;
  my $transid;
  my $trimal;
     open(ORTH,"<$ARGV[2]")||die "$!";
      while(<ORTH>){
          chomp $_;
          my @items=split/\t/,$_;
             $orthid=$items[0];
          my @orths=split/ /,$items[3];
          my @orths_sorted=sort @orths;
             open (ORTHPRO, ">$orthid.fasta");
             open (ORTHTRS, ">$orthid.cds");
             map{
                  ($species,$proid)=split/\|/,$_;  ###/correspondings between transcripts and proteins
                  $transid=$proid;
                  print ORTHPRO ">$species\n$pro{$proid}\n";
                  print ORTHTRS ">$species\n$cds{$transid}\n";
              }@orths_sorted;
             close ORTHPRO;
             close ORTHTRS;
          my $mafft="mafft --maxiterate 1000 --localpair --clustalout --quiet --thread 5 $orthid.fasta > $orthid.alignment" ;
             system($mafft);
         my $pal2nal="./pal2nal.v14/pal2nal.pl $orthid.alignment $orthid.cds -output clustal>$orthid.cds_align";
            system($pal2nal);
            unlink("$orthid.fasta");
            unlink("$orthid.cds");
            unlink("$orthid.alignment");
            if(not (-z "$orthid.cds_align")){
                 $trimal="./trimal/source/trimal -in $orthid.cds_align -out $orthid.cds_trimal -fasta";
                 system($trimal);
                 unlink("$orthid.cds_align");
             my $raxml="./standard-RAxML/raxmlHPC-PTHREADS -f ad -p 12345 -x 12345 -s $orthid.cds_trimal  -# 100 -m GTRGAMMAI -n $orthid.tre -k -T 30";
                system($raxml);
            }
     }
     close ORTH;
  my $mkdir="mkdir gene_trees";
     system($mkdir);
  my $move="mv RAxML_bootstrap.$orthid.tre gene_trees/";
     system($move); 
  my $cat="cat gene_trees/*.tre>raxml_bs.tre";
     system($cat);
  my $astral="java -jar ASTRAL/Astral/astral.5.6.2.jar -i raxml_bs.tre -o Lq_astral.tre";
     system($astral);

