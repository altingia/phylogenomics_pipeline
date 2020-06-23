##############Genome-wide phylogenomics analysis##############
=pod
   The scripts implemented the phylogenomics analysis with low copy genes
(here, we contracted the gene family size with copies less than 5)
references:
 1. https://github.com/chaoszhang/A-pro
=cut
#!/usr/bin/perl -w
  use strict;
  use warnings;
  die "USAGE:\n $0 <pro_sequences> <cds_sequences> <orthologroups_list>" unless(@ARGV==3);
############ read the protein and transcripts files################################
  open(PRO, "<$ARGV[0]")||die "$!";
  open(CDS, "<$ARGV[1]")||die "$!";
  my %pro;
  my %cds;
  my $name1; 
  my $seq1;
  my $name2; 
  my $seq2;
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
  my $transid;
  my @astral_trees;
  my $trimal;
     open(ORTH,"<$ARGV[2]")||die "$!";
      while(<ORTH>){
          chomp $_;
          my @items=split/\t/,$_;
             $orthid=$items[0];
          my @orths=split/ /,$items[3];
          my @orths_sorted=sort @orths;
          my %sp;
             open (ORTHPRO, ">$orthid.fasta");
             open (ORTHTRS, ">$orthid.cds");
             map{
                  ($species,$proid)=split/\|/,$_;  ###/correspondings between transcripts and proteins
                  $transid=$proid;
                  print ORTHPRO ">$proid\n$pro{$proid}\n";
                  print ORTHTRS ">$proid\n$cds{$transid}\n";
                  $sp{$proid}=$species;
              }@orths_sorted;
             close ORTHPRO;
             close ORTHTRS;
          my $mafft="mafft --maxiterate 1000 --localpair  --quiet --thread 5 $orthid.fasta > $orthid.alignment" ;
             system($mafft);
         my $pal2nal="./pal2nal.v14/pal2nal.pl $orthid.alignment $orthid.cds -output fasta>$orthid.cds_align";
            system($pal2nal);
            #unlink("$orthid.fasta");
           # unlink("$orthid.cds");
           # unlink("$orthid.alignment");
            if(not (-z "$orthid.cds_align")){
                 $trimal="./trimal/source/trimal -in $orthid.cds_align -out $orthid.cds_trimal -fasta";
                 system($trimal);
                 unlink("$orthid.cds_align");
             my $raxml="./standard-RAxML/raxmlHPC-PTHREADS -f ad -p 12345 -x 12345 -s $orthid.cds_trimal  -# 100 -m GTRGAMMAI -n $orthid.tre -k -T 30";
                system($raxml);
             my @besttree;
                open(BEST,"<RAxML_bestTree.$orthid.tre");
                while(<BEST>){
                  chomp $_;
                  push @besttree,$_;
                }
                close BEST;
             my $best=join '',@besttree;
                foreach my $key (keys %sp){
                 my $name=$sp{$key};
                    $best=~s/$key/$name/g;
                    push @astral_trees,$best;
               }
                @besttree=();
           }
} 
     close ORTH;
 my $temp=join '\n',@astral_trees;
    open (INPUT, "<$temp");
    open (OUTPUT,">astral-pro_input.tre");
    while(<INPUT>){
        chomp $_;
       print OUTPUT "$_\n";
    }
   close INPUT;
   close OUTPUT;
 my $astral_pro="java -D"java.library.path=./A-pro/ASTRAL-MP/lib" -jar /A-pro/ASTRAL-MP/astral.1.1.2.jar -i astral-pro_input.tre -o astral_output.tre";
    system($astral_pro);

