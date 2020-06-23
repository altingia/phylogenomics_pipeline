##############Genome-wide phylogenomics analysis##############
=pod
   The scripts implemented the phylogenomics analysis with  4dtv sites 
of one-to one orthologs between related species
=cut
###################alignment block################################
 #!/usr/bin/perl -w
  use strict;
  use warnings;
  use Bio::AlignIO;
  use Bio::AlignIO::clustalw;
  use Bio::AlignIO::fasta;
  use Bio::SeqIO;
  die "USAGE:\n $0 <pro_sequences> <cds_sequences> <orthologs_list><fourfold_sites_alignment>" unless(@ARGV==4);
############ read the protein and transcripts files################################
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
               my $seqname=join '_',($species,$orthid);
                  print ORTHPRO ">$seqname\n$pro{$proid}\n";
                  print ORTHTRS ">$seqname\n$cds{$transid}\n";
            }@orths_sorted;
    my $mafft="mafft --maxiterate 1000 --localpair --clustalout --quiet --thread 5 $orthid.fasta > $orthid.alignment" ;
       system($mafft);
    my $pal2nal="/home/wuwei/software/pal2nal.v14/pal2nal.pl $orthid.alignment $orthid.cds -output clustal>$orthid.cds_align";
       system($pal2nal);
       unlink("$orthid.fasta");
       unlink("$orthid.cds");
       unlink("$orthid.alignment");
    if(not (-z "$orthid.cds_align")){
              $trimal="/home/wuwei/software/trimal/source/trimal -in $orthid.cds_align -out $orthid.cds_trimal -fasta -gt 0.9 -cons 80";
              system($trimal);
              unlink("$orthid.cds_align");
    #############extraction of 4dtv sites###################################################################
         my %fourfoldsites=(
                       'TC'=> 'Ser',
                       'CT'=> 'Leu',
                       'CC'=> 'Pro',
                       'CG'=> 'Arg',
                       'AC'=> 'Thr',
                       'GT'=> 'Val',
                       'GC'=> 'Ala',
                       'GG'=> 'Gln',
                      );
         my %freq;
         my %seq;
         my @sp;
         my $trim_cds=Bio::SeqIO->new(-file=>"$orthid.cds_trimal",-format=>"fasta");
            while(my $seqs=$trim_cds->next_seq){
                 my $id=$seqs->id;
                    push @sp,$id;
                 my $seq=$seqs->seq;
                 my @seq=split//,$seq;
                 my $condon;
                    for(my $i=0;$i<@seq;$i+=3){
                           $condon=$seq[$i].$seq[$i+1];
                            if(exists $fourfoldsites{$condon}){
                               $seq{$id}{$i}=$seq[$i+2];
                               $freq{$i}++;
                          }
                   }
                     @seq=();
           }
            unlink("$orthid.cds_trimal"); 
           open (FOURFOLD,">$orthid.fourfold.align")||die "$!";
           foreach my $sp (@sp){
                   my $fourfoldseqs=();
                      foreach my $site (sort keys  %freq){
                          if($freq{$site}==($#sp+1)){
                           $fourfoldseqs.=$seq{$sp}{$site};
                         }
                   }
                  print FOURFOLD ">$sp\n$fourfoldseqs\n";
          }
           @sp=();
           push @orth_files,"$orthid.fourfold.align";
           close FOURFOLD;
       }
      else{
       next;
      }
    }
     my @orth_4dtv;
     my %fourfold_seqs;
     my %all;
     my $j=0;
     my (@loci,@species);
       open(ALLOCI,">$ARGV[3]")||die "$!";
 foreach my $file(@orth_files){
            open LOCI, "<$file" or die "$!";
            while(<LOCI>){
                  chomp $_;
                   if(/>(.*)/){
                   push @orth_4dtv,(split/_/,$1)[0];
                   }
                   else{
                   push @orth_4dtv,$_;
                   }
          }
            unlink("$file");
            close LOCI;
            push @loci,$j;
            %fourfold_seqs=@orth_4dtv;
            @species= keys %fourfold_seqs;
            @orth_4dtv=();
            foreach my $key(sort keys %fourfold_seqs){
               $all{$key}{$j}=$fourfold_seqs{$key};
            }
            %fourfold_seqs=();
                $j++;
   }
    my $concat;
       for my $sp (sort @species){
           for my $loci (@loci){
               $concat.=$all{$sp}{$loci};
           }
               print ALLOCI ">$sp\n$concat\n";
               $concat=();
      }
       close ALLOCI;
   my $jmodeltest="java -jar /home/wuwei/software/jmodeltest-2.1.10/jModelTest.jar -d all_species_fourfold.txt -o jmodeltest.out -i -f -g 4     -BIC -AIC -AICc -DT -v -a -w";
      system($jmodeltest); 
   my $raxml="/home/wuwei/software/standard-RAxML/raxmlHPC-PTHREADS -f ad  -p 12345 -x 12345 -s all_species_fourfold.txt -o atr  -# 1000 -m     GTRGAMMAI -n raxml_mecThreads";
     system($raxml);
