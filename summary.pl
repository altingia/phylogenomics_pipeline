#!/usr/bin/perl -w
use strict;
use warnings;
die "USAGE:\n $0 <groups.txt> <groups1.txt> <summary.txt>" unless (@ARGV==3);
open(IN, "<$ARGV[0]")||die "$!";
open(OUT1, ">$ARGV[1]")||die "$!";
open(OUT2, ">$ARGV[2]")||die "$!";

my @sp=("amtri","artha","cedem","eufer","gibil","lichi","lifor","metru","muacu","orsat","peame","potri","vivin");
   print OUT2 "OrthID\tamtri\tartha\tcedem\teufer\tgibil\tlichi\tlifor\tmetru\tmuacu\torsat\tpeame\tpotri\tvivim\n";
while(<IN>){
         chomp $_;
      my ($orthid,$orths)=split/: /,$_;
      my @orths=split/ /,$orths;
      my %sp;
      my @orths1;
         map { 
               my $sp=(split/\|/,$_)[0];
                  if($sp ne "pinig"){
                     $sp{$sp}++;
                     push @orths1,$_;
                  }
             }@orths;                
         @orths=();
      my $orths1=join "\t",@orths1;
         @orths1=();  
      my $sp_num=keys %sp;
      my $gene_num=0;
         foreach my $key (sort keys %sp){
                   $gene_num+=$sp{$key};
        }
         print OUT1 "$orthid\t$sp_num\t$gene_num\t$orths1\n";
      my @num;   
         foreach my $sp(@sp){
                 if(exists $sp{$sp}){
                   push @num,$sp{$sp};
                  }
                 else{
                   push @num,"0";
                 }
          }
       my $items=join "\t",($orthid,@num);
          print OUT2 "$items\n";
          @num=();
          %sp=();
  }

 close IN;
 close OUT1;
 close OUT2;          
