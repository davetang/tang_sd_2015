#!/usr/bin/env perl
#
# Download the strand information from http://www.well.ox.ac.uk/~wrayner/strand/
# Download the hg19 genome via the UCSC Genome Browser website
#
# Usage:
#
# check_hg19_base.pl ~/genome/hg19/hg19.fa Illumina_Omni2.5_genotypes_hg19.vcf > illumina_fixed.vcf
#

use strict;
use warnings;

my $usage = "Usage: $0 <genome.fa> <infile.vcf>\n";
my $genome = shift or die $usage;
my $infile = shift or die $usage;

my %genome = ();
warn "Storing $genome\n";
open(IN,'<',$genome) || die "Could not open $genome: $!\n";
my $chr = '';
while(<IN>){
   chomp;
   if (/^>(.*)$/){
      $chr = $1;
      $chr =~ s/^chr//;
   } else {
      $genome{$chr} .= $_;
   }
}
close(IN);
warn "Done storing $genome\n";

my $sf = '/home/san/dtang/project/wiluna/illumina/HumanOmni2.5-4v1_D-b37.strand';
my %strand = ();
warn "Storing $sf\n";
open(IN,'<',$sf) || die "Could not open $sf: $!\n";
while(<IN>){
   chomp;
   #kgp10000188     13      45106083        100     +       AG
   my ($id, $chr, $pos, $percent, $strand, $allele) = split(/\t/);
   $chr =~ s/^chr//;
   $strand{$chr}{$pos}->{'STRAND'} = $strand;
   $strand{$chr}{$pos}->{'ALLELE'} = $allele;
}
close(IN);
warn "Done storing $sf\n";

my $match = 0;
my $mismatch = 0;
my $missing_strand_info = 0;
my $strand_info = 0;

open(IN,'<',$infile) || die "Could not open $infile: $!\n";
ENTRY: while(<IN>){
   chomp;
   if (/^#/){
      print "$_\n";
      next;
   }

   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
   my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @rest) = split(/\t/);
   $chr =~ s/^chr//;

   # first things first, obtain the reference base at $pos according to $genome
   my $base = '';
   if (exists $genome{$chr}){
      $base = substr($genome{$chr}, $pos-1, 1);
      $base = uc($base);
   } else {
      die "Missing $chr in $genome\n";
   }

   # obtain information on the SNP
   my ($a, $b, $strand) = ('','','');
   if (exists $strand{$chr}{$pos}){
      ++$strand_info;
      $strand = $strand{$chr}{$pos}->{'STRAND'};
      my $allele = $strand{$chr}{$pos}->{'ALLELE'};
      $a = substr($allele, 0, 1);
      $b = substr($allele, 1, 1);
   } else {
      ++$missing_strand_info;
   }

   # on the positive strand
   if ($strand eq '+'){

      # two scenarios: $ref matches $base (great!)
      if ($ref eq $base){
         if ($alt eq '.'){
            if ($a eq $ref){
               $alt = $b;
            } elsif ($b eq $ref){
               $alt = $a;
            } else {
               next ENTRY;
            }
         } else {
            # great!
         }
         print join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @rest), "\n";
      }
      # and $ref does not match $base
      else {
         # is it because the ALT and REF are switched around?
         for (my $i = 0; $i < scalar(@rest); ++$i){
            if ($rest[$i] =~ /0\/0/){
               $rest[$i] = '1/1';
            } elsif ($rest[$i] =~ /1\/1/){
               $rest[$i] = '0/0';
            }
         }
         if ($alt eq $base){
            print join("\t", $chr, $pos, $id, $alt, $ref, $qual, $filter, $info, $format, @rest), "\n";
         }
         # switched around but the reference allele is missing
         elsif ($alt eq '.'){
            # obtain the alternative allele from the strand file
            # if it matches the genome reference
            if ($ref eq $a && $b eq $base){
               print join("\t", $chr, $pos, $id, $b, $a, $qual, $filter, $info, $format, @rest), "\n";
            } elsif ($ref eq $b && $a eq $base){
               print join("\t", $chr, $pos, $id, $a, $b, $qual, $filter, $info, $format, @rest), "\n";
            } else {
               # 27 cases
            }
         } else {
            # 7 cases
         }
      }

   # on the negative strand
   } elsif ($strand eq '-'){
      $base =~ tr/ACGT/TGCA/;
      if ($ref eq $base){
         if ($alt eq '.'){
            if ($a eq $ref){
               $alt = $b;
            } elsif ($b eq $ref){
               $alt = $a;
            } else {
               next ENTRY;
            }
         } else {
            # great!
         }
         $ref =~ tr/ACGT/TGCA/;
         $alt =~ tr/ACGT/TGCA/;
         print join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @rest), "\n";
      }
      # and $ref does not match $base
      else {
         # is it because the ALT and REF are switched around?
         for (my $i = 0; $i < scalar(@rest); ++$i){
            if ($rest[$i] =~ /0\/0/){
               $rest[$i] = '1/1';
            } elsif ($rest[$i] =~ /1\/1/){
               $rest[$i] = '0/0';
            }
         }
         if ($alt eq $base){
            $ref =~ tr/ACGT/TGCA/;
            $alt =~ tr/ACGT/TGCA/;
            print join("\t", $chr, $pos, $id, $alt, $ref, $qual, $filter, $info, $format, @rest), "\n";
         } elsif ($alt eq '.'){
            ++$match;
            # obtain the alternative allele from the strand file
            # if it matches the genome reference
            if ($ref eq $a && $b eq $base){
               $a =~ tr/ACGT/TGCA/;
               $b =~ tr/ACGT/TGCA/;
               print join("\t", $chr, $pos, $id, $b, $a, $qual, $filter, $info, $format, @rest), "\n";
            } elsif ($ref eq $b && $a eq $base){
               $a =~ tr/ACGT/TGCA/;
               $b =~ tr/ACGT/TGCA/;
               print join("\t", $chr, $pos, $id, $a, $b, $qual, $filter, $info, $format, @rest), "\n";
            } else {
               # 38 cases
            }
         } else {
            # 0 cases
         }
      }
   }
   # no strand information
   else {
      next ENTRY;
   }
}
close(IN);

#print "Match: $match\tMismatch: $mismatch\n";
warn("Strand info: $strand_info\nMissing strand info: $missing_strand_info\n");

exit(0);
