#!/usr/bin/env perl
#
# The MAP file was used to create a BED file of the SNP coordinates
# cat Illumina_Omni2.5_genotypes.map | awk '{OFS="\t"; print "chr"$1, $4-1, $4, $2, 0, "+"}' | gzip > omni.bed.gz
#
# These coordinates were lifted over to hg19 coordinates:
#
# liftOver omni.bed.gz ~/ucsc/hg18ToHg19.over.chain omni_hg19.bed omni_hg19_unmapped.out
# gzip omni_hg19.bed
#
# This script stores the lifted over coordinates from omni_hg19.bed.gz
# and reads in the hg18 VCF file, to create a hg19 VCF file
#
# Usage:
#
# hg18_to_hg19.pl > Illumina_Omni2.5_genotypes_hg19.vcf
#

use strict;
use warnings;

my $hg19 = 'omni_hg19.bed.gz';
my %hg19 = ();

open(IN,'-|',"gunzip -c $hg19") || die "Could not open $hg19: $!\n";
while(<IN>){
   chomp;
   #chr1    82153   82154   rs4477212       0       +
   my ($chr, $start, $end, $id, $score, $strand) = split(/\t/);
   $hg19{$id} = $end;
}
close(IN);

my $vcf = 'Illumina_Omni2.5_genotypes.vcf';
my $missing = 0;

open(IN,'<',$vcf) || die "Could not open $vcf: $!\n";
while(<IN>){
   chomp;
   if (/^#/){
      if (/^#CHROM/){
         s/\b\d+_//g;
         s/_/-/g;
         print "$_\n";
         next;
      } else {
         print "$_\n";
         next;
      }
   }
   my ($chr, $pos, $id, @rest) = split(/\t/);
   if (exists $hg19{$id}){
      print join("\t", $chr, $hg19{$id}, $id, @rest), "\n";
   } else {
      ++$missing;
      next;
   }
}
close(IN);

warn("$missing missing IDs\n");

exit(0);
