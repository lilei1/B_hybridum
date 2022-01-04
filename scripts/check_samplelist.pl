#!/usr/bin/perl
##by Li Lei, 2020-04-01, Berkeley;
#this is to check if the samples in the Gordon et al. 2020 in the vcf file and how many of them are not in the vcf file?
#You need 
use strict;
use warnings;
use Data::Dumper;

my ($file1, $file2) = @ARGV;

my %hash; #define a hash to store the positions for each SNP;
open(OUT,  "$file1") or die "Could not open $file1";
#my $header = <OUT>;
#print "$header";
foreach my $row (<OUT>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        #my @tt = split(/\.(\d+)/,$rtemp[0]);
        #print "$rtemp[0]\n";
         $hash{$rtemp[0]} = 0;
}
close (OUT);
#print Dumper(\%hash);

open(FILE,  "$file2") or die "Could not open $file2";
#my $header = <FILE>;
#print "$header";
foreach my $row (<FILE>){
        chomp $row;
        my @temp = split(/\t/,$row);
        if (exists $hash{$temp[0]}){
            print "$temp[0]\tYes\n"; 
        }
        else{
            print "$temp[0]\tNo\n";
        }
}
close (FILE);
