#!/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

#This script takes a vcf file, uses samtools to pull out the XRQ region, then blasts the HA412 genome and pulls out hits that are high enough.
#It also prevents multiple sites from getting the same location (due to mapping issues). 
#The script is adapted from Greg Owen at https://github.com/owensgl/wild_gwas_2018/edit/master/xrqpos2ha412pos_bwa.pl
#Cited by https://www.nature.com/articles/s41586-020-2467-6
#Here I need to convert the SNPs mapped back to B.distachyon back to B.stacei reference. 
my $input_file = $ARGV[0]; #Should be a gzvcf file
my $tmp_prefix = $ARGV[1];

#SLURM COMMANDS to load before running this script:
#system("module load samtools");
#system("module load bcftools");
#system("module load bwa");

my $ha412_ref = "/global/dna/projectdirs/RD/reseq_store/genomes/plants/Brachypodium_stacei/versions/1.0/Bstacei_316_v1.0.fa"; #stacei
my $xrq_ref = "/global/dna/projectdirs/RD/reseq_store/genomes/plants/Brachypodium_distachyon/versions/3.1/Bdistachyon_314_v3.0.fa";#distachyon
my $bwa = "bwa";
my $bcftools = "bcftools";
my $bases_surrounding=100; #Number of bases before and after the target site for blasting.
my $counter = 0;
my $ncores = 100; #Cores for fastq writing
my $ncores_bwa = 20; #Cores for BWA
my $min_mq = 40;

my $pm = new Parallel::ForkManager($ncores);

my @contig_list;
open FASTA, $ha412_ref;
while(<FASTA>){
  chomp;
  if ($_ =~ m/\>/){
    my $contig = $_;
    $contig =~ s/\>//g;
    $contig =~ s/ /, /g;
    push(@contig_list,$contig);
  }
}



#Turn SNPs into fastq of surrounding region
open(INPUT, "gunzip -c $input_file |");
open my $fastq, '>',  "$tmp_prefix.fastq";
my $first_line;
print STDERR "Converting SNP positions to fastq...\n";
while(<INPUT>){
  chomp;
  if ($_ =~ m/^#/){next;}
  my @a = split(/\t/,$_);
  $counter++;
  if ($counter > 1){$first_line++;}
  my $chr = $a[0];
  if ($chr eq "chr"){next;}
  $pm->start and next;
  my $pos = $a[1];
  my $start = $pos - $bases_surrounding;
  my $end = $pos + $bases_surrounding;
  if ($start < 0){next;}
  open CMD,'-|',"samtools faidx $xrq_ref $chr:$start-$end" or die $@;
  my $line;
  my $seq;
  while (defined($line=<CMD>)) {
    chomp($line);
    if ($line =~ m/>/){next;}
    $seq .= $line;
  }
  my $fastq_entry;
  if($first_line){
    $fastq_entry .= "\n";
  }else{
    $first_line++;
  }
  $fastq_entry .= "\@$chr:$pos:$start:$end\n";
  $fastq_entry .= "$seq\n";
  $fastq_entry .= "+\n";
  my $qual_length = ((2*$bases_surrounding)+1);
  foreach my $i (1..$qual_length){
    $fastq_entry .= "H";
  }
  print $fastq "$fastq_entry";
  if ($counter % 10000 == 0){print STDERR "Processed $chr\t$pos\n"}
  $pm->finish
}
$pm->wait_all_children;
close $fastq;
close INPUT;

print STDERR "Aligning SNP positions using BWA...\n";
#Align to HA412 and pull out SNP position from mapping.
system("$bwa mem -t $ncores_bwa $ha412_ref $tmp_prefix.fastq > $tmp_prefix.sam");
open SAM, "$tmp_prefix.sam";
$counter=0;
print STDERR "Pulling SNP positions from SAM file...\n";
my %loc_hash;
while(<SAM>){
  chomp;
  if ($_ =~ m/^@/){next;}
  $counter++;
  my @a = split(/\t/,$_);
  my @info = split(/:/,$a[0]);
  my $orig_chr = $info[0];
  my $orig_snp = $info[1];
  my $orig_start = $info[2];
  my $orig_end = $info[3];
  my $flag = $a[1];
  my $target_chr = $a[2];
  my $target_pos = $a[3];
  my $mq = $a[4];
  my $code = $a[5];
  if ($flag >= 2048){next;} #Skip secondary alignments
  if ($mq < $min_mq){next;} #Skip alignments with too poor mapping quality.
  #Try to figure out if there are soft clipping at the start for determining exact alignment location
  my @softclip = split(/S/, $code);
  my $clip_value = 0;
  if ($softclip[0] =~  /^[0-9]+$/ ){
    $clip_value = $softclip[0];
  }
  my $target_snp = $target_pos + $bases_surrounding - $clip_value;
  $loc_hash{$orig_chr.$orig_snp}{"chr"} = $target_chr;
  $loc_hash{$orig_chr.$orig_snp}{"pos"} = $target_snp;
  if ($counter % 10000 == 0){print STDERR "Processed sam $orig_chr\t$orig_snp\n"}
}

close SAM;
print STDERR "Outputting remapped VCF file...\n";
open(INPUT, "gunzip -c $input_file |");
open my $unsorted_vcf, '>', "$tmp_prefix.unsorted.vcf";
my @unaligned_array;
my %alignment_stats;
my %used_sites;
$counter=0;
my $first_contig_line;
my $printed_contig_lines;
while(<INPUT>){
  chomp;
  $counter++;
  if ($. == 1){print $unsorted_vcf "$_";next;}
  if ($_ =~ m/^##contig/){
    unless($first_contig_line){
      print $unsorted_vcf "\n##INFO=<ID=XRQ,Number=1,Type=String,Description='Bdistachyon_314_v3.0 genome location'";
      foreach my $contig (@contig_list){
        print $unsorted_vcf "\n##contig=<ID=$contig>";
      }
      print $unsorted_vcf "\n##source=Remapped to $ha412_ref using $bwa";
      $printed_contig_lines=1;
      $first_contig_line++;
    }
    next;
  }
  if ($_ =~ m/^#CHROM/){
    unless ($printed_contig_lines){
      print $unsorted_vcf "\n##INFO=<ID=Bdistachyon_314_v3.0,Number=1,Type=String,Description='Bdistachyon_314_v3.0 genome location'";
      foreach my $contig (@contig_list){
        print $unsorted_vcf "\n##contig=<ID=$contig>";
      }
      print $unsorted_vcf "\n##source=Remapped to $ha412_ref using $bwa";
      $printed_contig_lines=1;
    }
  }
  if ($_ =~ m/^#/){
    print $unsorted_vcf "\n$_";
    next;
  }
  
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  if ($loc_hash{$chr.$pos}{"chr"}){
    until (!$used_sites{$loc_hash{$chr.$pos}{'chr'}.$loc_hash{$chr.$pos}{'pos'}}){ #If the specific location was already used, shift the position by one.
       $loc_hash{$chr.$pos}{'pos'}++;
     }
     print $unsorted_vcf "\n$loc_hash{$chr.$pos}{'chr'}\t$loc_hash{$chr.$pos}{'pos'}";
     foreach my $i (2..6){
       print $unsorted_vcf "\t$a[$i]";
     }
     if ($a[7] ne '.'){
       print $unsorted_vcf "\t$a[7];";
     }else{
       print $unsorted_vcf "\t";
     } #Don't print info if its only .
     print $unsorted_vcf "XRQ=$chr.$pos";
     foreach my $i (8..$#a){
       print $unsorted_vcf "\t$a[$i]";
     }
     $alignment_stats{'aligned'}++;
     $used_sites{$loc_hash{$chr.$pos}{'chr'}.$loc_hash{$chr.$pos}{'pos'}}++;
  }else{
     $alignment_stats{'unaligned'}++;
     push(@unaligned_array,"$chr\t$pos");
  }
  if ($counter % 10000 == 0){print STDERR "Processed remapped vcf $chr\t$pos\n"}
}
close INPUT;
close $unsorted_vcf;

print STDERR "Outputting remapping stats...\n";
open my $alignment_stats, '>', "$tmp_prefix.Staceiconversionstats.txt";
unless($alignment_stats{'unaligned'}){
  $alignment_stats{'unaligned'} = 0;
}
unless($alignment_stats{'aligned'}){
  $alignment_stats{'aligned'} = 0;
}
my $percent_aligned = 100 * ($alignment_stats{'aligned'} / ($alignment_stats{'aligned'} + $alignment_stats{'unaligned'}));

print $alignment_stats "Percent_snps_remapped=$percent_aligned\n";
print $alignment_stats "SNPS_remapped=$alignment_stats{'aligned'}\n";
print $alignment_stats "SNPS_not_remapped=$alignment_stats{'unaligned'}\n";
print $alignment_stats "Unmapped sites";
foreach my $site (@unaligned_array){
  print $alignment_stats "\n$site";
}
my $sorted_vcf = "$tmp_prefix.remappedStacei.vcf.gz";
#Sort using BCFtools
system ("mkdir ./tmp_$tmp_prefix");

system("$bcftools sort -T ./tmp_$tmp_prefix $tmp_prefix.unsorted.vcf -O z > $sorted_vcf");
#Clean up temporary files
#system("rm $tmp_prefix.unsorted.vcf");
system("rm $tmp_prefix.fastq");
system("rm $tmp_prefix.sam"); 

