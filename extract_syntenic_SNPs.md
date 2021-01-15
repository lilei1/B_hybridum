# Objective: Pilar wants to extract the syntenic SNPs for constructing the phylogenetic tree.
  Josh developped some codes for filtering syntenic SNPs, and the link can be see here: https://github.com/jlevy44/JoshuaTree2/
  But it is hard for me to figure out. I debuged for a long while but could not make it work. Pilar need the file hurry, so I have to find other way to extract the syntenic SNPs.
  
  Here is the general idea about how I process the data: Joel Martin and I mapped the reads from distachyon, stacei, and hybridum samples back to Bd-21 reference, ABR114 reference, and the fake parents' reference (concatenate the Bd-21 reference and ABR114 reference) respectively, then called the variants with GATK. Then I've applied the same approach as (Loren Rieseberg's nature paper)[https://www.nature.com/articles/s41586-020-2467-6] and convert all of the SNPs relative to the distachyon reference/D-subgenome to the sylvaticum coordinate. By the way, for hybridum samples, I split them into two files of S-subgenome and D-subgenome.

The approach is that for each SNPs for D-subgenome or distachyon samples, we extract the 200bp context sequence surrounding the SNPs from the distachoyn reference genome. Then considering those 201 bp short sequences as reads, mapping them back to the stacei reference with bwa mem, then inferring the position for the target SNPs in the sylvaticum coordinate, and rewrite them into a single vcf file. Then I filtered out the discordant SNPs and merge all of the vcf files together.

Since John Lovell ran his pipeline GENESPACE for sylvaticum and distachyon for us and he figured out the syntenic blocks, I used his syntenic blocks data to do filtering and only keep SNPs in the syntenic blocks for the merged vcf file. Lastly, I convert the filtered vcf file into nexus file.


### File list:

```
#B. hybridum sample:
/global/homes/j/jieguo/Projects/reseq/Brachypodium_distachyon_redux_Bd21-3/set_hybridum/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz

```

### Step1: bcftools to split the hybridum samples into diatachyon and stacei part

```
#Only extract the B.distachyon samples
bcftools view /global/homes/j/jieguo/Projects/reseq/Brachypodium_distachyon_redux_Bd21-3/set_hybridum/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --regions Bdist.Bd1,Bdist.Bd2,Bdist.Bd3,Bdist.Bd4,Bdist.Bd5,Bdist.scaffold_12,Bdist.scaffold_14,Bdist.scaffold_135,Bdist.scaffold_180,Bdist.Bd1_centromere_containing_Bradi1g41430 -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/split_vcfs/Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz 

#extract the stacei samples from the hybridum file
sbatch -C haswell bcftools_split_stacei.sub

```