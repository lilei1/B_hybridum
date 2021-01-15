# Objective: Pilar wants to extract the syntenic SNPs for constructing the phylogenetic tree.
  Josh developped some codes for filtering syntenic SNPs, and the link can be see here: https://github.com/jlevy44/JoshuaTree2/
  But it is hard for me to figure out. I debuged for a long while but could not make it work. Pilar need the file hurry, so I have to find other way to extract the syntenic SNPs.
  
  Here is the general idea about how I process the data: Joel Martin and I mapped the reads from distachyon, stacei, and hybridum samples back to Bd-21 reference, ABR114 reference, and the fake parents' reference (concatenate the Bd-21 reference and ABR114 reference) respectively, then called the variants with GATK. Then I've applied the same approach as (Loren Rieseberg's nature paper)[https://www.nature.com/articles/s41586-020-2467-6] and convert all of the SNPs relative to the distachyon reference/D-subgenome to the sylvaticum coordinate. By the way, for hybridum samples, I split them into two files of S-subgenome and D-subgenome.

The approach is that for each SNPs for D-subgenome or distachyon samples, we extract the 200bp context sequence surrounding the SNPs from the distachoyn reference genome. Then considering those 201 bp short sequences as reads, mapping them back to the stacei reference with bwa mem, then inferring the position for the target SNPs in the sylvaticum coordinate, and rewrite them into a single vcf file. Then I filtered out the discordant SNPs and merge all of the vcf files together.

Since John Lovell ran his pipeline GENESPACE for sylvaticum and distachyon for us and he figured out the syntenic blocks, I used his syntenic blocks data to do filtering and only keep SNPs in the syntenic blocks for the merged vcf file. Lastly, I convert the filtered vcf file into nexus file.


### File list:

```
#B. hybridum samples

/global/homes/j/jieguo/Projects/reseq/Brachypodium_distachyon_redux_Bd21-3/set_hybridum/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz

#B.distachyon samples from pan genome project

/global/projectb/scratch/j_martin/brachy-thing/distachyon/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz 

#B. stacei samples

/global/projectb/sandbox/reseq/projects/working/Brachypodium_distachyon_redux_Bd21-3/set_stacei/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz

#distachyon reference

/global/dna/projectdirs/RD/reseq_store/genomes/plants/Brachypodium_distachyon/versions/3.1/Bdistachyon_314_v3.0.fa

#stacei reference

/global/dna/projectdirs/RD/reseq_store/genomes/plants/Brachypodium_stacei/versions/1.0/Bstacei_316_v1.0.fa


```

### Step1: bcftools to split the hybridum samples into diatachyon and stacei part

```
#Only extract the B.distachyon samples

bcftools view /global/homes/j/jieguo/Projects/reseq/Brachypodium_distachyon_redux_Bd21-3/set_hybridum/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --regions Bdist.Bd1,Bdist.Bd2,Bdist.Bd3,Bdist.Bd4,Bdist.Bd5,Bdist.scaffold_12,Bdist.scaffold_14,Bdist.scaffold_135,Bdist.scaffold_180,Bdist.Bd1_centromere_containing_Bradi1g41430 -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/split_vcfs/Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz 

#Change the header and chromosome ID:

##replace the chormosome id  the hybridum-diatachyon subgenome

zcat /global/u2/l/llei2019/cscratch/B_hybridum/split_vcfs/Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz | sed 's/Bdist.//g' >revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

#add the D_ to each genotype ID
grep "#" revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >head_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf 

vi head_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

add D_

#extract the non head
grep -v "#" revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >nohead_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

#recover the vcf file

cat head_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf nohead_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >biSNPs_only_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz


#extract the stacei samples from the hybridum file

sbatch -C haswell bcftools_split_stacei.sub

#Change the header and chromosome ID:

##replace the chormosome id  the hybridum-stacei subgenome

zcat /global/u2/l/llei2019/cscratch/B_hybridum/split_vcfs/Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz| sed 's/Bstacei.//g' >revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

#Then I need to change the genotypes with S_

grep "#" revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >head_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf 

vi revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >head_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

add S_

#extract the non head
grep -v "#" revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >nohead_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf 

#recover the vcf file

cat head_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf nohead_revised_Bstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf >revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf


```

### Step2: only keep biallelic SNPs for the converted vcf file with vcftools

```
#For distachyon vcf file

vcftools --gzvcf /global/projectb/scratch/j_martin/brachy-thing/distachyon/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /global/cscratch1/sd/llei2019/B_hybridum/remapped_stacei/biSNPs_only_distachyon_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

#zip the file:

gzip 'biSNPs_only_distachyon_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf'

#For distachyon sub-genome vcf file

vcftools --gzvcf /global/projectb/scratch/j_martin/brachy-thing/distachyon/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /global/cscratch1/sd/llei2019/B_hybridum/remapped_stacei/biSNPs_only_distachyon_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

gzip biSNPs_only_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz



#For stacei vcf file

vcftools --gzvcf /global/projectb/sandbox/reseq/projects/working/Brachypodium_distachyon_redux_Bd21-3/set_stacei/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /global/cscratch1/sd/llei2019/B_hybridum/remapped_stacei/biSNPs_Bstacei_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz

#For stacei sub-genome vcf file

vcftools --vcf /global/cscratch1/sd/llei2019/B_hybridum/split_vcfs/revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /global/cscratch1/sd/llei2019/B_hybridum/remapped_stacei/biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf

```

### Step3: Convert the distachyon reference to stacei reference coordinate:

```
#install the required software:

conda info -e

source activate /global/homes/l/llei2019/bscratch/software/my_work_en

conda install -c bioconda perl-parallel-forkmanager
 
#install bwa

conda install -c bioconda bwa

#Convert the distachyon samples with adapted (convert_reference.pl)[https://github.com/owensgl/wild_gwas_2018/blob/master/xrqpos2ha412pos_bwa.pl]

perl ~/Github/B_hybridum/scripts/convert_reference.pl ~/cscratch/B_hybridum/remapped_stacei/biSNPs_only_distachyon_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz Bdist

#Convert the distachyon sub-genome samples with adapted (convert_reference.pl)[https://github.com/owensgl/wild_gwas_2018/blob/master/xrqpos2ha412pos_bwa.pl]

perl ~/Github/B_hybridum/scripts/convert_reference.pl ~/cscratch/B_hybridum/remapped_stacei/biSNPs_only_revised_Bdist_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz Bdist_D_subgeno

```
### Step4: Since I tried to merge all of the processed files but with some discordance errors:

#Here is the file list

```
Distachyon_convert: Bdist.remappedStacei.distachyon_geno.vcf.gz
bcftools index Bdist.remappedStacei.distachyon_geno.vcf.gz

Distachyon subgenome: Bdist_D_subgeno.remappedStacei.vcf.gz
bcftools index Bdist_D_subgeno.remappedStacei.vcf.gz

Stacei subgenome: biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz

#Heads up: run bgzip need to use the conda enviroment "/global/homes/l/llei2019/bscratch/software/my_work_en"

bgzip biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf

bcftools index biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz

Stecei genome: biSNPs_Bstacei_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz.recode.vcf

```

#merge all of the vcf files using submit jobs!!!

```
sbatch -C haswell /global/projectb/scratch/llei2019/jobs/bcf_merge_hybridum.sub

```

##I first need figure out how many of those sites existing in the file!!!!
same position but different REF/ALT!!!!

Chr01   695     .       C       T       180.93  SnpGap
Chr01   695     .       G       C       601.92  SnpGap

### Step 5: Run vcftools for the files and figure out the concordance sites:

```
#figure put the distachyon and stacei problem
sbatch -C haswell ~/bscratch/jobs/run_vcftools_diff.site.sub 

#extract the condordant SNPs:

grep "B" Bdist_D_subgeno_remapped_vs_biSNPs_Bstacei_genotype.diff.sites_in_files|wc -l
96081

grep "B" Bdist_D_subgeno_remapped_vs_biSNPs_Bstacei_genotype.diff.sites_in_files|cut -f 1,2 >concordance_SNPs.list

#figure out he distachyon sub-genome and stacei problem

sbatch -C haswell ~/bscratch/jobs/run_vcftools_diff.site.sub 

#extract the condordant SNPs:
grep "B" Bdist.remapped_vs_biSNPs_Bstacei_genotype.diff.sites_in_files|cut -f 1,2>Bdist.remapped_vs_biSNPs_Bstacei.txt

##get the intersection of the two files!!!!!

awk 'NR==FNR { lines[$0]=1; next } $0 in lines' Bdist.remapped_vs_biSNPs_Bstacei.txt concordance_SNPs.list >common_concordance_SNPs.list

#create the bed file for the nonsyn SNPs

awk -F '\t' -v OFS='\t' '{ $(NF+1) = $2-1; print }' common_concordance_SNPs.list >common_concordance_SNPs.list.txt

paste  <(cut -f 1 common_concordance_SNPs.list.txt) <(cut -f 3 common_concordance_SNPs.list.txt) <(cut -f 2 common_concordance_SNPs.list.txt) >common_concordance_SNPs.list.bed

```

### Step 6: Only keep the concordant SNPs for all of the files with bcftool

```
bcftools view --regions-file /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs.list.bed /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/Bdist_D_subgeno.remappedStacei.vcf.gz -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs_Bdist_D_subgeno.remappedStacei.vcf.gz

bcftools view --regions-file /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs.list.bed /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/biSNPs_Bstacei_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz.recode.vcf.gz -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs_biSNPs_Bstacei_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz.recode.vcf.gz


bcftools view --regions-file /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs.list.bed /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/Bdist.remappedStacei.distachyon_geno.vcf.gz -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf.gz

bcftools view --regions-file /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs.list.bed /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs_biSNPs_revised_only_forBstacei_BH_genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.recode.vcf.gz


```

### Step 7: Merge all of the vcf files using submit jobs!!!

```
sbatch -C haswell ~/bscratch/jobs/bcf_merge_hybridum.sub 

```

### STEP 8: Cut the syntenic block from John Lovell's results and extract the SNPs in the syntenic region:

```
bcftools index common_concordance_SNPs_hybridum_stacei_distachyon.vcf.gz

bcftools view --regions-file /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/syntenic_block_Lovell.bed /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/common_concordance_SNPs_hybridum_stacei_distachyon.vcf.gz -O z -o /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/syntenic_common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf.gz

# count the SNPs
zgrep -v "#" syntenic_common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf.gz|wc -l
50013

#since I found there are stacei samples also mapping back to the fake parents, so I have to filter out the stacei samples from the hybridum file

vcftools --vcf /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/re_syntenic_common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf --remove /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/remove_indi.txt --recode --recode-INFO-all --out /global/u2/l/llei2019/cscratch/B_hybridum/remapped_stacei/filter_staceiInhyb_re_syntenic_common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf

```
### Step9 : convert the vcf file into the nexus format with (vcf2phylip.py)[https://github.com/edgardomortiz/vcf2phylip]

```
cd ~/bscratch/software/

git clone https://github.com/edgardomortiz/vcf2phylip.git

cd vcf2phylip/

python ~/bscratch/software/vcf2phylip/vcf2phylip.py --input filter_staceiInhyb_re_syntenic_common_concordance_SNPs_Bdist.remappedStacei.distachyon_geno.vcf.recode.vcf --phylip-disable --nexus 

```
Both final vcf and nexus files can be avaible in the google drive: https://drive.google.com/drive/u/1/folders/1S_ur1RGGQ_fBj8ZNnjsiGf0_sqTbpRLp
