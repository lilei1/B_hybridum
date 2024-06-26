# Objective: Pilar wants to extract the syntenic SNPs for constructing the phylogenetic tree.
I used a way (see [extract_syntenic_SNPs.md](https://github.com/lilei1/B_hybridum/blob/main/extract_syntenic_SNPs.md)) to do the syntenic SNP extraction but there are  some bugs. So I discussed with John and we came up a new idea (See [example picture](https://github.com/lilei1/B_hybridum/blob/main/pictures/example.png))!!!


### Step1: bcftools to split the hybridum samples into distachyon and stacei part

#Only extract the B.distachyon chromosomes

```
bcftools view /global/projectb/sandbox/reseq/projects/working/Brachypodium_distachyon_redux_Bd21-3/set_hybridum_more/genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.anno.vcf.gz --regions Bdist.Bd1,Bdist.Bd2,Bdist.Bd3,Bdist.Bd4,Bdist.Bd5,Bdist.scaffold_12,Bdist.scaffold_14,Bdist.scaffold_135,Bdist.scaffold_180,Bdist.Bd1_centromere_containing_Bradi1g41430 -O z -o  /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/raw_Bdist.vcf.gz

```

#Only extract the B.stacei chormosomes with [bcftools_split_stacei.sub](https://github.com/lilei1/B_hybridum/blob/main/jobs/bcftools_split_stacei.sub)!!

```
module load esslurm

sbatch -C skylake -A plant -q jgi_exvivo bcftools_split_stacei.sub
```

So now the raw stacei/S vcf file is:

`/global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/raw_Bstacei.vcf.gz`

the raw distachyon/D vcf file is:

` /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/raw_Bdist.vcf.gz`

### Step2: Filter out the eight stacei samples from the raw_Bdist.vcf.gz file 

BG19, BG26, BG31, Bsta5, IBD10, LP6.1, TE4.3, ABR114

The removed stacei sample file are stored [here](https://github.com/lilei1/B_hybridum/blob/main/files/removed_stacei_from_D.txt)
```
vcftools --gzvcf /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/raw_Bdist.vcf.gz --remove /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/removed_stacei_from_D.txt  --recode --recode-INFO-all --out purify_Bdist
```

### Step3: Filter out the distachyon samples from raw _Bstacei.vcf.gz

The removed distachyon sample file are stored [here](https://github.com/lilei1/B_hybridum/blob/main/files/removed_distachyon_from_S.txt)

```
vcftools --gzvcf /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/raw_Bstacei.vcf.gz --remove /global/u2/l/llei2019/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way/removed_distachyon_from_S.txt  --recode --recode-INFO-all --out purify_Bstacei
```

#Change the header and chromosome ID:

##replace the chromosome id  in the hybridum-stacei subgenome

```
sed 's/Bstacei.//g' purify_Bstacei.recode.vcf >revised_purify_Bstacei.recode.vcf

sed 's/Bdist.//g' purify_Bdist.recode.vcf >revised_purify_Bdist.recode.vcf

```

#Then I need to change the genotypes with S_

```
grep "#" revised_purify_Bstacei.recode.vcf >head_revised_purify_Bstacei.recode.vcf
grep -v "#" revised_purify_Bstacei.recode.vcf >Nohead_revised_purify_Bstacei.recode.vcf
```

```
vi head_revised_Bstacei_purify_hybridum.recode.vcf

#add S_ to the end of each accessions of hybridum
```

```
cat head_revised_purify_Bstacei.recode.vcf Nohead_revised_purify_Bstacei.recode.vcf >fix_header_revised_purify_Bstacei.recode.vcf 
```

#Handle distachyon sample

I need to change the genotypes with D_

```
grep "#" revised_purify_Bdist.recode.vcf >head_revised_purify_Bdist.recode.vcf

vi head_revised_Bdist_BH_hybridum.vcf.gz
add D_ to the end of the 

#extract the non head

grep -v "#" revised_purify_Bdist.recode.vcf  >nohead_revised_purify_Bdist.recode.vcf
 
#recover the vcf file
cat head_revised_purify_Bdist.recode.vcf nohead_revised_purify_Bdist.recode.vcf>fix_head_revised_purify_Bdist.recode.vcf

```


### Step 4: zip and index the vcf file

```
bgzip fix_header_revised_Bstacei_purify_hybridum.recode.vcf

bcftools index fix_header_revised_Bstacei_purify_hybridum.recode.vcf.gz

bgzip fix_header_revised_Bdist_purify_hybridum.vcf

bcftools index fix_header_revised_Bdist_purify_hybridum.vcf.gz


```

### Step5: Filtering

Since we think we need to do filtering and only extract the SNPs and 'PASS' one, so I can run vcf tools and only extract the PASS and the SNPs. 

Below description define the "PASS" with QUAL>=40  and INFO/QD >=5 according to JGI's SNP calling pipeline.

```
##FILTER=<ID=Filter1,Description="Set if true: QUAL<40 || FILTER!=\"PASS\" || INFO/QD <5">
##FILTER=<ID=SnpGap,Description="SNP within 10 bp of an indel">
##FILTER=<ID=IndelGap,Description="Indel within 3 bp of an indel">
##bcftools_filterVersion=1.9-181-g55b420d+htslib-1.9-162-g6e6cdfc-dirty
##bcftools_filterCommand=filter -m + -O z -g 10 -G 3 '-eQUAL<40 || FILTER!="PASS" || INFO/QD <5' -s + -o  genotype_gvcfs.f1.bf=g10-G3-Q40-QD5.vcf.gz genotype_gvcfs.filtered.vcf.gz; Date=Wed Sep 22 10:50:16 2021
````

#Extract the "PASS" one:
```
bcftools view -f PASS fix_header_revised_purify_Bstacei.recode.vcf > PASS_fix_header_revised_purify_Bstacei.vcf
bcftools view -f PASS fixed_header_revised_purify_Bdist.recode.vcf > PASS_fixed_header_revised_purify_Bdist.vcf
```

#Only extract the SNPs

```
vcftools --vcf PASS_fix_header_revised_purify_Bstacei.vcf --remove-indels  --recode --recode-INFO-all --out onlySNP_PASS_fix_header_revised_purify_Bstacei

vcftools --vcf PASS_fixed_header_revised_purify_Bdist.vcf --remove-indels  --recode --recode-INFO-all --out onlySNP_PASS_fixed_header_revised_purify_Bdist
```

### Step 6: Convert the vcf to hmp format 

Convert the vcf to hmp format with [run_vcf2hapmap_stacei.job](https://github.com/lilei1/B_hybridum/blob/main/jobs/run_vcf2hapmap_stacei.job) and [run_vcf2hapmap_dist.job](https://github.com/lilei1/B_hybridum/blob/main/jobs/run_vcf2hapmap_dist.job)

```

module load esslurm

sbatch -C skylake -A plant -q jgi_exvivo run_vcf2hapmap_stacei.job
Submitted batch job 2284053
(base) 

sbatch -C skylake -A plant -q jgi_exvivo run_vcf2hapmap_dist.job
Submitted batch job 2284260
(base) 

```

### Step 7: Extract the syntenic SNPs

Then extract the syntenic SNPs with [run_extract_syntenic_SNP.job](https://github.com/lilei1/B_hybridum/blob/main/jobs/run_extract_syntenic_SNP.job)

The python script run_extract_syntenic_SNP.job is [here](https://github.com/lilei1/B_hybridum/blob/main/jobs/run_extract_syntenic_SNP.job)
The lookup table of the position is from Ruben.

```
sbatch -C skylake -A plant -q jgi_exvivo run_extract_syntenic_SNP.job

Submitted batch job 2361256

#Check how many SNPs we have got:

llei2019@cori16 10:21:09 ~/cscratch/B_hybridum/Ruben_synteny/fake_parents_ref_way 
$ wc -l combined_D_dist_S_stacei.hmp.txt

 wc -l combined_D_dist_S_stacei.hmp.txt
52520013 combined_D_dist_S_stacei.hmp.txt

52,520,013 in total 
```

### Step 8: Convert to nexus format:

Then I need to think about convert the matrix into the nexus format with [run_convert_nexus.job](https://github.com/lilei1/B_hybridum/blob/main/jobs/run_convert_nexus.job)!!! The python script [convert_to_nexus.py](https://github.com/lilei1/B_hybridum/blob/main/scripts/convert_to_nexus.py)

```
sbatch -C skylake -A plant -q jgi_exvivo run_convert_nexus.job
```


