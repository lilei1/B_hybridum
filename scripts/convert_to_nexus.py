#!/usr/bin/env python
""" Convert the combined file into nexus formaat """
# by Li Lei with mentoring of Will, 10/04/2021
import argparse

def non_sample_headers (input_file_fh):
    input_file_fh.seek(0)
    headers = input_file_fh.readline().rstrip().split("\t")
    return headers[:6]

def get_merge_data(input_file_fh):
    input_file_fh.seek(0)
    headers = input_file_fh.readline().rstrip().split("\t")
    input_file_fh.seek(0)
    out = []
    for line in input_file_fh.readlines()[1:]:
        out.append({k:v for k, v in zip(headers, line.rstrip().split("\t"))})
    return out

# count the row of the file
def count_samples(merge_data):
    return len(merge_data[0])-6

#
def count_char(merge_fh,merge_data):
    sample_name = list(set(merge_data[0].keys()) - set(non_sample_headers(merge_fh)))[0]
    return len(get_sample_genotype(sample_name,merge_data))

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AG': 'R', 'CT': 'Y', 'CG': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M', 'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'N': 'N'}
def get_max_genotype_len(poly_status):
    poly_list = poly_status.split('/')
    alt_list = poly_list[1].split(',')
    possible_genotypes = alt_list + [poly_list[0]]
    return max(len(x) for x in possible_genotypes)

def get_haploid(genotype):
    if len(genotype) == 2:#only for SNPs
        if genotype[0] == genotype[1]:#homozygotes
            return genotype[0]
        return IUPAC[''.join(sorted(genotype))]
    index = int(len(genotype)/2)
    if genotype[:index] == genotype[index:]:
        return genotype[:index]
    return "?"

def get_haploid_with_gaps(genotype,poly_status1,poly_status2):
    haploid = get_haploid(genotype)
    max_len1 = get_max_genotype_len(poly_status1)
    max_len2 = get_max_genotype_len(poly_status2)
    max_len = max(max_len1,max_len2)
    return haploid + ("-" * (max_len - len(haploid)))

    #if #A/T,G,C
    #if all(len(x) == 1 for x in alt_list):#return true
    #    genotype
def get_sample_genotype(sample_name,merge_data): 
    return "".join([get_haploid_with_gaps(row[sample_name],row['Bd_SNPs'],row['St_SNPs']) for row in merge_data])
    
def get_all_genotypes(input_file_fh,output_fh):
    skip_headers = non_sample_headers(input_file_fh)
    merge_data = get_merge_data(input_file_fh)
    for sample_name in merge_data[0].keys():
        if sample_name in skip_headers:
            continue
        genotype = get_sample_genotype(sample_name,merge_data)
        print (f"\t{sample_name}\t{genotype}",file=output_fh)

def write_nexus_file(merge_fh,output_fh):
    merge_data = get_merge_data(merge_fh)
    print("#NEXUS",file=output_fh)
    print("",file=output_fh)
    print("BEGIN DATA;",file=output_fh)
    num_tax = count_samples(merge_data)
    num_char = count_char(merge_fh,merge_data)
    print(f"\tDIMENSIONS NTAX={num_tax} NCHAR={num_char};",file=output_fh)
    print("\tFORMAT DATATYPE=DNA MISSING=? GAP=- ;",file=output_fh)
    print("\tMATRIX",file=output_fh)
    get_all_genotypes(merge_fh,output_fh)
    print("",file=output_fh)
    print(";",file=output_fh)
    print("END;",file=output_fh)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert the nexus format")
    parser.add_argument(
        "merge_fh",
        type=argparse.FileType("r"),
        help="the merged file from extract_syntenic_snp.py",
    )
    parser.add_argument("output_fh", type=argparse.FileType("w"), help="The name of the nexus file")
    args = parser.parse_args()
    write_nexus_file(args.merge_fh,args.output_fh)