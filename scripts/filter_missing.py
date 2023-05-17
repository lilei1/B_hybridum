#!/usr/bin/env python
""" VCF Filtering """
# by Li Lei with mentoring of Will, 7/19/2021
import argparse


def exists_for_sample(genotype):
    """
    Returns True if the sample has defined genotype at both alleles
    Inputs:
        genotype: allelic status of each genotype
    """
    return genotype[0] != "." and genotype[2] != "."  # return true


def subset_passing(row, colnames, threshold):
    assert 0 <= threshold <= 1
    if len(colnames) == 0:
        raise ValueError("colnames cannot be an empty list")
    if threshold == 0:
        return True
    num_exists = 0
    for name in colnames:
        if exists_for_sample(row["SAMPLES"][name]):
            num_exists += 1
            if num_exists >= len(colnames) * threshold:
                return True
    return False


def missingness(row, groups, threshold):
    for sample_names in groups.values():
        if not subset_passing(row, sample_names, threshold):
            return False
    return True


def filter_vcf(file_handle, threshold, group_marker):
    for line in file_handle:
        line = line.rstrip()
        if line.startswith("##"):
            print(line)
            continue
        if line.startswith("#CHROM"):
            print(line)
            colnames = line.strip()[1:].split("\t")
            groups = {
                "D": colnames[colnames.index("FORMAT") + 1 : colnames.index(group_marker)],
                "S": colnames[colnames.index(group_marker) :],
            }
        else:
            tmp = line.strip().split("\t")
            row = {}
            samples = {}
            for i, name in enumerate(colnames):
                if i > colnames.index("FORMAT"):  # return the index of the colnames of "FORMAT"
                    samples[name] = tmp[i]
                else:
                    row[name] = tmp[i]
            row["SAMPLES"] = samples
            if missingness(row, groups, threshold):
                print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filtering of VCF files")
    parser.add_argument("vcf_file", type=argparse.FileType("r"), help="Name of VCF file to be filtered")
    parser.add_argument(
        "threshold",
        type=float,
        help="fraction of samples with defined genotype (in range 0 to 1)",
    )
    parser.add_argument(
        "group_marker",
        help="Name of left most column that defines start of second subset",
    )
    args = parser.parse_args()
    filter_vcf(args.vcf_file, args.threshold, args.group_marker)
