#!/usr/bin/env python
""" extract the syntenic SNPs based on the lookup table given by Ruben """
# by Li Lei with mentoring of Will, 8/23/2021
import argparse
import logging
import os
import sqlite3
import tempfile


logger = logging.getLogger(__name__)


def is_passing(fields):
    """Only keep the variants with the sixth column equal to PASS"""
    return fields[5] == "PASS"


def store(fields, storage):
    location = ",".join([fields[2], fields[3]])
    line = "\t".join(fields)
    cur = storage.cursor()  # define the transaction
    cur.execute("insert into hmp (location,line) values (?,?)", (location, line))
    storage.commit()  # end of the transaction


def initilize_storage():
    logger.debug("Initalizing database")
    temp_dir = getattr(os.environ, "CSCRATCH", "/tmp")
    connection = sqlite3.connect(tempfile.mkstemp(dir=temp_dir)[1])  # return the connection to the database
    cur = connection.cursor()  # define the transaction
    cur.execute("create table hmp (location,line)")  # create a table
    cur.execute("CREATE INDEX index_hmp_location ON hmp(location)")
    connection.commit()  # end of the transaction
    logger.debug("Database initialization was sucessfull")
    return connection  # allow us to send the connection to the databse later


def load_hmp(hmp, filter):
    storage = initilize_storage()
    logger.info("Starting to process HMP file %s", hmp.name)
    cur = storage.cursor()  # define the transaction
    for i, line in enumerate(hmp):
        logger.debug("From HMP file %s got line number %d containing %s", hmp.name, i, line)
        if i % 1000 == 0:
            logger.info("Processing line %d.", i)
        if i == 0:
            continue
        fields = line.rstrip().split("\t")
        if (filter and is_passing(fields)) or not filter:
            location = ",".join([fields[2], fields[3]])
            line = "\t".join(fields)
            cur.execute("insert into hmp (location,line) values (?,?)", (location, line))
            logger.debug(
                "From HMP file %s stored line number %d containing %s",
                hmp.name,
                i,
                line,
            )
    storage.commit()  # end of the transaction
    logger.info("Done processing HMP file %s", hmp.name)
    return storage


def retrieve(key_list, storage):
    assert isinstance(key_list, list)
    cur = storage.cursor()  # define the transaction
    result = list(
        cur.execute(
            f"select location, line from hmp where location IN ({','.join('?'*len(key_list))})", key_list
        )
    )
    storage.commit()  # end of the transaction
    out = {}
    for x in result:
        location = x[0]
        line = x[1]
        out[location] = line.split("\t")
    return out


def num_samples(storage):
    cur = storage.cursor()  # define the transaction
    result = list(cur.execute("select COUNT(*) from hmp"))
    storage.commit()  # end of the transaction
    return result[0][0]


def complementary(genotype):
    #cm_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    cm_dict = {"A": "A", "T": "T", "C": "C", "G": "G"}
    return "".join([cm_dict.get(x, x) for x in genotype])


def get_merged_header(hmp1_fh, hmp2_fh):
    hmp1_fh.seek(0)  # return the filehandle
    hmp2_fh.seek(0)
    hmp1_header = hmp1_fh.readline()
    hmp2_header = hmp2_fh.readline()
    hmp1_header_fileds = hmp1_header.rstrip().split("\t")
    hmp2_header_fileds = hmp2_header.rstrip().split("\t")
    assert len(hmp1_header_fileds[11:]) > 0
    assert len(hmp2_header_fileds[11:]) > 0
    return "\t".join(
        ["Bd_SNPs", "Bd_chr", "Bd_pos", "St_SNPs", "St_chr", "St_pos"]
        + hmp1_header_fileds[11:]
        + hmp2_header_fileds[11:]
    )

def get_acc_nb(hmp_fh):
    hmp_fh.seek(0)  # return the filehandle
    hmp_header = hmp_fh.readline()
    hmp_header_fileds = hmp_header.rstrip().split("\t")
    num_acc = len(hmp_header_fileds[11:])
    return num_acc

def merge(hmp1_fh, hmp2_fh, lookup_fh, merged_fh, filter):
    #breakpoint()
    hmp1 = load_hmp(hmp1_fh, filter)
    hmp2 = load_hmp(hmp2_fh, filter)
    num_acc1 = get_acc_nb(hmp1_fh)
    logger.info("HMP file %s has %d samples.", hmp1_fh.name, num_acc1)
    num_acc2 = get_acc_nb(hmp2_fh)
    logger.info("HMP file %s has %d samples.", hmp2_fh.name, num_acc2)
    print(get_merged_header(hmp1_fh, hmp2_fh), file=merged_fh)
    logger.info("Staring merge with lookup file %s.", lookup_fh.name)
    lines_to_process = []
    for i, line in enumerate(lookup_fh):
        if i == 0:
            continue
        if i % 1000 == 0:
            logger.info("Processing lookup table file line %d.", i)
        lines_to_process.append(line.split("\t"))
        if len(lines_to_process) > 10e3:
            merge_batch(hmp1, hmp2, num_acc1, num_acc2, lines_to_process, merged_fh)
            lines_to_process = []
    if len(lines_to_process) > 0:
        merge_batch(hmp1, hmp2, num_acc1, num_acc2,lines_to_process, merged_fh)


def merge_batch(hmp1, hmp2, num_acc1, num_acc2, lines_to_process, merged_fh):
    #breakpoint()
    hmp1_locations = []
    hmp2_locations = []
    for lookup_fields in lines_to_process:
        hmp1_locations.append(",".join(lookup_fields[0:2]))
        hmp2_locations.append(",".join(lookup_fields[4:6]))
    hmp1_dict = retrieve(hmp1_locations, hmp1)
    hmp2_dict = retrieve(hmp2_locations, hmp2)
    for lookup_fields in lines_to_process:
        transform = complementary if lookup_fields[3] == "reverse" else lambda x: x
        try:
            hmp1_fields = hmp1_dict[",".join(lookup_fields[0:2])]
        except KeyError:
            hmp1_fields = None
        try:
            hmp2_fields = hmp2_dict[",".join(lookup_fields[4:6])]
        except KeyError:
            hmp2_fields = None
        if hmp1_fields is None and hmp2_fields is None and lookup_fields[2] == transform(lookup_fields[6]):
            continue
        out_line = generate_output_line(hmp1_fields, num_acc1, hmp2_fields, num_acc2, lookup_fields)
        print(out_line, file=merged_fh)


def generate_output_line(hmp1_fields, num_acc1, hmp2_fields, num_acc2, lookup_fields):
    if hmp1_fields is None:
        hmp1_fields = (
            ["", lookup_fields[2] + "/" + lookup_fields[2]]
            + lookup_fields[0:2]
            + [""] * 7
            + [lookup_fields[2] + lookup_fields[2]] * num_acc1
        )
    if hmp2_fields is None:
        hmp2_fields = (
            ["", lookup_fields[6] + "/" + lookup_fields[6]]
            + lookup_fields[4:6]
            + [""] * 7
            + [lookup_fields[6] + lookup_fields[6]] * num_acc2
        )
    if lookup_fields[3] == "reverse":
        hmp2_fields = (
            [hmp2_fields[0]]
            + [complementary(hmp2_fields[1])]
            + hmp2_fields[2:11]
            + [complementary(x) for x in hmp2_fields[11:]]
        )
    out_fields = hmp1_fields[1:4] + hmp2_fields[1:4] + hmp1_fields[11:] + hmp2_fields[11:]
    assert len(out_fields) == 6 + num_acc1 + num_acc2
    return "\t".join(out_fields)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="merge hmp files")
    parser.add_argument(
        "hmp_file1",
        type=argparse.FileType("r"),
        help="The name of the hmp file1 to be merged",
    )
    parser.add_argument(
        "hmp_file2",
        type=argparse.FileType("r"),
        help="The name of the hmp file2 to be merged",
    )
    parser.add_argument("lookup_file", type=argparse.FileType("r"), help="The name of the lookup file")
    parser.add_argument("merged_file", type=argparse.FileType("w"), help="The name of the merged file")
    parser.add_argument(
        "-log",
        "--log",
        default="info",
        help=("Provide logging level. " "Example --log debug', default='info'"),
    )
    parser.add_argument("--filter", dest="filter", action="store_true", help="Filter out non-PASSing lines")
    parser.add_argument("--no-filter", dest="filter", action="store_false", help="Keep non-PASSing lines")
    parser.set_defaults(filter=True)

    args = parser.parse_args()
    levels = {
        "critical": logging.CRITICAL,
        "error": logging.ERROR,
        "warn": logging.WARNING,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    if args.log.lower() not in levels.keys():
        raise ValueError(f"log level given: {args.log}" f" -- must be one of: {' | '.join(levels.keys())}")
    logging.basicConfig(level=levels[args.log.lower()])
    merge(args.hmp_file1, args.hmp_file2, args.lookup_file, args.merged_file, args.filter)
