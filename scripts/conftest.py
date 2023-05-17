"""Define fixtures to be used to setup tests"""

from collections import UserList

import pytest


@pytest.fixture(name="hmp1_fh")
def fixture_hmp1_fh():
    fake_fh = UserList(
        [
            "rs# SNPalleles #CHROM POS Strand Genome_build Center ProtLSID assayLSID panelLSID QC_code Bhyb123-1\n",
            "rs5003416 T/C Bd2 58660304 4330.81 PASS unknown urn:unknown urn:unknown urn:unknown QC+ TT\n",
        ]
    )
    setattr(fake_fh, "name", "fake.tsv")
    setattr(fake_fh, "seek", lambda x: None)
    setattr(fake_fh, "readline", lambda: fake_fh[0])
    return fake_fh


@pytest.fixture(name="hmp2_fh")
def fixture_hmp2_fh():
    fake_fh = UserList(
        [
            "rs# SNPalleles #CHROM POS Strand Genome_build Center ProtLSID assayLSID panelLSID QC_code S288C\n",
            "rs1 T/C Chr01 360379 4330.81 PASS unknown urn:unknown urn:unknown urn:unknown QC+ TT\n",
        ]
    )
    setattr(fake_fh, "name", "fake.tsv")
    setattr(fake_fh, "seek", lambda x: None)
    setattr(fake_fh, "readline", lambda: fake_fh[0])
    return fake_fh


@pytest.fixture(name="lookup_fh")
def fixture_lookup_fh():
    fake_fh = UserList(["Bd2\t58660304\tT\tforward\tChr01\t360379\tT\t1\n"])
    setattr(fake_fh, "name", "fake.tsv")
    setattr(fake_fh, "seek", lambda x: None)
    setattr(fake_fh, "readline", lambda: fake_fh[0])
    return fake_fh

@pytest.fixture(name="merge_fh")
def fixture_merge_fh():
    fake_fh = UserList(["Bd_SNPs\tBd_chr\tBd_pos\tSt_SNPs\tSt_chr\tSt_pos\tRon2\tAL1A\tAL2D\tAL2E\tAL2F\t2B\n",
                        "T/T\tBd2\t58660303\tC/C\tChr01\t360378\tTT\tTT\tTT\tTT\tTT\tTT\n" 
    ])
    setattr(fake_fh, "name", "fake.tsv")
    setattr(fake_fh, "seek", lambda x: None)
    setattr(fake_fh, "readline", lambda: fake_fh[0])
    return fake_fh
