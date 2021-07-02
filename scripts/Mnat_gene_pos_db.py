#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Elias Dohmen"
__version__ = "0.1"
__email__ = "e.dohmen@wwu.de"
__institute__ = "IEB Münster"

import argparse
import re
import pandas as pd
import os
import sys
from collections import defaultdict


def gene_db_from_gff(gff):
    """
        reads the defined input gff and stores scaffold, gene name, start and stop position in the scaffold in a pandas dataframe
        :param gff: the input gff to read the information from
        :return: gene_db (pandas dataframe with genetic information)
        """

    d = {'contig': 'dummy', 'geneID': 'dummy', 'start': 1, 'end': 2, 'strand': '+'}
    gene_db = pd.DataFrame(data=d, index=[0])
    contig_lst = []
    geneID_lst = []
    start_lst = []
    end_lst = []
    strnd_lst = []


    with open(gff, 'r') as gin:
        for raw_line in gin:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()

                contig_lst.append(line[0])
                geneID_lst.append(re.search('ID=(\S+)',line[8]).group(1))
                start_lst.append(line[3])
                end_lst.append(line[4])
                strnd_lst.append(line[6])


    df1 = pd.DataFrame(
        {
            "contig": contig_lst,
            "geneID": geneID_lst,
            "start": start_lst,
            "end": end_lst,
            "strand": strnd_lst,
        },
    )

    gene_db = pd.concat([gene_db, df1])

    gene_db = gene_db[gene_db.contig != 'dummy']

    return gene_db


def main():
    """
    the main function of the script
    reads a gff and saves the created pandas dataframe to a csv file
    """

    gene_db = gene_db_from_gff("Macrotermes_natalensis_mRNA_sorted.gff")

    print("writing gene db...")
    gene_db.to_csv(r'gene_db2_pd.csv')

if __name__ == "__main__":
    main()