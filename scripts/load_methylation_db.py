#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Elias Dohmen"
__version__ = "0.1"
__email__ = "e.dohmen@uni-meunster.de"
__institute__ = "IEB Münster"

import argparse
import re
import pandas as pd
import os
import sys

def summarise(out_file, methylation_db):
    """
    creates a methylation summary and saves it to a file
    :param out_file: output file
    :param methylation_db: pandas db with methylation info
    :return: no return value
    """
    with open(out_file, 'w') as of:
        of.write("total sequenced CpGs: " + str(len(methylation_db)) + "\n")
        methylation_df = methylation_db.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs (contig+position duplicates removed): " + str(len(methylation_df)) + "\n")

        pval_df = methylation_db[methylation_db["p-val"] < 0.05]
        of.write("total methylated CpGs with P-Value < 0.05: " + str(len(pval_df)) + "\n")
        pval_dr = pval_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total methylated CpGs (contig+position duplicates removed) with P-Value < 0.05: " + str(len(pval_dr)) + "\n")

        # worker
        of.write("\nWorkers:\n")
        worker_df = methylation_db[methylation_db["caste"] == "worker"]
        of.write("total sequenced CpGs in workers: " + str(len(worker_df)) + "\n")
        worker_dr = worker_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in workers (duplicates removed): " + str(len(worker_dr)) + "\n")
        worker_pval = worker_df[worker_df["p-val"] < 0.05]
        of.write("total methylated CpGs in workers with P-Value < 0.05: " + str(len(worker_pval)) + "\n")
        worker_pval_dr = worker_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in workers with P-Value < 0.05 (duplicates removed): " + str(len(worker_pval_dr)) + "\n")

        # alate
        of.write("\nAlates:\n")
        alate_df = methylation_db[methylation_db["caste"] == "alate"]
        of.write("total sequenced CpGs in alates: " + str(len(alate_df)) + "\n")
        alate_dr = alate_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in alates (duplicates removed): " + str(len(alate_dr)) + "\n")
        alate_pval = alate_df[alate_df["p-val"] < 0.05]
        of.write("total methylated CpGs in alates with P-Value < 0.05: " + str(len(alate_pval)) + "\n")
        alate_pval_dr = alate_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in alates with P-Value < 0.05 (duplicates removed): " + str(len(alate_pval_dr)) + "\n")

        # queen
        of.write("\nQueens:\n")
        queen_df = methylation_db[methylation_db["caste"] == "queen"]
        of.write("total sequenced CpGs in queens: " + str(len(queen_df)) + "\n")
        queen_dr = queen_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in queens (duplicates removed): " + str(len(queen_dr)) + "\n")
        queen_pval = queen_df[queen_df["p-val"] < 0.05]
        of.write("total methylated CpGs in queens with P-Value < 0.05: " + str(len(queen_pval)) + "\n")
        queen_pval_dr = queen_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in queens with P-Value < 0.05 (duplicates removed): " + str(len(queen_pval_dr)) + "\n")

        # king
        of.write("\nKings:\n")
        king_df = methylation_db[methylation_db["caste"] == "king"]
        of.write("total sequenced CpGs in queens: " + str(len(queen_df)) + "\n")
        king_dr = king_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in kings (duplicates removed): " + str(len(king_dr)) + "\n")
        king_pval = king_df[king_df["p-val"] < 0.05]
        of.write("total methylated CpGs in kings with P-Value < 0.05: " + str(len(king_pval)) + "\n")
        king_pval_dr = king_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in kings with P-Value < 0.05 (duplicates removed): " + str(len(king_pval_dr)) + "\n")

        # CDS
        of.write("\nCDS:\n")
        cds_df = methylation_db[methylation_db["context"] == "cds"]
        of.write("total sequenced CpGs in CDS: " + str(len(cds_df)) + "\n")
        cds_dr = cds_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in CDS (duplicates removed): " + str(len(cds_dr)) + "\n")

        cds_pval = cds_df[cds_df["p-val"] < 0.05]
        of.write("total methylated CpGs in CDS with P-Value < 0.05: " + str(len(cds_pval)) + "\n")
        cds_pval_dr = cds_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in CDS with P-Value < 0.05 (duplicates removed): " + str(len(cds_pval_dr)) + "\n")

        # TEs
        of.write("\nTEs:\n")
        TE_df = methylation_db[methylation_db["context"] == "TE"]
        of.write("total sequenced CpGs in TEs: " + str(len(TE_df)) + "\n")
        TE_dr = TE_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in TEs (duplicates removed): " + str(len(TE_dr)) + "\n")

        TE_pval = TE_df[TE_df["p-val"] < 0.05]
        of.write("total methylated CpGs in TEs with P-Value < 0.05: " + str(len(TE_pval)) + "\n")
        TE_pval_dr = TE_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in TEs with P-Value < 0.05 (duplicates removed): " + str(len(TE_pval_dr)) + "\n")

        # introns
        of.write("\nIntrons:\n")
        intron_df = methylation_db[methylation_db["context"] == "intron"]
        of.write("total sequenced CpGs in introns: " + str(len(intron_df)) + "\n")
        intron_dr = intron_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in introns (duplicates removed): " + str(len(intron_dr)) + "\n")

        intron_pval = intron_df[intron_df["p-val"] < 0.05]
        of.write("total methylated CpGs in introns with P-Value < 0.05: " + str(len(intron_pval)) + "\n")
        intron_pval_dr = intron_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in introns with P-Value < 0.05 (duplicates removed): " + str(len(intron_pval_dr)) + "\n")

        # flanks
        of.write("\nFlanks:\n")
        flank_df = methylation_db[(methylation_db["context"] == "flank3") | (methylation_db["context"] == "flank5")]
        of.write("total sequenced CpGs in flanks: " + str(len(flank_df)) + "\n")
        flank_dr = flank_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flanks (duplicates removed): " + str(len(flank_dr)) + "\n")

        flank_pval = flank_df[flank_df["p-val"] < 0.05]
        of.write("total methylated CpGs in flanks with P-Value < 0.05: " + str(len(flank_pval)) + "\n")
        flank_pval_dr = flank_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flanks with P-Value < 0.05 (duplicates removed): " + str(len(flank_pval_dr)) + "\n")

        # flank3
        of.write("\nFlank3:\n")
        flank3_df = methylation_db[methylation_db["context"] == "flank3"]
        of.write("total sequenced CpGs in flank3: " + str(len(flank3_df)) + "\n")
        flank3_dr = flank3_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flank3 (duplicates removed): " + str(len(flank3_dr)) + "\n")

        flank3_pval = flank3_df[flank3_df["p-val"] < 0.05]
        of.write("total methylated CpGs in flank3 with P-Value < 0.05: " + str(len(flank3_pval)) + "\n")
        flank3_pval_dr = flank3_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flank3 with P-Value < 0.05 (duplicates removed): " + str(len(flank3_pval_dr)) + "\n")

        # flank5
        of.write("\nFlank5:\n")
        flank5_df = methylation_db[methylation_db["context"] == "flank5"]
        of.write("total sequenced CpGs in flank5: " + str(len(flank5_df)) + "\n")
        flank5_dr = flank5_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flank5 (duplicates removed): " + str(len(flank5_dr)) + "\n")                                                                                                                  
        
        flank5_pval = flank5_df[flank5_df["p-val"] < 0.05]
        of.write("total methylated CpGs in flank5 with P-Value < 0.05: " + str(len(flank5_pval)) + "\n")
        flank5_pval_dr = flank5_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in flank5 with P-Value < 0.05 (duplicates removed): " + str(len(flank5_pval_dr)) + "\n")

        # intergenic
        of.write("\nIntergenic:\n")
        intergenic_df = methylation_db[methylation_db["context"] == "intergenic"]
        of.write("total sequenced CpGs in intergenic: " + str(len(intergenic_df)) + "\n")
        intergenic_dr = intergenic_df.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in intergenic (duplicates removed): " + str(len(intergenic_dr)) + "\n")

        intergenic_pval = intergenic_df[intergenic_df["p-val"] < 0.05]
        of.write("total methylated CpGs in intergenic with P-Value < 0.05: " + str(len(intergenic_pval)) + "\n")
        intergenic_pval_dr = intergenic_pval.drop_duplicates(subset=['contig', 'position'], keep='last')
        of.write("total sequenced CpGs in intergenic with P-Value < 0.05 (duplicates removed): " + str(len(intergenic_pval_dr)) + "\n")

def main():
    print("load db...")
    methylation_db = pd.read_csv('methylation_db_pd.csv')
    # in jupyter notebook you can see the column names and first 10 rows of the database with
    #methylation_db.head(10)
    # or show infos about the database with 
    #methylation_db.info()


    print("summarise results...")
    summarise("summary_methylation_db.txt", methylation_db)
    sys.exit()

if __name__ == "__main__":
    main()
