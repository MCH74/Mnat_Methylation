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
from collections import defaultdict
from multiprocessing import Pool


def create_context_database(inputdir):
    """
    reads the context files (CDS, TE, intron, flank3, flank5) and stores start and end for every contig ID in a pandas database
    :param inputdir: directory with input files
    :return: context_db (pandas dataframe with genetic context information)
    """

    d = {'context': 'dummy', 'contig': 'dummy', 'start': 1, 'end': 2}
    context_db = pd.DataFrame(data=d, index=[0])
    context_lst = []
    contig_lst = []
    start_lst = []
    end_lst = []

    with open(os.path.join(inputdir, "CDS.gff")) as cds_in:
        for raw_line in cds_in:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                context_lst.append('cds')
                contig_lst.append(line[0])
                start_lst.append(int(line[1]))
                end_lst.append(int(line[2]))

    with open(os.path.join(inputdir, "TE.gff")) as intron_in:
        for raw_line in intron_in:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                context_lst.append('TE')
                contig_lst.append(line[0])
                start_lst.append(int(line[1]))
                end_lst.append(int(line[2]))

    with open(os.path.join(inputdir, "intron.gff")) as intron_in:
        for raw_line in intron_in:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                context_lst.append('intron')
                contig_lst.append(line[0])
                start_lst.append(int(line[1]))
                end_lst.append(int(line[2]))

    with open(os.path.join(inputdir, "flank3.gff")) as intron_in:
        for raw_line in intron_in:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                context_lst.append('flank3')
                contig_lst.append(line[0])
                start_lst.append(int(line[1]))
                end_lst.append(int(line[2]))

    with open(os.path.join(inputdir, "flank5.gff")) as intron_in:
        for raw_line in intron_in:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                context_lst.append('flank5')
                contig_lst.append(line[0])
                start_lst.append(int(line[1]))
                end_lst.append(int(line[2]))

    df1 = pd.DataFrame(
    {
        "context": context_lst,
        "contig": contig_lst,
        "start": start_lst,
        "end": end_lst,
    },
    )
    context_db = pd.concat([context_db, df1])

    context_db = context_db[context_db.context != 'dummy']
    return context_db

def read_methylation_file(inputdir, inp_file, context_db, caste_db, startpoint, endpoint):
    """
    reads methylation data into a pandas dataframe
    :param inputdir: directory with input file
    :param inp_file: input file
    :param context_db: pandas database with genetic context data
    :return: methylation_df: dataframe with methylation data
    """

    d = {'sample': 1, 'caste': "dummy", 'context': 'dummy', 'contig': 'dummy', 'position': 1, 'meth_percent': 0.0, 'methylated': 1, 'unmethylated': 1, 'p-val': 2}
    methylation_df = pd.DataFrame(data=d, index=[0])
    sample_lst = []
    caste_lst = []
    context_lst = []
    contig_lst = []
    pos_lst = []
    meth_perc_lst = []
    methylated_lst = []
    unmethylated_lst = []
    pval_lst = []
    sample = int(re.match("C_all_contexts_methylation_counts_(\d+)_", inp_file).group(1))
    caste = caste_db.at[sample, 'caste']

    raw_lines = []

    with open(os.path.join(inputdir, inp_file)) as meth_in:
        for i, raw_line in enumerate(meth_in):
            if i >= startpoint and i <= endpoint:
                try:
                    raw_lines.append(raw_line)
                except:
                    break
            if i > endpoint:
                break

    for raw_line in raw_lines:
        if raw_line[0] not in ('#', '\n'):
            line = raw_line.strip().split()
            contig = line[0]
            posi = int(line[1])

            condition1 = (context_db.contig == contig)

            #contig_df = context_db[context_db["contig"] == contig]
            contig_df = context_db[condition1]

            condition2 = (contig_df.start <= posi) & (contig_df.end >= posi)
            #specif_df = contig_df[(contig_df["start"] <= posi) & (contig_df["end"] >= posi)]
            specif_df = contig_df[condition2]

            if (len(specif_df.index) == 0):
                sample_lst.append(sample)
                caste_lst.append(caste)
                context_lst.append('intergenic')
                contig_lst.append(contig)
                pos_lst.append(posi)
                meth_perc_lst.append(float(line[3]))
                methylated_lst.append(int(line[4]))
                unmethylated_lst.append(int(line[5]))
                pval_lst.append(float(line[8]))
            else:
                for index, row in specif_df.iterrows():
                    sample_lst.append(sample)
                    caste_lst.append(caste)
                    context_lst.append(row["context"])
                    contig_lst.append(contig)
                    pos_lst.append(posi)
                    meth_perc_lst.append(float(line[3]))
                    methylated_lst.append(int(line[4]))
                    unmethylated_lst.append(int(line[5]))
                    pval_lst.append(float(line[8]))

    df1 = pd.DataFrame(
        {
            "sample": sample_lst,
            "caste": caste_lst,
            "context": context_lst,
            "contig": contig_lst,
            "position": pos_lst,
            "meth_percent": meth_perc_lst,
            "methylated": methylated_lst,
            "unmethylated": unmethylated_lst,
            "p-val": pval_lst,
        },
    )

    methylation_df = pd.concat([methylation_df, df1])

    methylation_df = methylation_df[methylation_df.context != 'dummy']

    return (methylation_df)

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

def multi_wrapper(t):
    inputdir, mfi, context_db, caste_db, startpoint, endpoint = t
    print("started one job on " + str(mfi))
    tdf = read_methylation_file(inputdir, mfi, context_db, caste_db, startpoint, endpoint)
    return tdf

def main():
    """
    The main function of the program.
    """

    parser = argparse.ArgumentParser(
        description='Reads methylation count files and summarises them to one file.')
    parser.add_argument("-i", "--inputdir", help="The directory with the *.txt_binom.bed_sorted sample and *.gff genetic context files (CDS, TE, intron, flank3, flank5).",
                        type=str, required=True)
    parser.add_argument("-o", "--output", help="The summary output file name (and file path).", type=str, required=True)
    args = parser.parse_args()

    caste_db = pd.DataFrame(
        {
            "sample": [1, 11, 16, 9, 14, 19, 10, 15, 5, 21, 22, 23],
            "caste": ["worker", "worker", "worker", "alate", "alate", "alate", "queen", "queen", "queen", "king", "king", "king"],
        },
    )
    caste_db = caste_db.set_index("sample")

    d = {'sample': 1, 'caste': "dummy", 'context': 'dummy', 'contig': 'dummy', 'position': 1, 'meth_percent': 0.0, 'methylated': 1, 'unmethylated': 1,
         'p-val': 2}
    methylation_db = pd.DataFrame(data=d, index=[0])

    print("start creating context database...")
    context_db = create_context_database(args.inputdir)
    print("finished creating context database...")

    methy_files = [mfi for mfi in os.listdir(args.inputdir) if mfi.endswith(".txt_binom.bed_sorted")]

    thread_arg_list = []
    for mfi in methy_files:
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 0, 500000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 500001, 1000000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 1000001, 1500000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 1500001, 2000000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 2000001, 2500000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 2500001, 3000000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 3000001, 3500000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 3500001, 4000000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 4000001, 4500000))
        thread_arg_list.append((args.inputdir, mfi, context_db, caste_db, 4500001, 5000000))


    p = Pool(80)
    mdf_results = p.map(multi_wrapper, thread_arg_list)
    #for mfi in methy_files:
    #    print("next file started: " + str(mfi))
    #    tdf = read_methylation_file(args.inputdir, mfi, context_db, caste_db)
    #    print("done reading, concatenating next...")
    #    methylation_db = pd.concat([methylation_db, tdf])
    #    print("done concatenating")
    #methylation_db = methylation_db[methylation_db.context != 'dummy']
    p.close()
    p.join()
    results_df = pd.concat(mdf_results)
    methylation_db = pd.concat([methylation_db, results_df])
    methylation_db = methylation_db[methylation_db.context != 'dummy']

    print("writing methylation db...")
    methylation_db.to_csv(r'methylation_db_pd.csv')

    print("summarise results...")
    summarise(args.output, methylation_db)
    sys.exit()

if __name__ == "__main__":
    main()
