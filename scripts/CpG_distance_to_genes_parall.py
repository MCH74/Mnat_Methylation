import re
import pandas as pd
import numpy as np
import math
import os
import sys
from multiprocessing import Pool

global gene_db

def apply_comp_function(x): return comp_function(x['contig'], x['position'])


def comp_function(contig, pos):
    global gene_db

    subgene_db = gene_db[gene_db['contig'] == contig]

    in_genes_db = subgene_db[(subgene_db['start'] < int(pos)) & (subgene_db['end'] > int(pos))]

    cond1 = subgene_db['start'] > int(pos)
    cond2 = subgene_db['start'] - int(pos) <= 50000
    cond3 = subgene_db['end'] < int(pos)
    cond4 = subgene_db['end'] - int(pos) >= -50000
    
    subgene1_db = subgene_db[cond1 | cond3]
    subgene_db = subgene_db[(cond1 & cond2) | (cond3 & cond4)]

    if (subgene_db.shape[0] == 0):
        if (in_genes_db.shape[0] > 0):
            if (subgene1_db.shape[0] > 0):
                return 'overlap'
            else:
                if (in_genes_db.shape[0] > 1):
                    return 'overlap'
                else:
                    return 0
        else:
            return 0

    elif (subgene_db.shape[0] > 1):
        return 'overlap'
    else:
        if (in_genes_db.shape[0] > 0):
            return 'overlap'
        else:
            if (subgene_db['start'].values[0] > int(pos)):
                res = math.ceil(abs((subgene_db['start'].values[0] - int(pos)) / 1000))
                if (subgene_db['strand'].values[0] == '+'):
                    return -res
                else:
                    return res
            else:
                res = math.ceil(abs((subgene_db['end'].values[0] - int(pos)) / 1000))
                if (subgene_db['strand'].values[0] == '+'):
                    return res
                else:
                    return -res


def multi_wrapper(t):
    sample, methylation_db, chnk = t
    print('job started on sample ' + str(sample) + ' at chunk ' + str(chnk))
    methylation_db['bin'] = methylation_db.apply(apply_comp_function, axis=1)
    return methylation_db

def main():
    global gene_db

    methylation_db = pd.read_csv('methylation_db_pd.csv', index_col=[0],
                                 dtype={'sample': int, 'caste': str, 'context': str, 'contig': str, 'position': int,
                                        'meth_percent': float, 'methylated': int, 'unmethylated': int, 'p-val': float})
    methylation_db = methylation_db.rename(columns={"p-val": "p_val"})

    gene_db = pd.read_csv('gene_db2_pd.csv', index_col=[0],
                          dtype={'contig': str, 'geneID': str, 'start': int, 'end': int, 'strand': str})

    d = {'sample': 1, 'caste': "dummy", 'context': 'dummy', 'contig': 'dummy', 'position': 1, 'meth_percent': 0.0, 'methylated': 1, 'unmethylated': 1,
             'p_val': 2.0, 'bin': "dummy"}
    methylation_bins_db = pd.DataFrame(data=d, index=[0])


    thread_arg_list = []
    grouped = methylation_db.groupby('sample')

    for sample,tdf in grouped:
        spl_lst = np.array_split(tdf, 7)
        for indx, tsdf in enumerate(spl_lst):
            chnk = indx+1
            thread_arg_list.append((sample, tsdf, chnk))

    p = Pool(70)
    mdf_results = p.map(multi_wrapper, thread_arg_list)
    p.close()
    p.join()

    results_df = pd.concat(mdf_results)
    methylation_bins_db = pd.concat([methylation_bins_db, results_df])
    methylation_bins_db = methylation_bins_db[methylation_bins_db.context != 'dummy']

    methylation_bins_db.to_csv(r'methylation_db_with_gene_dist_bins.csv')


if __name__ == "__main__":
    main()
