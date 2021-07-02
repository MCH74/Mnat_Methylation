#! /usr/bin/python3

import re

total_c_in_CpG_count = 0
total_c_in_gene = 0
total_c_in_intergenic = 0
total_c_in_exon = 0
total_c_in_intron = 0
total_c_in_UTR = 0

total_meth_count = 0
gene_meth_count = 0
intergenic_meth_count = 0

exon_meth_count = 0
intron_meth_count = 0
UTR_meth_count = 0

ambig_meth_count = 0


with open('CpG_context_methylation_counts_ambig.txt', 'r') as inp:
    for rawline in inp:
        if rawline[0] not in ('#','\n'):
            line = rawline.strip().split('\t')
            scaff = line[0]
            position = int(line[1])
            meth_perc = float(line[2])
            reads_meth = int(line[3])
            reads_unmeth = int(line[4])
            if len(line) > 5:
                annot = line[5]
            else:
                annot = ''

            if (reads_meth + reads_unmeth) >= 10:
                total_c_in_CpG_count += 1
                if re.search('mRNA', annot):
                    total_c_in_gene += 1
                    if re.search('CDS', annot):
                        total_c_in_exon += 1
                    elif re.search('intron', annot):
                        total_c_in_intron += 1
                    elif re.search('UTR', annot):
                        total_c_in_UTR += 1
                else:
                    total_c_in_intergenic += 1


            if (reads_meth + reads_unmeth) >= 10 and meth_perc >= 10.0 and reads_meth >= 5:
                if re.search('mRNA', annot):
                    gene_meth_count += 1
                    total_meth_count += 1
                    if re.search('CDS', annot):
                        exon_meth_count += 1
                    elif re.search('intron', annot):
                        intron_meth_count += 1
                    elif re.search('UTR', annot):
                        UTR_meth_count += 1

                elif re.search('ambiguous', annot):
                        ambig_meth_count += 1
                else:
                    intergenic_meth_count += 1
                    total_meth_count += 1

overall_meth_perc = float(float(total_meth_count)/float(total_c_in_CpG_count))*100

cpg_in_gene_meth_perc = float(float(gene_meth_count)/float(total_c_in_gene))*100
cpg_in_intergene_meth_perc = float(float(intergenic_meth_count)/float(total_c_in_intergenic))*100

gene_meth_perc = float(float(gene_meth_count)/float(total_meth_count))*100
intergene_meth_perc = float(float(intergenic_meth_count)/float(total_meth_count))*100

cpg_in_exon_perc = float(float(exon_meth_count)/float(total_c_in_exon))*100
cpg_in_intron_perc = float(float(intron_meth_count)/float(total_c_in_intron))*100
cpg_in_utr_perc = float(float(UTR_meth_count)/float(total_c_in_UTR))*100

exon_meth_perc = float(float(exon_meth_count)/float(gene_meth_count))*100
intron_meth_perc = float(float(intron_meth_count)/float(gene_meth_count))*100
utr_meth_perc = float(float(UTR_meth_count)/float(gene_meth_count))*100

ambig_meth_perc = float(float(ambig_meth_count)/float(total_c_in_CpG_count))*100

print('Total CpGs covered with at least 10 reads: ', total_c_in_CpG_count)

print('Methylated CpGs (of total CpGs): ', round(overall_meth_perc, 2), '%\n')

print('Methylated CpGs in genes (of total CpGs in genes): ', round(cpg_in_gene_meth_perc, 2), '%')
print('Methylated CpGs in intergenic regions (of total CpGs in intergenic regions): ', round(cpg_in_intergene_meth_perc, 2), '%\n')

print('Methylated CpGs in genes (of total CpGs): ', round(gene_meth_perc, 2), '%')
print('Methylated CpGs in intergenic regions (of total CpGs): ', round(intergene_meth_perc, 2), '%\n')

print('Methylated CpGs in exons (of total CpGs in exons): ', round(cpg_in_exon_perc, 2), '%')
print('Methylated CpGs in introns (of total CpGs in introns): ', round(cpg_in_intron_perc, 2), '%')
print('Methylated CpGs in UTRs (of total CpGs in UTRs): ', round(cpg_in_utr_perc, 2), '%\n')

print('Methylated CpGs in exons (of total methylated CpGs in genes): ', round(exon_meth_perc, 2), '%')
print('Methylated CpGs in introns (of total methylated CpGs in genes): ', round(intron_meth_perc, 2), '%')
print('Methylated CpGs in UTRs (of total methylated CpGs in genes): ', round(utr_meth_perc, 2), '%\n')

print('Methylated CpGs with ambiguous annotation (of total CpGs): ', round(ambig_meth_perc, 2), '%')