#! /usr/bin/env python

import re
import gzip
from collections import defaultdict

CpGs = defaultdict(dict)
annotation_dic = defaultdict(lambda: defaultdict(list))
final_file = []


with open('Macrotermes_natalensis.gff', 'r') as gff:
    for rawline in gff:
        if rawline[0] not in ('\n', '#'):
            line = rawline.strip().split('\t')
            scaff = line[0]
            annot = line[2]
            start = int(line[3])
            end = int(line[4])
            gid = ''
            if annot == 'mRNA':
                gid = re.search('ID=([^;]+)', line[8]).group(1)
            elif annot == 'CDS':
                gid = re.search('Parent=([^;]+)', line[8]).group(1)
            if (start, end, gid) not in annotation_dic[annot][scaff]:
                annotation_dic[annot][scaff].append((start, end, gid))


for scaf in annotation_dic['mRNA'].keys():
    for mRNAel in annotation_dic['mRNA'][scaf]:
        cds_list = []
        mRNAstart = mRNAel[0]
        mRNAend = mRNAel[1]
        for cds in annotation_dic['CDS'][scaf]:
            if cds[0] >= mRNAstart and cds[1] <= mRNAend and cds[2] == mRNAel[2]:
                cds_list.append(cds)
        cds_list.sort(key=lambda x: int(x[0]))
        if mRNAstart < cds_list[0][0]:
            if (int(mRNAstart), int(cds_list[0][0])-1) not in annotation_dic['UTR'][scaf]:
                annotation_dic['UTR'][scaf].append((int(mRNAstart), int(cds_list[0][0])-1))
        if mRNAend > cds_list[-1][1]:
            if (int(cds_list[-1][1])+1, int(mRNAend)) not in annotation_dic['UTR'][scaf]:
                annotation_dic['UTR'][scaf].append((int(cds_list[-1][1])+1, int(mRNAend)))
        for i in range(len(cds_list)-1):
            if (int(cds_list[i][1])+1, int(cds_list[i+1][0])-1) not in annotation_dic['intron'][scaf]:
                annotation_dic['intron'][scaf].append((int(cds_list[i][1])+1, int(cds_list[i+1][0])-1))


with open('CpG_OT_Old_queen_fat_body_c3_S1_L001_R1_001_trimmed_bismark_bt2.txt', 'r') as if1:
    for rawline in if1:
        if not re.match('Bismark methylation extractor', rawline):
            line = rawline.strip().split('\t')
            scaff = line[2]
            pos = line[3]
            methcall = line[4]

            try:
                CpGs[scaff][pos].append(methcall)
            except:
                CpGs[scaff][pos] = [methcall]


with open('CpG_OB_Old_queen_fat_body_c3_S1_L001_R1_001_trimmed_bismark_bt2.txt', 'r') as if2:
    for rawline in if2:
        if not re.match('Bismark methylation extractor', rawline):
            line = rawline.strip().split('\t')
            scaff = line[2]
            pos = line[3]
            methcall = line[4]

            try:
                CpGs[scaff][pos].append(methcall)
            except:
                CpGs[scaff][pos] = [methcall]


with open('Old_queen_fat_body_c3_S1_L001_R1_001_trimmed_bismark_bt2.bismark.cov', 'r') as inp:
    for rawline in inp:
        line = rawline.strip().split('\t')
        conti = line[0]
        posi = int(line[1])
        m_perc = line[3]
        m_count = line[4]
        um_count = line[5]

        if conti in CpGs and str(posi) in CpGs[conti]:
            annots = []
            for elem in annotation_dic['mRNA'][conti]:
                if posi >= elem[0] and posi <= elem[1]:
                    annots.append('mRNA')
            for elem in annotation_dic['CDS'][conti]:
                if posi >= elem[0] and posi <= elem[1]:
                    annots.append('CDS')
            for elem in annotation_dic['UTR'][conti]:
                if posi >= elem[0] and posi <= elem[1]:
                    annots.append('UTR')
            for elem in annotation_dic['intron'][conti]:
                if posi >= elem[0] and posi <= elem[1]:
                    annots.append('intron')
            if len(annots) > 2:
                annots = ['ambiguous']
            final_file.append('\t'.join((str(conti), str(posi), str(m_perc), str(m_count), str(um_count), str(','.join(annots)), '\n')))

with open('CpG_context_methylation_counts_ambig.txt', 'w') as outp:
    outp.write('# scaffold\tposition\tmethylation percentage\tcount methylated\tcount unmethylated\tannotation\n')
    for elem in final_file:
        outp.write(elem)
