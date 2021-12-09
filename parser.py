import os
import re
import sys
from collections import Counter

import pysam
import argparse
from more_itertools import pairwise


class GenePred(object):
    def __init__(self, file):
        self.file = file
        self.genes = {}
        self.gene_pos = []
        with open(self.file) as f:
            for i in f.readlines():
                line = i.strip().split()
                gene_name, exon_start, exon_end = line[11], line[8].strip(
                    ',').split(','), line[9].strip(',').split(',')
                chrom, strand, transcript = line[1], line[2], line[0]
                type = 'mRNA' if line[5] != line[6] else 'lncRNA'
                exons = [[int(x), int(y)]
                         for x, y in zip(exon_start, exon_end)]
                intron = [[x[1], y[0]]
                          for x, y in pairwise(exons)] if len(exons) > 1 else []
                self.genes.setdefault(gene_name, []).append(
                    [transcript, chrom, strand, type, exon_start, exon_end, exons, intron])
        for k in self.genes.keys():
            start = min([int(x[4][0]) for x in self.genes[k]])
            ends = max([int(x[5][-1]) for x in self.genes[k]])
            chrom = [x[1] for x in self.genes[k]]
            if len(set(chrom)) > 1:
                continue
            self.gene_pos.append([chrom[0], start, ends, k])


def parser(sam, SE=None):
    sam = pysam.AlignmentFile(sam)
    dic = {}
    res = []
    dic_segment = {}
    for i in sam.fetch():
        dic.setdefault(i.qname, {}).setdefault(
            'chrom', []).append(i.reference_name)
        if i.is_read1:
            dic.setdefault(i.qname, {}).setdefault('read1', []).append(i)
        elif i.is_read2:
            dic.setdefault(i.qname, {}).setdefault('read2', []).append(i)
        else:
            dic.setdefault(i.qname, {}).setdefault('read', []).append(i)
        dic_segment.setdefault(i.qname, []).append(i.query_alignment_length)
    for k in dic.keys():
        pos = []
        if len(list(set(dic[k]['chrom']))) == 1:
            lis = []
            if not SE:
                if len(dic[k]['read1']) == 2 and len(dic[k]['read2']) == 1:
                    lis = dic[k]['read1']
                    n = 'read2'
                elif len(dic[k]['read2']) == 2 and len(dic[k]['read1']) == 1:
                    lis = dic[k]['read2']
                    n = 'read1'
                else:
                    continue
            else:
                if len(dic[k]) == 2:
                    lis = dic[k]['read']
                    n = 'read'
            if (lis[0].qstart < lis[1].qstart and lis[0].reference_start > lis[1].reference_start):
                pos = [lis[1].reference_start, lis[0].reference_end]
            elif (lis[0].qstart > lis[1].qstart and lis[0].reference_start < lis[1].reference_start):
                pos = [lis[0].reference_start, lis[1].reference_end]
            if pos:
                if dic[k][n][0].reference_start >= pos[0] and dic[k][n][0].reference_end <= pos[1]:
                    res.append([dic[k]['chrom'][0], pos[0], pos[1], k])
    return res, dic_segment


def getInterect(res, gene_pos, out):
    o = open('%s_tmp.bed' % (out), 'w')
    for i in res:
        line = '\t'.join([str(x) for x in i])
        o.write('%s\n' % line)
    o.close()
    o = open('%s_gene_tmp.bed' % (out), 'w')
    for i in gene_pos:
        line = '\t'.join([str(x) for x in i])
        o.write('%s\n' % line)
    o.close()
    os.system('bedtools intersect -a %s -b %s -wao > %s_tmp.bed.anno' %
              ('%s_tmp.bed' % (out), '%s_gene_tmp.bed' % (out), out))
    return '%s_tmp.bed.anno' % out


def get_annotation(file, dic_segment, genes):
    total_res = []
    reads = {}
    with open(file) as f:
        for i in f.readlines():
            line = i.strip().split('\t')
            reads.setdefault(line[3], []).append(line)
    for r in reads.keys():
        n = 0
        for line in reads[r]:
            if n == 0:
                circ = line[:3]
                if line[5] == '-1':
                    if int(circ[2]) - int(circ[1]) <= 10000 and min(dic_segment[r]) >= 10:
                        total_res.append(
                            ['intergenic', '+', 'intergenic', circ[0], circ[1], circ[2], 0])
                        n += 1
                if int(line[1]) >= int(line[5]) and int(line[2]) <= int(line[6]):
                    gene = line[-2]
                    gene_exon_starts = []
                    gene_exon_ends = []
                    gene_introns = []
                    for t in genes[gene]:
                        gene_exon_starts += t[4]
                        gene_exon_ends += t[5]
                        gene_introns.append(t[7])
                    gene_exon_starts = list(set(gene_exon_starts))
                    gene_exon_starts.sort()
                    gene_exon_ends = list(set(gene_exon_ends))
                    gene_exon_ends.sort()
                    tag = [0, 0]
                    pos = [0, 0]
                    for es in gene_exon_starts:
                        if int(circ[1]) in range(int(es) - 5, int(es) + 5):
                            tag[0] = 1
                            pos[0] = int(es)
                    for ee in gene_exon_ends:
                        if int(circ[2]) in range(int(ee) - 5, int(ee) + 5):
                            tag[1] = 1
                            pos[1] = int(ee)
                    if tag == [1, 1]:
                        total_res.append(
                            [gene, genes[gene][0][2], 'exon', circ[0], circ[1], circ[2], 1])
                        n += 1
                    else:
                        if len(gene_introns) > 0:
                            for trans_intron in gene_introns:
                                if n == 0:
                                    for intron in trans_intron:
                                        if int(circ[1]) in range(intron[0] - 5, intron[0] + 5) and int(circ[2]) in range(intron[1] - 5, intron[1] + 5):
                                            total_res.append(
                                                [gene, genes[gene][0][2], 'intron', circ[0], circ[1], circ[2], 1])
                                            n += 1
                                        elif int(circ[1]) in range(intron[0] - 5, intron[0] + 5) and int(circ[2]) <=  intron[1] +5:
                                            total_res.append(
                                                [gene, genes[gene][0][2], 'intron', circ[0], circ[1], circ[2], 1])
                                            n += 1
                                        elif int(circ[1]) >=  intron[0] -5  and int(circ[2]) in range(intron[1] - 5, intron[1] + 5):
                                            total_res.append(
                                                [gene, genes[gene][0][2], 'intron', circ[0], circ[1], circ[2], 1])
                                            n += 1
    return total_res


def get_rc(s):
    dic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    a1 = ''
    for i in s:
        a1 += dic[i]
    return a1


def is_splice(x):
    start = re.search('AG', x[0])
    end = re.search('GT', x[1])
    if start and end:
        return True, start.span()[1], end.span()[0], '+'
    start = re.search('AC', x[0])
    end = re.search('CT', x[1])
    if start and end:
        return True, start.span()[1], end.span()[0], '+'
    start = re.search('TG', get_rc(x[0]))
    end = re.search('GA', get_rc(x[1]))
    if start and end:
        return True, start.span()[1], end.span()[0], '-'
    start = re.search('TC', get_rc(x[0]))
    end = re.search('CA', get_rc(x[1]))
    if start and end:
        return True, start.span()[1], end.span()[0], '-'
    return False, False, False, False


def check_splice_site(fasta, dic):
    fa = pysam.FastaFile(fasta)
    lis = []
    for i in dic:
        if i[-1] != 1:
            x = []
            chrom, start, end = i[3], i[4], i[5]
            if int(end) - int(start) <= 10000:
                x.append(fa.fetch(chrom, int(start) - 2, int(start) + 3))
                x.append(fa.fetch(chrom, int(end) - 2, int(end) + 3))
                state, start_l, end_l, strand = is_splice(x)
                if state:
                    if i[2] == 'intergenic':
                        i[1] = strand
                    lis.append(i[:-1])
        else:
            lis.append(i[:-1])
    return lis


def get_counts(lis, out):
    lis_counter = Counter(lis)
    o = open(out, 'w')
    for i in lis_counter.keys():
        line = i.strip().split()
        circ = '%s\t%s\t%s' % (line[3], line[4], line[5])
        gene = line[0]
        strand = line[1]
        circType = line[2]
        count = lis_counter[i]
        if circType == 'intergenic':
            if count >= 2:
                o.write('%s\t%s\n' %
                        (circ,str(count)))
                continue
            else:
                continue
        o.write('%s\t%s\n' %
                (circ, str(count)))
    o.close()


def main(chimericSam, gtf, out, genome, SE):

    test, dic_segment = parser(chimericSam, SE)
    os.system('gtfToGenePred -genePredExt -geneNameAsName2 %s %s.GenePred' % (gtf, out))
    a = GenePred('%s.GenePred'%out)
    anno = getInterect(test, a.gene_pos, out)
    b = get_annotation(anno, dic_segment, a.genes)
    c = check_splice_site(genome, b)
    d = ['\t'.join([str(y) for y in x]) for x in c]
    get_counts(d, out)
    os.system('rm %s %s %s %s' % ('%s_tmp.bed' %
                               (out), '%s_gene_tmp.bed' % (out), anno,'%s.GenePred'%out ))

    return out

def args():
    arg = argparse.ArgumentParser()
    arg.add_argument('-chimericSam',help='fastq1')
    arg.add_argument('-gtf',help='fastq2',default=None)
    arg.add_argument('-prefix',help='out prefix')
    arg.add_argument('-genome',help='genome fasta',type=str)
    arg.add_argument('-SE',action="store_true",help='is single end data')
    parser = arg.parse_args()
    return parser



