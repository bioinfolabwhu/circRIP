import os, sys, re
import pandas as pd
import numpy as np
import argparse
import pysam
import functools
import scipy.stats as st
from functools import reduce
from operator import itemgetter, attrgetter


class Interval(object):
    def __init__(self, x):
        self.start = x[0]
        self.end = x[1]

    def __str__(self):
        return '(%s, %s)' % (self.start, self.end)

    def intersect(self, other):
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if start < end:
            return Interval((start, end))
        return None

    def is_overlap(self, other):
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if start < end:
            return True
        return False

    def merge(self, other):
        if self.is_overlap(other):
            return Interval((min(self.start, other.start), max(self.end, other.end)))
        else:
            return None

    def to_list(self):
        return [self.start, self.end]


class GTFparser(object):
    def __init__(self, gtf):
        self.gene = {}
        self.transcript = {}
        with open(gtf) as f:
            for i in f.readlines():
                if i[0] != '#':
                    line = i.strip().split('\t')
                    if line[2] == 'gene':
                        gene_name = re.findall('gene_name "(.*?)"', i)[0]
                        gene_strand = line[6]
                        chr, start, end = line[0], int(line[3]) - 1, int(line[4])
                        self.gene.setdefault('gene_name', []).append(gene_name)
                        self.gene.setdefault('gene_strand', []).append(gene_strand)
                        self.gene.setdefault('chr', []).append(chr)
                        self.gene.setdefault('start', []).append(start)
                        self.gene.setdefault('end', []).append(end)

                    if line[2] in ['exon','UTR']:
                        exonStart, exonEnd = int(line[3]) - 1, int(line[4])
                        chr, strand = line[0], line[6]
                        geneName = re.findall('gene_name "(.*?)";', i)[0]
                        transcript_id = re.findall('transcript_id "(.*?)"; ', i)[0]
                        self.transcript.setdefault(geneName, {}).setdefault(transcript_id, []).append(
                            (exonStart, exonEnd))

        self.gene = pd.DataFrame(self.gene)

    def query(self, chr, start, end, inside):
        if inside:
            return self.gene.loc[
                ((self.gene['chr'] == chr) & (self.gene['start'] <= start) & (self.gene['end'] >= end)) ,]
        else:
            return self.gene.loc[
                ( (self.gene['chr'] == chr) & (self.gene['start'] <= start) & (
                    self.gene['end'] >= start)) | ( (self.gene['chr'] == chr) & (self.gene['start'] <= end) & (
                    self.gene['end'] >= end)),]

    def annotation(self, circchr, circstart, circend, inside):
        gquery = self.query(circchr, circstart, circend, inside)
        geneList = list(gquery['gene_name'])
        anno = []
        final_anno = []
        for g in geneList:
            for t in self.transcript[g].keys():
                exonInterval = [Interval(x) for x in self.transcript[g][t]]
                circExon = []
                n = 0
                for j in exonInterval:
                    tmp = Interval((circstart, circend)).intersect(j)
                    if tmp:
                        circExon.append(tmp.to_list())
                        n += 1
                if n > 0:
                    anno_tmp = [g, sorted(circExon, key=lambda x: x[0])]
                    anno.append(anno_tmp)
        if geneList != []:

            if anno == []:
                final_anno = [geneList[0], [[circstart, circend]], gquery.loc[gquery['gene_name']==geneList[0],'gene_strand'].iloc[0], 'intronic']
            else:
                for i in sorted(anno, key=lambda x: len(x[1]), reverse=True):
                    if i[-1][0][0] == circstart and i[-1][-1][1] == circend:
                        final_anno = i + [gquery.loc[gquery['gene_name']==i[0],'gene_strand'].iloc[0]] + ['exonic']
                        break
                if final_anno == []:
                    final_anno = sorted(anno, key=lambda x: len(x[1]), reverse=True)[0] + [gquery.loc[gquery['gene_name']==sorted(anno, key=lambda x: len(x[1]), reverse=True)[0][0],'gene_strand'].iloc[0]] + ['exonic+intron']
        else:

            final_anno = ['Unknow', [[circstart, circend]], '+', 'intergenic']
        if final_anno[1][0][0] != circstart:
            final_anno[1][0][0] = circstart
        if final_anno[1][-1][1] != circend:
            final_anno[1][-1][1] = circend
        return final_anno


class circRNA(object):
    def __init__(self, chr, start, end, redundant=5):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.redundant = int(redundant)

    def __str__(self):
        return self.chr + ':' + str(self.start) + '|' + str(self.end)

    def __eq__(self, other):
        if abs(self.start - other.start) <= self.redundant and abs(self.end - other.end) <= self.redundant \
                and self.chr == other.chr:
            return True
        return False

    def __hash__(self):
        return hash((self.chr, self.start, self.end))

    def add(self, other, attr, merge_type):
        if hasattr(self, attr) and hasattr(other, attr):
            if merge_type == 'sum':
                setattr(self, attr, getattr(self, attr) + getattr(other, attr))
            elif merge_type == 'collapse':
                setattr(self, attr, str(getattr(self, attr)) +
                        '|' + str(getattr(other, attr)))
            elif merge_type == 'distinct':
                self.add(other, attr, 'collapse')
                setattr(self, attr, '|'.join(
                    list(set(getattr(self, attr).split('|')))))

    def mergeCounts(self, other):
        if hasattr(self, 'counts') and hasattr(other, 'counts'):
            self.counts = [x + y for x, y in zip(self.counts, other.counts)]
            return self.counts
        return None

    def annotation(self, genes, inside):
        annoList = genes.annotation(self.chr, self.start, self.end, inside)
        # annoList[-1] = reduce(reduceList, annoList[-1])
        # if isinstance(annoList[-1], tuple):
        #     annoList[-1] = [x.to_list() for x in annoList[-1]]
        # else:
        #     annoList[-1] = [annoList[-1].to_list()]
        self.anno = annoList

    def getSequence(self, genome):
        sequence = ''
        dic = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N':'N'
        }
        if self.anno:
            for i in self.anno[1]:
                try:
                    sequence += genome.fetch(self.chr, i[0], i[1])
                except KeyError:
                    pass
            try:
                if self.anno[-2] == '+':
                    self.sequence = sequence
                else:
                    self.sequence = ''.join([dic[x] for x in sequence][::-1])
            except:
                print(self.chr, self.start, self.end)
                sys.exit()
            return self.sequence

    def getJunction(self, length):
        self.junction = self.sequence[-length:] + self.sequence[:length]
        return self.junction

def reduceList(x, y):
    if not isinstance(x, tuple):
        if x.is_overlap(y):
            return x.merge(y)
        else:
            return x, y
    else:
        if x[-1].is_overlap(y):
            return x[:-1] + (x[-1].merge(y),)
        else:
            return x + (y,)

def merge(lis):
    dic = {}
    dic_tmp = {}
    for circ in lis:
        dic_tmp.setdefault(str(circ), 0)
    lis = sorted(lis, key=lambda x: (x.chr, x.start))

    for i in range(len(lis)):
        if dic_tmp[str(lis[i])] == 0:
            dic_tmp[str(lis[i])] += 1
            dic.setdefault(str(lis[i]), lis[i])
            for j in range(i + 1, len(lis)):
                if lis[i] == lis[j]:
                    dic_tmp[str(lis[j])] += 1
                    dic[str(lis[i])].mergeCounts(lis[j])
                else:
                    if abs(lis[i].start - lis[j].start) > 5:
                        break
    return dic

def featureCount(gtf, ip_bam, s, p, input_bam, out, SE=False):
    if SE:
        cmd = 'featureCounts -a %s -o %s -t exon -g gene_name  -s %s --splitOnly  -T %s %s %s' \
              % (gtf, '%s.count' % out, s, p, input_bam, ip_bam)
    else:
        cmd = 'featureCounts -a %s -o %s -t exon -g gene_name  -s %s  --splitOnly -p -T %s %s %s' \
              % (gtf, '%s.count' % out, s, p, input_bam, ip_bam)
    print(cmd)
    os.system(cmd)
    return '%s.count' % out

def get_ctest_gene(ip_circ, input_circ, hostgene, host, R=1):
    gene = hostgene
    ip_input_counts = int(ip_circ) + int(input_circ)
    try:
        if gene in host.index:
            rate = R * int(host.loc[gene, 'ip']) / (int(host.loc[gene, 'input']) + R * int(host.loc[gene, 'ip']))
            ip_input_counts = int(ip_circ) + int(input_circ)
        else:
            rate = R * int(host.loc[:, 'ip'].sum()) / (
            int(host.loc[:, 'input'].sum()) + R * int(host.loc[:, 'ip'].sum()))
            ip_input_counts = int(ip_circ) + int(input_circ)
    except ZeroDivisionError:
        return 'NA'
    return ((1 - st.binom.cdf(int(ip_circ) + 1 - 1, ip_input_counts, rate)) + (
    1 - st.binom.cdf(int(ip_circ) - 1, ip_input_counts, rate))) / 2


def get_ctest_genome(ip_circ, input_circ, host, R=1):
    rate = R * int(host.loc[:, 'ip'].sum()) / (int(host.loc[:, 'input'].sum()) + R * int(host.loc[:, 'ip'].sum()))
    ip_input_counts = int(ip_circ) + int(input_circ)
    return ((1 - st.binom.cdf(int(ip_circ) + 1 - 1, ip_input_counts, rate)) + (
    1 - st.binom.cdf(int(ip_circ) - 1, ip_input_counts, rate))) / 2


def get_min(a, b):
    if a == 'NA' and b == 'NA':
        return 'NA'
    elif a != 'NA' and b == 'NA':
        return a
    elif a == 'NA' and b != 'NA':
        return b
    else:
        return min([float(a), float(b)])


def read_feature(file):
    data = pd.read_table(file, index_col=0, skiprows=1)
    data = data.iloc[:, [-1, -2]]
    data.columns = ['ip', 'input']
    return data


def getCircbed(circFile):
    data = pd.read_table(circFile, header=None)
    data.columns = ['chr', 'start', 'end',  'count']
    data.index = data['chr'] + ":" + data['start'].map(str) + "|" + data['end'].map(str)
    data = data.loc[:, ['count']]
    return data


def label(x, p=0.05, ratio=2):
    if x['pvalue'] == 'NA':
        return 'non-enriched'
    if float(x['ratio']) >= ratio and float(x['pvalue']) <= p:
        return 'enriched'
    else:
        return 'non-enriched'


def main(ip_circ, input_circ, gtf, ip_bam, input_bam, out, genomeFasta, s=2, p=30, SE=None):
    total = pd.concat([getCircbed(ip_circ), getCircbed(input_circ)], axis=1)
    total = total.fillna(0)
    genes = GTFparser(gtf)
    circRNAs = []
    for i in total.index:
        chr, start, end = re.split('[:|]', i)
        tmp_cir = circRNA(chr, start, end)
        tmp_cir.counts = list(total.loc[i, :])
        circRNAs.append(tmp_cir)
    circRNAs = merge(circRNAs)
    host = read_feature(featureCount(gtf, ip_bam, s, p, input_bam, out, SE))
    o = open(out, 'w')
    o2 = open('%s.junction.fa' % out, 'w')
    res = []
    genome = pysam.FastaFile(genomeFasta)
    for k in circRNAs.keys():
        circRNAs[k].annotation(genes,False)
        circRNAs[k].getSequence(genome)
        circRNAs[k].getJunction(50)
        p1 = get_ctest_gene(circRNAs[k].counts[0], circRNAs[k].counts[1], circRNAs[k].anno[0], host, R=1)
        p2 = get_ctest_genome(circRNAs[k].counts[0], circRNAs[k].counts[1], host, R=1)
        circRNAs[k].pvalue = get_min(p1, p2)
        if circRNAs[k].counts[1] / list(total.sum())[1] > 0:
            circRNAs[k].ratio = (circRNAs[k].counts[0] * 10 ** 6 / host['ip'].sum()) / (
            circRNAs[k].counts[1] * 10 ** 6 / host['input'].sum())
        else:
            circRNAs[k].ratio = np.inf
        circRNAs[k].label = 'non-enriched'
        if circRNAs[k].ratio >= 2 and circRNAs[k].pvalue <= 0.05:
            circRNAs[k].label = 'enriched'
        o.write('\t'.join([k, circRNAs[k].anno[0], str(int(circRNAs[k].counts[0])), \
                           str(int(circRNAs[k].counts[1])), str(circRNAs[k].counts[0] * 10 ** 6 / host['ip'].sum()),
                           str(circRNAs[k].counts[1] * 10 ** 6 / host['input'].sum()), str(circRNAs[k].ratio),
                           str(circRNAs[k].pvalue), circRNAs[k].label]) + '\n')
        res.append(
            [k, circRNAs[k].anno[0], str(int(circRNAs[k].counts[0])), \
             int(circRNAs[k].counts[1]), circRNAs[k].counts[0] * 10 ** 6 / host['ip'].sum(),
             circRNAs[k].counts[1] * 10 ** 6 / host['input'].sum(), circRNAs[k].ratio,
             circRNAs[k].pvalue, circRNAs[k].label]
        )
        o2.write('>%s\n' % k)
        o2.write('%s\n' % circRNAs[k].junction)

        with open(out + '.html', 'w') as f, open(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), 'template.html')) as tp:
            for i in tp.readlines():
                if '{{placeholder}}' in i:
                    i = 'var data = %s;\n' % str(res)
                f.write(i)

    o.close()


#main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
