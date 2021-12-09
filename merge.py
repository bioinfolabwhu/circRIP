import sys
import re
import argparse
import functools
from operator import itemgetter, attrgetter

class circRNA(object):
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

    def __eq__(self, other):
        if abs(self.start-other.start) <= 5 and abs(self.end - other.end) <= 5 and self.chr == other.chr:
            return True
        return False


    def __lt__(self, other):
        if self != other:
            if self.start < other:
                return True
        return False

    def __gt__(self, other):
        if self != other:
            if self.start > other:
                return True
        return False

    def add(self, other, attr, merge_type):
        if hasattr(self, attr) and hasattr(other, attr):
            if merge_type == 'sum':
                setattr(self, attr, getattr(self, attr)+getattr(other, attr))
            elif merge_type == 'collapse':
                setattr(self, attr, str(getattr(self, attr)) +
                        '|'+str(getattr(other, attr)))
            elif merge_type == 'distinct':
                self.add(other, attr, 'collapse')
                setattr(self, attr,  '|'.join(
                    list(set(getattr(self, attr).split('|')))))

    def __str__(self):
        return self.chr + '\t' + str(self.start) + '\t' + str(self.end)


def merge(lis):
    dic = {}
    tmp = []
    dic_tmp = {}
    for circ in lis:
        dic_tmp.setdefault(str(circ), 0)
    lis = sorted(lis, key=lambda x:(x.chr, x.start ))
    for i in range(len(lis)):
        if dic_tmp[str(lis[i])] == 0:
            dic_tmp[str(lis[i])] += 1
            dic.setdefault(str(lis[i]), lis[i])
            for j in range(i+1, len(lis)):
                if lis[i] == lis[j]:
                    dic_tmp[str(lis[j])] += 1
                    dic[str(lis[i])].add(lis[j], 'count', 'sum')
                else:
                    if  abs(lis[i].start - lis[j].start) > 5:
                        break
    return dic


def bowtie2(file):
    res = []
    with open(file) as f:
        for i in f.readlines()[1:]:
            line = i.strip().split()
            circ = line[0].split('|')[0]
            circ = re.split('[:-]',  circ)
            count = int(line[1])
            circins = circRNA(circ[0].strip(), int(circ[1]), int(circ[2]))
            circins.count = count
            res.append(circins)
    res = merge(res)
    return [res[x] for x in res.keys() ]

def star(file):
    res = []
    with open(file) as f:
        for i in f.readlines():
            line = i.strip().split()
            circ = line[:3]
            count = int(line[3])
            circins = circRNA(circ[0].strip(), int(circ[1]), int(circ[2]))
            circins.count = count
            res.append(circins)
    res = merge(res)
    return [res[x] for x in res.keys() ]

def main(bowtie2_res, star_Res, out):
    lis1 = bowtie2(bowtie2_res)
    lis2 = star(star_Res)
    res = merge(lis1 + lis2)
    o = open(out, 'w')
    for i in res.keys():
        o.write('%s\t%d\n' % (res[i], res[i].count))
    o.close()
    return out

def args():
    arg = argparse.ArgumentParser()
    arg.add_argument('-remapResult', help='remap result file')
    arg.add_argument('-fristscan', help='first scan file', default=None)
    arg.add_argument('-prefix', help='out prefix')
    parser = arg.parse_args()
    return parser
