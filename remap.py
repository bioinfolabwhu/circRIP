
import sys,os,re,subprocess,pysam
import pandas as pd
import argparse

def run_cmd(x):
    sys.stderr.write('running %s\n' % x)
    out = subprocess.run(
        x, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out.returncode != 0:
        sys.stderr.write('cmd running error\n')
        sys.stderr.write('args: %s\n' % out.args)
        sys.stderr.write('log: %s\n' % out.stdout)
        sys.stderr.write('log: %s\n' % out.stderr)
    else:
        sys.stderr.write('done\n')


def buildIndex(file, juncLen, prefix):
    o = open('%s.junction'%prefix, 'w')
    with open(file) as f:
        for i in f.readlines():
            if '>' == i[0]:
                o.write('%s\n'%i.strip())
            else:
                seq = i.strip()
                l = seq[:int(juncLen)]
                r = seq[-int(juncLen):]
                o.write('%s\n'%(r+l))
    o.close()
    run_cmd('bowtie2-build %s %s'%('%s.junction'%prefix,  prefix))
    return  prefix


def mapping(p, index,prefix, fq1,fq2, SE):
    if not SE:
        cmd = 'bowtie2  -p {p}  \
        -x  {index} \
        -1 {fq1} \
        -2 {fq2} \
        --end-to-end --no-unal -S {prefix}.sam  2>{prefix}.bowtie2.log'.format(p=p, index=index, fq1=fq1,fq2=fq2, prefix=prefix)
        run_cmd(cmd)
    else:
        cmd = 'bowtie2 -p {p} \
                -x  {index} \
                -U {fq1} \
                --end-to-end --no-unal -S {prefix}.sam  2>{prefix}.bowtie2.log'.format(p=p, index=index, fq1=fq1, prefix=prefix)
        run_cmd(cmd)
    return '%s.sam'%(prefix)


def get_df(sam, l, mymin):
    file =  pysam.AlignmentFile(sam)
    dic = {}
    for i in file.fetch():
        if not i.is_unmapped and i.get_tag('NM') <= 2:
            if i.reference_start <= l - mymin and i.reference_end >= l + mymin:
                if i.reference_name in dic.keys():
                    dic[i.reference_name] += 1
                else:
                    dic[i.reference_name] = 1
    data = pd.DataFrame([dic])
    return  data.T

def main(arg):
    if arg.buildIndex:
        index = buildIndex(arg.circFasta, arg.junLen,  arg.prefix)
    else:
        sam = mapping(arg.p, arg.index, arg.prefix, arg.fq1, arg.fq2)
        df = get_df(sam, int(arg.junLen), int(arg.minSpan))
        df.to_csv(arg.prefix, sep='\t')

def BuildIndex(circFasta, junLen, prefix):
    buildIndex(circFasta, junLen, prefix)

def JunctionMapping(fq1,fq2, index, p, prefix , minSpan, junLen,SE):
    sam = mapping(p, index, prefix, fq1, fq2,SE)
    df = get_df(sam, int(junLen), int(minSpan))
    df.to_csv(prefix, sep='\t')


def args():
    arg = argparse.ArgumentParser()
    arg.add_argument('-fq1',help='fastq1')
    arg.add_argument('-fq2',help='fastq2',default=None)
    arg.add_argument('-prefix',help='out prefix')
    arg.add_argument('-junLen',help='junction length',type=int)
    arg.add_argument('-minSpan',type=int,default=5)
    arg.add_argument('-circFasta',help='circRNA fasta',default=None)
    arg.add_argument('-index',default=None,type=str)
    arg.add_argument('-p', default=None, type=int)
    arg.add_argument('-buildIndex',action="store_true",help='build index mode')
    parser = arg.parse_args()
    return parser

