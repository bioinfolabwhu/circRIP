import sys
import argparse
import subprocess


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


def mapping(index, fastq1,fastq2, out, threads, chimSegmentMin):
    cmd = 'STAR --runThreadN {threads} \
        --genomeDir {index} \
        --readFilesIn {fastq1} {fastq2} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {out} \
        --outReadsUnmapped Fastx \
        --chimSegmentMin {chimSegmentMin} \
        --chimOutType Junctions SeparateSAMold  '.format(threads=threads, index=index, fastq1=fastq1, fastq2=fastq2, out=out, chimSegmentMin=chimSegmentMin)
    run_cmd(cmd)
    return '%sChimeric.out.sam'%out, '%sUnmapped.out.mate1'%out, '%sUnmapped.out.mate2'%out


def args():
    arg = argparse.ArgumentParser()
    arg.add_argument('-fq1',help='fastq1',default=' ')
    arg.add_argument('-fq2',help='fastq2',default=' ')
    arg.add_argument('-prefix',help='out prefix')
    arg.add_argument('-chimSegmentMin',help='junction length',type=int, default=10)
    arg.add_argument('-index',default=None,type=str)
    arg.add_argument('-p', default=None, type=int)
    parser = arg.parse_args()
    return parser

#arg = args()

#mapping(arg.index, arg.fq1, arg.fq2, arg.prefix, arg.p, arg.chimSegmentMin)
