import argparse,os
from starmapping import mapping
import parser
import remap
import ctest,merge

def circRNAIdentify(arg):
    chimeric, unmap1, unmqp2 = mapping(arg.star_index, arg.fq1, arg.fq2, arg.prefix, arg.p, arg.chimSegmentMin)

    parser.main(chimeric, arg.gtf, '%s_SplitAlignment.txt'%arg.prefix,arg.genome, arg.SE)

    remap.JunctionMapping(unmap1,unmqp2, arg.bowtie2_index, arg.p, \
                          '%s_PseudoReference.txt'%arg.prefix, arg.chimSegmentMin, arg.junLen,arg.SE)

    merge.main('%s_PseudoReference.txt'%arg.prefix,'%s_SplitAlignment.txt'%arg.prefix, arg.prefix)
    os.system('rm %s %s'%( '%s_PseudoReference.txt'%arg.prefix,'%s_SplitAlignment.txt'%arg.prefix ))

def EnrichedcircRNA(arg):
    ctest.main(arg.ip_circ, arg.input_circ, arg.gtf, arg.ip_bam, \
               arg.input_bam, arg.prefix,arg.G, arg.s, arg.p, arg.SE, arg.L, arg.FC, arg.pvalue)

def args():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    # circRNAIdentify
    circRNAIdentifyParm = sub.add_parser('circRNAIdentify', help='circRNA identification')
    circRNAIdentifyParm.add_argument('-fq1', help='Fastq1', default=' ')
    circRNAIdentifyParm.add_argument('-fq2', help='Fastq2', default=' ')
    circRNAIdentifyParm.add_argument('-prefix', help='Out prefix')
    circRNAIdentifyParm.add_argument('-chimSegmentMin', help='Minimum length of chimeric segment length, default value is 10', type=int, default=10)
    circRNAIdentifyParm.add_argument('-star_index',help='STAR genome index', default=None, type=str)
    circRNAIdentifyParm.add_argument('-bowtie2_index',help='known circRNA sequence bowtie2 index ', default=None, type=str)
    circRNAIdentifyParm.add_argument('-p',help='Threads, default value is 1',action="store", type=int,default=1)
    circRNAIdentifyParm.add_argument('-gtf', help='gtf format gene annotation file', default=None)
    circRNAIdentifyParm.add_argument('-genome', help='genome fasta file', type=str)
    circRNAIdentifyParm.add_argument('-SE', action="store_true", help='Single end mode')
    circRNAIdentifyParm.add_argument('-junLen', help='junction length', type=int)
    circRNAIdentifyParm.set_defaults(handle=circRNAIdentify)
    # EnrichedcircRNA
    EnrichedcircRNAParm = sub.add_parser('EnrichedcircRNA', help='Find IP-enriched circRNA')
    EnrichedcircRNAParm.add_argument('-ip_circ',help='circRNA count in IP')
    EnrichedcircRNAParm.add_argument('-input_circ',help='circRNA count in Input')
    EnrichedcircRNAParm.add_argument('-gtf',help='gtf format gene annotation')
    EnrichedcircRNAParm.add_argument('-ip_bam',help='BAM files of IP',type=str)
    EnrichedcircRNAParm.add_argument('-input_bam',type=str,help='BAM files of Input')
    EnrichedcircRNAParm.add_argument('-prefix',help='output prefix',type=str)
    EnrichedcircRNAParm.add_argument('-G',help='genome fasta')
    EnrichedcircRNAParm.add_argument('-s',help='strand specific: 0 (unstranded), 1 (stranded) and 2 (reversely stranded), default value is 0',default=0)
    EnrichedcircRNAParm.add_argument('-p',help='threads, default value is 1',action="store", type=int,default=1)
    EnrichedcircRNAParm.add_argument('-pvalue', help='pvalue threshold', action="store", type=float, default=0.05)
    EnrichedcircRNAParm.add_argument('-FC', help='IP/Input foldchange threshold', action="store", type=float, default=2)
    EnrichedcircRNAParm.add_argument('-L', help='lambda threshold', action="store", type=float, default=1)
    EnrichedcircRNAParm.add_argument('-SE',action="store_true",required=False,help='Single end mode',default=False)
    EnrichedcircRNAParm.set_defaults(handle=EnrichedcircRNA)

    return parser.parse_args()

arg = args()
arg.handle(arg)

