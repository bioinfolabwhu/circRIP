# circRIP
circRIP: An accurate tool for identifying circRNA-RBP interaction


Manual of circRIP v1.0.0

requisites

	pysam (0.14.1)
	pandas (1.0.3)
	featureCounts (v1.6.1)
	gtfToGenePred
	bedtools (v2.26.0)
	STAR (2.7.5a)
	bowtie2 (2.4.2)

step 1: Identify circRNA from IP/Input:

	usage:  python circRIP.py circRNAIdentify -h

	optional arguments:
    -h, --help            show this help message and exit
    -fq1 FQ1              Fastq1
    -fq2 FQ2              Fastq2
    -prefix PREFIX        Out prefix
    -chimSegmentMin CHIMSEGMENTMIN
                          Minimum length of chimeric segment length, default
                          value is 10
    -star_index STAR_INDEX
                          STAR genome index
    -bowtie2_index BOWTIE2_INDEX
                          known circRNA sequence bowtie2 index
    -p P                  Threads, default value is 1
    -gtf GTF              gtf format gene annotation file
    -genome GENOME        genome fasta file
    -SE                   Single end mode
    -junLen JUNLEN        junction length


	example:
    python circRIP.py  circRNAIdentify -fq1 test_ip_1.fq.gz -fq2 test_ip_2.fq.gz -star_index ~/starindex/ -bowtie2_index ~/bowtie2index/circBase -p 30   -gtf gencode.v32.annotation.gtf -genome GRCh38.p10.genome.fa -prefix  test_ip.circ -junLen 100

    python circRIP.py  circRNAIdentify -fq1 test_input_1.fq.gz -fq2 test_input_2.fq.gz -star_index ~/starindex/ -bowtie2_index ~/bowtie2index/circBase -p 30   -gtf gencode.v32.annotation.gtf -genome GRCh38.p10.genome.fa -prefix test_input.circ -junLen 100


  
  
step 2: Identify IP-enriched circRNA:

	usage: python circRIP.py EnrichedcircRNA -h
	
	optional arguments:
      -h, --help            show this help message and exit
      -ip_circ IP_CIRC      circRNA count in IP
      -input_circ INPUT_CIRC
                      circRNA count in Input
      -gtf GTF              gtf format gene annotation
      -ip_bam IP_BAM        BAM format files of IP
      -input_bam INPUT_BAM  BAM format files of Input
      -prefix PREFIX        output prefix
      -G G                  genome fasta
      -s S                  strand specific: 0 (unstranded), 1 (stranded) and 2
                      (reversely stranded), default value is 0
      -p P                  threads, default value is 1
      -SE                   Single end mode


	example: python circRIP.py EnrichedcircRNA -ip_circ test_ip.circ -input_circ test_input.circ -gtf gencode.v32.annotation.gtf -ip_bam test_ip.circAligned.out.bam -input_bam test_input.circAligned.out.bam -prefix test.final.out -G GRCh38.p10.genome.fa
