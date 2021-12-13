# circRIP

circRIP: An accurate tool for identifying circRNA-RBP interaction

## 1. requisites

    pysam (0.14.1)
    pandas (1.0.3)
    featureCounts (v1.6.1)
    gtfToGenePred
    bedtools (v2.26.0)
    STAR (2.7.5a)
    bowtie2 (2.4.2)

## 2. Identify circRNA from IP/Input
###     2.1 Identify circRNA using built-in methon in circRIP:

        circRIP identify circRNA by parsing STAR mapping results and re-mapping the unmapped reads to known circRNA back-splice junctions(BSJs)
    
        requisited input files:
	
            test_ip_1.fq.gz,test_ip_2.fq.gz,test_input_1.fq.gz,test_input_2.fq.gz  fastq files of IP/Input
            ~/starindex/    the STAR index
            ~/bowtie2index/circBase the bowtie2 index of circRNA junctions 
            gencode.v32.annotation.gtf  the GTF format gene annotation
            GRCh38.p10.genome.fa    the genome fasta file
            
        command lines:
	
		# IP
		python circRIP.py  circRNAIdentify -fq1 test_ip_1.fq.gz -fq2 test_ip_2.fq.gz \
		                   -star_index ~/starindex/ -bowtie2_index ~/bowtie2index/circBase \
		                   -p 30 -gtf gencode.v32.annotation.gtf -genome GRCh38.p10.genome.fa -prefix test_ip.circ -junLen 100
		# Input
		python circRIP.py  circRNAIdentify -fq1 test_input_1.fq.gz -fq2 test_input_2.fq.gz \
		                    -star_index ~/starindex/ -bowtie2_index ~/bowtie2index/circBase \
		                    -p 30 -gtf gencode.v32.annotation.gtf -genome GRCh38.p10.genome.fa -prefix test_input.circ -junLen 100


        
    results:
    
        the column of test_ip/input.circ is:
	
            1, Chromosome
            2, circRNA start coordinates
            3, circRNA end coordinates
            4, circRNA BSJs reads counts

###     2.2 Identify circRNA using other softwares:

     Users are free to use the circRNA recognition software, as long as the following conditions are met:
     
        1, transform the results to four columns format as described above.
        2, bam files of linear reads(to detecte forward-splice junctions, splice-aware aligners such as STAR should be used to perform reads mapping)
            
        
## 3. Identify IP-enriched circRNA:

circRIP identify IP-enriched circRNA by c-test:
    
    requisited input files: 
    
        test_ip.circ/test_input.circ    circRNAs in four columns format
        test_ip.circAligned.out.bam/test_input.circAligned.out.bam    bam files of linear reads
        gencode.v32.annotation.gtf    the GTF format gene annotation
        GRCh38.p10.genome.fa    the genome fasta file
        
    command lines:
       

	python circRIP.py EnrichedcircRNA -ip_circ test_ip.circ -input_circ test_input.circ -gtf gencode.v32.annotation.gtf \
	                  -ip_bam test_ip.circAligned.out.bam -input_bam test_input.circAligned.out.bam -prefix test.final.out \
			  -G GRCh38.p10.genome.fa


    results:
        circRIP generates two output files test.final.out and test.final.out.html, test.final.out is the txt format results and test.final.out.html is a web report page
        the column of test.final.out is:
            1, circRNA ID(Chromosome:start|end)
            2, circRNA host gene symbol
            3, circRNA BSJs in IP
            4, circRNA BSJs in Input
            5, circRNA CPM in IP
            6, circRNA CPM in Input
            7, IP/Input ratio
            8, P value
            9, circRNA status (enriched or non-enriched)
        
        test.final.out.html:
![image](https://github.com/bioinfolabwhu/imges/blob/main/demo.jpg)
            
