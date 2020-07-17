#!/bin/sh


cat /path/mapping/filelist.txt | sed 's/_1.fq.gz$//' | \
    xargs -n 1 -P 6 -I PREFIX \
    sh -c '
        lane_id=`basename PREFIX`
        sample=`dirname PREFIX`
        sample=`basename ${sample}`
        mkdir /path/mapping/${sample}
        
        echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
        
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"
        
        ## Align reads with BWA-MEM algorithm
        bwa mem -t 4 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
            /path/TAIR10_chr_all.fasta ${read1} ${read2} \
            > /path/mapping/${sample}/${sample}.TAIR10.bwa.sam \
            2> /path/mapping/${sample}/${sample}.TAIR10.bwa.log
        
        ## sort bam file
        java -jar -Djava.io.tmpdir=/home/temp /path/picard-tools-1.114/SortSam.jar \
            INPUT=/path/mapping/${sample}/${sample}.TAIR10.bwa.sam \
            OUTPUT=/path/mapping/${sample}/${sample}.TAIR10.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> /path/mapping/${sample}/${sample}.TAIR10.bwa.log 2>&1 && \
            rm -v /path/mapping/${sample}/${sample}.TAIR10.bwa.sam
        
        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar -Djava.io.tmpdir=/home/ryf/Data/temp /path/picard-tools-1.114/MarkDuplicates.jar \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=/path/mapping/${sample}/${sample}.TAIR10.bwa.sort.bam \
            OUTPUT=/path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam \
            METRICS_FILE=/path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.metrics \
            >> /path/mapping/${sample}/${sample}.TAIR10.bwa.log 2>&1 && \
            rm -v /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.bam
        
        ## index bam file
        samtools index /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /path/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /path/TAIR10_chr_all.fasta \
            -T RealignerTargetCreator -nt 4 \
            -o /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.realn.intervals \
            -I /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam \
            >> /path/mapping/${sample}/${sample}.TAIR10.bwa.log 2>&1
        
        java -jar /path/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /path/TAIR10_chr_all.fasta \
            -T IndelRealigner \
            -targetIntervals /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.realn.intervals \
            -o /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.realn.bam \
            -I /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam \
            >> /path/mapping/${sample}/${sample}.TAIR10.bwa.log 2>&1 && \
            rm -v /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam \
                  /path/mapping/${sample}/${sample}.TAIR10.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '
