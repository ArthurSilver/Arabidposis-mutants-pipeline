mkdir combined
mkdir hc_gvcf
mkdir ug_single
mkdir mapped
mkdir readcounts
mkdir reads
mkdir variants
## simulate reads, reads in following situations were removed (FLAG 3844):
##    read unmapped
##    not primary alignment
##    read fails platform/vendor quality checks
##    read is PCR or optical duplicate
##    supplementary alignment
##
sim_mutation_reads.pl --fasta /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
    --depth /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/Arabi_s545.TAIR10.bwa.sort.dedup.realn.hc.sampling_depths.csv \
    --random-size 1000 --samtools "-F 3844" --exclude mitochondria chloroplast --max-shared 7 \
    --group-file /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/groups1.txt \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.dat

## extract all reads cover the simulated regions in each sample
##
awk 'BEGIN{OFS="\t"} /^\#Tag/ || $1 == "MUT" {if(/\#Tag/){$2 = "#Chrom";} print;}' \
    /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.dat | \
    cut -f 2- | body sort -k1,1 -k2,2n \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.csv

awk 'BEGIN{OFS="\t";} {if(/#Chrom/) {print "##fileformat=VCFv4.1";
        print "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Simulated read depth\">";
        print "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Variant type\">";
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION";next;}
        print $1,$2,$4,$5,$6,"999\t.\tTYPE=snp\tSDP",$7;}' \
    /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.csv \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.vcf

awk '$1 ~ /SAM:/' /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.dat | \
    cut -f 2- \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.reads.sam

cat bamlist1.txt | \
    xargs -n 1 -P 4 -I BAM_FILE sh -c '
        sample=`basename BAM_FILE | cut -d"." -f1`
        
        echo "${sample}"
        
        extract_bam_pairs.pl --samtools "-F 3844" --extend 2000 --bam BAM_FILE \
            --input /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.csv \
            --patches /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.reads.sam \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sam \
            2> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.log
        
        samtools view -H BAM_FILE | \
            cat - /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sam \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.add.sam
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/ryf/Data/tmp -jar /opt/nfs/share/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.add.sam \
            OUTPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sort.log 2>&1
        
        samtools index /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sort.bam
        
        sam2fastq.pl -i /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sam \
            -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k \
            >> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.log 2>&1 && \
            rm -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.sam \
                  /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/${sample}_ex2k.add.sam
    '

find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads -name "*.fq" | \
    xargs -n 1 -P 12 -I {} gzip -v {}

## mapping and pre-processes

find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/reads/ -name "*_1.fq.gz" | \
    sed 's/_1.fq.gz$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | sed "s/_ex2k//"`
        lane_id=`basename PREFIX | sed "s/_ex2k//"`
        
        echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
        
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"
        
        ## Align reads with BWA-MEM algorithm
        bwa mem -t 1 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
           /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta ${read1} ${read2} \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sam \
            2> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.process.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/ryf/Data/tmp -jar /opt/nfs/share/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sam \
            OUTPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.process.log 2>&1 && \
            rm -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sam
        
        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar /opt/nfs/share/biosoft/picard-tools-1.114/MarkDuplicates.jar \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.bam \
            OUTPUT=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam \
            METRICS_FILE=/mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.metrics \
            >> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.process.log 2>&1 && \
            rm -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.bam
        
        ## index bam file
        samtools index /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam
        
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T RealignerTargetCreator -nt 1 \
            -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.realn.intervals \
            -I /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam \
            >> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.process.log 2>&1
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T IndelRealigner \
            -targetIntervals /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.realn.intervals \
            -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.realn.bam \
            -I /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam \
            >> /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.process.log 2>&1 && \
            rm -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam \
                  /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/${sample}.TAIR10.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '

## UnifiedGenotyper
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/ -name "*.realn.bam" | \
    xargs -n 1 -P 4 -I BAM_FILE \
    sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo "[`date`]: Start calling ${sample_id} ... "
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T UnifiedGenotyper -glm BOTH -nt 4 -stand_call_conf 30.0 -rf MappingQuality -mmq 20 \
            -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/${sample_id}.TAIR10.bwa.sort.dedup.realn.ug.vcf \
            -I BAM_FILE \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/${sample_id}.TAIR10.bwa.sort.dedup.realn.ug.log 2>&1
        
        echo "[`date`]: Finished calling ${sample_id} ... "
    '

find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/ -name "*.vcf" -print | \
    xargs -n 1 -P 4 -I {} bgzip {}
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/ -name "*.vcf.gz" -print | \
    xargs -n 1 -P 4 -I {} tabix -p vcf {}

## merge vcf files
bcftools merge -O v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/*.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz | \
    sed 's/\.\/\./0\/0/g' \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.vcf



## HaplotypeCaller
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/ -name "*.realn.bam" -print | \
    xargs -n 1 -P 2 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo "[`date`]: Start calling ${sample_id} ... "
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
            -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/${sample_id}.TAIR10.bwa.sort.dedup.realn.hc.gvcf \
            -I BAM_FILE \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/${sample_id}.TAIR10.bwa.sort.dedup.realn.hc.log 2>&1
        
        echo "[`date`]: Finished calling ${sample_id} ... "
    '

## Joint genotyping
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/ -name "*.hc.gvcf" -print | \
    xargs -n 1 -P 4 -I {} bgzip {}
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/ -name "*.hc.gvcf.gz" -print | \
    xargs -n 1 -P 4 -I {} tabix -p vcf {}

GVCF_FILEs=`find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/variants/ -name "*.hc.gvcf.gz" -print | \
    xargs -I GVCF_FILE echo -n "-V GVCF_FILE "`
java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
    -T GenotypeGVCFs -nt 4 -stand_call_conf 30.0 \
    -o /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.vcf \
    ${GVCF_FILEs} 2>& 1 | \
    tee /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.log


## check how many sites are callable
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.vcf \
    --combine-rows 0 1 --compare-rows 2 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.combined.vcf

map_records.pl --subject /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.combined.vcf \
    --query /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.csv \
    --rows1 0 1 --rows2 0 1 | cut -f 3,9 --complement | perl -ne 's/(.*)\s+.*?Combine=(.*?)\s+/$1\t$2\t/; print;' \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.called.csv

## Detect mutations


##
## Step1: Generate initial candidate targets
##

## UnifiedGenotyper
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.vcf \
    --quality 50 --rare-only 53 \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.vcf

cat /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.vcf | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.bed


## HaplotypeCaller
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.vcf \
    --quality 50 --rare-only 53 \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.vcf

cat /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.vcf | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.bed

## merge candidate target regions
cat /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.bed \
    /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.bed | \
    sort -k1,1 -k2,3n | bedtools merge -i - \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.bed

## count reads of all alleles in each strand
##
## Base Alignment Quality (BAQ) is a new concept deployed in samtools-0.1.9+.
## It aims to provide an efficient and effective way to rule out false SNPs caused by nearby INDELs.
## The default settings for mpileup is to filter reads with bitwise flag 0X704.
## So for pileup generation the following reads will not considered at all from the bam files:
##  1) 0x0400 (aka 1024 or "d") duplicate
##  2) 0x0200 (aka 512 or "f") failed QC
##  3) 0x0100 (aka 256 or "s") non primary alignment
##  4) 0x0004 (aka 4 or "u") unmapped
## Apply -A to use anomalous read pairs in mpileup, which are not used by default (requring r874+).
find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/mapped/ -name "*.TAIR10.bwa.sort.dedup.realn.bam" -print | \
    sed 's/.bam$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | cut -d"." -f1`
        
        echo "[`date`]: Start processing ${sample} ... "
        
        ## count in anomalous reads, mapping quality >= 20
        samtools mpileup -Ad100000 -q20 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.AR.readcounts
        
        
        ## count in anomalous reads, mapping quality >= 0
        samtools mpileup -Ad100000 -q0 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ0.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
            /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ0.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ0.AR.readcounts
        
        ## only use proper pairs, mapping quality >= 20
        samtools mpileup -d100000 -q20 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.NAR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.NAR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.${sample}.MQ20.NAR.readcounts
        
        echo "[`date`]: Finished processing ${sample}"
    '


## count in anomalous reads, MQ20
for f in `find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts \
    -name "Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.*.MQ20.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.AR.readcounts.list

## count in anomalous reads, MQ0
for f in `find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/ \
    -name "Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.*.MQ0.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ0.AR.readcounts.list


## count proper mapped reads only, MQ20
for f in `find /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/readcounts/ \
    -name "Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.*.MQ20.NAR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.NAR.readcounts.list

## UnifiedGenotyper

### with anomalous reads, MQ20
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.vcf

### with anomalous reads, MQ0
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.vcf

### without anomalous reads, MQ20
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.vcf


## HaplotypeCaller

### with anomalous reads, MQ20
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.vcf

### with anomalous reads, MQ0
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.vcf

### without anomalous reads, MQ20
fillVcfDepth.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.vcf \
    --list /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.vcf

## UnifiedGenotyper
##

## substitutions

### screen according to frequency

#### with anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 1 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.mut.s41c2t3d5m5.frequency.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.mut.s41c2t3d5m5.frequency.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.mut.s41c2t3d5m5.frequency.vcf

### screen according to topology
#### with anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.vcf \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf



##
## HaplotypeCaller
##

## substitutions

### screen according to frequency

#### with anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.mut.s41c2t3d5m5.frequency.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.mut.s41c2t3d5m5.frequency.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 41 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.mut.s41c2t3d5m5.frequency.vcf

### screen according to topology
### with anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/sample-group1.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/ibm3/wl/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta.masked --add-pass \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf

## Step4: Collect all cadidate mutations from different sources

##
## base substitution
##

## combine different readcount sets

### frequency

#### UnifiedGenotyper
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.mut.s41c2t3d5m5.frequency.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.mut.s41c2t3d5m5.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.mut.s41c2t3d5m5.frequency.vcf \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.mut.s41c2t3d5m5.frequency.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.mut.s41c2t3d5m5.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.mut.s41c2t3d5m5.frequency.vcf \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf


## combine different methods

## frequency
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf

#####topology

#### UnifiedGenotyper
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.mut.c2t3d5m5.topo.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.mut.c2t3d5m5.topo.combined.vcf

## combine different methods
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/hc_gvcf/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.hc.fq53.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/ug_single/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.ug.fq53.snp.mut.c2t3d5m5.topo.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.c2t3d5m5.topo.combined.vcf

## combine topology and frequency
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.s41c2t3d5m5.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.snp.mut.combined.vcf
######################################################################

## analysis of results

## compare results
vcf_process.pl --vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.snp.mut.combined.vcf \
    --secondary-vcf /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.simulated.vars.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag FPC --secondary-tag FNC --intersect-tag TPC \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.vcf


perl -ne 'next if (/^\#\#/); if (/\#/) {
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tMethods\tMut_Type\tFrequency\tMut_Set\tCategory\n"; next;}
    my @line = (split /\s+/); my $info = $line[7]; $info =~ /MA=(\w+)/; my $mut_allele = $1;
    my $fq_sum = 1; if($info =~ /Shared\=(\d+)/){$fq_sum = $1;}
    my $type = ($mut_allele eq $line[3]) ? "REF" : "ALT"; my $out_line = join "\t", @line;
    my $method = "UG_Single"; if ($info =~ /UG_Single\+HC_GVCF/){$method = "UG_Single+HC_GVCF";}elsif($info =~ /HC_GVCF/){$method = "HC_GVCF";}
    my @filters = (); if ($info =~ /Combine=AR;/) {push @filters, "AR";} if ($info =~ /Combine=NAR;/) {push @filters, "NAR";}
    if ($info =~ /Combine=MQ20;/) {push @filters, "MQ20";} if ($info =~ /Combine=MQ0;/) {push @filters, "MQ0";}
    if ($info !~ /NMISS=0/) {push @filters, "MISSING";} if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1;/)) {push @filters, "FPD";}
    if ($method eq "UG_Single") {push @filters, "UG";} if ($line[5] < 50) {push @filters, "LowQual";}
    my $set = "TP"; if ($info =~ /Combine=Grouped\+NonGrouped/){$set .= "+FQ";}elsif($info =~ /Combine=NonGrouped/){$set = "FQ";}
    if (scalar @filters == 0) {$set .= "(Confidence)";} else {my $filters = join ",", @filters; $set .= "($filters)";}
    my $cate = "TPC"; if ($info =~ /FPC/){$cate = "FPC";}elsif ($info =~ /FNC/){$cate = "FNC";}
    print "$out_line\t$method\t$type\t$fq_sum\t$set\t$cate\n";' \
    /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.vcf \
    > /mnt/san2/usr/ryf/arabidopsis/Col_mutant/01.processed/simulation/simu1/combined/Arabi_s54.simu1.TAIR10.bwa.sort.dedup.realn.fq53.snp.mut.csv
