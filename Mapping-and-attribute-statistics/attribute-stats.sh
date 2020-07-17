## idxstats

find /path/mapping/ -name "*.bam" | sed 's/.bam$//' | xargs -n 1 -P 4 -I PREFIX \
sh -c '
echo "processing PREFIX"

samtools idxstats PREFIX.bam > PREFIX.idxtemp 

echo "finish processing PREFIX"
'

grep "" /path/mapping/*/*.idxtemp \
>> /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.idxtemp.csv

## duplications
awk 'NR==7' /path/mapping/MSH5_B13_SL5/*.metrics \
> /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.metrics.csv

for infile in `ls /path/mapping/*/*.TAIR10.bwa.sort.dedup.metrics`;
do
    echo "${infile}"
    
    awk 'NR==8' ${infile} \
        >> /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.metrics.csv
done
##
## collect metrics
##
find /path/mapping -name "*.TAIR10.bwa.sort.dedup.realn.bam" | awk '!/Parents/' | \
    sed 's/.bam$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        echo "[`date`]: Start processing PREFIX"
        
        java -jar /path/picard-tools-1.114/CollectMultipleMetrics.jar \
            REFERENCE_SEQUENCE=/opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            PROGRAM=CollectAlignmentSummaryMetrics \
            PROGRAM=QualityScoreDistribution \
            PROGRAM=CollectInsertSizeMetrics \
            PROGRAM=MeanQualityByCycle \
            INPUT=PREFIX.bam \
            OUTPUT=PREFIX.stats \
            > PREFIX.stats.log 2>&1
        
        echo "[`date`]: Finished processing PREFIX"
    '


## alignment metrics
grep "PAIR" /path/mapping/*/*.TAIR10.bwa.sort.dedup.realn.stats.alignment_summary_metrics | \
    grep -v "CATEGORY" \
    >> /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.stats.alignment_summary_metrics.csv

## insert size
for infile in `ls /path/mapping/*/*.TAIR10.bwa.sort.dedup.realn.stats.insert_size_metrics`;
do
    echo "${infile}"
    
    sample_id=`basename ${infile} | cut -d"." -f1 -`
    
    awk -v id=${sample_id} 'BEGIN{OFS="\t"} NR==8 {print id,$0;}' ${infile} \
        >> /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.stats.insert_size_metrics.csv
done


## qulaity distributions
for infile in `ls /path/mapping/*/*.TAIR10.bwa.sort.dedup.realn.stats.quality_distribution_metrics`;
do
    name=`basename ${infile}`
    sample=`echo ${name} | cut -d"." -f1`
    
    echo "${sample}"
    
    sed "s/COUNT_OF_Q/${sample}/" ${infile} | awk '!/#/ && $1' > ${infile}.csv
done


count=0
for infile in `ls /path/mapping/*/*.TAIR10.bwa.sort.dedup.realn.stats.quality_distribution_metrics.csv`;
do
    count=`expr $count + 1`
    record=`expr $count - 1`
    join -j 1 /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.stats.quality_distribution_metrics.csv.${record} $infile \
    > /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.stats.quality_distribution_metrics.csv.${count}
    rm /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.stats.quality_distribution_metrics.csv.${record}
done
##
## check GC contents
##
find /path/mapping/ -name "*.TAIR10.bwa.sort.dedup.realn.bam" | awk '!/Parents/' | \
    sed 's/.bam$//' | xargs -n 1 -P 2 -I PREFIX \
    sh -c '
        echo "[`date`]: Start processing PREFIX"
        
        java -jar /path/picard-tools-1.114/CollectGcBiasMetrics.jar \
            REFERENCE_SEQUENCE=/opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            INPUT=PREFIX.bam \
            OUTPUT=PREFIX.GcBiasMetrics.csv \
            CHART_OUTPUT=PREFIX.GcBiasMetrics.pdf \
            SUMMARY_OUTPUT=PREFIX.GcBiasMetrics.summary.csv \
            > PREFIX.GcBiasMetrics.log 2>&1
        
        echo "[`date`]: Finished processing PREFIX"
    '




##
## depth of coverage
##
find /path/mapping/ -name "*.TAIR10.bwa.sort.dedup.realn.bam" | awk '!/Parents/' | \
    xargs -n 1 -P 8 -I BAM_FILE \
    sh -c '
        out_prefix=`echo BAM_FILE | sed "s/\.bam$//"`
        
        echo BAM_FILE
        
        /home/wl/bin/java -jar /path/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T DepthOfCoverage -nt 4 -rf BadCigar \
            -omitBaseOutput -omitIntervals -omitLocusTable --nBins 99 --start 1 --stop 100 -mmq 20 \
            -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -ct 40 -ct 50 -ct 100 \
            -o ${out_prefix}.MQ20.depth -I BAM_FILE > ${out_prefix}.MQ20.depth.log 2>&1   
'

cat /path/mapping/*/*.sample_summary | grep -hv "Total" | sort -h | uniq \
    > /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.MQ20.depth.sample_summary

cat /path/mapping/*/*.sample_statistics | sort -h | uniq \
    > /path/mapping/T0_set4.TAIR10.bwa.sort.dedup.realn.MQ20.depth.sample_statistics
