#####HC
cat bamlist.txt | sed 's/.bam$//' | xargs -n 1 -P 6 -I PREFIX \
sh -c '
sample=`basename PREFIX | cut -d"." -f1`

echo "start processing $sample"

mkdir /path/variants/HC/${sample}
       
java -jar /GenomeAnalysisTK.jar \
-R TAIR10_chr_all.fasta \
-T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
-o /path/variants/HC/$sample/${sample}.TAIR10.bwa.sort.dedup.realn.hc.gvcf \
-I PREFIX.bam \
> /path/variants/HC/$sample/${sample}.TAIR10.bwa.sort.dedup.realn.hc.log 2>&1

echo ">> Finished processing ${sample_id}"
'

ls /path/variants/HC/*/*.TAIR10.bwa.sort.dedup.realn.hc.gvcf | \
    xargs -n 1 -P 6 -I {} bgzip {}
ls /path/variants/HC/*/*.TAIR10.bwa.sort.dedup.realn.hc.gvcf.gz | \
    xargs -n 1 -P 6 -I {} tabix -p vcf {}

#####UG
cat bamlist.txt | sed 's/.bam$//' | xargs -n 1 -P 3 -I PREFIX \
sh -c '
sample=`basename PREFIX | cut -d"." -f1`

mkdir /path/variants/UG/${sample}

echo "start processing $sample"

java -jar /GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
-T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -nt 1 -rf MappingQuality -mmq 20 \
-o /path/variants/UG/$sample/${sample}.TAIR10.bwa.sort.dedup.realn.ug.vcf \
-I PREFIX.bam \
> /path/variants/UG/$sample/${sample}.TAIR10.bwa.sort.dedup.realn.ug.log 2>&1
        
echo ">> Finished processing ${sample_id}"
'

ls /path/variants/UG/*/*.TAIR10.bwa.sort.dedup.realn.ug.vcf | \
    xargs -n 1 -P 6 -I {} bgzip {}
ls /path/variants/UG/*/*.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz | \
    xargs -n 1 -P 6 -I {} tabix -p vcf {}
