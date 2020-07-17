##### merge single vcf file
bcftools merge -O v -l /path/snp/script/uglist.txt | \
sed 's/\.\/\./0\/0/g' | bgzip -c > /path/snp/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz && \
tabix -p vcf /path/snp/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz

Col17_Samples=`cat /path/snp/script/hclist.txt | \
    xargs -I GVCF_FILE echo -n "-V GVCF_FILE "`

/home/wl/bin/java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
    -T GenotypeGVCFs -nt 1 -stand_call_conf 30.0 \
    -o /path/snp/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf \
    ${Col17_Samples} 2>& 1 | \
    tee /path/snp/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.log

bgzip /path/snp/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf && \
    tabix -p vcf /path/snp/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz
    
######################################################################
## Step1: Generate initial candidate targets


##
## UnifiedGenotyper
##
bgzip -dc /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.bed

bgzip -dc /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /*/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
        my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq \
    > /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.bed


##
## HaplotypeCaller
##
bgzip -dc /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.bed

bgzip -dc /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /*/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
        my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq \
    > /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.bed


##
## merge candidate target regions
##
cat /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.bed \
    /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.bed | \
    sort -k1,1 -k2,3n | bedtools merge -i - \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.bed
######################################################################




######################################################################
## Step2: Count accurate allele depths for each locus and each sample

##
## SNV sites
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
find /path/Col_mutant -name "*.bam" -print | \
    sed 's/.bam$//' | xargs -n 1 -P 12 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | cut -d"." -f1`
        
        echo "[`date`]: Start processing ${sample} ... "
        
        ## count in anomalous reads, mapping quality >= 20
        samtools mpileup -Ad100000 -q20 -f /path/TAIR10_chr_all.fasta \
            -l /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.mpileup
        
        java -jar /home/wl/Data/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.readcounts
        
        ## count in anomalous reads, mapping quality >= 0
        samtools mpileup -Ad100000 -q0 -f /path/TAIR10_chr_all.fasta \
            -l /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.mpileup
        
        java -jar /home/wl/Data/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.readcounts
        
        ## only use proper pairs, mapping quality >= 20
        samtools mpileup -d100000 -q20 -f /path/TAIR10_chr_all.fasta \
            -l /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.mpileup
        
        java -jar /home/wl/Data/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/filename.TAIR10.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.readcounts
        
        echo "[`date`]: Finished processing ${sample}"
    '


## count in anomalous reads, mapping quality >= 20
for f in `find /path/readcounts/ \
    -name "filename.TAIR10.bwa.sort.dedup.realn.snp.*.MQ20.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.AR.readcounts.list

## count in anomalous reads, mapping quality >= 0
for f in `find /path/readcounts/ \
    -name "filename.TAIR10.bwa.sort.dedup.realn.snp.*.MQ0.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ0.AR.readcounts.list

## count proper mapped reads only, mapping quality >= 20
for f in `find /path/readcounts/ \
    -name "filename.TAIR10.bwa.sort.dedup.realn.snp.*.MQ20.NAR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.NAR.readcounts.list



## UnifiedGenotyper

### with anomalous reads, mapping quality >= 20
fillVcfDepth.pl --vcf /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.vcf.gz && \
    tabix -p vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.vcf.gz

### with anomalous reads, mapping quality >= 0
fillVcfDepth.pl --vcf /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.vcf.gz && \
    tabix -p vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.vcf.gz

### without anomalous reads, mapping quality >= 20
fillVcfDepth.pl --vcf /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.vcf.gz && \
    tabix -p vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.vcf.gz


## HaplotypeCaller

### with anomalous reads, mapping quality >= 20
fillVcfDepth.pl --vcf /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.vcf.gz && \
    tabix -p vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.vcf.gz

### with anomalous reads, mapping quality >= 0
fillVcfDepth.pl --vcf /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.vcf.gz && \
    tabix -p vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.vcf.gz

### without anomalous reads, mapping quality >= 20
fillVcfDepth.pl --vcf /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.vcf.gz \
    --list /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD | bgzip -c \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.vcf.gz && \
    tabix -p vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.vcf.gz



##
## Indels
##

## UnifiedGenotyper
Samples=`find /path/Col_mutant -name "*.bam" -print | \
    xargs -I BAM_FILE echo -n "-I BAM_FILE "`
java -jar /path/GenomeAnalysisTK.jar \
    -R /path/TAIR10_chr_all.fasta \
    -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -ip 50 \
    -L /path/ug/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.bed \
    -o /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf \
    ${Samples} 2>& 1 | \
    tee /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.log

bgzip /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf && \
    tabix -p vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf.gz


## HaplotypeCaller
Samples=`find /path/Col_mutant -name "*.bam" -print | \
    xargs -I BAM_FILE echo -n "-I BAM_FILE "`
java -jar /path/GenomeAnalysisTK.jar \
    -R /path
    /TAIR10_chr_all.fasta \
    -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -ip 50 \
    -L /path/hc/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.bed \
    -o /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.vcf \
    ${Samples} 2>& 1 | \
    tee /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.log

bgzip /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.vcf && \
    tabix -p vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.vcf.gz
######################################################################






######################################################################
## Step3: Screen out candidate mutations


##
## UnifiedGenotyper
##

## substitutions

### screen according to topology (by individuals)

#### with anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.c2t3d5m5.topology7.vcf

#### with anomalous reads, mapping quality >= 0
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.c2t3d5m5.topology7.vcf

#### without anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.c2t3d5m5.topology7.vcf


### screen according to topology (by mutants)

#### with anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups4.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.c2t3d5m5.topology4.vcf

#### with anomalous reads, mapping quality >= 0
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups4.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.c2t3d5m5.topology4.vcf

#### without anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups4.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.c2t3d5m5.topology4.vcf





### screen according to frequency

#### with anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.s132c2t1d5m0.frequency.vcf

#### with anomalous reads, mapping quality >= 0
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.s132c2t1d5m0.frequency.vcf

#### without anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.s132c2t1d5m0.frequency.vcf





## indels

### screen according to topology (by individuals)
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf.gz \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology7.vcf

### screen according to topology (by mutants)
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf.gz \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 \
    -g /path/filename.groups4.txt | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology4.vcf

### screen according to frequency
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.hc_multi.vcf.gz \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 1 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 0 --max-shared-freq 132 | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.s132c2t1d5m0.frequency.vcf




##
## HaplotypeCaller
##

## substitutions

### screen according to topology (by individuals)

#### with anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.c2t3d5m5.topology7.vcf

#### with anomalous reads, mapping quality >= 0
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.c2t3d5m5.topology7.vcf

#### without anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/filename.groups7.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.c2t3d5m5.topology7.vcf

### screen according to frequency

#### with anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.s132c2t1d5m0.frequency.vcf

#### with anomalous reads, mapping quality >= 0
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.s132c2t1d5m0.frequency.vcf

#### without anomalous reads, mapping quality >= 20
detect_mutations.pl -v /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.vcf.gz \
    --max-cmp-depth 2 --max-cmp-total 1 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 132 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.s132c2t1d5m0.frequency.vcf



## indels


### screen according to topology (by mutants)
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.vcf.gz \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 \
    -g /path/filename.groups4.txt | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.mut.c2t3d5m5.topology4.vcf

### screen according to frequency
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.hc_multi.vcf.gz \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 1 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 0 --max-shared-freq 132 | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /path/TAIR10_chr_all.fasta.masked --add-pass \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.mut.s132c2t1d5m0.frequency.vcf
######################################################################




######################################################################
## Step4: Collect all cadidate mutations from different sources

##
## base substitution
##


## combine different readcount sets

### topology (by individuals)

#### UnifiedGenotyper
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.c2t3d5m5.topology7.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.c2t3d5m5.topology7.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.c2t3d5m5.topology7.vcf \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topology7.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.c2t3d5m5.topology7.vcf \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.c2t3d5m5.topology7.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.c2t3d5m5.topology7.vcf \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topology7.combined.vcf

### topology (by mutants)

#### UnifiedGenotyper
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.c2t3d5m5.topology4.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.c2t3d5m5.topology4.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.c2t3d5m5.topology4.vcf \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topology4.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.c2t3d5m5.topology4.vcf \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.c2t3d5m5.topology4.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.c2t3d5m5.topology4.vcf \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topology4.combined.vcf


### frequency

#### UnifiedGenotyper
vcf_process.pl --vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.s132c2t1d5m0.frequency.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.s132c2t1d5m0.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.s132c2t1d5m0.frequency.vcf \
    > /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.s132c2t1d5m0.frequency.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.s132c2t1d5m0.frequency.vcf \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.s132c2t1d5m0.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.s132c2t1d5m0.frequency.vcf \
    > /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.s132c2t1d5m0.frequency.combined.vcf





## combine different methods

## topology (by individuals)
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topology7.combined.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topology7.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topology7.combined.vcf

## topology (by mutants)
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topology4.combined.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topology4.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topology4.combined.vcf

## frequency
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.snp.mut.s132c2t1d5m0.frequency.combined.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.snp.mut.s132c2t1d5m0.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.s132c2t1d5m0.frequency.combined.vcf


## combine topology and frequency
vcf_process.pl --vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topology7.combined.vcf \
    --secondary-vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topology4.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag TP7 --secondary-tag TP4 --intersect-tag "TP4+TP7" | \
    vcf_process.pl --vcf - \
    --secondary-vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.s132c2t1d5m0.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.combined.vcf

perl -ne 'next if (/^\#\#/);
    if (/\#/) {
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tLength\tMethods\tMut_Type\tFrequency\tMut_Set\n"; next;
    }
    my %total = (BRCA1_T2 => 20, BRCA1_T3 => 19, MRN_T20 => 14, MutLg_T26 => 30, MutLg_T27 => 23, Polf_T23 => 17, Polf_T24 => 11);
    my @line = (split /\s+/); my $info = $line[7]; $info =~ /MA=(\w+)/; my $mut_allele = $1;
    my $fq_sum = 1; if($info =~ /Shared\=(\d+)/){$fq_sum = $1;}
    my $ffqs = (/GRPFQ\=(.*?)\;/) ? $1 : 0; my @ffqs = (split /\,/, $ffqs);
    my $ffq_sum = 0; if (@ffqs > 0) {$ffq_sum += $_ for @ffqs;} $fq_sum += $ffq_sum;
    my %counts = (); while ($line[2] =~ m/(\w+\_T\d+)/g) {$counts{$1}++;}
    my @counts = ();
    for my $indv (sort keys %counts) {
        my $tag = $counts{$indv} < $total{$indv}? "partial" : "full";
        push @counts, "$indv:$counts{$indv}/$total{$indv}($tag)";
    }
    my $counts = join ",", @counts;
    my $type = ($mut_allele eq $line[3]) ? "REF" : "ALT"; my $out_line = join "\t", @line;
    my $method = "UG_Single";
    if ($info =~ /UG_Single\+HC_GVCF/){$method = "UG_Single+HC_GVCF";}elsif($info =~ /HC_GVCF/){$method = "HC_GVCF";}
    my @filters = (); if ($info =~ /Combine=AR;/) {push @filters, "AR";} if ($info =~ /Combine=NAR;/) {push @filters, "NAR";}
    if ($info =~ /Combine=MQ20;/) {push @filters, "MQ20";} if ($info =~ /Combine=MQ0;/) {push @filters, "MQ0";}
    if ($info !~ /NMISS=0/) {push @filters, "MISSING";} if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1;/)) {push @filters, "FPD";}
    if ($method eq "UG_Single") {push @filters, "UG";} if ($line[5] < 50) {push @filters, "LowQual";}
    my $tp_set = "TP7"; if ($info =~ /Combine=TP4\+TP7/){$tp_set = "TP4+TP7";}elsif($info =~ /Combine=TP4/){$tp_set = "TP4";}
    my $set = $tp_set;
    if ($info =~ /Combine=Grouped\+NonGrouped/){$set .= "+FQ";}elsif($info =~ /Combine=NonGrouped/){$set = "FQ";}
    if (scalar @filters == 0) {$set .= "(Confidence)";} else {my $filters = join ",", @filters; $set .= "($filters)";}
    print "$out_line\t0\t$method\t$type\t$fq_sum\[$counts\]\t$set\n";' \
    /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.combined.vcf \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.snp.mut.combined.csv





##
## indels
##

## combine different methods

## topology (by individuals)
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.mut.c2t3d5m5.topology7.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology7.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.c2t3d5m5.topology7.combined.vcf

## topology (by mutants)
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.mut.c2t3d5m5.topology4.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology4.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.c2t3d5m5.topology4.combined.vcf

## frequency
vcf_process.pl --vcf /path/hc_gvcf/filename.TAIR10.bwa.sort.dedup.realn.hc.indel.mut.s132c2t1d5m0.frequency.vcf \
    --secondary-vcf /path/ug_single/filename.TAIR10.bwa.sort.dedup.realn.ug.indel.mut.s132c2t1d5m0.frequency.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.s132c2t1d5m0.frequency.combined.vcf



## combine topology and frequency
vcf_process.pl --vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.c2t3d5m5.topology7.combined.vcf \
    --secondary-vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.c2t3d5m5.topology4.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag TP7 --secondary-tag TP4 --intersect-tag "TP4+TP7" | \
    vcf_process.pl --vcf - \
    --secondary-vcf /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.s132c2t1d5m0.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.combined.vcf


perl -ne 'next if (/^\#\#/);
    if (/\#/) {
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tLength\tMethods\tMut_Type\tFrequency\tMut_Set\n"; next;
    }
    my %total = (BRCA1_T2 => 20, BRCA1_T3 => 19, MRN_T20 => 14, MutLg_T26 => 30, MutLg_T27 => 23, Polf_T23 => 17, Polf_T24 => 11);
    my @line = (split /\s+/); my $info = $line[7]; $info =~ /MA=(\w+)/; my $mut_allele = $1;
    my $fq_sum = 1; if($info =~ /Shared\=(\d+)/){$fq_sum = $1;}
    my $ffqs = (/GRPFQ\=(.*?)\;/) ? $1 : 0; my @ffqs = (split /\,/, $ffqs);
    my $ffq_sum = 0; if (@ffqs > 0) {$ffq_sum += $_ for @ffqs;} $fq_sum += $ffq_sum;
    my %counts = (); while ($line[2] =~ m/(\w+\_T\d+)/g) {$counts{$1}++;}
    my @counts = ();
    for my $indv (sort keys %counts) {
        my $tag = $counts{$indv} < $total{$indv}? "partial" : "full";
        push @counts, "$indv:$counts{$indv}/$total{$indv}($tag)";
    }
    my $counts = join ",", @counts;
    my $type = ($mut_allele eq $line[3]) ? "REF" : "ALT"; my $out_line = join "\t", @line;
    my $method = "UG_Single";
    if ($info =~ /UG_Single\+HC_GVCF/){$method = "UG_Single+HC_GVCF";}elsif($info =~ /HC_GVCF/){$method = "HC_GVCF";}
    my @filters = (); if ($info =~ /Combine=AR;/) {push @filters, "AR";} if ($info =~ /Combine=NAR;/) {push @filters, "NAR";}
    if ($info =~ /Combine=MQ20;/) {push @filters, "MQ20";} if ($info =~ /Combine=MQ0;/) {push @filters, "MQ0";}
    if ($info !~ /NMISS=0/) {push @filters, "MISSING";} if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1;/)) {push @filters, "FPD";}
    if ($method eq "UG_Single") {push @filters, "UG";} if ($line[5] < 50) {push @filters, "LowQual";}
    my $tp_set = "TP7"; if ($info =~ /Combine=TP4\+TP7/){$tp_set = "TP4+TP7";}elsif($info =~ /Combine=TP4/){$tp_set = "TP4";}
    my $set = $tp_set;
    if ($info =~ /Combine=Grouped\+NonGrouped/){$set .= "+FQ";}elsif($info =~ /Combine=NonGrouped/){$set = "FQ";}
    if (scalar @filters == 0) {$set .= "(Confidence)";} else {my $filters = join ",", @filters; $set .= "($filters)";}
    print "$out_line\tINDEL\t$method\t$type\t$fq_sum\[$counts\]\t$set\n";' \
    /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.combined.vcf \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.indel.mut.combined.csv
######################################################################




######################################################################
## Step6: generate alignments and figures for manually inspections

for record in `cat /path/combined/filename.TAIR10.bwa.sort.dedup.realn.mut.candidate.set2.vcf | \
    perl -ne 'next unless (!/\#/); my ($chrom, $pos, $sample, $info) = (split /\s+/)[0,1,2,7];
    if($sample =~ /\;/) {$sample = "Shared"} $info =~ /MA=(\w+)/; print "$sample;$chrom:$pos#$1\n";'`;
do
    spec_sample=${record/;*}
    mutation=${record/*;}
    mut_base=${mutation/*\#}
    mutation=${mutation/\#*}
    chrom=${mutation/:*}
    mut_pos=${mutation/*:}

    start_pos=`echo "${mut_pos}-200" | bc`
    end_pos=`echo "${mut_pos}+200" | bc`

    echo "${spec_sample} ${chrom} ${mut_pos} ${mut_base}"
    
    if [[ -n /path/alignments/set2 ]]; then
        mkdir -pv /path/alignments/set2
    fi
    
    echo -e "${chrom} ${start_pos} ${end_pos}\n${chrom} ${start_pos} ${mut_pos}\n${chrom} ${mut_pos} ${end_pos}\n" | \
        fasta_process.pl --rows 0 1 2 --subset 1 2 --query - \
        --fasta /path/TAIR10_chr_all.fasta.masked \
        > /path/alignments/set2/${spec_sample}.${chrom}_${mut_pos}.fa

    for bam_file in `ls /path/Col_mutant/*.TAIR10.bwa.sort.dedup.realn.bam`;
    do
        sample=`basename ${bam_file} | cut -d"." -f1`
        samtools view -X -F 3844 ${bam_file} ${chrom}:${mut_pos}-${mut_pos} | \
            awk -v name=${sample} 'BEGIN {OFS = FS = "\t"}; {print ">"name"|"$2"|"$1"\n"$10;}' \
            >> /path/alignments/set2/${spec_sample}.${chrom}_${mut_pos}.fa
    done

    reference_align.pl -i /path/alignments/set2/${spec_sample}.${chrom}_${mut_pos}.fa \
        > /path/alignments/set2/${spec_sample}.${chrom}_${mut_pos}.aln.fasta && \
        rm -v /path/alignments/set2/${spec_sample}.${chrom}_${mut_pos}.fa
done




##
## manually inspection of bam alignments
##
echo "snapshotDirectory /path/alignments/set2" \
    > /path/combined/filename.TAIR10.bwa.sort.dedup.realn.mut.run_igv.set2.txt

for bam_file in `ls /path/Col_mutant/*.TAIR10.bwa.sort.dedup.realn.bam`;
do
    echo "load ${bam_file}" \
        >> /path/combined/filename.TAIR10.bwa.sort.dedup.realn.mut.run_igv.set2.txt
done

cat /path/combined/filename.TAIR10.bwa.sort.dedup.realn.mut.candidate.set2.vcf | \
    perl -ne 'next if (/^\#/); my ($chrom, $pos, $sample, $info) = (split /\s+/)[0,1,2,7];
    if($sample =~ /\;/) {$sample = "Shared"}
    my $start=$pos-55; my $end=$pos+55;
    print "goto $chrom:$start-$end\nsnapshot $sample\.$chrom\_", "$pos.ex55", ".part1.png\n";
    $start=$pos-300; $end=$pos+300;
    print "goto $chrom:$start-$end\nsnapshot $sample\.$chrom\_", "$pos.ex300", ".part1.png\n";
    $start=$pos-3000; $end=$pos+3000;
    print "goto $chrom:$start-$end\nsnapshot $sample\.$chrom\_", "$pos.ex3000", ".part1.png\n";' \
    >> /path/combined/filename.TAIR10.bwa.sort.dedup.realn.mut.run_igv.set2.txt
######################################################################

