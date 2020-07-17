for fastq in `find /path/Cleandata/ -name "*_1.fq.gz"`;
do
    echo ${fastq}
    gzip -dc ${fastq} | head -n 200001 | tail -n 1 | perl -ne 's/\@//; s/\:\d+\:\d+\s+\d\:N\:0/ /; print;'
done | paste - - > /path/T1_set6.lane_ids.csv

fastqc --noextract --format fastq --threads 6 find /path/Cleandata/*.fq.gz

## get total sequenced reads and GC contents
for qc in `find /path/Cleandata/ -name "*_fastqc.zip"`;
do
sample=`dirname ${qc}`
sample=`basename ${sample}`
echo -n "${sample}:"
unzip -p ${qc} */fastqc_data.txt | awk '/Filename/ || /Total Sequence/ || /^%GC/' | \
paste - - - | sed 's/\s/\t/g' | cut -f 2,5,7
done | sort -h | paste - - >/path/T1_set6.fastqc.csv

## get failed item
for qc in `find /path/Cleandata/ -name "*_fastqc.zip"`;
do
sample=`dirname ${qc}`
sample=`basename ${sample}`
echo "${sample}:"
unzip -p ${qc} */summary.txt | awk '/FAIL/' | sed 's/\s/\t/'
done >/path/T1_set6.fastqc.summary.csv
