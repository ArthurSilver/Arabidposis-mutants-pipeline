# Arabidposis-mutants-pipeline
## Mapping and attribute statistics 
### 1.fastqc.sh
Shell script of fastqc process to evaluate fastq file sequence quality <br>
### 2.mapping.sh
Shell script to align sequence file to TAIR10 reference <br>
### 3.attribute-stats.sh <br>
Shell script to evaluate bam files' attributes like depth,coverage,etc <br>
## Detect-mutation
### 1.callsnp.sh
Shell script to call mutation sites using GATK <br>
### 2.detect_mutation.sh
Shell script to screen out candidate snv and indel(insertion and deletion) sites <br>
## Simulation
### 1.simulation.sh
Shell script to simulate mutation detect pipeline to evaluate this method's FNR and callable sites ratio <br>
### 2.sim_mutation_reads.pl
Perl script to simulate reads with synthetic mutations use the method describled in Keightley et al.
