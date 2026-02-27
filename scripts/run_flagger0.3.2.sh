#!/bin/bash
#SBATCH --job-name=flagger
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=nol
#SBATCH --mem=400g 
date
## genome
matasm="gpT2T_mat_v1.0.genome.fna"
patasm="gpT2T_pat_v1.0.genome.fna"

export PATH=/Share/user/dux21/software/05.T2T_software/Winnowmap_add_I/bin/:$PATH
meryl count k=15 output merylDB <(cat "${matasm}" "${patasm}")
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -t 48 -W repetitive_k15.txt -ax map-pb -Y -L --eqx --cs -I8g <(cat "${matasm}" "${patasm}") gp_pb.fastq.gz |samtools view -hb | samtools sort -@20 > end.bam
cat "$matasm" "$patasm" > inteasm.fa
samtools faidx inteasm.fa
singularity exec `pwd`/secphase_v0.4.3.sif secphase --hifi -i end.bam -f inteasm.fa --threads 48 --outDir ./out1 --prefix gpT2T
singularity exec `pwd`/secphase_v0.4.3.sif correct_bam --primaryOnly -P ./out1/gpT2T.out.log -i end.bam -o end.corrected.bam
samtools index end.corrected.bam

## old pipeline using flagger v0.3.2
prefix="gpT2T"
assembly="inteasm.fa"
bam="end.corrected.bam"
samtools index end.bam
rustybam nucfreq -r $(sort -k2,2nr ${assembly}.fai | head -n 1 | awk 'BEGIN{FS="\t"} {printf $1 ":1-" $2}') $bam > ${prefix}.benchmark.nucfreq
expected_cov=$(awk '!/^#/{print $4+$5+$6+$7}' ${prefix}.benchmark.nucfreq | awk '{sum+=$1}END{print sum/NR}')
echo ${expected_cov}
samtools depth -aa -Q 0 $bam > ${prefix}.depth
singularity exec `pwd`/flagger_v0.3.2.sif depth2cov -d ${prefix}.depth -f $assembly.fai -o ${prefix}.cov
singularity exec `pwd`/flagger_v0.3.2.sif cov2counts -i ${prefix}.cov -o ${prefix}.counts
singularity exec `pwd`/flagger_v0.3.2.sif python3 /home/programs/src/fit_gmm.py --counts ${prefix}.counts --cov $expected_cov --output ${prefix}.table
singularity exec `pwd`/flagger_v0.3.2.sif find_blocks_from_table -c ${prefix}.cov -t ${prefix}.table -p ${prefix}
# correct false duplications
samtools depth -aa -Q 20 $bam > ${prefix}.high_mapq.depth
singularity exec `pwd`/flagger_v0.3.2.sif depth2cov -d ${prefix}.high_mapq.depth -f $assembly.fai -o ${prefix}.high_mapq.cov
cat ${prefix}.high_mapq.cov | awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= 5) {print contig"\t"$1-1"\t"$2}}' | bedtools merge -i - > ${prefix}.high_mapq.bed
bedtools subtract -a ${prefix}.duplicated.bed -b ${prefix}.high_mapq.bed > ${prefix}.dup_corrected.duplicated.bed
bedtools intersect -a ${prefix}.duplicated.bed -b ${prefix}.high_mapq.bed > ${prefix}.dup_to_hap.bed
cat ${prefix}.dup_to_hap.bed ${prefix}.haploid.bed | bedtools sort -i - | bedtools merge -i - > ${prefix}.dup_corrected.haploid.bed
date
classifications=('haploid' 'duplicated' 'collapsed' 'error')
for class in "${classifications[@]}"; do
	sort -k1,1 -k2,2n -k3,3n "${prefix}.${class}.bed" \
	| bedtools merge -i -  \
	| awk 'BEGIN{FS="\t"} {sum += $3-$2} END{if(sum) {print sum} else {print "0"}}'; done
