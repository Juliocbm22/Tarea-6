#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

echo "=========================================="
echo "4. Generar dataset"
echo "=========================================="
python3 scripts/generar_dataset.py

echo "=========================================="
echo "6. Simular reads con ART Illumina"
echo "=========================================="
art_illumina -ss HS25 -i haplotipos/hap1.fa -p -l 150 -f 10 -m 350 -s 50 -rs 101 -o fastq/hap1_
art_illumina -ss HS25 -i haplotipos/hap2.fa -p -l 150 -f 10 -m 350 -s 50 -rs 202 -o fastq/hap2_

echo "=========================================="
echo "7. Construir FASTQ de la muestra"
echo "=========================================="
cat fastq/hap1_1.fq fastq/hap2_1.fq > fastq/sample_R1.fastq
cat fastq/hap1_2.fq fastq/hap2_2.fq > fastq/sample_R2.fastq

echo "=========================================="
echo "8. Indexar y alinear con BWA"
echo "=========================================="
bwa index referencia/reference.fa
bwa mem -t 4 referencia/reference.fa fastq/sample_R1.fastq fastq/sample_R2.fastq > sam/sample.sam

echo "=========================================="
echo "9. Preparar el BAM"
echo "=========================================="
samtools view -b sam/sample.sam > bam/sample.bam
samtools sort bam/sample.bam -o bam/sample.sorted.bam
samtools index bam/sample.sorted.bam
samtools quickcheck bam/sample.sorted.bam
samtools flagstat bam/sample.sorted.bam > qc/flagstat.txt

echo "=========================================="
echo "10. faidx de la referencia"
echo "=========================================="
samtools faidx referencia/reference.fa

echo "=========================================="
echo "12 & 14. bcftools mpileup y call"
echo "=========================================="
REF="referencia/reference.fa"
BAM="bam/sample.sorted.bam"
PILEUP="resultados/sample.bcf"
RAW="resultados/raw.vcf"

bcftools mpileup -f "$REF" -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -Ou "$BAM" > "$PILEUP"
bcftools call -m -v -Ov -o "$RAW" "$PILEUP"

echo "=========================================="
echo "15. Pipe mpileup | call"
echo "=========================================="
bcftools mpileup -f "$REF" -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -Ou "$BAM" | bcftools call -m -v -Ov -o resultados/raw_pipe.vcf

echo "=========================================="
echo "20 & 21. bgzip e índice VCF"
echo "=========================================="
bgzip -f -c resultados/raw.vcf > resultados/raw.vcf.gz
bcftools index -f resultados/raw.vcf.gz

echo "=========================================="
echo "27. bcftools stats RAW"
echo "=========================================="
bcftools stats resultados/raw.vcf.gz > qc/stats_raw.txt

echo "=========================================="
echo "31 - 34. Normalizar llamadas"
echo "=========================================="
bcftools norm -m -any resultados/raw.vcf.gz -Oz -o resultados/split.vcf.gz
bcftools norm -f "$REF" -m -any resultados/raw.vcf.gz -Oz -o resultados/normalized.vcf.gz
bcftools sort resultados/normalized.vcf.gz -Oz -o resultados/normalized.sorted.vcf.gz
bcftools index -f resultados/normalized.sorted.vcf.gz

echo "=========================================="
echo "35. Normalizar TRUTH"
echo "=========================================="
bgzip -f -c truth/truth_variants.vcf > truth/truth_variants.vcf.gz
bcftools norm -f "$REF" -m -any truth/truth_variants.vcf.gz -Oz -o truth/truth.norm.vcf.gz
bcftools index -f truth/truth.norm.vcf.gz

echo "=========================================="
echo "36. bcftools isec (TRUTH vs CALLS)"
echo "=========================================="
rm -rf comparacion/isec
bcftools isec -c none truth/truth.norm.vcf.gz resultados/normalized.sorted.vcf.gz -p comparacion/isec

FN=$(grep -v '^#' comparacion/isec/0000.vcf | wc -l)
FP=$(grep -v '^#' comparacion/isec/0001.vcf | wc -l)
TP=$(grep -v '^#' comparacion/isec/0002.vcf | wc -l)

echo "TP = $TP"
echo "FP = $FP"
echo "FN = $FN"

echo "=========================================="
echo "44 - 49. BEDTools intersect"
echo "=========================================="
bedtools intersect -a resultados/normalized.sorted.vcf.gz -b bed/exons.bed -header > resultados/variants_exons.vcf
bedtools intersect -a resultados/normalized.sorted.vcf.gz -b bed/genes.bed -header > resultados/variants_genes.vcf
bedtools intersect -a resultados/normalized.sorted.vcf.gz -b bed/targets.bed -header > resultados/variants_targets.vcf
bedtools intersect -a resultados/normalized.sorted.vcf.gz -b bed/repeats.bed -header > resultados/variants_repeats.vcf
bedtools intersect -a resultados/normalized.sorted.vcf.gz -b bed/repeats.bed -v > resultados/non_repeat_variants.vcf

TOTAL=$(bcftools view -H resultados/normalized.sorted.vcf.gz | wc -l)
EXONES=$(grep -vc '^#' resultados/variants_exons.vcf || true)
GENES=$(grep -vc '^#' resultados/variants_genes.vcf || true)
TARGETS=$(grep -vc '^#' resultados/variants_targets.vcf || true)
REPEATS=$(grep -vc '^#' resultados/variants_repeats.vcf || true)
NON_REPEATS=$(grep -vc '^#' resultados/non_repeat_variants.vcf || true)

echo "Total = $TOTAL"
echo "Genes = $GENES"
echo "Exones = $EXONES"
echo "Targets = $TARGETS"
echo "Repetitivas = $REPEATS"
echo "No repetitivas = $NON_REPEATS"

echo "=========================================="
echo "50 & 51. Variantes a BED y bedtools closest"
echo "=========================================="
bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' resultados/normalized.sorted.vcf.gz > resultados/variants.bed
sort -k1,1 -k2,2n resultados/variants.bed > resultados/variants.sorted.bed
sort -k1,1 -k2,2n bed/genes.bed > bed/genes.sorted.bed
bedtools closest -a resultados/variants.sorted.bed -b bed/genes.sorted.bed -d > resultados/closest_gene.tsv

echo "=========================================="
echo "53. bedtools coverage"
echo "=========================================="
bedtools coverage -a bed/targets.bed -b bam/sample.sorted.bam > qc/target_coverage.tsv

echo "=========================================="
echo "55. bcftools concat"
echo "=========================================="
bcftools view -r chr1:1-250000 resultados/normalized.sorted.vcf.gz -Oz -o resultados/parte1.vcf.gz
bcftools view -r chr1:250001-500000 resultados/normalized.sorted.vcf.gz -Oz -o resultados/parte2.vcf.gz
bcftools index -f resultados/parte1.vcf.gz
bcftools index -f resultados/parte2.vcf.gz
bcftools concat resultados/parte1.vcf.gz resultados/parte2.vcf.gz -Oz -o resultados/reunido.vcf.gz

echo "=========================================="
echo "57. bcftools stats Normalized"
echo "=========================================="
bcftools stats resultados/normalized.sorted.vcf.gz > qc/stats_normalized.txt

echo "=========================================="
echo "PROCESO COMPLETO CON ÉXITO"
echo "=========================================="
