#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

REF="referencia/reference.fa"
BAM="bam/sample.sorted.bam"
PILEUP="resultados/sample.bcf"

echo "=== 12. bcftools mpileup ==="
bcftools mpileup -f "$REF" -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -Ou "$BAM" > "$PILEUP"

echo "=== 14. bcftools call ==="
RAW="resultados/raw.vcf"
bcftools call -m -v -Ov -o "$RAW" "$PILEUP"

echo "=== 15. Pipe mpileup | call ==="
bcftools mpileup -f "$REF" -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -Ou "$BAM" | bcftools call -m -v -Ov -o resultados/raw_pipe.vcf

echo "=== 20 & 21. bgzip and index ==="
bgzip -c resultados/raw.vcf > resultados/raw.vcf.gz
bcftools index -f resultados/raw.vcf.gz

echo "=== 22. Inspect view ==="
echo -n "Total raw variants: "
bcftools view -H resultados/raw.vcf.gz | wc -l
echo -n "Raw SNPs: "
bcftools view -H -v snps resultados/raw.vcf.gz | wc -l
echo -n "Raw INDELs: "
bcftools view -H -v indels resultados/raw.vcf.gz | wc -l

echo "=== 27. bcftools stats ==="
bcftools stats resultados/raw.vcf.gz > qc/stats_raw.txt
grep '^SN' qc/stats_raw.txt

echo "=== 31. Separar multialélicas ==="
bcftools norm -m -any resultados/raw.vcf.gz -Oz -o resultados/split.vcf.gz

echo "=== 32. Normalizar contra referencia ==="
bcftools norm -f referencia/reference.fa -m -any resultados/raw.vcf.gz -Oz -o resultados/normalized.vcf.gz

echo "=== 33 & 34. Ordenar e indexar normalizado ==="
bcftools sort resultados/normalized.vcf.gz -Oz -o resultados/normalized.sorted.vcf.gz
bcftools index -f resultados/normalized.sorted.vcf.gz

echo "=== 35. Normalizar TRUTH ==="
bgzip -c truth/truth_variants.vcf > truth/truth_variants.vcf.gz
bcftools norm -f referencia/reference.fa -m -any truth/truth_variants.vcf.gz -Oz -o truth/truth.norm.vcf.gz
bcftools index -f truth/truth.norm.vcf.gz

echo "=== 36. bcftools isec ==="
rm -rf comparacion/isec
bcftools isec -c none truth/truth.norm.vcf.gz resultados/normalized.sorted.vcf.gz -p comparacion/isec

FN=$(grep -v '^#' comparacion/isec/0000.vcf | wc -l)
FP=$(grep -v '^#' comparacion/isec/0001.vcf | wc -l)
TP=$(grep -v '^#' comparacion/isec/0002.vcf | wc -l)

echo "TP = $TP"
echo "FP = $FP"
echo "FN = $FN"

echo "=== 44 - 48. BEDTools intersect ==="
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

echo "=== 50 & 51. Variantes a BED y bedtools closest ==="
bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' resultados/normalized.sorted.vcf.gz > resultados/variants.bed
sort -k1,1 -k2,2n resultados/variants.bed > resultados/variants.sorted.bed
sort -k1,1 -k2,2n bed/genes.bed > bed/genes.sorted.bed
bedtools closest -a resultados/variants.sorted.bed -b bed/genes.sorted.bed -d > resultados/closest_gene.tsv

echo "=== 53. bedtools coverage ==="
bedtools coverage -a bed/targets.bed -b bam/sample.sorted.bam > qc/target_coverage.tsv

echo "=== 55. bcftools concat ==="
bcftools view -r chr1:1-250000 resultados/normalized.sorted.vcf.gz -Oz -o resultados/parte1.vcf.gz
bcftools view -r chr1:250001-500000 resultados/normalized.sorted.vcf.gz -Oz -o resultados/parte2.vcf.gz
bcftools index -f resultados/parte1.vcf.gz
bcftools index -f resultados/parte2.vcf.gz
bcftools concat resultados/parte1.vcf.gz resultados/parte2.vcf.gz -Oz -o resultados/reunido.vcf.gz

echo "=== 57. stats raw vs normalized ==="
bcftools stats resultados/raw.vcf.gz > qc/stats_raw.txt
bcftools stats resultados/normalized.sorted.vcf.gz > qc/stats_normalized.txt

echo "=== Pipeline completed successfully! ==="
