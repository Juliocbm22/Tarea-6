import random
import os

# ============================================================
# CONFIGURACIÓN Y PARÁMETROS
# ============================================================
random.seed(42)

LARGO = 500000
CHROM = "chr1"
BASES = ["A", "C", "G", "T"]

# Crear directorios necesarios
for d in ["referencia", "haplotipos", "truth", "bed", "fastq", "sam", "bam", "resultados", "qc", "scripts", "comparacion"]:
    os.makedirs(d, exist_ok=True)

# 1. GENERAR SECUENCIA DE REFERENCIA
secuencia_lista = [random.choice(BASES) for _ in range(LARGO)]

# Asegurar bases de los ejemplos del manual en posiciones 2531 y 4198 (1-based)
secuencia_lista[2530] = 'G'
secuencia_lista[4197] = 'C'

referencia = "".join(secuencia_lista)

# 2. GENERAR REGIONES ANOTADAS (BED)
genes = []
exons = []
targets = []
repeats = []

for i in range(20):
    g_start = 5000 + i * 24000
    g_end = g_start + 12000
    if g_end < LARGO:
        genes.append((CHROM, g_start, g_end, f"GEN_{i+1}"))
        for e in range(3):
            e_start = g_start + 1000 + e * 3500
            e_end = e_start + 800
            exons.append((CHROM, e_start, e_end, f"EXON_{i+1}_{e+1}"))
            t_start = max(0, e_start - 100)
            t_end = min(LARGO, e_end + 100)
            targets.append((CHROM, t_start, t_end, f"TARGET_{i+1}_{e+1}"))

for i in range(15):
    r_start = 18000 + i * 31000
    r_end = r_start + 1200
    if r_end < LARGO:
        repeats.append((CHROM, r_start, r_end, f"REPEAT_{i+1}"))

# 3. SELECCIONAR POSICIONES PARA VARIANTES (TOTAL: 290)
N_SNPS = 250
N_INS = 20
N_DEL = 20
TOTAL_VARS = N_SNPS + N_INS + N_DEL

# Fijar las dos primeras según la guía
posiciones = [2531, 4198]

# Asegurar que algunas caigan dentro de exones (~35)
exon_sample = random.sample(exons, 35)
for _, es, ee, _ in exon_sample:
    p = random.randint(es + 50, ee - 50)
    posiciones.append(p)

# Asegurar que algunas caigan dentro de repeticiones (~15)
for _, rs, re, _ in repeats:
    p = random.randint(rs + 50, re - 50)
    posiciones.append(p)

# Completar el resto de posiciones hasta 290 de forma espaciada
# Generar candidatos en todo el genoma
candidatos = list(range(6500, LARGO - 2000, 500))
random.shuffle(candidatos)

for c in candidatos:
    if len(posiciones) >= TOTAL_VARS:
        break
    # Chequear distancia mínima de 120 bp con las existentes
    if all(abs(c - p) >= 120 for p in posiciones):
        posiciones.append(c)

posiciones = sorted(posiciones[:TOTAL_VARS])

# Definir tipos de variantes
# 2 primeras son SNPs
tipos = ['SNP', 'SNP']
resto_tipos = ['SNP'] * (N_SNPS - 2) + ['INS'] * N_INS + ['DEL'] * N_DEL
random.shuffle(resto_tipos)
tipos = tipos + resto_tipos

# Construir variantes ordenadas
variantes = []
for idx, (pos, vtype) in enumerate(zip(posiciones, tipos)):
    var_id = f"VAR{idx+1}"
    if pos == 2531:
        ref_base = 'G'
        alt_base = 'A'
        gt = '0/1'
    elif pos == 4198:
        ref_base = 'C'
        alt_base = 'T'
        gt = '1/1'
    else:
        gt = '0/1' if random.random() < 0.55 else '1/1'
        if vtype == 'SNP':
            ref_base = referencia[pos - 1]
            opciones = [b for b in BASES if b != ref_base]
            alt_base = random.choice(opciones)
        elif vtype == 'INS':
            ref_base = referencia[pos - 1]
            ins_seq = "".join(random.choice(BASES) for _ in range(random.randint(1, 3)))
            alt_base = ref_base + ins_seq
        elif vtype == 'DEL':
            del_len = random.randint(1, 3)
            ref_base = referencia[pos - 1 : pos - 1 + del_len + 1]
            alt_base = ref_base[0]

    variantes.append({
        'pos': pos,
        'id': var_id,
        'ref': ref_base,
        'alt': alt_base,
        'type': vtype,
        'gt': gt
    })

# 4. CONSTRUIR HAPLOTIPOS
hap1_chars = []
hap2_chars = []

curr_pos = 1
for v in variantes:
    pos = v['pos']
    if pos > curr_pos:
        flank = referencia[curr_pos - 1 : pos - 1]
        hap1_chars.append(flank)
        hap2_chars.append(flank)
        curr_pos = pos

    gt = v['gt']
    vtype = v['type']

    if vtype == 'SNP':
        ref_b = v['ref']
        alt_b = v['alt']
        hap1_chars.append(alt_b if gt in ['0/1', '1/1'] else ref_b)
        hap2_chars.append(alt_b if gt == '1/1' else ref_b)
        curr_pos = pos + 1
    elif vtype == 'INS':
        ref_b = v['ref']
        alt_seq = v['alt']
        hap1_chars.append(alt_seq if gt in ['0/1', '1/1'] else ref_b)
        hap2_chars.append(alt_seq if gt == '1/1' else ref_b)
        curr_pos = pos + 1
    elif vtype == 'DEL':
        ref_seq = v['ref']
        alt_seq = v['alt']
        hap1_chars.append(alt_seq if gt in ['0/1', '1/1'] else ref_seq)
        hap2_chars.append(alt_seq if gt == '1/1' else ref_seq)
        curr_pos = pos + len(ref_seq)

if curr_pos <= LARGO:
    rest = referencia[curr_pos - 1 :]
    hap1_chars.append(rest)
    hap2_chars.append(rest)

hap1 = "".join(hap1_chars)
hap2 = "".join(hap2_chars)

# 5. FUNCION GUARDAR FASTA
def guardar_fasta(ruta, nombre, secuencia):
    with open(ruta, "w") as archivo:
        archivo.write(f">{nombre}\n")
        for i in range(0, len(secuencia), 80):
            archivo.write(secuencia[i:i + 80] + "\n")

# 6. GUARDAR REFERENCIA Y HAPLOTIPOS
guardar_fasta("referencia/reference.fa", CHROM, referencia)
guardar_fasta("haplotipos/hap1.fa", CHROM, hap1)
guardar_fasta("haplotipos/hap2.fa", CHROM, hap2)

# 7. GUARDAR VCF CON LA VERDAD CONOCIDA
with open("truth/truth_variants.vcf", "w") as archivo:
    archivo.write("##fileformat=VCFv4.2\n")
    archivo.write(f"##contig=<ID={CHROM},length={LARGO}>\n")
    archivo.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Synthetic variant type">\n')
    archivo.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    archivo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for v in variantes:
        archivo.write(f"{CHROM}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}\t.\tPASS\tTYPE={v['type']}\tGT\t{v['gt']}\n")

# 8. GUARDAR BEDs
def guardar_bed(ruta, intervalos):
    with open(ruta, "w") as f:
        for chrom, start, end, name in intervalos:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

guardar_bed("bed/genes.bed", genes)
guardar_bed("bed/exons.bed", exons)
guardar_bed("bed/targets.bed", targets)
guardar_bed("bed/repeats.bed", repeats)

# 9. REPORTE FINAL
n_snp = sum(1 for v in variantes if v['type'] == 'SNP')
n_ins = sum(1 for v in variantes if v['type'] == 'INS')
n_del = sum(1 for v in variantes if v['type'] == 'DEL')

print("Dataset sintético generado correctamente.")
print(f"Longitud referencia: {LARGO:,} bp")
print(f"SNPs: {n_snp}")
print(f"Inserciones: {n_ins}")
print(f"Deleciones: {n_del}")
print(f"Total variantes: {len(variantes)}")
