# Novel Protein Discovery Pipeline

All data are retrieved from IMG/M.

## Input files

### Proteins

Proteins are stored as a TSV file (`protein.tsv`).

| Column | Description |
|---------|-------------|
| 1 | Protein ID |
| 2 | Protein sequence |
| 3 | Start position |
| 4 | End position |

### Scaffolds

Scaffold information is stored as a TSV file (`scaffolds_info.tsv`).

| Column | Description |
|---------|-------------|
| 1 | Scaffold ID |
| 2 | Scaffold length |

---

# Filtering

To minimize the inclusion of fragmented or incomplete protein predictions while retaining potentially functional proteins, a four step filtering strategy was applied. The filtering criteria included:

1. Protein length
2. Scaffold length
3. Gene position within the scaffold
4. Presence of low complexity regions

---

## Step 1. Filter by protein length (≥35 aa)

```bash
awk -F '\t' '{print $1"\t"$2"\t"length($2)"\t"$3"\t"$4}' protein.tsv \
| awk -F '\t' '$3>=35 {print $1"\t"$2"\t"$4"\t"$5}' \
> filtered_step1_proteins.tsv
```

---

## Step 2. Filter by scaffold length (≥500 bp)

Keep scaffolds that are at least 500 bp long.

```bash
awk -F '\t' '$2>=500 {print $1}' scaffolds_info.tsv \
| sort -k1,1 \
> sorted_scaffolds_morethan500.tsv
```

Prepare the protein table for joining.

```bash
sed 's/|/\t/g' filtered_step1_proteins.tsv \
| awk -F '\t' '{print $1"|"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
| sort -k1,1 \
> sorted_filtered_step1_proteins.tsv
```

Retain only proteins located on scaffolds ≥500 bp.

```bash
join -t $'\t' -1 1 -2 1 \
sorted_scaffolds_morethan500.tsv \
sorted_filtered_step1_proteins.tsv \
> filtered_step2_proteins.tsv
```

---

## Step 3. Remove genes near scaffold edges

Genes beginning within 10 nucleotides of the scaffold start or ending within 10 nucleotides of the scaffold end are removed.

```bash
join -t $'\t' -1 1 -2 1 filtered_step2_proteins.tsv scaffolds_info.tsv \
| awk -F '\t' '$4 > 10 && $5 < $6 - 10 {print $1"\t"$2"\t"$3}' \
> filtered_step3_proteins.tsv
```

Convert the filtered proteins to FASTA.

```bash
awk '{print ">"$1"|"$2"\n"$3}' filtered_step3_proteins.tsv \
> filtered_step3_proteins.fa
```

---

## Step 4. Remove proteins with extensive low complexity regions

Mask low complexity regions using **tantan**.

```bash
tantan -p -x X filtered_step3_proteins.fa \
> masked_low_complexity_region.fa
```

Convert the masked FASTA back to TSV.

```bash
awk '/^>/ {
    if (seq) print name"\t"seq;
    name=substr($0,2);
    seq="";
    next
}
{
    seq=seq $0
}
END {
    if (seq) print name"\t"seq
}' masked_low_complexity_region.fa \
> masked_low_complexity_region.tsv
```

Proteins containing a run of ten consecutive masked residues (`XXXXXXXXXX`) are retained unchanged.

```bash
grep "XXXXXXXXXX" masked_low_complexity_region.tsv \
> keep_proteins1.tsv
```

For the remaining proteins, remove masked residues and retain only sequences that remain at least 35 amino acids long.

```bash
grep -v "XXXXXXXXXX" masked_low_complexity_region.tsv \
| sed 's/X//g' \
| awk 'length($2)>=35 {print $1"\t"$2}' \
> keep_proteins2.tsv
```

Combine the retained proteins.

```bash
cat keep_proteins1.tsv keep_proteins2.tsv \
> final_filtered_proteins.tsv
```

Generate the final FASTA file.

```bash
awk '{print ">"$1"\n"$2}' final_filtered_proteins.tsv \
> final_filtered_proteins.fa
```

---

# Searching

A protein was classified as **novel** if it showed no detectable matches to either:

* Pfam (v37)
* the IMG/M reference protein collection

---

# Search against Pfam

To facilitate analysis of the initial filtered dataset (~6 billion proteins), the FASTA file was split into approximately 6,000 smaller files. Each chunk was searched independently against the Pfam HMM database using `hmmsearch`.

```bash
hmmsearch \
    --cut_tc \
    --cpu 45 \
    --domtblout results_chunkx.domtblout \
    -o results_chunkx.hmmout \
    Pfam.hmm \
    chunkx.fa
```

After all searches completed, the `*.domtblout` files were merged into a single results table. Every protein appearing in this table was considered a Pfam hit.

## domtblout format

| Column | Description |
|---------|-------------|
| 1 | Target name (protein ID) |
| 2 | Target accession |
| 3 | Target length |
| 4 | Query name (Pfam domain) |
| 5 | Query accession |
| 6 | Query length |
| 7 | E-value |
| 8 | Full sequence bit score |
| 9 | Full sequence bias |
| 10 | Domain number |
| 11 | Total domains |
| 12 | Conditional E-value |
| 13 | Independent E-value |
| 14 | Domain bit score |
| 15 | Domain bias |
| 16 | HMM start |
| 17 | HMM end |
| 18 | Alignment start |
| 19 | Alignment end |
| 20 | Envelope start |
| 21 | Envelope end |
| 22 | Alignment accuracy |
| 23 | Target description |

## Collect Pfam hits

```bash
awk -F '\t' '{print $1}' final_hmmresults.domtblout \
| sort \
> protein_Pfam_hits.txt
```

## Remove Pfam hits

```bash
join -t $'\t' -v 1 -1 1 -2 1 \
final_filtered_proteins.tsv \
protein_Pfam_hits.txt \
| sort -k1,1 \
> pfam_novel_proteins.tsv
```

Generate FASTA.

```bash
awk '{print ">"$1"\n"$2}' pfam_novel_proteins.tsv \
> pfam_novel_proteins.fa
```

---

# Search against IMG/M reference proteins

Reference proteins are stored in four FASTA files:

* `bacteria.fa`
* `archaea.fa`
* `eukarya.fa`
* `viruses.fa`

Each database is searched independently.

## Create DIAMOND database

```bash
diamond makedb \
    --in input.fa \
    -d nr
```

## Search proteins

```bash
diamond blastp \
    -d nr.dmnd \
    -q pfam_novel_proteins.fa \
    -o matches.m8 \
    --outfmt 6 \
    qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore qlen slen
```

### Hit criteria

A protein is considered a reference match if it satisfies all of the following:

* Percent identity ≥30%
* Query coverage ≥70%
* Subject coverage ≥70%

Collect matching protein IDs.

```bash
awk '{
    if ($3>=30 &&
        (($8-$7)/$13)>=0.7 &&
        (($10-$9)/$14)>=0.7)
    print $1
}' output_table \
| sort -k1,1 \
> reference_hits_ids
```

Remove reference hits.

```bash
join -t $'\t' -v 1 -1 1 -2 1 \
pfam_novel_proteins.tsv \
reference_hits_ids \
> novel_proteins.tsv
```

Generate the final FASTA.

```bash
awk '{print ">"$1"\n"$2}' novel_proteins.tsv \
> novel_proteins.fa
```

---

# Clustering

Novel proteins are clustered using **MMseqs2 Linclust**.

Create the sequence database.

```bash
mmseqs createdb novel_proteins.fa seqDB
```

Run Linclust.

```bash
mmseqs linclust seqDB cluDB tmp \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0
```

Export cluster assignments.

```bash
mmseqs createtsv seqDB seqDB cluDB cluDB.tsv
```
