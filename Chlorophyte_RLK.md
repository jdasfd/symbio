# Chlorophyta RLK identification

## Preparation

Preparation process were recorded in the [RLK_identification.md](RLK_identification.md). All codes below were just the repeat among more algae data.

```
# script need python2
sudo apt install python2.7
```

```bash
mkdir -p ~/data/chlorophyta
cd ~/data/chlorophyta
mkdir info PROTEINS
```

## Select Chlorophyta proteins

### Genomes

All chlorophyta proteins were selected for covering the largest range of green algae. Some streptophyta proteins also included.

All proteins were saved into `~/data/chlorophyta/GENOMES`.

```bash
mkdir -p ~/data/chlorophyta/info
cd ~/data/chlorophyta

ls GENOMES > info/algae.lst


cat info/algae.lst | wc -l
#31
```

### Extract longest transcripts

```bash
cd ~/data/chlorophyta
mkdir PROTEINS

cat info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        python2.7 scripts/getLongestProteinFromGFF.py \
            GENOMES/{}/genome.pep GENOMES/{}/genome.gff \
            PROTEINS/{}.pep
    '
```

### Rename all protein files

```bash
cd ~/data/chlorophyta

```

## Identify RLK

### Count basic info of proteins

```bash
cd ~/data/chlorophyta
mkdir -p ~/data/chlorophyta/info/PRO_len

cat info/algae.lst |
    parallel -j 16 -k '
        faops size PROTEINS/{}.pep \
        > info/PRO_len/{}.length.tsv
    '

cat info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat info/PRO_len/{}.length.tsv |
            cut -f 1 |
            awk -v spe={} '\''{print ($0"\t"spe)}'\'' \
        >> info/pro_spe.tsv
    '
#457986
```

### Scan domains via `hmmscan`

```bash
cd ~/data/chlorophyta
mkdir -p ~/data/chlorophyta/DOMAIN/pfam

cat info/algae.lst |
    parallel -j 8 --keep-order '
        echo {}
        hmmscan --cpu 2 -E 0.1 --domE 0.1 -o DOMAIN/pfam/{}.txt \
            ../symbio/HMM/PFAM/Pfam-A.hmm PROTEINS/{}.pep
    '

# all query sequences
cat info/algae.lst |
    parallel -j 16 --keep-order '
        cat DOMAIN/pfam/{}.txt |
            grep -v "^#" |
            grep '^Query' |
            grep '\]$' |
            wc -l
        ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#457986
# all proteins were treated as query

cat info/algae.lst |
    parallel -j 16 --keep-order '
        echo "==> {}"
        perl ../symbio/scripts/hmm_results.pl -i DOMAIN/pfam/{}.txt \
        > DOMAIN/pfam/{}.tsv
    '

# check whether query seqs contained domains
cat info/algae.lst |
    parallel -j 12 --keep-order '
        cat DOMAIN/pfam/{}.tsv |
            tsv-select -H -f QUERY |
            tsv-uniq |
            wc -l
    ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#343391
#114595 proteins were failed in extracting domains
```

### Extract kinase domain

The reason of choosing kinase domains were recorded in the [RLK_identification.md](RLK_identification.md).

All 3 domains are:

- Pkinase
- PK_Tyr_Ser-Thr
- Pkinase_fungal

```bash
cd ~/data/chlorophyta/DOMAIN
mkdir -p KD

# all pkinase protein names
cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat pfam/{}.tsv |
            tsv-filter -H --or \
                --str-eq Domain:Pkinase \
                --str-eq Domain:Pkinase_fungal \
                --str-eq Domain:PK_Tyr_Ser-Thr |
            tsv-filter -H --le E_value:"1e-4" |
            tsv-select -H -f QUERY |
            sed 1d |
            tsv-uniq \
        > KD/{}.KD.lst
    '

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        if [ -f KD/{}.E_sort.KD.tsv ]; then
            rm KD/{}.E_sort.KD.tsv
        fi
        for acc in $(cat KD/{}.KD.lst)
        do
            cat pfam/{}.tsv |
                tsv-filter -H --str-eq QUERY:${acc} |
                sed 1d |
                tsv-select -f 1,2,3,4,5 |
                sort -gk 3,3 \
            >> KD/{}.E_sort.KD.tsv
        done
    '
```

### Filter all repeat domains

All domains were filtered according to the following discipline:

If the range of a domain overlapped with another one, the domain with less E-value will be kept.

```bash
cd ~/data/chlorophyta/DOMAIN
mkdir UNIQ

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat KD/{}.E_sort.KD.tsv |
            perl ../../symbio/scripts/domain_uniq.pl |
            tsv-select -f 1,4,5,2 \
        > UNIQ/{}.uniq.tsv
    '

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat UNIQ/{}.uniq.tsv |
            tsv-filter --or \
                --str-eq 4:Pkinase \
                --str-eq 4:Pkinase_fungal \
                --str-eq 4:PK_Tyr_Ser-Thr |
            tsv-select -f 1 |
            tsv-uniq \
        > UNIQ/{}.uniq.KD.lst
        cat UNIQ/{}.uniq.tsv |
            tsv-join -f UNIQ/{}.uniq.KD.lst -k 1 \
        > UNIQ/{}.uniq.KD.tsv
    '
```

### Identify TMD via `TMHMM2`

- Extract sequences containing kinase domain

```bash
cd ~/data/chlorophyta/DOMAIN
mkdir KD_SEQ

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        faops some ../PROTEINS/{}.pep \
            KD/{}.KD.lst KD_SEQ/{}.KD.fa
        NUM_KD=$(cat KD/{}.KD.lst | wc -l)
        NUM_FA=$(faops size KD_SEQ/{}.KD.fa | wc -l)
        if [[ "$NUM_KD" == "$NUM_FA" ]]; then
            echo "Extraction complete"
        else
            echo "Something wrong"
            rm KD_SEQ/{}.KD.fa
        fi
    '
```

- Upload file

Upload each PRO.fa to [TMHMM2](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) for prediction. Results could be recorded into the following files.

```bash
cd ~/data/chlorophyta/DOMAIN
mkdir TMD

cat ../info/algae.lst |
    parallel -j 1 -k '
        touch TMD/{}.tsv
    '

dos2unix TMD/*

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        perl ../../symbio/scripts/tmhmm_result.pl -i TMD/{}.tsv \
        > TMD/{}.TMD.tsv
    '

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat TMD/{}.TMD.tsv |
            cut -f 1 |
            tsv-uniq \
        > TMD/{}.TMD.lst
    '
```

### Combine results for RLK identification

RLK: N-ECD-TMD-KD-C

```bash
cd ~/data/chlorophyta/DOMAIN
mkdir COMBINE

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        for acc in $(cat TMD/{}.TMD.lst)
        do
            cat TMD/{}.TMD.tsv UNIQ/{}.uniq.KD.tsv |
                tsv-filter --str-eq 1:${acc} |
                sort -nk 2,2 \
            > COMBINE/{}.sort.KD_TMD.tsv
        done
    '

cat ../info/algae.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat COMBINE/{}.sort.KD_TMD.tsv |
            cut -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={} '\''{print spe"\t"$0}'\'' \
        > ../info/spe_KD.tsv
    '
```

```bash
sudo mount -t drvfs F: /mnt/f
ln -s /mnt/f/data/ ~/data
```
