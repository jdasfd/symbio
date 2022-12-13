# Chlorophyta RLK identification

## Preparation

Preparation process were recorded in the [RLK_identification.md](RLK_identification.md). All codes below were just the repeat among more algae data.

```bash
mkdir -p ~/data/chlorophyta
cd ~/data/chlorophyta
mkdir info PROTEINS
```

## Select Chlorophyta proteins

All chlorophyta proteins were selected for covering the largest range of green algae. Some streptophyta proteins also included.
All proteins were saved into `~/data/chlorophyta/PROTEINS`.

```bash
cd ~/data/chlorophyta
ls PROTEINS |
    sed 's/\.pep$//' \
    > info/algae.lst


cat algae.lst | wc -l
#31
```

## Identify RLK

### Scan domains via `hmmscan`

```bash
cd ~/data/chlorophyta
mkdir -p ~/data/chlorophyta/DOMAIN/pfam

cat algae.lst |
    parallel -j 8 --keep-order '
        echo {}
        hmmscan --cpu 2 -E 0.1 --domE 0.1 -o DOMAIN/pfam/{}.txt \
            ../symbio/HMM/PFAM/Pfam-A.hmm PROTEINS/{}.pep
    '

# all query sequences
ls DOMAIN/pfam/*.txt |
    parallel -j 16 --keep-order '
        cat {} |
            grep -v '^#' |
            grep '^Query' |
            grep '\]$' |
            wc -l
        ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#457986

ls DOMAIN/pfam/*.txt |
    parallel -j 16 --keep-order '
        echo "==> {/.}"
        perl ../symbio/scripts/hmm_results.pl -i {} \
        > DOMAIN/pfam/{/.}.pfam.tsv
    '

# check whether query seqs contained domains
ls DOMAIN/pfam/*.pfam.tsv |
    parallel -j 12 --keep-order '
        cat {} |
            tsv-select -H -f QUERY |
            tsv-uniq |
            wc -l
    ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#343391
```
