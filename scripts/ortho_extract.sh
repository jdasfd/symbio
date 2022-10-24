#!/usr/bin/env sh

usage () { echo "bash ortho_extract.sh <species_accession>" 1>&2; exit; }
[ $# -lt 1 ] && usage

S_ID=$1;

if ! grep -q -i $S_ID ~/data/symbio/Orthogroups/Orthogroups.tsv; then
    echo 1>&2 "[${S_ID}] is not identified" 1>&2; exit;
else
    cat ~/data/symbio/Orthogroups/Orthogroups.tsv | head -n 1 > ~/data/symbio/Orthogroups/${S_ID}.tsv
    cat ~/data/symbio/Orthogroups/Orthogroups.tsv | grep "${S_ID}" >> ~/data/symbio/Orthogroups/${S_ID}.tsv
fi

echo "==> Transposing"
cat ~/data/symbio/Orthogroups/${S_ID}.tsv |
    datamash transpose |
    mlr --itsv --ocsv cat \
    > ~/data/symbio/Orthogroups/${S_ID}.csv

if [[ ! -s ~/data/symbio/Orthogroups/${S_ID}.csv ]]; then
    echo 1>&2 "Orthogroups.csv extraction failed" 1>&2;
    rm ~/data/symbio/Orthogroups/${S_ID}.tsv && exit
else
    echo "Transposing succeed"
    rm ~/data/symbio/Orthogroups/${S_ID}.tsv
fi

echo "==> Processing END"
