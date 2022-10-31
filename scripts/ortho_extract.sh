#!/usr/bin/env sh

usage () { echo "bash ortho_extract.sh <species_accession>" 1>&2; exit; }
[ $# -lt 1 ] && usage

S_ID=$1;

if ! grep -q -i $S_ID $(pwd)/Orthogroups.tsv; then
    echo 1>&2 "[${S_ID}] is not identified" 1>&2; exit;
else
    cat $(pwd)/Orthogroups.tsv | head -n 1 > $(pwd)/${S_ID}.tsv
    cat $(pwd)/Orthogroups.tsv | grep "${S_ID}" >> $(pwd)/${S_ID}.tsv
fi

echo "==> Transposing"
cat $(pwd)/${S_ID}.tsv |
    datamash transpose |
    mlr --itsv --ocsv cat \
    > $(pwd)/${S_ID}.csv

if [[ ! -s $(pwd)/${S_ID}.csv ]]; then
    echo 1>&2 "Orthogroups.csv extraction failed" 1>&2;
    rm $(pwd)/${S_ID}.tsv && exit
else
    echo "Transposing succeed"
    rm $(pwd)/${S_ID}.tsv
fi

echo "==> Processing END"
