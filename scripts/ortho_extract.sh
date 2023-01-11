#!/usr/bin/env sh

usage () { echo "bash ortho_extract.sh <species_accession | ortho_num> <out_dir>" 1>&2; exit; }
[ $# -lt 1 ] && usage

S_ID=$1;
OUT_PATH=$2;

# path parameter exists or not
if [ ! -d ${OUT_PATH} ]; then
    echo "dir not found." 1>&2 && exit;
fi

# Orthogroups.tsv whether exists
if [ ! -f $(pwd)/Orthogroups.tsv ]; then
    echo "Orthogroups.tsv not found." 1>&2 && exit;
fi

if ! grep -q -i $S_ID $(pwd)/Orthogroups.tsv; then
    echo 1>&2 "[${S_ID}] is not identified" 1>&2; exit;
else
    cat $(pwd)/Orthogroups.tsv | head -n 1 > ${OUT_PATH}/${S_ID}.tsv
    cat $(pwd)/Orthogroups.tsv | grep "${S_ID}" >> ${OUT_PATH}/${S_ID}.tsv
fi

echo "==> Transposing"
cat ${OUT_PATH}/${S_ID}.tsv |
    datamash transpose |
    mlr --itsv --ocsv cat \
    > ${OUT_PATH}/${S_ID}.csv

if [[ ! -s ${OUT_PATH}/${S_ID}.csv ]]; then
    echo 1>&2 "Orthogroups.csv extraction failed for ${S_ID}." 1>&2;
    rm ${OUT_PATH}/${S_ID}.tsv && exit
else
    echo "Transposing succeed"
    rm ${OUT_PATH}/${S_ID}.tsv
fi

echo "==> Orthogroups extraction complete"
