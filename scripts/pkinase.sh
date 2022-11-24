#!/usr/bin/env sh

usage () { echo "bash pkinase.sh <tblout.txt>" 1>&2; exit; }
[ $# -lt 1 ] && usage

TBLOUT=$1;
SPECIES=$(echo ${TBLOUT} | sed 's/\.txt$//')

cat ${TBLOUT} |
    perl -nle '
        next if /^#/;
        @hmm = split (/\s+/, $_);
        print "$hmm[2]\t$hmm[0]\t$hmm[7]";
        ' \
    > ${SPECIES}.hmm.tmp

cat ${SPECIES}.hmm.tmp |
    tsv-filter --iregex 2:PKinase --le 3:"1e-4" |
    tsv-select -f 1 |
    tsv-uniq \
    > ${SPECIES}.name.tmp

cat ${SPECIES}.hmm.tmp |
    tsv-join -k 1 -f ${SPECIES}.name.tmp \
    > ${SPECIES}.pkinase.tsv

if [ -f ${SPECIES}.pkinase.tsv ]; then
    echo "Pkinase extraction complete"
    rm ${SPECIES}.hmm.tmp ${SPECIES}.name.tmp
else
    echo "Extraction failed"
fi
