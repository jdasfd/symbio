#!/usr/bin/env bash

usage () { echo "bash taxonomy.sh <species_name>" 1>&2; exit; }
[ $# -lt 1 ] && usage

txid=$(nwr info $1 --tsv | cut -f 1 | sed 1d)

nwr lineage $txid |
    grep -v '^no rank' |
    grep -v '^superkingdom' |
    grep -v '^species' |
    tsv-select -f 2 |
    paste -sd";"
