#!/usr/bin/env bash

usage () { echo "bash taxonomy.sh <species_name>" 1>&2; exit; }
[ $# -lt 1 ] && usage

# Check nwr dependencies
hash nwr 2>/dev/null || {
    echo >&2 "nwr is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install wang-q/tap/nwr"
}

# taxonomy only accept capital
species=$(echo $1 | sed 's/^\(.\)/\U\1/')

txid=$(nwr info ${species} --tsv | cut -f 1 | sed 1d)

nwr lineage $txid |
    grep -v '^no rank' |
    grep -v '^superkingdom' |
    grep -v '^species' |
    tsv-select -f 2 |
    paste -sd";"
