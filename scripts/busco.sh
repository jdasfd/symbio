#!/usr/bin/env sh

cd ~/data/symbio || exit

# GENOMES directories
if [ ! -d ./GENOMES ]; then
    echo 1>&2 "==> GENOMES unexist" 1>&2;
    exit;
fi

if [ ! -d ./BUSCO ]; then
    mkdir ~/data/symbio/BUSCO
    echo "==> mkdir BUSCO complete"
    echo
else
    echo "==> BUSCO is ready"
    echo
fi

# initialize result files
if [ ! -f ./BUSCO/busco_results.tsv ]; then
    echo -e "Species\tComplete\tSingle_copy\tMulti_copy\tMissing\tN_makers\tContigs\tN50" \
    > BUSCO/busco_results.tsv
    echo "==> initialize BUSCO results complete"
fi

# busco run for genomes
for file in $(ls GENOMES)
do
    # BUSCO for genomes
    if [ -f GENOMES/${file}/genome.fa ]; then
        # BUSCO processing
        echo "==>${file} BUSCO processing"
        busco -m genome -i GENOMES/${file}/genome.fa -o ${file} \
        --out_path ./BUSCO \
        -l viridiplantae_odb10 --quiet --cpu 10
        echo "==> BUSCO complete"
    fi

    if [ -f BUSCO/${file}/summary*.json ]; then
        cat BUSCO/${file}/summary*.json |
            jq -r -c '
            [.results.Complete,
            .results."Single copy",
            .results."Multi copy",
            .results.Missing,
            .results.n_markers,
            .results."Number of contigs",
            .results."Contigs N50"
            ] | @tsv
            ' |
            awk -v file=${file} '{print (file"\t"$0)}' \
            >> BUSCO/busco_results.tsv

            echo "==> BUSCO results summary complete"
            rm -rf BUSCO/${file}/*
    else
        echo "==> BUSCO results wrong, please check"
    fi

    rm *.log
done
