# RLK identification

## Preparation

Visit [Anaconda](https://www.anaconda.com/) to download install script.

- BUSCO

Visit [BUSCO userguide](https://busco.ezlab.org/busco_userguide.html#conda-package) for database.

```bash
# install anaconda
bash Anaconda3-2022.05-Linux-x86_64.sh
# enter to allow and complete installation
# start linux again, directly activating base

# avoid directly activating base when start linux
conda config --set auto_activate_base false
# conda install busco
conda create -n genome -c conda-forge -c bioconda busco=5.4.2
# activate env
conda activate genome
```

- OrthoFinder

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install orthofinder
conda install orthofinder
# primary_transcript.py to get the proteins from the primary transcripts
```

- HMMER

```bash
mkdir -p ~/data/symbio/HMM/PFAM
cd ~/data/symbio/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

hmmpress Pfam-A.hmm
```

- Other software

```bash
brew install dos2unix
```

## CDS or protein selection

### Genomes selection and CDS renaming

There are something to declare: 

- Multiple transcripts were not excluded here, which means different types of the species transcripts were not specified. Maybe the primary transcripts would be adopted afterwards.
- All genomes were put into directories named by `species`, then moved into directory `GENOMES`. Every directory basically had four files including genome, gff, cds and pep.
- Download directly from URL (no records in database)
  - fagopyrum_esculentum
  - allium_sativum
  - elaeis_guineensis
  - passiflora_edulis
  - daemonorops_jenkinsiana
  - sphagnum_magellanicum
  - adiantum_capillus-veneris (no longest gff)
  - ceratopteris_richardii, marchantia_paleacea
  - carica_papaya
  - taxus_chinensis
  - sphagnum_fallax
- Discard
  - punica_granatum (gff has problem)
  - picea_abies (ftp cannot be accessed)
  - phoenix_dactylifera

```bash
mkdir ~/data/symbio/info
cd ~/data/symbio/info

# manually pick some genomes into taxo_genom.xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f taxo_genom.xlsx |
    tsv-select -d, -f 1 |
    sed 1d |
    perl -nl -e '$_ = lc $_;print $_' \
    > genome.lst

wc -l genome.lst
#148 genome.lst

# process genomes
cd ~/data/symbio/GENOMES

for dir in $(ls)
do
    echo "==>{dir}"
    faops size ${dir}/genome.cds |
        cut -f 1 |
        awk -v spe=${dir} '{print ($0"\t"spe_$0)}' \
        > tmp
    faops replace ${dir}/genome.cds tmp ../CDS/${dir}.cds.fa
    rm tmp
done

# there would have a question: faops size will automatically stop at the first space
# so some cds files would be named incorrectly
# check whether changed successfully
cd ~/data/symbio

cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo "==>{}"
        faops size CDS/{}.cds.fa | cut -f 1 | head -n 1
    '
# manually check those cds done wrong

# manually change those species with problems
# following is an example for dealing with the `.||.||.` type cds
cd ~/data/symbio/GENOMES/<species>

# <species> means you should replace to the specific species
cat genome.cds |
    perl -nle '
        print uc $_ if /^[AGCTN]/i;
        print "><species>_$1" if /^>.+\|(.+?)\|\|-*\d\|\|CDS/;
    ' \
    > ../../CDS/<species>.cds.fa 

# finally check the renamed cds numbers
cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo "==>{}"
        CDS=$(cat CDS/{}.cds.fa | faops size stdin | wc -l)
        GEN=$(cat GENOMES/{}/genome.cds | faops size stdin | wc -l)
        if [[ $CDS == $GEN ]]; then
            echo check ok
        else
            echo check wrong
        fi
    '
# nothing wrong
```

### Protein processing

CDS always failed in grouping orthogroups correctly, so I chose protein files.

```bash
mkdir ~/data/symbio/PROTEIN
cd ~/data/symbio

cat info/genome.lst |
    parallel -j 6 --keep-order '
        if [ ! -f GENOMES/{}/genome.pep ]; then
            echo "==>{} without protein"
        else
            cp GENOMES/{}/genome.pep PROTEIN/{}.pep.fa
        fi
    '
# remove all species without pro seq in the genome.lst

conda activate genome

for file in $(ls PROTEIN)
do
    echo "==> ${file}"
    python /home/j/anaconda3/envs/genome/bin/primary_transcript.py PROTEIN/${file}
    echo
done

cd PROTEIN
ls primary_transcripts | wc -l
#148

mv primary_transcripts/ ..
cd ..

# make sure which file should change seq_id
# check seq_id could not be change
for file in $(ls primary_transcripts)
do
    id=$(cat primary_transcripts/${file} | grep '^>' | head -n 1)
    result=$(echo $id | grep '\s')
    if [[ $result != "" ]]; then
        echo "==>${file} name should change"
        echo $id
        echo
    fi
done
#persea_americana
#lemna_minor
#ceratophyllum_demersum
#cat primary_transcripts/<species>.pep.fa |
#    perl -nle '
#        print ">$1" if /^>.+?\|\|.+?\|\|.+?\|\|.+?\|\|(.+?)\|/;
#        print if /^[A-Z]/i;
#    ' > tmp && mv tmp primary_transcripts/<species>.pep.fa

for file in $(ls primary_transcripts)
do
    id=$(cat primary_transcripts/${file} | grep '^>' | head -n 1)
    result=$(echo $id | grep '\s')
    if [[ $result != "" ]]; then    
        cat primary_transcripts/${file} | perl -nle '
            print ">$1" if /^>(.+?)\s/;
            print if /^[A-Z]/i;
        ' \
        > tmp && mv tmp primary_transcripts/${file}
    fi
done
# run above code for sure again

# add species name before the id
for name in $(ls primary_transcripts | perl -p -e 's/\..+$//')
do
    faops size primary_transcripts/${name}.pep.fa |
        cut -f 1 |
        awk -v name=${name} '{print ($0"\t"name"_"$0)}' \
        > name.tmp
    faops replace primary_transcripts/${name}.pep.fa name.tmp tmp.fa \
    && mv tmp.fa primary_transcripts/${name}.pep.fa
done

# some protein files using . for termination codon
# change them to * for satisfying diamond
for file in $(ls primary_transcripts)
do
    SEQ=$(faops some -l 0 primary_transcripts/${file} \
            <(faops size primary_transcripts/${file}) \
            stdout | grep -v '^>')
    result=$(echo $SEQ | grep '\.')
    if [[ $result != "" ]]; then
        echo "==>${file} end is ."
    fi
done

# manually change them
for file in $(ls primary_transcripts)
do
    cat primary_transcripts/${file} |
        perl -nle 'if(/^>/){print;}else{$_ =~ s/\./\*/g; print $_;}' > tmp \
        && mv tmp primary_transcripts/${file}
done
# run above code for sure again

# upload to the workstation
rsync -avP primary_transcripts/ name@ip:jyq/data/symbio/
```

## BUSCO

Running BUSCO for checking the genomes or assemblies sequencing quality.

```bash
mkdir ~/data/symbio/BUSCO
cd ~/data/symbio

# check busco of every genome
bash scripts/busco.sh
```

## Identifying gene orthogroups via orthogroups

### OrthoFinder processing

```bash
# do following to transmit data to server
#rsync -avP /home/j/data/symbio/CDS name@ip:jyq/data/

# connect to server
# ssh name@ip

cd ~/jyq/data/symbio

# orthofinder in conda base env
conda activate

orthofinder -f ./primary_transcripts -og
# orthofinder [options] -f <dir>
# -og: Stop after inferring orthogroups

# cp results to local
#rsync -avP name@ip:jyq/data/symbio/primary_transcripts/OrthoFinder/Results_Oct20/Orthogroups ~/data/symbio/

# If hard limit, h > r already, then you just need to increase the soft limit:
ulimit -n 22004

# check the limit
ulimit -Sn

# continue after orthogroups inferred
orthofinder -fg Results_Oct20/ -M msa -X
# -fg <dir>: Start OrthoFinder from pre-computed orthogroups in <dir>
# -M <txt>: Method for gene tree inference. Options 'dendroblast' & 'msa' [Default = dendroblast]
# -X: Don't add species names to sequence IDs
```

### OrthoFinder results

Acquire results from the workstation using rsync and put in the dir `Orthogroups` (command not shown).

In Ceratophyllum demersum, there are proteins acquired from six box translation, which means that suspect repeat proteins should be removed from the results.

All protein names were recorded into a file named `info/rm_ceratophyllum_demersum.lst`.

```bash
cd ~/data/symbio/Orthogroups
mkdir groups

cat Orthogroups.tsv | cut -f 1 | sed 1d > ortho.lst

wc -l ortho.lst
#135923 ortho.lst
# Totally 135923 orthogroups were identified

# split each orthogroup into a tsv
cat ortho.lst |
    parallel -j 4 --keep-order '
        echo "==> {}"
        bash ../scripts/ortho_extract.sh {} ./groups/
    '

# here considering tsv-join for orthogroups identification
for ortho in $(cat ortho.lst)
do
    echo "==> ${ortho}"
    cat groups/${ortho}.csv |
        mlr --icsv --otsv cat |
        cut -f 2 |
        sed 1d |
        tr ',' '\n' |
        sed 's/^\s//' |
        sed '/^$/d' |
        awk -v ortho=${ortho} '{print ($0"\t"ortho)}' \
        >> group/all_pro_ortho.tsv
done

cat all_pro_ortho.tsv | wc -l
#5005513

cat all_pro_ortho.tsv |
    tsv-join -f ../info/rm_ceratophyllum_demersum.lst -k 1 -e -z \
    > tmp && mv tmp all_pro_ortho.tsv

cat all_pro_ortho.tsv | wc -l
#5005513

dos2unix Orthogroups.GeneCount.tsv

cat Orthogroups.GeneCount.tsv |
    tsv-select -H -f Orthogroup,Total \
    > ortho_gene.tsv
```

### Renaming and extracting

After OrthoFinder processing, those pep files contained `:` in their seq_ids will be converted to `_`. Protein extraction will fail if you do not change your seq_ids.

```bash
cd ~/data/symbio

rm info/rename_id.lst
# check those files whose seq_ids should be changed
for file in $(ls primary_transcripts)
do
    if grep -i -q ':' primary_transcripts/${file}; then
        echo "${file}" >> info/rename_id.lst
    fi
done
#anthoceros_angustus.pep.fa
#bonia_amplexicaulis.pep.fa
#calamus_simplicifolius.pep.fa
#cocos_nucifera.pep.fa
#daemonorops_jenkinsiana.pep.fa
#mesotaenium_endlicherianum.pep.fa
#olyra_latifolia.pep.fa
#phalaenopsis_equestris.pep.fa
#phyllostachys_edulis.pep.fa
#rhododendron_delavayi.pep.fa
#welwitschia_mirabilis.pep.fa

for file in $(cat info/rename_id.lst)
do
    cat primary_transcripts/${file} |
        perl -p -e 's/:/_/g' \
        > tmp && mv tmp primary_transcripts/${file}
done

# run again to check whether renaming succeed
for file in $(ls primary_transcripts)
do
    if grep -i -q ':' primary_transcripts/${file}; then
        echo "${file}"
    fi
done
# OK
```

## Identify RLK

RLK (receptor-like receptor) structures:

From N-terminal ECD (extracellular domain) - TMD (transmembrane domain) - KD (kinase domain) to C-terminal.

So domains should be formatted correctly for a true RLK.

### Identify domains among all proteins via `hmmscan`

```bash
cd ~/data/symbio
mkdir -p ~/data/symbio/DOMAIN/pfam

cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo {}
        hmmscan --cpu 2 -E 0.1 --domE 0.1 -o DOMAIN/pfam/{}.txt \
            ./HMM/PFAM/Pfam-A.hmm primary_transcripts/{}.pep.fa
    '

ls DOMAIN/pfam/*.txt |
    parallel -j 6 --keep-order '
        echo "==> {/.}"
        perl scripts/hmm_results.pl -i {} \
        > DOMAIN/pfam/{/.}.pfam.tsv
    '

ls DOMAIN/pfam/*.txt |
    parallel -j 12 --keep-order 'echo {/.}' \
    > DOMAIN/species.lst
```

### Extract kinase domain

All KDs were extracted according to `--iregex Domain:pk`, which means that other domains contained `pk` would also be included.

E-value cutoff <= 1e-4 for KDs - kinase domains are relatively more conserved than ECD.

```bash
cd ~/data/symbio/DOMAIN
mkdir -p ~/data/symbio/DOMAIN/KD

rm all_domains.tsv
# all domains were extracted
for file in $(ls pfam/*.pfam.tsv)
do
    cat ${file} |
        tsv-select -H -f Domain |
        sed 1d \
        >> all_domains.tsv
done

# basic info of domains identified
cat all_domains.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_domains.tsv

cat all_domains.tsv | wc -l
#19622
# all domains identified, including those domains scanned repeatedly from a location
# this will be considered later

# extract all domains contained pk with cutoff E-value <= 1e-4
# other domains would be included no matter E-value
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat pfam/{}.pfam.tsv |
            tsv-filter -H --iregex Domain:pk --le E_value:"1e-4" |
            sed 1d |
            tsv-select -f 1 |
            tsv-uniq \
            > KD/{}.lst
        cat pfam/{}.pfam.tsv |
            sed 1d |
            tsv-join -f KD/{/}.lst -k 1 |
            tsv-select -f 1,2,4,5,3 \
            > KD/{}.tsv
        rm KD/{}.lst
    '

rm all_KD.tsv
for file in $(ls KD/*.tsv)
do
    cat ${file} |
        tsv-select -f 2 |
        tsv-filter --iregex 1:pk \
        >> all_KD.tsv
done

cat all_KD.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_KD.tsv

rm all_KD_species.tsv
# all proteins with the pk domain among species
ls KD/*.tsv |
    parallel -j 12 --keep-order '
        cat {} |
            cut -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={/.} '\''{print (spe"\t"$0)}'\'' \
            >> all_KD_species.tsv
    '

# all proteins
cat all_KD_species.tsv |
    cut -f 2 |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#204137

# extract all protein sequences contained pk
mkdir -p ~/data/symbio/DOMAIN/seq_KD

ls KD/*.tsv |
    parallel -j 12 --keep-order '
        echo "==> {/.}"
        faops some -l 0 \
            ../primary_transcripts/{/.}.pep.fa \
            <(cat {} | cut -f 1 | tsv-uniq) \
            ./seq_KD/{/.}.fa
    '

# check extraction
wc -l ./seq_KD/*.fa | grep 'total' | perl -p -e 's/\s+(\d+).+$/$1\/2/' | bc
#204137
# extraction complete
```

### Predict transmembrane domain via `TMHMM2`

Manually upload files to [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) and copy results into a tsv file.

Although `selenium` could deal with this automatically, but there are always problems that I cannot solve. But still, a script will be saved in [scripts](https://github.com/jdasfd/symbio/tree/main/scripts). It will be completed later (maybe XD).

```bash
cd ~/data/symbio/DOMAIN

mkdir -p TMD COMBINE

cat species.lst |
    parallel -j 12 --keep-order '
        touch TMD/{}.tsv
    '
# manually doing this steps:
# upload pep.fa of each specie acquired from the previous step to TMHMM 2.0 websites
# copy all predicted results into a tsv named by species
# tsv-utils could only deal with those files in unix style (LF), convert them using dos2unix

dos2unix TMD/*.tsv

wc -l TMD/* | grep 'total'
#  204137 total
# alright

# retrieved potential pk proteins if they had at least a TMD
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        perl ../scripts/tmhmm_result.pl -i TMD/{}.tsv |
            tsv-select -f 1,4,2,3 \
            > TMD/{}.TMD.tsv
        '

cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat TMD/{}.tsv |
            cut -f 1 |
            tsv-join -f TMD/{}.TMD.tsv \
            -k 1 -e \
            > TMD/{}.not_TMD.lst
    '
```

### Extract all potential RLK

All proteins with characteristic that TMD appeared before KD (from N-terminal to C-terminal, sorted by start) were treated as potential targets of RLK.

After sorted by col3 (start pos), a script named `domain_order.pl` could just judge whether TMD showed before KD.

```bash
# Combine pfam results with TMD
# sort with the domain start pos
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat KD/{}.tsv |
            tsv-select -f 1,2,3,4 \
            > COMBINE/{}.tsv
        cat TMD/{}.TMD.tsv \
            >> COMBINE/{}.tsv
        cat COMBINE/{}.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            tsv-sort > COMBINE/{}.tmp.lst
        if [ -f COMBINE/{}.sort.tsv ]; then
            rm COMBINE/{}.sort.tsv
        fi
        for acc in $(cat COMBINE/{}.tmp.lst)
        do
            cat COMBINE/{}.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-sort -nk 3,3 \
                >> COMBINE/{}.sort.tsv
        done
        rm COMBINE/{}.tmp.lst
    '

# count for every specie
rm all_PK_pot_species.tsv
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat COMBINE/{}.sort.tsv |
            cut -f 1 |
            uniq |
            wc -l |
            awk -v spe={} '\''{print (spe"\t"$0)}'\'' \
            >> all_PK_pot_species.tsv
    '

# filter according to RLK structure
# that is TMD appeared before KD
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat COMBINE/{}.sort.tsv |
            perl ../scripts/domain_order.pl \
            > COMBINE/{}.PK.lst
        '

rm all_PK_true_species.tsv
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat COMBINE/{}.PK.lst |
            wc -l |
            awk -v spe={} '\''{print (spe"\t"$0)}'\'' \
            >> all_PK_true_species.tsv
    '

# combine two files
cat all_PK_pot_species.tsv |
    tsv-join -f all_PK_true_species.tsv -k 1 -a 2 |
    tsv-sort -k 1,1 |
    sed '1ispecies\tpot_PK\ttrue_PK' \
    > all_PK_species.tsv

rm all_PK_*_species.tsv
```

### Remove repeated domains

Because of redundant results acquired from `hmmscan` (using PFAM-A database with cutoff E-value 1e-1), a process should be adopted for reducing repeated domains.

The main principle - filtering according to E-value.

Protein domains of each protein were sorted by E-value (from less to more). The domains with smaller E-value were kept, and every repeated domains with larger E-value were removed.

```bash
cd ~/data/symbio/DOMAIN
mkdir -p UNIQ

# sort according to E_value from less to more
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        if [ -f UNIQ/{}.sort.tsv ]; then
            rm UNIQ/{}.sort.tsv
        fi
        for acc in $(cat COMBINE/{}.PK.lst)
        do
            cat pfam/{}.pfam.tsv |
                tsv-filter -H --str-eq QUERY:${acc} |
                sed 1d |
                tsv-select -f 1,2,3,4,5 |
                sort -gk 3,3 \
                >> UNIQ/{}.sort.tsv
        done
    '
# sort -g: use generic numerical value

cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat UNIQ/{}.sort.tsv |
            perl ../scripts/domain_uniq.pl |
            tsv-select -f 1,4,5,2 \
            > UNIQ/{}.uniq.tsv
    '

# combine TMD to uniq.tsv
cat species.lst |
    parallel -j 12 -keep-order '
        echo "==> {}"
        for acc in $(cat COMBINE/{}.PK.lst)
        do
            cat TMD/{}.TMD.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-select -f 1,3,4,2 \
                >> UNIQ/{}.uniq.tsv
        done
    '

# sort from start pos
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat UNIQ/{}.uniq.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            tsv-sort > UNIQ/{}.tmp.lst
        if [ -f UNIQ/{}.sort_uniq.tsv ]; then
            rm UNIQ/{}.sort_uniq.tsv
        fi
        for acc in $(cat UNIQ/{}.tmp.lst)
        do
            cat UNIQ/{}.uniq.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-sort -nk 2,2 \
                >> UNIQ/{}.sort_uniq.tsv
        done
        rm UNIQ/{}.tmp.lst
    '

rm all_KD_uniq.tsv
for species in $(cat species.lst)
do
    cat UNIQ/${species}.sort_uniq.tsv |
        tsv-select -f 4 |
        tsv-filter --iregex 1:pk \
        >> all_KD_uniq.tsv
done

cat all_KD_uniq.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_KD_uniq.tsv
```

### Pick pkinase domains only

All domains contained `pk` were collected. But RLK only contains the pkinase domain.

Those domains in `all_KD_uniq.tsv` were checked through [InterPro](https://www.ebi.ac.uk/interpro/), and true pkinase domain were kept.

- Pkinase
- PK_Tyr_Ser-Thr
- Pkinase_fungal

```bash
cd ~/data/symbio/DOMAIN
mkdir -p RLK

# all pkinase protein names
cat species.lst |
    parallel -j 12 --keep-order '
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-filter --or --str-eq 4:Pkinase \
                --str-eq 4:Pkinase_fungal \
                --str-eq 4:PK_Tyr_Ser-Thr |
            cut -f 1 |
            tsv-uniq \
            > RLK/{}.RLK.lst
    '

# all RLK that truly contained the pkinase domain
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-join -f RLK/{}.RLK.lst -k 1 \
            > RLK/{}.RLK.tsv
    '

echo -e "species\tRLK_num" > all_RLK_species.tsv
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat RLK/{}.RLK.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={} '\''{print (spe"\t"$0)}'\'' \
            >> all_RLK_species.tsv
    '

cat all_PK_species.tsv |
    tsv-join -H -f all_RLK_species.tsv -k species -a RLK_num \
    > all_count_species.tsv

perl ~/Scripts/fig_table/xlsx2csv.pl -f ../info/taxo_genom.xlsx |
    mlr --icsv --otsv cat |
    tsv-select -H -f SpeciesName,ComName |
    sed 1d |
    perl -nle 'print lc $_;' |
    sed '1ispecies\ttaxo' \
    > taxonomy.tsv

cat taxonomy.tsv |
    tsv-join -H -f all_count_species.tsv \
        -k species -a pot_PK,true_PK,RLK_num \
    > Species_result.tsv
```

### Count all ECDs

```bash
mkdir -p ~/data/symbio/group
cd ~/data/symbio

cat DOMAIN/RLK/*.RLK.lst > group/all_RLK.lst
cat group/all_RLK.lst | wc -l
#63888

# here considering tsv-join for orthogroups identification
for ortho in $(cat Orthogroups/ortho.lst)
do
    echo "==> ${ortho}"
    cat Orthogroups/groups/${ortho}.csv |
        mlr --icsv --otsv cat |
        cut -f 2 |
        sed 1d |
        tr ',' '\n' |
        sed 's/^\s//' |
        sed '/^$/d' |
        awk -v ortho=${ortho} '{print ($0"\t"ortho)}' \
        >> group/all_pro_ortho.tsv
done

cd ~/data/symbio/group

cat all_pro_ortho.tsv | wc -l

cat ../info/rm_ceratophyllum_demersum.lst |
    perl -nle '@array = split/\|/, $_; print "ceratophyllum_demersum_$array[8]";' \
    > rm_ceratophyllum_demersum.lst
#5005533

cat all_pro_ortho.tsv |
    tsv-join -f rm_ceratophyllum_demersum.lst -k 1 -e -z \
    > tmp && mv tmp all_pro_ortho.tsv

cat all_pro_ortho.tsv | wc -l
#5005513

cat all_RLK.lst |
    tsv-join -f all_pro_ortho.tsv -k 1 -a 2 \
    > RLK_ortho.tsv

cat RLK_ortho.tsv | wc -l
#63766

cat all_RLK.lst | tsv-join -f RLK_ortho.tsv -k 1 -e > RLK_without_ortho.lst

cat RLK_ortho.tsv | tsv-summarize -g 2 --count > ortho_count.tsv
wc -l ortho_count.tsv
#1486 ortho_count.tsv
```

### Picture all domains

```bash
cd ~/data/symbio/DOMAIN
mkdir -p PICTURE
mkdir -p SEQ

ls RLK/*.RLK.lst |
    sed 's/\.RLK\.lst$//' |
    parallel -j 12 -keep-order '
        echo "==> {/}"
        faops some -l 0 ../primary_transcripts/{/}.pep.fa \
            {}.RLK.lst SEQ/{/}.RLK.fa
        faops size SEQ/{/}.RLK.fa \
            > SEQ/{}.len.tsv
    '
```

## RNA-seq of AM symbiosis

```bash
mkdir ~/data/symbio/sra
cd ~/data/symbio/sra


```

```bash
for file in $(ls)
do
    hmmscan --cpu 8 -E 1e-4 --domE 1e-4 -o ${file}.txt \
    ../../symbio/HMM/PFAM/Pfam-A.hmm ${file}
done
```

```bash
perl ../scripts/tmhmm_result.pl -i TMD/actinidia_chinensis.tsv > actinidia_chinensis.RLK.lst

cat pfam/actinidia_chinensis.pfam.tsv | sed 1d | tsv-join -f actinidia_chinensis.RLK.lst -k 1 > actinidia_chinensis.RLK.tsv


```
